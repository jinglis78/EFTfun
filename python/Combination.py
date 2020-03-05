from ROOT import RooMultiVarGaussian, RooArgList, RooFormulaVar, RooArgSet, RooMinuit, TGraph, RooProdPdf, RooFit, RooGaussian, TFile, RooWorkspace, RooDataSet, RooStats, RooFit, RooCmdArg, TH2F, RooPoisson, TMatrixDSym
from ROOT.RooStats import ModelConfig,ProfileLikelihoodCalculator,HybridCalculator,FrequentistCalculator,HypoTestInverter,ProfileLikelihoodTestStat
import itertools,re
import Utils

class Combination:
    def __init__(self, theory=None, predictions=None,measurements=None,name='test',verbose=True,correlated_measurements=[]):

        self.name=name
        self.theory=theory
        self.verbose=verbose
        
        assert len(predictions)==len(measurements)
        for p,m in zip(predictions,measurements):
            if m.reco_level or m.getCov()==None:
                continue
            while m.getCov().Determinant()<10e-20:
                print 'Tiny covariance matrix determinant (',m.name,':',m.getCov().Determinant(),') rescaling measurement and prediction by 10'
                p.scale(10)
                m.scale(10)

        self.mvgaussians=[]
        self.poissons=[]
        self.corrected_predictions=[]
        self.data=RooArgSet()
        for prediction,measurement in zip(predictions,measurements):
            if measurement.reco_level:
                if measurement.name in correlated_measurements:
                    print 'Warning:',measurement.name,'is in correlated_measurements but this has no effect on reco level'
                for ibin,(p,n) in enumerate(zip(prediction.getPredictionFormulas(),measurement.central_values)):
                    self.data.add(n)
                    self.poissons+=[RooPoisson(measurement.name+"_bin"+str(ibin)+"_poisson",measurement.name+"_bin"+str(ibin)+"_poisson",n,p)]
            elif measurement.name not in correlated_measurements:
                muVec=RooArgList()
                for i,(p,c) in enumerate(zip(prediction.getPredictionFormulas(),measurement.getMeasurementCorrection())):
                    if c!=None:
                        self.corrected_predictions+=[RooFormulaVar("corrected_prediction_{}_{}".format(measurement.name,i),"@0-@1", RooArgList(p,c))]
                        muVec.add(self.corrected_predictions[-1])
                    else:
                        muVec.add(p)
                xVec=RooArgList()
                for x in measurement.central_values:
                    xVec.add(x)
                self.data.add(xVec)
                self.mvgaussians+=[RooMultiVarGaussian(measurement.name+"_mvgaussian",measurement.name+"_mvgaussian", xVec, muVec,measurement.getCov())]

        if len(correlated_measurements)>0:
            combined_muVec=RooArgList()
            combined_xVec=RooArgList()
            combined_nbins=sum([m.nbins for m in measurements if m.name in correlated_measurements])
            combined_cov=TMatrixDSym(combined_nbins)
            combined_ibin=0
            combined_measured=[]
            correlated_uncertainty_names=set()
            for m1 in measurements:
                for u in m1.uncertainties.values():
                    if u.uncorr or u.with_nuis_par:
                        continue
                    for m2 in measurements:
                        if m2==m1:
                            continue
                        if u.name in m2.uncertainties.keys() and not u.uncorr and not u.with_nuis_par:
                            correlated_uncertainty_names.add(u.name)
            correlated_uncertainties={}
            for u in correlated_uncertainty_names:
                correlated_uncertainties[u]=[]

            for prediction,measurement in zip(predictions,measurements):
                if measurement.name not in correlated_measurements:
                    continue
                for i,(p,c) in enumerate(zip(prediction.getPredictionFormulas(),measurement.getMeasurementCorrection())):
                    if c!=None:
                        self.corrected_predictions+=[RooFormulaVar("corrected_prediction_{}_{}".format(measurement.name,i),"@0-@1", RooArgList(p,c))]
                        combined_muVec.add(self.corrected_predictions[-1])
                    else:
                        combined_muVec.add(p)
                for x in measurement.central_values:
                    combined_xVec.add(x)
                combined_cov.SetSub(combined_ibin,combined_ibin,measurement.getCov(exclude=correlated_uncertainty_names))
                combined_ibin+=measurement.nbins
                for u in measurement.uncertainties.values():
                    if u.name in correlated_uncertainty_names:
                        correlated_uncertainties[u.name]+=[e*n for e,n in zip(u.getRelEffectSym(),measurement.getNominal())]
                combined_measured+=measurement.getNominal()

            for u in correlated_uncertainties.values():
                for i in range(combined_nbins):
                    for j in range(combined_nbins):
                        combined_cov[i,j]=combined_cov[i][j]+u[i]*u[j]
            self.data.add(combined_xVec)
            combined_corr=TMatrixDSym(combined_nbins)
            for i in range(combined_nbins):
                for j in range(combined_nbins):
                    combined_corr[i,j]=combined_cov[i][j]/(combined_cov[i][i]*combined_cov[j][j])**0.5
                    
            self.mvgaussians+=[RooMultiVarGaussian("combined_mvgaussian","combined_mvgaussian", combined_xVec, combined_muVec,combined_cov)]
                

        self.constraints=[]
        for theta in self.theory.thetas:
            constraint=RooGaussian("constraint_"+theta.GetName(),"constraint_"+theta.GetName(),theta ,RooFit.RooConst(.0) ,RooFit.RooConst(1.))
            self.constraints+=[constraint]
        pdfs=RooArgList()
        for mvgaussian in self.mvgaussians:
            pdfs.add(mvgaussian)
        for poisson in self.poissons:
            pdfs.add(poisson)
        for constraint in self.constraints:
            pdfs.add(constraint)
        self.pdf=RooProdPdf(self.name+'_pdf',self.name+'_pdf',pdfs)

        self.dataset=RooDataSet(self.name+'_data',self.name+'_data',RooArgSet(self.data))
        self.dataset.add(RooArgSet(self.data))
        self.nll=self.pdf.createNLL(self.dataset)

        self.minimizer=RooMinuit(self.nll)
        self.last_scan_range=None
        self.all_results=None
        self.robust=True


    def getFitResult(self):
        return self.minimizer.save()
               
    def profile1D(self,poi_name,scan_range=3,nsteps=101,inclOffset=False,extra_outfile=None):
        poi=self.theory.getCoeff(poi_name)
        result=[]
        print_level=-1
        # auto range
        if not isinstance(scan_range,tuple) and not isinstance(scan_range,list) :
            poiIsConst=poi.isConstant()
            poi.setConstant(False)
            self.minimize(minos=True)
            res=self.getFitResult()
            float_pars=res.floatParsFinal()
            for i in range(float_pars.getSize()):
                if float_pars[i].GetName()==poi_name:                    
                    scan_range=(scan_range*float_pars[i].getAsymErrorLo()+float_pars[i].getVal(),scan_range*float_pars[i].getAsymErrorHi()+float_pars[i].getVal())
            poi.setConstant(poiIsConst)
            self.last_scan_range=scan_range

        step=float(scan_range[1]-scan_range[0])/(nsteps-1)
        more_results={}
        percent=0
        for i in range(nsteps):
            if 100*i/nsteps>=percent:
                print percent,'%'
                percent+=10
            val=scan_range[0]+step*i
            poi.setVal(val)
            poiIsConst=poi.isConstant()
            poi.setConstant(True)
            strategy=0
            if self.robust:
                strategy=2
            if i==0:
                strategy=2
            result.append( (val,self.minimize(print_level=print_level,hesse=False,strategy=strategy)) )
            if extra_outfile!=None:
                res=self.getFitResult()
                float_pars=res.floatParsFinal()
                if i==0:
                    more_results[poi.GetName()]=[]
                    for ipar in range(float_pars.getSize()):
                        more_results[float_pars[ipar].GetName()]=[]
                more_results[poi.GetName()]+=[poi.getVal()]
                for ipar in range(float_pars.getSize()):
                    more_results[float_pars[ipar].GetName()]+=[float_pars[ipar].getVal()]
            poi.setConstant(poiIsConst)
        if extra_outfile!=None:
            self.all_results=more_results
            Utils.saveToTree(more_results,extra_outfile)
        
                        
        if not inclOffset:
            self.theory.reset()
            poiIsConst=poi.isConstant()
            poi.setConstant(False)
            minNLL=self.minimize(print_level=print_level)
            poi.setConstant(poiIsConst)
            result =[(r[0],r[1]-minNLL) for r in result]
        return result

    #TODO: unify with 1D
    def profileND(self,poi_names,scan_ranges=[5],nsteps=[21],inclOffset=False,extra_outfile=None):
        pois=[]
        for p in poi_names:
            pois+=[self.theory.getCoeff(p)]
        if len(nsteps)==1:
            nsteps=nsteps*len(poi_names)
        if len(scan_ranges)==1:
            scan_ranges=scan_ranges*len(poi_names)

        # auto range
        actual_scan_ranges=[]
        for the_poi,scan_range in zip(poi_names,scan_ranges):
            if not isinstance(scan_range,tuple) and not isinstance(scan_range,list) :
                poiIsConst=[]
                for p in pois:
                    poiIsConst+=[p.isConstant()]
                    p.setConstant(False)
                self.minimize()
                res=self.getFitResult()
                float_pars=res.floatParsFinal()
                for i in range(float_pars.getSize()):
                    if float_pars[i].GetName()==the_poi:
                        actual_scan_ranges+=[(-scan_range*float_pars[i].getError()+float_pars[i].getVal(),scan_range*float_pars[i].getError()+float_pars[i].getVal())]
                for isc,poi in zip(poiIsConst,pois):
                    poi.setConstant(isc)
            else:
                actual_scan_ranges+=[scan_range]           
        self.last_scan_range=actual_scan_ranges
        


        assert len(actual_scan_ranges)==len(poi_names)
        assert len(nsteps)==len(poi_names)

        scan_points=[]
        for sr,ns in zip(actual_scan_ranges,nsteps):
            scan_points.append([sr[0]+i*float(sr[1]-sr[0])/(ns-1) for i in range(ns)])
        grid=list(itertools.product(*scan_points))
        print_level=-1
        result=[]
        #TODO: don't store in memory
        more_results={}
        npoints=len(grid)
        print 'scanning',npoints,'points'
        percent=0
        for ipoint,gridpoint in enumerate(grid):
            if 100*ipoint/npoints>=percent:
                print percent,'%'
                percent+=10
            poiIsConst=[]
            for poi,val in zip(pois,gridpoint):
                poi.setVal(val)
                poiIsConst+=[poi.isConstant()]
                poi.setConstant(True)            
            strategy=0
            if self.robust:
                strategy=2
            if ipoint==0 or grid[ipoint][-1]<grid[ipoint-1][-1]:
                strategy=2
            result.append( gridpoint+(self.minimize(print_level=print_level,hesse=False,strategy=strategy),) )
            if extra_outfile!=None:
                res=self.getFitResult()
                float_pars=res.floatParsFinal()
                if ipoint==0:
                    for poi in pois:
                        more_results[poi.GetName()]=[]
                    for ipar in range(float_pars.getSize()):
                        more_results[float_pars[ipar].GetName()]=[]
                for poi in pois:
                    more_results[poi.GetName()]+=[poi.getVal()]
                for ipar in range(float_pars.getSize()):
                    more_results[float_pars[ipar].GetName()]+=[float_pars[ipar].getVal()]

            for pisc,poi in zip(poiIsConst,pois):
                poi.setConstant(pisc)
        if not inclOffset:
            self.theory.reset()
            for poi in pois:
                poi.setConstant(False)
            minNLL=self.minimize(print_level=print_level)
            for pisc,poi in zip(poiIsConst,pois):
                poi.setConstant(pisc)
            result =[r[:-1]+((r[-1]-minNLL),) for r in result]
        if extra_outfile!=None:
            Utils.saveToTree(more_results,extra_outfile)
            self.all_results=more_results
        return result

    def deactivateIrrelevantThetas(self):
        deactivated={}
        for c in self.theory.coeffs:
            if not c.isConstant() or c.getVal()!=0.:
                continue
            for theta in self.theory.thetas:
                if theta.isConstant():
                    continue
                if re.match(self.theory.makeMCstatName('.+',c.GetName(),'\d+'),theta.GetName()) or re.match(self.theory.makeMCstatName('.+',[c.GetName(),'.+'],'\d+'),theta.GetName()):
                    deactivated[theta]=theta.getVal()
                    theta.setVal(0.)
                    theta.setConstant()

        return deactivated


    def minimize(self,strategy=2,hesse=True,minos=False,print_level=None):
        deactivated=self.deactivateIrrelevantThetas()
        self.minimizer.setPrintLevel(1)
        if print_level!=None:
            self.minimizer.setPrintLevel(print_level)
        if not self.verbose or print_level<0:
            self.minimizer.setPrintLevel(-1)
            self.minimizer.setPrintEvalErrors(-1)
            self.minimizer.setVerbose(False)
        self.minimizer.setStrategy(strategy)
        self.minimizer.setErrorLevel(0.5)
        self.minimizer.migrad()
        if hesse:
            self.minimizer.hesse()
        if minos:
            self.minimizer.minos()
        for d in deactivated:
            d.setVal(deactivated[d])
        return self.nll.getVal()

    #TODO: other interval calculators, CLs.

    def fit(self):
        #TODO: print result even if quiet
        self.minimize(print_level=2,strategy=2,hesse=True)
        
    def calcProfileLikelihood(self,poi_names):
        if len(poi_names)>1:
            print 'can only do 1D intervals right now'
            return
        pois=RooArgSet()
        poiIsConst=[]
        for p in poi_names:            
            poi=self.theory.getCoeff(p)
            poiIsConst.append(poi.isConstant())
            poi.setConstant(False)
            pois.add(poi)
        pl = ProfileLikelihoodCalculator(self.dataset,self.pdf,pois)
        pl.SetConfidenceLevel(0.95)
        interval = pl.GetInterval()
        for poi in [pois.first()]:
            lowerLimit = interval.LowerLimit(poi)
            upperLimit = interval.UpperLimit(poi)
            print "95% CL interval on {0} is [{1},{2}]".format(poi.GetName(),lowerLimit,upperLimit)
            return lowerLimit,upperLimit


    def createWorkspace(self):
        self.workspace=RooWorkspace('w')
        getattr(self.workspace,'import')(self.pdf)
        getattr(self.workspace,'import')(self.dataset)

        modelconfig=ModelConfig(self.workspace)
        modelconfig.SetName(self.name+'_model')
        modelconfig.SetPdf(self.pdf)
        pois=RooArgSet()
        for c in self.theory.coeffs:
            pois.add(c)
        nps=RooArgSet()
        for t in self.theory.thetas:
            nps.add(t)

        modelconfig.SetParametersOfInterest(pois)
        modelconfig.SetObservables(RooArgSet(self.data))
        modelconfig.SetNuisanceParameters(nps)
        
        background_model =  modelconfig.Clone(self.name+'_model_background')
        oldval = pois.first().getVal()
        pois.first().setVal(0)
        background_model.SetSnapshot(pois)
        pois.first().setVal(oldval)
        getattr(self.workspace,'import')(modelconfig)
        getattr(self.workspace,'import')(background_model)

        return self.workspace

    def correlationHistPOIs(self):
        res=self.getFitResult()
        float_pars=res.floatParsFinal()
        cnames=self.theory.getCoeffNames()
        pois=[]
        for ipar in range(float_pars.getSize()):
            name=float_pars[ipar].GetName()
            if name in cnames:
                pois.append(name)
        if len(pois)<2:
            return None
        h=TH2F('correlationPOIs','correlation POIs',len(pois),0.,len(pois),len(pois),0,len(pois))
        for i,ipoi in enumerate(pois):
            for j,jpoi in enumerate(pois):
                h.SetBinContent(i+1,j+1,res.correlation(ipoi,jpoi))
        for i,ipoi in enumerate(pois):
            h.GetXaxis().SetBinLabel(i+1,ipoi)
            h.GetYaxis().SetBinLabel(i+1,ipoi)
        return h
        
            
        

    def pulls(self,thetas=True,pois=True,skipMCstat=True):
        pulls={}
        if self.verbose:
            print 'Minimizing for pulls'
        self.minimize(strategy=2,minos=True)
        res=self.getFitResult()
        float_pars=res.floatParsFinal()
        pars=[]
        if thetas:
            pars+=[t.GetName() for t in self.theory.thetas]
        if pois:
            pars+=[c.GetName() for c in self.theory.coeffs]
        for ipar in range(float_pars.getSize()):
            name=float_pars[ipar].GetName()
            if not name in pars:
                continue
            if skipMCstat and name.startswith('mcstat'):
                continue
            val=float_pars[ipar].getVal()
            hi=float_pars[ipar].getAsymErrorHi()
            lo=float_pars[ipar].getAsymErrorLo()            
            pulls[name]=(val,lo,hi)
        return pulls

    def pullEach(self,pois=True,thetas=False,skipMCstat=True):
        pulls={}
        pars=[]
        if thetas:
            pars+=[t for t in self.theory.thetas if not t.isConstant()]
        if pois:
            pars+=[c for c in self.theory.coeffs if not c.isConstant()]
        for thepar in pars:
            if skipMCstat and thepar.GetName().startswith('mcstat'):
                continue
            for par in pars:
                par.setVal(0.)
                par.setConstant()
            thepar.setConstant(False)
            if self.verbose:
                print 'Minimizing for pulls, with all POIs but',thepar.GetName(),'fixed at zero'
            self.minimize(strategy=2,minos=True)
            res=self.getFitResult()
            float_pars=res.floatParsFinal()
            for ipar in range(float_pars.getSize()):
                if float_pars[ipar].GetName()==thepar.GetName():
                    val=float_pars[ipar].getVal()
                    hi=float_pars[ipar].getAsymErrorHi()
                    lo=float_pars[ipar].getAsymErrorLo()            
                    pulls[thepar.GetName()]=(val,lo,hi)
                break
        for par in pars:
            par.setConstant(False)
        return pulls
