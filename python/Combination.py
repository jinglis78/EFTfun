from ROOT import RooMultiVarGaussian, RooArgList, RooFormulaVar, RooArgSet, RooMinuit, TGraph, RooProdPdf, RooFit, RooGaussian, TFile, RooWorkspace, RooDataSet, RooStats, RooFit, RooCmdArg, TH2F, RooPoisson, TMatrixDSym, TH1D, TH2D, TCanvas, TRandom3
from ROOT.RooStats import ModelConfig,ProfileLikelihoodCalculator,HybridCalculator,FrequentistCalculator,HypoTestInverter,ProfileLikelihoodTestStat
import itertools,re,math
import Utils

class Combination:
    def __init__(self, theory=None, predictions=None,measurements=None,name='test',verbose=True,correlated_measurements=[],robust=1):

        self.name=name
        self.theory=theory
        self.verbose=verbose

        self.predictions=predictions
        self.measurements=measurements

        self.robust=int(robust)
        
        self.maxreprofile=3
      
        assert len(predictions)==len(measurements)
        for p,m in zip(predictions,measurements):
            if m.reco_level or m.getCov()==None:
                continue
            det=m.getCov().Determinant()
            if det<=1e-200:
                raise Exception("Covariance degenerate")
            while det>10e20 or det==float('+inf'):
                print 'Huge covariance matrix determinant (',m.name,':',det,') rescaling measurement and prediction by 0.1'
                p.scale(0.1)
                m.scale(0.1)
                det=m.getCov(redo=True).Determinant()
                print 'New determinant:', det
            while det<10e-20:
                print 'Tiny covariance matrix determinant (',m.name,':',det,') rescaling measurement and prediction by 10'
                p.scale(10)
                m.scale(10)
                det=m.getCov(redo=True).Determinant()
                print 'New determinant:', det


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
                    self.poissons+=[RooPoisson(measurement.name+"_bin"+str(ibin)+"_poisson",measurement.name+"_bin"+str(ibin)+"_poisson",n,p,True)]
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
            print 'combining covariance matrices of:'
            for m in measurements:
                print "{} ({} bins)".format(m.name,m.nbins)
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
                found=[]
                for u in measurement.uncertainties.values():
                    if u.name in correlated_uncertainty_names:
                        correlated_uncertainties[u.name]+=[e*n for e,n in zip(u.getRelEffectSym(),measurement.getNominal())]
                        found+=[u.name]
                for uname in correlated_uncertainty_names:
                    if not uname in found:
                        correlated_uncertainties[uname]+=[0.]*measurement.nbins
                    
                combined_measured+=measurement.getNominal()

            for uname in correlated_uncertainties:
                u=correlated_uncertainties[uname]
                print "adding correlated uncertainty",uname
                print [x/n if n!=0 else 0 for x,n in zip(u,combined_measured)]
                for i in range(combined_nbins):
                    for j in range(combined_nbins):
                        combined_cov[i,j]=combined_cov[i][j]+u[i]*u[j]
            self.data.add(combined_xVec)
            combined_corr=TMatrixDSym(combined_nbins)
            for i in range(combined_nbins):
                for j in range(combined_nbins):
                    combined_corr[i,j]=combined_cov[i][j]/(combined_cov[i][i]*combined_cov[j][j])**0.5
#            mdraw=TH2D(combined_corr)
#            print combined_corr
#            c=TCanvas()
#            mdraw.Draw('COLZ')
#            c.SaveAs('fullcorr.root')
#            raw_input()
            self.mvgaussians+=[RooMultiVarGaussian("combined_mvgaussian","combined_mvgaussian", combined_xVec, combined_muVec,combined_cov)]
            self.combined_cov=combined_cov
        else:
            self.combined_cov=None
    

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
        self.nll=self.pdf.createNLL(self.dataset,RooFit.Optimize(False))

        self.minimizer=RooMinuit(self.nll)
        unused_pars=theory.getCoeffNames()
        for pred in predictions:
            for l in pred.linear_effects:
                if l in unused_pars:
                    unused_pars.remove(l)
            for q in pred.quadratic_effects:
                if q in unused_pars:
                    unused_pars.remove(q)
            for s in pred.scaling_pars:
                if s in unused_pars:
                    unused_pars.remove(s)
        for p in unused_pars:
            print 'setting unused parameter',p,'constant'
            theory.getCoeff(p).setConstant()
            
        
        self.last_scan_range=None
        self.all_results=None

    def getFitResult(self):
        return self.minimizer.save()
               
    def profile1D(self,poi_name,scan_range=3,nsteps=101,inclOffset=False,extra_outfile=None,nretries=0):
        nretry=nretries+1
        poi=self.theory.getCoeff(poi_name)
        if poi.isConstant():
            raise Exception("Cannot fit "+poi_name+", is the parameter implemented for this measurement?")

        result=[]
        print_level=-1
        # auto range
        if not isinstance(scan_range,tuple) and not isinstance(scan_range,list) :
            poiIsConst=poi.isConstant()
            poi.setConstant(False)
            self.minimize(minos=True)
            res=self.getFitResult()            
            if res.status()==0:
                float_pars=res.floatParsFinal()
                for i in range(float_pars.getSize()):
                    if float_pars[i].GetName()==poi_name:
                        actual_scan_range=(scan_range*float_pars[i].getAsymErrorLo()+float_pars[i].getVal(),scan_range*float_pars[i].getAsymErrorHi()+float_pars[i].getVal())
               
            else:
                print 'using profilelikelihood method to determine scan range'
                prof_range=self.calcProfileLikelihood([poi_name])
                actual_scan_range=prof_range[0]*scan_range/2,prof_range[1]*scan_range/2
            poi.setConstant(poiIsConst)
            self.last_scan_range=actual_scan_range
        else:
            actual_scan_range=scan_range
        step=float(actual_scan_range[1]-actual_scan_range[0])/(nsteps-1)
        more_results={}
        percent=0
        i_step=0
        steps=range(nsteps)[nsteps/2:]+list(reversed(range(nsteps)[:nsteps/2]))
        for i_step,i in enumerate(steps):
#            if i==nsteps/2-1:
            self.theory.reset()
            if 100*i_step/nsteps>=percent:
                print percent,'%'
                percent+=10
            val=actual_scan_range[0]+step*i
            poi.setVal(val)
            poiIsConst=poi.isConstant()
            poi.setConstant(True)
            strategy=0
            if self.robust:
                strategy=2
            if i_step==0:
                strategy=2
            if self.robust==3 and len(result)>0 and result[-1][1]<min(x[1] for x in result)+2:
                oldpars={}
                for c in self.theory.getNonConstCoeffs():
                    oldpars[c.GetName()]=c.getVal()
                tobeat=(result[-1][1]+1,oldpars)
            else:
                tobeat=None

            nll=self.minimize(print_level=print_level,hesse=False,strategy=strategy,tobeat=tobeat)
            if not nll is None:
                result.append( (val,nll) )
            if extra_outfile!=None:
                res=self.getFitResult()
                float_pars=res.floatParsFinal()
                if i_step==0:
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
        if nretry<self.maxreprofile and min(x[1] for x in result)<-0.01:
            print 'found better minimum, redoing scan with this starting point'            
            y0=999
            for x in result:
                if x[1]<y0:
                    x0=x[0]
                    y0=x[1]
            self.theory.reset()
            self.theory.setInitialVal(poi_name,x0)
            if not isinstance(scan_range,tuple) and not isinstance(scan_range,list) :
                xlow=(x0-result[0][0])
                ylow=(result[0][1]-y0)
                xhigh=(result[-1][0]-x0)
                yhigh=(result[-1][1]-y0)
                newlow=xlow*(5/ylow)**0.5 if ylow>0 else None
                newhigh=xhigh*(5/yhigh)**0.5 if ylow>0 else None
                if newlow is None and  not newhigh is None:
                    newlow=newhigh
                elif newhigh is None and not newlow is None:
                    newhigh=newlow
                elif newlow is None and newhigh is None:
                    newhigh=actual_scan_range[1]-actual_scan_range[0]
                    newlow=actual_scan_range[1]-actual_scan_range[0]
                new_range=(x0-newlow,x0+newhigh)
            else:
                new_range=scan_range
            return self.profile1D(poi_name,new_range,nsteps,inclOffset,extra_outfile,nretries=nretry)
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
            for r in result:
                minNLL=min(minNLL,r[-1])
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


    def minimize(self,strategy=2,hesse=True,minos=False,print_level=2,nshuffle=0,tobeat=None):
       
        deactivated=self.deactivateIrrelevantThetas()
        self.minimizer.setPrintLevel(1)
        if print_level!=None:
            self.minimizer.setPrintLevel(print_level)
        if not self.verbose or print_level<0:
            self.minimizer.setPrintLevel(-1)
            self.minimizer.setPrintEvalErrors(-1)
            self.minimizer.setNoWarn()
            self.minimizer.setVerbose(False)
        self.minimizer.setStrategy(strategy)
        self.minimizer.setErrorLevel(0.5)
        
        nan=False
        if math.isnan(self.nll.getVal()):
            nan=True
        for c in self.theory.getNonConstCoeffs():
            if math.isnan(c.getVal()):
                nan=True
        if nan:
            'something is nan, this is bad, how did we get here???'
            self.theory.dumpPars()
            print 'nll:',self.nll.getVal()
            print 'setting everything to zero, hope that helps...'
            for c in self.theory.getNonConstCoeffs():
                c.setVal(0.)
            return None

        #print 'doing Migrad'
        self.minimizer.migrad()
        if self.robust==2:
            res=self.getFitResult()
            if (res.status()!=0 or res.covQual()!=3):
                #print 'doing Simplex'
                self.minimizer.simplex()
                #print 'doing Migrad (2)'
                self.minimizer.migrad()
        if self.robust==3 or self.robust==4:
            res=self.getFitResult()
            bad=(res.status()!=0 or res.covQual()!=3 or tobeat is not None and self.nll.getVal()>tobeat[0])
            if nshuffle>0:
                print 'new nll after shuffle',self.nll.getVal()
            if (bad and nshuffle<20) or self.robust==4 and nshuffle<5:
                if nshuffle==0 and bad:
                    print "WARNING: trouble finding good minimum"
                #self.theory.dumpPars()
                print 'nll:',self.nll.getVal()
                print "Shuffling pars to obtain better min, attempt",nshuffle+1
                if tobeat is not None and not self.nll.getVal()<tobeat[0]:
                    for x in tobeat[1]:
                        self.theory.setVal(x,tobeat[1][x])
                oldpars=self.shuffle(nshuffle+(int(1e9*self.nll.getVal()) if self.nll.getVal()!=float("inf") else 0 ) )
                self.minimize(nshuffle=nshuffle+1,tobeat=(self.nll.getVal()+0.05,oldpars),strategy=strategy,hesse=hesse,minos=minos,print_level=print_level)
            if not tobeat is None and self.nll.getVal()>tobeat[0]:
                for x in tobeat[1]:
                    self.theory.setVal(x,tobeat[1][x])
                    self.minimizer.migrad()
            
        if hesse:
            #print 'doing Hesse'
            self.minimizer.hesse()

        if minos:
            #print 'doing Minos'
            self.minimizer.minos()
        for d in deactivated:
            d.setVal(deactivated[d])
        return self.nll.getVal()

    def shuffle(self,seed):
        #TODO: print result even if quiet
        rand=TRandom3(seed)
        oldpars={}
        newpars={}
        for c in self.theory.getNonConstCoeffs():
            old=c.getVal()
            oldpars[c.GetName()]=old
            nll0=self.nll.getVal()
            step0=0.1
            c.setVal(old+step0)
            nllup=self.nll.getVal()
            c.setVal(old-step0)
            nlldown=self.nll.getVal()
            c.setVal(old)
            dnll=abs(min(nllup-nll0,nlldown-nll0))
            if dnll==0:
                step=1.
            else:
                step=step0/dnll**0.5
            x=rand.Exp(step)*(1-rand.Integer(2)*2)
            #print c.GetName(),":",c.getVal(),'-->',c.getVal()+x
            newpars[c.GetName()]=c.getVal()+x
        for c in self.theory.getNonConstCoeffs():
            c.setVal(newpars[c.GetName()])
        return oldpars


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


    def createWorkspace(self,poi_names):
        assert poi_names=='all' or isinstance(poi_names,list)
        self.workspace=RooWorkspace('w')
        getattr(self.workspace,'import')(self.pdf)
        getattr(self.workspace,'import')(self.dataset)

        modelconfig=ModelConfig(self.workspace)
        modelconfig.SetName(self.name+'_model')
        modelconfig.SetPdf(self.pdf)
        pois=RooArgSet()
        nps=RooArgSet()
        for c in self.theory.coeffs:
            if poi_names=='all' or c.GetName() in poi_names:
                pois.add(c)
            elif not c.isConstant():
                nps.add(c)

        for t in self.theory.thetas:
            nps.add(t)

        modelconfig.SetParametersOfInterest(pois)
        modelconfig.SetObservables(RooArgSet(self.data))
        modelconfig.SetNuisanceParameters(nps)
        #TODO: set global observables
        
        #background_model =  modelconfig.Clone(self.name+'_model_background')
        #if pois.getSize()!=0:
        #    oldval = pois.first().getVal()
        #    pois.first().setVal(0)
        #    background_model.SetSnapshot(pois)
        #    pois.first().setVal(oldval)
        #    getattr(self.workspace,'import')(background_model)
        #getattr(self.workspace,'import')(modelconfig)


        return self.workspace

    def correlationHistPOIs(self):
        res=self.getFitResult()
        float_pars=res.floatParsFinal()
        cnames=self.theory.getCoeffNames()
        pois=[]
        errors=[]
        for ipar in range(float_pars.getSize()):
            name=float_pars[ipar].GetName()
            if name in cnames:
                pois.append(name)
                errors.append(float_pars[ipar].getError())
                
        if len(pois)<2:
            return None,None
        h=TH2F('correlationPOIs','correlation POIs',len(pois),0.,len(pois),len(pois),0,len(pois))
        hcov=TH2F('covariancePOIs','covariance POIs',len(pois),0.,len(pois),len(pois),0,len(pois))
        for i,ipoi in enumerate(pois):
            for j,jpoi in enumerate(pois):
                h.SetBinContent(i+1,j+1,res.correlation(ipoi,jpoi))
                hcov.SetBinContent(i+1,j+1,res.correlation(ipoi,jpoi)*errors[i]*errors[j])
        for i,ipoi in enumerate(pois):
            h.GetXaxis().SetBinLabel(i+1,ipoi)
            h.GetYaxis().SetBinLabel(i+1,ipoi)
            hcov.GetXaxis().SetBinLabel(i+1,ipoi)
            hcov.GetYaxis().SetBinLabel(i+1,ipoi)
        return h,hcov
        
            
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

    def PCA(self,groups=[]):
        import numpy as np
        scaling_pars=[]
        for p in self.predictions:
            for s in p.scaling_pars:
                scaling_pars.append(s)


        self.theory.setConstant()
        for s in scaling_pars:
            self.theory.setVal(s,1.)
            self.theory.setFloating(s)
        self.fit()
        res=self.getFitResult()
        float_pars=res.floatParsFinal()
        errs={}
        corr={}
        for ipar in range(float_pars.getSize()):
            name=float_pars[ipar].GetName()
            errs[name]=float_pars[ipar].getError()
            corr[name]={}
            for jpar in range(float_pars.getSize()):
                name2=float_pars[jpar].GetName()
                corr[name][name2]=res.correlation(name,name2)
        
        m_cov=[]
        for iparx,parx in enumerate(scaling_pars):
            m_cov.append([])
            for ipary,pary in enumerate(scaling_pars):
                m_cov[-1].append(errs[parx]*errs[pary]*corr[parx][pary])
        m_cov=np.matrix(m_cov)

        #print 'Covariance Matrix'
        #print m_cov
        c_inv=m_cov.getI()
        
        all_coeffs=[c for c in self.theory.getCoeffNames() if not c in scaling_pars]
        remaining_coeffs=list(all_coeffs)
        for g in groups:
            for c in g:
                remaining_coeffs.remove(c)

        results=[]
        for coeffs in groups+[remaining_coeffs]:
            m_para=[]
            for c in coeffs:
                m_para.append([])
                for p in self.predictions:
                    m_para[-1]+=[0. if x is None else x for x in ([e/n if not n==0 else 0 for e,n in zip(p.getEffectFloat(c),p.getNominal())] if not p.getEffectFloat(c) is None else [0.]*p.nbins) ]
            m_para=np.matrix(m_para)

            #print 'Parametrization Matrix'
            #print m_para
            a=m_para*c_inv*m_para.getT()
            w,v=np.linalg.eigh(a)
            res=zip(w,v.getT())
            res.sort(reverse=True, key=lambda x : x[0])
            results+=[(r[0],sorted(zip(r[1].tolist()[0],coeffs),key=lambda x: abs(x[0]),reverse=True)) for r in res]
            
        for ir,r in enumerate(results):
            print 'EV {}:'.format(ir),Utils.round2(r[0],3),'({})'.format(Utils.round2(1/r[0]**0.5,3) if r[0]>0 else 999999.)
            print ' + '.join(["{}*{}".format(round(x[0],3),x[1]) for x in r[1] if abs(x[0])>0.1]).replace('+ -','- ')
            print
        resmap=[]
        for r in results:
            d={}
            norm=0.
            vmax=0.
            vmin=0.
            for v,c in r[1]:
                if abs(v)<0.01:
                    v=0
                norm+=v**2
                d[c]=v
                vmax=max(vmax,v)
                vmin=min(vmin,v)
            norm=norm**0.5
            if abs(vmin)>abs(vmax):
                norm*=-1
            for c in d:
                d[c]=d[c]/norm
            resmap+=[(r[0],d)]
        print '#'+' '.join(all_coeffs)
        for i,r in enumerate(resmap):
            print 'c'+str(i)+' = ',
            for c in all_coeffs:
                if c in r[1]:
                    print r[1][c],
                else:
                    print 0.,
            print
        n=len(scaling_pars)
        h_cov=TH2D('h_total_cov','Total Covariance',n,0,n,n,0,n)
        for i,c in enumerate(scaling_pars):
            h_cov.GetXaxis().SetBinLabel(i+1,c)
            h_cov.GetYaxis().SetBinLabel(i+1,c)
        h_corr=h_cov.Clone('h_total_corr')
        h_corr.SetTitle('Total Correlation')
        for i in range(n):
            for j in range(n):
                h_cov.SetBinContent(i+1,j+1,m_cov[i,j])
                h_corr.SetBinContent(i+1,j+1,m_cov[i,j]/(m_cov[i,i]*m_cov[j,j])**0.5)
        n=len(all_coeffs)
        #h_para=TH2D('h_para','Parametrization Matrix',n,0,n,n,0,n)
        h_ev=TH2D('h_ev','Parametrization Matrix',n,0,n,n,0,n)
        h_ev_abs=TH2D('h_ev_abs','Parametrization Matrix (abs)',n,0,n,n,0,n)
        for iy,r in enumerate(resmap):
            for ix,c in enumerate(all_coeffs):
                h_ev.GetXaxis().SetBinLabel(ix+1,c)
                h_ev_abs.GetXaxis().SetBinLabel(ix+1,c)            
                if c in r[1]:
                    h_ev.SetBinContent(ix+1,iy+1,r[1][c])
                    h_ev_abs.SetBinContent(ix+1,iy+1,abs(r[1][c]))
            h_ev.GetYaxis().SetBinLabel(iy+1,'c_{'+str(iy)+'} '+'({})'.format(Utils.round2(1/r[0]**0.5,3) if r[0]>0 else 999999))
            h_ev_abs.GetYaxis().SetBinLabel(iy+1,'c_{'+str(iy)+'} '+'({})'.format(Utils.round2(1/r[0]**0.5,3) if r[0]>0 else 999999))

        ev_effects=[]
        for i,r in enumerate(resmap):
            ev_effects.append([])
            s=1/r[0]**0.5 if r[0]>0. else 0.
            if s>1000:
                s=0
            for p in self.predictions:
                effect=[0.]*p.nbins
                for c in all_coeffs:
                    if c in r[1] and p.getEffectFloat(c) is not None:
                        s1=r[1][c]
                        ceffect=[s1*s*e/n if not n==0 else 0 for e,n in zip(p.getEffectFloat(c),p.getNominal())]
                        effect=[x+y for x,y in zip(effect,ceffect)]
                ev_effects[-1].append(effect)

                        
        hs_effects=[]
        for i,ev_effect in enumerate(ev_effects):
            hs_effects.append([])
            for p,effect in zip(self.predictions,ev_effect):
                h=TH1D('h_effect{}_{}'.format(i,p.name),'effect of EV {} on {}'.format(i,p.name),len(effect),0,1)
                for ie,e in enumerate(effect):
                    h.SetBinContent(ie+1,e)
                hs_effects[-1].append(h)
        return h_cov,h_corr,h_ev,h_ev_abs,hs_effects
