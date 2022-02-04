from ROOT import RooRealVar, TMatrixDSym, TH1F, RooArgList, RooFormulaVar
from array import array
import Utils
import Uncertainty
import re

class Measurement:
    def __init__(self, name, nbins=None, measured=None, covariance=None, per_bin_width=False,binning=None,reco_level=False,theory=None,blind=None,binlabels=None):
        self.name=name
        self.nbins=nbins
        if measured!=None:
            self.setMeasured(measured)
        if covariance!=None:
            self.setCov(covariance)
        else:
            self.covMatrix=None
        self.per_bin_width=per_bin_width
        self.binning=binning
        self.binlabels=binlabels
        self.rescaled=1.
        self.uncertainties={}
        self.hepdata=None
        self.reco_level=reco_level
        self.theory=theory
        self.blind=blind
        

    def setMeasured(self,measured):
        if self.nbins==None:
            self.nbins=len(measured)
        assert len(measured)==self.nbins
        if self.blind is None:
            self.blind=self.nbins*[1]
        if self.binlabels is None:
            self.binlabels=range(self.nbins)
        self.central_values = []
        for i,(x_meas,b,binlabel) in enumerate(zip(measured,self.blind,self.binlabels)):
            x=RooRealVar(self.name+"_meas_"+str(binlabel),self.name+"_meas_"+str(binlabel),0,1e9)
            if b==1:
                x.setVal(x_meas)
            else:
                x.setVal(0)
            x.setConstant()
            self.central_values.append(x)

    def getMeasurementCorrectionString(self):
        need_to_correct=False
        for u in self.uncertainties.values():
            if u.with_nuis_par:
                need_to_correct=True
                break
        if not need_to_correct:
            return None
        
        formulas=[]
        for i,x in enumerate(self.central_values):
            form='-{}+{}'.format(x.getVal(),x.getVal())
            for u in self.uncertainties.values():
                if  u.with_nuis_par:
                    s=u.getStr(i)
                    if s!=None:
                        form+='*'+s
            formulas+=[form]
        return formulas

    def getMeasurementCorrection(self):
        formulas=self.getMeasurementCorrectionString()
        if formulas==None:
            return [None]*self.nbins
        self.measurement_correction=[]
        if self.binlabels is None:
            self.binlabels=range(self.nbins)
        for binlabel,form in zip(self.binlabels,formulas):
            params=RooArgList()
            for t in self.theory.thetas:
                if re.match('.*\\b'+t.GetName()+'\\b.*',form):
                    params.add(t)
            for c in self.central_values:
                if re.match('.*\\b'+c.GetName()+'\\b.*',form):
                    params.add(c)
            self.measurement_correction+=[RooFormulaVar("meas_correction_{}_{}".format(self.name,binlabel),"meas_correction_{}_{}".format(self.name,binlabel),form,params)]

        return self.measurement_correction

    def setCov(self,cov):
        if self.nbins==None:
            self.nbins=len(cov)

        assert len(cov)==self.nbins
        for row in cov:
            len(cov)==self.nbins

        self.covMatrix = TMatrixDSym(self.nbins)
        for i in range(self.nbins):
            for j in range(self.nbins):
                self.covMatrix[i,j]=cov[i][j]
                if i==j and cov[i][j]==0:
                    self.covMatrix[i,j]=1
                    print "Warning, covariance matrix element 0"
    def setMeasuredFromCfg(self,alldata):
        fromFile=False
        for e in ['.yoda','.yaml','.root']:
            if e in alldata:
                fromFile=True
                break
        if not fromFile:
            x=[float(x) for x in alldata.split()]
        else:
            x=[]
            for data in alldata.split():
                if data.split(':')[0].endswith('.yoda'):
                    spl=data.split(':')
                    assert len(spl)==2
                    x+=self.getMeasuredFromYoda(spl[0],spl[1])
                elif data.split(':')[0].endswith('.root'):
                    spl=data.split(':')
                    assert len(spl)==2
                    x+=self.getMeasuredFromRoot(spl[0],spl[1])
                elif data.split(':')[0].endswith('.yaml'):
                    spl=data.split(':')
                    if len(spl)==1:
                        spl.append(0)
                    assert len(spl)<=2
                    self.hepdata=data
                    x+=self.getMeasuredFromYaml(spl[0],int(spl[1]))
        self.setMeasured(x)

    def setCovFromCfg(self,data):
        if data.split(':')[0].endswith('.yoda'):
            spl=data.split(':')
            assert len(spl)==2
            self.setCov(self.getMatrixFromYoda(spl[0],spl[1]))
        elif data.split(':')[0].endswith('.root'):
            spl=data.split(':')
            assert len(spl)==2
            self.setCov(self.getMatrixFromRoot(spl[0],spl[1]))
        elif data.endswith('.yaml'):
            self.setCov(self.getMatrixFromYaml(data))
        else:
            cov=[]
            for y in data.split(','):
                cov+=[[float(x) for x in y.split()]]
            if not ',' in data:
                cov=[[float(data)]]

            self.setCov(cov)

    def setErrAndCorrFromCfg(self,data_err,data_corr):        
        if data_err.split(':')[0].endswith('.yoda'):
            import yoda
            filename_err,histoname_err=[s.strip() for s in data_err.split(':')]
            hist=yoda.read(filename_err)[histoname_err]
            if isinstance(hist,yoda.Scatter2D):
                errs=[(p.yErrs[0]+p.yErrs[1])/2 for p in hist.points]
            elif isinstance(hist,yoda.Histo1D):
                if self.per_bin_width:
                    errs=[b.heightErr for b in hist.bins]
                else:
                    errs=[b.areaErr for b in hist.bins]
            else:
                pass
                #TODO: throw error
                
        elif data_err.split(':')[0].endswith('.root'):
            import ROOT
            filename_err,histoname_err=[s.strip() for s in data.split(':')]
            f=ROOT.TFile(filename_err)
            histo=f.Get(histoname_err)            
            if self.per_bin_width:
                errs=[histo.GetBinError(iBin) for iBin in range(1,histo.GetNbinsX()+1)]
            else:
                errs=[histo.GetBinError(iBin)/histo.GetBinWidth(iBin) for iBin in range(1,histo.GetNbinsX()+1)]
            f.Close()
        elif data_err.split(':')[0].endswith('.yaml'):
            import yaml
            filename_err=data_err.split(':')[0]
            if ':' in data_err and data_err.split(':')[1].lower()!='all':
                label_err=data_err.split(':')[1]
                errs=Uncertainty.getErrsFromYaml(self,filename_err,label_err,rel=False)
                errs=[(abs(x[0])+abs(x[1]))/2 for x in errs]
            else:
                errs=Uncertainty.getErrsFromYaml(self,filename_err,err_label=None,rel=False)
#                f=open(filename_err.split(':')[0])
#                vals=yaml.load(f)['dependent_variables'][0]['values']
#                assert len(vals)==self.nbins
#                err_labels=[]
#                for y in vals:
#                    for error in y['errors']:
#                        err_labels+=[error['label']]
#                err_labels=list(set(err_labels))
#                total_errs2=[0.]*self.nbins
#                for label_err in err_labels:
#                    this_errs=Uncertainty.getErrsFromYaml(self,filename_err,label_err,rel=False)
#                    for ierr,err in enumerate(this_errs):
#                        total_errs2[ierr]+=(abs(err[0])+abs(err[1]))**2/4.
#                errs=[x**0.5 for x in total_errs2]
#                f.close()
            assert len(errs)==self.nbins

        else:
            errs=[float(x) for x in data_err.split()]

        if data_corr.split(':')[0].endswith('.yoda'):
            import yoda
            filename_corr,histoname_corr=[s.strip() for s in data_corr.split(':')]
            corr=self.getMatrixFromYoda(filename_corr,histoname_corr)
        elif data_corr.split(':')[0].endswith('.root'):
            import ROOT
            filename_corr,histoname_corr=[s.strip() for s in data.split(':')]
            corr=self.getMatrixFromRoot(filename_corr,histoname_corr)
        elif data_corr.split(':')[0].endswith('.yaml'):
            corr=self.getMatrixFromYaml(*(data_corr.split(':')))
        elif data_corr.lower()=='none':
            corr=[]
            for i in range(self.nbins):
                corr.append([1 if j==i else 0 for j in range(self.nbins)])
        else:
            corr=[]
            for y in data_corr.split(','):
                corr+=[[float(x) for x in y.split()]]
            if not ',' in data_corr:
                corr+=[float(data_corr)]
        assert len(corr)==self.nbins
        assert len(errs)==self.nbins
        cov=[]
        for i in range(self.nbins):
            cov+=[[0.]*self.nbins]
            for j in range(self.nbins):
                cov[i][j]=errs[i]*errs[j]*corr[i][j]
        self.setCov(cov)

    def addUncertaintyFromCfg(self,name,data):
        if self.reco_level:
            raise Exception("At reco level, the measurement has no uncertainty, it's a number -- the uncertainty is in the prediction")
        u=Uncertainty.createUncertaintyFromCfg(self,name,data,reco_level=False)
        self.addUncertainty(u)
            
    def addUncertainty(self,uncertainty):
        if uncertainty==None:
            print 'WARNING: empty uncertainty (likely not found in histo file)'
            return

        assert uncertainty.nbins==self.nbins
        if uncertainty.name in self.uncertainties:
            raise Exception('Added same uncertainty twice '+uncertainty.name)
        if uncertainty.with_nuis_par:
            self.theory.addTheta(uncertainty.name)
        self.uncertainties[uncertainty.name]=uncertainty

    def getMeasuredFromYoda(self,filename,name):
        import yoda
        
        hist=yoda.read(filename)[name]
        if isinstance(hist,yoda.Scatter2D):
            points=sorted([(p.x,p.y) for p in hist.points] , key=lambda p: p[0])
            measured=[p[1] for p in points]
        elif isinstance(hist,yoda.Histo1D):
            if self.per_bin_width:
                measured=[b.height for b in hist.bins]
            else:
                measured=[b.area for b in hist.bins]
        return measured

    def getMeasuredFromRoot(self,filename,name):
        import ROOT
        f=ROOT.TFile(filename)
        histo=f.Get(name)
        if self.per_bin_width:
            contents=[histo.GetBinContent(iBin)/histo.GetBinWidth(iBin) for iBin in range(1,histo.GetNbinsX()+1)]
        else:
            contents=[histo.GetBinContent(iBin) for iBin in range(1,histo.GetNbinsX()+1)]
        f.Close()
        return contents

    def getMeasuredFromYaml(self,filename,index):
        import yaml
        f=open(filename)
        data=yaml.load(f)
        vals=data['dependent_variables'][int(index)]['values']
        return [float(y['value']) for y in vals]

    def getMatrixFromYoda(self,filename,name):
        import yoda
        scatter=yoda.read(filename)[name]
        points=sorted([(p.x,p.y,p.z) for p in scatter.points] , key=lambda p: p[0])
        matrix=[]
        for i in range(self.nbins):
            matrix+=[[]]
            for j in range(self.nbins):
                matrix[-1]+=[points[i*self.nbins+j]]
            matrix[-1]=sorted([p for p in matrix[-1]] , key=lambda p: p[1])
        actual_matrix=[]
        for row in matrix:
            actual_matrix+=[[p[2] for p in row]]
        return actual_matrix

    def getMatrixFromRoot(self,filename,name):
        import ROOT
        f=ROOT.TFile(filename)
        histo=f.Get(name)
        matrix=[]
        for i in range(1,histo.GetNbinsX()+1):
            matrix+=[[]]
            for j in range(1,histo.GetNbinsY()+1):
                matrix[-1]+=[histo.GetBinContent(i,j)]
        return matrix

    def getMatrixFromYaml(self,filename,index=0):
        import yaml
        f=open(filename)
        data=yaml.load(f)
        vals=data['dependent_variables'][index]['values']
        assert len(vals)==self.nbins*self.nbins
        matrix=[]
        for i in range(self.nbins):
            matrix+=[[]]
            for j in range(self.nbins):
                matrix[-1]+=[float(vals[i*self.nbins+j]['value'])]
        return matrix

    def scale(self,s):
        if self.covMatrix!=None:
            self.covMatrix*=s*s  
        for x in self.central_values:
            x.setVal(x.getVal()*s)
        self.rescaled*=s

    def plotPrefit(self, outfile):
        if self.binning!=None:
            h=TH1F(self.name+'_data',self.name+'_data',self.nbins,array('f',self.binning))
        else:
            h=TH1F(self.name+'_data',self.name+'_data',self.nbins,0,1)
        errs=self.getErr(incl_np=True)
        for i in range(self.nbins):
            h.SetBinContent(i+1,self.central_values[i].getVal()/self.rescaled)
            h.SetBinError(i+1,errs[i]/self.rescaled)
        outfile.cd()
        h.Write()
        return h

    def plotPostfit(self, outfile,fitresult):
        if self.binning!=None:
            h=TH1F(self.name+'_data',self.name+'_data',self.nbins,array('f',self.binning))
        else:
            h=TH1F(self.name+'_data',self.name+'_data',self.nbins,0,1)
        errs=self.getErr(incl_np=False)
        for i,f in enumerate(self.getMeasurementCorrection()):
            h.SetBinContent(i+1,(self.central_values[i].getVal()+f.getVal())/self.rescaled)
            h.SetBinError(i+1,f.getPropagatedError(fitresult)/self.rescaled)
            h.SetBinError(i+1,(errs[i]**2/self.rescaled**2+h.GetBinError(i+1)**2)**0.5)
        outfile.cd()
        h.Write()
        return h

    def getCov(self,incl_np=False,exclude=[],redo=True):
        if incl_np==False and len(exclude)==0 and not redo:
            if hasattr(self,'bufferedCov'):
                return self.bufferedCov
        if self.reco_level:
            raise Exception("At reco level, the covariance matrix should not be used")
        if self.covMatrix:
            cov=TMatrixDSym(self.covMatrix)
        else:
            cov=TMatrixDSym(self.nbins)
        for u in self.uncertainties.values():            
            if u.already_in_cov:
                cov-=u.getCov(self.getMeas(),incl_np=True)
        for u in self.uncertainties.values():
            if u.name in exclude:
                continue
            if incl_np or not u.with_nuis_par:
                cov+=u.getCov(self.getMeas(),incl_np=incl_np)
        if incl_np==False and len(exclude)==0:
            self.bufferedCov=cov
        return cov

    def getCorr(self,incl_np=False):
        m=TMatrixDSym(self.getCov(incl_np=incl_np))
        var=[m[i][i] for i in range(self.nbins)]
        for i in range(self.nbins):
            for j in range(self.nbins):
                m[i][j]=m[i][j]/var[i]**0.5/var[j]**0.5 if var[i]!=0 and var[j]!=0 else 0
        return m

    def getRelErr(self,incl_np=False):
        if self.reco_level:
            return [1/n**0.5 for n in self.getMeas()]
        cov=self.getCov(incl_np=incl_np)
        return [(cov[i][i])**0.5/self.getNominal()[i] if self.getNominal()[i]!=0. else 0. for i in range(self.nbins)]

    def getErr(self,incl_np=False):
        if self.reco_level:
            return [n**0.5 for n in self.getMeas()]
        cov=self.getCov(incl_np=incl_np)
        return [(cov[i][i])**0.5 for i in range(self.nbins)]

    def getMeas(self):
        return [x.getVal() for x in self.central_values]

    def getNominal(self):
        return [x.getVal() for x in self.central_values]


    def rebin(self,rebin):
        import numpy as np
        
        if not isinstance(rebin,list):
            n=self.nbins//int(rebin)
            combine=[rebin]*n
            if sum(combine)<self.nbins:
                combine+=[self.nbins-rebin*n]
        else:
            combine=rebin

        assert sum(combine)==self.nbins
        
        self.binlabels=None
        self.blind=[1 if b!=0 else 0 for b in Utils.mergeBins(self.blind,combine)]
        self.binning=None
        print 'WARNING Rebinning: binlabels, blinding, and bin boundaries currently not properly updated'

        if not self.covMatrix is None:
            cov_np=[]
            for i in range(self.nbins):
                cov_np.append([])
                for j in range(self.nbins):
                    cov_np[-1].append(self.covMatrix[i,j])
            cov_np=np.matrix(cov_np)
            w,v=np.linalg.eigh(cov_np)
            res=zip(w,v.getT())
            cov_np_new=np.matrix([[0.]*len(combine)]*len(combine))
            for l,v in res:
                dummy_u=[l**0.5*x for x in v.tolist()[0]]
                dummy_u_new=Utils.mergeBins(dummy_u,combine)
                cov_np_new+=np.outer(dummy_u_new,dummy_u_new)
            self.covMatrix=TMatrixDSym(len(combine))
            for i in range(len(combine)):
                for j in range(len(combine)):
                    self.covMatrix[i,j]=cov_np_new[i,j]

        new_central_values=Utils.mergeBins(self.getNominal(),combine)
        for u in self.uncertainties.values():
            newEffectUp=Utils.mergeBins([x[1] for x in u.getEffect(self.getNominal())],combine)
            newEffectDown=Utils.mergeBins([x[0] for x in u.getEffect(self.getNominal())],combine)
            newRelEffectUp=[e/n for e,n in zip(newEffectUp,new_central_values)]
            newRelEffectDown=[e/n for e,n in zip(newEffectDown,new_central_values)]
            u.releffect=zip(newRelEffectDown,newRelEffectUp)
            u.nbins=len(combine)
                          
        self.nbins=len(combine)
        self.setMeasured(new_central_values)




