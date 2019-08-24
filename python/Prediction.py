from ROOT import RooFormulaVar, RooArgList, RooRealVar, TH1F
from math import log
from array import array
import os
import Utils

class Uncertainty:
    def __init__(self,nbins,affected,parname,releffect,interpol):
        self.nbins=nbins
        self.affected=affected
        self.parname=parname
        assert len(releffect)==self.nbins

        actual_releffect=[]
        for e in releffect:
            if isinstance(e,float):
                actual_releffect+=[(-e,e)]
            else:
                actual_releffect+=[e]
        assert len(actual_releffect)==self.nbins
        self.releffect=actual_releffect
        assert interpol in ['lin','exp']            
        self.interpol=interpol

        
    def getStr(self,iBin):
        if self.releffect[iBin]==0.:
            return None
        if isinstance(self.releffect[iBin],tuple) and self.releffect[iBin][0]==0 and self.releffect[iBin][1]==0:
            return None
        #TODO: define better interpolation functions
        if self.interpol=='exp':
            if abs(-self.releffect[iBin][0]/self.releffect[iBin][1]-1)<1e-2:
                return 'exp({0}*{1})'.format(self.releffect[iBin][1],self.parname)
            else:
                return '(exp(-{0}*{2})*({2}<0)+exp({1}*{2})*({2}>=0))'.format(self.releffect[iBin][0],self.releffect[iBin][1],self.parname) 
        if self.interpol=='lin':
            if abs(-self.releffect[iBin][0]/self.releffect[iBin][1]-1)<1e-2:
                return '(1+{0}*{1})'.format(self.releffect[iBin][1],self.parname)
            else:
                return '((1-{0}*{2})*({2}<0)+(1+{1}*{2})*({2}>=0))'.format(self.releffect[iBin][0],self.releffect[iBin][1],self.parname)

class Prediction:
    def __init__(self,name,theory,nbins=None, sm=None, per_bin_width=False,binning=None):
        self.default_range=(-1e9,1e9)
        self.name=name
        self.sm=sm
        self.theory=theory
        self.nbins=nbins
        self.uncertainties=[]
        self.normhisto=None
        self.per_bin_width=per_bin_width
        self.stat_unc_threshold=None
        self.extra_scaling=None
        self.binning=binning
        self.linear_effects={}
        self.quadratic_effects={}
        self.coeff_names=[]
        self.rescaled=1.

    
    def getLinearEffect(self,coeff_name,iBin):
        x=self.linear_effects[coeff_name][iBin]
        if not isinstance(x,tuple):
            effect=x
        else:
            effect=x[0]
            error=x[1]
            if self.stat_unc_threshold!=None:
                th=self.stat_unc_threshold
            elif self.theory.stat_unc_threshold!=None:
                th=self.theory.stat_unc_threshold
            if th>0 and effect!=0 and abs(error/effect)>th:
                self.addUncertainty(self.theory.makeMCstatName(self.name,coeff_name,iBin),[0.]*iBin+[error]+[0.]*(self.nbins-iBin-1),affected=coeff_name,interpol='lin')
        
        for u in self.uncertainties:
            if coeff_name in u.affected:
                ustr=u.getStr(iBin)
                if ustr!=None:
                    effect=str(effect)+'*'+ustr
        return effect


    def getQuadraticEffect(self,coeff_name1,coeff_name2,iBin):
        if coeff_name1<coeff_name2:
            c1=coeff_name1
            c2=coeff_name2
        else:
            c1=coeff_name2
            c2=coeff_name1
        if c1 in self.quadratic_effects and c2 in self.quadratic_effects[c1]:
            x=self.quadratic_effects[c1][c2][iBin]
        else:
            x=0.
        if not isinstance(x,tuple):
            effect=x
        else:
            effect=x[0]
            error=x[1]
            if self.stat_unc_threshold!=None:
                th=self.stat_unc_threshold
            elif self.theory.stat_unc_threshold!=None:
                th=self.theory.stat_unc_threshold
            if th>0 and effect!=0 and abs(error/effect)>th:
                self.addUncertainty('mcstat_'+self.name+'_'+coeff_name+'_quad_bin'+str(iBin),[0.]*iBin+[error]+[0.]*(self.nbins-iBin-1),affected=coeff_name,interpol='lin')

        for u in self.uncertainties:
            if c1+','+c2 in u.affected or c2+','+c1 in u.affected:
                ustr=u.getStr(iBin)
                if ustr!=None:
                    effect=str(effect)+'*'+ustr
        return effect
        

    def getPredictionStrings(self):
        eftprediction = []
        #TODO: tuple
        for iBin in range(self.nbins):

            formula=str(self.sm[iBin])
            for u in self.uncertainties:
                if 'sm' in u.affected:
                    ustr=u.getStr(iBin)
                    if ustr!=None:
                        formula+='*'+ustr

            for coeff_name in self.coeff_names:
                if coeff_name in self.linear_effects:
                    if self.getLinearEffect(coeff_name,iBin)!=0:
                        formula+="+{}*{}".format(coeff_name,self.getLinearEffect(coeff_name,iBin))

            for coeff_name in self.coeff_names:
                if coeff_name in self.quadratic_effects:
                    if self.getQuadraticEffect(coeff_name,coeff_name,iBin)!=0:
                        formula+="+{}*{}*{}".format(coeff_name,coeff_name,self.getQuadraticEffect(coeff_name,coeff_name,iBin))

            for coeff_name1 in self.coeff_names:
                if not coeff_name1 in self.quadratic_effects:
                    continue
                for coeff_name2 in self.coeff_names:
                    if not coeff_name2 in self.quadratic_effects:
                        continue
                    if coeff_name1==coeff_name2:
                        break
                    if self.getQuadraticEffect(coeff_name1,coeff_name2,iBin):
                        formula+="+{}*{}*{}".format(coeff_name1,coeff_name2,self.getQuadraticEffect(coeff_name1,coeff_name2,iBin))
            if self.extra_scaling!=None:
                formula="({})*({})".format(self.extra_scaling[iBin],formula)
            eftprediction+=[formula]

        sys_scaling=['']*self.nbins
        for u in self.uncertainties:
            if u.affected=='all':
                for iBin in range(self.nbins):
                    ustr=u.getStr(iBin)
                    if ustr!=None:
                        sys_scaling[iBin]+='*'+ustr

        predictions=['('+e+')'+s for e,s in zip(eftprediction,sys_scaling)]
        return predictions                   

    def getPredictionFormulas(self):
        self.predictions = []
        formulas=self.getPredictionStrings()
        for iBin,formula in enumerate(formulas):
            params=RooArgList()
            for c in self.theory.coeffs:
                params.add(c)
            for t in self.theory.thetas:
                params.add(t)
            #TODO: not allowed to be negative
            bin_prediction=RooFormulaVar("pred_{}_{}".format(self.name,iBin),"pred_{}_{}".format(self.name,iBin),formula,params)
            self.predictions.append(bin_prediction)
        return self.predictions

    def plotPrefit(self, outfile):
        if self.binning!=None:
            h=TH1F(self.name+'_prefit',self.name+'_prefit',self.nbins,array('f',self.binning))
        else:
            h=TH1F(self.name+'_prefit',self.name+'_prefit',self.nbins,0,1)
        for i,f in enumerate(self.getPredictionFormulas()):
            h.SetBinContent(i+1,f.getVal()/self.rescaled)
            uperror2=0.
            downerror2=0.
            for t in self.theory.thetas:
                if not t.isConstant():
                    oldval=t.getVal()
                    t.setVal(0)
                    nom=f.getVal()
                    t.setVal(1)
                    err1=f.getVal()-nom
                    t.setVal(-1)
                    err2=f.getVal()-nom
                    t.setVal(oldval)
                    if err1>0 and err1>err2:
                        uperror2+=err1**2
                    if err1<0 and err1<err2:
                        downerror2+=err1**2
                    if err2>0 and err2>err1:
                        uperror2+=err2**2
                    if err2<0 and err2<err1:
                        downerror2+=err2**2
            uperror=uperror2**0.5
            downerror=downerror2**0.5
            h.SetBinError(i+1,(uperror/2+downerror/2)/self.rescaled)
                   
        outfile.cd()
        h.Write()
        return h

        
    def plotPostfit(self, outfile, fitresult):
        hn=self.name+'_postfit'
        if self.binning!=None:
            h=TH1F(hn,hn,self.nbins,array('f',self.binning))
        else:
            h=TH1F(hn,hn,self.nbins,0,1)
        for i,f in enumerate(self.getPredictionFormulas()):
            h.SetBinContent(i+1,f.getVal()/self.rescaled)
            h.SetBinError(i+1,f.getPropagatedError(fitresult)/self.rescaled)
            
        outfile.cd()
        h.Write()
        return h

    def setSM(self,sm,rescale=1.):
        if self.nbins==None:
            self.nbins=len(sm)
        assert len(sm)==self.nbins
        self.sm=[x*rescale for x in sm]
        return True

    def setSMfromCfg(self,data):
        if Utils.isFile(data):
            spl=data.split(':')
            assert len(spl)==2
            return self.setSMfromFile(spl[0],spl[1])
        else:
            return self.setSM([float(x) for x in data.split()])
        
    def setSMfromFile(self,filename,histoname,rescale=1.):
        return self.setSM(self.readFile(filename,histoname)[0],rescale)

    def setNormHisto(self,normhisto,rescale=1.):
        if self.nbins==None:
            self.nbins=len(normhisto)
        assert len(normhisto)==self.nbins
        self.normhisto=[x*rescale for x in normhisto]        
        return True

    def setNormHistoFromCfg(self,data):
        if Utils.isFile(data):
            spl=data.split(':')
            assert len(spl)==2
            return self.setNormHistoFromFile(spl[0],spl[1])
        else:
            return self.setNormHisto([float(x) for x in data.split()])
        
    def setNormHistoFromFile(self,filename,histoname,rescale=1.):
        return self.setNormHisto(self.readFile(filename,histoname)[0],rescale)
        
    def setEffect(self,power,coeff_names,effect,uncertainty=None,rescale=1.,renorm=False,add=False,isrel=False):
        if uncertainty==None:
            uncertainty=[0.]*self.nbins
        if len(effect)==1:
            effect=effect*self.nbins
        if not isinstance(coeff_names,list):
            assert power==1
            coeff_names=[coeff_names]
        for c in coeff_names:
            if not c in self.coeff_names:
                self.coeff_names+=[c]        
        assert len(coeff_names) == power
        assert len(effect)==self.nbins
        assert len(uncertainty)==self.nbins

        actual_effect=[]
        for x,u in zip(effect,uncertainty):
            actual_effect+=[(x*rescale,u*rescale)]
        if isrel:
            actual_effect=[(s*e[0],s*e[1]) for s,e in zip(self.sm,actual_effect)]
        if renorm:
            actual_effect=[(s*(e[0]/n),s*(e[1]/n)) for s,e,n in zip(self.sm,actual_effect,self.normhisto)]

        if power==1:
            if add:
                self.linear_effects[coeff_names[0]]=[a+b for a,b in zip(self.linear_effects[coeff_names[0]],actual_effect)]
                #TODO: will also add a quadratic effect!
            else:
                self.linear_effects[coeff_names[0]]=actual_effect
        elif power==2:
            if coeff_names[0]<coeff_names[1]:
                c0,c1=coeff_names[0],coeff_names[1]
            else:
                c1,c0=coeff_names[0],coeff_names[1]
            if not c0 in self.quadratic_effects:
                self.quadratic_effects[c0]=dict()
            if add:
                self.quadratic_effects[c0][c1]=[a+b for a,b in zip(self.quadratic_effects[c0][c1],actual_effect)]
            else:
                self.quadratic_effects[c0][c1]=actual_effect
        else:
            pass #TODO: error
        return True
    
    def setEffectFromFile(self,power,coeff_names,filename,histoname,rescale=1.,renorm=False,add=False):
        if not os.path.exists(filename.strip()):
            return False
        else:
            actual_effect,actual_uncertainty=self.readFile(filename,histoname)
            self.setEffect(power,coeff_names,actual_effect,actual_uncertainty,rescale,renorm,add=add)
            return True

    def setEffectFromCfg(self,power,coeff_names,data,add=False):
        if '{' in data:
            options=set(o.lower() for o in data[data.find('{')+1:data.find('}')].split())
        else:
            options=set()
        data=data[data.find('}')+1:]
        isrel=False
        if 'abs' in options:
            isrel=False
            options.remove('abs')
        elif 'rel' in options:
            isrel=True
            options.remove('rel')
        renorm=None
        if 'renorm' in options:
            renorm=True
            options.remove('renorm')
        elif 'keepnorm' in options:
            renorm=False
            options.remove('keepnorm')
        if len(options)>0:
            raise Exception('unknow/incompatible options in effect '+data+' : '+' '.join(list(options)))

        actual_coeff_names=coeff_names.split()

        if Utils.isFile(data):
            filename,histoname=data.split(':')
            return self.setEffectFromFile(power,actual_coeff_names,filename,histoname,rescale=1.,renorm=renorm,add=add)
                
        else:
            effect=[float(d) for d in data.split()]
            return self.setEffect(power,actual_coeff_names,effect,uncertainty=None,rescale=1.,renorm=renorm,add=add,isrel=isrel)

    def setScalingFromCfg(self,data):
        if Utils.isFile(data):
            scaling=readFile(data)[0]
        else:
            scaling=data
            for c in ['+','-','*','/']:
                while True:
                    replaced=scaling.replace(c+' ',c).replace(' '+c,c)
                    if replaced==scaling:
                        break
                    else:
                        scaling=replaced
            scaling=scaling.split()
        self.setScaling(scaling)

    def setScaling(self,scaling):
        if len(scaling)==1:
            scaling=scaling*self.nbins
        assert len(scaling)==self.nbins
        self.extra_scaling=scaling

    def readFile(self,filename,histoname=None):
        if histoname==None:
            assert ":" in filename
            filename=filename.split(':')[0]
            histoname=filename.split(':')[1]
        filename=filename.strip()
        histoname=histoname.strip()
        if filename.lower().endswith('.root'):
            return self.readRoot(filename,histoname)
        elif filename.lower().endswith('.yoda'):
            return self.readYoda(filename,histoname)
        else:
            pass
        #TODO: error
        
    def readYoda(self,filename,histoname):
        import yoda
        histos=yoda.read(filename)
        hist=histos[histoname]
        if self.per_bin_width:
            contents=[b.height for b in hist.bins]
            uncertainty=[b.heightErr for b in hist.bins]
        else:
            contents= [b.area for b in hist.bins]
            uncertainty=[b.areaErr for b in hist.bins]
        return contents,uncertainty

    def readRoot(self,filename,histoname):
        import ROOT
        f=ROOT.TFile(filename)
        histo=f.Get(histoname)
        if self.per_bin_width:
            contents=[histo.GetBinContent(iBin)/histo.GetBinWidth(iBin) for iBin in range(1,histo.GetNbinsX()+1)]
            uncertainty=[histo.GetBinError(iBin)/histo.GetBinWidth(iBin) for iBin in range(1,histo.GetNbinsX()+1)]
        else:
            contents=[histo.GetBinContent(iBin) for iBin in range(1,histo.GetNbinsX()+1)]
            uncertainty=[histo.GetBinError(iBin) for iBin in range(1,histo.GetNbinsX()+1)]
        return contents,uncertainty
   
    def addUncertaintyFromCfg(self,name,data):
        if '{' in name:
            affected=set(o.strip() for o in name[name.find('{')+1:name.find('}')].split())
            if 'all' in affected:
                affected='all'
            thetaname=name[name.find('}')+1:].strip()
        else:
            thetaname=name
            affected='all'

        if '{' in data:
            options=set(o.lower() for o in data[data.find('{')+1:data.find('}')].split())
        else:
            options=[]
        data=data[data.find('}')+1:]
        interpol=None
        for o in options:
            if o in ['lin','exp']:
                interpol=o
                options.remove(o)
                break
        if 'abs' in options:
            isrel=False
            options.remove('abs')
        elif 'rel' in options:
            isrel=True
            options.remove('rel')
        renorm=None
        if 'renorm' in options:
            renorm=True
            options.remove('renorm')
        elif 'keepnorm' in options:
            renorm=False
            options.remove('keepnorm')
        if len(options)>0:
            raise Exception('unknow/incompatible options in uncertainty '+data+' : '+' '.join(list(options)))
        releffect=None
        if Utils.isFile(data):
            spl=data.split(',')
            assert len(spl)==2
            down_filename,down_histoname=spl[0].split(':')
            up_filename,up_histoname=spl[1].split(':')
            return self.addUncertaintyFromFile(thetaname,down_filename,up_filename,down_histoname,up_histoname,affected=affected,interpol=interpol,renorm=renorm)
        else:
            effect=[]
            for e in data.split():
                if ',' in e:
                    effect+=[ ( float(e.split(',')[0]) , float(e.split(',')[1]) ) ]
                else:
                    effect+=[ ( -float(e), float(e)) ]
            if len(effect)==1:
                effect=effect*self.nbins
            if isrel:
                releffect=effect
            else:
                releffect=[(e[0]/s,e[1]/s) for e,s in zip(effect,self.sm)]
            if renorm:
                releffect=[(e[0]*s/n,e[1]*s/n) for s,e,n in zip(self.sm,effect,self.normhisto)]
            return self.addUncertainty(thetaname,releffect=releffect,affected=affected,interpol=interpol)

    def addUncertaintyFromFile(self,name,down_filename,up_filename,down_histoname,up_histoname,affected='all',interpol='exp',renorm=False):
        #TODO: one-sided uncertainties
        #TODO: fix
        if down_filename==None:
            down=[0.]*self.nbins
        else:
            down=self.readFile(down_filename,down_histoname)[0]
        if up_filename==None:
            up=[0.]*self.nbins
        else:
            up=self.readFile(up_filename,up_histoname)[0]
        assert len(up)==self.nbins
        assert len(down)==self.nbins
        if renorm:
            up=[u/n-1 for u,n in zip(up,self.normhisto)]
            down=[d/n-1 for d,n in zip(down,self.normhisto)]
        else:
            up=[u/n-1 for u,n in zip(up,self.sm)]
            down=[d/n-1 for d,n in zip(down,self.sm)]
        releffect=zip(up,down)
        return self.addUncertainty(name,releffect,affected,interpol)

    def addUncertainty(self,parname,releffect,affected='all',interpol='exp'):
        if not isinstance(releffect,list):
            releffect=[releffect]*self.nbins
        assert len(releffect)==self.nbins
        self.uncertainties.append(Uncertainty(self.nbins,affected,parname,releffect,interpol))
        self.theory.addTheta(parname)
        return True

    def scale(self,s):
        self.sm=[x*s for x in self.sm]
        for k in self.linear_effects:
            for iBin in range(self.nbins):
                self.linear_effects[k][iBin]=(self.linear_effects[k][iBin][0]*s,self.linear_effects[k][iBin][1]*s)
        for k in self.quadratic_effects:
            for kk in self.quadratic_effects[k]:
                for iBin in range(self.nbins):
                    self.quadratic_effects[k][kk][iBin]=(self.quadratic_effects[k][kk][iBin][0]*s,self.quadratic_effects[k][kk][iBin][1]*s)
        if self.normhisto!=None:
            self.normhisto=[x*s for x in self.normhisto]
        self.rescaled*=s
