from ROOT import RooFormulaVar, RooArgList, RooRealVar, TH1F, TH2F
from math import log,sqrt
from array import array
import re
import os
import Utils
import Uncertainty
from Background import Background

class Prediction:
    def __init__(self,name,theory,nbins=None, sm=None, per_bin_width=False,binning=None,blind=None,binlabels=None):
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
        self.binlabels=binlabels
        self.blind=blind
        self.linear_effects={}
        self.quadratic_effects={}
        self.coeff_names=[]
        self.rescaled=1.
        self.backgrounds={}
        self.scaling_pars=[]

    #TODO: reduce duplication
    def getEffectFloat(self,coeffs,withUnc=False):
        if isinstance(coeffs,str):
            coeffs=[coeffs]
        if len(coeffs)==2:
            coeff_name1=coeffs[0]
            coeff_name2=coeffs[1]
            if coeff_name1<coeff_name2:
                c1=coeff_name1
                c2=coeff_name2
            else:
                c1=coeff_name2
                c2=coeff_name1
            if c1 in self.quadratic_effects and c2 in self.quadratic_effects[c1]:
                if not withUnc:
                    return [x[0] for x in self.quadratic_effects[c1][c2]]
                return self.quadratic_effects[c1][c2]
            else:
                return None
        elif len(coeffs)==1:
            c1=coeffs[0]
            if c1 in self.linear_effects:
                if not withUnc:
                    return [x[0] for x in self.linear_effects[c1]]
                return self.linear_effects[c1]
            else:
                return None
        else:
            return None
    
    def getLinearEffect(self,coeff_name,iBin,incl_syst=True):
        x=self.linear_effects[coeff_name][iBin]
        if not isinstance(x,tuple):
            effect=x
            mcstatname=None
        else:
            effect=x[0]
            error=x[1]
            if self.stat_unc_threshold!=None:
                th=self.stat_unc_threshold
            elif self.theory.stat_unc_threshold!=None:
                th=self.theory.stat_unc_threshold
            mcstatname=self.theory.makeMCstatName(self.name,coeff_name,iBin)
            if th>0 and effect!=0 and abs(error/effect)>th and not mcstatname in self.theory.getThetaNames():
                u=Uncertainty.Uncertainty(nbins=self.nbins,name=mcstatname,affected=coeff_name,interpol='lin',
                                releffect=[0.]*iBin+[abs(error/effect)]+[0.]*(self.nbins-iBin-1))
                self.addUncertainty(u)
            if self.getUnc(mcstatname):
                ustr=self.getUnc(mcstatname).getStr(iBin)
                effect=str(effect)+'*'+ustr                
        if incl_syst:
            for u in self.uncertainties:
                if u.name==mcstatname:
                    continue
                if coeff_name in u.affected:
                    ustr=u.getStr(iBin)
                    if ustr!=None:
                        effect=str(effect)+'*'+ustr
        return effect

    def getQuadraticEffectFloat(self,coeff_name1,coeff_name2,iBin):
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
        return effect


    def getQuadraticEffect(self,coeff_name1,coeff_name2,iBin,incl_syst=True):
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
            mcstatname=None
        else:
            effect=x[0]
            error=x[1]
            if self.stat_unc_threshold!=None:
                th=self.stat_unc_threshold
            elif self.theory.stat_unc_threshold!=None:
                th=self.theory.stat_unc_threshold
            mcstatname=self.theory.makeMCstatName(self.name,[coeff_name1,coeff_name2],iBin)
            if th>=0 and effect!=0 and abs(error/effect)>th and not mcstatname in self.theory.getThetaNames():
                u=Uncertainty.Uncertainty(nbins=self.nbins,name=mcstatname,affected=coeff_name1+','+coeff_name2,interpol='lin',
                                          releffect=[0.]*iBin+[abs(error/effect)]+[0.]*(self.nbins-iBin-1))
                self.addUncertainty(u)
            if self.getUnc(mcstatname):
                ustr=self.getUnc(mcstatname).getStr(iBin)
                effect=str(effect)+'*'+ustr
        if incl_syst:
            for u in self.uncertainties:
                if u.name==mcstatname:
                    continue
                if c1+','+c2 in u.affected or c2+','+c1 in u.affected:
                    ustr=u.getStr(iBin)
                    if ustr!=None:
                        effect=str(effect)+'*'+ustr
        return effect
        

    def getPredictionStrings(self):
        eftprediction = []

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
        for background in self.backgrounds.values():
            predictions=[p+'+'+b for p,b in zip(predictions,background.getPredictionStrings())]                
        return predictions                   

    def getPredictionFormulas(self):
        if self.binlabels is None:
            self.binlabels=range(self.nbins)
        self.predictions = []
        formulas=self.getPredictionStrings()
        for binlabel,formula in zip(self.binlabels,formulas):
            params=RooArgList()
            for c in self.theory.coeffs:
                if re.match('.*\\b'+c.GetName()+'\\b.*',formula):
                    params.add(c)
            for t in self.theory.thetas:
                if re.match('.*\\b'+t.GetName()+'\\b.*',formula):
                    params.add(t)
            bin_prediction=RooFormulaVar("pred_{}_{}".format(self.name,binlabel),"pred_{}_{}".format(self.name,binlabel),formula,params)
            self.predictions.append(bin_prediction)
        return self.predictions

    def plotPrefit(self, outfile,statonly=False):
        histoname=self.name+'_prefit'
        if statonly:
            histoname+='_statonly'
        if self.binning!=None:
            h=TH1F(histoname,histoname,self.nbins,array('f',self.binning))
        else:
            h=TH1F(histoname,histoname,self.nbins,0,1)
        for i,f in enumerate(self.getPredictionFormulas()):
            h.SetBinContent(i+1,f.getVal()/self.rescaled)
            uperror2=0.
            downerror2=0.            
            for t in self.theory.thetas:
                if statonly and not Utils.isMCstat(self.theory,t.GetName()):
                    continue
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

    def plotBkgPostfit(self, outfile, fitresult):
        hs=[]
        for b in self.backgrounds.values():
            hn=self.name+'_'+b.name+'_postfit'
            if self.binning!=None:
                h=TH1F(hn,hn,self.nbins,array('f',self.binning))
            else:
                h=TH1F(hn,hn,self.nbins,0,1)
            for i,f in enumerate(b.getPredictionStrings()):
                form=RooFormulaVar('form_'+self.name+'_'+b.name,f,self.theory.getThetaArgList())
                h.SetBinContent(i+1,form.getVal()/self.rescaled)
                h.SetBinError(i+1,form.getPropagatedError(fitresult)/self.rescaled)
            hs.append(h)           
            outfile.cd()
            h.Write()
        return hs


    def setSM(self,sm,rescale=1.):
        if self.nbins==None:
            self.nbins=len(sm)
        if self.blind==None:
            self.blind=self.nbins*[1]
        assert len(sm)==self.nbins
        assert len(self.blind)==self.nbins
        self.sm=[x*rescale if b==1 else 0 for x,b in zip(sm,self.blind)]
        return True

    def setSMfromCfg(self,data):
        if Utils.isFile(data):
            spl=data.split(':')
            assert len(spl)==2
            return self.setSMfromFile(spl[0],spl[1])
        else:
            return self.setSM([float(x) for x in data.split()])
        
    def setSMfromFile(self,filename,histoname,rescale=1.):
        return self.setSM(Utils.readFile(filename,histoname,per_bin_width=self.per_bin_width)[0],rescale)

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
            if data.lower()=='none' or data.lower()=='false' or data.lower()=='':
                self.normhisto=None
                return True
            return self.setNormHisto([float(x) for x in data.split()])
        
    def setNormHistoFromFile(self,filename,histoname,rescale=1.):
        return self.setNormHisto(Utils.readFile(filename,histoname,per_bin_width=self.per_bin_width)[0],rescale)
        
    def setEffect(self,power,coeff_names,effect,uncertainty=None,rescale=1.,renorm=False,add=False,isrel=False,corrections={}):
        if self.blind is None:
            self.blind==self.nbins*[1]
        effect=[e if b==1 else 0 for e,b in zip(effect,self.blind)]
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
            if self.normhisto==None:
                raise Exception('cannot use option "renorm", no histogram given to renorm prediction')
            else:
                actual_effect=[(s*(e[0]/n),s*(e[1]/n)) for s,e,n in zip(self.sm,actual_effect,self.normhisto)]
        if len(coeff_names)==1 and coeff_names[0] in corrections:
            corr=corrections[coeff_names[0]]
            actual_effect=[(e[0]+s*corr,e[1]) for s,e in zip(self.sm,actual_effect)]

        if power==1:
            if add and coeff_names[0] in self.linear_effects:
                self.linear_effects[coeff_names[0]]=[(a[0]+b[0],(a[1]**2+b[1]**2)**0.5) for a,b in zip(self.linear_effects[coeff_names[0]],actual_effect)]
            else:
                self.linear_effects[coeff_names[0]]=actual_effect
        elif power==2:
            if coeff_names[0]<coeff_names[1]:
                c0,c1=coeff_names[0],coeff_names[1]
            else:
                c1,c0=coeff_names[0],coeff_names[1]
            if not c0 in self.quadratic_effects:
                self.quadratic_effects[c0]=dict()
            if add and c1 in self.quadratic_effects[c0]:                
                self.quadratic_effects[c0][c1]=[(a[0]+b[0],(a[1]**2+b[1]**2)**0.5) for a,b in zip(self.quadratic_effects[c0][c1],actual_effect)]
            else:
                self.quadratic_effects[c0][c1]=actual_effect
        else:
            raise Exception('Only linear or quadratic effects possible')
        return True
    
    def setEffectFromFile(self,power,coeff_names,filename,histoname,rescale=1.,renorm=False,add=False,corrections={}):
        if not os.path.exists(filename.strip()):
            return False
        else:
            out=Utils.readFile(filename,histoname,per_bin_width=self.per_bin_width)
            if out is None:
                return False
            actual_effect,actual_uncertainty=out
            self.setEffect(power,coeff_names,actual_effect,actual_uncertainty,rescale,renorm,add=add,corrections=corrections)
            return True

    def setEffectFromCfg(self,power,coeff_names,data,add=False,corrections={}):
        actual_coeff_names=coeff_names.split()        
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
        
        if Utils.isFile(data):
            filename,histoname=data.split(':')
            if isrel:
                raise Exception('Effect from file shouldnt be relative')
            return self.setEffectFromFile(power,actual_coeff_names,filename,histoname,rescale=1.,renorm=renorm,add=add,corrections=corrections)
                
        else:
            effect=[float(d) for d in data.split()]
            return self.setEffect(power,actual_coeff_names,effect,uncertainty=None,rescale=1.,renorm=renorm,add=add,isrel=isrel,corrections=corrections)

    def setScalingFromCfg(self,data):
        assert '{' in data
        pars=list(set(d.lower() for d in data[data.find('{')+1:data.find('}')].split()))
        data=data[data.find('}')+1:]
        if Utils.isFile(data):
            scaling=Utils.readFile(data,per_bin_width=self.per_bin_width)[0]
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
        self.setScaling(scaling,pars)

    def setScaling(self,scaling,pars):
        if len(scaling)==1:
            scaling=scaling*self.nbins
        assert len(scaling)==self.nbins
        self.extra_scaling=scaling
        self.scaling_pars+=pars

    def addBackgroundFromCfg(self,name,data):
        if Utils.isFile(data):
            filename,histoname=data.split(':')
            values,errors=Utils.readFile(data,per_bin_width=self.per_bin_width)
        else:
            values=[float(d) for d in data.split()]
        self.addBackground(name,values)
        
    def addBackground(self,name,values):
        assert len(values)==self.nbins
        assert not name in self.backgrounds
        assert not name in self.theory.getCoeffNames()
        self.backgrounds[name]=Background(name,values)

    def addBkgUncertaintyFromCfg(self,name,data):
        bkgname=name[name.find('{')+1:name.find('}')].strip()
        assert len(bkgname)>0 and not ' ' in bkgname
        self.addUncertainty(Uncertainty.createUncertaintyFromCfg(self.backgrounds[bkgname],name,data,reco_level=True,is_bkg_unc=True,per_bin_width=self.per_bin_width,isBkgUnc=True))
        
    def addUncertaintyFromCfg(self,name,data,reco_level,is_bkg_unc=False):
        self.addUncertainty(Uncertainty.createUncertaintyFromCfg(self,name,data,reco_level=reco_level,is_bkg_unc=is_bkg_unc,per_bin_width=self.per_bin_width))
   
    def addUncertainty(self,uncertainty):
        if uncertainty==None:
            print 'WARNING: empty uncertainty (likely not found in histo file)'
            return
        if uncertainty.is_bkg_unc:
            found_bkg=False
            for a in uncertainty.affected:
                if a in self.backgrounds:
                    self.backgrounds[a].addUncertainty(uncertainty)
                    found_bkg=True
            if found_bkg:
                self.theory.addTheta(uncertainty.name)
            return found_bkg

        for u in self.uncertainties:
            if u.name==uncertainty.name:
                for a in uncertainty.affected:
                    if a in u.affected:
                        print 'WARNING: uncertainty with name',uncertainty.name,'affecting',a,'already exists, skipping'
                        return False
                    elif ',' in  a and a.split(',')[1]+','+a.split(',')[0] in u.affected:
                        print 'WARNING: uncertainty with name',uncertainty.name,'affecting',a,'already exists, skipping'
                        return False
        self.uncertainties.append(uncertainty)
        self.theory.addTheta(uncertainty.name)
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
        for b in self.backgrounds:
            self.backgrounds[b]=[x*s for x in self.backgrounds[b]]
        self.rescaled*=s

    def getNominal(self):
        return self.sm

    def getUnc(self,name):
        for u in self.uncertainties:
            if name==u.name:
                return u
        return None

    def plotScaling2D(self,outfile):
        par_names=list(self.quadratic_effects)
        if len(par_names)<1:
            return
        for iBin in range(self.nbins):
            h=TH2F(self.name+'_QuadEffects_bin'+str(iBin),'Quadratic effects',len(par_names),0.,len(par_names),len(par_names),0,len(par_names))
            hnorm=TH2F(self.name+'_QuadEffectsNorm_bin'+str(iBin),'Quadratic effects normalized',len(par_names),0.,len(par_names),len(par_names),0,len(par_names))
            for i,ipar in enumerate(par_names):
                for j,jpar in enumerate(par_names):
                    h.SetBinContent(i+1,j+1,self.getQuadraticEffectFloat(ipar,jpar,iBin))
                    if self.getQuadraticEffectFloat(jpar,jpar,iBin)>0 and self.getQuadraticEffectFloat(ipar,ipar,iBin)>0:
                        hnorm.SetBinContent(i+1,j+1,self.getQuadraticEffectFloat(ipar,jpar,iBin)/sqrt(self.getQuadraticEffectFloat(ipar,ipar,iBin)*self.getQuadraticEffectFloat(jpar,jpar,iBin)))
            for i,ipar in enumerate(par_names):
                h.GetXaxis().SetBinLabel(i+1,ipar)
                h.GetYaxis().SetBinLabel(i+1,ipar)
                hnorm.GetXaxis().SetBinLabel(i+1,ipar)
                hnorm.GetYaxis().SetBinLabel(i+1,ipar)
            outfile.cd()
            h.Write()
            hnorm.Write()

    def getQuadParaMatrix(self,ibin):
        import numpy as np
        m=[]
        for c1 in self.coeff_names:
            m.append([])
            for c2 in self.coeff_names:
                es=self.getEffectFloat([c1,c2])
                e=es[ibin] if es is not None else 0.
                if c1==c2:
                    m[-1].append(e)
                else:
                    m[-1].append(e/2)
        matrix=np.matrix(m)
        return matrix
    # M*old = new
    def getReParaMatrix(self,cs_new,cs_old,mapping):
        import numpy as np
        m=[]
        for cnew in cs_new:
            m.append([])
            for cold in cs_old:
                m[-1].append(mapping[cnew][cold])
        matrix=np.matrix(m)
        return matrix
                    
    def rePara(self,cs_new,mapping):
        print 'Warning, during reparametrization, statistical uncertainties are currently not propagated for quadratic effects'
       
        lin_old=self.linear_effects
        effects_matrices_old=[]
        for ibin in range(self.nbins):
            effects_matrices_old.append(self.getQuadParaMatrix(ibin))
        
        cs_old=self.coeff_names
        self.coeff_names=cs_new
            
        repara_matrix=self.getReParaMatrix(cs_new,cs_old,mapping)
        RI=repara_matrix.getI()

        self.linear_effects={}
        for inew,cnew in enumerate(cs_new):
            sum_effect=[0.]*self.nbins
            sum_err=[0.]*self.nbins
            for iold,cold in enumerate(cs_old):
                if not cold in lin_old:
                    continue
                scale=RI[iold,inew]
                sum_effect=[x+y[0]*scale for x,y in zip(sum_effect,lin_old[cold])]
                sum_err=[(x**2+(y[1]*scale)**2)**0.5 for x,y in zip(sum_err,lin_old[cold])]
            self.setEffect(1,cnew,sum_effect,sum_err)
            #self.linear_effects[c]=[(x,y) for x,y in zip(sum_effect,sum_err)]

        self.quadratic_effects={}
        effects_matrices_new=[RI.getT()*M*RI for M in effects_matrices_old]
        for i1,cnew1 in enumerate(cs_new):
            for i2,cnew2 in enumerate(cs_new):
                if i2<i1:
                    continue
                new_effect=[x[i1,i2] for x in effects_matrices_new]
                if not i1==i2:
                    new_effect=[2*x for x in new_effect]
                self.setEffect(2,[cnew1,cnew2],new_effect)
                    
                    
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
        self.binning=None
        print 'WARNING Rebinning: binlabels, blinding, bin boundaries, and backgrounds currently not properly updated'
        self.blind=[1 if b!=0 else 0 for b in Utils.mergeBins(self.blind,combine)]
        assert len(self.backgrounds)==0

        new_sm=Utils.mergeBins(self.sm,combine)
        for u in self.uncertainties:
            newEffectUp=Utils.mergeBins([x[1] for x in u.getEffect(self.getNominal())],combine)
            newEffectDown=Utils.mergeBins([x[0] for x in u.getEffect(self.getNominal())],combine)
            newRelEffectUp=[e/n for e,n in zip(newEffectUp,new_sm)]
            newRelEffectDown=[e/n for e,n in zip(newEffectDown,new_sm)]
            u.releffect=zip(newRelEffectDown,newRelEffectUp)
            u.nbins=len(combine)        

        for c in self.linear_effects:
            vals=Utils.mergeBins([x[0] for x in self.linear_effects[c]],combine)
            errs=Utils.mergeBins([x[1] for x in self.linear_effects[c]],combine)
            self.linear_effects[c]=zip(vals,errs)
        for c1 in self.quadratic_effects:
            for c2 in self.quadratic_effects[c1]:
                vals=Utils.mergeBins([x[0] for x in self.quadratic_effects[c1][c2]],combine)
                errs=Utils.mergeBins([x[1] for x in self.quadratic_effects[c1][c2]],combine)
                self.quadratic_effects[c1][c2]=zip(vals,errs)
        if not self.extra_scaling is None and isinstance(self.extra_scaling,list):
            self.extra_scaling=Utils.mergeBins(self.extra_scaling,combine,mode=0)

        if not self.normhisto is None:
            self.normhisto=Utils.mergeBins(self.sm,combine)
        self.sm=new_sm
        self.nbins=len(combine)


