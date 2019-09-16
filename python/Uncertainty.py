from ROOT import TMatrixDSym
import Utils

def createUncertaintyFromCfg(reference,name,data):
    if '{' in name:
        affected=set(o.strip() for o in name[name.find('{')+1:name.find('}')].split())
        if 'all' in affected:
            affected='all'
        thetaname=name[name.find('}')+1:].strip()
    else:
        thetaname=name
        affected='all'

    options,data=Utils.getOptionsAndData(data)
    interpol='lin'
    for o in options:
        if o in ['lin','exp']:
            interpol=o
            options.remove(o)
            break
    isrel=None
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
    uncorr=None
    if 'corr' in options:
        uncorr=False
        options.remove('corr')
    elif 'uncorr' in options:
        uncorr=True
        options.remove('uncorr')

    if len(options)>0:
        raise Exception('unknow/incompatible options in uncertainty '+data+' : '+' '.join(list(options)))
    releffect=None
    if Utils.isFile(data):
        spl=data.split(',')
        assert len(spl)==2
        down_filename,down_histoname=spl[0].split(':')
        up_filename,up_histoname=spl[1].split(':')
        if isrel!=None:
            print 'WARNING: abs/rel wont have an effect on uncertainties from histograms'
        return createUncertaintyFromFile(reference,thetaname,down_filename,up_filename,down_histoname,up_histoname,affected=affected,interpol=interpol,renorm=renorm,uncorr=uncorr)
    if Utils.isListOfNumbers(data,reference.nbins):
        effect=[]
        for e in data.split():
            if ',' in e:
                effect+=[ ( float(e.split(',')[0]) , float(e.split(',')[1]) ) ]
            else:
                effect+=[ ( -float(e), float(e)) ]
        if len(effect)==1:
            effect=effect*reference.nbins
        if isrel==True:
            releffect=effect
        elif isrel==False:
            releffect=[(e[0]/s,e[1]/s) for e,s in zip(effect,reference.getNominal())]
        else:
            raise Exception('Uncertainty needs to have option rel or abs')
        if renorm==True:
            if reference.normhisto==None:
                print 'WARNING: cannot use option "renorm", no histogram given to renorm prediction'
            else:
                releffect=[(e[0]*s/n,e[1]*s/n) for s,e,n in zip(reference.getNominal(),effect,reference.normhisto)]
        return createUncertainty(reference,thetaname,releffect=releffect,affected=affected,interpol=interpol,uncorr=uncorr)
    if isrel!=None:
        print 'WARNING: abs/rel wont have an effect on uncertainties from hepdata'
    if renorm!=None:
        print 'WARNING: renorm/keepnorm wont have an effect on uncertainties from hepdata'
    return createUncertaintyFromYaml(reference,thetaname,data,affected=affected,interpol=interpol,uncorr=uncorr)


def createUncertaintyFromYaml(reference,name,data,affected,interpol,uncorr):
    if ':' in data:
        filename,err_label=data.split(':')
    else:
        filename,err_label=reference.hepdata,data
    releffect=getErrsFromYaml(reference,filename,err_label,rel=True)
    return Uncertainty(reference.nbins,affected=affected,name=name,releffect=releffect,interpol=interpol,uncorr=uncorr)
    
def createUncertaintyFromFile(reference,name,down_filename,up_filename,down_histoname,up_histoname,affected='all',interpol='exp',renorm=False,uncorr=False):
    #TODO: one-sided uncertainties
    if down_filename==None:
        down=[0.]*reference.nbins
    else:
        down=Utils.readFile(down_filename,down_histoname)[0]
    if up_filename==None:
        up=[0.]*reference.nbins
    else:
        up=Utils.readFile(up_filename,up_histoname)[0]
    assert len(up)==reference.nbins
    assert len(down)==reference.nbins
    if renorm:
        up=[u/n-1 for u,n in zip(up,reference.normhisto)]
        down=[d/n-1 for d,n in zip(down,reference.normhisto)]
    else:
        up=[u/n-1 for u,n in zip(up,reference.getNominal())]
        down=[d/n-1 for d,n in zip(down,reference.getNominal())]
    releffect=zip(up,down)
    return createUncertainty(name,releffect,affected,interpol,uncorr)

def createUncertainty(reference,name,releffect,affected='all',interpol='exp',uncorr=False):
    if not isinstance(releffect,list):
        releffect=[releffect]*reference.nbins
    assert len(releffect)==reference.nbins
    return Uncertainty(reference.nbins,affected,name,releffect,interpol,uncorr)

def getErrsFromYaml(reference,filename,err_label,rel=True):
    import yaml
    f=open(filename)
    vals=yaml.load(f)['dependent_variables'][0]['values']
    assert len(vals)==reference.nbins
    errs=[]
    for y,meas in zip(vals,reference.getNominal()):
        for error in y['errors']:
            if error['label']==err_label:                    
                if 'symerror' in error:
                    if isinstance(error['symerror'],str) and '%' in error['symerror']:
                        err=float(error['symerror'].replace('%',''))/100.
                        if not rel:
                            err*=meas
                        errs.append((-err,err))
                    else:
                        err=float(error['symerror'])
                        if rel:
                            err/=meas
                        errs.append((-err,err))
                    break
                elif 'asymerror' in error:
                    if isinstance(error['asymerror']['minus'],str) and '%' in error['asymerror']['minus']:
                        asymerr=(float(error['asymerror']['minus'].replace('%',''))/100.,float(error['asymerror']['plus'].replace('%',''))/100.)
                        if not rel:
                            errs.append((asymerr[0]*meas,asymerr[1]*meas))
                        else:
                            errs.append(asymerr)
                    else:
                        asymerr=(float(error['asymerror']['minus']),float(error['asymerror']['plus']))
                        if rel:
                            errs.append((asymerr[0]/meas,asymerr[1]/meas))
                        else:
                            errs.append(asymerr)
                    break        
    f.close()
    return errs

class Uncertainty:
    def __init__(self,nbins,affected,name,releffect,interpol,uncorr=False):
        self.nbins=nbins
        self.affected=affected
        self.name=name
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
        self.uncorr=uncorr

        
    def getStr(self,iBin):
        if self.uncorr==True:
            raise Exception("Cannot profile uncorrelated uncertainty "+self.name)

        if self.releffect[iBin]==0.:
            return None
        if isinstance(self.releffect[iBin],tuple) and self.releffect[iBin][0]==0 and self.releffect[iBin][1]==0:
            return None
        #TODO: define better interpolation functions
        if self.interpol=='exp':
            if abs(-self.releffect[iBin][0]/self.releffect[iBin][1]-1)<1e-2:
                return 'exp({0}*{1})'.format(self.releffect[iBin][1],self.name)
            else:
                return '(exp(-{0}*{2})*({2}<0)+exp({1}*{2})*({2}>=0))'.format(self.releffect[iBin][0],self.releffect[iBin][1],self.name) 
        if self.interpol=='lin':
            if abs(-self.releffect[iBin][0]/self.releffect[iBin][1]-1)<1e-2:
                return '(1+{0}*{1})'.format(self.releffect[iBin][1],self.name)
            else:
                return '((1-{0}*{2})*({2}<0)+(1+{1}*{2})*({2}>=0))'.format(self.releffect[iBin][0],self.releffect[iBin][1],self.name)

    def getRelEffect(self):
        return self.releffect

    def getRelEffectSym(self):
        return [(x[0]-x[1])/2 for x in self.releffect]

    def getEffect(self,nominal):
        return [(n*x[0],n*x[1]) for x,n in zip(self.releffect,nominal)]

    def getEffectSym(self,nominal):
        return [(x[0]-x[1])/2 for x in self.getEffect(nominal)]

    def getCov(self,nominal):
        if self.interpol=='exp':
            print "WARNING: exponential interpolation wont affect covariance matrix of uncertainty",self.name
        cov = TMatrixDSym(self.nbins)
        if self.uncorr:
            for i in range(self.nbins):
                cov[i,i]=self.getEffectSym(nominal)[i]**2
        else:
            for i in range(self.nbins):
                for j in range(self.nbins):
                    cov[i,j]=self.getEffectSym(nominal)[i]*self.getEffectSym(nominal)[j]

        return cov
