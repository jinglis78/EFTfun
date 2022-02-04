from ROOT import TMatrixDSym
import Utils

def createUncertaintyFromCfg(reference,name,data,reco_level,is_bkg_unc=False,per_bin_width=True,isBkgUnc=False):
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
    with_nuis_par=reco_level
    if 'np' in options:
        options.remove('np')
        with_nuis_par=True
    already_in_cov=False
    if 'incov' in options:
        options.remove('incov')
        already_in_cov=True

    flip=False
    if 'flip' in options:
        options.remove('flip')
        flip=True
    sym=False
    if 'sym' in options:
        options.remove('sym')
        sym=True
    if len(options)>0:
        raise Exception('unknow/incompatible options in uncertainty '+data+' : '+' '.join(list(options)))
    releffect=None
    if Utils.isFile(data):
        spl=data.split(',')
        if len(spl)==1:
            up_filename,up_histoname=spl[0].split(':')
            down_filename,down_histoname=None,None
        else:
            if len(spl[0].strip())!=0:
                down_filename,down_histoname=spl[0].split(':')
            else:
                down_filename=None
                down_histoname=None
            if len(spl)>1 and len(spl[1].strip())!=0:
                up_filename,up_histoname=spl[1].split(':')
            else:
                up_filename=None
                up_histoname=None

        return createUncertaintyFromFile(reference,thetaname,down_filename,up_filename,down_histoname,up_histoname,affected=affected,interpol=interpol,renorm=renorm,uncorr=uncorr,with_nuis_par=with_nuis_par,already_in_cov=already_in_cov,is_bkg_unc=is_bkg_unc,flip=flip,sym=sym,per_bin_width=per_bin_width,isBkgUnc=isBkgUnc,isrel=isrel)
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
            releffect=[(e[0]/s,e[1]/s) if s!=0 else (0,0) for e,s in zip(effect,reference.getNominal())]
        else:
            raise Exception('Uncertainty needs to have option rel or abs')
        if renorm==True:
            if reference.normhisto==None:
                print 'WARNING: cannot use option "renorm", no histogram given to renorm prediction'
            else:
                releffect=[(e[0]*s/n,e[1]*s/n) for s,e,n in zip(reference.getNominal(),effect,reference.normhisto)]
        return createUncertainty(reference,thetaname,releffect=releffect,affected=affected,interpol=interpol,uncorr=uncorr,with_nuis_par=with_nuis_par,already_in_cov=already_in_cov,is_bkg_unc=is_bkg_unc,flip=flip,sym=sym)

    if isrel!=None:
        print 'WARNING: abs/rel wont have an effect on uncertainties from hepdata'
    if renorm!=None:
        print 'WARNING: renorm/keepnorm wont have an effect on uncertainties from hepdata'
    return createUncertaintyFromYaml(reference,thetaname,data,affected=affected,interpol=interpol,uncorr=uncorr,with_nuis_par=with_nuis_par,already_in_cov=already_in_cov,is_bkg_unc=is_bkg_unc,flip=flip)


def createUncertaintyFromYaml(reference,name,data,affected,interpol,uncorr,with_nuis_par=False,already_in_cov=False,is_bkg_unc=False,flip=False):
    return createUncertaintyFromYaml(reference,thetaname,data,affected=affected,interpol=interpol,uncorr=uncorr,with_nuis_par=with_nuis_par,already_in_cov=already_in_cov,is_bkg_unc=is_bkg_unc,flip=flip,sym=sym)


def createUncertaintyFromYaml(reference,name,data,affected,interpol,uncorr,with_nuis_par=False,already_in_cov=False,is_bkg_unc=False,flip=False,sym=False):
    if ':' in data:
        filename,err_label=data.split(':')
    else:
        filename,err_label=reference.hepdata,data
    releffect=getErrsFromYaml(reference,filename,err_label,rel=True)
    
    return Uncertainty(reference.nbins,affected=affected,name=name,releffect=releffect,interpol=interpol,uncorr=uncorr,with_nuis_par=with_nuis_par,already_in_cov=already_in_cov,is_bkg_unc=is_bkg_unc,flip=flip,sym=sym)
    
def createUncertaintyFromFile(reference,name,down_filename,up_filename,down_histoname,up_histoname,affected='all',interpol='exp',renorm=False,uncorr=False,with_nuis_par=False,already_in_cov=False,is_bkg_unc=False,flip=False,sym=False,per_bin_width=True,isBkgUnc=False,isrel=False):
    if not Utils.histoExists(down_filename,down_histoname) and not Utils.histoExists(up_filename,up_histoname): 
        print 'WARNING: did not find uncertainty in file:',':'.join([str(down_filename),str(down_histoname)]),'or',':'.join([str(up_filename),str(up_histoname)]),'skipping uncertainty...'
        return
    assert down_filename!=None or up_filename!=None
    if not affected=='all' and not len(affected)==1:
        raise Exception('Uncertainties from histograms can affect only affect the complete signal prediction ("all") or only "sm" or one specific operator')
    if isBkgUnc:
        nominal=reference.getNominal()
    elif affected=='all' or 'sm' in affected:
        if renorm:
            nominal=reference.normhisto
        else:
            nominal=reference.getNominal()
    else:       
        nominal=reference.getEffectFloat(list(affected)[0].split(','))
    if nominal==None:
        raise Exception('Could not find nominal effect for systematic affecting {}'.format(' '.join(list(affected))))
    pbw=per_bin_width
    if isrel:
        pbw=False    
    if up_filename!=None:
        up_variation=Utils.readFile(up_filename,up_histoname,per_bin_width=pbw)[0]
        if isrel:
            up_variation=[n*u+n for n,u in zip(nominal,up_variation)]
    if down_filename!=None:
        down_variation=Utils.readFile(down_filename,down_histoname,per_bin_width=pbw)[0]
        if isrel:
            down_variation=[n*d+n for n,d in zip(nominal,down_variation)]

    if down_filename==None:
        if sym:
            down=[2*n-x for n,x in zip(nominal, up_variation)]
        else:
            down=nominal
    else:
        down=down_variation
    if up_filename==None:
        if sym:
            up=[2*n-x for n,x in zip(nominal, down_variation)]
        else:
            up=nominal
    else:
        up=up_variation
    assert len(up)==reference.nbins
    assert len(down)==reference.nbins
    up=[u/n-1 if n!=0 else 0 for u,n in zip(up,nominal)]
    down=[d/n-1 if n!=0 else 0 for d,n in zip(down,nominal)]
    
    releffect=zip(down,up)

    return createUncertainty(reference,name,releffect,affected,interpol,uncorr,with_nuis_par=with_nuis_par,already_in_cov=already_in_cov,is_bkg_unc=is_bkg_unc,flip=flip,sym=sym)

def createUncertainty(reference,name,releffect,affected='all',interpol='exp',uncorr=False,with_nuis_par=False,already_in_cov=False,is_bkg_unc=False,flip=False,sym=False):
    if not isinstance(releffect,list):
        releffect=[releffect]*reference.nbins
    assert len(releffect)==reference.nbins
    return Uncertainty(reference.nbins,affected,name,releffect,interpol,uncorr,with_nuis_par=with_nuis_par,already_in_cov=already_in_cov,is_bkg_unc=is_bkg_unc,flip=flip,sym=sym)

def getErrsFromYaml(reference,filename,err_label,rel=True):
    import yaml
    f=open(filename)
    vals=yaml.load(f)['dependent_variables'][0]['values']
    assert len(vals)==reference.nbins
    errs=[]
    found=False
    if err_label==None:
        total_errs2=[]
    for y,meas in zip(vals,reference.getNominal()):
        total_err2=0.
        for error in y['errors']:               
            if err_label is None or error['label']==err_label:
                found=True
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
                if not err_label is None:
                    break
                else:
                    total_err2+=(abs(errs[-1][0])+abs(errs[-1][1]))**2/4.
        if err_label is None:
            total_errs2+=[total_err2]

    if not found:
        raise Exception("Did not find error label "+err_label+" in yaml")
    f.close()
    if err_label is None:
        return [x**0.5 for x in total_errs2]
    else:
        return errs

class Uncertainty:
    def __init__(self,nbins,affected,name,releffect,interpol,uncorr=False,with_nuis_par=False,already_in_cov=False,is_bkg_unc=False,flip=False,sym=False):
        self.nbins=nbins
        self.affected=affected
        self.name=name
        assert len(releffect)==self.nbins
        actual_releffect=[]
        for e in releffect:
            if isinstance(e,float):
                if flip:
                    actual_releffect+=[(e,-e)]
                else:
                    actual_releffect+=[(-e,e)]
            else:
                if flip:
                    actual_releffect+=[(e[1],e[0])]
                else:
                    actual_releffect+=[e]
        if sym:
            actual_releffect = [ ( -(e[1]-e[0])/2. ,(e[1]-e[0])/2. ) for e in actual_releffect]            
        assert len(actual_releffect)==self.nbins
        self.releffect=actual_releffect
        assert interpol in ['lin','exp']            
        self.interpol=interpol
        self.uncorr=uncorr
        self.with_nuis_par=with_nuis_par
        self.already_in_cov=already_in_cov
        self.is_bkg_unc=is_bkg_unc

    
        
    def getStr(self,iBin):
        if self.uncorr==True:
            raise Exception("Cannot profile uncorrelated uncertainty "+self.name)

        if self.releffect[iBin]==0.:
            return None
        if isinstance(self.releffect[iBin],tuple) and self.releffect[iBin][0]==0 and self.releffect[iBin][1]==0:
            return None
        #TODO: define better interpolation functions
        if self.interpol=='exp':
            if Utils.isSymmetric(self.releffect[iBin][0],self.releffect[iBin][1],1e-2):
                return 'exp({0}*{1})'.format(self.releffect[iBin][1],self.name)
            else:
                return '(exp(-{0}*{2})*({2}<0)+exp({1}*{2})*({2}>=0))'.format(self.releffect[iBin][0],self.releffect[iBin][1],self.name) 
        if self.interpol=='lin':
            if Utils.isSymmetric(self.releffect[iBin][0],self.releffect[iBin][1],1e-2):
                return '(1+{0}*{1})'.format(self.releffect[iBin][1],self.name)
            else:
                return '((1-{0}*{2})*({2}<0)+(1+{1}*{2})*({2}>=0))'.format(self.releffect[iBin][0],self.releffect[iBin][1],self.name)

    def getRelEffect(self):
        return self.releffect

    def getRelEffectSym(self):
        return [(x[1]-x[0])/2 for x in self.releffect]

    def getEffect(self,nominal):
        return [(n*x[0],n*x[1]) for x,n in zip(self.releffect,nominal)]

    def getEffectSym(self,nominal):
        return [(x[1]-x[0])/2 for x in self.getEffect(nominal)]

    def getCov(self,nominal,incl_np):
        if self.with_nuis_par and not incl_np:
#            print "WARNING: this uncertainty should be handled with a nuisance parameter:",self.name
            return TMatrixDSym(self.nbins)
        if self.interpol=='exp' and not incl_np:
            print "WARNING: exponential interpolation wont affect covariance matrix of uncertainty",self.name
        cov = TMatrixDSym(self.nbins)
        effect=self.getEffectSym(nominal)
        if self.uncorr:
            for i in range(self.nbins):
                cov[i,i]=effect[i]**2
        else:
            for i in range(self.nbins):
                for j in range(self.nbins):
                    cov[i,j]=effect[i]*effect[j]

        return cov
