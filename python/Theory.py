from ROOT import RooRealVar,RooFit,RooArgList
import re
import Utils

class Theory:
    def __init__(self,default_range=(-1e3,1e3),theta_range=(-4,4),stat_unc_threshold=0.3):
        self.default_range=default_range
        self.coeffs=[]
        self.thetas=[]
        self.pars={}
        self.initial_vals={}
        self.theta_range=theta_range
        self.stat_unc_threshold=stat_unc_threshold
        self.mapping=None

    def reset(self):
        for par in self.pars:
            if par in self.initial_vals:
                self.pars[par].setVal(self.initial_vals[par])
            else:
                self.pars[par].setVal(0.)

    def addCoeff(self,coeff_name,coeff_range=None):
        if coeff_range==None:
            coeff_range=self.default_range
        self.coeffs.append(RooRealVar(coeff_name,coeff_name,coeff_range[0],coeff_range[1]))
        self.pars[self.coeffs[-1].GetName()]=self.coeffs[-1]

    def addTheta(self,theta_name,theta_range=None):
        if theta_name in [t.GetName() for t in self.thetas]:
            return
        if theta_range==None:
            theta_range=self.theta_range
        self.thetas.append(RooRealVar(theta_name,theta_name,theta_range[0],theta_range[1]))
        self.pars[self.thetas[-1].GetName()]=self.thetas[-1]

    def setConstantFromCfg(self,constant,pois):
        if constant==None:
            pass
        elif constant.lower()=='all':
            self.setConstant()
#            for theta in self.thetas:
#                for c in self.coeffs:
#                    if c.GetName() in pois:
#                        continue
#                    if re.match(self.makeMCstatName('.+',c.GetName(),'\d+'),theta.GetName()):
#                        theta.setVal(0.)
#                        theta.setConstant()
        else:
            for c in constant.split(','):
                if '=' in c:
                    c,v=[x.strip() for x in c.split('=')]
                else:
                    c,v=c,0.
                self.setConstant(c,float(v))
                if float(v)!=0.:
                    self.initial_vals[c]=float(v)
#                for theta in self.thetas:
#                    if re.match(self.makeMCstatName('.+',c,'\d+'),theta.GetName()):
#                        theta.setVal(0.)
#                        theta.setConstant()
        for p in pois:
            self.setFloating(p)
                
        
    def setConstant(self,coeff=None,value=0.):
        if coeff==None:
            for c in self.coeffs:
                self.setConstant(c.GetName(),0.)
            return
        elif isinstance(coeff,dict):
            for key in coeff:
                self.setConstant(key,coeff[key])
            return
        else:
            for c in self.coeffs+self.thetas:
                if c.GetName() == coeff:
                    print 'setting',coeff,'constant to',value 
                    c.setVal(value)
                    c.setConstant()
                    return
        raise Exception('Did not find parameter: '+coeff)
        
    def setFloating(self,coeff=None):
        if coeff==None:
            for c in self.coeffs:
                self.setFloating(c.GetName())
        elif isinstance(coeff,list):
            for c in coeff:
                self.setFloating(c)
        else:
            for c in self.coeffs+self.thetas:
                if c.GetName() == coeff:
                    print 'setting',coeff,'floating'
                    c.setConstant(False)
                    for theta in self.thetas:
                        if re.match(self.makeMCstatName('.+',c.GetName(),'\d+'),theta.GetName()):
                            theta.setConstant(False)
                        elif re.match(self.makeMCstatName('.+',[c.GetName(),'.+'],'\d+'),theta.GetName()):
                            theta.setConstant(False)

        #TODO: drop constraint
    def useUncertainties(self,use=False):
        for theta in self.thetas:
            if use:
                theta.setConstant(False)
            else:
                theta.setVal(0.)
                theta.setConstant(True)
    def deactivateUncertainties(self,pattern=None):
        if pattern==None or len(pattern)==0:
            return
        for theta in self.thetas:
            if pattern=='all' or re.match(pattern,theta.GetName()):
                print 'deactivating uncertainty',theta.GetName()
                theta.setVal(0.)
                theta.setConstant(True)

    def activateUncertainties(self,pattern=None):
        if pattern==None:
            return
        for theta in self.thetas:
            if re.match(pattern,theta.GetName()):
                theta.setConstant(False)

    def getCoeffs(self):
        return self.coeffs

    def getCoeff(self,name):
        for c in self.coeffs+self.thetas:
            if c.GetName() == name:
                return c
        raise Exception('Did not find parameter '+name)
        
        
    def getConstantCoeffs(self):
        cs=[]
        for c in self.coeffs:
            if c.isConstant():
                cs.append(c)
        return cs

    def getNonConstCoeffs(self):
        cs=[]
        for c in self.coeffs:
            if not c.isConstant():
                cs.append(c)
        return cs

    
    def getConstantThetas(self):
        cs=[]
        for c in self.thetas:
            if c.isConstant():
                cs.append(c)
        return cs

    def getCoeffNames(self):
        return [c.GetName() for c in self.coeffs]

    def getThetaNames(self):
        return [c.GetName() for c in self.thetas]

    def getThetaArgList(self):
        l=RooArgList()
        for t in self.thetas:
            l.add(t)
        return l

    def getPar(self,par_name):
        if not par_name in self.pars:
            raise Exception('Did not find parameter: '+par_name)
        return self.pars[par_name]

    def getVal(self,par_name):
        if not par_name in self.pars:
            raise Exception('Did not find parameter: '+par_name)
        return self.pars[par_name].getVal()

    def setVal(self,par_name,val):
        if not par_name in self.pars:
            raise Exception('Did not find parameter: '+par_name)
        self.pars[par_name].setVal(val)

    def setInitialVal(self,par_name,val):
        self.setVal(par_name,val)
        self.initial_vals[par_name]=val
        print 'setting initival val for',par_name,'to',val
    

    def makeMCstatName(self,name,coeff_names,iBin):
        if isinstance(coeff_names,str):
             return 'mcstat_'+name+'_'+coeff_names+'_lin_bin'+str(iBin)
        elif isinstance(coeff_names,list):
            return 'mcstat_'+name+'_'+'_'.join(coeff_names)+'_quad_bin'+str(iBin)

    def dumpPars(self):
        for c in self.coeffs:
            print c.GetName(),c.getVal(),c.isConstant()
        for c in self.thetas:
            print c.GetName(),c.getVal(),c.isConstant()

    def rePara(self,cs_new,mapping):
        print 'reparametrizing:'
        print 'Note: parameter range will be:',self.default_range
        cs_old=self.getCoeffNames()
        self.coeffs=[]
        for cnew in cs_new:
            if cnew in cs_old:
                raise Exception('Need new parameter name for '+cnew)
            print cnew,'=',mapping[cnew]
            self.addCoeff(cnew)
        self.mapping=mapping


