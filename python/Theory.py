from ROOT import RooRealVar,RooFit,RooArgList,RooRandom
import re
import Utils
import logging

class Theory:
    def __init__(self,default_range=(-1e3,1e3),theta_range=(-5,5),stat_unc_threshold=-1):
        self.default_range=default_range
        self.coeffs=[]
        self.thetas=[]
        self.free_theta_names=[]
        self.global_obs={}
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

    def resetThetas(self):
        for par in self.getThetaNames():
            if par in self.initial_vals:
                self.pars[par].setVal(self.initial_vals[par])
            else:
                self.pars[par].setVal(0.)

                
    def addCoeff(self,coeff_name,coeff_range=None):
        if coeff_range==None:
            coeff_range=self.default_range
        self.coeffs.append(RooRealVar(coeff_name,coeff_name,coeff_range[0],coeff_range[1]))
        self.pars[self.coeffs[-1].GetName()]=self.coeffs[-1]

    def addTheta(self,theta_name,theta_range=None,free=False,multiplicity=0):
        if multiplicity!=0:
            for ibin in range(multiplicity):
                self.addTheta(theta_name+'_bin'+str(ibin),theta_range,free,multiplicity=0)
            return
        if theta_name in [t.GetName() for t in self.thetas]:
            if free and not theta_name in self.free_theta_names:
                raise Exception('Tried to create unconstrained uncertainty for constrained parameter '+theta_name)
            elif not free and theta_name in self.free_theta_names:
                raise Exception('Tried to create constrained uncertainty for unconstrained parameter '+theta_name)                
            return
        if theta_range==None:
            theta_range=self.theta_range
        if free:
            self.free_theta_names.append(theta_name)
        else:
            self.global_obs[theta_name]=RooRealVar(theta_name+'_in',theta_name+'_in',0,theta_range[0],theta_range[1])
            self.global_obs[theta_name].setConstant()
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
                    logging.info('setting {} constant to {}'.format(coeff,value) )
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
                    logging.info('setting '+coeff+' floating')
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
                logging.info('deactivating uncertainty '+theta.GetName())
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
    
    def getCoeffValues(self):
        res={}
        for c in self.coeffs:
            res[c.GetName()]=c.getVal()
        return res
    
    def getThetaValues(self):
        res={}
        for c in self.thetas:
            res[c.GetName()]=c.getVal()
        return res

    def getThetaArgList(self):
        l=RooArgList()
        for t in self.thetas:
            l.add(t)
        return l
    def getParArgList(self):
        l=RooArgList()
        for t in self.thetas:
            l.add(t)
        for c in self.coeffs:
            l.add(c)
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
        logging.info('setting initival val for {} to {}'.format(par_name,val))
    

    def makeMCstatName(self,name,coeff_names,iBin):
        if isinstance(coeff_names,str):
             return 'mcstat_'+name+'_'+coeff_names+'_lin_bin'+str(iBin)
        elif isinstance(coeff_names,list):
            return 'mcstat_'+name+'_'+'_'.join(coeff_names)+'_quad_bin'+str(iBin)

    def dumpPars(self):
        for c in self.coeffs:
            logging.info(' '.join([c.GetName(),str(c.getVal()),str(c.isConstant())]))
        for c in self.thetas:
            logging.info(' '.join([c.GetName(),str(c.getVal()),str(c.isConstant())]))

    def rePara(self,cs_new,mapping):
        logging.info('reparametrizing:')
        logging.info('Note: parameter range will be: '+str(self.default_range))
        cs_old=self.getCoeffNames()
        self.coeffs=[]
        for cnew in cs_new:
            if cnew in cs_old:
                raise Exception('Need new parameter name for '+cnew)
            logging.info(cnew+' = '+mapping[cnew])
            self.addCoeff(cnew)
        self.mapping=mapping

    def isConstant(self):
        for c in self.coeffs + self.thetas:
            if not c.isConstant():
                return False
        return True

    def shuffleGlobalObs(self,central_vals):
        for n,g in list(self.global_obs.items()):
            if not self.pars[n].isConstant():
                g.setVal(central_vals[n]+RooRandom.gaussian())

    def resetGlobalObs(self):
        for g in list(self.global_obs.values()):
            g.setVal(0.)

    def getGlobalObs(self):
        return list(self.global_obs.values())

