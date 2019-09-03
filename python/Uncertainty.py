from ROOT import TMatrixDSym

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
        cov = TMatrixDSym(self.nbins)
        if self.uncorr:
            for i in range(self.nbins):
                cov[i,i]=self.getEffectSym(nominal)[i]**2
        else:
            for i in range(self.nbins):
                for j in range(self.nbins):
                    cov[i,j]=self.getEffectSym(nominal)[i]*self.getEffectSym(nominal)[j]

        return cov
