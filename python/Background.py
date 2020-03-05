class Background:
    def __init__(self,name,values):
        self.name=name
        self.nominal=values
        self.nbins=len(values)
        self.uncertainties={}

    def getNominal(self):
        return self.nominal

    def addUncertainty(self,uncertainty):
        assert self.name in uncertainty.affected
        assert uncertainty.nbins==self.nbins
        assert uncertainty.with_nuis_par
        if uncertainty.name in self.uncertainties:
            raise Exception('Added same uncertainty twice '+uncertainty.name)
        self.uncertainties[uncertainty.name]=uncertainty
       
    def getPredictionStrings(self):
        ps=[]
        for iBin,x in enumerate(self.getNominal()):
            ps+=[str(x)]
            for u in self.uncertainties.values():
                ustr=u.getStr(iBin)
                if ustr!=None:
                    ps[-1]+="*"+u.getStr(iBin)
        return ps
