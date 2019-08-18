from ROOT import RooRealVar, TMatrixDSym, TH1F
from array import array

class Measurement:
    def __init__(self, name, nbins=None, measured=None, covariance=None, per_bin_width=False,binning=None):
        self.name=name
        self.nbins=nbins
        if measured!=None:
            self.setMeasured(measured)
        if covariance!=None:
            self.setCov(covariance)
        self.per_bin_width=per_bin_width
        self.binning=binning
        self.rescaled=1.
       

    def setMeasured(self,measured):
        if self.nbins==None:
            self.nbins=len(measured)
        assert len(measured)==self.nbins
        self.central_values = []
        for i,x_meas in enumerate(measured):
            x=RooRealVar(self.name+"_meas_"+str(i),self.name+"_meas_"+str(i),0,1e9)
            x.setVal(x_meas)
            x.setConstant()
            self.central_values.append(x)

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

    def setMeasuredFromCfg(self,data):
        if data.split(':')[0].endswith('.yoda'):
            spl=data.split(':')
            assert len(spl)==2
            self.setMeasuredFromYoda(spl[0],spl[1])
        elif data.split(':')[0].endswith('.root'):
            spl=data.split(':')
            assert len(spl)==2
            self.setMeasuredFromRoot(spl[0],spl[1])
        else:
            self.setMeasured([float(x) for x in data.split()])

    def setCovFromCfg(self,data):
        if data.split(':')[0].endswith('.yoda'):
            spl=data.split(':')
            assert len(spl)==2
            self.setCov(self.getMatrixFromYoda(spl[0],spl[1]))
        elif data.split(':')[0].endswith('.root'):
            spl=data.split(':')
            assert len(spl)==2
            self.setCov(self.getMatrixFromRoot(spl[0],spl[1]))
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
        else:
            corr=[]
            for y in data.split(','):
                corr+=[[float(x) for x in y.split()]]
            if not ',' in data:
                corr=[float(data)]

        assert len(corr)==self.nbins
        assert len(errs)==self.nbins
        cov=[]
        for i in range(self.nbins):
            cov+=[[0.]*self.nbins]
            for j in range(self.nbins):
                cov[i][j]=errs[i]*errs[j]*corr[i][j]
        self.setCov(cov)

    def setMeasuredFromYoda(self,filename,name):
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
        self.setMeasured(measured)

    def setMeasuredFromRoot(self,filename,name):
        import ROOT
        f=ROOT.TFile(filename)
        histo=f.Get(name)
        if self.per_bin_width:
            contents=[histo.GetBinContent(iBin)/histo.GetBinWidth(iBin) for iBin in range(1,histo.GetNbinsX()+1)]
        else:
            contents=[histo.GetBinContent(iBin) for iBin in range(1,histo.GetNbinsX()+1)]
        self.setMeasured(contents)
        f.Close()

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

    def scale(self,s):
        self.covMatrix*=s*s       
        for x in self.central_values:
            x.setVal(x.getVal()*s)
        self.rescaled*=s

    def plotPrefit(self, outfile):
        if self.binning!=None:
            h=TH1F(self.name+'_data',self.name+'_data',self.nbins,array('f',self.binning))
        else:
            h=TH1F(self.name+'_data',self.name+'_data',self.nbins,0,1)
        for i in range(self.nbins):
            h.SetBinContent(i+1,self.central_values[i].getVal()/self.rescaled)
            h.SetBinError(i+1,(self.covMatrix[i][i])**0.5/self.rescaled)
        outfile.cd()
        h.Write()
        return h
