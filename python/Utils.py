import math,re

#TODO: logger

def isFile(data):
    if data.split(':')[0].lower().endswith('.root') or data.split(':')[0].lower().endswith('.yoda'):
        return True
    return False

def saveToTree(results,outfile):
    from array import array
    from ROOT import TTree
    outfile.cd()
    t=TTree('tree','tree')
    arrays={}
    for r in results:
        arrays[r]=array( 'f', [ 0. ] )
        t.Branch( r, arrays[r], r+'/F' )
    for i in range(len(results[r])):
        for r in results:
            arrays[r][0]=results[r][i]
        t.Fill()
    t.Write()
    


def round2(number,digits=3):
    if number==0:
        return 0
    if math.log10(abs(number))>digits-1:
        return int(round(number))   
    return round(number, digits-int(math.floor(math.log10(abs(number))))-1)

def getOptionsAndData(data):
    if '{' in data:
        options=set(o.lower() for o in data[data.find('{')+1:data.find('}')].split())
    else:
        options=set()
    data=data[data.find('}')+1:].strip()
    return options,data


def readFile(filename,histoname=None,per_bin_width=True):
    if histoname==None:
        assert ":" in filename
        filename=filename.split(':')[0]
        histoname=filename.split(':')[1]
    filename=filename.strip()
    histoname=histoname.strip()
    if filename.lower().endswith('.root'):
        return readRoot(filename,histoname,per_bin_width)
    elif filename.lower().endswith('.yoda'):
        return readYoda(filename,histoname,per_bin_width)
    else:
        pass
    #TODO: error

def readYoda(filename,histoname,per_bin_width):
    import yoda
    histos=yoda.read(filename)
    hist=histos[histoname]
    if per_bin_width:
        contents=[b.height for b in hist.bins]
        uncertainty=[b.heightErr for b in hist.bins]
    else:
        contents= [b.area for b in hist.bins]
        uncertainty=[b.areaErr for b in hist.bins]
    return contents,uncertainty

def readRoot(filename,histoname,per_bin_width):
    import ROOT
    f=ROOT.TFile(filename)
    histo=f.Get(histoname)
    if not isinstance(histo,ROOT.TH1):
        raise Exception("Did not find histogram "+histoname+" in file "+filename)
    if per_bin_width:
        contents=[histo.GetBinContent(iBin)/histo.GetBinWidth(iBin) for iBin in range(1,histo.GetNbinsX()+1)]
        uncertainty=[histo.GetBinError(iBin)/histo.GetBinWidth(iBin) for iBin in range(1,histo.GetNbinsX()+1)]
    else:
        contents=[histo.GetBinContent(iBin) for iBin in range(1,histo.GetNbinsX()+1)]
        uncertainty=[histo.GetBinError(iBin) for iBin in range(1,histo.GetNbinsX()+1)]
    return contents,uncertainty

def isListOfNumbers(data,nbins):
    splt=data.split()
    if not len(splt)==1 and not len(splt)==nbins:
        return False
    for s in splt:
        numbers=s.split(',')
        for n in numbers:
            try:
                float(n)
            except ValueError:
                return False
    return True
    
    
def isMCstat(theory,name):
    if re.match(theory.makeMCstatName('.+','.*','\d+'),name)!=None:
        return True
    if re.match(theory.makeMCstatName('.+',['.*','.*'],'\d+'),name)!=None:
        return True
