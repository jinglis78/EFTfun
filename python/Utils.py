import math,re,os
from glob import glob
import logging

def isFile(data):
    if data.split(':')[0].lower().endswith('.root') or data.split(':')[0].lower().endswith('.yoda'):
        return True
    return False

def histoExists(f,h):
    if f==None or h==None:
        return False
    from ROOT import TFile
    f=TFile(f)
    htmp=f.Get(h)
    res=htmp!=None
    f.Close()
    return res

def saveToTree(results,outfile,name='tree'):
    from array import array
    from ROOT import TTree
    outfile.cd()
    t=TTree(name,name)
    arrays={}
    for r in results:
        arrays[r]=array( 'f', [ 0. ] )
        t.Branch( r, arrays[r], r+'/F' )
    for i in range(len(results[r])):
        for r in results:
            arrays[r][0]=results[r][i]
        t.Fill()
    t.Write()

def saveResultToTree(name,pois,result,toy_result,f):
    res_dict={}
    for ip,poi in enumerate(pois):
        res_dict[poi]=[r[ip] for r in result]
    res_dict['nll']=[r[-1] for r in result]
    if not toy_result is None:
        res_dict['pval']=[r[-2] for r in toy_result]
        res_dict['pval_err']=[r[-1] for r in toy_result]
    saveToTree(res_dict,f,'result_tree_'+name+'_'+'_'.join(pois))
  
def getResultFromFile(name,pois,resultsfiles):
    from ROOT import TFile
    for resultsfile in resultsfiles:
        if not resultsfile.endswith('.root'):
            raise Exception("Need rootfile input")
    logging.info('getting results from '+str(resultsfiles))
    result=[]
    toy_result=[]
    for fn in resultsfiles:
        f=TFile(fn)
        tree_name='result_tree_'+name+'_'+'_'.join(pois)
        t=f.Get(tree_name)
        if not 'pval' in [b.GetName() for b in t.GetListOfBranches()]:
            toy_result=None
        for e in t:
            result.append(tuple([getattr(e,p) for p in pois])+(e.nll,))
            if not toy_result is None:
                toy_result.append(tuple([getattr(e,p) for p in pois])+(e.pval,e.pval_err))
    return sorted(result),(sorted(toy_result) if not toy_result is None else None)

def getFiles(files):
    bad=False
    for fn in files.split():
        if not os.path.exists(fn):
            bad=True
            break
    if not bad:
        return files.split()
    for fn in glob(files):
        if not os.path.exists(fn):
            raise Exception("Can't find files "+str(files))
    return glob(files)

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
        filename,histoname=filename.split(':')[0],filename.split(':')[1]
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
    if isinstance(hist,yoda.Scatter2D):
        points=sorted([(p.x(),p.y()) for p in hist.points()] , key=lambda p: p[0])
        uncertainty=[(p.yErrs()[0]+p.yErrs()[1])/2 for p in hist.points()]
        contents=[p[1] for p in points]
        
    else:
        if per_bin_width:
            contents=[b.height() for b in hist.bins()]
            uncertainty=[b.heightErr() for b in hist.bins()]
        else:
            contents= [b.area() for b in hist.bins()]
            uncertainty=[b.areaErr() for b in hist.bins()]
    return contents,uncertainty

def readRoot(filename,histoname,per_bin_width):
    import ROOT
    f=ROOT.TFile(filename)
    histo=f.Get(histoname)
    if not isinstance(histo,ROOT.TH1):
        logging.warning("Did not find histogram "+histoname+" in file "+filename)
        return None
    if per_bin_width:
        contents=[histo.GetBinContent(iBin)/histo.GetBinWidth(iBin) for iBin in range(1,histo.GetNbinsX()+1)]
        uncertainty=[histo.GetBinError(iBin)/histo.GetBinWidth(iBin) for iBin in range(1,histo.GetNbinsX()+1)]
    else:
        contents=[histo.GetBinContent(iBin) for iBin in range(1,histo.GetNbinsX()+1)]
        uncertainty=[histo.GetBinError(iBin) for iBin in range(1,histo.GetNbinsX()+1)]
    return contents,uncertainty

def isListOfNumbers(data,nbins):
    splt=data.split()
    for s in splt:
        numbers=s.split(',')
        for n in numbers:
            try:
                float(n)
            except ValueError:
                return False
    if not len(splt)==1 and not len(splt)==nbins:
        raise Exception('List of numbers with wrong length ({}): '.format(len(splt))+data)
    return True
    
    
def isMCstat(theory,name):
    if re.match(theory.makeMCstatName('.+','.*','\d+'),name)!=None:
        return True
    if re.match(theory.makeMCstatName('.+',['.*','.*'],'\d+'),name)!=None:
        return True

def ask(question,option1=['yes','y', 'ye', ''],option2=['no','n']):
    print(question)
    o1 = set(option1)
    o2 = set(option2)
    choice = input().lower()
    if choice in o1:
        return True
    elif choice in o2:
        return False
    else:
        print("Please respond with "+','.join(option1)+" or "+','.join(option2))
        return ask(question,option1,option2)

def strip(s):
    return s.strip().replace('"','').replace("'",'')

def split(sin):
    dummy='SPLITTINGPROTECTION'
    s=sin.replace('\ ',dummy)
    spl=[x.replace(dummy,' ') for x in s.split()]
    return spl

def mergeBins(l,combine,mode=1):
    iold=0
    inew=0
    lnew=[0]*len(combine)
    for ic in combine:
       for i in range(ic):
           if mode==0:
               if lnew[inew]==0 or lnew[inew]==l[iold]:
                   lnew[inew]=l[iold]
               else:
                   assert False
           if mode==1:
               lnew[inew]=lnew[inew]+l[iold]
           if mode==2:
               lnew[inew]=(lnew[inew]**2+l[iold]**2)**0.5
           iold+=1
       inew+=1
    return lnew

def isSymmetric(up,down,epsilon=1e-2):
    if up==0 and down==0:
        return True
    if up==0 and down!=0 or up!=0 and down==0:
        return False
    if abs(up/down+1)<epsilon:
        return True

