import ROOT,sys
filename_both=sys.argv[1]
filename_1=sys.argv[2]
scaling_1=float(sys.argv[3])
filename_2=sys.argv[4]
scaling_2=float(sys.argv[5])
outname=sys.argv[6]
dirname=None
if len(sys.argv)>7:
    dirname=sys.argv[7]

    
f0=ROOT.TFile(filename_both)
subfiles=[ROOT.TFile(f) for f in [filename_1,filename_2]]
fout=ROOT.TFile(outname+'.root','RECREATE')
fout_unscaled=ROOT.TFile(outname+'_unscaled.root','RECREATE')
if dirname!=None:
    fout.mkdir(dirname)
    fout_unscaled.mkdir(dirname)

d=f0
if dirname!=None:
    d=f0.GetDirectory(dirname)
for histoname in [k.GetName() for k in d.GetListOfKeys()]:
    if dirname!=None:
        actual_histoname=histoname
        histoname=dirname+'/'+histoname
    h0=f0.Get(histoname).Clone()
    h0.SetLineColor(ROOT.kBlack)
    subhistos=[f.Get(histoname).Clone() for f in subfiles]
    reshisto=h0.Clone(actual_histoname)
    for h in subhistos:
        reshisto.Add(h,-1)
    for i in range(reshisto.GetNbinsX()):
        reshisto.SetBinError(i+1,0.)
    reshisto_unscaled=reshisto.Clone()
    reshisto.Scale(1./scaling_1/scaling_2)
    if dirname!=None:
        fout.cd(dirname)
    else:
        fout.cd()
    reshisto.Write()
    if dirname!=None:
        fout_unscaled.cd(dirname)
    else:
        fout_unscaled.cd()
    reshisto_unscaled.Write()

fout.Close()
fout_unscaled.Close()
f0.Close()
for f in subfiles:
    f.Close()
