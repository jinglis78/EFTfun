#!/usr/bin/env python
import sys,ROOT,re
from array import array
sys.path.insert(0, '../python')
sys.path.insert(0, 'python')
import Utils
goodcolors=[
    ROOT.TColor(9000,          0./255, 119./255, 187./255),
    ROOT.TColor(9001,         204./255,  51./255,  17./255),
    ROOT.TColor(9002,        34./255, 136./255,  51./255),
    ROOT.TColor(9003,      238./255, 119./255,  51./255),
    ROOT.TColor(9004,        187./255, 187./255, 187./255),
    ROOT.TColor(9005,      204./255, 187./255,  68./255),
    ROOT.TColor(9006,      170./255,  51./255, 119./255),
    ROOT.TColor(9007,         238./255, 102./255, 119./255),
    ROOT.TColor(9008,         51./255, 187./255, 238./255),
    ROOT.TColor(9009,          0./255, 153./255, 136./255),
    ROOT.TColor(9010,        110./255, 110./255, 110./255),
    ROOT.TColor(9011,         245/255., 66/255., 221/255.)
]

def canvasProfiles1D(name,filenames,xtitle,ytitle='#minus 2#Delta log(L)',linewidth=2,linestyle=1,outfolder='.',isNLL=True):
    files=[ROOT.TFile(fn) for fn in filenames]
    graphs=[]
    for i,fin in enumerate(files):
        for k in [x.GetName() for x in fin.GetListOfKeys()]:
            if re.match(pattern,k):
                g=fin.Get(k).Clone('g'+str(i))
                if isinstance(g,ROOT.TGraph):
                    graphs+=[g]
                    break
    f=ROOT.TFile(name+'.root','recreate')
    first=True
    for ig,g in enumerate(graphs):
        linecolor=9000+ig
        cl1=1
        cl2=3.841
        if first:
            c=ROOT.TCanvas(name,name,1024,1024 )
        c.SetRightMargin(0.02)
        c.SetTopMargin(0.08)
        c.SetLeftMargin(0.12)
        c.SetBottomMargin(0.12)
        c.SetTicks(1,1)
        xs=array('f',g.GetX())
        ys=array('f',g.GetY())
        imin=-1
        ymin=999.
        x0=0.
        for i,y in enumerate(ys):
            if y<ymin:
                ymin=y
                imin=i
                x0=xs[i]
        if isNLL:
            if imin>0 and ymin<0.01 and ymin>-0.001:
                g_inverse1=ROOT.TGraph(imin,ys[:imin],xs[:imin])
                g_inverse2=ROOT.TGraph(g.GetN()-imin,ys[imin:],xs[imin:])
                ci1=g_inverse1.Eval(cl1),g_inverse2.Eval(cl1)
                ci2=g_inverse1.Eval(cl2),g_inverse2.Eval(cl2)
        g.SetLineColor(linecolor)
        g.SetLineStyle(linestyle)
        g.SetLineWidth(linewidth)
        if first:
            g.Draw('AL')
        else:
            g.Draw('L')
        if isNLL :
            g.GetHistogram().GetYaxis().SetRangeUser(0.,g.GetHistogram().GetYaxis().GetXmax())
        g.GetHistogram().GetXaxis().SetTitleSize(0.05)
        g.GetHistogram().GetYaxis().SetTitleSize(0.05)
        g.GetHistogram().GetXaxis().SetLabelSize(0.04)
        g.GetHistogram().GetYaxis().SetLabelSize(0.04)
        g.GetHistogram().GetXaxis().SetTitleOffset(0.9)
        g.GetHistogram().GetYaxis().SetTitleOffset(1)
        g.GetHistogram().SetLineColor(ROOT.kBlack)
        g.GetHistogram().GetXaxis().SetTitle(xtitle)
        g.GetHistogram().GetYaxis().SetTitle(ytitle)
        if isNLL:
            xmax=g.GetHistogram().GetXaxis().GetXmax()
            xmin=g.GetHistogram().GetXaxis().GetXmin()
            cl1line=ROOT.TLine(xmin,cl1,xmax,cl1)
            cl2line=ROOT.TLine(xmin,cl2,xmax,cl2)
            cl1line.Draw()
            cl2line.Draw()

        if isNLL and imin>0 and ymin<0.01 and ymin>-0.001:
            ci2line1=ROOT.TLine(ci2[0],0,ci2[0],cl2)
            ci2line1.SetLineStyle(2)
            ci2line2=ROOT.TLine(ci2[1],0,ci2[1],cl2)
            ci2line2.SetLineStyle(2)
            ci1line1=ROOT.TLine(ci1[0],0,ci1[0],cl1)
            ci1line1.SetLineStyle(2)
            ci1line2=ROOT.TLine(ci1[1],0,ci1[1],cl1)
            ci1line2.SetLineStyle(2)
            ci1text=ROOT.TLatex(0.5,0.85,'Best fit: {}_{{ #minus{} }}^{{ #plus{} }}'.format(Utils.round2(x0,3),Utils.round2(abs(ci1[0]-x0),3),Utils.round2(abs(ci1[1]-x0),3)))
            ci1text.SetNDC()
            ci1text.SetTextSize(0.04)
            ci1text.SetTextFont(42)
            ci1text.SetTextAlign(21)
            ci2text=ROOT.TLatex(0.5,0.79,'95% CL: [{},{}]'.format(Utils.round2(ci2[0],2),Utils.round2(ci2[1],3)))
            ci2text.SetNDC()
            ci2text.SetTextSize(0.04)
            ci2text.SetTextFont(42)
            ci2text.SetTextAlign(21)
            ci1line1.Draw()
            ci1line2.Draw()
            ci2line1.Draw()
            ci2line2.Draw()
            ci1text.Draw()
            ci2text.Draw()
        first=False
    #g.Draw('AL')
    legend=ROOT.TLegend(0.6,0.65,0.9,0.9)
    legend.SetBorderSize(0)
    legend.SetLineStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.04)
    legend.SetFillStyle(0)
    for g in graphs:
        legend.AddEntry(g,g.GetName(),'L')
    legend.Draw()
    c.Draw()
    raw_input()
    f.Write()
    c.SaveAs(name+'.pdf')
    f.Close()

    return c

name=sys.argv[1]
pattern=sys.argv[2]
xtitle=sys.argv[3]
filenames=sys.argv[4:]
canvasProfiles1D(name,filenames,xtitle=xtitle,ytitle='#minus 2#Delta log(L)',linewidth=2,linestyle=1,outfolder='.',isNLL=False)

