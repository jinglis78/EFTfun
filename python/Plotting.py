import ROOT,math
from array import array
import Utils

def graphProfile1D(name,result,outfile=None):
    g=ROOT.TGraph()
    result=[r for r in result if not math.isnan(r[1]) and not math.isinf(r[1])]
    for i,r in enumerate(result):            
        g.SetPoint(i,r[0],r[1]*2)
    g.SetName(name)
    if outfile!=None:
        outfile.cd()
        g.Write()
    return g

def graphProfile2D(name,result,outfile=None):
    g=ROOT.TGraph2D()
    result=[r for r in result if not math.isnan(r[2]) and not math.isinf(r[2])]
    for i,r in enumerate(result):
        g.SetPoint(i,r[0],r[1],r[2]*2)
    g.SetName(name)
    if outfile!=None:
        outfile.cd()
        g.Write()

    return g


def canvasProfile1D(name,result,outfile,xtitle,ytitle='#minus 2#Delta log(L)',linecolor=ROOT.kRed+1,linewidth=2,linestyle=1,outfolder='.'):
    cl1=1
    cl2=3.841
    c=ROOT.TCanvas(name,name,1024,1024 )
    c.SetRightMargin(0.02)
    c.SetTopMargin(0.08)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.12)
    c.SetTicks(1,1)
    g=graphProfile1D(name,result)
    xs=array('f',g.GetX())
    ys=array('f',g.GetY())
    imin=-1
    ymin=999.
    xmin=0
    for i,y in enumerate(ys):
        if y<ymin:
            ymin=y
            imin=i
    if imin>-1 and ymin<0.01:
        g_inverse1=ROOT.TGraph(imin,ys[:imin],xs[:imin])
        g_inverse2=ROOT.TGraph(g.GetN()-imin,ys[imin:],xs[imin:])
        ci1=g_inverse1.Eval(cl1),g_inverse2.Eval(cl1)
        ci2=g_inverse1.Eval(cl2),g_inverse2.Eval(cl2)
    g.SetLineColor(linecolor)
    g.SetLineStyle(linestyle)
    g.SetLineWidth(linewidth)
    g.Draw('AL')
    g.GetHistogram().GetXaxis().SetTitleSize(0.05)
    g.GetHistogram().GetYaxis().SetTitleSize(0.05)
    g.GetHistogram().GetXaxis().SetLabelSize(0.04)
    g.GetHistogram().GetYaxis().SetLabelSize(0.04)
    g.GetHistogram().GetXaxis().SetTitleOffset(0.9)
    g.GetHistogram().GetYaxis().SetTitleOffset(1)
    g.GetHistogram().SetLineColor(ROOT.kBlack)
    g.GetHistogram().GetXaxis().SetTitle(xtitle)
    g.GetHistogram().GetYaxis().SetTitle(ytitle)
    xmax=g.GetHistogram().GetXaxis().GetXmax()
    xmin=g.GetHistogram().GetXaxis().GetXmin()
    cl1line=ROOT.TLine(xmin,cl1,xmax,cl1)
    cl2line=ROOT.TLine(xmin,cl2,xmax,cl2)
    cl1line.Draw()
    cl2line.Draw()

    if imin>-1 and ymin<0.01:
        ci2line1=ROOT.TLine(ci2[0],0,ci2[0],cl2)
        ci2line1.SetLineStyle(2)
        ci2line2=ROOT.TLine(ci2[1],0,ci2[1],cl2)
        ci2line2.SetLineStyle(2)
        ci1line1=ROOT.TLine(ci1[0],0,ci1[0],cl1)
        ci1line1.SetLineStyle(2)
        ci1line2=ROOT.TLine(ci1[1],0,ci1[1],cl1)
        ci1line2.SetLineStyle(2)
        ci1text=ROOT.TLatex(0.5,0.85,'68% CL.: [{},{}]'.format(Utils.round2(ci1[0]),Utils.round2(ci1[1])))
        ci1text.SetNDC()
        ci1text.SetTextSize(0.04)
        ci1text.SetTextFont(42)
        ci1text.SetTextAlign(21)
        ci2text=ROOT.TLatex(0.5,0.79,'95% CL: [{},{}]'.format(Utils.round2(ci2[0]),Utils.round2(ci2[1])))
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
    c.Draw()
    printAndSave(c,outfile,outfolder)

    return c

def canvasProfile2D(name,result,outfile,xtitle,ytitle,ztitle='#minus 2#Delta log(L)',linecolor=ROOT.kBlack,outfolder='.'):
    cl1=2.297
    cl2=5.991
    c=ROOT.TCanvas(name,name,1024,1024)
    c.SetRightMargin(0.18)
    c.SetTopMargin(0.12)
    c.SetLeftMargin(0.08)
    c.SetBottomMargin(0.12)
    c.SetTicks(1,1)
    pad1 = ROOT.TPad("pad1","",0,0,1,1)
    pad2 = ROOT.TPad("pad2","",0,0,1,1)
    for p in pad1,pad2:
        p.SetRightMargin(0.16)
        p.SetTopMargin(0.12)
        p.SetLeftMargin(0.1)
        p.SetBottomMargin(0.12)
    pad2.SetFillStyle(0)
    pad2.SetFillColor(0)
    pad2.SetFrameFillStyle(0)
    pad1.Draw()
    pad1.cd()
    g=graphProfile2D(name,result)
    g.Draw('COLZ')
    g.GetHistogram().SetTitle('')
    g.GetHistogram().GetXaxis().SetTitle(xtitle)
    g.GetHistogram().GetYaxis().SetTitle(ytitle)
    g.GetHistogram().GetZaxis().SetTitle(ztitle)
    g.GetHistogram().GetXaxis().SetTitleSize(0.05)
    g.GetHistogram().GetYaxis().SetTitleSize(0.05)
    g.GetHistogram().GetZaxis().SetTitleSize(0.05)
    g.GetHistogram().GetXaxis().SetLabelSize(0.04)
    #    g.GetHistogram().SetLineWidth(1)
    #    g.GetHistogram().SetLineColor(ROOT.kBlack)
    g.GetHistogram().GetYaxis().SetLabelSize(0.04)
    g.GetHistogram().GetXaxis().SetTitleOffset(0.8)
    g.GetHistogram().GetYaxis().SetTitleOffset(0.9)
    g.GetHistogram().GetZaxis().SetTitleOffset(0.9)
    g.GetHistogram().GetZaxis().SetRangeUser(0,10)
    h0=g.GetHistogram().Clone()
    h0.SetContour(20)
    h0.SetLineColor(linecolor)
    h0.SetLineWidth(2)
    h1=h0.Clone()
    h0.Draw('cont4z')
    pad1.Update()
    pad1.Modified()
    c.cd()
    pad2.Draw();
    pad2.cd();
    h1.SetContour(1)
    h1.SetContourLevel(0,cl1)
    h1.SetLineColor(linecolor)
    h1.SetLineWidth(3)
    h1.Draw('same cont3')
    h2=g.GetHistogram().Clone('h2')
    h2.SetContour(1)
    h2.SetContourLevel(0,cl2)
    h2.SetLineColor(linecolor)
    h2.SetLineStyle(2)
    h2.SetLineWidth(3)
    h2.Draw('cont3 same')
    pad2.Update()
    c.Draw()
    printAndSave(c,outfile,outfolder)
    return c


def canvasPostfit(name,data_prefit,pred_prefit,pred_postfit,outfile,outfolder):
    c=ROOT.TCanvas(name,name,1024,1024 )
    c.SetRightMargin(0.02)
    c.SetTopMargin(0.08)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.12)
    c.SetTicks(1,1)
    c.cd()
    ymax=0.
    for h in data_prefit,pred_prefit,pred_postfit:
        for i in range(h.GetNbinsX()):
            ymax=max(ymax,h.GetBinContent(i+1)+h.GetBinError(i+1))
    data_prefit=ROOT.TGraphErrors(data_prefit)
    pred_prefit=ROOT.TGraphErrors(pred_prefit)
    pred_postfit=ROOT.TGraphErrors(pred_postfit)
    pred_prefit.SetFillColor(ROOT.kBlue)
    pred_prefit.SetLineColor(ROOT.kWhite)
    pred_prefit.SetMarkerColor(ROOT.kBlue+1)
    pred_prefit.SetFillStyle(3354)
    pred_prefit.SetTitle('')
    pred_prefit.GetXaxis().SetTitle('')
    pred_prefit.SetMarkerStyle(21)
    pred_prefit.SetMarkerSize(1.3)
    pred_prefit.Draw('PA2')
    pred_prefit.GetYaxis().SetRangeUser(0,ymax*1.1)
    pred_postfit.SetFillColor(ROOT.kRed)
    pred_postfit.SetMarkerColor(ROOT.kRed+1)
    pred_postfit.SetLineColor(ROOT.kWhite)
    pred_postfit.SetFillStyle(3345)
    pred_postfit.SetMarkerStyle(22)
    pred_postfit.SetMarkerSize(1.6)
    pred_postfit.Draw('P2')
    data_prefit.SetLineColor(ROOT.kBlack)
    data_prefit.SetMarkerStyle(20)
    data_prefit.SetLineWidth(2)
    data_prefit.Draw('P')
    legend=ROOT.TLegend(0.7,0.65,0.9,0.9)
    legend.SetBorderSize(0)
    legend.SetLineStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.04)
    legend.SetFillStyle(0)
    legend.AddEntry(pred_prefit,'prefit','PF')
    legend.AddEntry(pred_postfit,'postfit','PF')
    legend.AddEntry(data_prefit,'measurement','PE')
    legend.Draw()

    printAndSave(c,outfile,outfolder)
    return c

def canvasPulls(name,pars,pulls,outfile,outfolder,pois=False,rescale=False):
    pars=[t for t in pars if t in pulls]
    if len(pars)==0:
        return ROOT.TCanvas(name)
    c=ROOT.TCanvas(name,name,int(50*(len(pars))),200 )
    c.SetRightMargin(0.05)
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.1)
    c.SetBottomMargin(0.15)
    c.SetTicks(1,1)
    h=ROOT.TH1F('dummy','',len(pars),-0.5,len(pars)-0.5)
    if pois and rescale:
        rescaled={}
        for p in pars:
            rescaled[p]=1.
            while pulls[p][2]-pulls[p][1]!=0 and pulls[p][2]-pulls[p][1]<20:
                pulls[p]=(pulls[p][0]*10,pulls[p][1]*10,pulls[p][2]*10)
                rescaled[p]=rescaled[p]*10
            while pulls[p][2]-pulls[p][1]!=0 and pulls[p][2]-pulls[p][1]>20:
                pulls[p]=(pulls[p][0]/10,pulls[p][1]/10,pulls[p][2]/10)
                rescaled[p]=rescaled[p]/10
                
    for i,p in enumerate(pars):
        h.GetXaxis().SetBinLabel(i+1,p)
        if rescale and rescaled[p]!=1.:
            h.GetXaxis().SetBinLabel(i+1,p+'#times'+str(Utils.round2(rescaled[p])))
    h.SetStats(0)
    h.SetLineColor(ROOT.kBlack)
    h.Draw()
    if not pois:
        h.GetYaxis().SetRangeUser(-3,3)
    elif rescale:
        h.GetYaxis().SetRangeUser(-15,15)
    else:
        ymin=0.
        ymax=0.
        for p in pulls:
            ymin=min(ymin,pulls[p][0]+pulls[p][1]*1.2)
            ymax=max(ymax,pulls[p][0]+pulls[p][2]*1.2)
        h.GetYaxis().SetRangeUser(ymin,ymax)
    if not pois:
        h.GetYaxis().SetTitle('pull')
    h.GetYaxis().CenterTitle()
#    h.GetYaxis().SetTitleSize(0.06)
#    h.GetYaxis().SetTitleOffset(0.5)
    h.GetXaxis().SetLabelSize(0.06)
    if not pois:
        outer=ROOT.TGraphAsymmErrors(h)
        inner=ROOT.TGraphAsymmErrors(h)
        for i in range(outer.GetN()):
            outer.SetPointError(i,1,1,2,2)
            inner.SetPointError(i,1,1,1,1)
        outer.SetFillColor(ROOT.kYellow)
        inner.SetFillColor(ROOT.kGreen)
        outer.Draw('2')
        inner.Draw('2')
    
    g=ROOT.TGraphAsymmErrors()
    g.SetMarkerStyle(20)
    g.SetMarkerSize(0.5)
    for i,t in enumerate(pars):
        g.SetPoint(i,i,pulls[t][0])
        g.SetPointError(i,0,0,abs(pulls[t][1]),abs(pulls[t][2]))
    g.Draw('sameP')
    h.Draw('same')
    h.Draw('axissame')
    printAndSave(c,outfile,outfolder)

def canvasCorrelation(name,h,outfile,outfolder):
    c=ROOT.TCanvas(name,name,1024,1024 )
    c.SetRightMargin(0.15)
    c.SetTopMargin(0.1)
    c.SetLeftMargin(0.1)
    c.SetBottomMargin(0.1)
    c.cd()
    h.SetStats(False)
    h.GetZaxis().SetRangeUser(-1,1)
    h.Draw('COLZ')
    printAndSave(c,outfile,outfolder)
    

def printAndSave(c,outfile,outfolder):
    outfile.cd()
    c.Write()
    c.SaveAs(outfolder+'/'+c.GetName()+'.pdf')
    c.SaveAs(outfolder+'/'+c.GetName()+'.png')
    
