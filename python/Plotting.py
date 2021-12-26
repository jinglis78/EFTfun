import ROOT,math
from ROOT import Math
from array import array
from math import log10
import Utils
import logging

def graphProfile1D(name,result,outfile=None,isNLL=True,makeProb=False):
    g=ROOT.TGraph()
    result=[r for r in result if not math.isnan(r[1]) and not math.isinf(r[1])]
    for i,r in enumerate(result):
        y=r[1]
        if makeProb:
            y=Math.chisquared_cdf_c(2*y,1)
        elif isNLL:
            y*=2
        g.SetPoint(i,r[0],y)
    g.SetName(name)
    if outfile!=None:
        outfile.cd()
        g.Write()
    return g

def graphProfile2D(name,result,outfile=None,isNLL=True,makeProb=False,conts1d=False):
    g=ROOT.TGraph2D()
    result=[r for r in result if not math.isnan(r[2]) and not math.isinf(r[2])]
    for i,r in enumerate(result):
        z=r[2]
        if makeProb:
            z=Math.chisquared_cdf_c(2*z,1 if conts1d else 2)
        elif isNLL:
            z*=2
        g.SetPoint(i,r[0],r[1],z)
    g.SetName(name)
    if outfile!=None:
        outfile.cd()
        g.Write()
    return g

def graphToys1D(name,result,outfile=None):
    g=ROOT.TGraphErrors()
    result=[r for r in result if not math.isnan(r[1]) and not math.isinf(r[1])]
    for i,r in enumerate(result):
        g.SetPoint(i,r[0],r[1])
        g.SetPointError(i,0,r[2])
    g.SetName(name)
    if outfile!=None:
        outfile.cd()
        g.Write()
    return g

def graphToys2D(name,result,outfile=None):
    g=ROOT.TGraph2DErrors()
    result=[r for r in result if not math.isnan(r[2]) and not math.isinf(r[2])]
    for i,r in enumerate(result):
        g.SetPoint(i,r[0],r[1],r[2])
        g.SetPointError(i,0,0,r[3])
    g.SetName(name)
    if outfile!=None:
        outfile.cd()
        g.Write()
    return g


def canvasProfile1D(name,result,outfile,xtitle,ytitle='#minus 2#Delta log(L)',linecolor=ROOT.kRed+1,linewidth=2,linestyle=1,outfolder='.',isNLL=True,toy_result=None,linecolor_toys=ROOT.kBlue+1,linewidth_toys=2,linestyle_toys=1):
    result=sorted(result)
    if not toy_result is None:
        toy_result=sorted(toy_result)
    cl1=1
    cl2=3.841
    pval1=1-0.6827
    pval2=0.05
    c=ROOT.TCanvas(name,name,1024,1024 )
    c.SetRightMargin(0.02)
    c.SetTopMargin(0.08)
    c.SetLeftMargin(0.12)
    c.SetBottomMargin(0.12)
    c.SetTicks(1,1)
    g2=None
    if toy_result:
        g=graphToys1D(name,toy_result)
        g2=graphProfile1D(name,result,isNLL=False,makeProb=True)
    else:
        g=graphProfile1D(name,result,isNLL=isNLL)
    xs=array('f',g.GetX())
    ys=array('f',g.GetY())
    imin=-1
    ymin=999.
    ymax=-1
    imax=-1
    x0=0.
    x1=0.
    for i,y in enumerate(ys):
        if y<ymin:
            ymin=y
            imin=i
            x0=xs[i]
        if y>ymax:
            ymax=y
            imax=i
            x1=xs[i]
    if isNLL:
        if imin>0 and ymin<0.01 and ymin>-0.001:
            g_inverse1=ROOT.TGraph(imin,ys[:imin],xs[:imin])
            g_inverse2=ROOT.TGraph(g.GetN()-imin,ys[imin:],xs[imin:])
            ci1=g_inverse1.Eval(cl1),g_inverse2.Eval(cl1)
            ci2=g_inverse1.Eval(cl2),g_inverse2.Eval(cl2)
    if not toy_result is None:
        if imax>0:
            g_inverse1=ROOT.TGraph(imax,ys[:imax],xs[:imax])
            g_inverse2=ROOT.TGraph(g.GetN()-imax,ys[imax:],xs[imax:])
            ci1=g_inverse1.Eval(pval1),g_inverse2.Eval(pval1)
            ci2=g_inverse1.Eval(pval2),g_inverse2.Eval(pval2)

    if not g2 is None:
        g2.SetLineColor(linecolor)
        g2.SetLineStyle(linestyle)
        g2.SetLineWidth(linewidth)
        g2.Draw('AL')

        g.SetLineColor(linecolor_toys)
        g.SetLineStyle(linestyle_toys)
        g.SetLineWidth(linewidth_toys)
        g.Draw('L')
    else:
        g.SetLineColor(linecolor)
        g.SetLineStyle(linestyle)
        g.SetLineWidth(linewidth)
        g.Draw('AL')
    if isNLL :
        g.GetHistogram().GetYaxis().SetRangeUser(0.,min(100,g.GetHistogram().GetYaxis().GetXmax()))
    if not toy_result is None:
        g2.GetHistogram().GetYaxis().SetRangeUser(0.,1.2)
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
    if not toy_result is None:
        xmax=g.GetHistogram().GetXaxis().GetXmax()
        xmin=g.GetHistogram().GetXaxis().GetXmin()
        cl1line=ROOT.TLine(xmin,pval1,xmax,pval1)
        cl2line=ROOT.TLine(xmin,pval2,xmax,pval2)
        cl1line.Draw()
        cl2line.Draw()
    if (isNLL and imin>0 and ymin<0.01 and ymin>-0.001) or (not toy_result is None and imax>0):
        ci2line1=ROOT.TLine(ci2[0],0,ci2[0],cl2 if isNLL else pval2)
        ci2line1.SetLineStyle(2)
        ci2line2=ROOT.TLine(ci2[1],0,ci2[1],cl2 if isNLL else pval2)
        ci2line2.SetLineStyle(2)
        ci1line1=ROOT.TLine(ci1[0],0,ci1[0],cl1 if isNLL else pval1)
        ci1line1.SetLineStyle(2)
        ci1line2=ROOT.TLine(ci1[1],0,ci1[1],cl1 if isNLL else pval1)
        ci1line2.SetLineStyle(2)
        if isNLL:
            err_low=abs(ci1[0]-x0)
            err_high=abs(ci1[1]-x0)
            central=round(x0,1-int(log10(err_low/2+err_high/2)))
            if central==0.:
                central=0
            ci1text=ROOT.TLatex(0.5,0.85,'Best fit: {}_{{ #minus{} }}^{{ #plus{} }}'.format(central,Utils.round2(err_low,2),Utils.round2(err_high,2)))
            ci1text.SetNDC()
            ci1text.SetTextSize(0.04)
            ci1text.SetTextFont(42)
            ci1text.SetTextAlign(21)
            ci1text.Draw()
        if not toy_result is None and imax>0:
            ci1text=ROOT.TLatex(0.35,0.85,'68% CL: [{},{}]'.format(Utils.round2(ci1[0],2),Utils.round2(ci1[1],2)))
            ci1text.SetNDC()
            ci1text.SetTextSize(0.04)
            ci1text.SetTextFont(42)
            ci1text.SetTextAlign(21)
            ci1text.Draw()
        if not toy_result is None:
            ci2text=ROOT.TLatex(0.35,0.79,'95% CL: [{},{}]'.format(Utils.round2(ci2[0],2),Utils.round2(ci2[1],2)))
        else:
            ci2text=ROOT.TLatex(0.5,0.79,'95% CL: [{},{}]'.format(Utils.round2(ci2[0],2),Utils.round2(ci2[1],2)))
        ci2text.SetNDC()
        ci2text.SetTextSize(0.04)
        ci2text.SetTextFont(42)
        ci2text.SetTextAlign(21)
        ci1line1.Draw()
        ci1line2.Draw()
        ci2line1.Draw()
        ci2line2.Draw()
        ci2text.Draw()
    if g2 is not None:
        legend=ROOT.TLegend(0.6,0.75,0.9,0.9)
        legend.SetBorderSize(0)
        legend.SetLineStyle(0)
        legend.SetTextFont(42)
        legend.SetTextSize(0.04)
        legend.SetFillStyle(0)
        legend.AddEntry(g,'Prob. from toys','LE')
        legend.AddEntry(g2,'Prob. asymptotic','L')
        legend.Draw()

    c.Draw()
    
    logging.info('68% CL{}: [{} , {}]'.format('' if toy_result is None else ' (toys)',ci1[0],ci1[1]))
    logging.info('95% CL{}: [{} , {}]'.format('' if toy_result is None else ' (toys)',ci2[0],ci2[1]))

    printAndSave(c,outfile,outfolder)
    if (not toy_result is None and imax>0) or isNLL and imin>0 and ymin<0.01 and ymin>-0.001:
        ftxt=open(outfolder+'/'+c.GetName()+'.txt','w')
        ftxt.write('{} {},{} {},{}'.format(x0,ci1[0],ci1[1],ci2[0],ci2[1]))

    if isNLL and ymin<-0.001:
        logging.warning("Fit was unable to find global minimum of likelihood, vary the parameter starting values (with '-s PAR=VAL')")
    return c

def canvasProfile2D(name,result,outfile,xtitle,ytitle,ztitle='#minus 2#Delta log(L)',linecolor=ROOT.kBlack,outfolder='.',isNLL=True,conts1d=False,toy_result=None,linecolor_toys=ROOT.kRed+1):
    ROOT.gStyle.SetPalette(ROOT.kThermometer);
    ROOT.TColor.InvertPalette();

    if conts1d:
        cl1=1
        cl2=3.841
    else:
        cl1=2.297
        cl2=5.991
    if not toy_result is None:
        cl1=1-0.6827
        cl2=0.05
    if not toy_result is None:
        c=ROOT.TCanvas(name,name,1280,1024)
        c.SetRightMargin(0.28)
        c.SetLeftMargin(0.1)
    else:
        c=ROOT.TCanvas(name,name,1024,1024)
        c.SetRightMargin(0.18)
        c.SetLeftMargin(0.08)
    c.SetTopMargin(0.12)
    c.SetBottomMargin(0.12)
    c.SetTicks(1,1)
    if not toy_result is None:
        g=graphToys2D(name,toy_result)
        g_asympt=graphProfile2D(name,result,makeProb=True,conts1d=conts1d)
    else:
        g=graphProfile2D(name,result)
        g.Draw('cont4')
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
    if not toy_result is None:
        g.Draw()
        g_asympt.Draw()
        g.GetHistogram().GetZaxis().SetRangeUser(0,1)
        h1=g.GetHistogram().Clone()
        h2=h1.Clone()
        h1.SetContour(1)
        h1.SetContourLevel(0,cl1)
        h1.SetLineColor(linecolor_toys)
        h1.SetLineWidth(3)
        h1.Draw('cont3')
        h2=g.GetHistogram().Clone('h2')
        h2.SetContour(1)
        h2.SetContourLevel(0,cl2)
        h2.SetLineColor(linecolor_toys)
        h2.SetLineStyle(2)
        h2.SetLineWidth(3)
        h2.Draw('cont3 same')
        hasympt1=g_asympt.GetHistogram().Clone()
        hasympt2=hasympt1.Clone()
        hasympt1.SetContour(1)
        hasympt1.SetContourLevel(0,cl1)
        hasympt1.SetLineColor(linecolor)
        hasympt1.SetLineWidth(3)
        hasympt1.Draw('cont3 same')
        hasympt2=hasympt1.Clone()
        hasympt2.SetContour(1)
        hasympt2.SetContourLevel(0,cl2)
        hasympt2.SetLineColor(linecolor)
        hasympt2.SetLineStyle(2)
        hasympt2.SetLineWidth(3)
        hasympt2.Draw('cont3 same')
        legend=ROOT.TLegend(0.72,0.6,0.99,0.88)
        legend.SetBorderSize(0)
        legend.SetLineStyle(0)
        legend.SetTextFont(42)
        legend.SetTextSize(0.04)
        legend.SetFillStyle(0)
        legend.AddEntry(h1,'68% CL toys','L')
        legend.AddEntry(h2,'95% CL toys','L')
        legend.AddEntry(hasympt1,'68% CL asymp.','L')
        legend.AddEntry(hasympt2,'95% CL asymp.','L')
        legend.Draw()

    else:
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


def canvasPostfit(name,data_prefit,data_postfit,pred_prefit,pred_postfit,outfile,outfolder):
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
    if data_postfit!=None:
        data_postfit=ROOT.TGraphErrors(data_postfit)
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
    pred_prefit.GetXaxis().SetTitle("Bins of "+name.replace("c_postfit_",""))
    pred_postfit.SetFillColor(ROOT.kRed)
    pred_postfit.SetMarkerColor(ROOT.kRed+1)
    pred_postfit.SetLineColor(ROOT.kWhite)
    pred_postfit.SetFillStyle(3345)
    pred_postfit.SetMarkerStyle(22)
    pred_postfit.SetMarkerSize(1.6)
    pred_postfit.Draw('P2')
    if data_postfit!=None:
        data_prefit.SetLineColor(ROOT.kGray+2)
        data_prefit.SetMarkerColor(ROOT.kGray+2)
        data_postfit.SetLineColor(ROOT.kBlack)
        data_postfit.SetMarkerColor(ROOT.kBlack)
        data_postfit.SetMarkerStyle(20)
        data_postfit.SetLineWidth(2)
        data_prefit.SetLineStyle(2)
    else:
        data_prefit.SetLineColor(ROOT.kBlack)
        data_prefit.SetMarkerColor(ROOT.kBlack)

    data_prefit.SetMarkerStyle(20)
    data_prefit.SetLineWidth(2)
    data_prefit.Draw('P')
    if data_postfit!=None:
        data_postfit.Draw('P')

    legend=ROOT.TLegend(0.6,0.65,0.9,0.9)
    legend.SetBorderSize(0)
    legend.SetLineStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.04)
    legend.SetFillStyle(0)
    legend.AddEntry(pred_prefit,'prefit','PF')
    legend.AddEntry(pred_postfit,'postfit','PF')
    legend.AddEntry(data_prefit,'measurement','PE')
    if data_postfit!=None:
        legend.AddEntry(data_postfit,'postfit measurement','PE')
    legend.Draw()

    printAndSave(c,outfile,outfolder,alsoLog=True)
    return c

def canvasPulls(name,pars,pulls,outfile,outfolder,pois=False,rescale=False):
    pars=[t for t in pars if t in pulls]
    if len(pars)==0:
        return ROOT.TCanvas(name)
    l=min(3,len(pars))
    c=ROOT.TCanvas(name,name,int(200*(l)),400 )
    c.SetRightMargin(0.05/l)
    c.SetTopMargin(0.05)
    c.SetLeftMargin(0.3/l)
    c.SetBottomMargin(0.15)
    c.SetTicks(1,1)
    h=ROOT.TH1F('dummy','',len(pars),-0.5,len(pars)-0.5)
    if pois and rescale:
        rescaled={}
        for p in pars:
            rescaled[p]=1.
            while pulls[p][2]-pulls[p][1]!=0 and pulls[p][2]-pulls[p][1]<2:
                pulls[p]=(pulls[p][0]*10,pulls[p][1]*10,pulls[p][2]*10)
                rescaled[p]=rescaled[p]*10
            while pulls[p][2]-pulls[p][1]!=0 and pulls[p][2]-pulls[p][1]>2:
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
#    elif rescale:
#        h.GetYaxis().SetRangeUser(-2.5,2.5)
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

def canvasCorrelation(name,h,outfile,outfolder,log=False):
    ROOT.gStyle.SetPalette(ROOT.kThermometer);
    ROOT.gStyle.SetPaintTextFormat("4.2f")
    c=ROOT.TCanvas(name,name,1024,1024 )
    c.SetRightMargin(0.15)
    c.SetTopMargin(0.1)
    c.SetLeftMargin(0.1)
    c.SetBottomMargin(0.1)
    c.cd()
    h.SetStats(False)
    if not log:
        h.GetZaxis().SetRangeUser(-1,1)
    if log:
        c.SetLogz()
    h.Draw('COLZTEXT')

    printAndSave(c,outfile,outfolder)

def canvasImpactsPage(pulls,impacts,poi,outfile,outfolder,max_impact,pars,number):
    n=len(pars)
    name='c_impacts_{}_{}'.format(poi,number)
    c=ROOT.TCanvas(name,name,800,100*n )
    scaletext=500./min(800,100*n)
    c.SetTopMargin(2./(n+3))
    c.SetBottomMargin(1./(n+3))
    c.SetRightMargin(0.05)
    c.SetLeftMargin(0.3)    
    c.cd()
    legend=ROOT.TLegend(0.1,1-(2./(n+3))*0.4,0.95,1)
    legend.SetBorderSize(0)
    legend.SetLineStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.06*scaletext)
    legend.SetFillStyle(0)
    legend.SetNColumns(3)

    hdown=ROOT.TH1D('h_down_impacts_{}_{}'.format(poi,number),'',n,0,n)
    hup=ROOT.TH1D('h_up_impacts_{}_{}'.format(poi,number),'',n,0,n)
    xmax=2.5
    relmax=0.8
    rescale=relmax*xmax/max_impact
    for i,imp in enumerate(pars):
        print imp,impacts[imp]
        hdown.SetBinContent(n-i,impacts[imp][0]*rescale)
        hdown.GetXaxis().SetBinLabel(n-i,imp)
        hup.SetBinContent(n-i,impacts[imp][1]*rescale)
        hup.GetXaxis().SetBinLabel(n-i,imp)
       
    hup.SetStats(0)
    hup.SetLabelSize(0.09*scaletext,'X')
    hup.SetLabelSize(0.06*scaletext,'Y')
    hup.SetTitleSize(0.06*scaletext,'XY')

    hup.GetYaxis().SetTitle('(#hat{#theta}-#theta_{0})')
    hup.GetYaxis().CenterTitle(1)
    hup.SetFillColor(ROOT.kRed)
    hdown.SetFillColor(ROOT.kBlue)
    hup.Draw('hbar')
    
    line1=ROOT.TLine(1,0,1,n)
    line1.Draw()
    line1.SetLineStyle(2)
    linem1=ROOT.TLine(-1,0,-1,n)
    linem1.SetLineStyle(2)
    linem1.Draw()

    hup.Draw('hbar,same')

    hup.SetBarWidth(0.33)
    hup.SetBarOffset(0.1)
    hup.GetYaxis().SetRangeUser(-xmax,xmax)
    hdown.SetBarWidth(0.33)
    hdown.SetBarOffset(0.57)
    hdown.Draw('hbar,same')
    g_pull=ROOT.TGraphAsymmErrors()
    for i in range(n):
        label=hup.GetXaxis().GetBinLabel(i+1)
        g_pull.SetPoint(i , pulls[label][0],hup.GetBinCenter(i+1))
        g_pull.SetPointError( i, abs(pulls[label][1]), abs(pulls[label][2]), 0, 0)
    g_pull.SetLineWidth(2)
    g_pull.SetMarkerStyle(21)
    g_pull.Draw('EP')

    line0=ROOT.TLine(0,0,0,n)
    line0.SetLineWidth(2)
    line0.Draw()

    
    axis = ROOT.TGaxis(-xmax, n, xmax, n,
                       -max_impact/relmax, max_impact/relmax, 510, "-")
    axis.SetLabelSize(0.06*scaletext)
    axis.SetLabelFont(42)
    axis.SetTitleSize(0.06*scaletext)
    axis.SetTitleFont(42)
    axis.SetLabelOffset(0.0)
    axis.SetTitle('#Delta '+poi)
    axis.CenterTitle(1)
    axis.Draw()

    legend.AddEntry(hup,'Impact #hat{#theta}+#Delta#hat{#theta}','F')
    legend.AddEntry(hdown,'Impact #hat{#theta}-#Delta#hat{#theta}','F')
    legend.AddEntry(g_pull,'#hat{#theta} and #pm#Delta#hat{#theta}','EP')
    legend.Draw()
    printAndSave(c,outfile,outfolder)
    
def canvasImpacts(pulls,impacts,poi,outfile,outfolder):
    ROOT.gStyle.SetHistMinimumZero(True)
    ROOT.gStyle.SetEndErrorSize(4)
    impacts_ordered=sorted([( max(abs(impacts[name][0]),abs(impacts[name][1])) , name) for name in impacts],reverse=True)
    max_impact=impacts_ordered[0][0]
    pars=[x[1] for x in impacts_ordered]
    per_page=10
    for n in range((len(pars)-1)//per_page+1):
        canvasImpactsPage(pulls,impacts,poi,outfile,outfolder,max_impact,pars[n*per_page:(n+1)*per_page],n)
                         
    

def printAndSave(c,outfile,outfolder,alsoLog=False):
    outfile.cd()
    c.Write()
    c.SaveAs(outfolder+'/'+c.GetName()+'.pdf')
    c.SaveAs(outfolder+'/'+c.GetName()+'.png')
    if alsoLog:
        c.SetLogy()
        c.SaveAs(outfolder+'/'+c.GetName()+'_log.pdf')
        c.SaveAs(outfolder+'/'+c.GetName()+'_log.png')
