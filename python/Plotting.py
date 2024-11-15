import ROOT,math
from ROOT import Math
from array import array
from math import log10
import Utils
import logging




def graphProfile1D(name,result,outfile=None,isNLL=True,makeProb=False,conts2d=False):
    g=ROOT.TGraph()
    result=[r for r in result if not math.isnan(r[1]) and not math.isinf(r[1])]
    for i,r in enumerate(result):
        y=r[1]
        if makeProb:
            y=Math.chisquared_cdf_c(2*y,2 if conts2d else 1)
        elif isNLL:
            y*=2
        g.SetPoint(i,r[0],y)
    g.SetName(name)
    if outfile!=None:
        outfile.cd()
        g.Write()
    return g


def make_line(formula, op, nu_val):
    # Create a dictionary to map variable names to values
    variables = {'cll1': 0, 'cHq3': 0, 'clq3': 0, 'cHl3': 0 }

    variables.update({
        'hist_sm': 0,
        'hist_cll1': 0,
        'hist_cHq3': 0,
        'hist_clq3': 0,
        'hist_cHl3': 0,
        'hist_cll1_quad': 0,
        'hist_cHq3_quad': 0,
        'hist_clq3_quad': 0,
        'hist_cHl3_quad': 0,
        'hist_cHl3_cll1_cross': 0,
        'hist_cHq3_clq3_cross': 0,
        'hist_cHq3_cHl3_cross': 0,
        'hist_cHq3_cll1_cross': 0,
        'hist_cll1_clq3_cross': 0,
        'hist_cHl3_clq3_cross': 0,
        'stat':0,
        'alphas':0,
        'scaleMu':0,
        'pdf':0,
        'ew':0
    })

    # Set the value of the relevant nuisance parameter
    variables["hist_"+op] = nu_val

    return variables

def canvasPredFormula(name, formula, data, uncertainties, pulls, outfile, outfolder):
    
    # Check for SetAtlasStyle and call it
    if hasattr(ROOT, "SetAtlasStyle"):
        ROOT.SetAtlasStyle()
    else:
        print("SetAtlasStyle is not defined in AtlasStyle.C")

    operators = {'cll1', 'cHq3', 'clq3', 'cHl3'} 

    # Dictionary of multigraphs

    multigraphs = {}

    min_result = {}
    max_result = {}

    for itr,op in enumerate(operators):
        print(op)

        # make multi graph
        mg = ROOT.TMultiGraph()
        
        mg.SetTitle("{};{}; Cross section [nb]".format(op, op))

        min_result[op] = 0
        max_result[op] = 1


        for nu_val in [-1, 0, 1]:
            print("nu_val = ", nu_val)
            print("op = ", op)

            results = []

            variables = make_line(formula, op, nu_val)

            graph = ROOT.TGraph()
            graph.GetXaxis().SetRangeUser(-0.2, 0.2)
            
            if op == "clq3":
                graph.GetXaxis().SetRangeUser(-0.001, 0.001)

            line_color = ROOT.kBlue if nu_val == 1 else ROOT.kRed if nu_val == -1 else ROOT.kBlack
            line_style = 2 if nu_val != 0 else 1

            graph.SetLineColor(line_color)
            graph.SetLineStyle(line_style)
            
            
            variables_sm = variables.copy()
            nominal = eval(formula, variables_sm)

            # Loop over x in the range -1 to 1 with 100 steps
            for x in range(-20, 21):
                variables_copy = variables.copy()
                variables_copy[op] = x / 20000  # Scale x to be in the range -1 to 1: 20
                if op != "clq3": 
                    variables_copy[op] = x / 200 # Scale x to be in the range -10 to 10

                print (formula)
                print (variables_copy)
                result = eval(formula, variables_copy)
                
                rel_result = (result - nominal) / nominal * 100
            
                #if x == 20:
                    #print ("x= ", x)
                    #print(formula)
                    #print("x/20 = ", x / 20.0)
                    #print("printing variables")
                    # print (variables_copy)
                    #print(result)
                results.append(result)
                if result < min_result[op]:
                    min_result[op] = rel_result
                if result > max_result[op]:
                    max_result[op] = rel_result

            #print ("min and max results: ", min_result, " and ", max_result)
            
            if min_result[op] < 0:
                graph.GetYaxis().SetRangeUser(min_result[op], max_result[op])
            else:
                graph.GetYaxis().SetRangeUser(min_result[op], max_result[op])

            graph.SetTitle("{};X title;Y title");

            for i, result in enumerate(results):
                print (i, result)
                if op != "clq3": 
                    print ("(i - 20) / 2. = ", (i - 20) / 2.)
                    graph.SetPoint(i, (i - 20) / 200., result)
                else: 
                    graph.SetPoint(i, (i - 20) / 20000., result)

            mg.Add(graph,"lp")


        multigraphs[op]= mg
    
    for op, mg in multigraphs.items():
        print (op)

        ci = pulls[op][0]
        ci_errors = (pulls[op][1], pulls[op][2])    

        # Create a new canvas for each multigraph
        c = ROOT.TCanvas(name+"_"+op, name+"_"+op, 1024, 1024)
        c.SetRightMargin(0.1)
        c.SetTopMargin(0.1)
        c.SetLeftMargin(0.15)
        c.SetBottomMargin(0.1)
        c.cd()


        data_points = ROOT.TGraphErrors()
        xval = 0.0001
        if op != "clq3":
            xval = 0.01
            
        data_points.SetPoint(0, xval, data["meas"])  # sm meas

        postfit_points = ROOT.TGraphErrors()
        postfit_points.SetPoint(0, ci, data["postfit"] )  # postfit

        pred_points = ROOT.TGraphErrors()
        pred_points.SetPoint(0, 0, data["sm"])  # sm pred

       

        print ("ci: ", ci)
        print ("data[op]: ", data[op], " +- ",uncertainties[op][1])
        print ("data['sm']: ", data['sm'], " +- ",uncertainties['sm'][1])
        print ("data['meas']: ", data['meas'], " +- ",uncertainties['meas'][1])
        print ("data['postfit']: ", data['postfit'], " +- ",uncertainties['postfit'][1])

        # Calculate maximum and minimum values
        max_value = max(data["sm"] + data[op] + uncertainties["postfit"][1] , data["postfit"] + uncertainties["postfit"][1])
        min_value = min(data["sm"] - data[op] - uncertainties["postfit"][1] ,  data["postfit"] - uncertainties["postfit"][1])

        # Set yMax and yMin
        # Final mass bin
        yMax = 0.3
        yMin = 0.1
        
        # Penultimate mass bin
        yMax = 1.4
        yMin = 0.8
    
        
        mg.SetMaximum(yMax)
        mg.SetMinimum(yMin)

        mg.Draw("a")
        
        pred_points.SetPointError(0, 0, uncertainties["sm"][1])  # Set errors for point 0
        #pred_points.SetPointError(1, 0, uncertainties[op][1] )  # Set errors for point 1
        

        data_points.SetPointError(0, 0,  uncertainties["meas"][1]) # Set errors for point 0

        postfit_points.SetPointError(0,0, uncertainties["postfit"][1])

        pred_points.SetFillColor(ROOT.kBlue)
        pred_points.SetMarkerColor(ROOT.kBlue+1)
        pred_points.SetLineColor(ROOT.kBlue+1)
        pred_points.SetFillStyle(3345)
        pred_points.SetMarkerStyle(20)
        pred_points.SetMarkerSize(2)
        pred_points.SetLineWidth(2)
        pred_points.Draw('P')

        data_points.SetFillColor(ROOT.kBlack)
        data_points.SetMarkerColor(ROOT.kBlack)
        data_points.SetLineColor(ROOT.kBlack)
        data_points.SetFillStyle(3345)
        data_points.SetMarkerStyle(20)
        data_points.SetMarkerSize(2)
        data_points.SetLineWidth(2)
        data_points.Draw('P')

        postfit_points.SetFillColor(ROOT.kRed)
        postfit_points.SetMarkerColor(ROOT.kRed+1)
        postfit_points.SetLineColor(ROOT.kRed+1)
        postfit_points.SetFillStyle(3345)
        postfit_points.SetMarkerStyle(22)
        postfit_points.SetMarkerSize(2)
        postfit_points.SetLineWidth(2)
        postfit_points.Draw('P')

        legend = ROOT.TLegend(0.65, 0.7, 0.9, 0.9)
        legend.SetTextSize(0.035)
        legend.SetFillStyle(0)
        legend.SetBorderSize(0)
        legend.AddEntry(data_points, "Meas.", "p")
        legend.AddEntry(pred_points, "Pred.", "p")
        legend.AddEntry(postfit_points, "Postfit", "p")


        legend.Draw("same")

        # Assuming ci, ci_errors[0], and ci_errors[1] are floating-point numbers
        sampleText = op + " = " + str(round(ci, 3)) + " [ +" + str(round(ci_errors[1], 3)) +  ", " + str(round(ci_errors[0], 3)) + " ]"

        if op == "clq3":    
            sampleText = op + " = " + str(round(ci, 6)) + " [ +" + str(round(ci_errors[1], 6)) +  ", " + str(round(ci_errors[0], 6)) + " ]"

        myText(0.2, 0.8, 1, 0.03, sampleText )

        c.Update()
        c.WaitPrimitive()
        printAndSave(c, outfile, outfolder)

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


def canvasProfile1D(name,result,outfile,xtitle,ytitle='#minus 2#Delta log(L)',linecolor=ROOT.kRed+1,linewidth=2,linestyle=1,outfolder='.',isNLL=True,toy_result=None,linecolor_toys=ROOT.kBlue+1,linewidth_toys=2,linestyle_toys=1,conts2d=False):
    result=sorted(result)
    if not toy_result is None:
        toy_result=sorted(toy_result)
    if conts2d:
        cl1=2.297
        cl2=5.991
    else:
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
        g2=graphProfile1D(name,result,isNLL=False,makeProb=True,conts2d=conts2d)
    else:
        g=graphProfile1D(name,result,isNLL=isNLL,conts2d=conts2d)
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

    if (isNLL and imin>0 and ymin<0.01 and ymin>-0.001) or (not toy_result is None and imax>0):
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
    
    SetAtlasStyle()
    

    # Setup Canvas 
    x1 = 0.04
    x2 = 0.96
    y1 = 0.30
    y2 = 0.96
    gap = 0.05
    
    aScale = (y2-y1)/(y1-gap-x1)

    c = ROOT.TCanvas(name,name,1000,1000)
    c.SetLogy()
    c.Draw()
    
    c.SetTicks(1, 1)

    c_2 = ROOT.TPad("c ratio", "c ratio", x1, x1+gap, x2, y1)
    c_2.Draw()
    c_2.SetLogy(0)
    c_2.SetTopMargin(0.02)
    c_2.SetBottomMargin(0.4) 
    c_2.SetFillStyle(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    

    c_1 = ROOT.TPad("c top", "c top", x1, y1, x2, y2)
    c_1.Draw()
    c_1.SetLogy(0)
    c_1.SetTopMargin(0.1)
    c_1.SetBottomMargin(0.02)
    c_1.SetFillStyle(0)
    


    print ("Canvas set up")

    c_1.cd()
    
    sampleText = r"pp \rightarrow l^{\pm} \nu"
    if "elmu" not in name and "el" in name:
        sampleText = r"pp \rightarrow e^{\pm} \nu"
    if "elmu" not in name and "mu" in name:
        sampleText = r"pp \rightarrow \mu^{\pm} \nu"
        

    ymax=0.
    for h in data_prefit,pred_prefit,pred_postfit:
        print (h.GetNbinsX())
        for i in range(h.GetNbinsX()):
            ymax=max(ymax,h.GetBinContent(i+1)+h.GetBinError(i+1))
        
    '''

    data_prefit_graph=ROOT.TGraphErrors(data_prefit)
    if data_postfit!=None:
        data_postfit_graph=ROOT.TGraphErrors(data_postfit)
    pred_postfit_graph=ROOT.TGraphErrors(pred_postfit)

    '''

    pred_prefit.SetFillColor(ROOT.kBlue)
    pred_prefit.SetLineColor(ROOT.kWhite)
    pred_prefit.SetMarkerColor(ROOT.kBlue+1)
    pred_prefit.SetLineWidth(2)
    pred_prefit.SetMarkerSize(1.6)
    pred_prefit.SetMarkerStyle(22)
    #pred_prefit.SetFillStyle(3354)
    pred_prefit.SetTitle('')
    pred_prefit.GetXaxis().SetTitle('')
    

    y_title = ""
    x_title = ""
    labels = []

    num_bins = pred_prefit.GetNbinsX()

    if num_bins < 15:
        y_title = r'd\sigma/dM [pb/GeV]'
        x_title = r"M_{T} [GeV]"
        if num_bins == 11:
            labels = ["200-250", "250-300", "300-350", "350-425", "425-500", "500-600", "600-750", "750-900", "900-1100", "1100-1400", "1400-2000"]
        elif num_bins == 12:
            labels = ["200-250", "250-300", "300-350", "350-425", "425-500", "500-600", "600-750", "750-900", "900-1100", "1100-1400", "1400-2000", "2000-5000"]

            
    else:
        y_title = r'd^{2}\sigma / dMd d\eta [pb/GeV]'
        x_title = r"|\eta^{lepton}| \otimes M_{T} [GeV]"

        if num_bins == 40:
            labels = ["200-250", "300-425", "425-600", "600-900", "900-2000"]
        elif num_bins == 44:
            labels = ["200-300", "300-425", "425-600", "600-900", "900-2000","2000-5000"]


    pred_prefit.GetYaxis().SetTitle(y_title)

    pred_prefit.GetXaxis().SetTitleSize(0.1)
    pred_prefit.GetXaxis().SetTitleOffset(0)
    pred_prefit.GetXaxis().SetLabelSize(0)
    pred_prefit.GetYaxis().SetTitleOffset(1.6)
    pred_prefit.GetYaxis().SetTitleSize(0.05)
    pred_prefit.GetYaxis().SetLabelSize(0.05)
    pred_prefit.SetMaximum(ymax*1.1)
    pred_prefit.SetStats(0)

    pred_prefit.DrawCopy("e")

    print (num_bins)

    #pred_prefit=ROOT.TGraphErrors(pred_prefit)

    pred_prefit.Draw('E2 same')
    pred_postfit.SetFillColor(ROOT.kRed)
    pred_postfit.SetMarkerColor(ROOT.kRed+1)
    pred_postfit.SetLineColor(ROOT.kWhite)
    #pred_postfit.SetFillStyle(3345)
    pred_postfit.SetMarkerStyle(22)
    pred_postfit.SetMarkerSize(1.6)
    pred_postfit.Draw('E2 same')
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
    data_prefit.Draw('PE same')
    if data_postfit!=None:
        data_postfit.Draw('PE same')

    leg = ROOT.TLegend(0.68,0.58,0.90,0.79)
    leg.SetTextSize(0.035)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)

    if  num_bins > 15:
        print("Plotting lines")
        line1 = ROOT.TLine(0.3, 0, 0.3, 0.08)
        line2 = ROOT.TLine(0.6, 0, 0.6, 0.08)
        line3 = ROOT.TLine(0.75, 0, 0.75, 0.08)
        line4 = ROOT.TLine(0.9, 0, 0.9, 0.08)

        line1.SetLineWidth(2)  # Set line thickness to 2 pixels
        line2.SetLineWidth(2)
        line3.SetLineWidth(2)
        line4.SetLineWidth(2)

        line1.Draw("same")
        line2.Draw("e")
        line3.Draw("p")
        line4.Draw("p")

    leg.AddEntry(pred_prefit,'Prefit pred.','pe')
    leg.AddEntry(pred_postfit,'Postfit pred.','pe')
    leg.AddEntry(data_prefit,'Measurement','pe')
    if data_postfit!=None:
        leg.AddEntry(data_postfit,'Postfit meas.','PE')
    leg.Draw("same")

    ATLASLabel(0.50, 0.83, "Internal")
    myText(0.5, 0.75, 1, 0.042, sampleText )

    c_2.cd()

    prefit_ratio = pred_prefit.Clone()
    prefit_ratio.Divide(data_prefit)
    prefit_ratio.GetYaxis().SetNdivisions(505)
    prefit_ratio.GetYaxis().CenterTitle()
    prefit_ratio.SetFillColor(ROOT.kBlue)
    prefit_ratio.SetLineColor(ROOT.kWhite)
    prefit_ratio.SetMarkerColor(ROOT.kBlue+1)
    prefit_ratio.SetTitle('')

    num_points = prefit_ratio.GetNbinsX()

    # New section: Shift the prefit points slightly to the left
    x_shift = -0.02  # Shift amount to move the prefit points left
    shifted_prefit = ROOT.TGraphErrors(num_points)  # Create a new TGraphErrors to hold shifted points
    for i in range(num_points):
        x_value = prefit_ratio.GetXaxis().GetBinCenter(i+1)
        y_value = prefit_ratio.GetBinContent(i+1)
        y_error = prefit_ratio.GetBinError(i+1)
        shifted_prefit.SetPoint(i, x_value + x_shift, y_value)  # Shift x-value
        shifted_prefit.SetPointError(i, 0, y_error)

    # Update the graph to reflect the changes
    print(aScale*pred_prefit.GetXaxis().GetTitleSize())

    shifted_prefit.SetMarkerStyle(21)
    shifted_prefit.SetMarkerSize(1.3)
    shifted_prefit.SetFillColor(ROOT.kBlue)
    shifted_prefit.SetMarkerColor(ROOT.kBlue+1)
    shifted_prefit.SetLineColor(ROOT.kWhite)

    prefit_ratio.SetMarkerStyle(21)
    prefit_ratio.SetMarkerSize(1.3)
    prefit_ratio.GetYaxis().SetTitleSize(0.10)
    prefit_ratio.GetYaxis().CenterTitle()
    prefit_ratio.GetYaxis().SetTitleOffset(1.6)
    prefit_ratio.GetYaxis().SetLabelSize(0.08)
    prefit_ratio.GetXaxis().SetLabelSize(0)
    prefit_ratio.SetLineWidth(5)

    if data_postfit != None:
        postfit_ratio = pred_postfit.Clone("postfit_ratio")
        postfit_ratio.Divide(data_postfit) 
        postfit_ratio = ROOT.TGraphErrors(postfit_ratio)
        postfit_ratio.SetFillColor(ROOT.kRed)
        postfit_ratio.SetLineColor(ROOT.kWhite)
        postfit_ratio.SetMarkerColor(ROOT.kRed+1)
        postfit_ratio.SetTitle('')
        postfit_ratio.SetMarkerStyle(21)
        postfit_ratio.SetMarkerSize(1.3)
        
        data_ratio = data_postfit.Clone("data_ratio")
        data_ratio.Divide(data_postfit) 
        #data_ratio = ROOT.TGraphErrors(data_ratio)
        data_ratio.SetFillColor(ROOT.kBlack)
        data_ratio.SetLineColor(ROOT.kWhite)
        data_ratio.SetMarkerColor(ROOT.kBlack)
        data_ratio.SetTitle('')
        data_ratio.SetMarkerStyle(21)
        data_ratio.SetMarkerSize(1.3)
        
    else:
        postfit_ratio = pred_postfit.Clone("postfit_ratio")
        postfit_ratio.Divide(data_prefit)
        postfit_ratio = ROOT.TGraphErrors(postfit_ratio)
        postfit_ratio.SetFillColor(ROOT.kRed)
        postfit_ratio.SetLineColor(ROOT.kWhite)
        postfit_ratio.SetMarkerColor(ROOT.kRed+1)
        postfit_ratio.SetTitle('')
        postfit_ratio.SetMarkerStyle(21)
        postfit_ratio.SetMarkerSize(1.3)
        postfit_ratio.SetLineWidth(5)
        
        data_ratio = data_prefit.Clone("data_ratio")
        data_ratio.Divide(data_prefit) 
        #data_ratio = ROOT.TGraphErrors(data_ratio)
        data_ratio.SetFillColor(ROOT.kBlack)
        data_ratio.SetLineColor(ROOT.kWhite)
        data_ratio.SetMarkerColor(ROOT.kBlack)
        data_ratio.SetTitle('')
        data_ratio.SetMarkerStyle(21)
        data_ratio.SetMarkerSize(1.3)
        
    # New section: Shift the data points slightly to the left
   
    shifted_data = ROOT.TGraphErrors(num_points)  # Create a new TGraphErrors to hold shifted points
    for i in range(num_points):
        x_value = data_ratio.GetXaxis().GetBinCenter(i+1)
        y_value = data_ratio.GetBinContent(i+1)
        y_error = data_ratio.GetBinError(i+1)
        shifted_data.SetPoint(i, x_value + 0.02, y_value)  # Shift x-value
        shifted_data.SetPointError(i, 0, y_error)

    # Update the graph to reflect the changes
    print(aScale*pred_prefit.GetXaxis().GetTitleSize())

    shifted_data.SetMarkerStyle(21)
    shifted_data.SetMarkerSize(1.3)
    shifted_data.SetFillColor(ROOT.kBlack)
    shifted_data.SetMarkerColor(ROOT.kBlack)
    shifted_data.SetLineColor(ROOT.kWhite)

    line = ROOT.TLine(prefit_ratio.GetXaxis().GetXmin(), 1.0, prefit_ratio.GetXaxis().GetXmax(), 1.0)
    line.SetLineWidth(2)
    line.SetLineStyle(2)
    line.SetLineColor(ROOT.kBlack)

    ratio = pred_prefit.Clone()
    ratio.Divide(pred_prefit)
    ratio.GetYaxis().SetNdivisions(505)
    ratio.GetYaxis().CenterTitle()
    ratio.GetYaxis().SetTickLength(0.02)
    ratio.SetMinimum(0.5)
    ratio.SetMaximum(1.5)
    ratio.SetTitle("")
    ratio.SetMarkerSize(0)
    ratio.SetLineWidth(0)
    ratio.GetYaxis().SetTitle("Pred. / Meas.")
    ratio.GetXaxis().SetTitleOffset(5)  
    ratio.GetXaxis().SetTitleSize(0.1)
    ratio.GetXaxis().SetTickLength(0.1)
    ratio.GetXaxis().SetTitleOffset(-0.04)
    ratio.GetYaxis().SetTitleOffset(0.6)
    ratio.GetXaxis().SetLabelSize(0.1)
    ratio.GetXaxis().SetLabelOffset(0.04)
    ratio.GetYaxis().SetTitleSize(0.12)
    ratio.GetYaxis().SetLabelSize(0.1)
    ratio.SetStats(0)

    if num_bins < 15:
        for i, label in enumerate(labels):
            ratio.GetXaxis().SetBinLabel(i+1, label)
    else:
        ratio.GetXaxis().SetLabelSize(0)

    ratio.DrawCopy("e")
    
    

    # Modified part: Use shifted prefit ratio instead of the original
    shifted_prefit.Draw("same PE2")
    postfit_ratio.Draw("same PE2")
    shifted_data.Draw("same PE2")
    line.Draw("same")
    
    if num_bins < 15:
        myText(0.85, 0.02, 1, 0.15, x_title)
    else:   
        myText(0.7, 0.02, 1, 0.15, x_title)

        # Define the axis labels and the corresponding bin ranges
        labels = ["200-300", "300-425", "425-600", "600-900", "900-2000"]
        bin_ranges = [0.24, 0.48, 0.65, 0.77, 0.865]  # Bin ranges where labels change

        for i in range(len(labels)):
            xoffset =  bin_ranges[i]
            text = labels[i]
            myText(xoffset, 0.2, 1, 0.1, text)

    c_1.Update()
    c_2.Update()
    c.Update()

    printAndSavePads(c, c_1, c_2, outfile, outfolder, alsoLog=True)
    return c
  
def canvasRatio(name,data_prefit,data_postfit,pred_prefit,pred_postfit,outfile,outfolder):

    # Check for SetAtlasStyle and call it
    if hasattr(ROOT, "SetAtlasStyle"):
        ROOT.SetAtlasStyle()
    else:
        print("SetAtlasStyle is not defined in AtlasStyle.C")

    # Setup Canvas 
    x1 = 0.04
    x2 = 0.96
    y1 = 0.30
    y2 = 0.96
    gap = 0.05
    
    aScale = (y2-y1)/(y1-gap-x1)

    c = ROOT.TCanvas(name,name,1000,1000)
    c.SetLogy()
    c.Draw()

    c_2 = ROOT.TPad("c ratio", "c ratio", x1, x1+gap, x2, y1)
    c_2.Draw()
    c_2.SetLogy(0)
    c_2.SetTopMargin(0.02)
    c_2.SetBottomMargin(0.4) 
    c_2.SetFillStyle(0)

    c_1 = ROOT.TPad("c top", "c top", x1, y1, x2, y2)
    c_1.Draw()
    c_1.SetLogy(0)
    c_1.SetTopMargin(0.1)
    c_1.SetBottomMargin(0.02)
    c_1.SetFillStyle(0)


    print ("Canvas set up")

    c_1.cd()

    sampleText = r"pp \rightarrow l^{\pm} \nu"
    if "elmu" not in name and "el" in name:
        sampleText = r"pp \rightarrow e^{\pm} \nu"
    if "elmu" not in name and "mu" in name:
        sampleText = r"pp \rightarrow \mu^{\pm} \nu"

    ymax=0.
    for h in data_prefit,pred_prefit,pred_postfit:
        print (h.GetNbinsX())
        for i in range(h.GetNbinsX()):
            ymax=max(ymax,h.GetBinContent(i+1)+h.GetBinError(i+1))
        


    pred_prefit.SetFillColor(ROOT.kBlue)
    pred_prefit.SetLineColor(ROOT.kWhite)
    pred_prefit.SetMarkerColor(ROOT.kBlue+1)
    pred_prefit.SetLineWidth(2)
    pred_prefit.SetMarkerSize(1.6)
    pred_prefit.SetMarkerStyle(22)
    #pred_prefit.SetFillStyle(3354)
    pred_prefit.SetTitle('')
    pred_prefit.GetXaxis().SetTitle('')
    

    y_title = ""
    x_title = ""
    labels = []

    num_bins = pred_prefit.GetNbinsX()

    if num_bins < 15:
        y_title = r'd\sigma/dM [pb/GeV]'
        x_title = r"M_{T} [GeV]"
        if num_bins == 11:
            labels = ["200-250", "250-300", "300-350", "350-425", "425-500", "500-600", "600-750", "750-900", "900-1100", "1100-1400", "1400-2000"]
        elif num_bins == 12:
            labels = ["200-250", "250-300", "300-350", "350-425", "425-500", "500-600", "600-750", "750-900", "900-1100", "1100-1400", "1400-2000", "2000-5000"]

            
    else:
        y_title = r'd^{2}\sigma / dMd d\eta [pb/GeV]'
        x_title = r"|\eta^{lepton}| \otimes M_{T} [GeV]"

        if num_bins == 40:
            labels = ["200-250", "300-425", "425-600", "600-900", "900-2000"]
        elif num_bins == 44:
            labels = ["200-300", "300-425", "425-600", "600-900", "900-2000","2000-5000"]


    pred_prefit.GetYaxis().SetTitle(y_title)

    pred_prefit.GetXaxis().SetTitleSize(0.1)
    pred_prefit.GetXaxis().SetTitleOffset(0)
    pred_prefit.GetXaxis().SetLabelSize(0)
    pred_prefit.GetYaxis().SetTitleOffset(1.6)
    pred_prefit.GetYaxis().SetTitleSize(0.05)
    pred_prefit.GetYaxis().SetLabelSize(0.05)
    pred_prefit.SetMaximum(ymax*1.1)

    pred_prefit.DrawCopy("e")

    print (num_bins)

    #pred_prefit=ROOT.TGraphErrors(pred_prefit)

    pred_prefit.Draw('E2 same')
    pred_postfit.SetFillColor(ROOT.kRed)
    pred_postfit.SetMarkerColor(ROOT.kRed+1)
    pred_postfit.SetLineColor(ROOT.kWhite)
    #pred_postfit.SetFillStyle(3345)
    pred_postfit.SetMarkerStyle(22)
    pred_postfit.SetMarkerSize(1.6)
    pred_postfit.Draw('E2 same')
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
    data_prefit.Draw('PE same')
    if data_postfit!=None:
        data_postfit.Draw('PE same')

    leg = ROOT.TLegend(0.68,0.58,0.90,0.79)
    leg.SetTextSize(0.035)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)

    if  num_bins > 15:
        print("Plotting lines")
        line1 = ROOT.TLine(0.3, 0, 0.3, 0.15)
        line2 = ROOT.TLine(0.6, 0, 0.6, 0.15)
        line3 = ROOT.TLine(0.75, 0, 0.75, 0.15)
        line4 = ROOT.TLine(0.9, 0, 0.9, 0.15)

        line1.SetLineWidth(2)  # Set line thickness to 2 pixels
        line2.SetLineWidth(2)
        line3.SetLineWidth(2)
        line4.SetLineWidth(2)

        line1.Draw("same")
        line2.Draw("e")
        line3.Draw("p")
        line4.Draw("p")

    leg.AddEntry(pred_prefit,'Prefit pred.','pe')
    leg.AddEntry(pred_postfit,'Postfit pred.','pe')
    leg.AddEntry(data_prefit,'Measurement','pe')
    if data_postfit!=None:
        leg.AddEntry(data_postfit,'postfit measurement','PE')
    leg.Draw("same")

    ATLASLabel(0.50, 0.83, "Work in Progress")
    myText(0.5, 0.75, 1, 0.042, sampleText )

    c_2.cd()

    prefit_ratio = pred_prefit.Clone()
    prefit_ratio.Divide(data_prefit)
    prefit_ratio.GetYaxis().SetNdivisions(505)
    prefit_ratio.GetYaxis().CenterTitle()
    #prefit_ratio = ROOT.TGraphErrors(prefit_ratio)
    prefit_ratio.SetFillColor(ROOT.kBlue)
    prefit_ratio.SetLineColor(ROOT.kWhite)
    prefit_ratio.SetMarkerColor(ROOT.kBlue+1)
    #prefit_ratio.SetFillStyle(3345)
    prefit_ratio.SetTitle('')


    num_points = prefit_ratio.GetNbinsX()


    
    # Update the graph to reflect the changes
    print(aScale*pred_prefit.GetXaxis().GetTitleSize())


    prefit_ratio.SetMarkerStyle(21)
    prefit_ratio.SetMarkerSize(1.3)
    prefit_ratio.GetYaxis().SetTitleSize(0.10)
    prefit_ratio.GetYaxis().CenterTitle()
    prefit_ratio.GetYaxis().SetTitleOffset(0.5)
    prefit_ratio.GetYaxis().SetLabelSize(0.08)
    prefit_ratio.GetXaxis().SetLabelSize(0)
    prefit_ratio.SetLineWidth(5)
    if data_postfit!=None:
        postfit_ratio = pred_postfit.Clone("postfit_ratio")
        postfit_ratio.Divide(data_postfit) 
        postfit_ratio = ROOT.TGraphErrors(postfit_ratio)
        postfit_ratio.SetFillColor(ROOT.kRed)
        postfit_ratio.SetLineColor(ROOT.kWhite)
        postfit_ratio.SetMarkerColor(ROOT.kRed+1)
        #postfit_ratio.SetFillStyle(3354)
        postfit_ratio.SetTitle('')
        postfit_ratio.SetMarkerStyle(21)
        postfit_ratio.SetMarkerSize(1.3)
    else:
        postfit_ratio = pred_postfit.Clone("postfit_ratio")
        postfit_ratio.Divide(data_prefit)
        postfit_ratio = ROOT.TGraphErrors(postfit_ratio)
        postfit_ratio.SetFillColor(ROOT.kRed)
        postfit_ratio.SetLineColor(ROOT.kWhite)
        postfit_ratio.SetMarkerColor(ROOT.kRed+1)
        #postfit_ratio.SetFillStyle(3354)
        postfit_ratio.SetTitle('')
        postfit_ratio.SetMarkerStyle(21)
        postfit_ratio.SetMarkerSize(1.3)
        postfit_ratio.SetLineWidth(5)

    line = ROOT.TLine(prefit_ratio.GetXaxis().GetXmin(), 1.0, prefit_ratio.GetXaxis().GetXmax(), 1.0)
    line.SetLineWidth(2)
    line.SetLineStyle(2)
    line.SetLineColor(ROOT.kBlack)


    hnumer = pred_prefit.Clone()
    hdenom = data_prefit.Clone()

    hnumer.Divide(hnumer, hdenom, 1.0, 1.0, "B" )

    ratio = pred_prefit.Clone()
    ratio.Divide(pred_prefit)
    ratio.GetYaxis().SetNdivisions(505)
    ratio.GetYaxis().CenterTitle()
    ratio.GetYaxis().SetTickLength(0.02)
    ratio.GetYaxis().SetTitleSize(0.1)
    ratio.SetMinimum(0.8)
    ratio.SetMaximum(1.2)
    ratio.SetTitle("")
    ratio.SetMarkerSize(0)
    ratio.SetLineWidth(0)
    ratio.GetYaxis().SetTitle("Pred./Meas.")
    #ratio.GetXaxis().SetTitle(x_title)
    ratio.GetXaxis().SetTitleOffset(5)  
    ratio.GetXaxis().SetTitleSize(0.1)
    ratio.GetXaxis().SetTickLength(0.1)
    ratio.GetXaxis().SetTitleOffset(0)
    ratio.GetXaxis().SetLabelSize(0.1)
    ratio.GetXaxis().SetLabelOffset(0.04)
    ratio.GetYaxis().SetTitleSize(0.1)
    ratio.GetYaxis().SetLabelSize(0.1)

    if num_bins < 15:
        for i, label in enumerate(labels):
            ratio.GetXaxis().SetBinLabel(i+1, label)
    else:
        ratio.GetXaxis().SetLabelSize(0)

    fitline = ROOT.TF1("fitline","[0]", 0, ratio.GetXaxis().GetXmax())

    hnumer.Fit("fitline")
    fitresult = hnumer.GetFunction("fitline")
    fitresult.SetLineStyle(1)
    fitresult.SetLineColor(8)
    chiSquare = fitresult.GetChisquare()
    n_dof = fitresult.GetNDF()
    YChi = fitresult.GetParameter(0)
    ErrorChi = fitresult.GetParError(0)
    leg.AddEntry(fitresult,"y = "+"{:.3f}".format(YChi) + " +/- " + "{:.3f}".format(ErrorChi),'e')


    ratio.DrawCopy("e")

    prefit_ratio.Draw("same PE2")
    postfit_ratio.Draw("same PE2")
    line.Draw("same")
    fitresult.Draw("same")
    if num_bins < 15:
        myText(0.85, 0.02, 1, 0.15, x_title)
    else:   
        myText(0.7, 0.02, 1, 0.15, x_title)

        # Define the axis labels and the corresponding bin ranges
        labels = ["200-300", "300-425", "425-600", "600-900", "900-2000"]

        bin_ranges = [0.22, 0.46, 0.65, 0.78, 0.88]  # Bin ranges where labels change

        for i in range(len(labels)):
            xoffset =  bin_ranges[i]
            text = labels[i]
            myText(xoffset, 0.2, 1, 0.1, text)

    c_1.cd()

    text = r"y = Fit to Prefit/Meas."
    myText(0.64, 0.53, 1, 0.044 , text)

    myText(0.66, 0.47, 1, 0.042,  r'\chi^{2}' + ' / dof = '+"{:.1f}".format(chiSquare) + " / " +str(n_dof))


    c_1.Update()
    c_2.Update()
    c.Update()

    printAndSavePads(c, c_1, c_2, outfile,outfolder,alsoLog=True)
    return c


def canvasPulls(name,pars,pulls,outfile,outfolder,pois=False,rescale=False):
    # Check for SetAtlasStyle and call it
    if hasattr(ROOT, "SetAtlasStyle"):
        ROOT.SetAtlasStyle()
    else:
        print("SetAtlasStyle is not defined in AtlasStyle.C")
    
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
    else:
        ymin=0.
        ymax=0.
        for p in pulls:
            ymin=min(ymin,pulls[p][0]+pulls[p][1]*1.2)
            ymax=max(ymax,pulls[p][0]+pulls[p][2]*1.2)
        h.GetYaxis().SetRangeUser(ymin,ymax)
    if not pois:
        h.GetYaxis().SetTitle('pull')
    #h.GetYaxis().SetRangeUser(-1,1)
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
    g.SetMinimum(-1.0)
    g.SetMaximum(1.0)
    for i,t in enumerate(pars):
        g.SetPoint(i,i,pulls[t][0])
        g.SetPointError(i,0,0,abs(pulls[t][1]),abs(pulls[t][2]))
    g.Draw('sameP')
    h.Draw('same')
    h.Draw('axissame')
    printAndSave(c,outfile,outfolder)

def canvasCorrelation(name,h,outfile,outfolder,log=False):
    # Check for SetAtlasStyle and call it
    if hasattr(ROOT, "SetAtlasStyle"):
        ROOT.SetAtlasStyle()
    else:
        print("SetAtlasStyle is not defined in AtlasStyle.C")
    
    
    ROOT.gStyle.SetPalette(ROOT.kThermometer);
    ROOT.gStyle.SetPaintTextFormat("4.2f")
    c=ROOT.TCanvas(name,name,1024,1024 )
    c.SetRightMargin(0.15)
    c.SetTopMargin(0.15)
    c.SetLeftMargin(0.15)
    c.SetBottomMargin(0.15)
    c.cd()
    h.SetStats(False)
    if not log:
        h.GetZaxis().SetRangeUser(-1,1)
    if log:
        c.SetLogz()
    h.Draw('COLZTEXT')

    ATLASLabel(0.15, 0.86, "Internal")
    sampleText = r"pp \rightarrow l^{\pm} \nu"
    if "elmu" not in outfolder and "el" in outfolder:
        sampleText = r"pp \rightarrow e^{\pm} \nu"
    if "elmu" not in outfolder and "mu" in outfolder:
        sampleText = r"pp \rightarrow \mu^{\pm} \nu"
        
    print ("outfile: ", outfolder)
        
    myText(0.6, 0.86, 1, 0.042, sampleText )
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
        print(imp,impacts[imp])
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

def printAndSavePads(c,pad1,pad2,outfile,outfolder,alsoLog=False):
    outfile.cd()
    c.Write()
    c.SaveAs(outfolder+'/'+c.GetName()+'.pdf')
    c.SaveAs(outfolder+'/'+c.GetName()+'.png')
    if alsoLog:
        pad1.SetLogy()
        c.SaveAs(outfolder+'/'+c.GetName()+'_log.pdf')
        c.SaveAs(outfolder+'/'+c.GetName()+'_log.png')


def ATLASLabel(x, y, text=None, color=1):
    l = ROOT.TLatex()
    l.SetNDC()
    l.SetTextFont(72)
    l.SetTextColor(color)

    delx = 0.115 * 696 * ROOT.gPad.GetWh() / (472 * ROOT.gPad.GetWw())

    l.DrawLatex(x, y, "ATLAS")
    if text:
        p = ROOT.TLatex()
        p.SetNDC()
        p.SetTextFont(42)
        p.SetTextColor(color)
        p.DrawLatex(x + delx, y, text)

def ATLASLabelOld(x, y, Preliminary=False, color=1):
    l = ROOT.TLatex()
    l.SetNDC()
    l.SetTextFont(72)
    l.SetTextColor(color)
    l.DrawLatex(x, y, "ATLAS")
    if Preliminary:
        p = ROOT.TLatex()
        p.SetNDC()
        p.SetTextFont(42)
        p.SetTextColor(color)
        p.DrawLatex(x + 0.115, y, "Preliminary")

def ATLASVersion(version=None, x=0.88, y=0.975, color=1):
    if version:
        versionString = "Version {}".format(version)
        l = ROOT.TLatex()
        l.SetTextAlign(22)
        l.SetTextSize(0.04)
        l.SetNDC()
        l.SetTextFont(72)
        l.SetTextColor(color)
        l.DrawLatex(x, y, versionString)

def SetAtlasStyle():
    atlasStyle = AtlasStyle()
    ROOT.gROOT.SetStyle("ATLAS")
    ROOT.gROOT.ForceStyle()
    print("\nApplying ATLAS style settings...\n")

def myText(x,y,color,size, text):

    #Double_t tsize=0.05
    l = ROOT.TLatex()
    
    l.SetTextSize(size)
    l.SetNDC()
    l.SetTextColor(color)
    l.DrawLatex(x,y,text)



def AtlasStyle():
    atlasStyle = ROOT.TStyle("ATLAS", "Atlas style")

    # Use plain black on white colors
    icol = 0  # WHITE
    atlasStyle.SetFrameBorderMode(icol)
    atlasStyle.SetFrameFillColor(icol)
    atlasStyle.SetCanvasBorderMode(icol)
    atlasStyle.SetCanvasColor(icol)
    atlasStyle.SetPadBorderMode(icol)
    atlasStyle.SetPadColor(icol)
    atlasStyle.SetStatColor(icol)

    # Set the paper & margin sizes
    atlasStyle.SetPaperSize(20, 26)

    # Set margin sizes
    atlasStyle.SetPadTopMargin(0.05)
    atlasStyle.SetPadRightMargin(0.05)
    atlasStyle.SetPadBottomMargin(0.16)
    atlasStyle.SetPadLeftMargin(0.16)

    # Set title offsets (for axis label)
    atlasStyle.SetTitleXOffset(1.4)
    atlasStyle.SetTitleYOffset(1.4)

    # Use large fonts
    font = 42  # Helvetica
    tsize = 0.05
    atlasStyle.SetTextFont(font)
    atlasStyle.SetTextSize(tsize)
    atlasStyle.SetLabelFont(font, "x")
    atlasStyle.SetTitleFont(font, "x")
    atlasStyle.SetLabelFont(font, "y")
    atlasStyle.SetTitleFont(font, "y")
    atlasStyle.SetLabelFont(font, "z")
    atlasStyle.SetTitleFont(font, "z")
    atlasStyle.SetLabelSize(tsize, "x")
    atlasStyle.SetTitleSize(tsize, "x")
    atlasStyle.SetLabelSize(tsize, "y")
    atlasStyle.SetTitleSize(tsize, "y")
    atlasStyle.SetLabelSize(tsize, "z")
    atlasStyle.SetTitleSize(tsize, "z")

    # Use bold lines and markers
    atlasStyle.SetMarkerStyle(20)
    atlasStyle.SetMarkerSize(1.2)
    atlasStyle.SetHistLineWidth(2)
    atlasStyle.SetLineStyleString(2, "[12 12]")  # postscript dashes

    # Get rid of X error bars (as recommended in ATLAS figure guidelines)
    atlasStyle.SetErrorX(0.0001)
    # Get rid of error bar caps
    atlasStyle.SetEndErrorSize(0.)

    # Do not display any of the standard histogram decorations
    atlasStyle.SetOptTitle(0)
    atlasStyle.SetOptStat(0)
    atlasStyle.SetOptFit(0)

    # Put tick marks on top and RHS of plots
    atlasStyle.SetPadTickX(1)
    atlasStyle.SetPadTickY(1)

    return atlasStyle

def ATLAS_LABEL(x, y, color=1):
    l = ROOT.TLatex()
    l.SetNDC()
    l.SetTextFont(72)
    l.SetTextColor(color)
    l.DrawLatex(x, y, "ATLAS")

def myTGraphErrorsDivide(g1, g2):
    n1 = g1.GetN()
    n2 = g2.GetN()

    if n1 != n2:
        print("Vectors do not have the same number of entries!")
        return None

    g3 = ROOT.TGraphErrors()

    for i in range(n1):
        x1, y1 = ROOT.Double(0), ROOT.Double(0)
        g1.GetPoint(i, x1, y1)

        x2, y2 = ROOT.Double(0), ROOT.Double(0)
        g2.GetPoint(i, x2, y2)

        if x1 != x2:
            print("x1 and x2 do not match!")
            continue

        dx1 = g1.GetErrorX(i)
        dy1 = g1.GetErrorY(i) / y1 if y1 != 0 else 0
        dy2 = g2.GetErrorY(i) / y2 if y2 != 0 else 0

        g3.SetPoint(i, x1, y1 / y2 if y2 != 0 else 0)
        g3.SetPointError(i, dx1, math.sqrt(dy1 ** 2 + dy2 ** 2) * (y1 / y2) if y2 != 0 else 0)

    return g3

