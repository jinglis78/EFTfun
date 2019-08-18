#!/usr/bin/env python
import sys
sys.path.insert(0, '../python')
sys.path.insert(0, 'python')
from ConfigParser import SafeConfigParser
from Measurement import Measurement
from Theory import Theory
from Prediction import Prediction
from Combination import Combination
from ROOT import TFile,RooWorkspace,gROOT,TCanvas
import Utils,Plotting
from optparse import OptionParser
import re,os
gROOT.SetBatch(True)

#TODO: pulls, nicer plots

def main():    
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="config",
                      help="Read input config FILE. Mandatory.", metavar="FILE")
    parser.add_option("-o", "--outname", dest="outname",default='output',
                      help="Write output to OUTNAME.root and in folder OUTNAME/", metavar="OUTNAME")
    parser.add_option("-m", "--method",
                      dest="method",default="fit",
                      help="Method to be run, can be one of fit,scan,profilelikelihood,workspace,histo (fit - perform a maximum likelihood fit,  - perform a profile likelihood scan of one or multiple POIs, profilelikelihood - run the RooStats ProfileLikelihoodCalculator to calculate confidece intervals, workspace - write out workspace and ModelConfig, histo - just write out histograms with given parameter values e.g. for validation of the parametrisation)"
    )
    parser.add_option("-p", "--poi",dest='poi',default=None,
                      help="Parameters of interest, comma separated. In the case of a maximum likelihood fit, POIs will always be floating.")
    parser.add_option("-q", "--quiet",
                      action="store_false", dest="verbose", default=True,
                      help="Print less.")
    parser.add_option("-r", "--range",
                      dest="scanrange",default='3',metavar='SCANRANGE',
                      help="Scan range, either given as min:max or one number. If only one number, scan will be performed in range (-SCANRANGE*error+central,central+SCANRANGE*error). For multiple POIs, multiple scan ranges can be given, separated by commas. Only relevant for method 'scan'.")
    parser.add_option("-c", "--constant", dest="constant", metavar="CONSTANTPARS",
                      help="Parameters that are set constant ('all' for all coefficients except nuisance parameters). Can also be set to specific values with par=value. Can both be nuisance parameters and POIs. Comma separated.")
    parser.add_option("-d", "--deactivate-uncertainties",metavar='PATTERN',
                      dest="deactivate_uncertainties", default=None,
                      help="Deactivate systematics that match the pattern given.")
    parser.add_option("-n", "--nsteps",dest="nsteps",default='200',
                      help="Number of steps in scan (only relevant for method 'scan').")
    parser.add_option("-t", "--ntoys",dest="ntoys",default=1000,
                      help="Number of toys (currently not used).")
    parser.add_option("-a", "--asimov",
                      action="store_true", dest="asimov", default=False,
                      help="Run on asimov dataset (i.e. exchange data with sm prediction)?")
    parser.add_option("-l", "--linear",
                      action="store_true", dest="linear", default=False,
                      help="Consider only effects linear in EFT coefficients?")
    parser.add_option("-s", "--save",
                      action="store_true", dest="save", default=False,
                      help="Save all floating parameter values during scan in TTree?")
    parser.add_option("-f", "--force",
                      action="store_true", dest="force", default=False,
                      help="Run without user interaction?")

    options,args=parser.parse_args()

    if  options.method==None:
        print 'need method'
        parser.print_help()
    elif  options.config==None or not os.path.exists(options.config):
        print 'need config'
        parser.print_help()
    elif (options.method=='scan' or options.method=='profilelikelihood' ) and options.poi==None:
        print 'need poi'
        parser.print_help()
    else:
        fit(options.config,options.outname,options.method,options.poi,options.scanrange,options.nsteps,options.ntoys,options.asimov,options.linear,options.deactivate_uncertainties,options.constant,options.verbose,options.save,options.force)


def fit(config,outname,method,poi,scanrange,nsteps,ntoys,asimov,linear,deactivate_uncertainties,constant,verbose,save,force):
    if verbose:
        print '--------------------------------------'
        print 'running fit with following parameters:'
        print 'config',config
        print 'outname',outname
        print 'poi',poi
        print 'scanrange',scanrange
        print 'nsteps',nsteps
        print 'ntoys',ntoys
        print 'asimov',asimov
        print 'linear',linear
        print 'deactivate_uncertainties',deactivate_uncertainties
        print 'constant',constant
        print 'verbose',verbose
        print 'save',save
        print '--------------------------------------'
        
    main_cfg_name=config
    main_cfg = SafeConfigParser()
    main_cfg.optionxform = str
    main_cfg_name=Utils.expandCfg(main_cfg_name)
    main_cfg.read(main_cfg_name)

    theory=Theory()
    if main_cfg.has_section('Ranges'):
        ranges = dict(x for x in main_cfg.items('Ranges'))
    else:
        ranges = {}
    for c in main_cfg.get('General','coeffs').split():
        if c in ranges:
            rng=[float(x) for x in ranges[c].split()]
        else:
            rng=[float(x) for x in main_cfg.get('General','default_range').split()]
        assert len(rng)==2
        theory.addCoeff(c,rng)
        
    measurements=[]
    predictions=[]
    for cfg_name in main_cfg.get('General','measurement_configs').split():
        cfg=SafeConfigParser()
        cfg.optionxform = str
        if not os.path.exists(cfg_name):
            cfg_name_orig=cfg_name
            cfg_name='/'.join(main_cfg_name.split('/')[:-1])+'/'+cfg_name
        if not os.path.exists(cfg_name):
            raise Exception("Could not find path "+cfg_name_orig+" or "+cfg_name)
        cfg.read([main_cfg_name,Utils.expandCfg(cfg_name)])

        per_bin_width=False
        if cfg.has_option('Measurement','per_bin_width'):
            per_bin_width=cfg.getboolean('Measurement','per_bin_width')
        binning=None
        if cfg.has_option('Measurement','binning'):
            binning=[float(x) for x in cfg.get('Measurement','binning').split()]


        measurement=Measurement(cfg.get('Measurement','name'),per_bin_width=per_bin_width,binning=binning)
        assert not measurement.name in [m.name for m in measurements]
        
        if asimov:
            measurement.setMeasuredFromCfg(cfg.get('Measurement','sm'))
        else:
            measurement.setMeasuredFromCfg(cfg.get('Measurement','measured'))

        if cfg.has_option('Measurement','covariance'):
            measurement.setCovFromCfg(cfg.get('Measurement','covariance'))
        elif cfg.has_option('Measurement','correlation') and cfg.has_option('Measurement','uncertainty'):
            measurement.setErrAndCorrFromCfg(cfg.get('Measurement','uncertainty'),cfg.get('Measurement','correlation'))
        else:
            raise Exception("Found neither covariance nor uncertainty and correlation for measurement in "+cfg_name)
        prediction=Prediction(cfg.get('Measurement','name'),theory,per_bin_width=per_bin_width,binning=binning)
        assert prediction.name not in [p.name for p in predictions]
        
        prediction.setSMfromCfg(cfg.get('Measurement','sm'))
        if cfg.has_option('Measurement','renorm sm'):        
            prediction.setNormHistoFromCfg(cfg.get('Measurement','renorm sm'))
        if cfg.has_option('Measurement','scaling_factor'):
            prediction.setScalingFromCfg(cfg.get('Measurement','scaling_factor'))
        if cfg.has_option('Measurement','stat_unc_threshold'):
            prediction.stat_unc_threshold=float(cfg.get('Measurement','stat_unc_threshold'))
        if cfg.has_section('Theo Systs'):
            for theo_syst_name,theo_syst in cfg.items('Theo Systs'):
                prediction.addUncertaintyFromCfg(theo_syst_name,theo_syst)
        if cfg.has_section('Linear Effects'):
            for name,effect in cfg.items('Linear Effects'):
                if not prediction.setEffectFromCfg(1,name,effect):
                    print "WARNING: Could not find linear effect for "+name+" -- assuming it is zero"
        if not linear and cfg.has_section('Quadratic Effects'):
            names_set=[]
            for name,effect in cfg.items('Quadratic Effects'):
                if prediction.setEffectFromCfg(2,name,effect):
                    names_set+=[' '.join(name.split())]
            for name,effect in cfg.items('Quadratic Effects'):
                reverse=name.split()[1]+' '+name.split()[0]
                if not name in names_set and not reverse in names_set:
                    print "WARNING: Could not find quadratic effect for "+name+" or "+reverse+" -- assuming it is zero"
                        
        measurements+=[measurement]
        predictions+=[prediction]
                    
    if verbose:
        print 'Ingredients ------------------------------------------------'
        for measurement,prediction in zip(measurements,predictions):
            print 'Measurement:',measurement.name
            print 'Measured:', ' '.join([str(x.getVal()) for x in measurement.central_values])

            print 'Covariance matrix:'
            measurement.covMatrix.Print()

            print 'Prediction formulas:'
            for s in prediction.getPredictionStrings():
                print 
                print s
            print
            print '------------------------------------------------'

    if not force:
        raw_input('Does the above look right (check before we are getting spammed by RooFit and Minuit)? Press Enter to continue.')

    if not os.path.exists(outname):
        os.mkdir(outname)
    f=TFile(outname+'.root','RECREATE')
    combination=Combination(theory,predictions,measurements,main_cfg.get('General','name'),verbose=verbose) 

    if poi==None:
        pois=[]
    else:
        pois=poi.split(',')

    theory.deactivateUncertainties(deactivate_uncertainties)
    theory.setConstantFromCfg(constant,pois)
                    
    if method=='workspace':
        w=combination.createWorkspace()
        w.Write()
        f.Close()

    if method=='histo':
        for p in predictions:
            p.plotPrefit(f)
        for m in measurements:
            m.plotPrefit(f)
        f.Close()

    if method=='scan':
        extra_outfile=None
        if save:
            extra_outfile=f
        if len(pois)==1:
            poi=pois[0]
            if ':' in str(scanrange):
                actual_scanrange=[float(x) for x in scanrange.split(':')]
            else:        
                actual_scanrange=float(scanrange)
            prof=combination.profile1D(poi,actual_scanrange,int(nsteps),extra_outfile=extra_outfile)
            Plotting.graphProfile1D('g_scan_'+poi,prof,f)
            Plotting.canvasProfile1D('c_scan_'+poi,prof,f,poi,outfolder=outname)
        else:
            #TODO
            scanranges=scanrange.split(',')
            actual_scanranges=[]
            for s in scanranges:
                if ':' in s:
                    actual_scanranges.append([float(x) for x in s.split(':')])
                else:
                    actual_scanranges.append(float(s))
            nsteps=[int(n) for n in nsteps.split(',')]
            prof=combination.profileND(pois,actual_scanranges,nsteps,extra_outfile=extra_outfile)
            if len(pois)==2:
                Plotting.graphProfile2D('scan_'+'_'.join(pois),prof,f)
                Plotting.canvasProfile2D('scan_'+'_'.join(pois),prof,f,pois[0],pois[1],outfolder=outname)

    if method=='fit':
        if verbose:
            print 'doing fit'
        for p in pois:
            theory.setFloating(p)
        data_histos_prefit=[]
        pred_histos_prefit=[]
        for p,m in zip(predictions,measurements):
            pred_histos_prefit+=[p.plotPrefit(f)]
            data_histos_prefit+=[m.plotPrefit(f)]
        combination.minimize()
        corr=combination.getFitResult().correlationHist()
        corrPOIs=combination.correlationHistPOIs()
        corr.Write()
        corrPOIs.Write()
        Plotting.canvasCorrelation('c_correlation',corr,f,outname)
        Plotting.canvasCorrelation('c_correlationPOIs',corrPOIs,f,outname)
        
        for p,m,h_p_pre,h_m_pre in zip(predictions,measurements,pred_histos_prefit,data_histos_prefit):
            h_p_post=p.plotPostfit(f,combination.getFitResult())
            Plotting.canvasPostfit('c_postfit_'+m.name,h_m_pre,h_p_pre,h_p_post,outfile=f,outfolder=outname)

    if method=='pulls':
        pulls=combination.pulls()
        theta_names=[]
        for t in theory.thetas:
            if t.GetName() in pulls:
                print t.GetName(),pulls[t.GetName()]
                theta_names+=[t.GetName()]
        Plotting.canvasPulls('c_pulls',theta_names,pulls,outfile=f,outfolder=outname)
            

    if method=='profilelikelihood':
        combination.calcProfileLikelihood(pois)

    f.Close()
if __name__== "__main__":
    main()



