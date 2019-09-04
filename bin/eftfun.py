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
import re,os,time
gROOT.SetBatch(True)

def main():    
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="config",
                      help="Read input config FILE. Mandatory.", metavar="FILE")
    parser.add_option("-o", "--outname", dest="outname",default=None,
                      help="Write output to OUTNAME.root and in folder OUTNAME/", metavar="OUTNAME")
    parser.add_option("-m", "--method",
                      dest="method",default="fit",
                      help="Method to be run, can be one of fit,scan,profilelikelihood,workspace,histo (fit - perform a maximum likelihood fit with MIGRAD+HESSE,  - perform a profile likelihood scan of one or multiple POIs, pulls - calculate nuisance parameter pulls and best fit values and errors with MINOS, profilelikelihood - run the RooStats ProfileLikelihoodCalculator to calculate confidece intervals, workspace - write out workspace and ModelConfig, histo - just write out histograms with given parameter values e.g. for validation of the parametrisation)"
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
                      help="Parameters that are set constant ('all' for all coefficients except nuisance parameters). Can be both nuisance parameters and POIs. Comma separated.")
    parser.add_option("-s", "--set", dest="initial", metavar="PARAMETERSETTING",
                      help="Set parameters to specific (initial) values with par1=x,par2=y,... Can be both nuisance parameters and POIs.")
    parser.add_option("-d", "--deactivate-uncertainties",metavar='PATTERN',
                      dest="deactivate_uncertainties", default=None,
                      help="Deactivate systematics (i.e. set them constant) that match the pattern given.")
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
    parser.add_option("-w", "--write_more",
                      action="store_true", dest="write_more", default=False,
                      help="Save all floating parameter values during scan in TTree? Will also plot non-linear parameter correlations.")
    parser.add_option("-f", "--force",
                      action="store_true", dest="force", default=False,
                      help="Run without user interaction?")
    parser.add_option("--override",
                      dest="override", default=None,
                      help="Override option in config. Syntax: Configname1:Section1.Option1=Value1,Configname2:Section2.Option2=Value2,... . Configname can also be 'measurements' or 'main' as shortcut for all measurments or the main cfg.")

    #TODO: "robust" option
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
        fit(options.config,options.outname,options.method,options.poi,options.scanrange,options.nsteps,options.ntoys,options.asimov,options.linear,options.deactivate_uncertainties,options.constant,options.initial,options.verbose,options.write_more,options.force,options.override)


def fit(config,outname,method,poi,scanrange,nsteps,ntoys,asimov,linear,deactivate_uncertainties,constant,initial,verbose,write_more,force,override):
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
        print 'write_more',write_more
        print 'override',override
        print '--------------------------------------'
    #TODO: make sure all options are used
    main_cfg_name=config
    SafeConfigParser.getAndRm = Utils.getAndRm
    SafeConfigParser.getbooleanAndRm = Utils.getbooleanAndRm
    SafeConfigParser.itemsAndRm = Utils.itemsAndRm
    main_cfg = SafeConfigParser()
    main_cfg.optionxform = str
    if override!=None:
        overrides=override.split(',')
    else:
        overrides=[]
    overrides_main=[]
    for ovr in overrides:
        if ovr.split(':')[0]=='main' or ovr.split(':')[0]==main_cfg_name:
            overrides_main+=[ovr.split(':')[1]]

    main_cfg_name=Utils.expandCfg(main_cfg_name,overrides_main)
    main_cfg.read(main_cfg_name)

    theory=Theory()
    if main_cfg.has_section('Ranges'):
        ranges = dict(x for x in main_cfg.itemsAndRm('Ranges'))
    else:
        ranges = {}
    for c in main_cfg.getAndRm('General','coeffs').split():
        if c in ranges:
            rng=[float(x) for x in ranges[c].split()]
        else:
            rng=[float(x) for x in main_cfg.getAndRm('General','default_range').split()]
        assert len(rng)==2
        theory.addCoeff(c,rng)
        
    measurements=[]
    predictions=[]
    for cfg_name in main_cfg.getAndRm('General','measurement_configs').split():
        cfg=SafeConfigParser()
        cfg.optionxform = str
        if not os.path.exists(cfg_name):
            cfg_name_orig=cfg_name
            cfg_name='/'.join(main_cfg_name.split('/')[:-1])+'/'+cfg_name
        if not os.path.exists(cfg_name):
            raise Exception("Could not find path "+cfg_name_orig+" or "+cfg_name)
        overrides_measurement=[]
        for ovr in overrides:
            if ovr.split(':')[0]=='measurements' or ovr.split(':')[0]==cfg_name:
                overrides_measurement+=[ovr.split(':')[1]]

        cfg.read(Utils.expandCfg(cfg_name,overrides_measurement))

        per_bin_width=False
        if cfg.has_option('Measurement','per_bin_width'):
            per_bin_width=cfg.getbooleanAndRm('Measurement','per_bin_width')
        binning=None
        if cfg.has_option('Measurement','binning'):
            binning=[float(x) for x in cfg.getAndRm('Measurement','binning').split()]


        measurement=Measurement(cfg.getAndRm('Measurement','name'),per_bin_width=per_bin_width,binning=binning)
        assert not measurement.name in [m.name for m in measurements]
        
        if asimov:
            measurement.setMeasuredFromCfg(cfg.getAndRm('Measurement','sm'))
        else:
            measurement.setMeasuredFromCfg(cfg.getAndRm('Measurement','measured'))

        if cfg.has_option('Measurement','covariance'):
            measurement.setCovFromCfg(cfg.getAndRm('Measurement','covariance'))
        elif cfg.has_option('Measurement','correlation') and cfg.has_option('Measurement','uncertainty'):
            measurement.setErrAndCorrFromCfg(cfg.getAndRm('Measurement','uncertainty'),cfg.getAndRm('Measurement','correlation'))
        else:
            raise Exception("Found neither covariance nor uncertainty and correlation for measurement in "+cfg_name)

        if cfg.has_section('Exp Systs'):
            for exp_syst_name,exp_syst in cfg.itemsAndRm('Exp Systs'):
                measurement.addUncertaintyFromCfg(exp_syst_name,exp_syst)
                
        prediction=Prediction(measurement.name,theory,per_bin_width=per_bin_width,binning=binning)
        assert prediction.name not in [p.name for p in predictions]
        
        prediction.setSMfromCfg(cfg.getAndRm('Measurement','sm'))
        if cfg.has_option('Measurement','renorm sm'):        
            prediction.setNormHistoFromCfg(cfg.getAndRm('Measurement','renorm sm'))
        if cfg.has_option('Measurement','scaling_factor'):
            prediction.setScalingFromCfg(cfg.getAndRm('Measurement','scaling_factor'))
        if cfg.has_option('Measurement','stat_unc_threshold'):
            prediction.stat_unc_threshold=float(cfg.getAndRm('Measurement','stat_unc_threshold'))
        if cfg.has_section('Theo Systs'):
            for theo_syst_name,theo_syst in cfg.itemsAndRm('Theo Systs'):
                prediction.addUncertaintyFromCfg(theo_syst_name,theo_syst)
        if cfg.has_section('Linear Effects'):
            for name,effect in cfg.itemsAndRm('Linear Effects'):
                if not prediction.setEffectFromCfg(1,name,effect):
                    print "WARNING: Could not find linear effect for "+name+" -- assuming it is zero"
        if not linear and cfg.has_section('Quadratic Effects'):
            names_set=[]
            for name,effect in cfg.items('Quadratic Effects'):
                if prediction.setEffectFromCfg(2,name,effect):
                    names_set+=[' '.join(name.split())]
            for name,effect in cfg.itemsAndRm('Quadratic Effects'):
                reverse=name.split()[1]+' '+name.split()[0]
                if not name in names_set and not reverse in names_set:
                    print "WARNING: Could not find quadratic effect for "+name+" or "+reverse+" -- assuming it is zero"

        for s in cfg.sections():
            if len(cfg.options(s))>0:
                raise Exception('Unused option '+','.join(cfg.options(s))+' in section '+s)
        
        measurements+=[measurement]
        predictions+=[prediction]
                    
    if verbose:
        print 'Ingredients ------------------------------------------------'
        for measurement,prediction in zip(measurements,predictions):
            print 'Measurement:',measurement.name
            print 'Measured:', ' '.join([str(x.getVal()) for x in measurement.central_values])
            print 'Uncertainty:', ' '.join(["{:.2f}% ".format(100*measurement.getCov()[i][i]**0.5/measurement.central_values[i].getVal()) for i in range(measurement.nbins)])

            print 'Covariance matrix:'
            measurement.getCov().Print()

            print 'Prediction formulas:'
            for i,s in enumerate(prediction.getPredictionStrings()):
                print
                print 'bin {})'.format(i)
                print s
            print
            print '------------------------------------------------'
    main_name=main_cfg.getAndRm('General','name')
    if outname==None:
        outname=main_name
    outname='results/'+outname
    if not os.path.exists(outname):
        os.makedirs(outname)
    f=TFile(outname+'.root','RECREATE')
    combination=Combination(theory,predictions,measurements,main_name,verbose=verbose) 

    if poi==None:
        pois=[]
    else:
        pois=poi.split(',')

    theory.setConstantFromCfg(constant,pois)
    theory.deactivateUncertainties(deactivate_uncertainties)
    if initial!=None:
        for x in initial.split(','):
            par=x.split("=")[0].strip()
            val=float(x.split("=")[1].strip())
            theory.setInitialVal(par,val)

    for s in main_cfg.sections():
        if len(main_cfg.options(s))>0:
            raise Exception('Unused option '+','.join(main_cfg.options(s))+' in section '+s)


    if not force:
        raw_input('Does the above look right (check before we are getting spammed by RooFit and Minuit)? Press Enter to continue.')

    start_time = time.time()
                    
    if method=='workspace':
        w=combination.createWorkspace()
        w.Write()
        f.Close()

    elif method=='histo':
        for p in predictions:
            p.plotPrefit(f)
        for m in measurements:
            m.plotPrefit(f)
        f.Close()

    elif method=='scan':
        extra_outfile=None
        if write_more:
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
            if write_more:
                for par in combination.all_results:
                    if par==poi:
                        continue
                    prof2=[(prof[i][0],y) for i,y in enumerate(combination.all_results[par])]
                    Plotting.graphProfile1D('g_scan_'+par,prof2,f,isNLL=False)
                    Plotting.canvasProfile1D('c_scan_'+par,prof2,f,poi,par,outfolder=outname,isNLL=False)
        else:
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

    elif method=='fit':
        if verbose:
            print 'doing fit'
        #for p in pois:
        #    theory.setFloating(p)
        data_histos_prefit=[]
        pred_histos_prefit=[]
        for p,m in zip(predictions,measurements):
            pred_histos_prefit+=[p.plotPrefit(f)]
            data_histos_prefit+=[m.plotPrefit(f)]
        combination.fit()
        corr=combination.getFitResult().correlationHist()
        corrPOIs=combination.correlationHistPOIs()
        corr.Write()
        corrPOIs.Write()
        Plotting.canvasCorrelation('c_correlation',corr,f,outname)
        Plotting.canvasCorrelation('c_correlationPOIs',corrPOIs,f,outname)
        
        for p,m,h_p_pre,h_m_pre in zip(predictions,measurements,pred_histos_prefit,data_histos_prefit):
            h_p_post=p.plotPostfit(f,combination.getFitResult())
            Plotting.canvasPostfit('c_postfit_'+m.name,h_m_pre,h_p_pre,h_p_post,outfile=f,outfolder=outname)

    elif method=='pulls':
        #for p in pois:
        #    theory.setFloating(p)        
        pulls=combination.pulls()
        pullsEach=combination.pullEach()
        print 'Results when floating all (MINOS errors)'
        for t in theory.coeffs+theory.thetas:
            if t.GetName() in pulls:
                print '{} = {}_{{ {} }}^{{ {} }}'.format(t.GetName(),*[ Utils.round2(p) for p in pulls[t.GetName()]])
        if len(pullsEach)>0:
            print 'Results when setting all but one poi to zero (MINOS errors)'
        for t in theory.coeffs:
            if t.GetName() in pullsEach:
                print '{} = {}_{{ {} }}^{{ {} }}'.format(t.GetName(),*[ Utils.round2(p) for p in pullsEach[t.GetName()]])

        Plotting.canvasPulls('c_pulls',theory.getThetaNames(),pulls,outfile=f,outfolder=outname)
        Plotting.canvasPulls('c_bestfit_all',theory.getCoeffNames(),pulls,outfile=f,outfolder=outname,pois=True,rescale=True)
        Plotting.canvasPulls('c_bestfit_each',theory.getCoeffNames(),pullsEach,outfile=f,outfolder=outname,pois=True,rescale=True)
            

    elif method=='profilelikelihood':
        combination.calcProfileLikelihood(pois)

    #TODO: implement hypotest inverter
        
    else:
        raise Excpetion("unkown method "+method)
    
    f.Close()

    print 'Method took',time.time() - start_time,'seconds'
    
if __name__== "__main__":
    main()



