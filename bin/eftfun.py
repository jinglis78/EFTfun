#!/usr/bin/env python
import sys
sys.path.insert(0, '../python')
sys.path.insert(0, 'python')
import MyConfigParser
from Measurement import Measurement
from Theory import Theory
from Prediction import Prediction
from Combination import Combination
from ROOT import TFile,RooWorkspace,gROOT,TCanvas
import Utils,Plotting
from optparse import OptionParser
import re,os,time
from subprocess import call
gROOT.SetBatch(True)

def main():    
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="config",
                      help="Read input config FILE. Mandatory.", metavar="FILE")
    parser.add_option("-o", "--outname", dest="outname",default=None,
                      help="Write output to OUTNAME.root and in folder OUTNAME/", metavar="OUTNAME")
    parser.add_option("-m", "--method",
                      dest="method",default="fit",
                      help="Method to be run, can be one of fit,scan,profilelikelihood,workspace,histo (fit - perform a maximum likelihood fit with MIGRAD+HESSE,  - perform a profile likelihood scan of one or multiple POIs, pulls - calculate nuisance parameter pulls and best fit values and errors with MINOS, profilelikelihood - run the RooStats ProfileLikelihoodCalculator to calculate confidece intervals, pca - run principal component analysis for linear fit, workspace - write out workspace and ModelConfig, histo - just write out histograms with given parameter values e.g. for validation of the parametrisation)"
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
    parser.add_option("--quad",
                      action="store_true", dest="quad", default=False,
                      help="Consider only effects quadratic in EFT coefficients?")
    parser.add_option("-w", "--write_more",
                      action="store_true", dest="write_more", default=False,
                      help="Save all floating parameter values during scan in TTree? Will also plot non-linear parameter correlations.")
    parser.add_option("--cls",
                      action="store_true", dest="cls", default=False,
                      help="Do CLs limit (only for frequentist/hybrid/asymptotic calculator)/")
    parser.add_option("-f", "--force",
                      action="store_true", dest="force", default=False,
                      help="Run without user interaction?")
    parser.add_option("--override",
                      dest="override", default=None,
                      help="Override option in config. Syntax: Configname1:Section1.Option1=Value1,Configname2:Section2.Option2=Value2,... . Configname can also be 'measurements' or 'main' as shortcut for all measurments or the main cfg.")
    parser.add_option("--measurements",
                      dest="measurements", default='all',
                      help="Comma separated list of measurements to be included (default: include all measurements)")
    parser.add_option("--repara",
                      action="store_true", dest="repara", default=False,
                      help="Use reparametrization given in config?")
    parser.add_option("-g", "--groups",dest="groups",default=None,
                      help="Parameter groups for PCA, e.g. c1,c2:c2,c3,c4:c5.")
    parser.add_option("--robust",
                      dest="robust", default=1,
                      help="Robustness. Best to be kept at 1 or 2. Should the fit be complicated or fail too often use 3. 0 for faster fits.")
    parser.add_option("--offset",
                      action="store_true", dest="offset", default=False,
                      help="Don't subtract min NLL")
    parser.add_option("--implicitscan",
                      dest="implicitscan", default=None,
                      help="Profile variable manually in case it is too complicated for minuit (e.g. due to local minima), syntax: paramter:steps:min:max")
    parser.add_option("--conts1d",
                      action="store_true", dest="conts1d", default=False,
                      help="Use contours corresponding to 1D chi2 also 2D fits")


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
        fit(options.config,options.outname,options.method,options.poi,options.scanrange,options.nsteps,options.ntoys,options.asimov,options.linear,options.quad,options.deactivate_uncertainties,options.constant,options.initial,options.verbose,options.write_more,options.cls,options.force,options.override,options.measurements,options.repara,options.groups,options.robust,options.offset,options.implicitscan,options.conts1d)


def fit(config,outname,method,poi,scanrange,nsteps,ntoys,asimov,linear,quad,deactivate_uncertainties,constant,initial,verbose,write_more,cls,force,override,measurements_included,repara,groups,robust,offset,implicitscan,conts1d):
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
        print 'quad',quad
        print 'deactivate_uncertainties',deactivate_uncertainties
        print 'constant',constant
        print 'initial',initial
        print 'verbose',verbose
        print 'write_more',write_more
        print 'cls',cls
        print 'override',override
        print 'measurements',measurements_included
        print 'repara',repara
        print 'groups',groups
        print 'robust',robust
        print 'offset',offset
        print 'implicitscan',implicitscan
        print 'conts1d',conts1d
        print '--------------------------------------'

    main_cfg_name=config
    main_cfg = MyConfigParser.MyConfigParser()
    main_cfg.optionxform = str

    if override!=None:
        overrides=override.split(',')
    else:
        overrides=[]
    overrides_main=[]
    for ovr in overrides:
        if ovr.split(':')[0]=='main' or ovr.split(':')[0]==main_cfg_name:
            overrides_main+=[':'.join(ovr.split(':')[1:])]
   
    main_cfg_name=MyConfigParser.expandCfg(main_cfg_name,overrides_main)
    main_cfg.read(main_cfg_name)

    theory=Theory()
    if main_cfg.has_section('Ranges'):
        ranges = dict(x for x in main_cfg.itemsAndRm('Ranges'))
    else:
        ranges = {}
    default_range=None
    if main_cfg.has_option('General','default_range'):
        default_range=main_cfg.getAndRm('General','default_range')
    if main_cfg.has_option('General','coeffs'):
        for c in main_cfg.getAndRm('General','coeffs').split():
            if c in ranges:
                rng=[float(x) for x in ranges[c].split()]
            elif default_range!=None:
                rng=[float(x) for x in default_range.split()]
            else:
                raise Exception("No range given for coefficient "+c+" and no default range set")
            assert len(rng)==2
            theory.addCoeff(c,rng)
    correlated_measurements=[]
    if main_cfg.has_option('General','correlated_measurements'):
        correlated_measurements=main_cfg.getAndRm('General','correlated_measurements').split()
    c_mapping={}
    cs_new=[]
    if main_cfg.has_section('Reparameterization'):
        for cnew,para in main_cfg.itemsAndRm('Reparameterization'):
            cs_new.append(cnew)
            c_mapping[cnew]={}
            for cold,x in zip(theory.getCoeffNames(),[float(x) for x in para.split()]):
                c_mapping[cnew][cold]=x
    measurements=[]
    predictions=[]
    for cfg_name in main_cfg.getAndRm('General','measurement_configs').split():
        cfg=MyConfigParser.MyConfigParser()
        cfg.optionxform = str
        if not os.path.exists(cfg_name):
            cfg_name_orig=cfg_name
            cfg_name='/'.join(main_cfg_name.split('/')[:-1])+'/'+cfg_name
        if not os.path.exists(cfg_name):
            raise Exception("Could not find path "+cfg_name_orig+" or "+cfg_name)
        overrides_measurement=[]
        for ovr in overrides:
            if ovr.split(':')[0]=='measurements' or ovr.split(':')[0]==cfg_name:
                overrides_measurement+=[':'.join(ovr.split(':')[1:])]
        cfg.read(MyConfigParser.expandCfg(cfg_name,overrides_measurement))        
        measurement_name=cfg.getAndRm('Measurement','name')
        if not measurements_included=='all' and not measurement_name in [x.strip() for x in measurements_included.split(',')]:
            continue
        reco_level=False
        if cfg.has_option('Measurement','reco_level'):
            reco_level=cfg.getbooleanAndRm('Measurement','reco_level')
        per_bin_width=False
        if cfg.has_option('Measurement','per_bin_width'):
            per_bin_width=cfg.getbooleanAndRm('Measurement','per_bin_width')
        binning=None
        if cfg.has_option('Measurement','binning'):
            binning=[float(x) for x in cfg.getAndRm('Measurement','binning').split()]
        binlabels=None
        if cfg.has_option('Measurement','binlabels'):
            binlabels=cfg.getAndRm('Measurement','binlabels').split()
        blind=None
        if cfg.has_option('Measurement','blind'):
            blind=[int(x) for x in cfg.getAndRm('Measurement','blind').split()]
        if verbose:
            print "creating prediction for", measurement_name
        prediction=Prediction(measurement_name,theory,per_bin_width=per_bin_width,binning=binning,blind=blind,binlabels=binlabels)
        assert prediction.name not in [p.name for p in predictions]
        
        prediction.setSMfromCfg(cfg.getAndRm('Measurement','sm'))
        if cfg.has_option('Measurement','renorm_sm'):        
            prediction.setNormHistoFromCfg(cfg.getAndRm('Measurement','renorm_sm'))
        if cfg.has_option('Measurement','scaling_factor'):
            prediction.setScalingFromCfg(cfg.getAndRm('Measurement','scaling_factor'))
        if cfg.has_option('Measurement','stat_unc_threshold'):
            prediction.stat_unc_threshold=float(cfg.getAndRm('Measurement','stat_unc_threshold'))
        rebin=None
        if cfg.has_option('Measurement','rebin'):
            rebin=cfg.getAndRm('Measurement','rebin')

        #linear_corrections={}
        #if cfg.has_section('Linear Corrections'):
        #    for name,effect in cfg.itemsAndRm('Linear Corrections'):
        #        if name not in theory.getCoeffNames():
        #            continue
        #        else:                    
        #            linear_corrections[name]=float(effect)
            
        if cfg.has_section('Linear Effects'):
            if quad:
                cfg.remove_section('Linear Effects')
            else:
                for name,effect in cfg.itemsAndRm('Linear Effects'):
                    if name not in theory.getCoeffNames():
                        continue
                    if not prediction.setEffectFromCfg(1,name,effect):
                        print "WARNING: Could not find linear effect for "+name+" and measurement "+measurement_name+" -- assuming it is zero"

        if cfg.has_section('Linear Corrections'):
            if quad:
                cfg.remove_section('Linear Corrections')
            else:
                for name,effect in cfg.itemsAndRm('Linear Corrections'):
                    if name not in theory.getCoeffNames():
                        continue
                    if not prediction.setEffectFromCfg(1,name,effect,add=True):
                        print "WARNING: Could not find linear correction for "+name+" and measurement "+measurement_name+" -- assuming it is zero"

        if cfg.has_section('Quadratic Effects'):
            if linear:
                cfg.remove_section('Quadratic Effects')
            else:
                names_set=[]
                for name,effect in cfg.items('Quadratic Effects'):
                    if not len(name.split())==2:
                        raise Exception('Quadratic effects need to be denoted with two coefficients (e.g. c1 c1 or c1 c2)and not like this: '+name)
                    if name.split()[0] not in theory.getCoeffNames() or  name.split()[1] not in theory.getCoeffNames():
                        continue                    
                    if prediction.setEffectFromCfg(2,name,effect):
                        names_set+=[' '.join(name.split())]
                for name,effect in cfg.itemsAndRm('Quadratic Effects'):
                    if name.split()[0] not in theory.getCoeffNames() or  name.split()[1] not in theory.getCoeffNames():
                        continue
                    reverse=name.split()[1]+' '+name.split()[0]
                    if not name in names_set and not reverse in names_set:
                        print "WARNING: Could not find quadratic effect for "+name+" or "+reverse+" and measurement "+measurement_name+" -- assuming it is zero"
        if cfg.has_section('Prediction Uncertainties'):
            for theo_syst_name,theo_syst in cfg.itemsAndRm('Prediction Uncertainties'):
                prediction.addUncertaintyFromCfg(theo_syst_name,theo_syst,reco_level=reco_level)

        if cfg.has_section('Backgrounds'):
            for name,data in cfg.itemsAndRm('Backgrounds'):
                assert reco_level
                prediction.addBackgroundFromCfg(name,data)
            if cfg.has_section('Background Uncertainties'):
                for bkg_syst_name,bkg_syst in cfg.itemsAndRm('Background Uncertainties'):
                    assert reco_level
                    prediction.addBkgUncertaintyFromCfg(bkg_syst_name,bkg_syst)
                    
        if verbose:
            print "creating measurement for", measurement_name
        measurement=Measurement(measurement_name,per_bin_width=per_bin_width,binning=binning,reco_level=reco_level,theory=theory,blind=blind,binlabels=binlabels)
        assert not measurement.name in [m.name for m in measurements]
        if asimov:
            measurement.setMeasuredFromCfg(cfg.getAndRm('Measurement','measured'))
            asmimov_data=prediction.getNominal()
            for b in prediction.backgrounds.values():
                asmimov_data=[x+y for x,y in zip(asmimov_data,b.getNominal())]
            measurement.setMeasured(asmimov_data)
            
        else:
            measurement.setMeasuredFromCfg(cfg.getAndRm('Measurement','measured'))
        if cfg.has_option('Measurement','covariance'):
            measurement.setCovFromCfg(cfg.getAndRm('Measurement','covariance'))
        elif cfg.has_option('Measurement','correlation') and cfg.has_option('Measurement','correlated_uncertainty'):
            measurement.setErrAndCorrFromCfg(cfg.getAndRm('Measurement','correlated_uncertainty'),cfg.getAndRm('Measurement','correlation'))
        elif not measurement.reco_level:
            print 'WARNING: neither covariance nor correlation+uncertainty in config, covariance will be built from sources only, will likely mean uncorrelated statistical uncertainties'
        if cfg.has_section('Measurement Uncertainties'):
            if verbose:
                print "adding uncertainties for",measurement_name
            for exp_syst_name,exp_syst in cfg.itemsAndRm('Measurement Uncertainties'):
                measurement.addUncertaintyFromCfg(exp_syst_name,exp_syst)
                
        for s in cfg.sections():
            if len(cfg.options(s))>0:
                raise Exception('Unused option '+','.join(cfg.options(s))+' in section '+s)

        if not rebin is None:
            if len(rebin.split())==1:
                combinebins=int(rebin)
            else:
                combinebins=[int(x) for x in rebin.split()]
            measurement.rebin(combinebins)
            prediction.rebin(combinebins)

        if method=='pca':
            print 'Modify model for PCA: add extra scaling factors'
            scaling=[]
            for i in range(prediction.nbins):
                s='mu_{}_{}'.format(prediction.name,i)
                scaling.append(s)
                theory.addCoeff(s,(-100,100))
            prediction.setScaling(scaling,scaling)

            
        measurements+=[measurement]
        predictions+=[prediction]
    if repara and c_mapping!=None:
        theory.rePara(cs_new,c_mapping)
        for p in predictions:
            p.rePara(cs_new,c_mapping)
    if verbose:
        print 'Ingredients ------------------------------------------------'
        for measurement,prediction in zip(measurements,predictions):
            print 'Measurement:',measurement.name
            print
            print 'Measured:', ' '.join([str(x.getVal()) for x in measurement.central_values])
            print
            if measurement.binlabels and isinstance(measurement.binlabels[0],str):
                print "Labels: "+' '.join(measurement.binlabels)
                print
            if measurement.getMeasurementCorrectionString()!=None:
                print 'Measurement correction formulas:'
                for i,s in enumerate(measurement.getMeasurementCorrectionString()):
                    print
                    print 'bin {}:'.format(i)
                    print s.replace(')+',')\n+').replace('))','))\n')
            if not measurement.reco_level:
                print 'Uncertainty (excluding nuisance parameters):', ' '.join(["{:.2f}% ".format(100.*x) for x in measurement.getRelErr()])
                print
                print 'Correlation matrix:'
                measurement.getCorr().Print()

            print 'Prediction formulas:'
            for i,s in enumerate(prediction.getPredictionStrings()):
                print
                print 'bin {}:'.format(i)
                print s.replace(')+',')\n+').replace('))','))\n')
            print
            print '------------------------------------------------'
    main_name=main_cfg.getAndRm('General','name')
    if outname==None:
        outname=main_name
    outname='results/'+outname
    if not os.path.exists(outname):
        os.makedirs(outname)
    if os.path.exists(outname+'.root'): # and (force or Utils.ask("Outfile {0}.root already exists, update (u) or recreate (r)? You can also ctrl+c and set a different name with the -o flag.".format(outname),['u'],['r'])):        
        f=TFile(outname+'.root','UPDATE')
        print 'WARNING: updating',outname+'.root'
    else:
        f=TFile(outname+'.root','RECREATE')

    combination=Combination(theory,predictions,measurements,main_name,correlated_measurements=correlated_measurements,verbose=verbose,robust=robust) 

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
        if len(pois)==0:
            pois='all'
        w=combination.createWorkspace(pois)
        w.Write()

    elif method=='histo':
        for p in predictions:
            p.plotPrefit(f)
            p.plotPrefit(f,statonly=True)
            p.plotScaling2D(f)
        for m in measurements:
            m.plotPrefit(f)

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
            if implicitscan is None:
                prof=combination.profile1D(poi,actual_scanrange,int(nsteps),extra_outfile=extra_outfile,inclOffset=offset)
            else:
                profs=[]
                impl_par,impl_n,impl_min,impl_max=implicitscan.split(':')
                print 'repeating scan with',impl_par,'set from',impl_min,'to',impl_max,'in',impl_n,'steps'
                for impl_x in [float(impl_min)+float(x)*(float(impl_max)-float(impl_min))/(int(impl_n)-1) for x in range(int(impl_n))]:
                    theory.setConstant(impl_par,impl_x)
                    theory.setInitialVal(impl_par,impl_x)
                    profs.append(combination.profile1D(poi,actual_scanrange,int(nsteps),extra_outfile=extra_outfile,inclOffset=True))
                prof=[]
                total_min=1e9
                for pointresults in map(list, zip(*profs)):
                    best=1e9
                    for x in pointresults:
                        best=min(x[1],best)
                    total_min=min(best,total_min)
                    prof.append((pointresults[0][0],best))
                prof=[(x[0],x[1]-total_min) for x in prof]
            Plotting.graphProfile1D('g_scan_'+poi,prof,f)
            Plotting.canvasProfile1D('c_scan_'+poi,prof,f,poi,outfolder=outname)
            if write_more:
                for par in combination.all_results:
                    if par==poi:
                        continue
                    prof2=[(prof[i][0],y) for i,y in enumerate(combination.all_results[par])]
                    Plotting.graphProfile1D('g_scan_'+'_'.join([poi,par]),prof2,f,isNLL=False)
                    Plotting.canvasProfile1D('c_scan_'+'_'.join([poi,par]),prof2,f,poi,par,outfolder=outname,isNLL=False)
        else:
            scanranges=scanrange.split(',')
            actual_scanranges=[]
            for s in scanranges:
                if ':' in s:
                    actual_scanranges.append([float(x) for x in s.split(':')])
                else:
                    actual_scanranges.append(float(s))
            nsteps=[int(n) for n in nsteps.split(',')]
            if implicitscan is None:
                prof=combination.profileND(pois,actual_scanranges,nsteps,extra_outfile=extra_outfile,inclOffset=offset)
            else:
                profs=[]
                impl_par,impl_n,impl_min,impl_max=implicitscan.split(':')
                print 'repeating scan with',impl_par,'set from',impl_min,'to',impl_max,'in',impl_n,'steps'
                for impl_x in [float(impl_min)+float(x)*(float(impl_max)-float(impl_min))/(int(impl_n)-1) for x in range(int(impl_n))]:
                    theory.setConstant(impl_par,impl_x)
                    theory.setInitialVal(impl_par,impl_x)
                    profs.append(combination.profileND(pois,actual_scanranges,nsteps,extra_outfile=extra_outfile,inclOffset=True))
                prof=[]
                total_min=1e9
                for pointresults in map(list, zip(*profs)):
                    best=1e9
                    for x in pointresults:
                        best=min(x[-1],best)
                    total_min=min(best,total_min)
                    prof.append(pointresults[0][:-1]+(best,))
                prof=[x[:-1]+((x[-1]-total_min),) for x in prof]
            if len(pois)==2:
                Plotting.graphProfile2D('g_scan_'+'_'.join(pois),prof,f)
                Plotting.canvasProfile2D('c_scan_'+'_'.join(pois),prof,f,pois[0],pois[1],outfolder=outname,conts1d=conts1d)
    elif method=='fit':
        if verbose:
            print 'doing fit'
        data_histos_prefit=[]
        pred_histos_prefit=[]
        for p,m in zip(predictions,measurements):
            pred_histos_prefit+=[p.plotPrefit(f)]
            data_histos_prefit+=[m.plotPrefit(f)]
        combination.fit()
        corr=combination.getFitResult().correlationHist()        
        corr.Write()
        Plotting.canvasCorrelation('c_correlation',corr,f,outname)
        corrPOIs,covPOIs=combination.correlationHistPOIs()
        if corrPOIs:
            corrPOIs.Write()
            covPOIs.Write()
            Plotting.canvasCorrelation('c_correlationPOIs',corrPOIs,f,outname)
            Plotting.canvasCorrelation('c_covariancePOIs',covPOIs,f,outname,log=True)
        
        for p,m,h_p_pre,h_m_pre in zip(predictions,measurements,pred_histos_prefit,data_histos_prefit):
            h_p_post=p.plotPostfit(f,combination.getFitResult())
            hs_bkg_post=p.plotBkgPostfit(f,combination.getFitResult())
            h_m_post=None
            if m.getMeasurementCorrectionString()!=None:
                h_m_post=m.plotPostfit(f,combination.getFitResult())
            Plotting.canvasPostfit('c_postfit_'+m.name,h_m_pre,h_m_post,h_p_pre,h_p_post,outfile=f,outfolder=outname)

    elif method=='pulls':
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
        limit=combination.calcProfileLikelihood(pois)
        ftxt=open(outname+'/'+'interval_pl_'+pois[0]+'.txt','w')
        ftxt.write(','.join([str(x) for x in limit]))
        ftxt.close()

    elif method=='pca':
        if groups==None:
            grouping=[]
        else:
            grouping=[g.split(',') for g in groups.split(":")]
        h_cov,h_corr,h_ev,h_ev_abs,hs_effects=combination.PCA(grouping)
        for h in h_cov,h_corr,h_ev,h_ev_abs:
            h.Write()
        for ev in hs_effects:
            for h in ev:
                h.Write()

    elif method=='asymptotic' or method=='hybrid' or method=='frequentist':
        raw_input('WARNING, the method '+method+' is experimental and work in progress, press enter to proceed')
        if len(pois)!=1:
            print "Need exactly one POI to run asymptotic/hybrid/frequentist limits"
        w=combination.createWorkspace(pois)
        w.Write()
        f.Close()

        fn=outname+'.root'
        ws='w'
        model=combination.name+'_model'
        model_bkg=combination.name+'_model_background'
        data=combination.name+'_data'
        if method=='asymptotic':
            calculator=2
        if method=='hybrid':
            calculator=1
        if method=='frequentist':
            calculator=0
        teststat=2
        docls='false'
        if cls:
            docls='true'
        if not ':' in str(scanrange) or len(scanrange.split(':'))!=2:
            print 'WARNING, need to provide a scan range with "-r xmin:max"'
            return
        xmin=scanrange.split(':')[0]
        xmax=scanrange.split(':')[1]
        rootcmd='root -l scripts/StandardHypoTestInvDemo.C(\"{0}\",\"{1}\",\"{2}\",\"{3}\",\"{4}\",{5},{6},{7},{8},{9},{10},{11})'.format(
            fn,ws,model,model_bkg,data,calculator,teststat,docls,nsteps,xmin,xmax,ntoys)
        print 'Running:',rootcmd
        print 
        call(rootcmd.split())
        
    else:
        raise Exception("unkown method "+method)
    
    f.Close()

    print 'Method took',time.time() - start_time,'seconds'
    
if __name__== "__main__":
    main()



