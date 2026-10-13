# 1. General Information
EFT Fitter for UNfolded measurements -- now also fitting reco measurements.

For help:
```
./bin/eftfun.py
```

Prerequisites: python2.7, ROOT6 with RooFit (ideally ROOT 6.24 or newer), optionally YODA (typically installed with [Rivet](https://rivet.hepforge.org/trac/wiki/GettingStarted))

Config templates with explanations and examples:
```
configs/TopCfgTemplate.cfg
configs/MeasurementCfgTemplate.cfg
```

[Tutorial](https://gitlab.cern.ch/eft-tools/eft-fun/blob/master/EFTfit.md) (might be outdated in parts)

Fitting Examples:
```
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m fit -c all -o example1
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m fit -c all -p cW,cHDD,cHWB,cll1,cHl3 -l -o example2
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m fit -c all -p cW,cHDD,cHWB,cll1,cHl3 -l -d ".*mcstat.*" -o example3
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m pulls -c all -o example4
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m pulls -c all -p cW,cHDD,cHWB,cll1,cHl3 -l -o example5
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m scan -p cW -c all -o example6
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m scan -p cW -c all -w -o example7
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m scan -p cW -c all -d all -o example8
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m scan -p cW -c all -r 5 -o example9
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m scan -p cW -c all -r -10:15 -n 500 -o example10
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m scan -p cW -c all -a -o example11
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m scan -p cW,cll1 -c all -l -o example12
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m profilelikelihood -p cW -c all -o example13
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m workspace -l -o example14
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m workspace -l --override measurements:Measurement.stat_unc_threshold=-1 -o example15
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m pca -l -o example16
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m scan -p cW -c all -n 10 -t 400 -o example17
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m scan -p cW -c all -n 10 -t 400 -o example18 --njobs 8
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m impacts -c all -p cW -o example19


# 2. Performing CCDY Fits

Setting up the environment:
```
setupATLAS
source ROOTstart.sh 
```
Running the fit, e.g:
```
source scripts/CCDY_u35_1d_22-11-24.sh
```
The running script should perform six fits: linear and linear+quadratic electron, muon, and combined.
These should be saved in e.g `/results/u35_22-11-24/`, under six different sub directories:
`mT_elmu` , `mT_el`, `mT_mu`, `mT_elmu_lin`, `mT_el_lin`, `mT_mu_lin`. 

The `mt_elmu` fit should contain a `--doBigFit` argument, in order to produce the cross-section vs Wilson coefficient plots. 

The fitpulls methods should have `pullsEachNP` enabled in order for the NP pulls to be produced when varying each Wilson coefficient one at a time. The prefit-postfit comparison and pulls on Wilson coefficients (`pullsEach`) should always be produced, along with the correlation matrix for the Wilson coefficients.

# 3. Git Workflow

Once you have made a change, perform the following commands to add the change and push it to the CCDY branch:
```
git status
git add .
git commit -m "Message describing change"
git push origin CCDY
```