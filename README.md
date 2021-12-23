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
```