EFT Fitter for UNfolded measurements.

For help:
```
./bin/eftfun.py
```

Prerequisites: python2.7, ROOT6 with RooFit, YODA (typically installed with Rivet)

On lxplus:
```
setupATLAS
asetup 21.6.6,AthGeneration
source setupRivet.sh
lsetup "root 6.14.04-x86_64-slc6-gcc62-opt"
```


Examples:
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
```