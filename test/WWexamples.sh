./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m fit -c all -o example1 -f
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m fit -c all -p cW,cHDD,cHWB,cll1,cHl3 -l -o example2 -f
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m fit -c all -p cW,cHDD,cHWB,cll1,cHl3 -l -d ".*mcstat.*" -o example3 -f
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m pulls -c all -o example4 -f
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m pulls -c all -p cW,cHDD,cHWB,cll1,cHl3 -l -o example5 -f
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m scan -p cW -c all -o example6 -f
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m scan -p cW -c all -w -o example7 -f
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m scan -p cW -c all -d all -o example8 -f
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m scan -p cW -c all -r 5 -o example9 -f
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m scan -p cW -c all -r -10:15 -n 500 -o example10 -f
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m scan -p cW -c all -a -o example11 -f
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m scan -p cW,cll1 -c all -l -o example12 -f
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m profilelikelihood -p cW -c all -o example13 -f
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m workspace -l -o example14 -f
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m workspace -l --override measurements:Measurement.stat_unc_threshold=-1 -o example15 -f
