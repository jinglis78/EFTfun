./bin/eftfun.py -i configs/WW_WZ_combi/combi_EWd6.cfg -c all -p cWWWL2 -m scan -f -r -3.5:3.5  -s cWWWL2=1 -o combi_all_EWdim6
./bin/eftfun.py -i configs/WW_WZ_combi/combi_EWd6.cfg -c all -p cWL2 -m scan -f -r -8:5   -o combi_all_EWdim6
./bin/eftfun.py -i configs/WW_WZ_combi/combi_EWd6.cfg -c all -p cBL2 -m scan -f -r -20:20  -o combi_all_EWdim6

./bin/eftfun.py -i configs/WW_WZ_combi/ww_EWd6.cfg -c all -p cWWWL2 -m scan -f -r -3.5:3.5  -s cWWWL2=-1 -o ww_all_EWdim6
./bin/eftfun.py -i configs/WW_WZ_combi/ww_EWd6.cfg -c all -p cWL2 -m scan -f -r -8:5   -o ww_all_EWdim6
./bin/eftfun.py -i configs/WW_WZ_combi/ww_EWd6.cfg -c all -p cBL2 -m scan -f -r -20:20  -o ww_all_EWdim6

./bin/eftfun.py -i configs/WW_WZ_combi/wz_EWd6.cfg -c all -p cWWWL2 -m scan -f -r -3.5:3.5  -s cWWWL2=1 -o wz_all_EWdim6
./bin/eftfun.py -i configs/WW_WZ_combi/wz_EWd6.cfg -c all -p cWL2 -m scan -f -r -8:5   -o wz_all_EWdim6
./bin/eftfun.py -i configs/WW_WZ_combi/wz_EWd6.cfg -c all -p cBL2 -m scan -f -r -160:200  -o wz_all_EWdim6

./bin/eftfun.py -i configs/WW_WZ_combi/combi_EWd6.cfg -c all -p cWWWL2 -m scan -f -r -40:40  -s cWWWL2=1 -o combi_lin_EWdim6 -l
./bin/eftfun.py -i configs/WW_WZ_combi/combi_EWd6.cfg -c all -p cWL2 -m scan -f -r -5:5   -o combi_lin_EWdim6 -l
./bin/eftfun.py -i configs/WW_WZ_combi/combi_EWd6.cfg -c all -p cBL2 -m scan -f -r -120:100  -o combi_lin_EWdim6 -l

./bin/eftfun.py -i configs/WW_WZ_combi/ww_EWd6.cfg -c all -p cWWWL2 -m scan -f -r -160:0  -s cWWWL2=-1 -o ww_lin_EWdim6 -l
./bin/eftfun.py -i configs/WW_WZ_combi/ww_EWd6.cfg -c all -p cWL2 -m scan -f -r -15:8   -o ww_lin_EWdim6 -l
./bin/eftfun.py -i configs/WW_WZ_combi/ww_EWd6.cfg -c all -p cBL2 -m scan -f -r -120:100  -o ww_lin_EWdim6 -l

./bin/eftfun.py -i configs/WW_WZ_combi/wz_EWd6.cfg -c all -p cWWWL2 -m scan -f -r -30:50  -s cWWWL2=1 -o wz_lin_EWdim6 -l
./bin/eftfun.py -i configs/WW_WZ_combi/wz_EWd6.cfg -c all -p cWL2 -m scan -f -r -4:5   -o wz_lin_EWdim6 -l
./bin/eftfun.py -i configs/WW_WZ_combi/wz_EWd6.cfg -c all -p cBL2 -m scan -f -r -400:400  -o wz_lin_EWdim6 -l
