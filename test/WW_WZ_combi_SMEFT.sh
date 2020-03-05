
./bin/eftfun.py -i configs/WW_WZ_combi/combi.cfg -c all -m fit -f -o combi_fit_const
./bin/eftfun.py -i configs/WW_WZ_combi/wz.cfg -c all  -m fit -f -o wz_fit_const
./bin/eftfun.py -i configs/WW_WZ_combi/ww.cfg -c all  -m fit -f -o ww_fit_const

./bin/eftfun.py -i configs/WW_WZ_combi/combi.cfg -c cHq1,cHq3,cHl1,cHl3,cHe,cHd,cHu,clq1,clq3 -m fit -f -o combi_fit_aTGC
./bin/eftfun.py -i configs/WW_WZ_combi/wz.cfg -c cHq1,cHq3,cHl1,cHl3,cHe,cHd,cHu,clq1,clq3  -m fit -f -o wz_fit_aTGC
./bin/eftfun.py -i configs/WW_WZ_combi/ww.cfg -c cHq1,cHq3,cHl1,cHl3,cHe,cHd,cHu,clq1,clq3  -m fit -f -o ww_fit_aTGC

./bin/eftfun.py -i configs/WW_WZ_combi/combi.cfg -c cHl1,cHl3,cHe,clq1,clq3 -m fit -f -o combi_fit_aTGCVqq
./bin/eftfun.py -i configs/WW_WZ_combi/wz.cfg -c cHl1,cHl3,cHe,clq1,clq3  -m fit -f -o wz_fit_aTGCVqq
./bin/eftfun.py -i configs/WW_WZ_combi/ww.cfg -c cHl1,cHl3,cHe,clq1,clq3  -m fit -f -o ww_fit_aTGCVqq

./bin/eftfun.py -i configs/WW_WZ_combi/combi.cfg -c clq1,clq3 -m fit -f -o combi_fit_aTGCVqqVll
./bin/eftfun.py -i configs/WW_WZ_combi/wz.cfg -c clq1,clq3  -m fit -f -o wz_fit_aTGCVqqVll
./bin/eftfun.py -i configs/WW_WZ_combi/ww.cfg -c clq1,clq3  -m fit -f -o ww_fit_aTGCVqqVll


./bin/eftfun.py -i configs/WW_WZ_combi/combi.cfg -m fit -f -o combi_fit -l
./bin/eftfun.py -i configs/WW_WZ_combi/wz.cfg  -m fit -f -o wz_fit -l 
./bin/eftfun.py -i configs/WW_WZ_combi/ww.cfg  -m fit -f -o ww_fit -l


./bin/eftfun.py -i configs/WW_WZ_combi/combi.cfg -c all -m fit -f -o combi_pulls_const

# cW cHDD cHWB cHq1 cHq3 cHl1 cHl3 cHd cHu cll1;

for i in cW cHDD cHWB cHq1 cHq3 cHl1 cHl3 cHe cHd cHu cll1 clq1 clq3 ;
do
    ./bin/eftfun.py -i configs/WW_WZ_combi/combi.cfg -c all -p $i -m scan -f -o combi_quad
    ./bin/eftfun.py -i configs/WW_WZ_combi/wz.cfg -c all -p $i -m scan -f -o wz_quad
    ./bin/eftfun.py -i configs/WW_WZ_combi/ww.cfg -c all -p $i -m scan -f -o ww_quad
done
./bin/eftfun.py -i configs/WW_WZ_combi/combi.cfg -c all -p cW -m scan -f -o combi_quad -r -0.25:0.25 -s cW=0.1
./bin/eftfun.py -i configs/WW_WZ_combi/wz.cfg -c all -p cW -m scan -f -o wz_quad -r -0.25:0.25 -s cW=-0.1
./bin/eftfun.py -i configs/WW_WZ_combi/ww.cfg -c all -p cW -m scan -f -o ww_quad -r -0.25:0.25

./bin/eftfun.py -i configs/WW_WZ_combi/combi.cfg -c all -p cHq3 -m scan -f -o combi_quad -r -0.2:0.15 
./bin/eftfun.py -i configs/WW_WZ_combi/wz.cfg -c all -p cHq3 -m scan -f -o wz_quad -r -0.2:0.15 
./bin/eftfun.py -i configs/WW_WZ_combi/ww.cfg -c all -p cHq3 -m scan -f -o ww_quad -r -0.4:0.4

for i in cW cHDD cHWB cHq1 cHq3 cHl1 cHl3 cHe cHd cHu cll1 clq1 clq3 ;
do
    ./bin/eftfun.py -i configs/WW_WZ_combi/combi.cfg -c all -p $i -m scan -f -o combi_lin -l -r 3
    ./bin/eftfun.py -i configs/WW_WZ_combi/wz.cfg -c all -p $i -m scan -f -o wz_lin -l -r 3
    ./bin/eftfun.py -i configs/WW_WZ_combi/ww.cfg -c all -p $i -m scan -f -o ww_lin -l -r 3
done

for i in cW cHDD cHWB cHq1 cHq3 cHl1 cHl3 cHe cHd cHu cll1 clq1 clq3 
do
    ./bin/eftfun.py -i configs/WW_WZ_combi/combi.cfg -c cHq1,cHq3,cHl1,cHl3,cHd,cHu,clq1,clq3 -p $i -m scan -f -o combi_profATGC -l -r 3
    ./bin/eftfun.py -i configs/WW_WZ_combi/wz.cfg -c cHq1,cHq3,cHl1,cHl3,cHd,cHu,clq1,clq3 -p $i -m scan -f -o wz_profATGC -l -r 3
    ./bin/eftfun.py -i configs/WW_WZ_combi/ww.cfg -c cHq1,cHq3,cHl1,cHl3,cHd,cHu,clq1,clq3 -p $i -m scan -f -o ww_profATGC -l -r 3
done

for i in cW cHDD cHWB cHq1 cHq3 cHl1 cHl3 cHe cHd cHu cll1 clq1 clq3 
do
    ./bin/eftfun.py -i configs/WW_WZ_combi/combi.cfg -c cHq1,cHq3,cHl1,cHl3,cHDD,cHWB,cll1 -p $i -m scan -f -o combi_profVqq -l -r 3
    ./bin/eftfun.py -i configs/WW_WZ_combi/wz.cfg -c cHq1,cHq3,cHl1,cHl3,cHDD,cHWB,cll1 -p $i -m scan -f -o wz_profVqq -l -r 3
    ./bin/eftfun.py -i configs/WW_WZ_combi/ww.cfg -c cHq1,cHq3,cHl1,cHl3,cHDD,cHWB,cll1 -p $i -m scan -f -o ww_profVqq -l -r 3
done

for i in cW cHDD cHWB cHq1 cHq3 cHl1 cHl3 cHe cHd cHu cll1 clq1 clq3 
do
    ./bin/eftfun.py -i configs/WW_WZ_combi/combi.cfg -c cHDD,cHWB,cHq1,cHq3,cHd,cHu,cll1,clq1,clq3 -p $i -m scan -f -o combi_profVll -l -r 3
    ./bin/eftfun.py -i configs/WW_WZ_combi/wz.cfg -c cHDD,cHWB,cHq1,cHq3,cHd,cHu,cll1,clq1,clq3 -p $i -m scan -f -o wz_profVll -l -r 3
    ./bin/eftfun.py -i configs/WW_WZ_combi/ww.cfg -c cHDD,cHWB,cHq1,cHq3,cHd,cHu,cll1,clq1,clq3 -p $i -m scan -f -o ww_profVll -l -r 3
done


for i in cW cHDD cHWB cll1 cHq1 cHq3;
do
    for j in cW cHDD cHWB cll1 cHq1 cHq3;
    do
	if [ $i != $j ]
	then
	    ./bin/eftfun.py -i configs/WW_WZ_combi/combi.cfg -c all -p $i,$j -m scan -f -o combi_lin2D -l -r 3 -n 100
	    ./bin/eftfun.py -i configs/WW_WZ_combi/wz.cfg -c all -p $i,$j -m scan -f -o wz_lin2D -l -r 3
	    ./bin/eftfun.py -i configs/WW_WZ_combi/ww.cfg -c all -p $i,$j -m scan -f -o ww_lin2D -l -r 3
	fi
    done
done
