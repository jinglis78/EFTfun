#cW cHbox cHDD cHG cHW cHB cHWB cHl1 cHl3 cHe cHj1 cHj3 cHu cHd cll cll1 clj1 clj3 cee ceu ced cle clu cld cje
for i in cW cHj3 cHDD cHWB cHl1 cHl3 cHe cll1 cHj1 cHu cHd clj1 clj3 clu cld ceu ced cje cjj11 cjj18 cjj31 cjj38 cuu1 cuu8 cdd1 cdd8 cud1 cud8 cju1 cju8 cjd1 cjd8 cG
do
    for m in combi WW VBFZ WZ ZZ 
    do
	./bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_data_lin  -l -c all --robust 3
	./bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_data_quad  -c all  --robust 3
	
    done
done

./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan -p cHj3 -r -0.3:0.1 -f -c all --robust 3 -o combi_data_quad
./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan -p cje -r -0.25:0.25 -f -c all --robust 3 -o combi_data_quad
./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan -p cju1 -r -3:3 -s cju1=1 -f -c all --robust 3 -o combi_data_quad
./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan -p cju8 -r -8:8  -f -c all --robust 3 -o combi_data_quad
./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan -p cld -r -0.4:0.2 -f -c all --robust 3 -o combi_data_quad
./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan -p clj1 -r -0.3:0.3 -f -c all --robust 3 -o combi_data_quad
./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan -p clu -r -0.3:0.3 -f -c all --robust 3 -o combi_data_quad
./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan -p cuu1 -r -4:3 -f -c all --robust 3 -o combi_data_quad
./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan -p cld -r -0.4:0.4 -f -c all --robust 3 -o combi_data_quad
./bin/eftfun.py -i configs/EW_combi/ZZ.cfg -m scan -p cje -r -0.25:0.25 -f -c all --robust 3 -o ZZ_data_quad
./bin/eftfun.py -i configs/EW_combi/ZZ.cfg -m scan -p ced -r -0.6:0.6 -f -c all --robust 3 -o ZZ_data_quad
./bin/eftfun.py -i configs/EW_combi/ZZ.cfg -m scan -p clj1 -r -0.2:0.3 -f -c all --robust 3 -o ZZ_data_quad
./bin/eftfun.py -i configs/EW_combi/VBFZ.cfg -m scan -p cHd -r -20:20 -f -c all --robust 3 -o VBFZ_data_quad
./bin/eftfun.py -i configs/EW_combi/VBFZ.cfg -m scan -p cjd1 -r -4:4 -s cjd1=-2 -f -c all --robust 3 -o VBFZ_data_quad
./bin/eftfun.py -i configs/EW_combi/VBFZ.cfg -m scan -p cju8 -r -8:8 -f -c all --robust 3 -o VBFZ_data_quad
./bin/eftfun.py -i configs/EW_combi/VBFZ.cfg -m scan -p cud1 -r -8:8 -s cud1=2 -f -c all --robust 3 -o VBFZ_data_quad
./bin/eftfun.py -i configs/EW_combi/VBFZ.cfg -m scan -p cud8 -r -12:12 -f -c all --robust 3 -o VBFZ_data_quad
./bin/eftfun.py -i configs/EW_combi/VBFZ.cfg -m scan -p cG -r -3.5:3.5 -s cG=-1 -f -c all --robust 3 -o VBFZ_data_quad
./bin/eftfun.py -i configs/EW_combi/VBFZ.cfg -m scan -p cuu1 -r -4:5 -f -c all --robust 3 -o VBFZ_data_quad
./bin/eftfun.py -i configs/EW_combi/WW.cfg -m scan -p clj1 -r -0.3:0.3 -f -c all --robust 3 -o WW_data_quad
./bin/eftfun.py -i configs/EW_combi/WW.cfg -m scan -p clj3 -r -0.4:0.2 -f -c all --robust 3 -o WW_data_quad
./bin/eftfun.py -i configs/EW_combi/WZ.cfg -m scan -p clj1 -f -c all --robust 3 -o WW_data_quad
./bin/eftfun.py -i configs/EW_combi/WZ.cfg -m scan -p clj3 -f -c all --robust 3 -o WW_data_quad
