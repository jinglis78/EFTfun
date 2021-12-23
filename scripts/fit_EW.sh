const='cHbox,cll,cll1,clj1,clj3,cee,ceu,ced,cle,clu,cld,cje'
const2='cHbox,cll,cll1,clj1,clj3,cee,ceu,ced,cle,clu,cld,cje,cHW,cHB,cHe,cHu,cHd'
#  WW WWj VBFZ WZ ZZ
#cW cHbox cHDD cHG cHW cHB cHWB cHl1 cHl3 cHe cHj1 cHj3 cHu cHd cll cll1 clj1 clj3 cee ceu ced cle clu cld cje
for m in WW VBFZ WZ ZZ combi 
do
    for i in cW cHDD cHWB cHl1 cHl3 cHe cll1 cHj1 cHj3 cHu cHd clj1 clj3 clu cld ceu ced cje cee cll cle cjj11 cjj18 cjj31 cjj38 cuu1 cuu8 cdd1 cdd8 cud1 cud8 cju1 cju8 cjd1 cjd8 cG
    do
	./bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_lin -a -l -c all
	./bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -s ${i}=0.005 -c all
	#./bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_prof_non4f -a -l -c $const
	#./bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_prof_some -a -l -c $const2
	
    done
done

m=WW
i=cW
./bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s cW=0.05 -r -0.3:0.3
m=WZ
./bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s cW=0.05 -r -0.3:0.3
m=ZZ
i=cHG
./bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.001 -r -0.02:0.02
m=WW
./bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.001 -r -0.2:0.2
i=cHW
./bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=1 -r -30:30
m=ZZ
./bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=1 -r -3:3

m=WW
i=cW
./bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s cW=0.05 -r -0.3:0.3
m=WZ
./bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s cW=0.05 -r -0.3:0.3
m=VBFZ
./bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.1 -r -0.4:0.4
m=combi
./bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.1 -r -0.3:0.3
m=ZZ
i=cHG
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.001 -r -0.02:0.02
m=combi
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.001 -r -0.02:0.02
m=WW
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.001 -r -0.2:0.2
i=cHW
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=1 -r -30:30
m=ZZ
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=1 -r -3:3
m=combi
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=1 -r -3:3
i=cHWB
m=ZZ
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.1 -r -0.3:0.3
m=WW
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.1 -r -6:6
i=cHl1
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=1 -r -100:100
i=cHj1
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.1 -r -0.4:0.4
i=cHj3
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.1 -r -0.5:0.5
m=WZ
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.1 -r -0.3:0.3
m=combi
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.1 -r -0.3:0.3
m=WW
i=cHu
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.1 -r -1:1
m=combi
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.1 -r -1:1
i=cHd
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.1 -r -1:1
m=WW
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.1 -r -1:1
m=WZ
i=cll
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.1 -r -300:300
i=cll1
m=VBFZ
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.1 -r -5:5
i=clj1
m=WW
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.05 -r -1:1
m=WZ
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.05 -r -5:5
m=ZZ
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.05 -r -1:1
i=clj3
m=WW
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.05 -r -1:1
m=WZ
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.05 -r -5:5
m=ZZ
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.05 -r -0.1:0.1
m=combi
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.05 -r -0.1:0.1
i=ceu
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.05 -r -1:1
i=cle
m=WZ
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.05 -r -30:30
m=WW
i=clu
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.05 -r -1:1
i=cld
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.05 -r -4:4
m=ZZ
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.05 -r -2:2
m=WZ
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.05 -r -1:1
m=combi
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.05 -r -1:1
i=cje
bin/eftfun.py -i configs/EW_combi/${m}.cfg -m scan -p $i -f -o ${m}_quad -a -c all -s $i=0.05 -r -5:5
