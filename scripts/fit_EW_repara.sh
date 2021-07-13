const=cVff8,c2q2l4,c2q2l5,c2q2l6,c4q1,c4q2,c4q3,c4q4,c4q5,c4q6,c4q7,c4q8,c4q9,c4q10,c4q11,c4q12,c4q13,c4q14
fit="c3W cHq3 cVff0 cVff1 cVff2 cVff3 cVff4 cVff5 cVff6 cVff7 c2q2l0 c2q2l1 c2q2l2 c2q2l3 c4q0"
all="c3W cHq3 cVff0 cVff1 cVff2 cVff3 cVff4 cVff5 cVff6 cVff7 c2q2l0 c2q2l1 c2q2l2 c2q2l3 c4q0 cVff8 c2q2l4 c2q2l5 c2q2l6 c4l0 c4l1 c4l2 c4q1 c4q2 c4q3 c4q4 c4q5 c4q6 c4q7 c4q8 c4q9 c4q10 c4q11 c4q12 c4q13 c4q14"


for i in $fit
do
    ./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan -a -p $i --repara -o repara_prof_quad -f -c $const,cVff7 -n 50 --implicitscan cVff7:5:-5:5 --robust 3
    ./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan -p $i --repara -o repara_data_prof_quad -f -c $const,cVff7 -n 50 --implicitscan cVff7:5:-5:5 --robust 3
done


for i in $all
do
    ./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan -a -c all -p $i --repara -l -o repara_lin -f --robust 3
    ./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan -a -c all -p $i --repara -o repara_quad -f --robust 3
done

for i in $all
do
    ./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan -c all -p $i --repara -l -o repara_data_lin -f --robust 3
    ./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan -c all -p $i --repara -o repara_data_quad -f --robust 3
done


for i in $fit
do
    ./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan -a -p $i --repara -l -o repara_prof_lin -f -c $const --robust 3
    ./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan -p $i --repara -l -o repara_data_prof_lin -f -c $const --robust 3
done


./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan --repara -f -c $const -o repara_data_prof_quad -p cHq3 -r -0.35:0.1 --robust 3 -n 100 --implicitscan cVff7:30:-5:5
./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan --repara -f -c $const -o repara_data_prof_quad -p c3W -r -0.2:0.2 --robust 3 -n 100 --implicitscan cVff7:30:-5:5
./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan --repara -f -c $const -o repara_data_prof_quad -p cVff0 -r -0.25:0.35 --robust 3 -n 100 --implicitscan cVff7:30:-5:5
./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan --repara -f -c $const -o repara_data_prof_quad -p cVff1 -r -0.65:0.25 --robust 3 -n 100 --implicitscan cVff7:30:-5:5
./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan --repara -f -c $const -o repara_data_prof_quad -p cVff2 -r -1:1.2 --robust 3 -n 100 --implicitscan cVff7:30:-5:5
./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan --repara -f -c $const -o repara_data_prof_quad -p c2q2l0 -r -0.3:0.1 --robust 3 -n 100 --implicitscan cVff7:30:-5:5
./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan --repara -f -c $const -o repara_data_prof_quad -p c2q2l1 -r -0.4:0.4 --robust 3 -n 100 --implicitscan cVff7:30:-5:5
./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan --repara -f -c $const -o repara_data_prof_quad -p c4q0 -r -0.2:0.15 --robust 3 -n 100 --implicitscan cVff7:30:-5:5

./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan --repara -f -c $const -o repara_data_prof_quad -p cVff5 -r -2.5:1.2 --robust 3 -n 100 --implicitscan cVff7:30:-5:5

./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan --repara -f -c $const -o repara_prof_quad -p cHq3 -r -0.3:0.3 --robust 3 -n 100 -a --implicitscan cVff7:30:-5:5
./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan --repara -f -c $const -o repara_prof_quad -p c3W -r -0.2:0.2 --robust 3 -n 100 -a --implicitscan cVff7:30:-5:5
./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan --repara -f -c $const -o repara_prof_quad -p cVff0 -r -0.3:0.3 --robust 3 -n 100 -a --implicitscan cVff7:30:-5:5
./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan --repara -f -c $const -o repara_prof_quad -p cVff1 -r -0.4:0.4 --robust 3 -n 100 -a --implicitscan cVff7:30:-5:5
./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan --repara -f -c $const -o repara_prof_quad -p cVff2 -r -1:1 --robust 3 -n 100 -a --implicitscan cVff7:30:-5:5
./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan --repara -f -c $const -o repara_prof_quad -p c2q2l0 -r -0.3:0.3 --robust 3 -n 100 -a --implicitscan cVff7:30:-5:5
./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan --repara -f -c $const -o repara_prof_quad -p c2q2l1 -r -0.4:0.4 --robust 3 -n 100 -a --implicitscan cVff7:30:-5:5
./bin/eftfun.py -i configs/EW_combi/combi.cfg -m scan --repara -f -c $const -o repara_prof_quad -p c4q0 -r -0.15:0.15 --robust 3 -n 100 -a --implicitscan cVff7:30:-5:5

