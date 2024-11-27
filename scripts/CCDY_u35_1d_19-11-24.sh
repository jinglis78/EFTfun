./bin/eftfun.py -i ./configs/CCDY_u35_1d_19-11-24/Tutorial_elmu.cfg -o ./results/u35_19-11-24/mT_elmu -m fitpulls > log_elmu.txt
./bin/eftfun.py -i ./configs/CCDY_u35_1d_19-11-24/Tutorial_el.cfg -o ./results/u35_19-11-24/mT_el -m fitpulls > log_el.txt
./bin/eftfun.py -i ./configs/CCDY_u35_1d_19-11-24/Tutorial_mu.cfg -o ./results/u35_19-11-24/mT_mu -m fitpulls > log_mu.txt

./bin/eftfun.py -i ./configs/CCDY_u35_1d_19-11-24/Tutorial_elmu.cfg -o ./results/u35_19-11-24/mT_elmu_lin -m fitpulls --linear > log_elmu_lin.txt
./bin/eftfun.py -i ./configs/CCDY_u35_1d_19-11-24/Tutorial_el.cfg -o ./results/u35_19-11-24/mT_el_lin -m fitpulls --linear > log_el_lin.txt
./bin/eftfun.py -i ./configs/CCDY_u35_1d_19-11-24/Tutorial_mu.cfg -o ./results/u35_19-11-24/mT_mu_lin -m fitpulls --linear > log_mu_lin.txt

