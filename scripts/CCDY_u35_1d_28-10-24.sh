./bin/eftfun.py -i ./configs/CCDY_u35_1d_clipping/Tutorial_mu.cfg -o ./results/u35_28-10-24/mT_mu -m fitpulls > log_mu.txt
./bin/eftfun.py -i ./configs/CCDY_u35_1d_clipping/Tutorial_el.cfg -o ./results/u35_28-10-24/mT_el -m fitpulls > log_el.txt
./bin/eftfun.py -i ./configs/CCDY_u35_1d_clipping/Tutorial_elmu.cfg -o ./results/u35_28-10-24/mT_elmu -m fitpulls > log_el.txt

./bin/eftfun.py -i ./configs/CCDY_u35_1d_clipping/Tutorial_mu.cfg -o ./results/u35_28-10-24/mT_mu_lin -m fitpulls --linear > log_mu.txt
./bin/eftfun.py -i ./configs/CCDY_u35_1d_clipping/Tutorial_el.cfg -o ./results/u35_28-10-24/mT_el_lin -m fitpulls --linear > log_el.txt
./bin/eftfun.py -i ./configs/CCDY_u35_1d_clipping/Tutorial_elmu.cfg -o ./results/u35_28-10-24/mT_elmu_lin -m fitpulls --linear > log_el.txt
