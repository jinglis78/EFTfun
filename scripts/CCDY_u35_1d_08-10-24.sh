#./bin/eftfun.py -i ./configs/CCDY_u35_1d_08-10-24/Tutorial_elmu.cfg -o ./results/CCDY_top1d_tests/mT_comb_u35_08-10-24 -m fitpulls > log_elmu.txt
#./bin/eftfun.py -i ./configs/CCDY_u35_1d_08-10-24/Tutorial_el.cfg -o ./results/CCDY_top1d_tests/mT_el_u35_08-10-24 -m fitpulls > log_el.txt
./bin/eftfun.py -i ./configs/CCDY_u35_1d_08-10-24/Tutorial_mu.cfg -o ./results/CCDY_top1d_tests/mT_mu_u35_08-10-24 -m fitpulls > log_mu.txt

#./bin/eftfun.py -i ./configs/CCDY_u35_1d_08-10-24/Tutorial_elmu_lin.cfg -o ./results/CCDY_top1d_tests/mT_comb_lin_u35_08-10-24 --linear -m fitpulls > log_lin_elmu.txt
#./bin/eftfun.py -i ./configs/CCDY_u35_1d_08-10-24/Tutorial_mu_lin.cfg -o ./results/CCDY_top1d_tests/mT_mu_lin_u35_08-10-24 --linear -m fitpulls > log_lin_mu.txt
#./bin/eftfun.py -i ./configs/CCDY_u35_1d_08-10-24/Tutorial_el_lin.cfg -o ./results/CCDY_top1d_tests/mT_el_lin_u35_08-10-24 --linear -m fitpulls > log_lin_el.txt