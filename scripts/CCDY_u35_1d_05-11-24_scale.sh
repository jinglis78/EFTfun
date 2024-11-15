./bin/eftfun.py -i ./configs/CCDY_u35_1d_05-11-24_scale/Tutorial_elmu_env.cfg -o ./results/u35_05-11-24_scale/mT_elmu_lin_env -m fitpulls --linear > log_elmu.txt
./bin/eftfun.py -i ./configs/CCDY_u35_1d_05-11-24_scale/Tutorial_el_env.cfg -o ./results/u35_05-11-24_scale/mT_el_lin_env -m fitpulls --linear > log_el.txt
./bin/eftfun.py -i ./configs/CCDY_u35_1d_05-11-24_scale/Tutorial_mu_env.cfg -o ./results/u35_05-11-24_scale/mT_mu_lin_env -m fitpulls --linear > log_emu.txt

./bin/eftfun.py -i ./configs/CCDY_u35_1d_05-11-24_scale/Tutorial_elmu_orig.cfg -o ./results/u35_05-11-24_scale/mT_elmu_lin_orig -m fitpulls --linear > log_el_mulin.txt
./bin/eftfun.py -i ./configs/CCDY_u35_1d_05-11-24_scale/Tutorial_el_orig.cfg -o ./results/u35_05-11-24_scale/mT_el_lin_orig -m fitpulls --linear > log_el.txt
./bin/eftfun.py -i ./configs/CCDY_u35_1d_05-11-24_scale/Tutorial_mu_orig.cfg -o ./results/u35_05-11-24_scale/mT_mu_lin_orig -m fitpulls --linear > log_emu.txt
