for i in {1..12}; do
    # Copy the original file to a new file
    ./bin/eftfun.py -i ./configs/CCDY_u35_1d_clipping/Tutorial_elmu_${i}.cfg -o ./results/u35_clipping/mT_elmu_lin_${i} -m fitpulls --linear > log_clipping.txt
done
