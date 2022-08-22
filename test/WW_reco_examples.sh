# Needs a recent root version with RooFit! I use this code on my laptop with root 6.22, should also work on lxplus with lsetup root

# setup
git clone ssh://git@gitlab.cern.ch:7999/eft-tools/eft-fun.git -b WW2022
cd eft-fun

# These config files are the most important thing to study, they contain many additional comments
# configs/WW_reco/WW_reco.cfg
# configs/WW_reco/WW_incl_htlep.cfg

# Additional help:
./bin/eftfun.py -h
# Fits can be performed like below

# maximum likelihood fit (output are postfit plots and correlation matrixes)
./bin/eftfun.py -i configs/WW_reco/WW_reco.cfg -o WW_reco_test_fit -m fit
# profile likelihood scan of mu (-f skips the question of whether everything looks ok)
./bin/eftfun.py -i configs/WW_reco/WW_reco.cfg -o WW_reco_test_scan -m scan -p mu -f
# this is a fit without sytematics:
./bin/eftfun.py -i configs/WW_reco/WW_reco.cfg -o WW_reco_test_scan_nosyst -m scan -p mu -d all -f
# very important, optimize without data, using the asimov data set
./bin/eftfun.py -i configs/WW_reco/WW_reco.cfg -o WW_reco_test_scan_asimov -m scan -p mu -a -f
# sometimes we might want to set the fit range by hand
./bin/eftfun.py -i configs/WW_reco/WW_reco.cfg -o WW_reco_test_scan_fixrange -m scan -p mu -r 0:0.2 -f
# nuisance parameter pulls (and fit results for parameters of interest)
./bin/eftfun.py -i configs/WW_reco/WW_reco.cfg -o WW_reco_test_scan -m pulls -f
# systematics impact plot (takes a while and prints out a lot of stuff...)
./bin/eftfun.py -i configs/WW_reco/WW_reco.cfg -o WW_reco_test_impacts -m impacts -p mu -f
# calculate confidence intervals with toys (takes ages! even for 10 steps and 200 toys each -- this uses 4 cores, better would be cluster (needs to be implemented)
# follow the progress in results/WW_reco_test_scan_toys/job_?/log
./bin/eftfun.py -i configs/WW_reco/WW_reco.cfg -o WW_reco_test_scan_toys -m scan -p mu -r 0:0.2 -n 10 -t 200 --njobs 4 -f


# once eft effects are implemented 1D scans can be performed like this (keeping all other parameters constant):
./bin/eftfun.py -i configs/WW_reco/WW_reco.cfg -o WW_reco_scan -m scan -p c1 -c all
# 2D scans can be performed like this (keeping all other parameters constant, reducing the number of scan points):
./bin/eftfun.py -i configs/WW_reco/WW_reco.cfg -o WW_reco_scan2d -m scan -p c1,c2 -c all -n 20
# we can also profile over other parameters:
./bin/eftfun.py -i configs/WW_reco/WW_reco.cfg -o WW_reco_scan_profiled -m scan -p c1
# typically we need to keep at least a few constant when profiling...
./bin/eftfun.py -i configs/WW_reco/WW_reco.cfg -o WW_reco_scan_profilesome -m scan -p c1 -c c5,c7
