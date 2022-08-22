# setup:
git clone ssh://git@gitlab.cern.ch:7999/eft-tools/eft-fun.git -b WW2022
cd eft-fun

# maximum likelihood fit
./bin/eftfun.py -i configs/WW_reco/WW_reco.cfg -o WW_reco_test_fit -m fit
# profile likelihood scan of mu (-f skips the question of whether everything looks ok)
./bin/eftfun.py -i configs/WW_reco/WW_reco.cfg -o WW_reco_test_scan -m scan -p mu -f
# this is a fit without sytematics:
./bin/eftfun.py -i configs/WW_reco/WW_reco.cfg -o WW_reco_test_scan_nosyst -m scan -p mu -d all
# very important, optimize without data, using the asimov data set
./bin/eftfun.py -i configs/WW_reco/WW_reco.cfg -o WW_reco_test_scan_asimov -m scan -p mu -a
# sometimes we might want to set the fit range by hand
./bin/eftfun.py -i configs/WW_reco/WW_reco.cfg -o WW_reco_test_scan_fixrange -m scan -p mu -r -1.0:1.15
# systematics impact plot (takes a while!)
./bin/eftfun.py -i configs/WW_reco/WW_reco.cfg -o WW_reco_test_impacts -m impacts -p mu -f


# once eft effects are implemented 1D scans can be performed like this (keeping all other parameters constant):
./bin/eftfun.py -i configs/WW_reco/WW_reco.cfg -o WW_reco_scan -m scan -p c1 -c all
# 2D scans can be performed like this (keeping all other parameters constant, reducing the number of scan points):
./bin/eftfun.py -i configs/WW_reco/WW_reco.cfg -o WW_reco_scan2d -m scan -p c1,c2 -c all -n 20
# we can also profile over other parameters:
./bin/eftfun.py -i configs/WW_reco/WW_reco.cfg -o WW_reco_scan_profiled -m scan -p c1
# typically we need to keep at least a few constant when profiling...
./bin/eftfun.py -i configs/WW_reco/WW_reco.cfg -o WW_reco_scan_profilesome -m scan -p c1 -c c5,c7
