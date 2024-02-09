./bin/eftfun.py -i configs/ATLAS_boosted_ttbar_13TeV/config.cfg -o boosted_ttbar_workspace_linear_asimov -l -a -m workspace -f
./bin/eftfun.py -i configs/ATLAS_boosted_ttbar_13TeV/config.cfg -o boosted_ttbar_workspace_linear -l -m workspace -f
./bin/eftfun.py -i configs/ATLAS_boosted_ttbar_13TeV/config.cfg -o boosted_ttbar_workspace -m workspace -f
./bin/eftfun.py -i configs/ATLAS_boosted_ttbar_13TeV/config.cfg -o boosted_ttbar_workspace_asimov -a -m workspace -f

./bin/eftfun.py -i configs/ATLAS_boosted_ttbar_13TeV/config.cfg -o boosted_ttbar_workspace_mu_linear_asimov -l -a -m muworkspace -f
./bin/eftfun.py -i configs/ATLAS_boosted_ttbar_13TeV/config.cfg -o boosted_ttbar_workspace_mu_linear -l -m muworkspace -f

./bin/eftfun.py -i configs/ATLAS_boosted_ttbar_13TeV/config.cfg -o boosted_ttbar_workspace_para -l -m printpara -f
