./bin/eftfun.py -i configs/EW_combi/combi_top.cfg -o EW_combi_workspace_linear_asimov -l -a -m workspace -f
./bin/eftfun.py -i configs/EW_combi/combi_top.cfg -o EW_combi_workspace_linear -l -m workspace -f
./bin/eftfun.py -i configs/EW_combi/combi_top.cfg -o EW_combi_workspace -m workspace -f
./bin/eftfun.py -i configs/EW_combi/combi_top.cfg -o EW_combi_workspace_asimov -a -m workspace -f

./bin/eftfun.py -i configs/EW_combi/combi_top.cfg -o EW_combi_workspace_mu_linear_asimov -l -a -m muworkspace -f
./bin/eftfun.py -i configs/EW_combi/combi_top.cfg -o EW_combi_workspace_mu_linear -l -m muworkspace -f

./bin/eftfun.py -i configs/EW_combi/combi_top.cfg -o EW_combi_workspace_para -l -m printpara -f
