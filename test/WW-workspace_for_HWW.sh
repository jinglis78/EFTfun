./bin/eftfun.py -i configs/ATLAS_WW_2019/profileSMEFTexample.cfg --override "main:General.coeffs=mu_ww,measurements:Measurement.scaling_factor=mu_ww,measurements:Prediction Uncertainties.{all} ww_norm= {rel} 0,measurements:Prediction Uncertainties.{all} ww_ewk= {rel} 0" -m workspace -f -o SMWW_workspace_mu_ww_noEFT

./bin/eftfun.py -i configs/ATLAS_WW_2019/profileSMEFTexample.cfg --override "main:General.coeffs=mu_ww,measurements:Measurement.scaling_factor=mu_ww,measurements:Prediction Uncertainties.{all} ww_norm= {rel} 0,measurements:Prediction Uncertainties.{all} ww_ewk= {rel} 0" -m workspace -f -a -o SMWW_workspace_mu_ww_noEFT_asimov

./bin/eftfun.py -i configs/ATLAS_WW_2019/profileSMEFTexample.cfg --override "main:General.coeffs=mu_ww cW cHG cHbox cHDD cHWB cHW cHq1 cHq3 cHl1 cHl3 cHd cHu cll1 cld clu clq1 clq3,measurements:Measurement.scaling_factor=mu_ww,measurements:Prediction Uncertainties.{all} ww_norm= {rel} 0,measurements:Prediction Uncertainties.{all} ww_ewk= {rel} 0" -m workspace -f -o SMWW_workspace_mu_ww_EFT

./bin/eftfun.py -i configs/ATLAS_WW_2019/profileSMEFTexample.cfg --override "main:General.coeffs=mu_ww cW cHG cHbox cHDD cHWB cHW cHq1 cHq3 cHl1 cHl3 cHd cHu cll1 cld clu clq1 clq3,measurements:Measurement.scaling_factor=mu_ww,measurements:Prediction Uncertainties.{all} ww_norm= {rel} 0,measurements:Prediction Uncertainties.{all} ww_ewk= {rel} 0" -m workspace -f -a -o SMWW_workspace_mu_ww_EFT_asimov

./bin/eftfun.py -i configs/ATLAS_WW_2019/profileSMEFTexample.cfg --override "main:General.coeffs=,measurements:Prediction Uncertainties.{all} ww_ewk= {rel} 0" -m workspace -f -o SMWW_workspace_ww_norm_noEFT

./bin/eftfun.py -i configs/ATLAS_WW_2019/profileSMEFTexample.cfg --override "main:General.coeffs=,measurements:Prediction Uncertainties.{all} ww_ewk= {rel} 0" -m workspace -f -a -o SMWW_workspace_ww_norm_noEFT_asimov

./bin/eftfun.py -i configs/ATLAS_WW_2019/profileSMEFTexample.cfg --override "measurements:Prediction Uncertainties.{all} ww_ewk= {rel} 0" -m workspace -f -o SMWW_workspace_ww_norm_EFT

./bin/eftfun.py -i configs/ATLAS_WW_2019/profileSMEFTexample.cfg --override "measurements:Prediction Uncertainties.{all} ww_ewk= {rel} 0" -m workspace -f -a -o SMWW_workspace_ww_norm_EFT_asimov
