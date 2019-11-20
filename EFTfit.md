# 3. EFT Fit

In this part of the tutorial, we will perform an EFT fit using histograms like the ones generated in the previous exercise (see [here](https://gitlab.cern.ch/eft-tools/smeft-jos/blob/master/EFTMCGeneration.md) for the MC generation exercise).
The EFT fit preparation consists of the following steps:

1. Setup the fitting code
1. Download HEP data and include it in fit config
1. Add Standard Model prediction
1. Add experimental uncertainties
1. Add theory uncertainties
1. Add EFT predictions
1. Perform fits
1. Reco measurement examples

## 3.1 Setting Up the Fitting Tool

### 3.1.1 Installation

The tool can be found here: 
<https://gitlab.cern.ch/eft-tools/eft-fun>
and downloaded with
```
git clone https://gitlab.cern.ch/eft-tools/eft-fun.git
# We will always assume you work from the base directory
cd eft-fun
```

Setup instructions can be found in the [README](https://gitlab.cern.ch/eft-tools/eft-fun/blob/master/README.md).
If possible, I would run it on a laptop, which requires ROOT with python bindings and RooFit installed.
Alternatively, it can also be run on lxplus, after setting up an appropriate release (see [README](https://gitlab.cern.ch/eft-tools/eft-fun/blob/master/README.md)).
TODO: We will also try to provide a Docker container.

To try out the code (just to check that it doesn't crash), run one of the examples from the [README](https://gitlab.cern.ch/eft-tools/eft-fun/blob/master/README.md), e.g.:
```
./bin/eftfun.py -i configs/ATLAS_WW_2019/SMEFTexample.cfg -m fit -c all -o example1
```

### 3.1.2 Introduction to the Code

This fitting tool can be used to perform EFT fits with unfolded measurements. In principle other fits of unfolded measurements can be performed as well (as they are typically simpler, with less POIs and often only linear parameter dependence). It is not planned to support folded fits (not because it is impossible but for more complicated problems a tried and tested tool like Histfitter might be more appropriate). 

The tool will implement a parameterization of EFT effects as a quadratic function, allows to include theory uncertainties on predictions, and has some additional convenience features. The preferred input format is measurements in the form of HEPdata and ROOT or Yoda histograms for theory predictions and theory uncertainties. Multiple measurements can be combined (but so far the correlation of experimental uncertainties is not supported, this is work in progress).

The fit model is built using config files.
The default statistical model is given by a multivariate Gaussian that represents measurement uncertainties, which are typically correlated between bins. The theory prediction fit to this multivariate Gaussian can contain nuisance parameters, constraint with normal distributions.

There are multiple functions to perform fits, parameter scans, statistical tests, and to debug the fit.

**It is by no means mandatory to perform SM EFT fits with this tool**. It is just simple enough to be used in this tutorial.
The statistical model created by this tool is not terribly complex (a printout is essentially human-readable), it should be possible to create a similar model from scratch if such an approach is preferred.

The tool has been validated for a limited number of use cases. Fit results from [arXiv:1808.06577 [hep-ph]](https://arxiv.org/abs/1808.06577) and [arXiv:1905.04242 [hep-ex]](https://arxiv.org/abs/1905.04242) have been reproduced. The tool is still work in progress. If you intend to use it for an analysis [let me know](mailto:hannes.mildner@cern.ch). I am happy to debug the code or add additional features on request.

### 3.1.3 Setup Config Files

The starting point for the config of this tutorial should be `configs/TopCfgTemplate.cfg` and `configs/MeasurementCfgTemplate.cfg`.
You can also have a look at the configs in `configs/ATLAS_WW_2019`, which are two EFT fits of the leading lepton pT distribution in WW production (using different models and assumptions).
Start by copying the templates:
```
mkdir configs/ATLAS_high_mass_DY_8TeV
cp configs/TopCfgTemplate.cfg configs/ATLAS_high_mass_DY_8TeV/Tutorial.cfg
cp configs/MeasurementCfgTemplate.cfg configs/ATLAS_high_mass_DY_8TeV/M_ee.cfg
```

In the `Tutorial.cfg` set a name, list the relevant coefficients you found in the first part of the tutorial, and set maximal parameter ranges to +-1000 (we assume that during MC generation Lambda was fixed to 1 TeV, so this corresponds to 1000/TeV^2, which is a pretty strong coupling). The config should look similar to this:
```
[General]
name = Tutorial
coeffs = c1 c2 c3 # actually the parameters from step one, case sensitive
measurement_configs = M_ee.cfg 
[Ranges]
c1 = -1000 1000
c2 = -1000 1000
c3 = -1000 1000
```
Instead of using the `[Ranges]` section, you could also set the `default_range`, which will be used if a parameter is given no range.
In `M_ee.cfg` you can set a measurement name, e.g.
```
name = M_ee
```
the rest will be filled in the following steps.

## 3.2 Add HEPdata

### 3.2.1 Download HEPdata

Go to <https://www.hepdata.net/record/ins1467454> and check out the HEPdata.
Download all in yaml format and put it in the hepdata subdirectory. From the command line:
```
wget https://www.hepdata.net/download/submission/ins1467454/1/yaml -O hepdata/HEPData-ins1467454-v1-yaml.tar.gz
tar xf hepdata/HEPData-ins1467454-v1-yaml.tar.gz -C hepdata/
```
We are interested in Table 18, you can already study it on the [webpage](https://www.hepdata.net/record/ins1467454?version=1&table=Table18).

### 3.2.2 Add HEPdata to Config

Set the central values of Table 18 as `measured` in your measurement config:
```
measured =  hepdata/HEPData-ins1467454-v1-yaml/Table_18.yaml
```
If you don't have the result of a measurement of interest in HEPdata --> prepare HEPdata ;) . If that shouldn't be possible, you can also add the cross section as a list of numbers, with one number per bin:
```
measured = 0.220630372493 0.100382409178 ...
```

The next three entries in the config, `covariance`, `correlation`, and `correlated_uncertainty` can be deleted.
In any case only `covariance` or the combination of `correlation`  and `correlated_uncertainty` should be used.
If the full covariance matrix of a measurement is published, this can be used directly in the fit. The combination of `correlation`  and `correlated_uncertainty` serves the same purpose.
For our measurement, the covariance matrix is not published, it will be constructed from the published uncertainty sources.
An alternative use of `covariance` / `correlation`+`correlated_uncertainty` is to encode the statistical uncertainties, which will be correlated between bins for an unfolded measurement. For simplicity, we will assume this correlation is negligible for this measurement.

## 3.3 Add SM Prediction

### 3.3.1 Best Prediction

First we will add a precise SM prediction, a fixed order NNLO QCD + NLO EWK calculation, that also contains photon-induced (PI) production, together with a born -> dressed correction, all taken from [the publication](https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PAPERS/STDM-2014-06/).
```
sm = 0.232434 0.106692 0.0506748 0.0264052 0.0145314 0.00831072 0.00472926 0.002126 0.000694617 0.000173075 0.000031465 0.000004004
```

### 3.3.2 LO Prediction

Our EFT predictions are calculated at leading order. We will make the assumption that QCD and EWK corrections to the SM prediction (differentially in mll) can also be applied to the prediction containing the effects of dimension six operators.
To implement this, we add the MadGraph LO SM histogram as shown below.
A histogram as the one generated in the previous exercise, but with more statistics, can be found in the `histos/ATLAS_high_mass_DY_8TeV` folder.
Add it to the configuration:
```
renorm_sm = histos/ATLAS_high_mass_DY_8TeV/Rivet_run_cll1_0_NPlq1.root:ATLAS_2016_I9999999/d18-x01-y01
# or, if you work with the YODA histograms
#renorm_sm = histos/ATLAS_high_mass_DY_8TeV/Rivet_run_cll1_0_NPlq1.yoda:/ATLAS_2016_I9999999/d18-x01-y01
```
We can then rescale effects by the sm/renorm_sm ratio (see below).
This can be, but surely isn't always, a good approach. One option, which we will implement, is to include the higher order corrections on the central SM+EFT prediction, but keep the scale uncertainties of the LO prediction.

### 3.3.3 More Setup

Finally, we add the following more cryptic options in the `[Measurement]` section
```
per_bin_width = True
stat_unc_threshold = 0.1
```
`per_bin_width = True` means that we don't use the bin content but bin content/width for all histogram inputs.
We need to do this since this is the way the cross section is published.
If the relative statistical uncertainty in a bin of a template corresponding to an EFT effect is above `stat_unc_threshold`, a nuisance parameter will be added that parametrizes this uncertainty in the fit (it is assumed that SM predictions are available with sufficient statistics and no such feature is implemented for them).
The `scaling_factor` option can be removed.
One can add a `binning` but this will only affect plots and isn't always helpful.

## 3.4 Add Experimental Uncertainties
We will skip to the `[Measurement Uncertainties]` section now.
By looking at the `Table_18.yaml` file or the [HEPdata page](https://www.hepdata.net/record/ins1467454?version=1&table=Table18) you will see that the measurement has fifteen uncertainty sources. All uncertainties affecting the measurement, including theory uncertainties on background or unfolding, will be classified as measurement uncertainties.

We will let the fitting code construct a covariance matrix from these uncertainties. To do so, it needs to know which uncertainties are correlated and which are uncorrelated between bins (uncertainties that are partially correlated need to be tackled by adding the covariance matrix as described in section 3.2).
Add uncertainties in the following fashion:
```
stat =  {uncorr} stat
trig_cor =  {corr} sys,trig cor
trig_unc =  {uncorr} sys,trig unc
reco_cor = {corr} sys,reco cor
id_cor = {corr} sys,id cor
iso_cor = {corr} sys,iso cor
iso_unc = {uncorr} sys,iso unc
eres_cor = {corr} sys,Eres cor
escale_cor = {corr} sys,Escale cor
mult_cor = {corr} sys,mult. cor
mult_unc = {uncorr} sys,mult. unc
... (four more to add)
```
The left-hand side name needs to be unique but doesn't serve a specific purpose so far. In the future it could be used to correlate uncertainties of different measurements.
On the right hand side it needs to be specified in braces whether an uncertainty is correlated or uncorrelated and the label, which can be found in the yaml file or on the HEPdata page, has to be written.

Instead of just using the label (e.g. `sys,trig cor`) one can also use the full path `hepdata/HEPData-ins1467454-v1-yaml/Table_18.yaml:sys,trig cor` (this is necessary if the central value is not taken from hepdata). One can also write uncertainties as a list of numbers. In that case one needs to specify whether the uncertainty is relative or absolute, e.g.:
```
dy_ee_stat =  {uncorr rel} 0.01 0.02 ...
```
There is also the option to not include the uncertainty in the covariance matrix but to introduce a nuisance parameter in the fit that shifts the measured cross section, with the option 'np'. This nuisance parameter can then be correlated, also with reco level fits.

## 3.5 Add Theory Uncertainties
Now we will add theory uncertainties, i.e. uncertainties that affect predictions for dilepton production, be it SM DY or due to dimension-six operator effects.
These will be included in the fit as nuisance parameters. This is done in the `[Prediction Uncertainties]` section.

In the publication, uncertainties on the SM prediction are estimated to be

Uncertainty | Relative effect (symmetrized)
--- | ---
PDF | 0.014 0.014 0.015 0.016 0.016 0.017 0.018 0.020 0.022 0.025 0.028 0.031
alpha_S |  0.009 0.008 0.008 0.007 0.007 0.006 0.006 0.005 0.004 0.002 0.001 0.003
scale | 0.005 0.005 0.004 0.003 0.003 0.003 0.002 0.003 0.003 0.004 0.004 0.005
P.I. |0.005 0.008 0.011 0.015 0.019 0.022 0.026 0.031 0.040 0.051 0.071 0.0112

The implementation of theory uncertainties in the config file works similar to the implementation of experimental uncertainties. In the config file we now additionally need to specify the affected parts of the prediction on the left hand side. It can be `{all}`, `{sm}`, an EFT coefficent `{c1}` or a pair `{c1,c2}` if the uncertainty affects only certain linear or quadratic effects.
In contrast to the experimental uncertainties, the left-hand-side name is of relevance: it corresponds to the nuisance parameter name. This fact can be used to correlate theory uncertainties affecting different processes (or even different measurements).
On the right hand side the following options are possible:

- `rel`/`abs`: Relative or linear effect.
- `lin`/`exp`: The type of interpolation and extrapolation for the effect. Linear might be more intuitive, exponential can be more stable (it can't become negative) and might be a better choice for large rate uncertainties. Nuisance parameters have Gaussian constraints but using exponential interpolation + Gaussian constraint should be the same as linear interpolation + LogNormal constraint.
- `keepnorm`/`renorm`: Apply the sm/renorm_sm correction to the effect? 

Fill the section in the following fashion assuming that alpha_S has the same effect on all parts of the prediction but the PI uncertainty only affects the SM:
```
{all} dy_pdf = {rel lin keepnorm} 0.014 0.014 0.015 0.016 0.016 0.017 0.018 0.020 0.022 0.025 0.028 0.031
{sm} dy_scale = {rel lin keepnorm} 0.005 0.005 0.004 0.003 0.003 0.003 0.002 0.003 0.003 0.004 0.004 0.005
... (two more to add)
```
For good measure, we will also add the following uncertainties (these are dummy uncertainties and don't correspond to what one would use in a serious fit)

- A 1% flat uncertainties affecting everything, to account for uncertainties due to non-perturbative effects (remember we compare a fixed order prediction with a fudge factor to particle level)
- Two scale uncertainties, one of 5% on all linear and one of 10% on all quadratic EFT effects. All linear effects (and all quadratic effects) will be affected correlated but there is no correlation between linear and quadratic. This scheme is extremely simplistic and ad hoc, this is a part where more study and theory input would be needed. Even though it wont be a big effect: the uncertainty should be rescaled by sm/renorm_sm, as we do the same with the EFT effects.

```
{all} dy_np = {rel lin keepnorm} 0.01
{cll1} dy_4f_operator_scale_lin = {rel lin renorm} 0.05
{clq1} dy_4f_operator_scale_lin = {rel lin renorm} 0.05
{clq3} dy_4f_operator_scale_lin = {rel lin renorm} 0.05
{cll1,cll1} dy_4f_operator_scale_sq = {rel lin renorm} 0.1
{clq1,clq1} dy_4f_operator_scale_sq = {rel lin renorm} 0.1
{clq3,clq3} dy_4f_operator_scale_sq = {rel lin renorm} 0.1
{cll1,clq1} dy_4f_operator_scale_sq = {rel lin renorm} 0.1
{cll1,clq3} dy_4f_operator_scale_sq = {rel lin renorm} 0.1
{clq1,clq3} dy_4f_operator_scale_sq = {rel lin renorm} 0.1
```
This is a case where wildcards, documented in `configs/TopCfgTemplate.cfg`, could also be useful to compactify the notation, especially once more operators are involved.

## 3.6 Add EFT Predictions
Finally, we will add the EFT templates generated in the previous exercise.
There are already higher stats templates prepared for you in the `histos/ATLAS_high_mass_DY_8TeV/` folder, the naming convention should be mostly clear from the previous exercise.

### 3.6.1 Linear Effects
Add the linear EFT effects (the ones generated with `NP^2==1`) in the `[Linear Effects]` section in the following fashion.
```
cll1 = {renorm abs} histos/ATLAS_high_mass_DY_8TeV/Rivet_run_cll1_1_NPSQeq1.root:ATLAS_2016_I9999999/d18-x01-y01
clq1 = {renorm abs} histos/ATLAS_high_mass_DY_8TeV/Rivet_run_clq1_1_NPSQeq1.root:ATLAS_2016_I9999999/d18-x01-y01
... (one more to add)
```
Again there is the option to do `renorm` or `keepnorm` (we settled for `renorm`) and `abs` or `rel` (which makes more sense if effects are given as numbers instead of histograms).

### 3.6.1 Quadratic Effects
Quadratic effects are added in the same way, only that two coefficients are on the left hand side:
```
cll1 cll1 = {renorm abs} histos/ATLAS_high_mass_DY_8TeV/Rivet_run_cll1_1_NPeq1.root:ATLAS_2016_I9999999/d18-x01-y01
clq1 clq1 = {renorm abs} histos/ATLAS_high_mass_DY_8TeV/Rivet_run_clq1_1_NPeq1.root:ATLAS_2016_I9999999/d18-x01-y01
clq3 clq3 = {renorm abs} histos/ATLAS_high_mass_DY_8TeV/Rivet_run_clq3_1_NPeq1.root:ATLAS_2016_I9999999/d18-x01-y01
cll1 clq1 = {renorm abs} histos/ATLAS_high_mass_DY_8TeV/Rivet_run_cll1_1_clq1_1_cross.root:ATLAS_2016_I9999999/d18-x01-y01
cll1 clq3 = {renorm abs} histos/ATLAS_high_mass_DY_8TeV/Rivet_run_cll1_1_clq3_1_cross.root:ATLAS_2016_I9999999/d18-x01-y01
... (one more to add)
```
Only one of `c1 c2` , `c2 c1` need to be added. However, if you add both the code will automatically skip the non-existing hitogram.

## 3.7 Perform Fits
Now that the model is built we can run fits and statistical tests.
The fitter prints usage instruction when calling:
```
./bin/eftfun.py
```

First, run a maximum likelihood fit, with all EFT parameters fixed to zero:
```
./bin/eftfun.py -i configs/ATLAS_high_mass_DY_8TeV/Tutorial.cfg -m fit -c all
```
The fitter will ask you to check the statistical model.
You should do this and you can also use the opportunity to try to understand the model that is built.
Once the fit is run, the `results` folder will contain post-fit plots and a correlation matrix.

Then, fit one parameter at a time or even all of them:
```
./bin/eftfun.py -i configs/ATLAS_high_mass_DY_8TeV/Tutorial.cfg -m fit -c all -p clq3
./bin/eftfun.py -i configs/ATLAS_high_mass_DY_8TeV/Tutorial.cfg -m fit
```

Additional information on nuisance parameters and POIs can be gathered from pull plots:
```
./bin/eftfun.py -i configs/ATLAS_high_mass_DY_8TeV/Tutorial.cfg -m pulls -c all
./bin/eftfun.py -i configs/ATLAS_high_mass_DY_8TeV/Tutorial.cfg -m pulls 
```

Estimate a confidence interval on one parameter with the Profilelikelihood method and do the same once without theory uncertainties and once in a linearised approach (where quadratic EFT coefficients are dropped):
```
./bin/eftfun.py -i configs/ATLAS_high_mass_DY_8TeV/Tutorial.cfg -m profilelikelihood -c all -p clq3
./bin/eftfun.py -i configs/ATLAS_high_mass_DY_8TeV/Tutorial.cfg -m profilelikelihood -c all -p clq3 -d all
./bin/eftfun.py -i configs/ATLAS_high_mass_DY_8TeV/Tutorial.cfg -m profilelikelihood -c all -p clq3 -l
```

More information about the structure of the minimum can be gathered by running a scan. Even more info will be stored with the `-w` flag (checkout the results). Also run an Asimov fit with `-a`:
```
./bin/eftfun.py -i configs/ATLAS_high_mass_DY_8TeV/Tutorial.cfg -m scan -c all -p clq3
./bin/eftfun.py -i configs/ATLAS_high_mass_DY_8TeV/Tutorial.cfg -m scan -c all -p clq3 -w
./bin/eftfun.py -i configs/ATLAS_high_mass_DY_8TeV/Tutorial.cfg -m scan -c all -p clq3 -a
```

Now run a 2D scan to better understand parameter correlations (this takes a bit longer):
```
./bin/eftfun.py -i configs/ATLAS_high_mass_DY_8TeV/Tutorial.cfg -m scan -c all -p clq3,cll1
./bin/eftfun.py -i configs/ATLAS_high_mass_DY_8TeV/Tutorial.cfg -m scan -c all -p clq3,clq1
```
The last one looks interesting, let's adapt the scan range:
```
./bin/eftfun.py -i configs/ATLAS_high_mass_DY_8TeV/Tutorial.cfg -m scan -c all -p clq3,clq1 -r -0.12:0,-0.1:0.15
```

Finally, one can create a RooStats workspace and model config with the following command.
```
./bin/eftfun.py -i configs/ATLAS_high_mass_DY_8TeV/Tutorial.cfg -m workspace
```
You can print its contents with
```
root -l results/Tutorial.root
RooWorkspace* w = (RooWorkspace*)_file0->Get("w")
w->Print()
```
The workspace could be combined with other measurements (also fits of folded data) or used to calculate limits with frequentist methods, e.g. using [this](https://root.cern.ch/doc/v614/StandardHypoTestInvDemo_8C.html) script.

It is also planned to add functions to implement frequentist and hybrid limits using the appropriate RooStats calculator directly to this fitting tool.

## 3.8 Reco measurement examples (not part of tutorial)

A reco level measurement has to be declared with
```
reco_level = True
```
Measured and predicted values are then given by the number of observed events instead of measured cross section. There is no covariance matrix. There is a new category 'Backgrounds' All uncertainties (experimental and theoretical) are 'Prediction Uncertainties' (if they affect the sm prediction) or 'Background Uncertainties'.
See below for an example of a reco level measurement:
```
[Measurement]
name = reco_measurement1
measured = 1000 1000
sm = 850 880
reco_level = True

[Linear Effects]
c1 = 1 0

[Backgrounds]
b1 = 100 100

[Background Uncertainties]
{b1} b1_norm = {rel} 0.3 0.3
```
and the equivalent unfolded measurement:
```
[Measurement]
name = unfold_measurement1
measured = 0.9 0.9
sm = 0.85 0.88
covariance = 0.001 0, 0 0.001

[Linear Effects]
c1 = 0.001 0

[Measurement Uncertainties]
{all} b1_norm = {rel flip np} 0.03 0.03
```