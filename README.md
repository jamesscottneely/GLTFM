# GLTFM
This repository contains the code for the Generalized Long Term Fault Memory (GLTFM) model described in Neely et al. *[Submitted]*. GLTFM allows users to calculate large earthquake probabilities. GLTFM incorporates the impact of residual strain after an earthquake in these calculations. Users input a paleoseismic record with the inter-event times of past large earthquakes and GLTFM calculates earthquake probabilities for the current quiescent period. If you have any questions/comments/concerns, please reach out to me at jneely@uchicago.edu

## Installation
Because GLTFM runs on Python, we recommend using either Anaconda (https://www.anaconda.com/download/) or Miniconda (https://docs.conda.io/en/latest/miniconda.html) to create and manage an environment to run the GLTFM package  After downloading this repository, create a conda environment with the necessary dependencies using the ```GLTFM_env.yml``` file provided. Using the Terminal window, navigate to the GLTFM repository and run the following ```conda env create -f GLTFM_env.yml``` to create the environment. See the Examples section below to test the succesful installation of the code. NOTE: This package has only been tested on a macOS environment and is not guaranteed to work in Windows.

## Input Files
*PaleoRecord.txt*: [Required] This file contains one column with the observed earthquake inter-event times (in years) in chronological order. For example if earthquakes occured in 1900, 1950, 1962, and 2000, this file would contain the following inter-event times:
```
50
12
38
```
Note: There is always one fewer inter-event times than earthquakes.

*ParamFile.txt*: [Required] This file contains the input parameters needed to run GLTFM. This file contains 10 input parameters:

```currentYr```: The current year of interest for the probability calculations. 

```priorEQ```: The year when the most recent earthquake in the paleoseismic record occured. 

```T_Forecast```: Window of time (in years) of interest for earthquake probability forecasts. For example, the probability of an earthquake in the next 30 years.

```yrsForecast```: Number of years after the most recent earthquake to calculate the probabilities.

```z0_flag```: Flag for the initial strain $Z_0$ (in years) after the first earthquake in the paleoseismic record. Either set a value (e.g. ```10```) to fix it or set to ```Free``` to have it estimated via maximum likelihood.

```R_flag```: Flag for the strain drop $R$ (in years) after an earthquake in the paleoseismic record. This flag can be set to a specific value (e.g. ```150```), estimated via maximum likelihood by setting to ```Free```, or estimated using tectonic information. To estimate R using tectonic information--using the Wells and Coppersmith 1994 equation, you must include the estimated earthquake magnitude, fault geometry (SS, NO, R, All) and slip accumulation rate (in mm/yr). The input should be within brackets and look like ```[7.8 SS 50]```. Future implementations will allow for variable values for $R$.

```incOpenInt```: This is a Y or N flag to indicate whether the quiescent period after the most recent earthquake should be included in the maximum likelihood estimation. Usually set to Y as the open interval can provide additional information.

```model```: This is the base probability model used for GLTFM. The model name must correspond to one of the 2-parameter continuous probability models in ```scipy.stats```. As discussed in the manuscript we prefer using the Weibull model which in ```scipy.stats``` is ```weibull_min```.

```iterFit```: This is a Y or N flag to indicate whether the free parameters should be iteratively solved for (as described in Neely et al. *submitted*) or if they should be solved for all at once.

```params```: This list ```[shape scale]``` is for the shape and scale parameters of GLTFM's underlying probability distribution. Either a specific value should be specified or set to ```Free``` to estimate via maximum likelihood. For example ```[Free Free]``` or ```[2.0 Free]```. See ```scipy.stats``` for a description of the shape and scale parameters for the different probability distributions.

*PaleoRecord_R.txt*: [In development] This file contains one column with the estimated earthquake $R$ values in chronological order. Each earthquake needs an assigned $R$ value. This allows for variable $R$ values for the different earthquakes. ```R_flag``` must be set to ```Variable``` to use this feature. NOTE: This functionality is still under development and has not been tested fully.

## GLTFM Package Structure and Key Functions
The GLTFM package consists of four modules. Below is a brief description of each module and some of the functions that they contain. 

```readInputs.py``` processes the input files for GLTFM. It contains  ```paleoRead``` which loads the paleoseismic record and ```paramRead``` which reads and processes the *ParamFile.txt* file. It also contains ```WC_equation``` which converts tectonic information into an estimate of $R$ using the Wells and Coppersmith (1994) equations.

```fittingGLTFM.py``` contains the functions that calculate earthquake probabilities using both the GLTFM model and existing renewal models (such as Lognormal or Brownian Passage Time). ```fit_GLTFM``` fits the GLTFM model to the supplied paleoseismic record and returns the best fitting parameters ($R$, $Z$, shape, scale depending on what was specified in *ParamFile.txt*). For specific $R$, initial strain $Z_0$, and inter-event times, ```z_val``` returns the estimated residual strain (in years) after each earthquake. This module also contains the ```GLTFM_gen``` class which when provided with residual strain $Z$, shape and scale parameters, and GLTFM's base probability model (usually set to weibull_min) has the functionality to estimate probability density functions [```GLTFM_gen.pdf()```], cumulative distribution functions [```GLTFM_gen.cdf()```], instantaneous hazard rates [```GLTFM_gen.haz_z()```], and T-year forecasts (30, 50, etc.) [```GLTFM_gen.Tforecast()```]. This module also contains ```fit_Renewal``` which fits the specified renewal model to the paleoseismic record.

```plotFunctions.py``` contains the ```plotHaz``` function which plots the GLTFM hazard rate history curve for the paleoseismic record. The plot also indicates the amount of residual strain (in years) after each each earthquake.

```uncertainty.py``` contains the functions to estimate the uncertainty of the forecast models. Forecast uncertainty is estimated by generating synthetic paleoseismic records with set parameters and length, estimating the parameters for the synthetic record, and then applying the estimated parameters to the real paleoseismic record. ```GLTFM_uncertainty``` performs this uncertainty analysis using ```sim_GLTFM``` which generates synthetic paleoseismic records. ```Renewal_uncertainty``` and ```sim_Renewal``` peform a similar function for the renewal models.

## Examples
We have provided the ```Applying_GLTFM.py``` script to help users implement GLTFM. This script takes a paleoseismic record, calculates earthquake probabilities and uncertainties, and returns several plots and files of interest. This script has been designed to show the range of calculations available and can be updated accordingly to the user's needs. It also compares GLTFM to the existing lognormal, Brownian Passage Time, and Weibull models. In the ```exampleData``` folder, we have provided the San Andreas Fault Pallet Creek paleoseismic record (Scharer et al., 2011) and a sample parameter file for testing. To run the script, please change ```dataPath```, ```include_uncertainty```, ```seed_val```, and ```num_sims``` fields as indicated in the script.

Here is a brief summary of the output files and plots generated by the ```Applying_GLTFM.py``` script:

### Output Files
```R_out.txt```: $R$ value for each earthquake in the paleoseismic sequence in chronological order.

```z_out.txt```: Residual strain (in years) after each in the paleoseismic sequence in chronological order.

```model_out.txt```: For GLTFM and the selected renewal models, this file indicates the best fitting shape and scale parameters, the corresponding loglikelihood value found in the maximum likelihood estimation, and the ##yrsForecast indicating the probability of an earthquake occuring in the next ## years after the ```currentYr```.

```PDF_out.txt```: Inter-event time probability density function estimates for the current quiescent period for GLTFM and renewal models.

```CDF_out.txt```: Inter-event time cumulative distribution function estimates for the current quiescent period for GLTFM and renewal models.

```yrsForecast_out.txt```: Probability of an earthquake occuring in the next ## years for the current quiescent period for GLTFM and renewal models.

```GLTFM_UncertaintyParams_out.txt```: GLTFM $R$, $Z_0$, shape parameter, and scale parameter estimates from each of the simulated catalogs in the uncertainty analysis. Note that $Z_0$ is not intrinsic to the underlying parameters so it estimated for the real catalog based on $R$, shape parameter, and scale parameter.

```GLTFM_UncertaintyPDF_out.txt```: Estimated GLTFM inter-event time PDFs for the current quiescent period when the $R$, $Z_0$, shape parameter, and scale parameter estimates from ```GLTFM_UncertaintyParams_out.txt``` are applied to the real paleoseismic record.

```GLTFM_UncertaintyCDF_out.txt```: Estimated GLTFM inter-event time CDFs for the current quiescent period when the $R$, $Z_0$, shape parameter, and scale parameter estimates from ```GLTFM_UncertaintyParams_out.txt``` are applied to the real paleoseismic record.

```GLTFM_Uncertainty##yrsForecast_out.txt```: Estimated GLTFM probability of an earthquake occuring in the next ## years for the current quiescent period  when the $R$, $Z_0$, shape parameter, and scale parameter estimates from ```GLTFM_UncertaintyParams_out.txt``` are applied to the real paleoseismic record.

```GLTFM_Uncertainty##yrsForecastToday_out.txt```: Estimated GLTFM probability of an earthquake occuring in the next ## years after the ```currentYr``` when the $R$, $Z_0$, shape parameter, and scale parameter estimates from ```GLTFM_UncertaintyParams_out.txt``` are applied to the real paleoseismic record.

```_UncertaintyParams_out.txt```: Renewal model shape parameter, and scale parameter estimates from each of the simulated catalogs in the uncertainty analysis.

```_UncertaintyPDF_out.txt```: Estimated renewal model inter-event time PDFs for the current quiescent period when the shape parameter and scale parameter estimates from ```_UncertaintyParams_out.txt``` are applied to the real paleoseismic record.

```_Uncertainty##yrsForecast_out.txt```: Estimated renewal model probability of an earthquake occuring in the next ## years for the current quiescent period when the renewal shape parameter and scale parameter estimates from ```_UncertaintyParams_out.txt``` are applied to the real paleoseismic record.

```_Uncertainty##yrsForecastToday_out.txt```: Estimated renewal model probability of an earthquake occuring in the next ## years after the ```currentYr``` when the renewal shape parameter and scale parameter estimates from ```_UncertaintyParams_out.txt``` are applied to the real paleoseismic record.

### Output Plots
```GLTFM_fitted.pdf```: A three-paneled plot showing (A) GLTFM's estimated hazard rate history for the paleoseimic history with best fitting shape and scale parameters and the $R$ value and residual strain for each earthquake. Panel (B) shows the estimated inter-event time probability density functions for the current quiescent period for different modelss. The mean, standard deviation, and loglikelihood are indicated. Panel (C) is the estimated ##-year forecast (the probability of an earthquake occuring in the next ## years) for the current quiescent period. The probability for ```currentYr``` is indicated. 

```Uncertainty##yrsForecastToday_KDE.pdf```: A kernel density plot of the probability of an earthquake occuring in the next ## years after the ```currentYr``` when the simulation shape parameter and scale parameter estimates (GLTFM and renewal) and $R$, and $Z_0$ (GLTFM) are applied to the real paleoseismic record.

```Uncertainty##yrsForecastToday_Histograms.pdf```: A histogram plot of the probability of an earthquake occuring in the next ## years after the ```currentYr``` when the simulation shape parameter and scale parameter estimates (GLTFM and renewal) and $R$, and $Z_0$ (GLTFM) are applied to the real paleoseismic record.

```_Uncertainty_Plots.pdf```: Estimated inter-event time probability density functions for the current quiescent period based on simulated catalogs (Panel A). Estimated probability of an earthquake occuring in the next ## years for the current quiescent period based on simulated catalogs (Panel B). Remaining panels -- distribution of parameters from simulated catalogs.



