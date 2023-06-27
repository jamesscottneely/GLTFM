# GLTFM
Generalized Long Term Fault Memory Model. 

## Introduction
This repository contains the code for the Generalized Long Term Fault Memory (GLTFM) model described in Neely et al. *Submitted*. GLTFM allows users to calculate large earthquake probabilities that incorporate the impact of residual strain after an earthquake. Users input a paleoseismic record with the inter-event times of past large earthquakes and GLTFM calculates earthquake probabilities for the current quiescent period.

## Installation
TO BE COMPLETED

## Input Files

*PaleoRecord.txt*: [Required] This file contains one column with the observed earthquake inter-event times (in years) in ascending order. For example if earthquakes occured in 1900, 1950, 1962, and 2000, this file would contain the following inter-event times:
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

```model```: This is the base probability to use for GLTFM. The model name must correspond to one of the 2-parameter continuous probability models in ```scipy.stats```. As discussed in the manuscript we prefer using the Weibull model which in ```scipy.stats``` is ```weibull_min```.

```iterFit```: This is a Y or N flag to indicate whether the free parameters should be iteratively solved for (as described in Neely et al. *submitted*) or if they should be solved for all at once.

```params```: This list is for the shape and scale parameters of GLTFM's underlying probability distribution. Either a specific value should be specified or set to ```Free``` to estimate via maximum likelihood. For example ```[Free Free]``` or ```[2.0 Free]```. See ```scipy.stats``` for a description of the shape and scale parameters for the different probability distributions.



