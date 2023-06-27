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

```R_flag```: Flag for the strain drop $R$ (in years) after an earthquake in the paleoseismic record. This flag can be set to a specific value (e.g. ```150```), estimated via maximum likelihood by setting to ```Free```, or estimated using tectonic information. To estimate R using

```incOpenInt```:

```model```:

```iterFit```:

```params```:

