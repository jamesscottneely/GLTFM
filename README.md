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

```T_Forecast```: Window of time (in years) of interest for earthquake probability forecasts.

```yrsForecast```:

```z0_flag```:

```R_flag```:

```incOpenInt```:

```model```:

```iterFit```:

```params```:

