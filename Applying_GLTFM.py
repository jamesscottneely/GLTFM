#!/usr/bin/env python3		

############### Import packages
import GLTFM
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import matplotlib.transforms as transforms
import pandas as pd
import sys

############################  USER DEFINED INPUTS
dataPath='ENTER PATH HERE' # Enter path to directory that contains input records
include_uncertainty = True # Set Boolean (True or False) if you want to include uncertainties
seed_val= 128231839 # Define a seed value for random number generator. Only used if include_uncertainty = TRUE
num_sims = 10 # Enter the number of simulations for uncertainty analysis. Only used if include_uncertainty = TRUE
############################  


#### Read in inter-event times
eq_inter = GLTFM.paleoRead(dataPath+'/PaleoRecord.txt')

#### Dictionary of renewal models to compare GLTFM to. 'name' is the scipy.stats.rv_continuous name
renewal_mod = {'Lognormal':{'name':'lognorm'},'BPT':{'name':'invgauss'},'Weibull':{'name':'weibull_min'}}

#### Read in the necessary data
paramDic = GLTFM.paramRead(dataPath)
q_period,eqtimes = GLTFM.eqTimes(eq_inter,paramDic['priorEQ'],paramDic['currentYr'])

#### Create array of time values for probability calculation
xVals = np.arange(paramDic['yrsForecast']+1)

##### Find the best fitting GLTFM parameters
R_out,z_out,P0_out,P1_out,loglikeGLTFM = GLTFM.fit_GLTFM(u=eq_inter,q_period=q_period,paramDic = paramDic)

#### Calculate GLTFM Probability Distribution for current quiescent period
rv_GLTFM = GLTFM.GLTFM_gen(z=z_out[-1],model=paramDic['model'],P0=P0_out,P1=P1_out)
GLTFM_PDF = rv_GLTFM.pdf(xVals) # GLTFM probability density function after most recent earthquake
GLTFM_CDF =  rv_GLTFM.cdf(xVals)  # GLTFM cumulative distribution function after most recent earthquake
GLTFM_Haz = rv_GLTFM.haz_z(xVals) # GLTFM hazard rate function after most recent earthquake
GLTFM_Tforecast = rv_GLTFM.Tforecast(xVals,paramDic['T_Forecast']) # GLTFM T-year forecast after most recent earthquake
GLTFM_Tforecast_today = rv_GLTFM.Tforecast(q_period,paramDic['T_Forecast'])# GLTFM T-year forecast after most recent earthquake for time of interest



### Apply standard probability models to data
for k in renewal_mod:
	#### Find best fitting renewal model parameters
	bestParam,loglike = GLTFM.fit_Renewal(distribution=renewal_mod[k]['name'],u=eq_inter,q_period=q_period,incOpenInt=paramDic['incOpenInt'])
	renewal_mod[k]['bestParams'] = bestParam
	renewal_mod[k]['loglike'] = loglike
	#### Calculate PDF, CDF, Hazard, and T-year forecast for renewal models
	rv_renew = getattr(stats,renewal_mod[k]['name'])(*bestParam) ### Generate random variable
	renewal_mod[k]['PDF'] = rv_renew.pdf(xVals)
	renewal_mod[k]['CDF'] = rv_renew.cdf(xVals)
	renewal_mod[k]['Haz'] = rv_renew.pdf(xVals)/rv_renew.sf(xVals)
	renewal_mod[k]['T_Forecast'] = (rv_renew.cdf(xVals+paramDic['T_Forecast'])-rv_renew.cdf(xVals))/rv_renew.sf(xVals)
	renewal_mod[k]['TForecast_Today'] = (rv_renew.cdf(q_period+paramDic['T_Forecast'])-rv_renew.cdf(q_period))/rv_renew.sf(q_period)
	renewal_mod[k]['mean'] = rv_renew.mean()
	renewal_mod[k]['std'] = rv_renew.std()
	

#### Save data to files
np.savetxt(dataPath+'/R_out.txt',R_out,fmt='%.4f') # Save R values to file. One for each earthquake
np.savetxt(dataPath+'/z_out.txt',z_out,fmt='%.4f') # Save Z (residual strain) values to file. One for each earthquake
pdf_param_out = [['GLTFM',P0_out,P1_out,loglikeGLTFM,GLTFM_Tforecast_today]]+ [[k,renewal_mod[k]['bestParams'][0],renewal_mod[k]['bestParams'][2],renewal_mod[k]['loglike'],renewal_mod[k]['TForecast_Today']] for k in renewal_mod] # List of parameter values
pd.DataFrame(pdf_param_out,columns=['Model','Shape','Scale','Loglikelihood',"{:.0f}yrsForecast".format(paramDic['T_Forecast'])]).to_csv(dataPath+'/model_out.txt',index=False) # Save best parameters for each model 
pdf_vals_out = [GLTFM_PDF]+ [renewal_mod[k]['PDF'] for k in renewal_mod] # List of PDFs
pd.DataFrame(pdf_vals_out).transpose().to_csv(dataPath+'/PDF_out.txt',index=False,header=['GLTFM']+[k for k in renewal_mod]) # Save PDF for each model 
cdf_vals_out = [GLTFM_CDF]+ [renewal_mod[k]['CDF'] for k in renewal_mod] # List of CDFs
pd.DataFrame(cdf_vals_out).transpose().to_csv(dataPath+'/CDF_out.txt',index=False,header=['GLTFM']+[k for k in renewal_mod]) # Save CDF for each model 
Tfore_vals_out = [GLTFM_Tforecast]+ [renewal_mod[k]['T_Forecast'] for k in renewal_mod] # List of T_forecasts
pd.DataFrame(Tfore_vals_out).transpose().to_csv(dataPath+"/{:.0f}yrsForecast_out.txt".format(paramDic['T_Forecast']),index=False,header=['GLTFM']+[k for k in renewal_mod]) # Save T_forecast for each model 


### Create a plot of GLTFM's hazard rate history (including residual strain), plot PDFs and the T-year forecasts
axd = plt.figure(figsize=(12,6)).subplot_mosaic(
"""
AAAA
BBCC
""")
#### Plot hazard rate curve for paleoseismic record
GLTFM.plotHaz(axd['A'],z_out,eqtimes,eq_inter,R_out,P0_out,P1_out,paramDic,q_period)
#### Plot PDFs
for k in renewal_mod:
	style_dic = GLTFM.plot_styles(renewal_mod[k]['name'])
	label_data = "({:.0f}".format(renewal_mod[k]['mean']) + r"$\pm$" + "{:.0f}".format(renewal_mod[k]['std']) + ", {:.1f})".format(renewal_mod[k]['loglike'])
	labelstr = style_dic['name'] +" " +label_data #
	axd['B'].plot(xVals,renewal_mod[k]['PDF'],color=style_dic['color'],label=labelstr,lw=1,ls=style_dic['ls'])
	label_data = "({:.2f})".format(renewal_mod[k]['TForecast_Today'])
	labelstr = style_dic['name'] +" " +label_data #
	axd['C'].plot(xVals,renewal_mod[k]['T_Forecast'],color=style_dic['color'],label=labelstr,lw=1)

style_dic = GLTFM.plot_styles('GLTFM')
label_data = "({:.0f}".format(rv_GLTFM.mean()) + r"$\pm$" + "{:.0f}".format(rv_GLTFM.std()) + ", {:.1f})".format(loglikeGLTFM)
labelstr = 'GLTFM' +" " +label_data #
axd['B'].plot(xVals,renewal_mod[k]['PDF'],color=style_dic['color'],label=labelstr,lw=1,ls=style_dic['ls'])

label_data = "({:.2f})".format(GLTFM_Tforecast_today)
labelstr = 'GLTFM' +" " +label_data 
axd['C'].plot(xVals,GLTFM_Tforecast,color=style_dic['color'],label=labelstr,lw=1,ls=style_dic['ls'])


titlestr = 'Estimated Hazard Rate History'
axd['A'].set_title(titlestr)

axd['B'].set_title('Estimated Current Inter-Event Time PDF')
axd['B'].set_xlabel('Years Since Most Recent Earthquake')
axd['B'].set_ylabel('Probability Density')
axd['B'].set_ylim(ymin=0)
axd['B'].set_xlim([0,paramDic['yrsForecast']])

tstr = 'Estimated ' + "{:.0f}".format(paramDic['T_Forecast']) + '-Year Forecast'
axd['C'].set_title(tstr)
axd['C'].set_xlabel('Years Since Most Recent Earthquake')
axd['C'].set_ylabel('Probability')
axd['C'].set_ylim(ymin=0)
axd['C'].set_xlim([0,paramDic['yrsForecast']])

#### Annotate current year
axd['B'].axvline(q_period,color='grey',ls='--')
trans = transforms.blended_transform_factory(
    axd['B'].transData, axd['B'].transAxes)
axd['B'].text(q_period, 0.01, paramDic['currentYr'],va ='bottom', rotation = 270,transform=trans)
trans = transforms.blended_transform_factory(
    axd['C'].transData, axd['C'].transAxes)
axd['C'].text(q_period, 0.01, paramDic['currentYr'], va ='bottom',rotation = 270,transform=trans)

lx = 0
ly=1.05
axd['A'].text(lx,ly,'(A)',ha='center',va='bottom',transform=axd['A'].transAxes)
axd['B'].text(lx,ly,'(B)',ha='center',va='bottom',transform=axd['B'].transAxes)
axd['C'].text(lx,ly,'(C)',ha='center',va='bottom',transform=axd['C'].transAxes)
	

axd['C'].axvline(q_period,color='grey',ls='--')
axd['B'].legend(title=r'Model ($\mu\pm\sigma$,LL)')
axd['C'].legend(title='Model (Pr. Today)')
plotFile = dataPath + '/GLTFM_fitted.pdf'
plt.tight_layout()
plt.savefig(plotFile)
plt.close()

if not include_uncertainty: # Stop here if uncertainty analysis unwanted
	sys.exit()

###### ###### ###### ###### ###### ###### ###### ###### ###### 
###### Uncertainty testing
###### ###### ###### ###### ###### ###### ###### ###### ###### 
sub_len = (paramDic['currentYr'] - eqtimes[0]) # Set the length of synthetic catalogs. Here we assume
rng = np.random.default_rng(seed=seed_val) # Set the seed value for the random number generator


#### Estimate GLTFM uncertainty
param_array,PDF_array,TForecast_array,TForecast_Today_array = GLTFM.GLTFM_uncertainty(num_sims,paramDic,P0_out,P1_out,R_out,eq_inter,q_period,sub_len,rng) #### Estimate GLTFM uncertainty

#### Save GLTFM Uncertainty data
np.savetxt(dataPath+'/GLTFM_UncertaintyParams_out.txt',param_array,header='R z0 Shape Scale',comments='',fmt='%.8f') # Save R,Z0,Shape and Scale parameters
np.savetxt(dataPath+'/GLTFM_UncertaintyPDF_out.txt',PDF_array,fmt='%.8f') # Save simulated PDFs
np.savetxt(dataPath+"/GLTFM_Uncertainty{:.0f}yrsForecast_out.txt".format(paramDic['T_Forecast']),TForecast_array,fmt='%.8f') # Save simulated PDFs
np.savetxt(dataPath+"/GLTFM_Uncertainty{:.0f}yrsForecastToday_out.txt".format(paramDic['T_Forecast']),TForecast_Today_array,fmt='%.8f') # Save simulated PDFs


## Estimate renewal model uncertainty and save data
for k in renewal_mod:
	# Renewal models
	param_array_Renewal,PDF_array_Renewal,TForecast_array_Renewal,TForecast_Today_array_Renewal = GLTFM.Renewal_uncertainty(num_sims,paramDic,renewal_mod[k]['bestParams'],renewal_mod[k]['name'],q_period,sub_len,rng)
	np.savetxt(dataPath+'/'+k+'_UncertaintyParams_out.txt',param_array_Renewal,header='Shape Loc Scale',comments='',fmt='%.8f') # Save Shape and Scale parameters
	np.savetxt(dataPath+'/'+k+'_UncertaintyPDF_out.txt',PDF_array_Renewal,fmt='%.8f') # Save simulated PDFs
	np.savetxt(dataPath+'/'+k+"_Uncertainty{:.0f}yrsForecast_out.txt".format(paramDic['T_Forecast']),TForecast_array_Renewal,fmt='%.8f') # Save simulated PDFs
	np.savetxt(dataPath+'/'+k+"_Uncertainty{:.0f}yrsForecastToday_out.txt".format(paramDic['T_Forecast']),TForecast_Today_array_Renewal,fmt='%.8f') # Save simulated PDFs


######### Generate Figures to show uncertainty in forecasts. The following plotting codes
#### have been written so that they can be run independently of the simulations runs above
#### (i.e. you can run the simulations once and then modify the visualization codes as need be in a separate script)



#### Plot Kernel Density Estimates of EQ Probability in next N years based on simulations
axd = plt.figure(figsize=(4,4)).subplot_mosaic(
"""
A
""")
### Plot renewal models
for k in renewal_mod:
	filename = dataPath+'/'+k+"_Uncertainty{:.0f}yrsForecastToday_out.txt".format(paramDic['T_Forecast'])
	data = np.loadtxt(filename) # load data
	style_dic = GLTFM.plot_styles(renewal_mod[k]['name'])
	density = stats.gaussian_kde(data)
	xs = np.linspace(0,1,1000)
	axd['A'].plot(xs,density(xs),color=style_dic['color'],label = k,ls=style_dic['ls'])
	
### Plot GLTFM model
filename = dataPath+"/GLTFM_Uncertainty{:.0f}yrsForecastToday_out.txt".format(paramDic['T_Forecast'])
data = np.loadtxt(filename)
style_dic = GLTFM.plot_styles('GLTFM')
density = stats.gaussian_kde(data)
xs = np.linspace(0,1,1000)
axd['A'].plot(xs,density(xs),color=style_dic['color'],label = 'GLTFM')

### Annotate
axd['A'].set_xlabel('Estimated ' + "{:.0f}".format(paramDic['T_Forecast']) + '-Year Forecast')
axd['A'].set_title('Forecast Uncertainty')
axd['A'].set_xlim([0,1])	
axd['A'].set_ylim(ymin=0)
axd['A'].legend()


### save kernel density plot
plotFile = dataPath+"/Uncertainty{:.0f}yrsForecastToday_KDE.pdf".format(paramDic['T_Forecast'])
plt.tight_layout()
plt.savefig(plotFile)
plt.close()

#### Individual histograms of EQ Probability in next N years based on simulations. Same data
#### as KDE plot but now shown in histogram form
axd = plt.figure(figsize=(8,8)).subplot_mosaic(
"""
AB
CD
""")
ax_l = ['A','B','C','D']
idx = 0
for k in renewal_mod:
	filename = dataPath+'/'+k+"_Uncertainty{:.0f}yrsForecastToday_out.txt".format(paramDic['T_Forecast'])
	data = np.loadtxt(filename)
	style_dic = GLTFM.plot_styles(renewal_mod[k]['name'])
	axd[ax_l[idx]].hist(data,bins=100,range=(0,1),color=style_dic['color'])
	axd[ax_l[idx]].set_xlim([0,1])
	axd[ax_l[idx]].set_title(k)
	axd[ax_l[idx]].set_ylabel("Number of Simulations")
	axd[ax_l[idx]].set_xlabel('Estimated ' + "{:.0f}".format(paramDic['T_Forecast']) + '-Year Forecast')

	idx +=1 
filename = dataPath+"/GLTFM_Uncertainty{:.0f}yrsForecastToday_out.txt".format(paramDic['T_Forecast'])
data = np.loadtxt(filename)
style_dic = GLTFM.plot_styles('GLTFM')
axd[ax_l[idx]].hist(data,bins=100,range=(0,1),color=style_dic['color'])
axd[ax_l[idx]].set_ylabel("Number of Simulations")
axd[ax_l[idx]].set_xlabel('Estimated ' + "{:.0f}".format(paramDic['T_Forecast']) + '-Year Forecast')
axd[ax_l[idx]].set_xlim([0,1])
axd[ax_l[idx]].set_title('GLTFM')
lx = -.01
ly = 1.01
axd['A'].text(lx,ly,'(A)',ha='center',va='bottom',transform=axd['A'].transAxes)
axd['B'].text(lx,ly,'(B)',ha='center',va='bottom',transform=axd['B'].transAxes)
axd['C'].text(lx,ly,'(C)',ha='center',va='bottom',transform=axd['C'].transAxes)
axd['D'].text(lx,ly,'(D)',ha='center',va='bottom',transform=axd['D'].transAxes)
### Save histograms
plotFile = dataPath+"/Uncertainty{:.0f}yrsForecastToday_Histograms.pdf".format(paramDic['T_Forecast'])
plt.tight_layout()
plt.savefig(plotFile)
plt.close()


#### PLot Uncertainty of current PDFs, N-year forecasts, and PDF paramater uncertainties
for k in renewal_mod:
	### Read in parameters
	paramfile = dataPath+'/'+k+"_UncertaintyParams_out.txt"
	dataParam = pd.read_csv(paramfile,header=0,sep=" ")
	#### Read in PDFs
	PDFfile = dataPath+'/'+k+"_UncertaintyPDF_out.txt"
	dataPDF = np.loadtxt(PDFfile)	
	#### Read in N-years forecast
	filename = dataPath+'/'+k+"_Uncertainty{:.0f}yrsForecast_out.txt".format(paramDic['T_Forecast'])
	data =  np.loadtxt(filename)
	#### Create Plots
	axd = plt.figure(figsize=(8,8)).subplot_mosaic(
		"""
		AB
		CD
		""")
	#### Plot PDF uncertainty
	alpha = .1
	style_dic = GLTFM.plot_styles(renewal_mod[k]['name'])
	axd['A'].plot(np.tile(xVals,(np.shape(dataPDF)[1],1)).T,dataPDF,color=style_dic['color'],alpha=alpha)

	axd['B'].plot(np.tile(xVals,(np.shape(data)[1],1)).T,data,color=style_dic['color'],alpha=alpha)

	axd['A'].set_title('Estimated Current Inter-Event Time PDF')
	axd['A'].set_xlabel('Years Since Most Recent Earthquake')
	axd['A'].set_ylabel('Probability Density')
	axd['A'].set_ylim(ymin=0)
	axd['A'].set_xlim([0,paramDic['yrsForecast']])
	tstr = 'Estimated ' + "{:.0f}".format(paramDic['T_Forecast']) + '-Year Forecast'
	axd['B'].set_title(tstr)
	axd['B'].set_xlabel('Years Since Most Recent Earthquake')
	axd['B'].set_ylabel('Probability')
	axd['B'].set_ylim(ymin=0)
	axd['B'].set_xlim([0,paramDic['yrsForecast']])

	#### Annotate current year
	axd['A'].axvline(q_period,color='grey',ls='--')
	trans = transforms.blended_transform_factory(
	    axd['A'].transData, axd['A'].transAxes)
	axd['A'].text(q_period, 0.01, paramDic['currentYr'],va ='bottom', rotation = 270,transform=trans)
	axd['B'].axvline(q_period,color='grey',ls='--')
	trans = transforms.blended_transform_factory(
 	   axd['B'].transData, axd['B'].transAxes)
	axd['B'].text(q_period, 0.01, paramDic['currentYr'], va ='bottom',rotation = 270,transform=trans)
	
	#### Add the parameter histograms
	style_dic = GLTFM.plot_styles(renewal_mod[k]['name'])
	axd['C'].hist(dataParam['Shape'],color=style_dic['color'])
	axd['C'].set_xlabel("Shape Parameter")
	axd['C'].set_ylabel('Number of Simulations')
	axd['D'].hist(dataParam['Scale'],color=style_dic['color'])
	axd['D'].set_xlabel("Scale Parameter")
	axd['D'].set_ylabel('Number of Simulations')

	### Save Uncertainty plots
	plotFile = dataPath+'/'+k+"_Uncertainty_Plots.pdf"
	plt.tight_layout()
	plt.savefig(plotFile)
	plt.close()

### Plot the GLTFM Forecast and parameter uncertainties
 ### Read in parameters
paramfile = dataPath+"/GLTFM_UncertaintyParams_out.txt"
dataParam = pd.read_csv(paramfile,header=0,sep=" ")
#### Read in PDFs
PDFfile = dataPath+"/GLTFM_UncertaintyPDF_out.txt"
dataPDF = np.loadtxt(PDFfile)	
#### Read in N-years forecast
filename = dataPath+"/GLTFM_Uncertainty{:.0f}yrsForecast_out.txt".format(paramDic['T_Forecast'])
data =  np.loadtxt(filename)
#### Create Plots
axd = plt.figure(figsize=(8,12)).subplot_mosaic(
	"""
	AB
	CD
	EF
	""")
#### Plot PDF uncertainty
alpha = .1
style_dic = GLTFM.plot_styles('GLTFM')
axd['A'].plot(np.tile(xVals,(np.shape(dataPDF)[1],1)).T,dataPDF,color=style_dic['color'],alpha=alpha)
axd['B'].plot(np.tile(xVals,(np.shape(data)[1],1)).T,data,color=style_dic['color'],alpha=alpha)
axd['A'].set_title('Estimated Current Inter-Event Time PDF')
axd['A'].set_xlabel('Years Since Most Recent Earthquake')
axd['A'].set_ylabel('Probability Density')
axd['A'].set_ylim(ymin=0)
axd['A'].set_xlim([0,paramDic['yrsForecast']])
tstr = 'Estimated ' + "{:.0f}".format(paramDic['T_Forecast']) + '-Year Forecast'
axd['B'].set_title(tstr)
axd['B'].set_xlabel('Years Since Most Recent Earthquake')
axd['B'].set_ylabel('Probability')
axd['B'].set_ylim(ymin=0)
axd['B'].set_xlim([0,paramDic['yrsForecast']])
#### Annotate current year
axd['A'].axvline(q_period,color='grey',ls='--')
trans = transforms.blended_transform_factory(
    axd['A'].transData, axd['A'].transAxes)
axd['A'].text(q_period, 0.01, paramDic['currentYr'],va ='bottom', rotation = 270,transform=trans)
axd['B'].axvline(q_period,color='grey',ls='--')
trans = transforms.blended_transform_factory(
   axd['B'].transData, axd['B'].transAxes)
axd['B'].text(q_period, 0.01, paramDic['currentYr'], va ='bottom',rotation = 270,transform=trans)
#### Add the parameter histograms
axd['C'].hist(dataParam['Shape'],color=style_dic['color'])
axd['C'].set_xlabel("Shape Parameter")
axd['C'].set_ylabel('Number of Simulations')
axd['D'].hist(dataParam['Scale'],color=style_dic['color'])
axd['D'].set_xlabel("Scale Parameter")
axd['D'].set_ylabel('Number of Simulations')
axd['E'].hist(dataParam['R'],color=style_dic['color'])
axd['E'].set_xlabel(r'$R$')
axd['E'].set_ylabel('Number of Simulations')
axd['F'].hist(dataParam['z0'],color=style_dic['color'])
axd['F'].set_xlabel(r'$Z_0$')
axd['F'].set_ylabel('Number of Simulations')


### Save Uncertainty plots
plotFile = dataPath+"/GLTFM_Uncertainty_Plots.pdf"
plt.tight_layout()
plt.savefig(plotFile)
plt.close()




