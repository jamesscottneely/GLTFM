#!/usr/bin/env python3		

############### Import packages
import GLTFM
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as stats
import matplotlib.transforms as transforms



def plot_styles(model):
	#####
	# Function to retrieve plotting formats 
	# Input:
	# model: name of probability model
	# Output:
	# styles: dictionary for the different models
	#####
	styles = {'GLTFM':{'color':"#ffa600",'ls':'-','name':'GLTFM'},'expon':{'color':'#003f5c','ls':'-','name':'Exponential'},'lognorm':{'color':'#7a5195','ls':'--','name':'Lognormal'},'invgauss':{'color':'#ef5675','ls':'-','name':'BPT'},'weibull_min':{'color':'#21a2de','ls':'-','name':'Weibull'}} # line styles
	return styles[model]

def plotHaz(axd,z_out,eqtimes,eq_inter,R_out,P0_out,P1_out,paramDic,q_period):
	#####
	# Function to plot the hazard rate history of the paleoseismic record
	# Input: 
	# axd: matplotlib axis for plotting
	# z_out: array of z_i values
	# eqtimes: array of earthquake times
	# eq_inter: array of earthquake inter-event times
	# R_out: array of R values for earthquake drop
	# P0_out: probability shape parameter value
	# P1_out: probability scale parameter value
	# paramDic: dictionary of fitting parameters
	# q_period: current quiescent period
	# Output: plot of hazard rate history with residual strain indicated after each earthquake
	#####
	style_dic = plot_styles('GLTFM')
	timeVals = np.arange(max(eq_inter)+1)# Time increments for forecast
	rv_GLTFM_Haz = GLTFM.GLTFM_gen(z=z_out.reshape(len(z_out),1),model=paramDic['model'],P0=P0_out,P1=P1_out)
	Haz_hist = rv_GLTFM_Haz.haz_z(np.tile(timeVals,(len(z_out),1))) # GLTFM hazard rate function after each  earthquake
	### Plot hazard curve
	maxHaz = [] # List for tracking max hazard of each inter-event time
	for eq_idx in range(len(eqtimes)-1):
		axd.plot(eqtimes[eq_idx]+timeVals[:int(eq_inter[eq_idx])+1],Haz_hist[eq_idx,0:int(eq_inter[eq_idx])+1],color=style_dic['color'],ls=':') # Plot hazard curve
		maxHaz.append(Haz_hist[eq_idx,int(eq_inter[eq_idx])])
		Ztextstr =  "{:.0f}".format(z_out[eq_idx])
		Rtextstr = "{:.0f}".format(R_out[eq_idx])
		axd.text(eqtimes[eq_idx],Haz_hist[eq_idx,0],Ztextstr,va='top',ha='center',color='black',fontsize=8) #plot Zi string
		axd.text(eqtimes[eq_idx],maxHaz[eq_idx-1],Rtextstr,va='top',ha='left',color='black',fontsize=8,rotation=270) #plot R string

	### Annotate final earthquake
	Ztextstr =  "{:.0f}".format(z_out[-1])
	Rtextstr = "{:.0f}".format(R_out[-1])
	axd.text(eqtimes[-1],Haz_hist[-1,0],Ztextstr,va='top',ha='center',color='black',fontsize=8) #plot Zi string
	axd.text(eqtimes[-1],maxHaz[-1],Rtextstr,va='top',ha='left',color='black',fontsize=8,rotation=270) #plot R string
	axd.text(eqtimes[-1],(maxHaz[-1]-Haz_hist[-1,0])/2,paramDic['priorEQ'],va='center',fontweight='bold')

	### Plot hazard curve
	axd.plot(eqtimes[-1]+timeVals[:int(q_period)+1],Haz_hist[-1,:int(q_period)+1],color=style_dic['color'],ls=':') # Plot current quiescence period
	### Plot earthquake times
	axd.vlines(eqtimes[0],Haz_hist[0,0],maxHaz[0],color=style_dic['color']) # Plot first earthquake
	axd.vlines(eqtimes[1:],Haz_hist[1:,0],maxHaz,color=style_dic['color']) # Plot rest of earthquakes

	ptext = "Base Model: " + plot_styles(paramDic['model'])['name'] + " "+r"Shape: {:.1f} ".format(P0_out) + r"Scale: {:.1f}".format(P1_out) + "\n" + "Paleoseismic Record Mean: {:.0f}, SD: {:.0f}, CV: {:.1f} ".format(np.mean(eq_inter),np.std(eq_inter),np.std(eq_inter)/np.mean(eq_inter)) 
	axd.text(.5,.9,ptext,va='top',ha='center',color='black',transform=axd.transAxes) #include resdiual text

	axd.set_xlabel('Year')
	axd.set_ylabel('Hazard Rate')
	trans = transforms.blended_transform_factory(
    	axd.transData, axd.transAxes)
	axd.text(paramDic['currentYr'], 0.01, paramDic['currentYr'],va ='bottom', rotation = 270,transform=trans)
	axd.axvline(q_period+max(eqtimes),color='grey',ls='--')
	return axd

