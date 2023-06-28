############### Import packages
import numpy as np
import scipy.stats as stats
import GLTFM
from scipy.integrate import quad
import copy

def sim_GLTFM(model,P0,P1,R,sub_len,rng):
	#### 
	# Function to simulate GLTFM catalogs
	####

	### Generat synthetic catalog
	if isinstance(R, np.ndarray):
		if len(np.unique(R)) > 1:
			print ('Warning: Simulations require only 1 R value. Largest R selected for simulation.')
		R = max(R)
	rv = getattr(stats,model)(P0,0,scale=P1)
	u = rv.rvs(size=1,random_state=rng)[0] # Generate random inter-event before first earthquake
	z0 = max(0,u-R) # Corresponding initial start strain
	z_list = [z0] # Initialize residual strain list
	t_list = [0] # Initialize EQ time list
	u_list = [] # Initialize EQ inter-event times list
	t = 0  # Initialize time variable
	z = z0 # Initialize z
	while t < sub_len: # Generate synthetic catalog 
		rand_GLTFM = GLTFM.GLTFM_gen_rand(z=z,P0=P0,P1=P1,model=model)
		u = rand_GLTFM.rvs(size=1,random_state=rng)[0]
		z = max(0,u+z-R) # Calculate residual strain after next earthquake
		u_list.append(u) # Append values
		z_list.append(z) # Append values
		t = t+u # Track time
		t_list.append(t) # Append values
	q_period_sim = sub_len - t_list[-2] # Simulated current quiescent period
	return np.array(z_list[:-1]),np.array(t_list[:-1]),np.array(u_list[:-1]),q_period_sim,z_list[-1],t_list[-1],u_list[-1]
	
def sim_GLTFM_test(u_obs,q_period,u_next,paramDic):		
	#### 
	# Function to test how well model performs on synthetic records
	####
	## Fit GLTFM
	R_out,z_out,P0_out,P1_out,loglikeGLTFM = GLTFM.fit_GLTFM(u=u_obs,q_period=q_period,paramDic = paramDic)
	## Find expected GLTFM
	rand_GLTFM = GLTFM.GLTFM_gen_rand(z=z_out[-1],P0=P0_out,P1=P1_out,model=paramDic['model'])
	expected_GLTFM_u = rand_GLTFM.mean()	
	expectedSD_GLTFM_u = rand_GLTFM.std()	
	### Calculate scores
	ig_score = 	-1*np.log(rand_GLTFM.pdf(u_next))
	mom_score = -1*(((u_next-expected_GLTFM_u)/expectedSD_GLTFM_u)**2)-np.log(expectedSD_GLTFM_u**2)
	propLin_score = quad(rand_GLTFM.pdfsq,0,np.inf)[0] - 2*rand_GLTFM.pdf(u_next)
	return ig_score,mom_score,propLin_score	
	
def GLTFM_uncertainty(num_sims,paramDic,P0,P1,R,u_real,q_period_real,sub_len,rng):
	#### 
	# Function to estimate uncertainty of GLTFM forecasts
	# Inputs: 
	# num_sims: number of simulations to run
	# paramDic: input parameter dictionary
	# P0: best estimate for the shape parameter
	# P1: best estimate for the scale parameter
	# R: value for R - earthquake drop
	# u_real: paleoseismic record of interest
	# q_period_real: current quiescent period for paleoseismic record of interest
	# sub_len: length of subcatalogs to generate - usually set to length of paleoseismic record of interest
	# rng: random number generator
	# Outputs: 
	# param_array: array of the best fitting R, Z0, Shape, and Scale parameters from simulation. Note that if any are fixed initially, they will be fixed here. Also Z0 is estimated from the paleoseismic record of interest not the simulated record
	# PDF_array: array of estimated inter-event time PDFs for the current quiescent period when simulated paramaters applied to paleoseismic record of interest
	# TForecast_array: array of estimated earthquake probabilities in the next ## years for the current quiescent period when simulated paramaters applied to paleoseismic record of interest
	# TForecast_Today_array: array of estimated earthquake probabilities in the next ## years for current moment in time when simulated paramaters applied to paleoseismic record of interest
	####
	xVals = np.arange(0,paramDic['yrsForecast']+1)
	param_array = np.zeros((num_sims,4)) # Initialize paramater array for subcatalogs
	PDF_array = np.zeros((len(xVals),num_sims))
	TForecast_array = np.zeros((len(xVals),num_sims))
	TForecast_Today_array = np.zeros((num_sims))
	idx = 0
	for n in range(num_sims):
		z,t,u,q_period_sim,z_next,t_next,u_next = sim_GLTFM(paramDic['model'],P0,P1,R,sub_len,rng) # Simulate catalog
		R_out,z_out,P0_out,P1_out,loglikeGLTFM = GLTFM.fit_GLTFM(u=u,q_period=q_period_sim,paramDic = paramDic) # Estimate parameters from simulated catalog
		tempDic = copy.deepcopy(paramDic)  # Create a temporary dictionary for the loop
		tempDic['iterFit'] = 'N' # All parameters are accounted for except z0 so iterative approach no longer needed in simulation.
		tempDic['modelParams'] = [P0_out,P1_out] # Set model params for calculation
		tempDic['R'] = R_out[0] # Set R to simulated value
		R_out_fit,z_out_fit,P0_out_fit,P1_out_fit,loglikeGLTFM_fit = GLTFM.fit_GLTFM(u=u_real,q_period=q_period_real,paramDic = tempDic)# Use fitted parameters from simulate catalog to estimate Zi's
		rv_GLTFM_temp =GLTFM.GLTFM_gen(z=z_out_fit[-1],model=tempDic['model'],P0=P0_out,P1=P1_out) # Construct GLTFM from ranom variable
		GLTFM_PDF_temp = rv_GLTFM_temp.pdf(xVals) # GLTFM probability density function after most recent earthquake
		GLTFM_Tforecast_temp = rv_GLTFM_temp.Tforecast(xVals,tempDic['T_Forecast']) # GLTFM T-year forecast after most recent earthquake
		PDF_array[:,idx] = GLTFM_PDF_temp
		TForecast_array[:,idx] = GLTFM_Tforecast_temp	
		TForecast_Today_array[idx] = rv_GLTFM_temp.Tforecast(q_period_real,tempDic['T_Forecast']) 
		param_array[idx,:] = [R_out[0],z_out_fit[0],P0_out,P1_out]
		idx+=1	
	return param_array,PDF_array,TForecast_array,TForecast_Today_array

def sim_Renewal(model,params,sub_len,rng):
	#### 
	# Function to simulate renewal catalogs
	####

	### Generate synthetic catalog
	rv = getattr(stats,model)(*params)
	t_list = [0] # Initialize EQ time list
	u_list = [] # Initialize EQ inter-event times list
	t = 0  # Initialize time variable
	while t < sub_len: # Generate synthetic catalog 
		u = rv.rvs(size=1,random_state=rng)[0]
		u_list.append(u) # Append values
		t = t+u # Track time
		t_list.append(t) # Append values
	q_period_sim = sub_len - t_list[-2] # Simulated current quiescent period
	return np.array(t_list[:-1]),np.array(u_list[:-1]),q_period_sim,t_list[-1],u_list[-1]
	

	
def sim_Renewal_test(u_obs,q_period,u_next,model,paramDic):		
	#### 
	# Function to test how well model performs on synthetic records
	####	
	#### Fit renewal model to data
	bestParams,loglike = GLTFM.fit_Renewal(distribution=model,u=u_obs,q_period=q_period,incOpenInt=paramDic['incOpenInt'])
	#### Create random variable for best fitting params
	rv_renew = getattr(stats,model)(*bestParams) ### Generate random variable
	rv_sq = renew_sq(bestParams[0],bestParams[2],model) # Generat for sqauring
	### Calculate expected time until next earthquake
	expected_Renew_u = rv_renew.mean()	
	expectedSD_Renew_u = rv_renew.std()	
	### Calculate scores
	ig_score = 	-1*np.log(rv_renew.pdf(u_next))
	mom_score = -1*(((u_next-expected_Renew_u)/expectedSD_Renew_u)**2)-np.log(expectedSD_Renew_u**2)
	propLin_score = quad(rv_sq.pdfsq_f,0,np.inf)[0] - 2*rv_renew.pdf(u_next)
	return ig_score,mom_score,propLin_score

def Renewal_uncertainty(num_sims,paramDic,params,model,q_period_real,sub_len,rng):
	#### 
	# Function to estimate uncertainty of GLTFM forecasts
	####
	xVals = np.arange(0,paramDic['yrsForecast']+1)
	if model == 'expon':
		param_array = np.zeros((num_sims,2)) # Initialize paramater array for subcatalogs
	else:
		param_array = np.zeros((num_sims,3)) # Initialize paramater array for subcatalogs
	PDF_array = np.zeros((len(xVals),num_sims))
	TForecast_array = np.zeros((len(xVals),num_sims))
	TForecast_Today_array = np.zeros((num_sims))
	idx = 0
	for n in range(num_sims):
		t,u,q_period_sim,t_next,u_next = sim_Renewal(model,params,sub_len,rng)# Simulate catalog
		bestParams,loglike = GLTFM.fit_Renewal(distribution=model,u=u,q_period=q_period_sim,incOpenInt=paramDic['incOpenInt']) # Estimate parameters from simulated catalog
		rv_renew = getattr(stats,model)(*bestParams) ### Generate random variable
		Renewal_PDF_temp = rv_renew.pdf(xVals) # renewal probability density function after most recent earthquake
		Renewal_Tforecast_temp = (rv_renew.cdf(xVals+paramDic['T_Forecast'])-rv_renew.cdf(xVals))/rv_renew.sf(xVals)
		PDF_array[:,idx] = Renewal_PDF_temp
		TForecast_array[:,idx] = Renewal_Tforecast_temp	
		TForecast_Today_array[idx] = (rv_renew.cdf(q_period_real+paramDic['T_Forecast'])-rv_renew.cdf(q_period_real))/rv_renew.sf(q_period_real)
		param_array[idx,:] = bestParams
		idx+=1	
	return param_array,PDF_array,TForecast_array,TForecast_Today_array	


class renew_sq(stats.rv_continuous):
	def __init__(self,P0,P1,model,*args,**kwargs):
		super().__init__(*args, **kwargs)
        # init function
		self.P0 = P0
		self.P1 = P1
		self.model = model	
		
	def pdfsq_f(self, u):
		rv_mod = getattr(stats,self.model)(self.P0,0,scale=self.P1) # Define base model
		PDFsq =  rv_mod.pdf(u)**2
		return PDFsq

	
