import numpy as np
import scipy.optimize as opt
import math
import scipy.stats as stats
import scipy.interpolate as interpolate
import sys


def eqTimes(u,priorEQ,currentYr):
	#####
	# Function converts eq inter-event data to specific dates based on when
	# most recent prior earthquake occured (priorEQ) and the current year.
	# Returns current quiescent time (q_period) and earthquake dates
	#####
	q_period = currentYr - priorEQ # Calculate current quiescence time
	eqtimes = np.flip(currentYr-np.cumsum(np.flip(np.concatenate((u,[q_period]))))) 
	return q_period,eqtimes

def fit_GLTFM(u,q_period,paramDic):
	##### 
	# Wrapper function to find the best fitting GLTFM parameters
	# Inputs:
	# u: Array of inter-event times
	# q_period: length of current quiescent period
	# paramDic: Dictionary of input fitting paramaters
	# Outputs:
	# R_out: Array of R values (earthquake drop size) for each earthquake
	# z_out: Array of Z values (residual strain) after each earthquake
	# P0_out: Best fitting paramater (0) of GLTFM PDF
	# P1_outL Best fitting paramater (1) of GLTFM PDF
	# -1*res.fun: Loglikelihood of fitted GLTFM PDF for observed inter-events and quiescence
	#####	
	#### Generate initial guesses for the parameters of interest - only used if that parameter is free
	z0_init = 0 # Initial z0 guess
	Zbounds = [0,2*np.max(u)] # Bounds for inversion to find Z0	
	R_init = np.mean(u)*1.5  # Initial R drop estimate. Only used if R not set based on tectonic information
	Rbounds = [0,np.max(u)] # Bounds for R. Only used if R not set based on tectonic information
	modelParams0_init, loc1, modelParams1_init = getattr(stats,paramDic['model']).fit(u,floc=0) # Initial guess for PDF parameters
	P1bounds = [0,np.inf] # Bounds for PDF scale parameter

	#### Inversion related parameters
	method = 'L-BFGS-B' # inversion method
	options = {'maxiter':1000} # set number of iterations for inversion
	## Inversion specific initialization parameters based on baseline PDF for GLTFM. Recommend using weibull_min
	if paramDic['model'] == 'weibull_min':
		P0bounds = [1,10] # Min=1 so hazard rate curve never decreases with time
	else:
		P0bounds = [0,10] # All other PDFs can have shape paramater as low as 0
	#### Handling the R value 
	if paramDic['R_flag'] =='Variable': # If setting different R for each earthquakes
		try:
			paramDic['R_flag'] = GLTFM.RRead(paramDic['dataPath']+'/PaleoRecord_R.txt') # Read in R values
		except:
			print("ParamFile indicated that R_flag set to Variable but PaleoRecord_R.txt containing R's for each earthquake not found.")
			sys.exit(1) 
	### One step
	if paramDic['iterFit'] == 'N':
		fit_mask = [paramDic['R_flag'] =='Free',paramDic['z0_flag'] =='Free',paramDic['modelParams'][0]=='Free',paramDic['modelParams'][1]=='Free'] # Create mask for params to fit
		init_vals = np.array([R_init,z0_init,modelParams0_init,modelParams1_init]) # Array of initial guesses
		params = init_vals[fit_mask] # These are the paramaters that are being fit
		boundsInit = np.array([Rbounds,Zbounds,P0bounds,P1bounds]) # Bounds used for fitting
		bounds = boundsInit[fit_mask] # Select the appropriate bounds
		#### Minimize function and find best params
		res = opt.minimize(neg_logLike_GLTFM, params, args=(paramDic['model'],u,paramDic['incOpenInt'],q_period,paramDic['R_flag'],paramDic['z0_flag'],paramDic['modelParams'][0],paramDic['modelParams'][1],fit_mask), method = method,bounds=bounds,options=options)
	#### Itertatively
	elif paramDic['iterFit'] == 'Y':
		bestParam,loglike = fit_Renewal(distribution=paramDic['model'],u=u,q_period=q_period,incOpenInt=paramDic['incOpenInt']) # First fit a renewal model
		modelParams = [bestParam[0],bestParam[2]] # Extract initial shape and scale parameters
		fit_mask = [paramDic['R_flag'] =='Free',paramDic['z0_flag'] =='Free',modelParams[0]=='Free',modelParams[1]=='Free'] # Create mask for params to fit. modelParams always set to values based on above.
		init_vals = np.array([R_init,z0_init,modelParams0_init,modelParams1_init]) # Array of initial guesses 
		params = init_vals[fit_mask] # These are the paramaters that are being fit
		boundsInit = np.array([Rbounds,Zbounds,P0bounds,P1bounds])  # Bounds used for fitting
		bounds = boundsInit[fit_mask] # Select the appropriate bounds
		#### Now fit GLTFM to find Z0
		res = opt.minimize(neg_logLike_GLTFM, params, args=(paramDic['model'],u,paramDic['incOpenInt'],q_period,paramDic['R_flag'],paramDic['z0_flag'],modelParams[0],modelParams[1],fit_mask), method = method,bounds=bounds,options=options)
		#### Update the u's to account for initial residual strain using initial Z0 estimate
		u_alt = u + z_val(u,paramDic['R_flag'],res.x[0])[:-1] # Temporary u calculation
		#### Now refit the renewal model with the lengthened inter-event times incorporating residual strain times
		bestParam,loglike = fit_Renewal(distribution=paramDic['model'],u=u_alt,q_period=q_period+z_val(u,paramDic['R_flag'],res.x[0])[-1],incOpenInt=paramDic['incOpenInt'])
		modelParams = [bestParam[0],bestParam[2]] # Extract final shape and scale parameters	
		#### Now re-fit GLTFM to find final Z0
		res = opt.minimize(neg_logLike_GLTFM, params, args=(paramDic['model'],u,paramDic['incOpenInt'],q_period,paramDic['R_flag'],paramDic['z0_flag'],modelParams[0],modelParams[1],fit_mask), method = method,bounds=bounds,options=options)
		paramDic['modelParams'] = modelParams # Update modelParams in dictionary to be this best fitting values
	### Extract fitted paramaters from fitted values
	loc_idx = 0
	# Extract best fitting R value
	if fit_mask[0]: # if being fit
		R_out = res.x[loc_idx]
		loc_idx += 1
	else:
		R_out = paramDic['R_flag']
	if isinstance(R_out,float) or isinstance(R_out,int):
		R_out = np.repeat(R_out,len(u)+1)
	# Extract z0_out
	if fit_mask[1]: # if being fit
		z0_out = res.x[loc_idx]
		loc_idx += 1
	else:
		z0_out = paramDic['z0_flag']
	z_out = z_val(u,R_out,z0_out)	
	# Extract P0 (shape parameter)
	if fit_mask[2]: # if being fit
		P0_out = res.x[loc_idx]
		loc_idx += 1
	else:
		P0_out = paramDic['modelParams'][0]
	# Extract P1 (scale parameter)
	if fit_mask[3]: # if being fit
		P1_out = res.x[loc_idx]
	else:
		P1_out = paramDic['modelParams'][1]
	return R_out,z_out,P0_out,P1_out,-1*res.fun


def neg_logLike_GLTFM(params,mod_name,u,incOpenInt,q_period,R,z0,P0,P1,fit_mask):
	#####
	# Function to find estimates for R,z0,P0,P1
	####
	loc_idx = 0
	### Determine the variables of interest for fitting
	if fit_mask[0]:
		R = params[loc_idx]
		loc_idx += 1
	if fit_mask[1]:
		z0 = params[loc_idx]
		loc_idx += 1
	if fit_mask[2]:
		P0 = params[loc_idx]
		loc_idx += 1
	if fit_mask[3]:
		P1 = params[loc_idx]
	z = z_val(u,R,z0) # Calculate z values
	val = neg_loglike_Gen(u,q_period,z,mod_name,P0,P1,incOpenInt)
	return val


def neg_loglike_Gen(u,q_period,z,mod_name,P0,P1,incOpenInt):
	#####
	# Function to find estimates for R,z0,P0,P1 depending on whether open or closed interval included
	####
	GLTFM_rv_closed = GLTFM_gen(z=z[:-1],model=mod_name,P0=P0,P1=P1) # Create GLTFM PDFs for closed intervals
	if incOpenInt == 'Y':
		# Create GLTFM PDF for open interval
		GLTFM_rv_open = GLTFM_gen(z=z[-1],model=mod_name,P0=P0,P1=P1)
		# Calculate loglikelihood for open and closed periods
		val = np.sum(np.log(GLTFM_rv_closed.pdf(u)))+np.log(1-GLTFM_rv_open.cdf(q_period))
	elif incOpenInt == 'N':
		# Calculate loglikelihood for closed periods
		val = np.sum(np.log(GLTFM_rv_closed.pdf(u)))
	return -1*val



#### Fit renewal models depending on 
def fit_Renewal(*args,distribution,u,q_period,incOpenInt):
	#####
	# Wrapper function for fitting standard probability models to observed paleoseismic
	# record. Leverages existing Scipy modules.
	#####
	params = getattr(stats,distribution).fit(u,floc=float(0)) # Initial parameter guess uses closed intervals only
	if distribution == 'expon': # If exponential distribution
		res = opt.minimize(neg_logLike_renewal, [params[1]], args=(distribution,incOpenInt,u,q_period), method = 'L-BFGS-B') # Perform optimization 
		bestparams = [0,res.x[0]]	
	else: # For 2 paramater models
		res = opt.minimize(neg_logLike_renewal, [params[0],params[2]], args=(distribution,incOpenInt,u,q_period), method = 'L-BFGS-B') # Perform optimization 
		bestparams = [res.x[0],0,res.x[1]] # Add location shift of 0 back in 
	loglike = -res.fun
	return bestparams,loglike
	
def neg_logLike_renewal(params,distribution,incOpenInt,u,q_period):
	#####
	# Minimization function for commonly used probability models
	#####	
	if distribution =='expon':
		rv_mod = getattr(stats,distribution)(0,params[0]) # Get PDF
	else: 
		rv_mod = getattr(stats,distribution)(params[0],0,params[1]) # Get PDF
	if incOpenInt == 'Y':
		val = np.sum(np.log(rv_mod.pdf(u)))+np.log(1-rv_mod.cdf(q_period))
	elif incOpenInt == 'N':
		val = np.sum(np.log(rv_mod.pdf(u)))
	return -1*val	

def z_val(u,R,z0):
	### Function 
	### u: inter-event time(s)
	### r: earthquake drop size in probability
	### z0: residual strain (in years) after first earthquake 
	zs = np.zeros(len(u))
	z = z0
	if isinstance(R,float): # Convert single R to array with length equal to num. eqs.
		R = np.repeat(R,len(u)+1)
	for idx in range(len(u)): # loop through inter-event times
		zs[idx] = max(0,z+u[idx]-R[idx+1]) # Calculate y
		z = zs[idx] 
	return np.concatenate(([z0],zs))
	

	

class GLTFM_gen(stats.rv_continuous):
	#### Define a new class for GLTFM probability distributions based on the scipy.stats.rv_continous class
    
    # define init with z (residual strain), model (base), P0 (param 0) and P1 (param 1)
	def __init__(self, z,P0,P1,model,*args,**kwargs):
		super().__init__(*args, **kwargs)
        # init function
		self.z = z #
		self.P0 = P0
		self.P1 = P1
		self.model = model		
		
	# pdf function for residual strain
	def _pdf(self, u):
		rv_mod = getattr(stats,self.model)(self.P0,0,scale=self.P1) # Define base model
		PDF_z = rv_mod.pdf(u+self.z)/(1-rv_mod.cdf(self.z))
		return PDF_z
	
	# cdf function for residual strain
	def _cdf(self,u):
		rv_mod = getattr(stats,self.model)(self.P0,0,scale=self.P1) # Define base model
		CDF_z = (rv_mod.cdf(u+self.z)-rv_mod.cdf(self.z))/(1-rv_mod.cdf(self.z))
		return CDF_z
	
	# Hazard rate function for residual strain	
	def haz_z(self,u):
		rv_mod = getattr(stats,self.model)(self.P0,0,scale=self.P1) # Define base model
		Haz_z = rv_mod.pdf(u+self.z)/(1-rv_mod.cdf(u+self.z))
		return Haz_z
	
	# T-year forecast function for residual strain
	def Tforecast(self,u,T_Forecast):
		#####
		# Function calculates the T-year forecast for the GLTFM model. Aka probability of an
		# earthquake in the next T-years.
		#####
		rv_mod = getattr(stats,self.model)(self.P0,0,scale=self.P1) # Define base model
		forecast = (rv_mod.cdf(u+self.z+T_Forecast)-rv_mod.cdf(u+self.z))/(1-rv_mod.cdf(u+self.z))
		return forecast     
	
	# Used for scoring models only	
	def pdfsq(self, u):
		rv_mod = getattr(stats,self.model)(self.P0,0,scale=self.P1) # Define base model
		PDFsq_z =(rv_mod.pdf(u+self.z)/(1-rv_mod.cdf(self.z)))**2
		return PDFsq_z
		
class GLTFM_gen_rand(GLTFM_gen):
	#### Define a new sub-class for GLTFM specifically for sampling random numbers 
	def __init__(self,*args,**kwargs):
		super().__init__(*args, **kwargs)
		self.ppf_func = self.create_ppf()
		
	def create_ppf(self):
		xs = np.arange(0, 1000001,1)
		func_ppf = interpolate.interp1d(self.cdf(xs), xs, fill_value='extrapolate')
		return func_ppf	
		
   # inverse cdf function
	def _ppf(self, u):
 		return self.ppf_func(u)  
	
