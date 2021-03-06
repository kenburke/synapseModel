import numpy as np
from .math_funcs import hill, diffuse

def _sim_CaV_opening(params, no_stims, no_trials, no_syn, text_display = False):
    '''
    ####################################
    # Simulate Calcium Channel Opening #
    ####################################
    Note that subsequent calcium concentration modelling depends on model detail level:
        -phenom facil = [Ca] ~ cav currents
        -core model = [Ca] ~ cav_currents convolved with exponential decay
        -diffusion model = [Ca] ~ cav_currents diffused along distance and decayed
    '''

    phenom_facil = params['phenom_facil']       # True if phenomenological facilitation
    diffusion_model = params['diffusion_model']       # True if modeling diffusion and synt7

    if text_display:
        print("Simulating Calcium Channel Opening...")

    cav_successes = np.zeros(np.array([no_stims,no_trials,no_syn]).astype(int))
    cav_currents = np.zeros(np.array([no_stims,no_trials,no_syn]).astype(int))

    for i in range(params["num_cav"]*params["num_cav_ratio"]):
        cav_successes = np.random.uniform(size=(no_stims,no_trials,no_syn)) < params["cav_p_open"]
        cav_currents += cav_successes*params["cav_i"]/params["num_cav_ratio"]

    if phenom_facil:
        # calcium concentration defined only by CaV currents
        Ca_t = np.zeros(np.array(cav_currents.shape))
        Ca_t += cav_currents
    
    elif diffusion_model:
        if text_display:
            print("Simulating [Ca] Diffusion...")
            
        # define decay time constant
        ca_kernel = np.exp(-params["stim_int"]*np.arange(no_stims)/params["ca_decay"])
        
        # define initial calcium concentration at given distance
        Ca_t = np.zeros(np.array(cav_currents.shape))
        diffusion_params = [
            'r_cav_ves',
            'cav_i',
            'd_ca',
            'k_on',
            'k_d',
            'B_total',
            'ca_basal'
            ]
            
        ordered_params = [params[x] for x in diffusion_params]
        
        # set distance to be matrix of same size as cav_currents (stim_num,trial,no_syn)
        distances = np.zeros(np.array(cav_currents.shape))
        distances += params['r_cav_ves']
        
        # NOTE: IF YOU WANT TO MAKE DISTANCES RANDOM, here is easiest place to do it
        # put distance and cav_currents matrices into ordered_params at correct location
        ordered_params[0] = distances
        ordered_params[1] = cav_currents
        
        # NOTE: This assumes params['cav_i'] in units of pA (unlike core model)
        ca_initial = diffuse(*ordered_params) + params['ca_basal'] # in micromolar
        
        # generate [Ca](stim_num,trial) by convolving ca_kernel with ca_intial
        # crop for no_stim length
        # this mimics decay due to pumps / other extrusion mechanisms
        
        Ca_t = np.apply_along_axis(lambda m: np.convolve(m,ca_kernel), axis=0, arr=ca_initial)
        Ca_t = Ca_t[0:no_stims,:,:]

    else:
        # define exponential decay as [Ca] kernel
        ca_kernel = np.exp(-params["stim_int"]*np.arange(no_stims)/params["ca_decay"])
    
        # generate [Ca](stim_num,trial) by convolving with cav_currents
        # crop for no_stim length

        Ca_t = np.apply_along_axis(lambda m: np.convolve(m,ca_kernel), axis=0, arr=cav_currents)
        Ca_t = Ca_t[0:no_stims,:,:]
    if diffusion_model:
        return (cav_successes,cav_currents,ca_initial,ca_kernel,Ca_t)
    else:        
        return (cav_successes,cav_currents,ca_kernel,Ca_t)

def _sim_vesicle_release(params,Ca_t,cav_successes,text_display = False):
    '''
    #########################################
    # Simulate Ca-Dependent Vesicle Release #
    #########################################
    Note that facilitation has different mechanisms:
        -phenom model = facilitation through nonlinear scaling
        -core model = facilitation through shift along release function
        -diffusion model = facilitation through secondary high-affinity hill function (e.g. Synt7)
    '''
    
    if text_display:
        print("Simulating [Ca]-Dependent Vesicle Release...")

    phenom_facil = params['phenom_facil']       # True if phenomenological facilitation
    diffusion_model = params['diffusion_model']       # True if modeling diffusion and synt7

    # apply hill function to obtain probability of vesicle release
    p_v = hill(Ca_t,S=1,ec50=params["ca_ec50"],n=params["ca_coop"])

    # extract values corresponding to action-potential timepoints
    # multiplied by CaV opening success/failure (to prevent vesicle
    # release due purely to residual calcium).
    # Also multiply by a scaling factor to indicate maximum release probability

    if phenom_facil:
        # facilitate according to a simple parameter 'A' if previous success
        p_v[1,:,:] *= 1 + params['phenom_param'] * (1 - p_v[1,:,:]) * cav_successes[0,:,:]
        corrected_p = p_v*cav_successes*params["vesicle_prox"]
        
    elif diffusion_model:
        # facilitate according to secondary hill function
        # NOTE: because of decay, second pulse is likely facilitated already
        # in a way that depends on first hill function
        # FIRST: scale according to 
        corrected_p = p_v*cav_successes*params["vesicle_prox"]
        corrected_p *= hill(Ca_t,S=params["S_facil"],ec50=params["ca_ec50_facil"],n=params["ca_coop_facil"]) + 1.0
        
    else:
        corrected_p = p_v*cav_successes*params["vesicle_prox"]
        
    # then randomly sample to generate quantal events
    p_v_successes = (np.random.uniform(size=corrected_p.shape) < corrected_p)*1

    return(p_v,corrected_p,p_v_successes)

def _sim_vesicle_depletion(params,p_v_successes,no_stims,text_display = False):
    '''
    ################################
    # Simulate Vesicular Depletion #
    ################################
    (if off, just return p_v_successes as quantal_content_per_syn)
    '''

    if(params["depletion_on"]):

        if text_display:
            print("Simulating Vesicular Depletion...")

        # Now, release_successes must be multiplied by whether or not a vesicle is present
        # So we must model the readily-releasable pool and its depletion / replenishment
        # To do this, we need a matrix of the shape of pool_size that indicates whether
        # the RRP is occupied (modeled as only 1 vesicle with exponential recovery):

        pool_size = np.zeros(p_v_successes.shape,dtype=int) + 1

        pool_recovery_fraction = 1-np.exp(-params["stim_int"]/params["pool_tau"]) # frac recovered every stim interval

        release_successes = np.zeros(p_v_successes.shape,dtype=int)

        # Now we iterate through each stimulus (in parallel) and check if the empty pools
        # recover (doing a flat fraction according to pool_tau generates exponential recovery
        # on average)

        for i in range(no_stims):
            # first calculate whether a vesicle is present
            # for pool_size == 0, check rand against exponential recovery
            pool_zero_inds = np.array(np.where(pool_size[i,:,:]==0))
            cols,chunks = (pool_zero_inds[0,:],pool_zero_inds[1,:])
            rows = np.zeros(cols.shape,dtype=int)+i
            pool_size[rows,cols,chunks] = (np.random.uniform(size=chunks.size) < pool_recovery_fraction ) * 1

            # now check if there was a successful vesicle release (pool_size * p_v_successes)
            release_successes[i,:,:] = pool_size[i,:,:]*p_v_successes[i,:,:]

            # if you're still within range, for all successes, set next stim's pool size to zero
            if i+1<no_stims:
                release_inds = np.array(np.where((p_v_successes[i,:,:]==1)&(pool_size[i,:,:]==1)))
                cols,chunks = (release_inds[0,:],release_inds[1,:])
                rows = np.zeros(cols.shape,dtype=int)+i+1     # next stimulus
                pool_size[rows,cols,chunks] = 0
                # for all inds where the pool_size was 0, nothing (they get another chance to reload)

        quantal_content_per_syn = release_successes
    else:
        quantal_content_per_syn = p_v_successes

    return quantal_content_per_syn

def _sim_ampa_responses(params,quantal_content_per_syn,text_display = False):
    '''
    ###########################
    # Simulate AMPA Responses #
    ###########################

    NOTE: If you wanted to integrate into NEURON, here is where you would do it.
    Take "quantal_content_per_syn" and feed as lookup table for a NEURON
    simulation with trains of M stimuli, over N trials and O synapses
    for a quantal_content_per_syn of shape MxNxO = (no_stims,no_trials,no_syn)

    Either that, or incorporate the math found in this file into the NEURON synapse directly
    '''
    if text_display:
        print("Simulating AMPA Responses...")

    # obtain total quantal content per trial (sum across synapses)
    # in order to obtain total EPSC

    quantal_content = np.sum(quantal_content_per_syn,axis=2)

    # quantify bulk EPSC and individual synapse EPSC
    epsc = quantal_content*params["quantal_size"]
    epsc_per_syn = quantal_content_per_syn*params["quantal_size"]

    epsc_ave = np.mean(epsc,axis=1)

    return(quantal_content,epsc,epsc_per_syn,epsc_ave)
