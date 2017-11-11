import os
import numpy as np
import math
import matplotlib.pyplot as plt
from .io import load_input_pickle, save_input_pickle, load_output_pickle, save_output_pickle
from .check import check_params
from .math_funcs import ampa, hill

class Simulation:
    """
    A class for one iteration of model
    """

    def __init__(self, name=None, params=None, params_from_file = False):

        if name:
            self.name = name
        else:
            self.name = math.floor(np.random.uniform()*1e5)
            
        print("name : "+str(self.name))
        
        if params:
            self.params = params
            
        else:
            if params_from_file:
                self.params = load_input_pickle(params)
                
            else:
                #Define default params
                self.params = {
                                                ### General Simulation Params

                    "fs" : 2e4,                 # per second
                    "sweep_length" : 1,         # seconds
                    "num_trials" : 300,         # number of simulation trials to run
                    "num_stim" : 2,             # number of AP stimuli to run
                    "stim_int" : 0.05,          # interval of AP stimuli (sec)
                    "stim1_time" : 0.1,         # time of first stimulus (sec)

                                                ### Calcium Channel Params

                    "num_syn" : 100,            # number of synapses
                    "num_cav" : 1,              # number of voltage-gated calcium channels in complex
                                                # NOTE: [Ca] comes from all channels, as if
                                                # in a cluster (but still open independently)
                    "cav_p_open" : 0.83,        # probability of CaV opening per AP
                    "cav_i" : 1,                # single CaV current per AP (arbitrary units)
                    "ca_decay" : 0.05,          # exponential decay time constant for [ca]
                    "num_cav_ratio" : 2,        # mult num_cav and divide cav_i by this

                                                ### Transmission / Hill Function Params
                    "vesicle_prox" : 0.25,      # Probability your vesicle is nearby the CaV complex
                    "ca_coop" : 3.72,           # Hill Function param N for calcium ion cooperativity
                                                # in vesicular release, see Scimemi & Diamond
                    "ca_ec50" : 0.7,            # Hill Function param EC_50 for calcium concentration
                                                # for half-maximal vesicle release probability                

                                                ### Vesicular Depletion Params

                    #probably need num vesicles in RRP? and to adjust exponential resetting func?
                    "depletion_on" : False,     # Turn on depletion?
                    "pool_tau" : 1.00,          # Readily-Releasable Pool recovery time constant.
                                                # Here, we model the RRP as a binary variable, either
                                                # occupied by one vesicle or not. Thus, we use this
                                                # to define a probability of occupancy [0,1)

                                                ### Postsynaptic Params

                    "quantal_size" : -10,       # single AMPA current per vesicle (pA)
                    "a_tau" : (0.001,0.005)     # fast and slow AMPA current time constants

                                                # NOTE: EC_50 and [Ca] decay constant were constrained 
                                                #       to approximate observed changes in EPSC amplitude
                                                #       by SKF and Baclofen, given 40% reduction in
                                                #       I_cal = N_cav * p_open * cav_i       
                    }
        
        self.default_runs = []
        print("")
        print("Running Model with Default Parameters...")
        self.default_runs.append(self.run_default())


    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return self.name
        
    def run_default(self):
        sim_run = self._runModel()
        return sim_run
    
    def run_modulation(self,parameter="cav_p_open",mod_range=[(x+1)*2/20 for x in range(20)]):
        """
        run a bunch of simulations, modulating one parameter each time
        store output in list
        """
        sim_runs = []
 
        print("Running Modulation of "+parameter+" for range:")
        print(mod_range)
        print("")

        for x in range(len(mod_range)):
            print("Run #"+str(x+1)+" with "+parameter+" at "+str(mod_range[x])+" of default")
            alt_params = self.params.copy()
            alt_params[parameter] *= mod_range[x]
            sim_runs.append(self._runModel(params=alt_params))
        
        print("")
        print("----")
        print("Done with Modulation of "+parameter)
        print("----")

        return (sim_runs,mod_range)

    def run_analysis(self,sim_runs=None):
        """
        takes list of simulation_run objects
        """
        
        if sim_runs is None:
            sim_runs = self.default_runs
            
        print("Running Analysis:")
        print("")
        
        amp = np.zeros(len(sim_runs))
        ppr = np.zeros(len(sim_runs))
        cv_invsq = np.zeros(len(sim_runs))
        mean_epsc = []
        
        for i in range(len(sim_runs)):
            print("Analyzing Run #{0}".format(i+1))
            curr = sim_runs[i]
            amp[i] = curr.epsc_ave[0]
            ppr[i] = curr.epsc_ave[1]/curr.epsc_ave[0]
            cv_invsq[i] = (np.mean(curr.epsc[1,:]**2))/(np.var(curr.epsc[1,:]))

            print("Recreating Mean EPSC from Run #{0}".format(i+1))
            mean_epsc.append(self.plot_epsc_trace(sim_run=curr,plot=False))
#             plt.plot(curr.time,mean_epsc[i])

        beep = lambda x: os.system("echo '\a';sleep 0.5;" * x)
        beep(1)
#         plt.ylabel('Membrane Current (pA)')
#         plt.xlabel('Time, seconds')
#         plt.title('Mean EPSC (across trials)')
#         plt.show()
        
        return (amp,ppr,cv_invsq)
        
    def plot_epsc_trace(self,sim_run=None,plot=True,average=True):
        """
        plots epscs (either all of them, or average)
        """
        
        if sim_run is None:
            sim_run = self.default_runs[0]
        
        no_samples = float(sim_run.params["fs"])*sim_run.params["sweep_length"]
        no_trials = sim_run.params["num_trials"]

        epsc_activity = np.zeros(np.array([no_samples,no_trials]).astype(int))

        cav_inds = np.array(np.where(sim_run.quantal_content))
        big_rows, small_rows, cols = (
            sim_run.ap_inds[cav_inds[0,:]],
            cav_inds[0,:],
            cav_inds[1,:],
            )

        epsc_activity[big_rows.astype(int),cols.astype(int)] = \
            sim_run.quantal_content[small_rows.astype(int),cols.astype(int)]

        # define AMPA kernel
        ampa_kernel = ampa(
            sim_run.time,
            quantal_size = sim_run.params["quantal_size"],
            tau = sim_run.params["a_tau"]
            )

        if average:
            inds = np.ix_(np.arange(epsc_activity.shape[0]),np.arange(epsc_activity.shape[1]))
            Vm_t = np.apply_along_axis(lambda m: np.convolve(m,ampa_kernel), axis=0, \
                arr=np.mean(epsc_activity[inds],axis=1))

            # Vm_t is of shape (no_samples,num_traces)        
            Vm_t = Vm_t[0:len(sim_run.time)]
            if plot:
                plt.plot(sim_run.time,mean_epsc)
                plt.ylabel('Membrane Current (pA)')
                plt.xlabel('Time, seconds')
                plt.title('Mean EPSC (across trials)')
                plt.show()
        else:
            inds = np.ix_(np.arange(epsc_activity.shape[0]),np.arange(epsc_activity.shape[1]))
            Vm_t = np.apply_along_axis(lambda m: np.convolve(m,ampa_kernel), axis=0, \
                arr=epsc_activity[inds])

            # Vm_t is of shape (no_samples,num_traces)        
            Vm_t = Vm_t[0:len(sim_run.time),:]
            if plot:
                plt.plot(sim_run.time,Vm_t)
                plt.ylabel('Membrane Current (pA)')
                plt.xlabel('Time, seconds')
                plt.title('EPSC (multiple trials)')
                plt.show()

        return Vm_t

    def plot_hill_func(self,sim_run=None,trace=0,synapse=0,average=False):
        """
        to plot where a certain synapse and trace is on the hill function
        (only one trace/synapse at a time)
        """
    
        if sim_run is None:
            sim_run = self.default_runs[0]
        
        cav_hits = sim_run.Ca_t[:,trace,synapse]
        
        p_v_func = hill(np.arange(200)/100.,S=1,ec50=sim_run.params["ca_ec50"],n=sim_run.params["ca_coop"])
        plt.plot(np.arange(200)/100.,p_v_func)
        for i in range(len(cav_hits)):
            plt.plot((cav_hits[i],cav_hits[i]),(0,1))
        plt.show()
    
    def plot_I_ca_trace(self,sim_run=None,trace=0,synapse=0,average=False):
        """
        to plot multiple traces or synapses, put in tuple
        """
        
        if sim_run is None:
            sim_run = self.default_runs[0]
        
        trace = np.array([trace]).flatten()
        synapse = np.array([synapse]).flatten()
        
        ## GETTING [Ca](t) (instead of [Ca](stim_num))
        # define exponential decay as [Ca] kernel
        ca_kernel = np.exp(-sim_run.time/sim_run.params["ca_decay"])

        # generate [Ca](t,trial) by convolving against cav_activity
        # crop for sweep length (note, time of peaks = ap_inds)
        no_samples = float(sim_run.params["fs"])*sim_run.params["sweep_length"]
        no_trials = sim_run.params["num_trials"]
        no_syn = sim_run.params["num_syn"]
        
        cav_activity = np.zeros(np.array([no_samples,no_trials,no_syn]).astype(int))
        
        cav_inds = np.array(np.where(sim_run.cav_openings))
        big_rows, small_rows, cols, chunks = (
            sim_run.ap_inds[cav_inds[0,:]],
            cav_inds[0,:],
            cav_inds[1,:],
            cav_inds[2,:]
            )
            
        cav_activity[big_rows.astype(int),cols.astype(int),chunks.astype(int)] = \
            sim_run.cav_openings[small_rows.astype(int),cols.astype(int),chunks.astype(int)]
        
        Ca_t = np.apply_along_axis(lambda m: np.convolve(m,ca_kernel), axis=0, \
            arr=cav_activity[np.ix_(np.arange(cav_activity.shape[0]),trace,synapse)])

        # Ca_t is of shape (no_samples,num_traces,num_synapses)        
        Ca_t = Ca_t[0:len(sim_run.time),:,:]
        
        if average:
            Ca_per_syn = np.mean(Ca_t,axis=1)
            plt.plot(sim_run.time,Ca_per_syn)
            plt.ylabel('[Ca], arbitrary units')
            plt.xlabel('Time, seconds')
            plt.title('Mean [Ca] per bouton (across trials)')
            plt.show()
            Ca_per_trace = np.mean(Ca_t,axis=2)
            plt.plot(sim_run.time,Ca_per_trace)
            plt.ylabel('[Ca], arbitrary units')
            plt.xlabel('Time, seconds')
            plt.title('Mean [Ca] per trial (across boutons)')
            plt.show()
        else:
            plt.plot(sim_run.time,Ca_t.reshape(cav_activity.shape[0],Ca_t.shape[1]*Ca_t.shape[2]))
            plt.ylabel('[Ca], arbitrary units')
            plt.xlabel('Time, seconds')
            plt.title('[Ca] for multiple trials')
            plt.show()

    
    def _runModel(self,params=None):
        """
        simple model of a synapse

        INPUT:
        -dictionary of parameters for simulation 

        OUTPUT: 
        -"simulation" object that contains the following data:
            - time
            - ap_times
            - ap_inds
            - cav_openings
            - ca_kernel
            - Ca_t
            - p_v_successes
            - quantal_content_per_syn
            - epsc_per_syn
            - quantal_content
            - epsc
    
            (also, the following if Vesicle Depletion is on:)
            (-
            (-
        """
        if params is None:
            params = self.params

        check_params(params)    

        FS = float(params["fs"])
        dt = 1./FS
        no_samples = FS*params["sweep_length"]
        no_trials = params["num_trials"]
        no_stims = params["num_stim"]
        no_syn = params["num_syn"]
        time = np.arange(no_samples)*dt

        # array of AP times    
        ap_times = np.zeros(math.floor(no_samples))
        ap_inds = np.zeros(no_stims)

        for i in range(no_stims):
            curr_stim_time = params["stim1_time"] + i * params["stim_int"]
            curr_ind = math.floor(curr_stim_time*FS)
            try:
                ap_times[curr_ind] = 1
                ap_inds[i] = curr_ind
            except IndexError as err:
                print("INDEXERROR: {0}".format(err))
                print("Stimulation parameters may exceed length of sweep")
                return

        ####################################
        # Simulate Calcium Channel Opening #
        ####################################

    #     print("Simulating Calcium Channel Opening...")

        cav_openings = np.zeros(np.array([no_stims,no_trials,no_syn]).astype(int))
        cav_successes = np.zeros(np.array([no_stims,no_trials,no_syn]).astype(int))

        for i in range(params["num_cav"]*params["num_cav_ratio"]):
            cav_successes = np.random.uniform(size=(no_stims,no_trials,no_syn)) < params["cav_p_open"] 
            cav_openings += cav_successes*params["cav_i"]/params["num_cav_ratio"]    

        # define exponential decay as [Ca] kernel
        ca_kernel = np.exp(-params["stim_int"]*np.arange(no_stims)/params["ca_decay"])

        # generate [Ca](stim_num,trial) by convolving with cav_openings
        # crop for no_stim length

        Ca_t = np.apply_along_axis(lambda m: np.convolve(m,ca_kernel), axis=0, arr=cav_openings)
        Ca_t = Ca_t[0:no_stims,:,:]

        #########################################
        # Simulate Ca-Dependent Vesicle Release #
        #########################################

    #     print("Simulating [Ca]-Dependent Vesicle Release...")

        # apply hill function to obtain probability of vesicle release
        p_v = hill(Ca_t,S=1,ec50=params["ca_ec50"],n=params["ca_coop"])

        # extract values corresponding to action-potential timepoints
        # multiplied by CaV opening success/failure (to prevent vesicle
        # release due purely to residual calcium).
        # Also multiply by a scaling factor that indicates the probability
        # that the vesicle is nearby the calcium channel cluster

        corrected_p = p_v*cav_successes*params["vesicle_prox"]

        # then randomly sample to generate quantal events
        p_v_successes = (np.random.uniform(size=corrected_p.shape) < corrected_p)*1


        if(params["depletion_on"]):

            ################################
            # Simulate Vesicular Depletion #
            ################################
    
    #         print("Simulating Vesicular Depletion...")
    
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


        ###########################
        # Simulate AMPA Responses #
        ###########################    

    #     print("Simulating AMPA Responses...")

        # obtain total quantal content per trial (sum across synapses)
        # in order to obtain total EPSC

        quantal_content = np.sum(quantal_content_per_syn,axis=2)

        # quantify bulk EPSC and individual synapse EPSC
        epsc = quantal_content*params["quantal_size"]
        epsc_per_syn = quantal_content_per_syn*params["quantal_size"]

        epsc_ave = np.mean(epsc,axis=1)

        #####################
        # Packaging Results #
        #####################   

    #   print("Packaging Results....")

        data = {
            "params" : params,
            "time" : time,
            "ap_times" : ap_times,
            "ap_inds" : ap_inds,
            "cav_openings" : cav_openings,
            "ca_kernel" : ca_kernel,
            "Ca_t" : Ca_t,
            "p_v_successes" : p_v_successes,
            "quantal_content_per_syn" : quantal_content_per_syn,
            "epsc_per_syn" : epsc_per_syn,
            "quantal_content" : quantal_content,
            "epsc" : epsc,
            "epsc_ave" : epsc_ave
            }

        sim_run = simulation_run(data)

        return sim_run
        

class simulation_run:
    """
    A simple class for storing data (and parameters) from a simulation run
    """

    def __init__(self, data):
        for key,value in data.items():
            setattr(self,key,value)
    
    def update(self,newdata):
        for key,value in newdata:
            setattr(self,key,value)
