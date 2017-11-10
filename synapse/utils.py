from .model import runModel, ampa, hill
from matplotlib import pyplot as plt
import numpy as np
import math
import os


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
        if params is None:
            params = self.params
        return runModel(params)

#     data = {
#         "params" : params,
#         "time" : time,
#         "ap_times" : ap_times,
#         "ap_inds" : ap_inds,
#         "cav_openings" : cav_openings,
#         "ca_kernel" : ca_kernel,
#         "Ca_t" : Ca_t,
#         "p_v_successes" : p_v_successes,
#         "quantal_content_per_syn" : quantal_content_per_syn,
#         "epsc_per_syn" : epsc_per_syn,
#         "quantal_content" : quantal_content,
#         "epsc" : epsc,
#         "epsc_ave" : epsc_ave
#         }
