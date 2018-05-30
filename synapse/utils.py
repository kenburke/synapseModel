import os
import numpy as np
import math
import matplotlib.pyplot as plt
import datetime
from .i_o import load_input_pickle, save_session_pickle, get_user_params, save_output_plot
from .check import check_params
from .math_funcs import ampa, hill, diffuse
from .model import _sim_CaV_opening, _sim_vesicle_release, _sim_vesicle_depletion, _sim_ampa_responses

class Simulation:
    """
    A class for one iteration of model

    Major attributes:
        - params:
            -current set of parameters if you call _runModel (without specifying new parameters)
        - default_runs
            -list of "simulation_run" objs that contain results runs of a simulation, including:
                -sim_run.params, same format as Simulation.params, for replication

    """

    def __init__(self, name=None, params=None, params_from_file=False, params_from_user=False):
        """Init by loading param file and running one simulation"""

        print("")
        if name:
            self._name = name
        else:
            self._name = input("Simulation Name : ")

        print("Name : "+str(self._name))

        self.plot_path = os.getcwd()+'/session/'+self._name+'_plots/'
        try:
            os.mkdir(self.plot_path)
        except (FileExistsError):
            beep = lambda x: os.system("echo '\a';sleep 0.5;" * x)
            beep(1)
            print("WARNING : FOLDER PATH ALREADY EXISTS")
            print(self.plot_path)
            print("WRITING OVER")
            for fn in os.listdir(self.plot_path):
                os.remove(self.plot_path+fn)

        if params:
            self.params = params
        else:
            if params_from_file:
                self.params = load_input_pickle(params_from_file)
            elif params_from_user:
                self.params = get_user_params()
            else:
                #Define default params
                self.params = load_input_pickle('default')

        self.default_runs = []      # array of simulation runs with default parameters
        self.mod_runs = []          # array of tuples that contain 0) a list of simulation runs
                                    # and 1) a dictionary clarifying which parameter was given
                                    # which value for each run. (for convenience, can also
                                    # determine by comparing the simulation_run.params
                                    # directly


        print("Running Model with Default Parameters...")
        self.run_default()
        print("")

    def __str__(self):
        """Print out name, current params and stored simulation runs"""

        print("")
        s = "NAME : "+self._name+"\n\n"
        s += "PARAMS :"
        print(s)

        for key, val in self.params.items():
            l = (21-len(key))//7
            print("{0}".format(key)+"\t"*l+":\t{0}".format(val))

        s = "\nRuns stored in DEFAULT_RUNS = "+str(len(self.default_runs))
        print(s)

        s = "\nRuns stored in MOD_RUNS = "+str(len(self.mod_runs))
        print(s)

        return ""

    def save(self,fn=None):
        """
        Save this model into .pkl file
        """

        if (fn==None):
            fn = self._name

        print("")
        print("Saving Simulation Object into session/"+fn+".pkl")
        return save_session_pickle(self,fn)
        print("----")

    def replicate(self,simulation_run):
        """
        Replicate (i.e. rerun) a simulation with the same params (but different random seed)
        WARNING: DOES NOT CURRENTLY REPLICATE NONUNIFORM SYNAPSES
        
        INPUT: simulation_run object, already completed
        OUTPUT: simulation_run object from another simulation with same params as input
        """

        return self._runModel(params=simulation_run.params)

    def run_default(self):
        """
        Runs model with current default parameters (stored in self.params)
        INPUT: none
        OUTPUT: simulation_run object
        """
        sim_run = self._runModel()

        self.default_runs.append(sim_run)

        return sim_run

    def run_modulation(self,parameter="cav_p_open",mod_range=[(x+1)/20 for x in range(20)],synapse_distr=False):
        """
        run a bunch of simulations, modulating one parameter each time
        synapse_distr should be a tuple from i_o.get_synapse_range()
        store output in tuple containing list of simulation run objects and mod_range
        """
        sim_runs = []

        print("Running Modulation of "+parameter+" for range:")
        print(mod_range)
        print("")

        for x in range(len(mod_range)):
            print("Run #"+str(x+1)+" with "+parameter+" at "+str(mod_range[x])+" of default")
            alt_params = self.params.copy()
            alt_params[parameter] *= mod_range[x]
            sim_runs.append(self._runModel(params=alt_params,nonuniform_parameter=synapse_distr))

        print("")
        print("----")
        print("Storing sim_runs in mod_runs as\n(sim_runs, mod_dict)\nwhere mod_dict = {parameter:mod_range}\n")
        self.mod_runs.append((sim_runs,{parameter:mod_range}))
        print("Done with Modulation of "+parameter)
        print("----")

        return (sim_runs,{parameter:mod_range})

    def run_analysis(self,sim_runs=None,notify=False):
        """
        takes list of simulation_run objects
        runs a bunch of basic analyses (PPR, CV, etc) and saves plots
        """

        if sim_runs is None:
            sim_runs = self.default_runs

        print("Analysis Start:")
        print("")

        amp = np.zeros(len(sim_runs))
        ppr = np.zeros(len(sim_runs))
        cv_invsq = np.zeros(len(sim_runs))
        mean_epsc = []

        for i in range(len(sim_runs)):
            print("Analyzing Run #{0}".format(i+1))
            curr = sim_runs[i]
            amp[i] = curr.data["epsc_ave"][0]
            ppr[i] = curr.data["epsc_ave"][1]/curr.data["epsc_ave"][0]
            cv_invsq[i] = (np.mean(curr.data["epsc"][1,:]**2))/(np.var(curr.data["epsc"][1,:]))

            print("Recreating Mean EPSC from Run #{0}".format(i+1))
            mean_epsc.append(self.plot_epsc_trace(sim_run=curr,plot=False))

        if notify:
            beep = lambda x: os.system("echo '\a';sleep 0.5;" * x)
            beep(2)

        print("")
        print("Saving Plots to "+self.plot_path)

        plt.plot(amp/amp[-1],ppr/ppr[-1])
        plt.xlabel('EPSC Amplitude (norm.)')
        plt.ylabel('Paired-Pulse Ratio (norm.)')
        plt.title('PPR vs. Amplitude')
        plt.axis([0.0,1.1,0.0,3.0])
        save_output_plot(self.plot_path,'PPR_v_Amp')
        plt.close()

        plt.plot(amp/amp[-1],cv_invsq/cv_invsq[-1])
        plt.xlabel('EPSC Amplitude (norm.)')
        plt.ylabel('C.V.^-2 (norm.)')
        plt.title('C.V.^-2 vs. Amplitude')
        plt.axis([0.0,1.1,0.0,3.0])
        save_output_plot(self.plot_path,'CV_v_Amp')
        plt.close()

        for i in range(len(mean_epsc)):
            plt.plot(curr.data['time'],mean_epsc[i])
        plt.xlabel('Time (seconds)')
        plt.ylabel('Mean EPSC (pA)')
        plt.title('Mean EPSC across conditions')
        save_output_plot(self.plot_path,'EPSC_v_Time')
        plt.close()

        return (amp,ppr,cv_invsq,mean_epsc)

    def plot_epsc_trace(self,sim_run=None,plot=True,average=True):
        """
        plots epscs (either all of them, or average)
        """

        if sim_run is None:
            sim_run = self.default_runs[0]

        no_samples = float(sim_run.params["fs"])*sim_run.params["sweep_length"]
        no_trials = sim_run.params["num_trials"]

        epsc_activity = np.zeros(np.array([no_samples,no_trials]).astype(int))

        cav_inds = np.array(np.where(sim_run.data["quantal_content"]))
        big_rows, small_rows, cols = (
            sim_run.data["ap_inds"][cav_inds[0,:]],
            cav_inds[0,:],
            cav_inds[1,:],
            )

        epsc_activity[big_rows.astype(int),cols.astype(int)] = \
            sim_run.data["quantal_content"][small_rows.astype(int),cols.astype(int)]

        # define AMPA kernel
        ampa_kernel = ampa(
            sim_run.data["time"],
            quantal_size = sim_run.params["quantal_size"],
            tau = sim_run.params["a_tau"]
            )

        if average:
            inds = np.ix_(np.arange(epsc_activity.shape[0]),np.arange(epsc_activity.shape[1]))
            Vm_t = np.apply_along_axis(lambda m: np.convolve(m,ampa_kernel), axis=0, \
                arr=np.mean(epsc_activity[inds],axis=1))

            # Vm_t is of shape (no_samples,num_traces)
            Vm_t = Vm_t[0:len(sim_run.data["time"])]
            if plot:
                plt.plot(sim_run.data["time"],Vm_t)
                plt.ylabel('Membrane Current (pA)')
                plt.xlabel('Time, seconds')
                plt.title('Mean EPSC (across trials)')
                plt.show()
        else:
            inds = np.ix_(np.arange(epsc_activity.shape[0]),np.arange(epsc_activity.shape[1]))
            Vm_t = np.apply_along_axis(lambda m: np.convolve(m,ampa_kernel), axis=0, \
                arr=epsc_activity[inds])

            # Vm_t is of shape (no_samples,num_traces)
            Vm_t = Vm_t[0:len(sim_run.data["time"]),:]
            if plot:
                plt.plot(sim_run.data["time"],Vm_t)
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

        cav_hits = sim_run.data["Ca_t"][:,trace,synapse]

        p_v_func = hill(np.arange(200)/100.,S=1,ec50=sim_run.params["ca_ec50"],n=sim_run.params["ca_coop"])
        plt.plot(np.arange(200)/100.,p_v_func)
        for i in range(len(cav_hits)):
            plt.plot((cav_hits[i],cav_hits[i]),(0,1))
        plt.ylabel('Probbility of Vesicle Release')
        plt.xlabel('Calcium Concentration (arb. units)')
        plt.title('Location of [Ca] response on Hill Function for sequential APs')
        plt.show()

    def plot_I_ca_trace(self,sim_run=None,trace=0,synapse=0,average=False):
        """
        to plot multiple traces or synapses, put in tuple
        """

        if sim_run is None:
            sim_run = self.default_runs[0]
        
        diffusion_model = sim_run.params["diffusion_model"]
        
        trace = np.array([trace]).flatten()
        synapse = np.array([synapse]).flatten()

        ## GETTING [Ca](t) (instead of [Ca](stim_num))
        # define exponential decay as [Ca] kernel
        ca_kernel = np.exp(-sim_run.data["time"]/sim_run.params["ca_decay"])

        # generate [Ca](t,trial) by convolving against cav_activity
        # crop for sweep length (note, time of peaks = ap_inds)
        no_samples = float(sim_run.params["fs"])*sim_run.params["sweep_length"]
        no_trials = sim_run.params["num_trials"]
        no_syn = sim_run.params["num_syn"]

        cav_activity = np.zeros(np.array([no_samples,no_trials,no_syn]).astype(int))

        cav_inds = np.array(np.where(sim_run.data["cav_currents"]))
        big_rows, small_rows, cols, chunks = (
            sim_run.data["ap_inds"][cav_inds[0,:]],
            cav_inds[0,:],
            cav_inds[1,:],
            cav_inds[2,:]
            )

        cav_activity[big_rows.astype(int),cols.astype(int),chunks.astype(int)] = \
            sim_run.data["cav_currents"][small_rows.astype(int),cols.astype(int),chunks.astype(int)]
        
        if diffusion_model:
            distances = np.zeros(np.array(cav_activity.shape))
            distances += sim_run.params['r_cav_ves']
                
            diffusion_params = [
                'r_cav_ves',
                'cav_i',
                'd_ca',
                'k_on',
                'k_d',
                'B_total',
                'ca_basal'
                ]
                
            ordered_params = [sim_run.params[x] for x in diffusion_params]            
            ordered_params[0] = distances
            ordered_params[1] = cav_activity
            ca_initial = diffuse(*ordered_params) + sim_run.params['ca_basal'] # in micromolar
            Ca_t = np.apply_along_axis(lambda m: np.convolve(m,ca_kernel), axis=0, \
                arr=ca_initial[np.ix_(np.arange(ca_initial.shape[0]),trace,synapse)])
        else:
            Ca_t = np.apply_along_axis(lambda m: np.convolve(m,ca_kernel), axis=0, \
                arr=cav_activity[np.ix_(np.arange(cav_activity.shape[0]),trace,synapse)])

        # Ca_t is of shape (no_samples,num_traces,num_synapses)
        Ca_t = Ca_t[0:len(sim_run.data["time"]),:,:]

        if average:
            Ca_per_syn = np.mean(Ca_t,axis=1)
            plt.plot(sim_run.data["time"],Ca_per_syn)
            plt.ylabel('[Ca], arbitrary units')
            plt.xlabel('Time, seconds')
            plt.title('Mean Functional [Ca] per bouton (across trials)')
            plt.show()
            Ca_per_trace = np.mean(Ca_t,axis=2)
            plt.plot(sim_run.data["time"],Ca_per_trace)
            plt.ylabel('[Ca], arbitrary units')
            plt.xlabel('Time, seconds')
            plt.title('Mean Functional [Ca] per trial (across boutons)')
            plt.show()
        else:
            plt.plot(sim_run.data["time"],Ca_t.reshape(cav_activity.shape[0],Ca_t.shape[1]*Ca_t.shape[2]))
            plt.ylabel('[Ca], arbitrary units')
            plt.xlabel('Time, seconds')
            plt.title('[Ca] for multiple trials')
            plt.show()


    def _runModel(self,params=None,text_display=False,nonuniform_parameter=False):
        """
        simple model of a synapse

        INPUT:
        -dictionary of parameters for simulation
        -text_display toggle
        -toggle to distribute a parameter across synapses
            --> If tuple is passed, it is assumed to be output of get_synapse_range()

        OUTPUT:
        -"simulation" object that contains the following data:
            - params <- dictionary of the format of Simulation.params, for reproducing experiments
            - time
            - ap_times
            - ap_inds
            - cav_currents
            - ca_kernel
            - Ca_t
            - p_v_successes
            - quantal_content_per_syn
            - epsc_per_syn
            - quantal_content
            - epsc
        """

        if params is None:
            params = self.params
            
        check_params(params)

        if nonuniform_parameter:
            syn_param,syn_range = nonuniform_parameter

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
        if nonuniform_parameter:
            #NOTE: As warned, rounding-up error if no_syn not multiple of len(syn_range)
            r_ves_key = {
                1 : 15,
                2 : 23.3,
                3 : 29.2,
                4 : 33.8,
                5 : 37.6,
                }
            syn_set_size = int(no_syn / (len(syn_range))) + 1
            for x in range(len(syn_range)):
                params[syn_param] = syn_range[x]
                params['r_cav_ves'] = r_ves_key[syn_range[x]]
                calcium_temp = _sim_CaV_opening(
                    params, no_stims, no_trials, syn_set_size, text_display=text_display)
                if (x == 0):
                    if params["diffusion_model"]:
                        cav_successes,cav_currents,ca_initial,ca_kernel,Ca_t = calcium_temp
                    else:
                        cav_successes,cav_currents,ca_kernel,Ca_t = calcium_temp
                else:
                    cav_successes = np.dstack((cav_successes,calcium_temp[0]))
                    cav_currents = np.dstack((cav_currents,calcium_temp[1]))
                    ca_initial = np.dstack((ca_initial,calcium_temp[2]))
                    Ca_t = np.dstack((Ca_t,calcium_temp[4]))
                    
        else:
            if params["diffusion_model"]:
                cav_successes,cav_currents,ca_initial,ca_kernel,Ca_t = _sim_CaV_opening(
                    params, no_stims, no_trials, no_syn, text_display=text_display)
            else:
                cav_successes,cav_currents,ca_kernel,Ca_t = _sim_CaV_opening(
                    params, no_stims, no_trials, no_syn, text_display=text_display)
              
        #########################################
        # Simulate Ca-Dependent Vesicle Release #
        #########################################
        p_v,corrected_p,p_v_successes = _sim_vesicle_release(
            params,Ca_t,cav_successes,text_display=text_display)

        ################################
        # Simulate Vesicular Depletion #
        ################################
        if(params["depletion_on"]):
            quantal_content_per_syn = _sim_vesicle_depletion(
                params,p_v_successes,no_stims,text_display=text_display)
        else:
            quantal_content_per_syn = p_v_successes

        ###########################
        # Simulate AMPA Responses #
        ###########################
        quantal_content,epsc,epsc_per_syn,epsc_ave = _sim_ampa_responses(
            params,quantal_content_per_syn,text_display=text_display)

        #####################
        # Packaging Results #
        #####################

        if text_display:
            print("Packaging Results....")
        
        if params["diffusion_model"]:
            data = {
                "time" : time,
                "ap_times" : ap_times,
                "ap_inds" : ap_inds,
                "cav_successes" : cav_successes,
                "cav_currents" : cav_currents,
                "ca_kernel" : ca_kernel,
                "ca_initial" : ca_initial,
                "Ca_t" : Ca_t,
                "p_v_successes" : p_v_successes,
                "quantal_content_per_syn" : quantal_content_per_syn,
                "epsc_per_syn" : epsc_per_syn,
                "quantal_content" : quantal_content,
                "epsc" : epsc,
                "epsc_ave" : epsc_ave
                }
        else:
            data = {
                "time" : time,
                "ap_times" : ap_times,
                "ap_inds" : ap_inds,
                "cav_successes" : cav_successes,
                "cav_currents" : cav_currents,
                "ca_kernel" : ca_kernel,
                "Ca_t" : Ca_t,
                "p_v_successes" : p_v_successes,
                "quantal_content_per_syn" : quantal_content_per_syn,
                "epsc_per_syn" : epsc_per_syn,
                "quantal_content" : quantal_content,
                "epsc" : epsc,
                "epsc_ave" : epsc_ave
                }


        sim_run = simulation_run(params,data)

        return sim_run


class simulation_run:
    """
    A simple class for storing data (and parameters) from a simulation run
    """

    def __init__(self, params, data):
        self.params = params
        self.data = data

    def __str__(self):
        print(self.__repr__())
        print("\tPARAMS...")
        for key, val in self.params.items():
            l = (21-len(key))//7
            print("\t{0}".format(key)+"\t"*l+":\t{0}".format(val))

        print("\n\tDATA...")
        for key, val in self.data.items():
            l = (21-len(key))//7
            print("\t{0}".format(key))

        return ""
