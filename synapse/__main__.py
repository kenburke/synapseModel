import sys
import numpy as np
import matplotlib.pyplot as plt
from .i_o import load_input_pickle, get_user_params, get_user_modulation, dumpclean, save_output_plot, get_synapse_range
from .utils import Simulation

print("")

# Some quick stuff to make sure the program is called correctly
if len(sys.argv) != 3:
    print("Usage: python -m synapse [-L| -M] <session_name>")
    sys.exit(0)
    
# Choose which effect to model

if sys.argv[1][0:2] == '-L':
    print("Option: Load simulation parameter file")
    input_fn = input("Insert input filename from input folder > ")
    print("")
    SIM = Simulation(name=sys.argv[2],params_from_file = input_fn)
    print("")
    print("Running simulation with params from "+input_fn)

elif sys.argv[1][0:2] == '-M':
    print("Option: Manual insertion of parameters")
    print("")   
    SIM = Simulation(name=sys.argv[2],params_from_user = True)
    print("")   
    print("Running simulation with user-definied params")
    
else:
    print("Usage: python -m synapse [-L| -M] <session_name>")
    sys.exit(0)

# Show user what parameters they're working with
print(SIM)
print("----")   

# get the parameter for modulation and the range of modulation, then run
mod_param,mod_range = get_user_modulation()
syn_param,syn_range = get_synapse_range()

outer_param_first = 'num_cav'
outer_range_first = [1]
base_val_first = SIM.params[outer_param_first]
outer_param_second = 'r_cav_ves'
outer_range_second = [15]
base_val_second = SIM.params[outer_param_second]

outer_amp = []
outer_ppr = []
outer_cv_invsq = []
outer_epsc= []

for x in range(len(outer_range_first)):
    print("OUTER RUN #"+str(x+1)+" with "+outer_param_first+" set to "+str(outer_range_first[x])+" and "+outer_param_second+" set to "+str(outer_range_second[x]))
    SIM.params[outer_param_first] = outer_range_first[x]
    SIM.params[outer_param_second] = outer_range_second[x]
    SIM.run_modulation(parameter=mod_param,mod_range=mod_range,synapse_distr=(syn_param,syn_range))
    
    print("")
    print("----")   
    print("")
    print('Running Analysis on most recent modulation run, with following params : ')
    print("OUTER PARAM "+outer_param_first+" = "+str(outer_range_first[x]))
    print("OUTER PARAM "+outer_param_second+" = "+str(outer_range_second[x]))
    print("INNER PARAM:")
    print(SIM.mod_runs[-1][1])
    print("")
    
    # notify just pings you if it takes a long time
    amp,ppr,cv_invsq,mean_epsc = SIM.run_analysis(sim_runs=SIM.mod_runs[-1][0],notify=True)
    print("----")   
    outer_amp.append(amp)
    outer_ppr.append(ppr)
    outer_cv_invsq.append(cv_invsq)
    outer_epsc.append(mean_epsc)
    
    for i in range(len(mean_epsc)):
        plt.plot(SIM.default_runs[0].data['time'],mean_epsc[i])
    plt.xlabel('Time (seconds)')
    plt.ylabel('Mean EPSC (pA)')
    plt.title('Mean EPSC across conditions')
    save_output_plot(SIM.plot_path,'EPSC_v_Time_'+str(x))
    plt.close()

print("Saving Plots to "+SIM.plot_path)

for x in range(len(outer_ppr)):
    amp,ppr = outer_amp[x], outer_ppr[x]
    plt.plot(amp/amp[-1],ppr/ppr[-1])

plt.xlabel('EPSC Amplitude (norm.)')
plt.ylabel('Paired-Pulse Ratio (norm.)')
plt.title('PPR vs. Amplitude')
plt.axis([0.0,1.1,0.0,3.0])
save_output_plot(SIM.plot_path,'PPR_v_Amp')
plt.close()

for x in range(len(outer_ppr)):
    amp,cv_invsq = outer_amp[x], outer_cv_invsq[x]
    plt.plot(amp/amp[-1],cv_invsq/cv_invsq[-1])

plt.xlabel('EPSC Amplitude (norm.)')
plt.ylabel('C.V.^-2 (norm.)')
plt.title('C.V.^-2 vs. Amplitude')
plt.axis([0.0,1.1,0.0,3.0])
save_output_plot(SIM.plot_path,'CV_v_Amp')
plt.close()

# now plot Cd / Mn style Amp and PPR plots as a function of r_cav_ves at 0.6 modulation
ppr_slice_unmod = []
amp_slice_unmod = []
ppr_slice_mod = []
amp_slice_mod = []
ppr_slice_norm = []
amp_slice_norm = []

for x in range(len(outer_ppr)):
    amp, ppr = outer_amp[x], outer_ppr[x]
    amp_slice_unmod.append(amp[-1]) # unmodulated amp
    amp_slice_mod.append(amp[11]) # modulated amp
    amp_slice_norm.append(amp[11]/amp[-1]) # normalized change in amp
    ppr_slice_unmod.append(ppr[-1])
    ppr_slice_mod.append(ppr[11])
    ppr_slice_norm.append(ppr[11]/ppr[-1])

plt.plot(outer_range_second, amp_slice_unmod)
plt.plot(outer_range_second, amp_slice_mod)
plt.xlabel('Distance between CaV Cluster and Vesicle (nm)')
plt.ylabel('EPSC Amplitude (pA)')
plt.title('EPSC Amp. vs Distance')
plt.axis([0.0,60,0.0,-300])
save_output_plot(SIM.plot_path,'Amp_vs_Distance')
plt.close()

plt.plot(outer_range_second, ppr_slice_unmod)
plt.plot(outer_range_second, ppr_slice_mod)
plt.xlabel('Distance between CaV Cluster and Vesicle (nm)')
plt.ylabel('Paired-Pulse Ratio (raw)')
plt.title('PPR vs Distance')
plt.axis([0.0,60,0.0,3.0])
save_output_plot(SIM.plot_path,'PPR_vs_Distance')
plt.close()

plt.plot(outer_range_second, amp_slice_norm)
plt.plot(outer_range_second, ppr_slice_norm)
plt.xlabel('Distance between CaV Cluster and Vesicle (nm)')
plt.ylabel('EPSC Amp or PPR (norm.)')
plt.title('PPR and Amp. vs Distance')
plt.axis([0.0,60,0.0,3.0])
save_output_plot(SIM.plot_path,'PPR_Amp_vs_Distance')
plt.close()


SIM.params[outer_param_first] = base_val_first
SIM.params[outer_param_second] = base_val_second


#print("")
#print("----")   
#print("")
#SIM.run_modulation(parameter=mod_param,mod_range=mod_range)

## Run most recent modulation analysis
#print('Running Analysis on most recent modulation run, with following params : ')
#print(SIM.mod_runs[-1][1])
#print("")
## notify just pings you if it takes a long time
#amp,ppr,cv_invsq,mean_epsc = SIM.run_analysis(sim_runs=SIM.mod_runs[-1][0],notify=True)

#print("----")   
#print("")
#SIM.save()
#print("")
#######

#print("--------------------")   
#print("")
#print("")
#print("Done")
#print("")
#print("If you would like to work more with this session from command line python,")
#print("i_o.load_session_pickle('{0}')".format(sys.argv[2]))
#print("")