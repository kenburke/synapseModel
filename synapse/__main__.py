import sys
import numpy as np
import matplotlib.pyplot as plt
from .i_o import load_input_pickle, get_user_params, get_user_modulation, dumpclean
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

outer_param = 'ca_decay'
#outer_range = [0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.5]
outer_range = [0.05]
base_val = SIM.params[outer_param]

outer_amp = []
outer_ppr = []
outer_cv_invsq = []
outer_epsc= []

for x in range(len(outer_range)):
    print("OUTER RUN #"+str(x+1)+" with "+outer_param+" set to "+str(outer_range[x]))
    SIM.params[outer_param] = outer_range[x]
    SIM.run_modulation(parameter=mod_param,mod_range=mod_range)
    
    print("")
    print("----")   
    print("")
    print('Running Analysis on most recent modulation run, with following params : ')
    print("OUTER PARAM "+outer_param+" = "+str(outer_range[x]))
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
    i_o.save_output_plot(SIM.plot_path,'EPSC_v_Time_'+str(x))
    plt.close()


print("Saving Plots to "+SIM.plot_path)

for x in range(len(outer_ppr)):
    amp,ppr = outer_amp[x], outer_ppr[x]
    plt.plot(amp/amp[-1],ppr/ppr[-1])

plt.xlabel('EPSC Amplitude (norm.)')
plt.ylabel('Paired-Pulse Ratio (norm.)')
plt.title('PPR vs. Amplitude')
plt.axis([0.0,1.1,0.0,3.0])
i_o.save_output_plot(SIM.plot_path,'PPR_v_Amp')
plt.close()

for x in range(len(outer_ppr)):
    amp,cv_invsq = outer_amp[x], outer_cv_invsq[x]
    plt.plot(amp/amp[-1],cv_invsq/cv_invsq[-1])

plt.xlabel('EPSC Amplitude (norm.)')
plt.ylabel('C.V.^-2 (norm.)')
plt.title('C.V.^-2 vs. Amplitude')
plt.axis([0.0,1.1,0.0,3.0])
i_o.save_output_plot(SIM.plot_path,'CV_v_Amp')
plt.close()

SIM.params[outer_param] = base_val


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