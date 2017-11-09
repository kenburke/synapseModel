import sys
import numpy as np
import matplotlib.pyplot as plt
from .io import write_results, load_input_pickle, get_user_params, dumpclean
from .utils import Simulation
from .model import runModel

print("")
print("")

# Some quick stuff to make sure the program is called correctly
if len(sys.argv) != 3:
    print("Usage: python -m synapse [-L| -M] <output file>")
    sys.exit(0)
    
# Choose which effect to model

if sys.argv[1][0:2] == '-L':
    print("Option: Load simulation parameter file")
    input_fn = input("Insert input filename from input/:")
    params = load_input_pickle(input_fn)
    print("----")   
    print("Running simulation with params from "+input_fn+" :")

elif sys.argv[1][0:2] == '-M':
    print("Option: Manual insertion of parameters")
    params = get_user_params()
    print("----")   
    print("Running simulation with user-definied params:")

elif sys.argv[1][0:2] == '-K':
    print("Option: Running with defaults")
    
else:
    print("Usage: python -m synapse [-L| -M] <output file>")
    sys.exit(0)


print("")
print("----")
print("Initializing Simulation")
print("")
SIM = Simulation()
print("----")   
print("")
num_cav_ratio_mod = [1,2,3,5,7,10,20,50,100]
mod_range = (np.arange(10)+1)*0.1
cv_runs = np.zeros((len(mod_range),len(num_cav_ratio_mod)))
ppr_runs = np.zeros((len(mod_range),len(num_cav_ratio_mod)))
amp_runs = np.zeros((len(mod_range),len(num_cav_ratio_mod)))

for i in range(len(num_cav_ratio_mod)):
    print("")
    print("num_cav_ratio = {0}".format(num_cav_ratio_mod[i]))   
    print("")
    SIM.params["num_cav_ratio"] = num_cav_ratio_mod[i]
    SIM.mod_runs,mod = SIM.run_modulation(parameter="cav_p_open",mod_range=mod_range)
    print("")
    print("")
    print("----")   
    r = len(mod)-1
    amp,ppr,cv = SIM.run_analysis(SIM.mod_runs)
    print("PPR = {0}".format(ppr))
    print("")
    print("----")   
    for j in range(len(mod_range)):
        ppr_runs[j,i] = ppr[j]
        amp_runs[j,i] = amp[j]
        cv_runs[j,i] = cv[j]

for i in range(len(num_cav_ratio_mod)):
    plt.plot(amp_runs[:,i]/amp_runs[r,i],ppr_runs[:,i]/ppr_runs[r,i])
plt.show()  
  
for i in range(len(num_cav_ratio_mod)):
    plt.plot(amp_runs[:,i]/amp_runs[r,i],cv_runs[:,i]/cv_runs[r,i])
plt.show()    

print(plt.rcParams['axes.prop_cycle'].by_key()['color'])


# for j in range(len(mod_range)):
#     ppr_runs[j,i] = ppr[j]
#     amp_runs[j,i] = amp[j]
#     cv_runs[j,i] = cv[j]


# plt.plot(num_cav_ratio_mod,ppr_runs[0,:])
# plt.plot(num_cav_ratio_mod,ppr_runs[1,:])
# plt.show()
# plt.plot(amp_runs[0,:],ppr_runs[0,:])
# plt.plot(amp_runs[1,:],ppr_runs[1,:])
# plt.show()
# plt.plot(amp_runs[0,:]/amp_runs[1,:],ppr_runs[0,:]/ppr_runs[1,:])
# plt.show()
# plt.plot(num_cav_ratio_mod,ppr_runs[0,:]/ppr_runs[1,:])
# plt.show()
# 
# for i in range(len(mod)):
#     print(mod[i])
#     print(amp[i]/amp[r])
#     print(ppr[i]/ppr[r])
#     print(cv[i]/cv[r])
#     print("")

# plt.plot(amp/amp[r],ppr/ppr[r])
# plt.show()    
# plt.plot(amp/amp[r],cv/cv[r])
# plt.show()    

# print("Plotting Basline Traces...")
# SIM.plot_I_ca_trace(SIM.mod_runs[1],trace=(0,1,2,3),synapse=(0,1,2,3,4,5,6,7,8),average=False)
# print("...")
# SIM.plot_epsc_trace(SIM.mod_runs[1],trace=np.arange(20),average=False)
# print("...")
# SIM.plot_epsc_trace(SIM.mod_runs[9],trace=np.arange(50),average=True)
# print("Plotting Modulation Traces...")
# SIM.plot_I_ca_trace(SIM.mod_runs[0],trace=(0,1,2,3),synapse=(0,1,2,3,4,5,6,7,8),average=False)
# print("...")
# SIM.plot_epsc_trace(SIM.mod_runs[0],trace=np.arange(20),average=False)
# print("...")
# for i in range(9):
#     print("Plotting Modulation Level {0}...".format(i*0.1+0.1))
#     SIM.plot_epsc_trace(SIM.mod_runs[i],trace=np.arange(50),average=True)

# print("Average Quantal Content")
# print(np.mean(SIM.mod_runs[1].quantal_content,axis=1))
# print(np.mean(SIM.mod_runs[0].quantal_content,axis=1))
# print("Average P_Vesicle per Action Potential (should be between 0.1 and 0.2)")
# print(np.mean(SIM.mod_runs[1].quantal_content,axis=1)/50)
# print(np.mean(SIM.mod_runs[0].quantal_content,axis=1)/50)
# print("Plotting Baseline Hill...")
# SIM.plot_hill_func(sim_run=SIM.mod_runs[1])
# print("Plotting Modulation Hill...")
# SIM.plot_hill_func(sim_run=SIM.mod_runs[0])

print("----")   
print("Done")
print("Writing Output to "+sys.argv[2])
#write_results(sim, sys.argv[2])
print("")
print("Done")