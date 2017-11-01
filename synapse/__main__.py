import sys
from .io import write_results, load_input_pickle, dumpclean
from .model import runModel

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

else:
    print("Usage: python -m synapse [-P| -N| -I] <input directory> <output file>")
    sys.exit(0)


print("")
dumpclean(params)
print("----")   
sim = runModel(params)
print("----")   
print("")
print("Writing Output to"+sys.argv[2])
write_results(sim, sys.argv[2])
print("")
print("Done")