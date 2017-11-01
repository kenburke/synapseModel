import sys
from .io import write_results
from .model import runModel

# Some quick stuff to make sure the program is called correctly
if len(sys.argv) < 4:
    print("Usage: python -m synapse [-P| -N| -I] <input directory> <output file>")
    sys.exit(0)
    
# Choose which effect to model

if sys.argv[1][0:2] == '-P':
    print("Modeling Suppression of Open Probability, P"
    traces = runModel(sys.argv[2])
    print("Writing Output to"+sys.argv[3])
    write_results(sys.argv[3], traces)

elif sys.argv[1][0:2] == '-N':
    print("Modeling Suppression of Number of Channels, N"
    traces = runModel(sys.argv[2])
    print("Writing Output to"+sys.argv[3])
    write_results(sys.argv[3], traces)


elif sys.argv[1][0:2] == '-I':
    print("Modeling Suppression of Single-Channel Current, I"
    traces = runModel(sys.argv[2])
    print("Writing Output to"+sys.argv[3])
    write_results(sys.argv[3], traces)

else:
    print("Usage: python -m synapse [-P| -N| -I] <input directory> <output file>")
    sys.exit(0)
