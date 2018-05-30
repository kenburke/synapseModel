import pickle
import os
from matplotlib import pyplot as plt

def load_input_pickle(name):
    with open(os.getcwd()+'/input/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)
def save_session_pickle(obj, name):
    with open(os.getcwd()+'/session/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def save_input_pickle(obj, name):
    with open(os.getcwd()+'/input/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_session_pickle(name):
    with open(os.getcwd()+'/session/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)

def save_output_plot(file_loc,filename):
    plt.savefig(file_loc+filename+'.svg', format='svg')

def get_user_params():

    params = load_input_pickle('default')

    print("DEFAULT PARAMS :")

    for key, val in params.items():
        l = (21-len(key))//7
        print("{0}".format(key)+"\t"*l+":\t{0}".format(val))

    print("\nNow we will enter our values")

    inters = [
        'num_syn',
        'num_cav',
        'num_stim',
        'num_trials',
        'num_cav_ratio'
    ]

    bool_keys = [
        'depletion_on',
        'phenom_facil',
        'diffusion_model'
    ]

    for k,v in params.items():

        if k in bool_keys:
            val = input(k+" (use 0/1 for False/True) = ")
            val = bool(int(val))
        elif k=='a_tau':
            val1 = input(k+", first (fast) = ")
            val2 = input(k+", second (slow) = ")
            val = (float(val1),float(val2))
        else:
            val = float(input(k+" = "))

        if k in inters:
            val = int(val)

        params[k] = val

    print("-----")
    print("")
    print("NEW PARAMS :")

    for key, val in params.items():
        l = (21-len(key))//7
        print("{0}".format(key)+"\t"*l+":\t{0}".format(val))

    print("")
    name = input("Enter Filename (no extension) > ")
    save_input_pickle(params,name)
    print("")

    return params

def get_user_modulation():
    '''
    gets relevant input for Simulation.run_modulation:
    - which parameter to modulate
    - what range to modulate it over
    '''

    mod_param = input('Enter parameter name for range modulation > ').lower()
    reg_list = input('Regularly-spaced range (R) or manual list of values (L)? > ').lower()

    if reg_list == 'r':

        lower = float(input('Enter lower value > '))
        upper = float(input('Enter upper value > '))
        length = float(input('Enter length of list > '))

        mod_range = [lower + x*(upper-lower)/(length-1) for x in range(int(length))]

    elif reg_list == 'l':
        mod_range = []

        while True:
            print(mod_range)
            val = input('Enter next value (q to quit) > ')
            if val.lower() == 'q':
                break
            else:
                mod_range.append(float(val))

    return (mod_param,mod_range)

def get_synapse_range():
    '''
    gets relevant input for Simulation._runModel with nonuniform_parameter turned on:
    - which parameter to modulate
    - what range to modulate it over
    '''

    mod_param = input('Enter parameter name to distribute over synapses > ').lower()
    reg_list = input('Regularly-spaced range (R) or manual list of values (L)? > ').lower()

    if reg_list == 'r':

        lower = float(input('Enter lower value > '))
        upper = float(input('Enter upper value > '))
        print("WARNING : For best results ensure synapse number is a multiple of this range...")
        length = float(input('Enter length of list > '))

        mod_range = [lower + x*(upper-lower)/(length-1) for x in range(int(length))]

    elif reg_list == 'l':
        mod_range = []

        while True:
            print(mod_range)
            val = input('Enter next value (q to quit) > ')
            if val.lower() == 'q':
                break
            else:
                if mod_param=='num_cav_ratio' or mod_param=="num_cav":
                    mod_range.append(int(val))
                else:
                    mod_range.append(float(val))

    return (mod_param,mod_range)

def dumpclean(obj):
    if type(obj) == dict:
        for k, v in obj.items():
            if hasattr(v, '__iter__'):
                print(k)
                dumpclean(v)
            else:
                print('%s : %s'.format(k, v))
    elif type(obj) == list:
        for v in obj:
            if hasattr(v, '__iter__'):
                dumpclean(v)
            else:
                print(v)
    else:
        print(obj)
