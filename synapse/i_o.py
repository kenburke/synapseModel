import pickle

def load_input_pickle(name):
    with open('input/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)

def save_output_pickle(obj, name):
    with open('output/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def save_input_pickle(obj, name):
    with open('input/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_output_pickle(name):
    with open('output/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)

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

    for k,v in params.items():

        if k=='depletion_on':
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
