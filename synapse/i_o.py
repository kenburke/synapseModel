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

def write_results(simulation, output_file):
    """
    Input: 
        - "simulation" object from model.runModel
    """

def get_user_params():
    return

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
