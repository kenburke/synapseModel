import pickle

def save_pickle(obj, name):
    with open('input/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_pickle(name):
    with open('input/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)
