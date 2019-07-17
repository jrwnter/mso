import os
import pickle

data_dir = os.path.dirname(os.path.realpath(__file__))

def get_file(name):
    full_path = os.path.join(data_dir, name)
    if not os.path.exists(full_path):
        return None
    else:
        return full_path


def load_predict_model_from_pkl(filename):
    ''' The pickled models are assumed to have a predict method'''
    pkl = get_file(filename)
    if pkl:
        with open(pkl, 'rb') as pkl:
            return pickle.load(pkl).predict
    else:
        return None

def load_transform_model_from_pkl(filename):
    ''' The pickled models are assumed to have a predict method '''
    pkl = get_file(filename)
    if pkl:
        with open(pkl, 'rb') as pkl:
            return pickle.load(pkl).transform
    else:
        return None
