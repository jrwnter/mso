'''
Author: your name
Date: 2020-08-31 13:28:56
LastEditTime: 2020-09-24 19:16:41
LastEditors: Please set LastEditors
Description: In User Settings Edit
FilePath: \mso\objectives\emb_functions.py
'''
"""
Modeule with scoring functions that take molecular CDDD embeddings (positions of the particles in the swarm) as input.
"""
from scipy.spatial.distance import cdist
from mso.data import load_predict_model_from_pkl

bace_score_512 = load_predict_model_from_pkl('bace_classifier.pkl')
egfr_score_512 = load_predict_model_from_pkl('egfr_classifier.pkl')


def check_valid_smiles(swarm):
    def calculate_score(func):
        @wraps(func)
        def wrapper(swarm, *args, **kwargs):
            bo = np.array([Chem.MolFromSmiles(smi) is not None for smi in swarm.smiles])
            score = func(swarm.x, *args, **kwargs)
            score = np.where(bo, score, -100)
            return score
        return wrapper
    return calculate_score


@check_valid_smiles(swarm=None)
def distance_score(x, target, metric="cosine"):
    """
    Function tha calculates the distance between an input molecular embedding and a target molecular embedding.
    :param x: input molecular embedding
    :param target: target molecular embedding
    :param metric: The metric used by scipy.spatial.distance.cdist to compute the distance in space.
    :return: The distance between input and target.
    """
    score = cdist(x, target, metric).flatten()
    return score