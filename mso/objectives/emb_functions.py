"""
Modeule with scoring functions that take molecular CDDD embeddings (positions of the particles in the swarm) as input.
"""
from scipy.spatial.distance import cdist
from mso.data import load_predict_model_from_pkl

bace_score_512 = load_predict_model_from_pkl('bace_classifier.pkl')
egfr_score_512 = load_predict_model_from_pkl('egfr_classifier.pkl')


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