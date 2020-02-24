"""
Module with scoring functions that take RDKit mol objects as input for scoring.
"""
import warnings
from mso.data import data_dir
import os
import pandas as pd
import numpy as np
from functools import wraps
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
from rdkit.Chem.Descriptors import qed, MolLogP
from rdkit import DataStructs


smarts = pd.read_csv(os.path.join(data_dir, "sure_chembl_alerts.txt"), header=None, sep='\t')[1].tolist()
alert_mols = [Chem.MolFromSmarts(smart) for smart in smarts if Chem.MolFromSmarts(smart) is not None]
    
def check_valid_mol(func):
    """
    Decorator function that checks if a mol object is None (resulting from a non-processable SMILES string)
    :param func: the function to decorate.
    :return: The decorated function.
    """
    @wraps(func)
    def wrapper(mol, *args, **kwargs):
        if mol is not None:
            return func(mol, *args, **kwargs)
        else:
            return 0
    return wrapper

@check_valid_mol
def qed_score(mol):
    """
    Quantitative Drug Likeness (QED)
    :param mol: input mol
    :return: score
    """
    try:
        score = qed(mol)
    except :
        score = 0
    return score

@check_valid_mol
def tan_sim(mol, ref_smiles):
    """
    Calculates the Tanimoto similarity between a input molecule and a refence molecule (SMILES) based on ECFP4
    fingerprints.
    :param mol: input molecule
    :param ref_smiles: reference molecule as SMILES.
    :return: The Tanimoto similarity
    """
    ref_mol = Chem.MolFromSmiles(ref_smiles)
    fp_query = AllChem.GetMorganFingerprint(mol, 2)
    fp_ref = AllChem.GetMorganFingerprint(ref_mol, 2)
    sim = DataStructs.TanimotoSimilarity(fp_ref, fp_query)
    return sim

@check_valid_mol
def substructure_match_score(mol, query, kind="any"):
    """
    :param mol: input molecule
    :param query: A list or a single SMARTS pattern the query is checked against.
    :param kind: "any": input should match one of the queries.  "all": should match all.
    :return: 1 if it matches, 0 if not.
    """
    if not isinstance(query, list):
        query = [query]
    if kind == "any":
        match = np.any([mol.HasSubstructMatch(sub) for sub in query])
    elif kind == "all":
        match = np.all([mol.HasSubstructMatch(sub) for sub in query])
    else:
        raise ValueError("use kind == any or all")

    if match:
        score = 1
    else:
        score = 0
    return score

@check_valid_mol
def sa_score(mol):
    """
    Synthetic acceptability score as proposed by Ertel et al..
    """
    try:
        score = sascorer.calculateScore(mol)
    except:
        score = 0
    return score

@check_valid_mol
def logp_score(mol):
    """
    crippen logP
    """
    score = Chem.Crippen.MolLogP(mol)
    return score

@check_valid_mol
def penalized_logp_score(mol, alpha=1):
    """
    penalized logP score as defined by .
    """
    score = reward_penalized_log_p(mol)
    return alpha * score

@check_valid_mol
def heavy_atom_count(mol):
    """
    Number of heavy atoms in molecule
    """
    hac = Chem.Descriptors.HeavyAtomCount(mol)
    return hac


@check_valid_mol
def molecular_weight(mol):
    """molecular weight"""
    mw = Chem.Descriptors.MolWt(mol)
    return mw


@check_valid_mol
def penalize_long_aliphatic_chains(mol, min_members):
    """
    Score that is 0 for molecules with aliphatic chains longer than min_members.
    """
    query = Chem.MolFromSmarts("[AR0]" + "~[AR0]"*(min_members - 1))
    if mol.HasSubstructMatch(query):
        score = 0
    else:
        score = 1
    return score

@check_valid_mol
def penalize_macrocycles(mol):
    """ 0 for molecules with macrocycles."""
    score = 1
    ri = mol.GetRingInfo()
    for x in ri.AtomRings():
        if len(x) > 8:
            score = 0
            break
    return score

@check_valid_mol
def tox_alert(mol):
    """
    0 if a molecule matches a structural alert as defined by the included list from surechembl.
    """
    if np.any([mol.HasSubstructMatch(alert) for alert in alert_mols]):
        score = 0
    else:
        score = 1
    return score

try:
    import cdddswarm.data.sascorer
    import networkx as nx
    @check_valid_mol
    def reward_penalized_log_p(mol):
        """
        Reward that consists of log p penalized by SA and # long cycles,
        as described in (Kusner et al. 2017). Scores are normalized based on the
        statistics of 250k_rndm_zinc_drugs_clean.smi dataset
        Code taken from implementation of:
        You, Jiaxuan, et al. "Graph Convolutional Policy Network for Goal-Directed
        Molecular Graph Generation." arXiv preprint arXiv:1806.02473 (2018).
        https://github.com/bowenliu16/rl_graph_generation
        """
        # normalization constants, statistics from 250k_rndm_zinc_drugs_clean.smi
        logP_mean = 2.4570953396190123
        logP_std = 1.434324401111988
        SA_mean = -3.0525811293166134
        SA_std = 0.8335207024513095
        cycle_mean = -0.0485696876403053
        cycle_std = 0.2860212110245455

        try:
            log_p = MolLogP(mol)
        except ValueError:
            return 0
        try:
            SA = -sascorer.calculateScore(mol)
        except ZeroDivisionError:
            return 0

        # cycle score
        cycle_list = nx.cycle_basis(nx.Graph(
            Chem.rdmolops.GetAdjacencyMatrix(mol)))
        if len(cycle_list) == 0:
            cycle_length = 0
        else:
            cycle_length = max([len(j) for j in cycle_list])
        if cycle_length <= 6:
            cycle_length = 0
        else:
            cycle_length = cycle_length - 6
        cycle_score = -cycle_length

        normalized_log_p = (log_p - logP_mean) / logP_std
        normalized_SA = (SA - SA_mean) / SA_std
        normalized_cycle = (cycle_score - cycle_mean) / cycle_std

        return normalized_log_p + normalized_SA + normalized_cycle
except:
    warnings.warn("failed to load reward_penalized_log_p score. Consider installing package networkx")
    reward_penalized_log_p = None

fps = np.load(os.path.join(data_dir, "chembl_fps.npy"), allow_pickle=True).item()
@check_valid_mol
def has_chembl_substruct(mol):
    """0 for molecuels with substructures (ECFP2 that occur less often than 5 times in ChEMBL."""
    fp_query = AllChem.GetMorganFingerprint(mol, 1, useCounts=False)
    if np.any([bit not in fps for bit in fp_query.GetNonzeroElements().keys()]):
        return 0
    else:
        return 1