# Molecular Swarm Optimization (MSO)

Implementation of the method proposed in the paper "Efficient Multi-Objective Molecular Optimization in a Continuous Latent Space" by Robin Winter, Floriane Montanari, Andreas Steffen, Hans Briem, Frank Noé and Djork-Arné Clevert.<sup>1</sup>

### Dependencies
- python 3
- numpy
- rdkit
- scipy
- [cddd](https://github.com/jrwnter/cddd)

### Installing
```
python setup.py install
```
### Example
As a first simple experiment we will optimize a query molecule with respect to the drug likeness score (QED Bickerton et al.). We will start the optimization from a simple benzene molecule which has a QED score of 0.31.
```python
from mso.optimizer import BasePSOptimizer
from mso.objectives.scoring import ScoringFunction
from mso.objectives.mol_functions import qed_score, heavy_atom_count
from cddd.inference import InferenceModel
infer_model = InferenceModel() # The CDDD inference model used to encode/decode molecular SMILES strings to/from the CDDD space
init_smiles = "c1ccccc1" # SMILES representation of benzene
scoring_functions = [ScoringFunction(qed_score, "qed", is_mol_func=True)]
```
After loading some packages and defning the inference model, starting point and the objective functoin for the optimization we can create an instance of the Particle Swarm Optimizer with only one swarm and 100 particles.
```python
opt = BasePSOptimizer.from_query(
    init_smiles=init_smiles,
    num_part=100,
    num_swarms=1,
    inference_model=infer_model,
    scoring_functions=scoring_functions)
```
Now we can just run the optimization for a few steps.
```python
opt.run(20)
```
The best results are summarized in opt.best_solutions . The optimization history (best solution in each step and swarm) is summarized in opt.best_fitness_history . Most of the times the optimizer should be able to find a solution with a score higher than 0.9 already after a few steps.



### References
[1] Chemical Science, 2019, DOI: 10.1039/C9SC01928F https://pubs.rsc.org/en/content/articlelanding/2019/SC/C9SC01928F#!divAbstract

