"""
Module that defines the ScoringFunction class.
"""
import numpy as np
from scipy.interpolate import interp1d

DEFAULT_DESIRABILITY = [{"x": 0.0, "y": 0.0}, {"x": 1.0, "y": 1.0}]

class ScoringFunction:
    """
    Class that handles the integration of functions used to evaluate the particles/molecules
    in the particle swarm.
    """
    def __init__(self, func, name, description=None, desirability=None, truncate_left=True,
                 truncate_right=True, weight=100, is_mol_func=False):
        """
        :param func: A function that takes either a single RDKit mol object as input or an array
            of particle positions (num_particles, ndim) in the CDDD space as input and outputs a
            single score or an array of num_particles scores respectively. Scoring functions with
            additional arguments should be defined as partial.
        :param name: A unique Name of the scoring function. Used for bookkeeping.
        :param description: A description of the scoring function.
        :param desirability: A list of dictionaries where each dictionary {"x": x, "y": y} defines
            a point on the desirability curve used to scale the output of the scoring function into
            the range [0, 1]. If None, a default desirability curve is used which is linear in the
            range [0, 1].
        :param truncate_left: Flag whether the desirability is truncated on the left side (lowest
            defined x), thus set to the same const value for all smaller x or linearly extapolated.
        :param truncate_right: Flag whether the desirability is truncated on the right side (highest
            defined x), thus set to the same const value for all higher x or linearly extrapolated.
        :param weight: The weight of the scoring function in the combined (weighted average) score
            in a multi-objective optimization.
        :param is_mol_func: Flag that defines if the scoring function expects a RDKit mol object
            (True) or an array of particle positions (False).
        """

        self.func = func
        self.name = name
        self.description = description
        self.weight = weight
        self.is_mol_func = is_mol_func
        self._desirability = desirability or DEFAULT_DESIRABILITY
        self.desirability_function = self._create_desirability_function(
            self._desirability,
            truncate_left=truncate_left,
            truncate_right=truncate_right)

    def _create_desirability_function(self, desirability, truncate_left=True, truncate_right=True):
        """
        Method that returns a function that calculates the desirability score for a given input
        unscaled score. Linearly interpolates between points provided.
        :param desirability: List of dictionaries that define points that lie on the
            desirability curve.
        :param truncate_left: Flag whether the desirability is truncated on the left side
            (lowest defined x), thus set to the same const value for all smaller x or
            linearly extrapolated.
        :param truncate_right: Flag whether the desirability is truncated on the right side
            (highest defined x), thus  set to the same const value for all higher x or linearly
            extrapolated.
        :return: A function that calculates the desirability score for a input unscaled score.
        """
        x = [point['x'] for point in desirability]
        y = [point['y'] for point in desirability]
        assert len(x) == len(y)
        if truncate_left:
            x = [x[0] - 1] + x
            y = [y[0]] + y
        if truncate_right:
            x.append(x[-1] + 1)
            y.append(y[-1])
        return interp1d(x, y, fill_value='extrapolate')

    def __call__(self, input):
        """
        Calling a ScoringFunction instance evaluates the scoring function and rescales the scores
        with respect to the desirability scaling and the weight.
        :param input: Either a RDKit mol object or an array of particle positions
            (num_particles, ndim) in the CDDD space.
        :return:
            unscaled_scores: The unscaled output of the scoring function call.
            scaled_scores: The unscaled score scaled with respect to the desirability curve and
                multiplied by the weight of the function.
            desirability_scores: The unscaled score scaled only with respect to the desirability
                curve.
        """
        if self.is_mol_func:
            unscaled_scores = np.array([self.func(mol) for mol in input])
        else:
            unscaled_scores = self.func(input)
        desirability_scores = self.desirability_function(unscaled_scores)
        scaled_scores = desirability_scores * self.weight

        return unscaled_scores, scaled_scores, desirability_scores

    @classmethod
    def from_dict(cls, dictionary):
        """
        Classmethod to create a ScoringFunction instance from a dictionary defining its parameters.
        :param dictionary: A Dictionary defining the ScoringFunction parameters.
        :return: A ScoringFunction instance.
        """
        name = dictionary['name']
        func = dictionary['function']
        description = dictionary['description']
        desirability = dictionary.get('desirability', None)
        weight = dictionary.get('weight', 100)
        is_mol_func = dictionary.get('is_mol_func', True)
        return cls(func=func,
                   name=name,
                   description=description,
                   desirability=desirability,
                   weight=weight,
                   is_mol_func=is_mol_func)
    def to_dict(self):
        """
        Classmethod to write out a ScoringFunction instance parameters to a dictionary.
        :return: A Dictionary with the parameters of the ScoringFunction instance.
        """
        return {'name': self.name,
                'description': self.description,
                'desirabilty': self._desirability,
                'weight': self.weight}

    def __repr__(self):
        return 'mso.objective.ScoringFunction name={} desirability={} weight={}'.format(
            self.name,
            self._desirability,
            self.weight)
