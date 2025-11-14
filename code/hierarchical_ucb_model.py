import numpy as np

from ucb import UCB
from settings import *
from mutation_utils import mutation_at


class HierarchicalUCBmodel:
    def __init__(self, cxcr4sequence_raw):
        self._num_tested = 0
        self._nucleotidesUCB = np.array([UCB(n_arms=Settings().get(SettingsKeys.AMINO_ACIDS_NUM),
                                             exclude_actions=[Settings().get(SettingsKeys.AMINO_ACIDS).index(cxcr4sequence_raw[i])])
                                         for i in range(len(cxcr4sequence_raw))])

    def mutate(self, cxcr4sequence):
        position = -1
        max_value = -np.inf
        # update UCB's totals
        for i in range(len(self._nucleotidesUCB)):
            self._nucleotidesUCB[i].total_time = self._num_tested
            this_ucb_value = self._nucleotidesUCB[i].get_untested_or_best_ucb_value()
            if this_ucb_value > max_value:
                max_value = this_ucb_value
                position = i

        new = Settings().get(SettingsKeys.AMINO_ACIDS)[self._nucleotidesUCB[position].pick_action()]
        return mutation_at(cxcr4sequence, position, new)

    def get_ucb_value_at_pos(self, pos):
        return self._nucleotidesUCB[pos].get_best_ucb_value()

    @property
    def length(self):
        return len(self._nucleotidesUCB)

    def get_model_at_pos(self, pos):
        return self._nucleotidesUCB[pos]

    def register_reward(self, position, new, reward):
        self._num_tested += 1
        acid_index = Settings().get(SettingsKeys.AMINO_ACIDS).index(new)
        self._nucleotidesUCB[position].register_reward(acid_index, reward)

    def __str__(self):
        ucbs_str = ', '.join(f"{i} : {str(self._nucleotidesUCB[i])}\n" for i in range(len(self._nucleotidesUCB)))
        return f"Hierarchical UCB model({ucbs_str})"
