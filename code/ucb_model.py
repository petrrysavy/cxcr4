import numpy as np

from ucb import UCB
from settings import *
from mutation_utils import mutation_at


class UCBmodel:
    def __init__(self, sequence_len, cxcr4sequence_raw):
        self._positionUCB = UCB(sequence_len)
        self._nucleotidesUCB = np.array([UCB(n_arms=Settings().get(SettingsKeys.AMINO_ACIDS_NUM)) for _ in range(sequence_len)])

    def mutate(self, cxcr4sequence):
        position = self._positionUCB.pick_action()
        new = Settings().get(SettingsKeys.AMINO_ACIDS)[self._nucleotidesUCB[position].pick_action()]
        return mutation_at(cxcr4sequence, position, new)

    def register_reward(self, position, new, reward):
        self._positionUCB.register_reward(position, reward)
        acid_index = Settings().get(SettingsKeys.AMINO_ACIDS).index(new)
        self._nucleotidesUCB[position].register_reward(acid_index, reward)