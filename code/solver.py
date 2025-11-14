from abc import ABC, abstractmethod

from settings import *
from ioutils import load_single_fasta_sequence

class BaseSolver(ABC):
    def __init__(self, reward):
        self._settings = Settings()
        self._cxcr4seq = load_single_fasta_sequence(self._settings.get_file(SettingsKeys.CXCR4_RAW_SEQ_FILE))
        self._reward = reward

    @property
    def cxcr4seq(self):
        return self._cxcr4seq

    @property
    def reward(self):
        return self._reward

    @abstractmethod
    def solve(self):
        pass
