import numpy as np
from settings import *


def mutation_at(sequence, position, new):
    sequence = list(sequence)
    original = str(sequence[position])
    sequence[position] = new
    new = str(sequence[position])
    return "".join(sequence), position, original, new


def random_mutation_at(sequence, position):
    return mutation_at(sequence, position, np.random.choice(Settings().get(SettingsKeys.AMINO_ACIDS)))


def random_position(sequence):
    return np.random.randint(0, len(sequence))


def hamming_distance(seq1, seq2):
    """Calculate the Hamming distance between two sequences. Coded by ChatGPT."""
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of equal length")

    return sum(el1 != el2 for el1, el2 in zip(seq1, seq2))