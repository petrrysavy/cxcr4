import os
from enum import Enum


class SettingsKeys(Enum):
    PROTEIN_NAME = 'PROTEIN_NAME'
    POSITIVE_CONNECTION_PROTEIN = 'POSITIVE_CONNECTION_PROTEIN'
    NEGATIVE_CONNECTION_PROTEIN = 'NEGATIVE_CONNECTION_PROTEIN'
    BASE_DIR = 'BASE_DIR'
    INPUT_SEQUENCES_DIR = 'INPUT_SEQUENCES_DIR'
    RESULTS_DIR = 'RESULTS_DIR'
    PLOTS_DIR = 'PLOTS_DIR'
    SPEED_PPI_DIR = 'SPEED_PPI_DIR'
    HHBLITS_COMMAND = 'HHBLITS_COMMAND'
    PREDICTIONS_CACHE_DIR = "PREDICTIONS_CACHE_DIR"  # relative to base dir
    POSITIVE_PREDICTIONS_CACHE_FILE = 'POSITIVE_PREDICTIONS_CACHE_FILE'  # relative to base dir/predictions cache dir
    NEGATIVE_PREDICTIONS_CACHE_FILE = 'NEGATIVE_PREDICTIONS_CACHE_FILE'
    RAW_SEQUENCES_DIR = "RAW_SEQUENCES_DIR"  # relative to base dir
    POSITIVE_RAW_SEQ_FILE = 'POSITIVE_RAW_SEQ_FILE'  # relative to base dir/raw sequences dir
    NEGATIVE_RAW_SEQ_FILE = 'NEGATIVE_RAW_SEQ_FILE'
    CXCR4_RAW_SEQ_FILE = "CXCR4_RAW_SEQ_FILE"
    CACHE_BACKUP_THRESHOLD = 'CACHE_BACKUP_THRESHOLD'
    AMINO_ACIDS_DICT = 'AMINO_ACIDS_DICT'
    AMINO_ACIDS = 'AMINO_ACIDS'
    AMINO_ACIDS_NUM = 'AMINO_ACIDS_NUM'
    DEBUG = 'DEBUG'
    HARD_EVAL_REPEATS = "HARD_EVAL_REPEATS"

    # CD4 - experiment
    CD4_RAW_SEQ_FILE = "CD4_RAW_SEQ_FILE"
    CD4_SERINE_POSITION = "CD4_SERINE_POSITION"
    OKCT4_RAW_SEQ_FILE_1 = "OKCT4_RAW_SEQ_FILE_1"
    HLADP_RAW_SEQ_FILE = "HLADP_RAW_SEQ_FILE"
    HLADQB1_RAW_SEQ_FILE = "HLADQB1_RAW_SEQ_FILE"
    HLADRB1_RAW_SEQ_FILE = "HLADRB1_RAW_SEQ_FILE"
    OKCT4_CACHE_FILE_1 = "OKCT4_CACHE_FILE_1"
    HLADP_CACHE_FILE = "HLADP_CACHE_FILE"
    HLADQB1_CACHE_FILE = "HLADQB1_CACHE_FILE"
    HLADRB1_CACHE_FILE = "HLADRB1_CACHE_FILE"
    CD4_PROTEIN = "CD4_PROTEIN"
    OKCT4A_PROTEIN = "OKCT4A_PROTEIN"
    HLADP_PROTEIN = "HLADP_PROTEIN"
    HLADQB1_PROTEIN = "HLADQB1_PROTEIN"
    HLADRB1_PROTEIN = "HLADRB1_PROTEIN"


class Settings:
    """
    This class works as settings for the application. It is a singleton class. Initial version coded by ChatGPT.
    """
    _instance = None
    _settings = {}

    def __new__(cls, **kwargs):
        if cls._instance is None:
            cls._instance = super(Settings, cls).__new__(cls)
            cls._instance._load_defaults()
            # Update with user-provided settings on first creation
            for key, value in kwargs.items():
                cls._instance.set(key, value)
        return cls._instance

    def _load_defaults(self):
        """Load default settings"""
        self.set(SettingsKeys.BASE_DIR, "/" + "home" + os.sep + "rysavpe1" + os.sep + "it")

        self.set(SettingsKeys.SPEED_PPI_DIR, "SpeedPPI")
        self.set(SettingsKeys.HHBLITS_COMMAND,
                 self.get(SettingsKeys.BASE_DIR) + os.sep + "hhsuite" + os.sep + "bin" +os.sep + "hhblits")

        self.set(SettingsKeys.PREDICTIONS_CACHE_DIR, "predictionscache")
        self.set(SettingsKeys.RAW_SEQUENCES_DIR, "rawsequences")
        self.set(SettingsKeys.INPUT_SEQUENCES_DIR, "inputtmp")
        self.set(SettingsKeys.RESULTS_DIR, "results")
        self.set(SettingsKeys.PLOTS_DIR, "plots")

        self.set(SettingsKeys.POSITIVE_PREDICTIONS_CACHE_FILE, "sdf1cache.csv")
        self.set(SettingsKeys.NEGATIVE_PREDICTIONS_CACHE_FILE, "gp120cache.csv")

        self.set(SettingsKeys.POSITIVE_RAW_SEQ_FILE, "SDF-1.fasta")
        self.set(SettingsKeys.NEGATIVE_RAW_SEQ_FILE, "GP120.fasta")
        self.set(SettingsKeys.CXCR4_RAW_SEQ_FILE, "CXCR4.fasta")

        self.set(SettingsKeys.PROTEIN_NAME, "CXCR4")
        self.set(SettingsKeys.POSITIVE_CONNECTION_PROTEIN, "SDF-1")
        self.set(SettingsKeys.NEGATIVE_CONNECTION_PROTEIN, "GP120")
        self.set(SettingsKeys.AMINO_ACIDS_DICT, {
            'A': 'Alanine',
            'C': 'Cysteine',
            'D': 'Aspartic Acid',
            'E': 'Glutamic Acid',
            'F': 'Phenylalanine',
            'G': 'Glycine',
            'H': 'Histidine',
            'I': 'Isoleucine',
            'K': 'Lysine',
            'L': 'Leucine',
            'M': 'Methionine',
            'N': 'Asparagine',
            'P': 'Proline',
            'Q': 'Glutamine',
            'R': 'Arginine',
            'S': 'Serine',
            'T': 'Threonine',
            'V': 'Valine',
            'W': 'Tryptophan',
            'Y': 'Tyrosine'
        })
        self.set(SettingsKeys.AMINO_ACIDS, list(self.get(SettingsKeys.AMINO_ACIDS_DICT).keys()))
        self.set(SettingsKeys.AMINO_ACIDS_NUM, len(self.get(SettingsKeys.AMINO_ACIDS)))
        self.set(SettingsKeys.DEBUG, True)
        self.set(SettingsKeys.CACHE_BACKUP_THRESHOLD, 1)  # though implemented, it is probably better not to waste a single prediction

        self.set(SettingsKeys.HARD_EVAL_REPEATS, 50)

        # CD4 experiment
        self.set(SettingsKeys.CD4_RAW_SEQ_FILE, "CD4.fasta")
        self.set(SettingsKeys.CD4_SERINE_POSITION, 81)
        self.set(SettingsKeys.OKCT4_RAW_SEQ_FILE_1, "OKCT4AVH.fasta")
        self.set(SettingsKeys.HLADP_RAW_SEQ_FILE, "HLA-DP.fasta")
        self.set(SettingsKeys.HLADQB1_RAW_SEQ_FILE, "HLA-DQB1.fasta")
        self.set(SettingsKeys.HLADRB1_RAW_SEQ_FILE, "HLA-DRB1.fasta")
        self.set(SettingsKeys.OKCT4_CACHE_FILE_1, "okct4acache.csv")
        self.set(SettingsKeys.HLADP_CACHE_FILE, "hladpcache.csv")
        self.set(SettingsKeys.HLADQB1_CACHE_FILE, "hladqb1cache.csv")
        self.set(SettingsKeys.HLADRB1_CACHE_FILE, "hladrb1cache.csv")
        self.set(SettingsKeys.CD4_PROTEIN, "CD4")
        self.set(SettingsKeys.OKCT4A_PROTEIN, "OKCT4AVH")
        self.set(SettingsKeys.HLADP_PROTEIN, "HLA-DP")
        self.set(SettingsKeys.HLADQB1_PROTEIN, "HLA-DQB1")
        self.set(SettingsKeys.HLADRB1_PROTEIN, "HLA-DRB1")

    def set(self, key: SettingsKeys, value):
        """Set a setting value using enum as key."""
        self._settings[key] = value

    def get(self, key: SettingsKeys, default=None):
        """Retrieve a setting value or return default if not found."""
        return self._settings.get(key, default)

    def get_file(self, key: SettingsKeys):
        if key in (SettingsKeys.PREDICTIONS_CACHE_DIR, SettingsKeys.RAW_SEQUENCES_DIR, SettingsKeys.INPUT_SEQUENCES_DIR, \
                 SettingsKeys.RESULTS_DIR, SettingsKeys.SPEED_PPI_DIR, SettingsKeys.PLOTS_DIR):
            return self.get(SettingsKeys.BASE_DIR) + os.sep + self.get(key)
        elif key in (SettingsKeys.POSITIVE_PREDICTIONS_CACHE_FILE, SettingsKeys.NEGATIVE_PREDICTIONS_CACHE_FILE,
                     SettingsKeys.OKCT4_CACHE_FILE_1, SettingsKeys.HLADP_CACHE_FILE, SettingsKeys.HLADQB1_CACHE_FILE,
                     SettingsKeys.HLADRB1_CACHE_FILE):
            return self.get_file(SettingsKeys.PREDICTIONS_CACHE_DIR) + os.sep + self.get(key)
        elif key in (SettingsKeys.POSITIVE_RAW_SEQ_FILE, SettingsKeys.NEGATIVE_RAW_SEQ_FILE, SettingsKeys.CXCR4_RAW_SEQ_FILE,
                     SettingsKeys.CD4_RAW_SEQ_FILE, SettingsKeys.OKCT4_RAW_SEQ_FILE_1, SettingsKeys.HLADP_RAW_SEQ_FILE,
                     SettingsKeys.HLADQB1_RAW_SEQ_FILE, SettingsKeys.HLADRB1_RAW_SEQ_FILE):
            return self.get_file(SettingsKeys.RAW_SEQUENCES_DIR) + os.sep + self.get(key)

    def load_from_env(self, key: SettingsKeys, default=None):
        """Load a setting from environment variables."""
        return os.getenv(key.value, self.get(key, default))
