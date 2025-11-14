import os

import numpy as np
import pandas as pd
from datetime import datetime
import subprocess
import time

from predictions_cache import PredictionsCache
from run_command import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class InteractionEvaluator:
    def __init__(self, fixed_sequence_file, input_sequences_dir, results_directory, predictions_cache: PredictionsCache,
                 other_protein_name, force=False, modified_protein_name=Settings().get(SettingsKeys.PROTEIN_NAME)):
        self._fixed_sequence_file = fixed_sequence_file
        self._predictions_cache = predictions_cache
        self._results_directory = results_directory
        self._input_sequences_dir = input_sequences_dir
        self._other_protein_name = other_protein_name
        self._force = force
        self._modified_protein_name = modified_protein_name

    def score(self, cxcr4_sequence, force=None, run_command=RunCommand.DEFAULT, timestamp=None, info=""):
        force = force if force is not None else self._force
        if not force or (isinstance(force, int) and not isinstance(force, bool)):  # watch out - bool is subclass of int
            results = self._predictions_cache.search_sequence(cxcr4_sequence)
            if results and isinstance(force, int) and len(results) > force:
                return np.random.choice(results)

        # we do not use cached results, so call the command

        # first, store the current_sequence, if it is not already stored
        if not timestamp:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        cxcr4_file = self._input_sequences_dir + os.sep + f"{self._modified_protein_name}-{timestamp}.fasta"
        if not os.path.exists(cxcr4_file):
            record = SeqRecord(Seq(cxcr4_sequence), id="1", description=f"{self._modified_protein_name}_"+timestamp)
            SeqIO.write([record], cxcr4_file, "fasta")

        # now run the command
        command = run_command.format_command(cxcr4_file, self._fixed_sequence_file, self._results_directory)
        start_time = time.time()
        try:
            result = subprocess.run(command, capture_output=True, text=True, check=True, cwd=Settings().get_file(SettingsKeys.SPEED_PPI_DIR))
            end_time = time.time()
        except subprocess.CalledProcessError as e:
            print("Command failed with return code:", e.returncode)
            print("Error output:", e.stderr)
            return None

        # locate the output file
        resultfile = self._results_directory + os.sep + f"{self._modified_protein_name}-{timestamp}-{self._other_protein_name}_metrics.csv"
        values = pd.read_csv(resultfile).iloc[0]
        self._predictions_cache.add_row(cxcr4_sequence, values["num_contacts"], values["avg_if_plddt"], values["pdockq"],
                                        run_command.value, end_time-start_time, timestamp, info=str(info))
        return values["pdockq"]

    @property
    def cache(self):
        return self._predictions_cache





