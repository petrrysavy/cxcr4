import pandas as pd
import os
import shutil
from datetime import datetime

from settings import *
import numpy as np


class PredictionsCache:
    """This class serves as a cache for predictions. The original version was coded by ChatGPT using this propmt:

    <i>
    Hi, write me a python class that will serve as a cache for a csv file and which will store a dataframe.
    The first column of the csv will be named cxcr4 and will contain a sequence. Next columns will be named
    num_contacts,avg_if_plddt, pdockq, and setting. num_contacts is integer, the next two will be double,
    setting is string. The class needs to contain method for searching for a cxcr4 sequence, which then returns
    the pdfckq value list (sequence does not need to be unique). Next method will allow for adding a new row
    into the data frame. Then, there will be methods for iterating over cxcr4 and pdockq pairs. Once the class
    is instantiated, the underlying file is read from hard drive. After 10 new records, the original file is
    copied to a one with extension ".backup" with a timestamp, and the original contents as well as the new records
    are stored in the original file.
    </i>
    """

    def __init__(self, filepath, backup_threshold=Settings().get(SettingsKeys.CACHE_BACKUP_THRESHOLD)):
        self._filepath = filepath
        self._df = None
        self._load_data()
        self._new_records = 0  # Track new records added
        self._backup_threshold = backup_threshold  # Backup after 10 new records

    def search_sequence(self, cxcr4_sequence):
        """Search for a cxcr4 sequence and return the pdockq value(s)."""
        result = self._df[self._df['cxcr4'] == cxcr4_sequence]['pdockq']
        return result.tolist() if not result.empty else []

    def _load_data(self):
        # Load existing data or create empty DataFrame if file doesn't exist
        if os.path.exists(self._filepath):
            self._df = pd.read_csv(self._filepath)
        else:
            self._df = pd.DataFrame(columns=["cxcr4", "num_contacts", "avg_if_plddt", "pdockq", "setting", "time", "timestamp", "info"])

        # Ensure correct data types
        self._df["num_contacts"] = self._df["num_contacts"].astype(int)
        self._df["avg_if_plddt"] = self._df["avg_if_plddt"].astype(float)
        self._df["pdockq"] = self._df["pdockq"].astype(float)
        self._df["time"] = self._df["time"].astype(int)

    def add_row(self, cxcr4, num_contacts, avg_if_plddt, pdockq, setting, time, timestamp, info=""):
        """Add a new row to the DataFrame."""
        new_row = {
            'cxcr4': cxcr4,
            'num_contacts': num_contacts,
            'avg_if_plddt': avg_if_plddt,
            'pdockq': pdockq,
            'setting': setting,
            'time': time,
            'timestamp': timestamp,
            'info': info
        }
        self._df = self._df.append(new_row, ignore_index=True)
        self._new_records += 1

        # Check if backup is needed
        if self._new_records >= self._backup_threshold:
            self._backup_and_save()

    def iterate_cxcr4_pdockq(self):
        """Iterate over cxcr4 and pdockq pairs."""
        for i, row in self._df[['cxcr4', 'pdockq']].iterrows():
            yield row['cxcr4'], row['pdockq'], i

    def iterate_cxcr4_pdockq_timestamp(self):
        for i, row in self._df[['cxcr4', 'pdockq', 'timestamp']].iterrows():
            yield row['cxcr4'], row['pdockq'], row['timestamp'], i

    def search_sequence_by_timestamp(self, timestamp):
        """Search for a cxcr4 sequence by timestamp and return the index."""
        #matches = self._df.loc[self._df['timestamp'] == timestamp]
        matching = self._df['timestamp'] == timestamp
        idx = np.where(matching)
        return idx[0][0] if len(idx[0]) > 0 else None

    def ith_entry(self, i):
        return self._df.iloc[i]['cxcr4'], self._df.iloc[i]['pdockq']

    def ith_row_full(self, i):
        return self._df.iloc[i]

    def _backup_and_save(self):
        """Backup the original file and save new records."""
        # Create a backup with a timestamp
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        if os.path.exists(self._filepath):
            backup_filepath = f"{self._filepath}.backup_{timestamp}"
            shutil.copy(self._filepath, backup_filepath)

        # Save the current DataFrame to the original file
        # modified manually - load first before storing ...
        new_records = self._df.tail(self._backup_threshold)
        self._load_data()
        self._df = self._df.append(new_records, ignore_index=True)
        # some records might be lost if two threads call this method concurrently, but that's ok

        self._df.to_csv(self._filepath, index=False)

        # Reset the new records counter
        self._new_records = 0
