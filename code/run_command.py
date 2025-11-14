from enum import Enum
import os
from settings import *


class RunCommand(Enum):
    DEFAULT = ("UNIREF30-23-DEF", Settings().get_file(SettingsKeys.SPEED_PPI_DIR) + os.sep + "predict_single.sh")
    ORIGINAL = ("UNICLUST30-18-DEF", Settings().get_file(SettingsKeys.SPEED_PPI_DIR) + os.sep + "predict_single_orig.sh")

    def __init__(self, setting_name, command):
        self._value_ = setting_name
        self._command = command

    @property
    def value(self):
        return self._value_

    @property
    def command(self):
        return self._command

    def format_command(self, first_sequence, second_sequence, results_directory):
        return ["bash", self._command, first_sequence, second_sequence, Settings().get(SettingsKeys.HHBLITS_COMMAND), results_directory + os.sep]