from datetime import datetime

from abc import ABC, abstractmethod
from interaction_evaluator import *
from predictions_cache import *


class BaseReward(ABC):
    def __init__(self):
        s = Settings()
        positive_cache = PredictionsCache(s.get_file(SettingsKeys.POSITIVE_PREDICTIONS_CACHE_FILE))
        self.positive_evaluator = InteractionEvaluator(s.get_file(SettingsKeys.POSITIVE_RAW_SEQ_FILE),
                                                       s.get_file(SettingsKeys.INPUT_SEQUENCES_DIR),
                                                       s.get_file(SettingsKeys.RESULTS_DIR),
                                                       positive_cache,
                                                       s.get(SettingsKeys.POSITIVE_CONNECTION_PROTEIN))
        negative_cache = PredictionsCache(s.get_file(SettingsKeys.NEGATIVE_PREDICTIONS_CACHE_FILE))
        self.negative_evaluator = InteractionEvaluator(s.get_file(SettingsKeys.NEGATIVE_RAW_SEQ_FILE),
                                                       s.get_file(SettingsKeys.INPUT_SEQUENCES_DIR),
                                                       s.get_file(SettingsKeys.RESULTS_DIR),
                                                       negative_cache,
                                                       s.get(SettingsKeys.NEGATIVE_CONNECTION_PROTEIN))

    def reward(self, cxcr4_sequence, info="", force=None):
        # provide timestamp, so that we do not run duplicate alignment (=small speedup)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        positive_score = self.positive_evaluator.score(cxcr4_sequence, timestamp=timestamp, info=info, force=force)
        negative_score = self.negative_evaluator.score(cxcr4_sequence, timestamp=timestamp, info=info, force=force)
        return self.combine_results(positive_score, negative_score)

    @abstractmethod
    def combine_results(self, positive_score, negative_score):
        pass


class DifferenceReward(BaseReward):
    def combine_results(self, positive_score, negative_score):
        return positive_score - negative_score
