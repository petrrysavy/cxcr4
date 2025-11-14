from solver import BaseSolver
from hierarchical_ucb_model import HierarchicalUCBmodel
from mutation_utils import hamming_distance
from settings import *


# bounded depth reinforcement learning solver
class BoundedLDHRLSolver(BaseSolver):
    def __init__(self, reward):
        super().__init__(reward)

    def iterate_cache(self):
        for cxcr4, pdockQpositive, timestamp, i in self.reward.positive_evaluator.cache.iterate_cxcr4_pdockq_timestamp():
            j = self.reward.negative_evaluator.cache.search_sequence_by_timestamp(timestamp)

            if j is None:
                continue

            cxcr4_negative, pdockQnegative = self.reward.negative_evaluator.cache.ith_entry(j)

            #assert cxcr4 == cxcr4_negative, f"The positive and negative cachce differ at {i}-th position." # FALSe
            if cxcr4 != cxcr4_negative:
                print(f"mismatched timestamps at {cxcr4}, {cxcr4_negative}, {timestamp}")
                continue

            if hamming_distance(self.cxcr4seq, cxcr4) > 5:
                continue

            yield cxcr4, pdockQpositive, pdockQnegative, self.reward.combine_results(pdockQpositive, pdockQnegative)
    def load_from_cache(self) -> HierarchicalUCBmodel:
        model = HierarchicalUCBmodel(self.cxcr4seq)

        for cxcr4, pdockQpositive, pdockQnegative, reward in self.iterate_cache():
            # enumerate all different positions
            for position, (original, new) in enumerate(zip(self.cxcr4seq, cxcr4)):
                if original != new:
                    model.register_reward(position, new, reward)

        print(model)
        return model

    def solve(self, model=None):
        if model is None: model = HierarchicalUCBmodel(self.cxcr4seq)
        while True:
            noimprove = 0
            depth = 0
            best = None
            bestScore = float('-inf')
            curr = float('-inf')
            current = self.cxcr4seq
            while noimprove < 2 and depth < 5:
                current, position, original, new = model.mutate(current)

                val = self.reward.reward(current, info={"position": position, "original": original, "new": new,
                                                        "type": "pointmutation"})
                model.register_reward(position, new, val)

                if val > bestScore:
                    best = current
                    bestScore = val

                if val > curr:
                    curr = val
                    noimprove = 0
                else:
                    noimprove += 1

                depth += 1

            print("finished episode")
            print(best)

        return best
