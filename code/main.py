from bounded_LD_HRL_solver import BoundedLDHRLSolver
from scoring_function import DifferenceReward


def main():
    solver = BoundedLDHRLSolver(DifferenceReward())
    model = solver.load_from_cache()
    solver.solve(model=model)


if __name__ == '__main__':
    main()
