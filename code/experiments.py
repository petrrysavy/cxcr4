import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

import ast
from sklearn.decomposition import PCA

from jinja2 import Template
from scipy.stats import binomtest

from bounded_LD_HRL_solver import BoundedLDHRLSolver
from scoring_function import DifferenceReward
from settings import Settings, SettingsKeys
from ioutils import load_single_fasta_sequence


def plot_per_acid(data, position, directory, native):
    labels = list(data.keys())
    means = [val[0] for val in data.values()]
    stds = [val[1] for val in data.values()]

    x_pos = np.arange(len(labels))
    plt.bar(x_pos, means, yerr=stds, capsize=5, alpha=0.7, color='skyblue')
    plt.xticks(x_pos, labels, rotation=45)
    plt.xlabel("Amino Acid (a)")
    plt.ylabel("Q(a) +- stdev")
    plt.title(f"Q({position+1}, a) [native acid is {native}]")

    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()

    plt.savefig(directory + os.sep + f"CXCR4-acid-qs-{position+1}.png", dpi=300)

    plt.close()


def plot_per_acid_pdock(data, position, directory, native, protein):
    labels = list(data.keys())
    means = [val[0] for val in data.values()]
    stds = [val[1] for val in data.values()]

    x_pos = np.arange(len(labels))
    plt.bar(x_pos, means, yerr=stds, capsize=5, alpha=0.7, color='skyblue')
    plt.xticks(x_pos, labels, rotation=45)
    plt.xlabel("Amino Acid (a)")
    plt.ylabel(f"{protein} pDockQ +- stdev")
    plt.title(f"{protein} pDockQ [native acid is {native}]")

    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()

    plt.savefig(directory + os.sep + f"CXCR4-{protein}-pdockq-{position+1}.png", dpi=300)

    plt.close()


def plot_counts(data, position, directory, native):
    labels = list(data.keys())
    counts = [val for val in data.values()]
    x_pos = np.arange(len(labels))
    plt.bar(x_pos, counts, capsize=5, alpha=0.7, color='skyblue')
    plt.xticks(x_pos, labels, rotation=45)
    plt.xlabel("Amino Acid (a)")
    plt.ylabel("n(p,a)")
    plt.title(f"n({position+1}, a) [native acid is {native}]")

    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()

    plt.savefig(directory + os.sep + f"CXCR4-acid-ns-{position+1}.png", dpi=300)

    plt.close()


def plot_per_position(means, directory):
    x_pos = np.arange(len(means)) + 1
    plt.bar(x_pos, means, capsize=5, alpha=0.7, color='skyblue')
    plt.xlabel("Position (1 to 353)")
    plt.ylabel("max_a Q(p,a)")
    plt.title(f"Q-values per position")

    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()

    plt.savefig(directory + os.sep + f"CXCR4-position-qs.png", dpi=300)

    plt.close()


def draw_ucb_criterion():
    solver = BoundedLDHRLSolver(DifferenceReward())
    model = solver.load_from_cache()
    print("loaded from cache")

    AMINO_ACIDS = Settings().get(SettingsKeys.AMINO_ACIDS)
    AMINO_ACIDS_NUM = Settings().get(SettingsKeys.AMINO_ACIDS_NUM)
    CXCR4SEQ = load_single_fasta_sequence(Settings().get_file(SettingsKeys.CXCR4_RAW_SEQ_FILE))

    seqLen = model.length
    ucb_vals = np.zeros(seqLen)
    for position in range(seqLen):
        ucb_submodel = model.get_model_at_pos(position)
        ucb_values_submodel = {AMINO_ACIDS[i]: (ucb_submodel.q_est(i), ucb_submodel.stdev(i)) for i in range(AMINO_ACIDS_NUM)}
        ucb_values_submodel[CXCR4SEQ[position]] = (np.nan, np.nan)
        ucb_vals[position] = ucb_submodel.get_greedy_q_value()
        nsa = {AMINO_ACIDS[i]: (ucb_submodel.n_trials(i)) for i in range(AMINO_ACIDS_NUM)}

        plot_per_acid(ucb_values_submodel, position, Settings().get_file(SettingsKeys.PLOTS_DIR) + os.sep + "ucb_scores", CXCR4SEQ[position])
        plot_counts(nsa, position, Settings().get_file(SettingsKeys.PLOTS_DIR) + os.sep + "ucb_scores", CXCR4SEQ[position])

    plot_per_position(ucb_vals, Settings().get_file(SettingsKeys.PLOTS_DIR) + os.sep + "ucb_scores")

    top_50_indices = sorted(range(len(ucb_vals)), key=lambda i: ucb_vals[i], reverse=True)[:50]
    print(f"Best indices: {[i + 1 for i in top_50_indices]}")

    print(f"Best UCB value: {np.max(ucb_vals)}")
    print(f"Happening at: {np.argmax(ucb_vals)+1}")

    # Create list of (index, value) pairs, 1-based index
    indexed = [(i + 1, val) for i, val in enumerate(ucb_vals)]
    # Sort by value descending
    indexed.sort(key=lambda x: x[1], reverse=True)
    # Write to file
    with open(Settings().get_file(SettingsKeys.PLOTS_DIR) + os.sep + "ucb_scores" + os.sep + "per_position.csv", "w") as f:
        for idx, val in indexed:
            f.write(f"{idx}, {val}, {Settings().get(SettingsKeys.AMINO_ACIDS)[model.get_model_at_pos(idx - 1).get_best_arm()]}\n")
    

def format_pval_latex(p, threshold=0.01, digits=2):
    if p > threshold:
        return f"${p:.{digits}f}$"
    exponent = int(f"{p:.0e}".split('e')[1])
    base = p / 10**exponent
    return r"$" + f"{base:.{digits}f} \\cdot 10^{{{exponent}}}" + r"$"


def print_top_k_sequences():
    CXCR4SEQ = load_single_fasta_sequence(Settings().get_file(SettingsKeys.CXCR4_RAW_SEQ_FILE))

    solver = BoundedLDHRLSolver(DifferenceReward())
    df = pd.DataFrame(solver.iterate_cache(), columns=['sequence', 'pdockpos', 'pdockneg', 'reward'])
    df = df.sort_values(by="reward", ascending=False)

    latex_template = r"""
    \begin{tabular}{cccc}
    \hline
    Mutations & $\frac{\text{pDockQ}}{\text{SDF-1}}$ & $\frac{\text{pDockQ}}{\text{GP120}}$ & $r(s)$ \\
    \hline
    {% for item in data %}
    {{ item.substitutions }} & {{ "{:.3f}".format(item.pdockpos) }} & {{ "{:.3f}".format(item.pdockneg) }} & {{ "{:.3f}".format(item.reward) }} \\
    {% endfor %}
    \hline
    \end{tabular}
    """

    def describe_mutations(original: str, new: str) -> str:
        changes = []
        for i, (orig_char, new_char) in enumerate(zip(original, new), start=1):
            if orig_char != new_char:
                #changes.append(f"${i}: {orig_char} \\rightarrow {new_char}$")
                changes.append(f"{orig_char}{i}{new_char}")
        return ", ".join(changes)

    df['substitutions'] = df.apply(lambda row: describe_mutations(CXCR4SEQ, row['sequence']), axis=1)

    TOP_K = [10, 25, 100]

    template = Template(latex_template)
    for num in TOP_K:
        rendered = template.render(data=df.iloc[:num].to_dict(orient="records"))
        with open(Settings().get_file(SettingsKeys.PLOTS_DIR) + os.sep + "sorted" + os.sep + f"top-{num}.tex", "w") as f:
            f.write(rendered)
    #df.to_excel(Settings().get_file(SettingsKeys.PLOTS_DIR) + os.sep + "sorted" + os.sep + f"results.xlsx", index=False)
    df.to_csv(Settings().get_file(SettingsKeys.PLOTS_DIR) + os.sep + "sorted" + os.sep + f"results.csv", index=False)

    REGIONS = [("N-terminal", 1, 34), ("ECL1", 100, 104), ("ECL2", 174, 192), ("ECL3", 267, 273)]  # by https://pmc.ncbi.nlm.nih.gov/articles/PMC3074590/
    HOTSPOTS = [116,97,98,99,100,101,113,116,117,166,169,171,178,180,182,183,184,185,187,188,189,191,192,193,196,199,200,202,203,
                272,26,31,32,34,265,272,279,285,29,33,266,275,276,278,281,27,30,38,41,267,268,277,282,288,290,259,273,280,284,259,273,284,287,291,292,28]
    HOTSPOTS = set([i-1 for i in HOTSPOTS])
    for alternative in ["greater", "less"]:
        rownames = [i[0] for i in REGIONS] + ["other", "hotspots"]
        regionsdf = pd.DataFrame(np.zeros((len(REGIONS) + 2, len(TOP_K))), index=rownames, columns=[f"$k={i}$" for i in TOP_K]).astype(int)
        hotspots = np.zeros(len(TOP_K))

        regionsdf['start'] =  [r[1]  for r in REGIONS] + ["-", "-"]
        regionsdf['end'] = [r[2]  for r in REGIONS] + ["-", "-"]
        for num in TOP_K:
            for seq in df.iloc[:num]["sequence"]:
                for i, (orig_char, new_char) in enumerate(zip(CXCR4SEQ, seq), start=1):
                    if orig_char != new_char:
                        inregion = False
                        for region, start, end in REGIONS:
                            start -= 1
                            end -= 1  # to natural ordering!
                            if i >= start and i <= end:
                                regionsdf.loc[region, f"$k={num}$"] += 1
                                inregion = True
                        if i in HOTSPOTS:
                            regionsdf.loc["hotspots", f"$k={num}$"] += 1
                        if not inregion:
                            regionsdf.loc["other", f"$k={num}$"] += 1
        columns = ["start", "end"]
        for num in TOP_K:
            regionsdf[f"$\\frac{{p\\mathrm{{-value}}}}{{k={num}}}$"] = "-"
            for region, start, end in REGIONS:
                #pval = hypergeom.sf(regionsdf.loc[region, f"$k={num}$"] - 1, 352, end-start+1, sum(regionsdf[f"$k={num}$"]))
                pval = binomtest(regionsdf.loc[region, f"$k={num}$"], sum(regionsdf[f"$k={num}$"][:-1]), (end - start + 1) / 352, alternative=alternative).pvalue
                regionsdf.loc[region, f"$p\\mathrm{{-value}}_{{k={num}}}$"] = format_pval_latex(pval)
            pval = binomtest(regionsdf.loc["hotspots", f"$k={num}$"], sum(regionsdf[f"$k={num}$"][:-1]), (len(hotspots)) / 352, alternative=alternative).pvalue
            regionsdf.loc["hotspots", f"$p\\mathrm{{-value}}_{{k={num}}}$"] = format_pval_latex(pval)
            columns += [f"$k={num}$", f"$p\\mathrm{{-value}}_{{k={num}}}$"]
        regionsdf = regionsdf[columns]
        latex = regionsdf.to_latex(index=True, column_format="lcc|cc|cc|cc")
        with open(Settings().get_file(SettingsKeys.PLOTS_DIR) + os.sep + "sorted" + os.sep + ("regions.tex" if alternative == "greater" else "regions-lower.tex"), "w") as f:
            f.write(latex)




# this method is work of Adéla Drahokoupilová
def draw_pca():
    solver = BoundedLDHRLSolver(DifferenceReward())
    df = pd.DataFrame(solver.iterate_cache(), columns=['sequence', 'pdockpos', 'pdockneg', 'reward'])
    print(df)

    PROPERTIES = {
        "A": [1.28, 0.05, 1.00, 0.31, 6.11, 0.42, 0.23],
        "G": [0.00, 0.00, 0.00, 0.00, 6.07, 0.13, 0.15],
        "V": [3.67, 0.14, 3.00, 1.22, 6.02, 0.27, 0.49],
        "L": [2.59, 0.19, 4.00, 1.70, 6.04, 0.39, 0.31],
        "I": [4.19, 0.19, 4.00, 1.80, 6.04, 0.30, 0.45],
        "F": [2.94, 0.29, 5.89, 1.79, 5.67, 0.30, 0.38],
        "Y": [2.94, 0.30, 6.47, 0.96, 5.66, 0.25, 0.41],
        "W": [3.21, 0.41, 8.08, 2.25, 5.94, 0.32, 0.42],
        "T": [3.03, 0.11, 2.60, 0.26, 5.60, 0.21, 0.36],
        "S": [1.31, 0.06, 1.60, -0.04, 5.70, 0.20, 0.28],
        "R": [2.34, 0.29, 6.13, -1.01, 10.74, 0.36, 0.25],
        "K": [1.89, 0.22, 4.77, -0.99, 9.99, 0.32, 0.27],
        "H": [2.99, 0.23, 4., 0.13, 7.69, 0.27, 0.30],
        "D": [1.60, 0.11, 2.78, -0.77, 2.95, 0.25, 0.20],
        "E": [1.56, 0.15, 3.78, -0.64, 3.09, 0.42, 0.21],
        "N": [1.60, 0.13, 2.95, -0.60, 6.52, 0.21, 0.22],
        "Q": [1.56, 0.18, 3.95, -0.22, 5.65, 0.36, 0.25],
        "M": [2.35, 0.22, 4.43, 1.23, 5.71, 0.38, 0.32],
        "P": [2.67, 0.00, 2.72, 0.72, 6.80, 0.13, 0.34],
        "C": [1.77, 0.13, 2.43, 1.54, 6.35, 0.17, 0.41]
    }

    def encode_sequence(seq, max_length=352):
        encoded = [PROPERTIES[aa] for aa in seq if aa in PROPERTIES]
        encoded = np.array(encoded)
        encoded.tolist()
        return np.array(encoded)

    df["encoded"] = np.array(df["sequence"].apply(lambda seq: encode_sequence(seq) if isinstance(seq, str) else None))

    encodings = []
    for item in df["encoded"]:
        encodings.append(item.flatten())

    X = np.stack(encodings)
    pca = PCA(n_components=3)
    X2D = pca.fit_transform(X)

    explained_var = pca.explained_variance_ratio_[:2]
    print(f"Explained variance: PC1 = {explained_var[0] * 100:.1f}%, PC2 = {explained_var[1] * 100:.1f}%")

    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(X2D[:, 0], X2D[:, 1], c=df["pdockpos"], cmap='viridis', alpha=1)
    plt.colorbar(scatter, label="pDockQ")
    plt.title("PCA of embedded CXCR4 sequences colored by SDF-1 pDockQ")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.grid(True)
    plt.savefig(Settings().get_file(SettingsKeys.PLOTS_DIR) + os.sep + "pca" + os.sep + f"CXCR4-SDF1.png", dpi=300)
    plt.show()

    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(X2D[:, 0], X2D[:, 1], c=df["pdockneg"], cmap='viridis', alpha=1)
    plt.colorbar(scatter, label="pDockQ")
    plt.title("PCA of embedded CXCR4 sequences colored by GP120 pDockQ")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.grid(True)
    plt.savefig(Settings().get_file(SettingsKeys.PLOTS_DIR) + os.sep + "pca" + os.sep + f"CXCR4-GP120.png", dpi=300)
    plt.show()

    # End of Adéla's code
    plt.figure(figsize=(8, 6))
    plt.bar(range(1, len(pca.explained_variance_ratio_) + 1), pca.explained_variance_ratio_)
    plt.xlabel("Principal Component")
    plt.ylabel("Explained Variance Ratio")
    plt.title("Explained Variance by Component")
    plt.savefig(Settings().get_file(SettingsKeys.PLOTS_DIR) + os.sep + "pca" + os.sep + f"pca-variance.png", dpi=300)
    plt.show()



def main():
    print_top_k_sequences()
    draw_ucb_criterion()
    draw_pca()


if __name__ == '__main__':
    main()

