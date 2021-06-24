#!/python3
"""

Functions to plot epitope scores

"""
from typing import List, Tuple
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import seaborn as sns


def plotEpitopeScores(epitope_data,
                      MutationVector=None,
                      EpitopeVector=None,
                      DSSPVector=None,
                      output_filepath=None,
                      curve_labels=None,
                      method_name="Bepipred 2.0",
                      Verbose=False):

    sns.set_style("darkgrid")

    if Verbose:
        print(epitope_data)
        print(EpitopeVector)

    plt.clf()
    fig, ax = plt.subplots()
    fig.dpi = 600
    fig.set_figwidth(9.6)
    matplotlib.rcParams["savefig.dpi"] = fig.dpi

    Colormap = matplotlib.cm.Set1

    # Plot each sequence score line;
    for e, epitopes in enumerate(epitope_data):
        if isinstance(epitopes, list):
            raise
            # FIXME: TRULY NECESSARY?
            X = list(range(len(epitopes)))
            Y = list(map(lambda x: float(x["Score"]), epitopes))
        else:
            X = list(range(epitopes.shape[0]))
            Y = epitopes["Score"]

        A = np.trapz(Y, dx=1)

        if curve_labels is not None:
            label = f"{curve_labels[e]}: {round(A, 2)}"
            label = curve_labels[e]
        else:
            label = None

        ax.plot(
            X,
            Y,
            label=label,
            color=Colormap(e),
            alpha=0.6,
            # linestyle='dashed'
        )

    plt.ylim(min(Y) - 0.05, max(Y) + 0.02 * len(epitope_data) + 0.05)

    # plt.xlabel("Amino Acid Position")
    plt.xlabel("Posição do Aminoácido", fontsize=16)
    # plt.ylabel("BepiPred Score")
    plt.ylabel(f"Pontuação no {method_name}", fontsize=16)

    # Plot mutation letters;
    if MutationVector:
        for m, Mut in enumerate(MutationVector):
            if Mut:
                if Verbose:
                    print(Mut)

                if DSSPVector is not None:
                    secondary_structure_code = list(DSSPVector)[m][2]

                    dsspY = max(Y) - 0.02
                    ax.text(
                        m,
                        dsspY,
                        s=secondary_structure_code,
                        color="black"
                    )

                for T, Track in enumerate(reversed(epitope_data)):
                    Letter = Track["Residue"][m]

                    mutY = max(Y) + 0.03 * T

                    ax.text(
                        m,
                        mutY,
                        s=Letter,
                        color=Colormap(len(epitope_data) - T - 1),
                        fontsize=14
                    )

    # Plot continuous limit line;
    if EpitopeVector:
        line_segments = simplify_epitope_vector(EpitopeVector)
        if Verbose:
            print(line_segments)
        for (frm, to, E) in line_segments:
            if E:
                # STYLE = "dashed" if E else None
                ax.plot(
                    range(frm, to),
                    [0.502] * (to - frm),
                    color="red",
                    linewidth=1
                )

    # PLOT BASE THRESHOLD LINE
    ax.plot(X, [0.5] * len(X), color="red", linewidth=1)

    if curve_labels:
        ax.legend(loc=8)

    if output_filepath is None:
        plt.show()
    else:
        plt.savefig(output_filepath)


def simplify_epitope_vector(EpitopeVector) -> List[Tuple[int, int, bool]]:
    """

    Converts the vector of length L containing overlapped epitope counts
    for the sequence of length L into a list of segments.



    """
    output: List[Tuple[int, int, bool]] = []
    current = None

    for e, Epitope in enumerate(EpitopeVector):
        if current is None:
            current = (e, 0, Epitope > 0)

        else:
            V = Epitope > 0
            if V != current[2] or e == len(EpitopeVector) - 1:
                (c, _, v) = current
                output.append((c, e - 1, v))
                current = (e, 0, V)

    return output
