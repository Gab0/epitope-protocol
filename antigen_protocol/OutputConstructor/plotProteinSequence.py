#!/bin/python

"""

This module converts a polypeptide sequence with mutations into an
informative plot figure.

"""

from typing import List

from collections import namedtuple
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


class MutationTextPositionAlternator():
    """

    Helper class: this object keeps track of the position
    of the alternative AA labels on the text.

    """
    IDX = 0
    X_POS = None
    LoadedBoxes: List[np.ndarray] = []

    def __init__(self,
                 SequenceSegments: List[str],
                 multiplier: float = 1,
                 hclearance=16):
        self.hclearance = hclearance

        assert multiplier > 0

        self.VariableNLabelRows = [
            self.calculate_label_nrows_from_horizontal_clearance(SEG)
            for SEG in SequenceSegments
        ]

        self.Positions = [
            W * multiplier
            for W in range(1, len(self.VariableNLabelRows) + 1)
        ]

    def get_offset(self, X_POS):
        if self.X_POS is None:
            self.X_POS = X_POS

        else:
            DeltaX = abs(X_POS - self.X_POS)
            if DeltaX > self.hclearance:
                self.X_POS = X_POS
                self.IDX = 0
            else:
                self.X_POS = X_POS
                self.IDX += 1
                if self.IDX >= len(self.Positions):
                    self.IDX = 0

        return self.Positions[self.IDX]

    def calculate_label_nrows_from_horizontal_clearance(self, var_residues):
        nrows = 0
        for v, _ in enumerate(var_residues):
            segment_nrows = sum(var_residues[v: v + self.hclearance])
            nrows = max(nrows, segment_nrows)

        return nrows

    def coord_violates_box(self, X, Y):
        clearance = 2

        def inside(V, L, U):
            return L - clearance < V < U + clearance

        for B in self.LoadedBoxes:
            Xv = inside(X, B[0][0], B[1][0])
            Yv = inside(Y, B[0][1], B[1][1])

            if Xv and Yv:
                return True

        return False

    def reboot(self):
        self.IDX = 0


def separate_sequence_segments(MutationVector, XSIZE):

    LabelSegments = []
    for i, MUT in enumerate(MutationVector):
        L = (i // XSIZE) + 1
        if L > len(LabelSegments):
            LabelSegments.append([])

        V = len(list(MUT.variations.keys())) if MUT else 0
        LabelSegments[-1].append(V)

    return LabelSegments


def calculate_optionalXSIZE(SeqLength, bounds=(80, 140)):
    def score(XSIZE, SeqLength):
        W = (SeqLength / XSIZE) % 1
        b = min(abs(1 - W), W) * 10
        c = round(SeqLength / XSIZE)

        return b + c

    scores = {
        b: score(b, SeqLength)
        for b in range(bounds[0], bounds[1] + 1)
    }

    XSIZERank = sorted(scores.keys(), key=lambda k: scores[k])

    return XSIZERank[0]


# -- Mutation functions;
def get_mutation_color(mut):
    W = fetch_percentages(mut)

    MAX = max([w[1] for w in W])

    if mut.SpecialSwap:
        return "Special"
    if MAX < 90:
        return "Common"
    if MAX < 75:
        return "Uncommon"
    if MAX > 10:
        return "Rare"

    return "Unknown."


def fetch_percentages(mut):
    return sorted(
        [
            (altAA, 100 * len(mut.variations[altAA]) / mut.get_nbseq())
            for altAA in mut.variations.keys()
        ],
        key=lambda w: w[-1],
        reverse=True
    )


def fetchVariationText(AA, mut, idx):
    def buildLine(AA, pct):
        return "%s (%.2f%%)" % (AA, pct)

    W = fetch_percentages(mut)

    return ["%i" % idx] + [
        buildLine(*w)
        for w in W
    ]


def calculate_anchor(BBox, top=False):
    BBoxV = BBox.get_points()
    aX = BBoxV.T[0].sum() / 2

    if top:
        aY = max(BBoxV.T[1])
    else:
        aY = min(BBoxV.T[1])

    return (aX, aY)


def drawBoxClassic(ax, Texts, MUT_TPA, pointerXY, hfont, FontSize):
    (pointerX, pointerY) = pointerXY
    for Text in Texts:
        # -- Draw variant label;
        variantPositionY =\
            pointerY + MUT_TPA.get_offset(pointerX)

        ax.text(
            pointerX,
            variantPositionY,
            Text,
            **hfont,
            fontsize=0.7 * FontSize,
            rotation=0
        )

    return variantPositionY


def PadBBox(BBox, padx, pady):
    hX = padx // 2
    hY = pady // 2

    BB = BBox.get_points()

    Padding = [
        [-hX, -hY],
        [hX, hY]
    ]

    BBox.set_points(BB + Padding)
    return BBox


def AxisGetConstantOffset(fig, coord):
    """

    Transforms a coordinate as absolute fraction of the figure
    and returns it as axis values.

    """
    _coord = tuple((fig.get_size_inches()[i] / (7 * fig.get_dpi())) * coord[i] for i in range(2))
    return fig.dpi_scale_trans.transform(_coord)


def drawBoxOrganized(ax, Texts, MUT_TPA, pointer, hfont, FontSize):
    """

    DEPRECATED method of drawing information about a SNP right above
    its position on the AA sequence.

    """

    _, pointerY = pointer
    _, pdY = AxisGetConstantOffset(ax.get_figure(), (0, 0.0014))
    print(pdY)
    variantPositionY = pointerY

    BBoxStyle = dict(facecolor='silver',
                     edgecolor='brown',
                     boxstyle='round',
                     clip_on=True)

    posX0 = posX = 1  # max(1, pointerX - 32)

    posY = variantPositionY

    HT = False
    while True:
        if posX > ax.get_xlim()[-1]:
            posX = posX0
            HT = True

        if posY > ax.get_ylim()[-1]:
            m = "Invalid position for label."
            print(m)
            # raise RuntimeError(m)

        T = ax.text(
            posX,
            posY,
            "\n".join(Texts),
            bbox=BBoxStyle,
            fontsize=0.7 * FontSize,
            **hfont,
            zorder=10
        )

        rawBBox = get_object_bbox(ax, T)
        BBox = PadBBox(rawBBox, 3, 1.4)

        Conflict = False
        for Other in MUT_TPA.LoadedBoxes:
            if BBox.overlaps(Other):
                Conflict = True
                break

        if not Conflict:
            MUT_TPA.LoadedBoxes.append(BBox)
            break

        T.remove()

        if not HT:
            posX += 0.3
        else:
            W = ax.transAxes.transform((0, 0.3))
            _, dY = ax.transData.transform(W)
            posY += dY

    return posY


def get_object_bbox(ax, obj):
    renderer = ax.get_figure().canvas.get_renderer()
    BBox = obj.get_tightbbox(renderer)

    return BBox.transformed(ax.transData.inverted())


def getYPointer(fig, r, MUT_TPA, YRowSpacer=0.04):
    """

    Most used method of calculating the Y pointer
    based on a given line of the aminoacid sequence.

    """
    if r < 0:
        return 0
    try:
        NLBL = MUT_TPA.VariableNLabelRows[r]
    except IndexError:
        NLBL = 5

    # YRow = max(2, NLBL // 5) * 0.6 # 0.25
    # YRow = NLBL * (VariantLabelSpacer * 1.05)
    #_, YRow = AxisGetConstantOffset(fig, (0, 0.07))

    return YRowSpacer + getYPointer(fig, r - 1, MUT_TPA)


def ShowStatistic(s):
    if isinstance(s, int):
        return "%s" % s

    if isinstance(s, float):
        return "%.2f" % s

    raise AttributeError("Invalid statistics.")


SequenceParameters = namedtuple("SequenceParameters",
                                "XSIZE ModelBounds YSpacer")


def InitializeStatistics():
    Statistics = {
        "#AASNP": 0,
        "#AASNPinsideEpitopes": 0,
        "#NonSynMutRatio": 0.0,
        "#AGCR": 0.0,
    }

    PrintableStatistics = {
        "#AASNP": r"#$\Delta$AA",
        "#AASNPinsideEpitopes": r"#$\Delta$AA $\in$ Epitopes",
        "#NonSynMutRatio": "Non Synonymous Mutation Ratio",
        "#AGCR": "AA Group Change Ratio",
        "AASNP/SEQ_NB": "Ratio: $\\Delta$AA / #Ind"
    }
    return Statistics, PrintableStatistics


def ShowSequence(Sequence: str,
                 AlnSize: int,
                 Identifier=None,
                 MutationVector=None,
                 EpitopeVector=None,
                 ModelBounds=None,
                 XSIZE=80,
                 OutputFilePath=None):

    if Sequence.endswith("*"):
        print("Trimming trailing * from AA sequence.")
        Sequence = Sequence.strip("*")

    drawBox = False
    if XSIZE == 0:
        XSIZE = calculate_optionalXSIZE(len(Sequence))
        print("XSIZE calculated as %i for sequence of length %i." %
              (XSIZE, len(Sequence)))

    fig, ax = plt.subplots()

    VariantLabelSpacer = 0.20

    SequenceSegments = separate_sequence_segments(MutationVector, XSIZE)

    MUT_TPA = MutationTextPositionAlternator(
        SequenceSegments,
        multiplier=VariantLabelSpacer
    )

    YSIZE = getYPointer(fig, len(Sequence) // XSIZE, MUT_TPA) + 0.5

    print()
    print("Drawing for %s" % Identifier)
    print("YSIZE = %.2f" % YSIZE)
    print("Variable NRows: %s" % MUT_TPA.VariableNLabelRows)

    AxisBounds = (XSIZE, YSIZE)

    FontSize = 7

    fig.dpi = 600

    matplotlib.rcParams["savefig.dpi"] = fig.dpi
    hfont = {
        'fontname': 'Monospace',
        "color": "black",
        "weight": None
    }

    ModelSpanIndexes: List[int] = []
    if ModelBounds is not None:
        for f, t in ModelBounds:
            ModelSpanIndexes += list(range(f, t))

    # -- Statistics;
    Statistics, PrintableStatistics = InitializeStatistics()

    ax.axis([0, AxisBounds[0], 0, AxisBounds[1]])

    # -- ITERATE AMINOACID RESIDUES;
    for i, AA in enumerate(Sequence):
        AAFont = hfont.copy()
        pointerX = i % AxisBounds[0]
        lineY = i // AxisBounds[0]

        if not pointerX:
            MUT_TPA.reboot()

        pointerY = AxisBounds[1] - getYPointer(fig, lineY, MUT_TPA)
        assert pointerY > 0

        # -- Process Mutation;
        Mut = MutationVector[i]
        inMut = MutationVector is not None and Mut

        MutationHighlightColor = {
            "Rare": "yellowgreen",
            "Uncommon": "darkorange",
            "Common": "peachpuff",
            "Special": "cyan"
        }

        MutationHighlightLegend = {
            "Special": "Residue Group Swap",
            "Common": "Common mutation",
            "Uncommon": "Uncommon mutation",
            "Rare": "Rare mutation"
        }

        MutationHighlightLegend = {
            "Special": "Troca do grupo do resíduo",
            "Common": "Troca comum",
            "Uncommon": "Troca incomum",
            "Rare": "Troca rara"
        }

        MutationDegree =\
            get_mutation_color(MutationVector[i]) if Mut else "Common"

        if MutationDegree == "Special":
            Statistics["#AGCR"] += 1

        AABox = dict(
            facecolor=MutationHighlightColor[MutationDegree],
            edgecolor='none',
            alpha=1,
            pad=0.0
        ) if inMut else None

        # -- Evaluate epitope positions;
        if EpitopeVector is not None:
            InEpitope = EpitopeVector[i]
            if InEpitope:
                if inMut:
                    Statistics["#AASNPinsideEpitopes"] += 1
                AAFont["weight"] = "bold"
                AAFont["color"] = "purple"

        if ModelSpanIndexes:
            if i in ModelSpanIndexes:
                _, offY = AxisGetConstantOffset(fig, (0, 0.015))
                ax.text(pointerX, pointerY - offY, chr(8212), **hfont)

        # -- Draw aminoacid label;
        AAText = ax.text(
            pointerX,
            pointerY,
            AA,
            bbox=AABox,
            **AAFont,
            fontsize=FontSize,
        )

        # -- Draw current Sequence index (periodically);
        if not (i + 1) % 20:
            _, nbpY = AxisGetConstantOffset(fig, (0, 0.025))

            ax.text(
                pointerX,
                pointerY + nbpY,
                str(i + 1),
                fontsize=4
            )

        TextAnchor = calculate_anchor(get_object_bbox(ax, AAText), top=True)

        # -- Evaluate mutations;
        if inMut:
            mut = MutationVector[i]
            if mut:
                inMut = True
                Statistics["#AASNP"] += 1

                Texts = fetchVariationText(AA, mut, i + 1)
                # drawBoxClassic(ax, Texts, MUT_TPA,
                #               (pointerX, pointerY), hfont, FontSize)

            # -- Draw highlight box around aminoacid label;
            if drawBox:
                drawBoxOrganized(
                    ax,
                    Texts,
                    MUT_TPA,
                    TextAnchor,
                    hfont,
                    FontSize
                )

    # -- Calculate additional statistics;
    Statistics["#NonSynMutRatio"] = Statistics["#AASNP"] / len(Sequence)
    Statistics["#AGCR"] = Statistics["#AGCR"] / len(Sequence)
    Statistics["AASNP/SEQ_NB"] = Statistics["#AASNP"] / AlnSize

    # -- Plot statistics;
    InfoText = []
    if Identifier is not None:
        InfoText.append(Identifier)

    for _, Stat in enumerate(Statistics.keys()):
        message = "%s: %s" % (PrintableStatistics[Stat],
                              ShowStatistic(Statistics[Stat]))
        InfoText.append(message)

    YSEP = 0.12
    _, summary_Yspacer = AxisGetConstantOffset(fig, (0, YSEP))
    InfoY = pointerY - summary_Yspacer

    if False:
        InfoText = ax.text(
            0,
            InfoY,
            "\n".join(InfoText),
            fontsize=0.8 * FontSize
        )

    # -- Create plot legend;
    LegendRefs = []
    for MutName, MutColor in MutationHighlightColor.items():
        p = mpatches.Patch(
            color=MutColor,
            label=MutationHighlightLegend[MutName]
        )
        LegendRefs.append(p)

    # -- Define plot legend position;
    _, lgY = AxisGetConstantOffset(fig, (0, 0.016))

    print(lgY)
    XL = ax.get_xlim()
    print(XL)

    LEGBBOX_COORD = (XSIZE - 5, (InfoY + pointerY) * 0.536)

    LEGBBOX_ABS = ax.transData.transform(LEGBBOX_COORD)
    LEGBBOX = ax.transAxes.inverted().transform(LEGBBOX_ABS)
    print("LEGBc %s" % (LEGBBOX_COORD,))
    print("LEGBf %s" % LEGBBOX)

    ax.legend(
        handles=LegendRefs,
        fontsize=0.8 * FontSize,
        loc="upper right",
        bbox_to_anchor=LEGBBOX
    )

    # -- Disable axis X and Y rulers;
    ax.axis('off')

    # plt.tight_layout()

    ShowPlot(OutputFilePath)


def ShowPlot(OutputFilePath):
    if OutputFilePath is None:
        plt.show()

    else:
        plt.savefig(
            OutputFilePath,
            transparent=True,
            bbox_inches='tight',
            pad_inches=1
        )
