from pathlib import Path
import os
import tempfile

os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", os.path.join(tempfile.gettempdir(), "fontconfig"))
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "images"
DPI = 600

COLORS = {
    "input": "#EAF4F4",
    "process": "#FFF3D6",
    "analysis": "#EDE7F6",
    "output": "#E5EEF9",
    "accent": "#1F5F6B",
    "accent2": "#B8323B",
    "text": "#1D2733",
    "muted": "#52616B",
    "line": "#2F4858",
}


plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 13,
        "axes.titlesize": 17,
        "axes.labelsize": 14,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.dpi": 160,
    }
)


def box(ax, xy, text, width=2.15, height=0.78, color="process", fontsize=13):
    x, y = xy
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.055,rounding_size=0.10",
        linewidth=1.6,
        edgecolor=COLORS["line"],
        facecolor=COLORS[color],
    )
    ax.add_patch(patch)
    ax.text(
        x + width / 2,
        y + height / 2,
        text,
        ha="center",
        va="center",
        fontsize=fontsize,
        color=COLORS["text"],
        linespacing=1.16,
    )
    return (x, y, width, height)


def arrow(ax, start, end):
    sx, sy, sw, sh = start
    ex, ey, ew, eh = end
    patch = FancyArrowPatch(
        (sx + sw, sy + sh / 2),
        (ex, ey + eh / 2),
        arrowstyle="-|>",
        mutation_scale=16,
        linewidth=1.7,
        color=COLORS["line"],
        shrinkA=8,
        shrinkB=8,
    )
    ax.add_patch(patch)


def arrow_down(ax, start, end):
    sx, sy, sw, sh = start
    ex, ey, ew, eh = end
    patch = FancyArrowPatch(
        (sx + sw / 2, sy),
        (ex + ew / 2, ey + eh),
        arrowstyle="-|>",
        mutation_scale=16,
        linewidth=1.7,
        color=COLORS["line"],
        shrinkA=8,
        shrinkB=8,
    )
    ax.add_patch(patch)


def arrow_points(ax, start_xy, end_xy, connectionstyle="arc3,rad=0"):
    patch = FancyArrowPatch(
        start_xy,
        end_xy,
        arrowstyle="-|>",
        mutation_scale=16,
        linewidth=1.7,
        color=COLORS["line"],
        shrinkA=6,
        shrinkB=6,
        connectionstyle=connectionstyle,
    )
    ax.add_patch(patch)


def setup(figsize=(11, 4.6)):
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    ax.set_axis_off()
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    return fig, ax


def save(fig, name):
    path = OUT / name
    fig.savefig(path, dpi=DPI, bbox_inches="tight", facecolor="white", pad_inches=0.08)
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight", facecolor="white", pad_inches=0.08)
    plt.close(fig)


def pipeline_overview():
    fig, ax = setup((16.4, 5.0))
    ax.set_xlim(0, 16.4)
    ax.set_ylim(0, 4.8)

    b1 = box(ax, (0.35, 3.15), "FASTQ\nreads", color="input")
    b2 = box(ax, (0.35, 1.95), "Sample\nsheet", color="input")
    b3 = box(ax, (0.35, 0.75), "Reference\nfiles", color="input")
    qc = box(ax, (3.0, 3.15), "Read QC\nFastQC + MultiQC")
    quant = box(ax, (5.7, 1.95), "Quantification\nSalmon or STAR")
    counts = box(ax, (8.35, 1.95), "Gene count\nmatrix", color="output")
    de = box(ax, (11.15, 1.95), "Differential\nexpression", color="analysis")
    interp = box(ax, (13.9, 1.95), "Plots,\nenrichment,\ninterpretation", color="output", width=2.05, fontsize=11.5)

    arrow(ax, b1, qc)
    arrow(ax, qc, quant)
    arrow(ax, b2, quant)
    arrow(ax, b3, quant)
    arrow(ax, quant, counts)
    arrow(ax, counts, de)
    arrow(ax, de, interp)

    save(fig, "09-rnaseq-pipeline-overview.png")


def quantification_routes():
    fig, ax = setup((16.2, 6.0))
    ax.set_xlim(0, 16.2)
    ax.set_ylim(0, 5.8)

    fastq = box(ax, (0.4, 2.55), "FASTQ\nreads", color="input")
    salmon = box(ax, (3.05, 3.9), "Salmon\ntranscriptome\nquantification", height=0.84)
    quant = box(ax, (5.75, 3.9), "quant.sf\nper sample", color="output")
    txi = box(ax, (8.35, 3.9), "tximport\ngene summary")

    align = box(ax, (3.05, 1.0), "STAR/HISAT2\ngenome\nalignment", height=0.84)
    bam = box(ax, (5.75, 1.0), "BAM files", color="output")
    fc = box(ax, (8.35, 1.0), "featureCounts\ngene counting")
    counts = box(ax, (11.35, 2.55), "Gene count\nmatrix", color="output")
    stats = box(ax, (14.05, 2.55), "DESeq2\nedgeR\nlimma-voom", color="analysis", width=1.55, fontsize=11.5)

    arrow(ax, fastq, salmon)
    arrow(ax, salmon, quant)
    arrow(ax, quant, txi)
    arrow_points(
        ax,
        (txi[0] + txi[2], txi[1] + txi[3] / 2),
        (counts[0], counts[1] + counts[3] * 0.68),
        connectionstyle="arc3,rad=-0.08",
    )
    arrow(ax, fastq, align)
    arrow(ax, align, bam)
    arrow(ax, bam, fc)
    arrow_points(
        ax,
        (fc[0] + fc[2], fc[1] + fc[3] / 2),
        (counts[0], counts[1] + counts[3] * 0.32),
        connectionstyle="arc3,rad=0.08",
    )
    arrow(ax, counts, stats)

    ax.text(3.05, 5.15, "Transcriptome route", fontsize=14, weight="bold", color=COLORS["accent"])
    ax.text(3.05, 0.25, "Genome alignment route", fontsize=14, weight="bold", color=COLORS["accent"])

    save(fig, "09-rnaseq-quantification-routes.png")


def three_inputs():
    fig, ax = setup((12.2, 4.9))
    ax.set_xlim(0, 12.2)
    ax.set_ylim(0, 4.9)

    ref = box(ax, (0.55, 3.35), "Reference\nfiles", color="input")
    reads = box(ax, (0.55, 2.1), "FASTQ\nreads", color="input")
    sheet = box(ax, (0.55, 0.85), "Sample\nsheet", color="input")

    q1 = box(ax, (3.4, 3.35), "What can\nreads be\nassigned to?")
    q2 = box(ax, (3.4, 2.1), "What was\nsequenced,\nand how well?")
    q3 = box(ax, (3.4, 0.85), "What should\nthe model\nknow?")

    matrix = box(ax, (6.55, 2.1), "Reliable\ncount matrix", color="output")
    model = box(ax, (9.45, 2.1), "Interpretable\nresults", color="analysis")

    arrow(ax, ref, q1)
    arrow(ax, reads, q2)
    arrow(ax, sheet, q3)
    arrow(ax, q1, matrix)
    arrow(ax, q2, matrix)
    arrow(ax, q3, matrix)
    arrow(ax, matrix, model)

    save(fig, "09-rnaseq-three-inputs.png")


def deseq2_flow():
    fig, ax = setup((12.4, 5.8))
    ax.set_xlim(0, 12.4)
    ax.set_ylim(0, 5.8)

    counts = box(ax, (0.55, 3.8), "Count\nmatrix", color="input")
    meta = box(ax, (0.55, 1.25), "Sample\nmetadata", color="input")
    dds = box(ax, (3.35, 2.5), "DESeqDataSet\nmodel design")
    transform = box(ax, (6.0, 3.8), "vst/rlog\ntransformed\ncounts", color="analysis")
    qc = box(ax, (8.9, 3.8), "PCA,\nclustering,\nheatmaps", color="output")
    model = box(ax, (6.0, 1.25), "Negative\nbinomial\nmodel", color="analysis")
    results = box(ax, (8.9, 1.25), "results()\nlfcShrink()", color="output")

    arrow(ax, counts, dds)
    arrow(ax, meta, dds)
    arrow(ax, dds, transform)
    arrow(ax, transform, qc)
    arrow(ax, dds, model)
    arrow(ax, model, results)

    ax.text(
        6.05,
        5.05,
        "exploratory visualization",
        fontsize=13,
        weight="bold",
        color=COLORS["accent"],
    )
    ax.text(
        6.05,
        0.45,
        "formal differential expression testing",
        fontsize=13,
        weight="bold",
        color=COLORS["accent"],
    )

    save(fig, "09-rnaseq-deseq2-flow.png")


def qc_decision_tree():
    fig, ax = setup((12.8, 6.0))
    ax.set_xlim(0, 12.8)
    ax.set_ylim(0, 6.0)

    multiqc = box(ax, (0.55, 2.75), "MultiQC\nsummary", color="input")
    reads = box(ax, (3.35, 4.75), "Uneven\nread counts?")
    adapters = box(ax, (3.35, 2.75), "Adapters or\nlow-quality ends?")
    stranded = box(ax, (3.35, 0.75), "Strandedness\nmismatch?")

    action1 = box(ax, (6.45, 4.75), "Check library\nnotes and\nsample swaps", color="output")
    action2 = box(ax, (6.45, 2.75), "Trim or rerun\nwith careful\nread handling", color="output")
    action3 = box(ax, (6.45, 0.75), "Correct Salmon\nor featureCounts\nsettings", color="output")
    rerun = box(ax, (10.0, 2.75), "Rerun QC\nbefore DE", color="analysis")

    arrow(ax, multiqc, reads)
    arrow(ax, multiqc, adapters)
    arrow(ax, multiqc, stranded)
    arrow(ax, reads, action1)
    arrow(ax, adapters, action2)
    arrow(ax, stranded, action3)
    arrow(ax, action2, rerun)
    arrow_points(
        ax,
        (action1[0] + action1[2], action1[1] + action1[3] / 2),
        (rerun[0], rerun[1] + rerun[3] * 0.72),
        connectionstyle="arc3,rad=-0.10",
    )
    arrow_points(
        ax,
        (action3[0] + action3[2], action3[1] + action3[3] / 2),
        (rerun[0], rerun[1] + rerun[3] * 0.28),
        connectionstyle="arc3,rad=0.10",
    )

    save(fig, "09-rnaseq-qc-decision-tree.png")


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    pipeline_overview()
    three_inputs()
    quantification_routes()
    deseq2_flow()
    qc_decision_tree()
