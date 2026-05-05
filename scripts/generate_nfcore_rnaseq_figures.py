from pathlib import Path
import os
import tempfile

os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", os.path.join(tempfile.gettempdir(), "fontconfig"))
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "images"
DPI = 450

COLORS = {
    "input": "#EAF4F4",
    "process": "#FFF3D6",
    "route": "#EDE7F6",
    "output": "#E5EEF9",
    "qc": "#F8DDD8",
    "accent": "#1F5F6B",
    "text": "#1D2733",
}


def setup(figsize=(12, 5)):
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.set_axis_off()
    return fig, ax


def box(ax, x, y, text, w=2.0, h=0.72, kind="process", fs=11):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.04,rounding_size=0.08",
        facecolor=COLORS[kind],
        edgecolor=COLORS["accent"],
        linewidth=1.45,
    )
    ax.add_patch(patch)
    ax.text(
        x + w / 2,
        y + h / 2,
        text,
        ha="center",
        va="center",
        fontsize=fs,
        color=COLORS["text"],
        linespacing=1.15,
    )
    return (x, y, w, h)


def left(b):
    return (b[0], b[1] + b[3] / 2)


def right(b):
    return (b[0] + b[2], b[1] + b[3] / 2)


def top(b):
    return (b[0] + b[2] / 2, b[1] + b[3])


def bottom(b):
    return (b[0] + b[2] / 2, b[1])


def right_at(b, frac):
    return (b[0] + b[2], b[1] + b[3] * frac)


def left_at(b, frac):
    return (b[0], b[1] + b[3] * frac)


def arrow_xy(ax, start, end, curved=0.0):
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=15,
            color=COLORS["accent"],
            linewidth=1.6,
            shrinkA=6,
            shrinkB=6,
            connectionstyle=f"arc3,rad={curved}",
        )
    )


def arrow(ax, a, b, curved=0.0):
    ax.add_patch(
        FancyArrowPatch(
            (a[0] + a[2], a[1] + a[3] / 2),
            (b[0], b[1] + b[3] / 2),
            arrowstyle="-|>",
            mutation_scale=15,
            color=COLORS["accent"],
            linewidth=1.6,
            shrinkA=5,
            shrinkB=5,
            connectionstyle=f"arc3,rad={curved}",
        )
    )


def arrow_down(ax, a, b):
    ax.add_patch(
        FancyArrowPatch(
            (a[0] + a[2] / 2, a[1]),
            (b[0] + b[2] / 2, b[1] + b[3]),
            arrowstyle="-|>",
            mutation_scale=15,
            color=COLORS["accent"],
            linewidth=1.6,
            shrinkA=5,
            shrinkB=5,
        )
    )


def save(fig, name):
    path = OUT / name
    fig.savefig(path, dpi=DPI, bbox_inches="tight", facecolor="white", pad_inches=0.08)
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight", facecolor="white", pad_inches=0.08)
    plt.close(fig)


def overview():
    fig, ax = setup((15, 5.2))
    ax.set_xlim(0, 15)
    ax.set_ylim(0, 5.2)

    samples = box(ax, 0.35, 3.75, "samplesheet.csv", kind="input", w=2.2)
    reads = box(ax, 0.35, 2.45, "FASTQ files", kind="input", w=2.2)
    refs = box(ax, 0.35, 1.15, "FASTA +\nGTF/GFF", kind="input", w=2.2)

    check = box(ax, 3.25, 3.05, "Check inputs\nmerge lanes", kind="qc", w=2.15)
    qc = box(ax, 5.85, 3.05, "FastQC\ntrim/filter\nrRNA optional", kind="qc", w=2.25)
    route = box(ax, 8.55, 2.15, "Choose analysis route\nalignment or\npseudoalignment", kind="route", w=2.55, h=0.95)
    matrices = box(ax, 11.85, 2.95, "Count and\nabundance\nmatrices", kind="output", w=2.35, h=0.9)
    multiqc = box(ax, 11.85, 1.15, "MultiQC\nHTML report", kind="output", w=2.35, h=0.9)

    arrow_xy(ax, right(samples), left(check), curved=-0.08)
    arrow_xy(ax, right(reads), left(check))
    arrow_xy(ax, right(refs), left(route), curved=0.05)
    arrow(ax, check, qc)
    arrow_xy(ax, right(qc), left(route), curved=-0.08)
    arrow_xy(ax, right(route), left(matrices), curved=-0.08)
    arrow_xy(ax, right(route), left(multiqc), curved=0.08)

    ax.text(12.0, 4.15, "Use for DESeq2,\nedgeR, limma", fontsize=11, color=COLORS["text"])
    ax.text(12.0, 0.45, "Start QC review here", fontsize=11, color=COLORS["text"])

    save(fig, "15-nfcore-rnaseq-overview.png")


def route_choice():
    fig, ax = setup((15.2, 6.1))
    ax.set_xlim(0, 15.2)
    ax.set_ylim(0, 6.1)

    start = box(ax, 0.35, 2.65, "Choose an\nanalysis route", kind="input", w=2.15)

    routes = [
        box(ax, 3.35, 4.85, "--aligner\nstar_salmon", kind="route", w=2.55),
        box(ax, 3.35, 3.45, "--aligner\nstar_rsem", kind="route", w=2.55),
        box(ax, 3.35, 2.05, "--aligner\nhisat2", kind="route", w=2.55),
        box(ax, 3.35, 0.55, "--pseudo_aligner\nsalmon/kallisto\n--skip_alignment", kind="route", w=2.55, h=0.88),
    ]
    outputs = [
        box(ax, 6.6, 4.85, "STAR BAMs\nSalmon counts", kind="output", w=2.6),
        box(ax, 6.6, 3.45, "STAR BAMs\nRSEM counts", kind="output", w=2.6),
        box(ax, 6.6, 2.05, "HISAT2 BAMs\nQC only\nno counts", kind="output", w=2.6, h=0.88),
        box(ax, 6.6, 0.55, "Fast transcript\nquantification\nfewer alignment QC plots", kind="output", w=2.6, h=0.88),
    ]
    guidance = [
        box(ax, 10.2, 4.85, "Best first choice\nfor most bulk\nRNA-seq projects", kind="qc", w=3.0, h=0.88),
        box(ax, 10.2, 2.05, "Use when STAR\nmemory is too high", kind="qc", w=3.0, h=0.78),
    ]

    for source_frac, b in zip([0.92, 0.65, 0.35, 0.08], routes):
        arrow_xy(ax, right_at(start, source_frac), left(b), curved=0.0)
    for a, b in zip(routes, outputs):
        arrow(ax, a, b)
    arrow(ax, outputs[0], guidance[0])
    arrow(ax, outputs[2], guidance[1])

    save(fig, "15-nfcore-rnaseq-routes.png")


def samplesheet():
    fig, ax = setup((15, 5.9))
    ax.set_xlim(0, 15)
    ax.set_ylim(0, 5.9)

    headers = [
        box(ax, 0.5, 4.55, "sample", kind="input", w=2.1),
        box(ax, 3.1, 4.55, "fastq_1", kind="input", w=2.1),
        box(ax, 5.7, 4.55, "fastq_2", kind="input", w=2.1),
        box(ax, 8.3, 4.55, "strandedness", kind="input", w=2.25),
        box(ax, 11.2, 4.55, "optional\nmetadata", kind="input", w=2.25),
    ]
    details = [
        box(ax, 0.35, 2.9, "Same name for\nextra lanes of\nthe same sample", kind="process", w=2.4, h=0.95),
        box(ax, 2.95, 2.9, "Read 1 file\nmust be gzipped\nFASTQ", kind="process", w=2.4, h=0.95),
        box(ax, 5.55, 2.9, "Read 2 file\nblank for\nsingle-end", kind="process", w=2.4, h=0.95),
        box(ax, 8.25, 2.9, "unstranded,\nforward,\nreverse, or auto", kind="process", w=2.45, h=0.95),
        box(ax, 11.05, 2.9, "seq_platform,\nseq_center,\nBAM paths", kind="process", w=2.45, h=0.95),
    ]
    for a, b in zip(headers, details):
        arrow_down(ax, a, b)

    checked = box(
        ax,
        1.25,
        0.55,
        "A clean samplesheet prevents most painful pipeline failures",
        kind="qc",
        w=11.6,
        h=0.78,
    )
    target_y = checked[1] + checked[3]
    for i, d in enumerate(details):
        target_x = checked[0] + checked[2] * (i + 1) / (len(details) + 1)
        arrow_xy(ax, bottom(d), (target_x, target_y), curved=0.0)

    save(fig, "15-nfcore-rnaseq-samplesheet.png")


def outputs():
    fig, ax = setup((15.4, 6.2))
    ax.set_xlim(0, 15)
    ax.set_ylim(0, 6.2)

    root = box(ax, 0.4, 2.65, "results/", kind="input", w=1.85)
    folders = [
        box(ax, 3.0, 4.85, "fastqc/\ntrimgalore/\nfastp/", kind="qc", w=2.35, h=0.85),
        box(ax, 3.0, 3.55, "star_salmon/\nstar_rsem/\nhisat2/", kind="output", w=2.35, h=0.85),
        box(ax, 3.0, 2.25, "rseqc/\nqualimap/\ndupradar/", kind="qc", w=2.35, h=0.85),
        box(ax, 3.0, 0.95, "multiqc/", kind="output", w=2.35, h=0.85),
    ]

    count = box(ax, 6.4, 3.8, "merged gene\ncounts / TPM", kind="output", w=2.55, h=0.85)
    se = box(ax, 10.0, 3.8, "Summarized\nExperiment RDS", kind="output", w=2.55, h=0.85)
    logs = box(ax, 6.4, 2.55, "alignment logs\nquantification logs", kind="process", w=2.55, h=0.85)
    html = box(ax, 6.4, 1.15, "multiqc_report.html\nstart reading here", kind="qc", w=2.9, h=0.9)
    pipe = box(ax, 10.35, 0.25, "pipeline_info/\nversions, trace,\nrun reports", kind="output", w=2.65, h=1.0)

    for b in folders:
        arrow_xy(ax, right(root), left(b), curved=0.0)
    arrow(ax, folders[1], count)
    arrow(ax, count, se)
    arrow(ax, folders[1], logs)
    arrow(ax, folders[3], html)
    arrow_xy(ax, right(folders[3]), left(pipe), curved=0.18)

    ax.text(10.0, 2.65, "Use matrices for\nDESeq2, edgeR,\nand limma.", fontsize=11, color=COLORS["text"])
    ax.text(
        10.15,
        1.75,
        "Use MultiQC to decide\nwhether the run is\ntrustworthy.",
        fontsize=11,
        color=COLORS["text"],
    )

    save(fig, "15-nfcore-rnaseq-outputs.png")


def alignment_vs_pseudo():
    fig, ax = setup((15, 6.1))
    ax.set_xlim(0, 15)
    ax.set_ylim(0, 6.1)

    question = box(ax, 0.35, 2.65, "What do you\nneed from the\nanalysis?", kind="input", w=2.35, h=0.95)

    q1 = box(ax, 3.4, 4.85, "Routine gene\nexpression?", kind="process", w=2.45, h=0.8)
    q4 = box(ax, 3.4, 3.55, "Many samples\nor limited\ncompute?", kind="process", w=2.45, h=0.85)
    q2 = box(ax, 3.4, 2.15, "Need IGV,\nsplice junctions,\nintrons, fusions?", kind="process", w=2.45, h=0.95)
    q3 = box(ax, 3.4, 0.8, "Poor or unusual\nannotation?", kind="process", w=2.45, h=0.8)

    pseudo = box(ax, 6.8, 4.15, "Pseudoalignment\nSalmon-only", kind="route", w=2.75, h=0.9)
    align = box(ax, 6.8, 1.7, "Genome alignment\nSTAR or HISAT2", kind="route", w=2.75, h=0.9)
    hybrid = box(ax, 10.5, 2.45, "Best beginner\ndefault:\nstar_salmon", kind="qc", w=2.8, h=0.95)
    pseudo_note = box(ax, 10.5, 4.35, "Fast, small outputs\nno BAM by default", kind="output", w=2.8, h=0.85)
    align_note = box(ax, 10.5, 0.75, "Richer QC\nBAM files\nmore storage", kind="output", w=2.8, h=0.9)

    for source_frac, b in zip([0.92, 0.65, 0.35, 0.08], [q1, q2, q3, q4]):
        arrow_xy(ax, right_at(question, source_frac), left(b), curved=0.0)
    arrow_xy(ax, right(q1), left_at(pseudo, 0.72), curved=-0.05)
    arrow_xy(ax, right(q4), left_at(pseudo, 0.28), curved=0.05)
    arrow_xy(ax, right(q2), left_at(align, 0.72), curved=0.04)
    arrow_xy(ax, right(q3), left_at(align, 0.28), curved=-0.04)
    arrow(ax, pseudo, pseudo_note)
    arrow_xy(ax, right(align), left(hybrid), curved=-0.05)
    arrow_xy(ax, right(align), left(align_note), curved=0.08)

    save(fig, "15-nfcore-rnaseq-alignment-vs-pseudo.png")


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    overview()
    route_choice()
    samplesheet()
    outputs()
    alignment_vs_pseudo()
