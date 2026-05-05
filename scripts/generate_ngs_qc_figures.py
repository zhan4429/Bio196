from pathlib import Path
import os
import tempfile

os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", os.path.join(tempfile.gettempdir(), "fontconfig"))

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "images"
DPI = 600

COLORS = {
    "green": "#DDF2D1",
    "yellow": "#FFF0B3",
    "red": "#F6C6C3",
    "blue": "#2364AA",
    "orange": "#D77A22",
    "teal": "#1F7A7A",
    "purple": "#6A4C93",
    "text": "#1D2733",
    "muted": "#52616B",
    "line": "#2F4858",
    "input": "#EAF4F4",
    "process": "#FFF3D6",
    "output": "#E5EEF9",
    "analysis": "#EDE7F6",
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


def save(fig, name):
    OUT.mkdir(exist_ok=True)
    path = OUT / name
    fig.savefig(path, dpi=DPI, bbox_inches="tight", facecolor="white", pad_inches=0.08)
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight", facecolor="white", pad_inches=0.08)
    plt.close(fig)


def setup(figsize):
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    return fig, ax


def per_base_quality():
    fig, ax = setup((11.5, 5.6))
    x = np.arange(1, 151)
    good = 37.5 - 0.014 * x - 1.8 / (1 + np.exp(-(x - 132) / 8))
    poor = 36 - 0.035 * x - 15 / (1 + np.exp(-(x - 102) / 11))

    ax.axhspan(28, 42, color=COLORS["green"], alpha=0.86, zorder=0)
    ax.axhspan(20, 28, color=COLORS["yellow"], alpha=0.88, zorder=0)
    ax.axhspan(0, 20, color=COLORS["red"], alpha=0.82, zorder=0)

    ax.plot(x, good, color=COLORS["blue"], linewidth=3, label="Good library")
    ax.plot(x, poor, color=COLORS["orange"], linewidth=3, label="Poor 3' tail")
    ax.fill_between(x, good - 2.1, good + 1.4, color=COLORS["blue"], alpha=0.16)
    ax.fill_between(x, poor - 3.4, poor + 2.2, color=COLORS["orange"], alpha=0.18)

    ax.text(3, 35.5, "pass zone", color=COLORS["muted"], fontsize=12, weight="bold")
    ax.text(3, 24.1, "warning zone", color=COLORS["muted"], fontsize=12, weight="bold")
    ax.text(3, 12.0, "fail zone", color=COLORS["muted"], fontsize=12, weight="bold")
    ax.annotate(
        "quality often declines\nnear the 3' end",
        xy=(126, poor[125]),
        xytext=(86, 16.5),
        arrowprops=dict(arrowstyle="-|>", color=COLORS["line"], linewidth=1.3),
        color=COLORS["text"],
        fontsize=12,
    )

    ax.set_title("Per-base sequence quality")
    ax.set_xlabel("Position in read")
    ax.set_ylabel("Phred quality score")
    ax.set_xlim(1, 150)
    ax.set_ylim(0, 42)
    ax.grid(axis="y", color="white", linewidth=1.1)
    ax.legend(frameon=False, loc="lower left")
    save(fig, "09-ngs-qc-per-base-quality.png")


def adapter_content():
    fig, ax = setup((11.2, 5.4))
    x = np.arange(1, 151)
    clean = 0.2 + 0.12 * np.sin(x / 10)
    contaminated = 55 / (1 + np.exp(-(x - 72) / 9))
    light = 8 / (1 + np.exp(-(x - 118) / 10))

    ax.plot(x, clean, color=COLORS["teal"], linewidth=3, label="Clean library")
    ax.plot(x, light, color=COLORS["blue"], linewidth=3, label="Minor read-through")
    ax.plot(x, contaminated, color=COLORS["orange"], linewidth=3, label="Short inserts")
    ax.axhline(5, color=COLORS["line"], linewidth=1.2, linestyle="--")
    ax.text(4, 6.5, "investigate / trim threshold", color=COLORS["muted"], fontsize=12)
    ax.annotate(
        "adapter signal starts near\ninsert length",
        xy=(74, 28),
        xytext=(24, 38),
        arrowprops=dict(arrowstyle="-|>", color=COLORS["line"], linewidth=1.3),
        color=COLORS["text"],
        fontsize=12,
    )

    ax.set_title("Adapter content across read positions")
    ax.set_xlabel("Position in read")
    ax.set_ylabel("Reads with adapter sequence (%)")
    ax.set_xlim(1, 150)
    ax.set_ylim(0, 62)
    ax.grid(axis="y", color="#E8ECEF", linewidth=0.8)
    ax.legend(frameon=False, loc="upper left")
    save(fig, "09-ngs-qc-adapter-content.png")


def gc_content():
    fig, ax = setup((10.8, 5.4))
    gc = np.linspace(0, 100, 600)
    expected = np.exp(-0.5 * ((gc - 42) / 8.5) ** 2)
    observed_clean = np.exp(-0.5 * ((gc - 43) / 9.0) ** 2)
    observed_contam = 0.82 * np.exp(-0.5 * ((gc - 42) / 8.5) ** 2) + 0.34 * np.exp(
        -0.5 * ((gc - 63) / 5.0) ** 2
    )

    expected /= expected.max()
    observed_clean /= observed_clean.max()
    observed_contam /= observed_contam.max()

    ax.plot(gc, expected, color=COLORS["muted"], linewidth=2.4, linestyle="--", label="Expected organism")
    ax.plot(gc, observed_clean, color=COLORS["blue"], linewidth=3, label="Clean sample")
    ax.plot(gc, observed_contam, color=COLORS["orange"], linewidth=3, label="Possible contamination")
    ax.fill_between(gc, 0, observed_contam, color=COLORS["orange"], alpha=0.10)
    ax.annotate(
        "secondary peak",
        xy=(63, observed_contam[np.argmin(abs(gc - 63))]),
        xytext=(70, 0.72),
        arrowprops=dict(arrowstyle="-|>", color=COLORS["line"], linewidth=1.3),
        color=COLORS["text"],
        fontsize=12,
    )

    ax.set_title("Per-sequence GC content")
    ax.set_xlabel("GC content per read (%)")
    ax.set_ylabel("Relative number of reads")
    ax.set_xlim(15, 85)
    ax.set_ylim(0, 1.12)
    ax.grid(axis="y", color="#E8ECEF", linewidth=0.8)
    ax.legend(frameon=False, loc="upper left")
    save(fig, "09-ngs-qc-gc-content.png")


def multiqc_outlier():
    fig, ax = setup((11.8, 5.8))
    samples = ["ctrl_1", "ctrl_2", "ctrl_3", "trt_1", "trt_2", "trt_3"]
    q30 = np.array([91, 89, 61, 88, 90, 87])
    reads = np.array([24, 26, 3.2, 25, 23, 27])
    dup = np.array([28, 31, 64, 29, 34, 30])
    metrics = [q30, reads, dup]
    metric_names = ["Q30 bases", "Reads (M)", "Duplication"]
    thresholds = [80, 20, 50]
    directions = ["low", "low", "high"]

    ax.set_xlim(0, len(samples))
    ax.set_ylim(0, len(metrics))
    ax.set_axis_off()

    for row, values in enumerate(metrics):
        y = len(metrics) - row - 1
        for col, value in enumerate(values):
            if directions[row] == "low":
                bad = value < thresholds[row]
            else:
                bad = value > thresholds[row]
            color = COLORS["red"] if bad else COLORS["green"]
            ax.add_patch(Rectangle((col, y), 1, 1, facecolor=color, edgecolor="white", linewidth=2))
            label = f"{value:.0f}%" if row != 1 else f"{value:.1f}"
            ax.text(col + 0.5, y + 0.5, label, ha="center", va="center", fontsize=13, color=COLORS["text"])

    for col, sample in enumerate(samples):
        ax.text(col + 0.5, 3.15, sample, ha="center", va="bottom", fontsize=12, rotation=35, color=COLORS["text"])
    for row, name in enumerate(metric_names):
        ax.text(-0.18, len(metrics) - row - 0.5, name, ha="right", va="center", fontsize=13, color=COLORS["text"])

    ax.annotate(
        "one sample is consistently abnormal",
        xy=(2.5, 1.5),
        xytext=(3.75, 2.65),
        arrowprops=dict(arrowstyle="-|>", color=COLORS["line"], linewidth=1.4),
        color=COLORS["text"],
        fontsize=13,
    )
    ax.text(0, -0.38, "MultiQC-style summary: scan rows for samples that break the pattern.", fontsize=12, color=COLORS["muted"])
    ax.set_title("Spotting an outlier sample across QC metrics", loc="left", y=1.08, fontsize=17, weight="bold")
    save(fig, "09-ngs-qc-multiqc-outlier.png")


def box(ax, xy, text, color="process", width=2.05, height=0.74):
    x, y = xy
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.055,rounding_size=0.08",
        linewidth=1.5,
        edgecolor=COLORS["line"],
        facecolor=COLORS[color],
    )
    ax.add_patch(patch)
    ax.text(x + width / 2, y + height / 2, text, ha="center", va="center", fontsize=12.5, color=COLORS["text"])
    return x, y, width, height


def arrow(ax, start, end):
    sx, sy, sw, sh = start
    ex, ey, ew, eh = end
    ax.add_patch(
        FancyArrowPatch(
            (sx + sw, sy + sh / 2),
            (ex, ey + eh / 2),
            arrowstyle="-|>",
            mutation_scale=15,
            linewidth=1.6,
            color=COLORS["line"],
            shrinkA=7,
            shrinkB=7,
        )
    )


def qc_workflow():
    fig, ax = setup((13.6, 4.8))
    ax.set_axis_off()
    ax.set_xlim(0, 13.6)
    ax.set_ylim(0, 4.6)

    raw = box(ax, (0.35, 2.0), "Raw\nFASTQ", color="input")
    fastqc = box(ax, (2.85, 2.0), "FastQC")
    multiqc = box(ax, (5.25, 2.0), "MultiQC\nsummary")
    decision = box(ax, (7.75, 2.0), "Interpret\nmodules", color="analysis")
    trim = box(ax, (10.25, 3.0), "Trim/filter\nfastp", color="process")
    align = box(ax, (10.25, 1.0), "Proceed to\nalignment", color="output")

    for a, b in [(raw, fastqc), (fastqc, multiqc), (multiqc, decision)]:
        arrow(ax, a, b)
    arrow(ax, decision, trim)
    arrow(ax, decision, align)

    ax.text(10.25, 4.02, "adapter, low quality,\ntoo many Ns", fontsize=11.5, color=COLORS["muted"])
    ax.text(10.25, 0.44, "QC acceptable\nor expected warnings", fontsize=11.5, color=COLORS["muted"])
    ax.set_title("Raw-read QC decision workflow", loc="left", fontsize=17, weight="bold")
    save(fig, "09-ngs-qc-workflow.png")


if __name__ == "__main__":
    per_base_quality()
    adapter_content()
    gc_content()
    multiqc_outlier()
    qc_workflow()
