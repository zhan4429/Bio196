"""Generate the genomic coordinate systems figure for Chapter 2."""

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "images"

# Base colours for nucleotides — vivid enough to read, muted enough for a textbook
BASE_COLORS = {
    "A": "#E74C3C",  # red
    "T": "#3498DB",  # blue
    "G": "#2ECC71",  # green
    "C": "#F39C12",  # amber
}
HIGHLIGHT_COLORS = {
    "A": "#C0392B",
    "T": "#2980B9",
    "G": "#27AE60",
    "C": "#E67E22",
}
DIM_ALPHA = 0.30  # opacity for bases outside the highlighted region
HIGHLIGHT_ALPHA = 1.0

TEXT = "#1F2933"
MUTED = "#6B7280"
ACCENT = "#2F5D62"
BED_COLOR = "#C0392B"   # red for BED annotation
GFF_COLOR = "#27AE60"   # green for GFF annotation


def draw_base_row(ax, y, sequence, numbering_start, highlight_start, highlight_end,
                  box_w=0.72, box_h=0.72, gap=0.04, x_offset=0.8):
    """Draw a row of nucleotide boxes with position numbers underneath.

    highlight_start/end are indices into the sequence (0-based python index)
    that should be visually highlighted.
    """
    for i, base in enumerate(sequence):
        x = x_offset + i * (box_w + gap)
        in_region = highlight_start <= i < highlight_end

        # Box
        alpha = HIGHLIGHT_ALPHA if in_region else DIM_ALPHA
        fc = HIGHLIGHT_COLORS[base] if in_region else BASE_COLORS[base]
        patch = FancyBboxPatch(
            (x, y), box_w, box_h,
            boxstyle="round,pad=0.02,rounding_size=0.08",
            linewidth=1.2 if in_region else 0.6,
            edgecolor="#FFFFFF" if in_region else "#D1D5DB",
            facecolor=fc,
            alpha=alpha,
        )
        ax.add_patch(patch)

        # Base letter
        ax.text(
            x + box_w / 2, y + box_h / 2, base,
            ha="center", va="center",
            fontsize=13, weight="bold",
            color="white",
            alpha=HIGHLIGHT_ALPHA if in_region else 0.55,
        )

        # Position number below
        pos_label = str(numbering_start + i)
        ax.text(
            x + box_w / 2, y - 0.25, pos_label,
            ha="center", va="center",
            fontsize=8.5,
            color=TEXT if in_region else MUTED,
            weight="bold" if in_region else "normal",
        )


def draw_bracket(ax, y, x_start, x_end, color, label, box_w=0.72, gap=0.04,
                 x_offset=0.8, bracket_y_offset=-0.50):
    """Draw a bracket below the position numbers with a label."""
    left = x_offset + x_start * (box_w + gap) + 0.02
    right = x_offset + (x_end - 1) * (box_w + gap) + box_w - 0.02
    by = y + bracket_y_offset
    tick = 0.08

    # Bracket lines
    lw = 2.0
    ax.plot([left, left], [by + tick, by], color=color, lw=lw, solid_capstyle="round")
    ax.plot([left, right], [by, by], color=color, lw=lw, solid_capstyle="round")
    ax.plot([right, right], [by, by + tick], color=color, lw=lw, solid_capstyle="round")

    # Label centered below bracket
    ax.text(
        (left + right) / 2, by - 0.22, label,
        ha="center", va="center",
        fontsize=9, weight="bold", color=color,
    )


def coordinate_systems():
    sequence = list("ATGCTAGCTAGCTAGCTA")  # 18 bases
    n = len(sequence)

    # The highlighted region: 0-based indices 4..9 (python slice [4:10])
    # 0-based half-open: [4, 10)  → positions 4,5,6,7,8,9
    # 1-based closed:    [5, 10]  → positions 5,6,7,8,9,10
    hl_start = 4   # python index of first highlighted base
    hl_end = 10    # python index past last highlighted base

    fig_w = 16
    fig_h = 7.8
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.set_axis_off()
    ax.set_xlim(0, fig_w)
    ax.set_ylim(0, fig_h)

    # --- Title ---
    ax.text(
        fig_w / 2, 7.35,
        "Genomic Coordinate Systems: Same Region, Different Conventions",
        ha="center", va="center",
        fontsize=18, weight="bold", color=TEXT,
    )

    # =====================================================================
    # Top section: 0-based, half-open
    # =====================================================================
    top_label_y = 6.55
    top_row_y = 5.55

    ax.text(
        fig_w / 2, top_label_y,
        "0-based, half-open  (BED, BAM)   →   region [4, 10)  covers 6 bases: TAGCTA",
        ha="center", va="center",
        fontsize=11.5, color=MUTED,
    )

    draw_base_row(ax, top_row_y, sequence, numbering_start=0,
                  highlight_start=hl_start, highlight_end=hl_end)

    draw_bracket(ax, top_row_y, hl_start, hl_end, BED_COLOR,
                 "BED:  chrom  4  10   (start=4, end=10)")

    # =====================================================================
    # Bottom section: 1-based, closed
    # =====================================================================
    bot_label_y = 3.65
    bot_row_y = 2.65

    ax.text(
        fig_w / 2, bot_label_y,
        "1-based, closed  (GFF, GTF, VCF, SAM)   →   region [5, 10]  covers 6 bases: TAGCTA",
        ha="center", va="center",
        fontsize=11.5, color=MUTED,
    )

    draw_base_row(ax, bot_row_y, sequence, numbering_start=1,
                  highlight_start=hl_start, highlight_end=hl_end)

    draw_bracket(ax, bot_row_y, hl_start, hl_end, GFF_COLOR,
                 "GFF/VCF:  seqid  5  10   (start=5, end=10)")

    # =====================================================================
    # Connecting note
    # =====================================================================
    ax.text(
        fig_w / 2, 1.45,
        "Both notations describe the same 6 bases.  "
        "Converting:  1-based start = 0-based start + 1  ;  end is the same integer.",
        ha="center", va="center",
        fontsize=10.5, color=ACCENT, weight="bold",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#EAF4F4", edgecolor=ACCENT,
                  linewidth=1.2),
    )

    # --- Save ---
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / "02-coordinate-systems.png", dpi=300, bbox_inches="tight",
                facecolor="white")
    plt.close(fig)
    print(f"Saved {OUT / '02-coordinate-systems.png'}")


if __name__ == "__main__":
    coordinate_systems()
