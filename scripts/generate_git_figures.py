from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Circle


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "images"

COLORS = {
    "work": "#E8F3F1",
    "stage": "#F6E7CB",
    "repo": "#DDE8F7",
    "branch": "#E8E1F5",
    "danger": "#F7D9D5",
    "maybe": "#F4E3B2",
    "accent": "#2F5D62",
    "text": "#1F2933",
    "muted": "#6B7280",
    "line": "#334E52",
    "white": "#FFFFFF",
}


def setup(figsize=(13, 5.5)):
    fig, ax = plt.subplots(figsize=figsize)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.set_axis_off()
    return fig, ax


def rounded(ax, x, y, w, h, text, fill="work", fs=13, weight="normal"):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.05,rounding_size=0.12",
        linewidth=1.6,
        edgecolor=COLORS["accent"],
        facecolor=COLORS[fill],
    )
    ax.add_patch(patch)
    ax.text(
        x + w / 2,
        y + h / 2,
        text,
        ha="center",
        va="center",
        fontsize=fs,
        weight=weight,
        color=COLORS["text"],
        linespacing=1.15,
    )
    return (x, y, w, h)


def arrow(ax, a, b, label=None, curve=0.0, yoffset=0.0):
    sx, sy, sw, sh = a
    ex, ey, ew, eh = b
    start = (sx + sw, sy + sh / 2 + yoffset)
    end = (ex, ey + eh / 2 + yoffset)
    patch = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=16,
        linewidth=1.8,
        color=COLORS["line"],
        shrinkA=7,
        shrinkB=7,
        connectionstyle=f"arc3,rad={curve}",
    )
    ax.add_patch(patch)
    if label:
        mx = (start[0] + end[0]) / 2
        my = (start[1] + end[1]) / 2
        ax.text(
            mx,
            my + 0.35,
            label,
            ha="center",
            va="center",
            fontsize=11,
            color=COLORS["accent"],
            weight="bold",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none"),
        )


def save(fig, name):
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / name, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def git_three_areas():
    fig, ax = setup((13.5, 5.4))
    ax.set_xlim(0, 13.5)
    ax.set_ylim(0, 5.4)

    ax.text(0.4, 5.05, "Git's three areas", fontsize=20, weight="bold", color=COLORS["text"])
    ax.text(
        0.4,
        4.66,
        "Edit files, choose what belongs in the next snapshot, then save it permanently.",
        fontsize=12,
        color=COLORS["muted"],
    )

    work = rounded(ax, 0.5, 2.45, 3.1, 1.35, "Working directory\nfiles you are editing", "work", 13, "bold")
    stage = rounded(ax, 5.05, 2.45, 3.1, 1.35, "Staging area\nselected changes", "stage", 13, "bold")
    repo = rounded(ax, 9.6, 2.45, 3.1, 1.35, "Repository\ncommit history", "repo", 13, "bold")

    arrow(ax, work, stage, "git add")
    arrow(ax, stage, repo, "git commit")

    unstaged = rounded(ax, 0.85, 0.75, 2.4, 0.75, "modified file", "white", 11)
    staged = rounded(ax, 5.4, 0.75, 2.4, 0.75, "ready for commit", "white", 11)
    commit = rounded(ax, 9.95, 0.75, 2.4, 0.75, "snapshot", "white", 11)
    arrow(ax, unstaged, staged, "stage", yoffset=0)
    arrow(ax, staged, commit, "save", yoffset=0)

    ax.text(0.55, 1.85, "You can edit freely here.", fontsize=11, color=COLORS["muted"])
    ax.text(4.85, 1.85, "This is your packing list.", fontsize=11, color=COLORS["muted"])
    ax.text(9.45, 1.85, "Commits are recoverable checkpoints.", fontsize=11, color=COLORS["muted"])

    save(fig, "10-git-three-areas.png")


def git_branching():
    fig, ax = setup((13.5, 5.2))
    ax.set_xlim(0, 13.5)
    ax.set_ylim(0, 5.2)

    ax.text(0.4, 4.85, "Branching lets experiments stay separate", fontsize=20, weight="bold", color=COLORS["text"])
    ax.text(
        0.4,
        4.46,
        "Try a new idea on a branch. Merge it only when it is ready.",
        fontsize=12,
        color=COLORS["muted"],
    )

    main_y = 2.25
    feat_y = 3.45
    xs = [1.0, 2.7, 4.4, 6.1, 7.8, 9.5, 11.2]
    labels = ["A", "B", "C", "F", "G", "H", "I"]
    for x, lab in zip(xs, labels):
        circ = Circle((x, main_y), 0.24, facecolor=COLORS["repo"], edgecolor=COLORS["accent"], lw=1.5)
        ax.add_patch(circ)
        ax.text(x, main_y, lab, ha="center", va="center", fontsize=12, weight="bold", color=COLORS["text"])
    for x1, x2 in zip(xs[:-1], xs[1:]):
        ax.plot([x1 + 0.24, x2 - 0.24], [main_y, main_y], color=COLORS["line"], lw=2)

    feat_xs = [4.4, 6.1, 7.8]
    feat_labels = ["D", "E", "E'"]
    for x, lab in zip(feat_xs, feat_labels):
        circ = Circle((x, feat_y), 0.24, facecolor=COLORS["branch"], edgecolor=COLORS["accent"], lw=1.5)
        ax.add_patch(circ)
        ax.text(x, feat_y, lab, ha="center", va="center", fontsize=12, weight="bold", color=COLORS["text"])
    ax.plot([2.7 + 0.22, 4.4 - 0.24], [main_y + 0.05, feat_y - 0.05], color=COLORS["line"], lw=2)
    ax.plot([4.4 + 0.24, 6.1 - 0.24], [feat_y, feat_y], color=COLORS["line"], lw=2)
    ax.plot([6.1 + 0.24, 7.8 - 0.24], [feat_y, feat_y], color=COLORS["line"], lw=2)
    ax.plot([7.8 + 0.24, 9.5 - 0.24], [feat_y - 0.05, main_y + 0.05], color=COLORS["line"], lw=2)

    ax.text(0.7, 1.55, "main branch", fontsize=12, weight="bold", color=COLORS["accent"])
    ax.text(4.5, 4.02, "feature branch", fontsize=12, weight="bold", color=COLORS["accent"])
    ax.text(3.05, 3.15, "git switch -c", fontsize=10, color=COLORS["muted"], rotation=31)
    ax.text(8.3, 3.1, "git merge", fontsize=10, color=COLORS["muted"], rotation=-31)
    ax.text(9.6, 1.58, "main now contains the finished experiment", fontsize=11, color=COLORS["muted"])

    save(fig, "10-git-branching.png")


def project_structure():
    fig, ax = setup((11.5, 8.2))
    ax.set_xlim(0, 11.5)
    ax.set_ylim(0, 8.2)

    ax.text(0.4, 7.78, "What belongs in Git?", fontsize=20, weight="bold", color=COLORS["text"])
    ax.text(0.4, 7.38, "Track small handwritten files. Keep large data and reproducible outputs elsewhere.", fontsize=12, color=COLORS["muted"])

    col1 = rounded(ax, 0.55, 6.1, 3.1, 0.72, "Track in Git", "work", 13, "bold")
    col2 = rounded(ax, 4.2, 6.1, 3.1, 0.72, "Do not track", "danger", 13, "bold")
    col3 = rounded(ax, 7.85, 6.1, 3.1, 0.72, "Maybe track", "maybe", 13, "bold")

    track = [
        "scripts/*.R, *.py, *.sh",
        "Snakefile, nextflow.config",
        "environment.yml",
        "samplesheet.csv",
        "README.md",
        ".gitignore",
    ]
    no = [
        "FASTQ files",
        "BAM, SAM, CRAM",
        "reference genomes",
        "large result folders",
        "temporary work/",
        "secrets and keys",
    ]
    maybe = [
        "final small tables",
        "publication figures",
        "small BED/GTF subsets",
        "reports/*.qmd",
        "metadata dictionaries",
        "tiny example data",
    ]

    def items(x, y, vals):
        for i, val in enumerate(vals):
            rounded(ax, x, y - i * 0.72, 3.1, 0.48, val, "white", 10.5)

    items(0.55, 5.35, track)
    items(4.2, 5.35, no)
    items(7.85, 5.35, maybe)

    rounded(
        ax,
        1.25,
        0.45,
        9.0,
        0.75,
        "Rule of thumb:\nlarge, generated, private, or downloadable files stay out of Git.",
        "repo",
        11,
        "bold",
    )

    save(fig, "10-project-structure.png")


def github_workflow():
    fig, ax = setup((13.5, 7.0))
    ax.set_xlim(0, 13.5)
    ax.set_ylim(0, 7.0)

    ax.text(0.4, 6.62, "Fork and pull request workflow", fontsize=20, weight="bold", color=COLORS["text"])
    ax.text(0.4, 6.23, "A safe way to contribute changes without editing the original project directly.", fontsize=12, color=COLORS["muted"])

    boxes = [
        rounded(ax, 0.45, 4.45, 1.85, 0.82, "Original\nrepository", "repo", 11, "bold"),
        rounded(ax, 2.95, 4.45, 1.85, 0.82, "Fork to\nyour account", "work", 11, "bold"),
        rounded(ax, 5.45, 4.45, 1.85, 0.82, "Clone\nlocally", "work", 11, "bold"),
        rounded(ax, 7.95, 4.45, 1.85, 0.82, "Create\nbranch", "branch", 11, "bold"),
        rounded(ax, 10.45, 4.45, 1.85, 0.82, "Commit\nchanges", "stage", 11, "bold"),
    ]
    for a, b in zip(boxes[:-1], boxes[1:]):
        arrow(ax, a, b)

    push = rounded(ax, 10.45, 3.05, 1.85, 0.82, "Push to\nyour fork", "stage", 11, "bold")
    pr = rounded(ax, 10.45, 1.55, 1.85, 0.82, "Open pull\nrequest", "branch", 11, "bold")
    review = rounded(ax, 6.85, 1.55, 2.15, 0.82, "Review and\nchanges", "maybe", 11, "bold")
    merge = rounded(ax, 2.95, 1.55, 2.15, 0.82, "Merge into\noriginal", "repo", 11, "bold")

    # Route the lower workflow around boxes so connector lines do not cross labels.
    ax.add_patch(
        FancyArrowPatch(
            (11.38, 4.45),
            (11.38, 3.90),
            arrowstyle="-|>",
            mutation_scale=16,
            linewidth=1.8,
            color=COLORS["line"],
        )
    )
    ax.add_patch(
        FancyArrowPatch(
            (11.38, 3.05),
            (11.38, 2.40),
            arrowstyle="-|>",
            mutation_scale=16,
            linewidth=1.8,
            color=COLORS["line"],
        )
    )
    for a, b in [(pr, review), (review, merge)]:
        ax.add_patch(
            FancyArrowPatch(
                (a[0], a[1] + a[3] / 2),
                (b[0] + b[2], b[1] + b[3] / 2),
                arrowstyle="-|>",
                mutation_scale=16,
                linewidth=1.8,
                color=COLORS["line"],
                shrinkA=7,
                shrinkB=7,
            )
        )

    ax.text(0.55, 3.55, "You work in your copy.", fontsize=11, color=COLORS["muted"])
    ax.text(6.15, 2.75, "The maintainer reviews before merging.", fontsize=11, color=COLORS["muted"])
    ax.text(3.0, 5.65, "GitHub", fontsize=12, weight="bold", color=COLORS["accent"])
    ax.text(6.0, 5.65, "Your computer", fontsize=12, weight="bold", color=COLORS["accent"])

    save(fig, "10-github-workflow.png")


def undo_safety():
    fig, ax = setup((13.5, 7.8))
    ax.set_xlim(0, 13.5)
    ax.set_ylim(0, 7.8)

    ax.text(0.45, 7.35, "Undoing changes safely", fontsize=20, weight="bold", color=COLORS["text"])
    ax.text(
        0.45,
        6.96,
        "First ask where the change is. Then choose the least destructive command.",
        fontsize=12,
        color=COLORS["muted"],
    )

    header_y = 5.82
    rounded(ax, 0.55, header_y, 3.15, 0.46, "Situation", "repo", 10.5, "bold")
    rounded(ax, 4.35, header_y, 3.35, 0.46, "Command", "repo", 10.5, "bold")
    rounded(ax, 8.35, header_y, 4.55, 0.46, "What it does", "repo", 10.5, "bold")

    rows = [
        (
            4.80,
            "Already pushed?",
            "git revert <hash>",
            "Safest for public history.\nCreates a new commit that undoes the old one.",
            "work",
        ),
        (
            3.55,
            "Staged by mistake?",
            "git restore --staged <file>",
            "Removes the file from the staging area.\nYour edits stay in the working directory.",
            "stage",
        ),
        (
            2.30,
            "Unstaged edits unwanted?",
            "git restore <file>",
            "Discards local edits in that file.\nUse only when you do not need those edits.",
            "maybe",
        ),
        (
            1.05,
            "Need to move HEAD?",
            "git reset --soft HEAD~1",
            "Undoes the last commit but keeps changes staged.\nUsually safer than hard reset.",
            "repo",
        ),
    ]

    for y, problem, command, note, fill in rows:
        rounded(ax, 0.55, y, 3.15, 0.72, problem, fill, 10.0, "bold")
        rounded(ax, 4.35, y, 3.35, 0.72, command, "white", 9.6, "bold")
        rounded(ax, 8.35, y - 0.03, 4.55, 0.78, note, "white", 8.8)
        ax.add_patch(
            FancyArrowPatch(
                (3.73, y + 0.36),
                (4.32, y + 0.36),
                arrowstyle="-|>",
                mutation_scale=13,
                linewidth=1.4,
                color=COLORS["line"],
                shrinkA=3,
                shrinkB=3,
            )
        )
        ax.add_patch(
            FancyArrowPatch(
                (7.73, y + 0.36),
                (8.32, y + 0.36),
                arrowstyle="-|>",
                mutation_scale=13,
                linewidth=1.4,
                color=COLORS["line"],
                shrinkA=3,
                shrinkB=3,
            )
        )

    rounded(
        ax,
        0.55,
        0.18,
        12.35,
        0.58,
        "Last resort: git reset --hard discards local changes. Use it only when you are completely sure.",
        "danger",
        9.0,
        "bold",
    )

    save(fig, "10-git-undo-safety.png")


if __name__ == "__main__":
    git_three_areas()
    git_branching()
    project_structure()
    github_workflow()
    undo_safety()
