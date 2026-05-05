from pathlib import Path
import os
import tempfile

os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", os.path.join(tempfile.gettempdir(), "fontconfig"))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "images"
DPI = 450

COLORS = {
    "input": "#EAF4F4",
    "process": "#FFF3D6",
    "model": "#EDE7F6",
    "output": "#E5EEF9",
    "warning": "#F8DDD8",
    "good": "#DDEFD8",
    "accent": "#1F5F6B",
    "orange": "#C96F3A",
    "green": "#4B8B55",
    "rose": "#B45A63",
    "gold": "#C49A32",
    "text": "#1D2733",
    "muted": "#6B7280",
    "grid": "#D9DEE6",
}


def setup(figsize=(12, 5)):
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.set_axis_off()
    return fig, ax


def save(fig, name):
    path = OUT / name
    fig.savefig(path, dpi=DPI, bbox_inches="tight", facecolor="white", pad_inches=0.08)
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight", facecolor="white", pad_inches=0.08)
    plt.close(fig)


def box(ax, x, y, text, w=2.2, h=0.8, kind="process", fs=11):
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


def right(b):
    return (b[0] + b[2], b[1] + b[3] / 2)


def left(b):
    return (b[0], b[1] + b[3] / 2)


def top(b):
    return (b[0] + b[2] / 2, b[1] + b[3])


def bottom(b):
    return (b[0] + b[2] / 2, b[1])


def left_at(b, frac):
    return (b[0], b[1] + b[3] * frac)


def right_at(b, frac):
    return (b[0] + b[2], b[1] + b[3] * frac)


def arrow_xy(ax, start, end, curved=0.0, color=None):
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=15,
            color=color or COLORS["accent"],
            linewidth=1.6,
            shrinkA=6,
            shrinkB=6,
            connectionstyle=f"arc3,rad={curved}",
        )
    )


def arrow(ax, a, b, curved=0.0):
    arrow_xy(ax, right(a), left(b), curved=curved)


def ml_workflow():
    fig, ax = setup((15.4, 5.6))
    ax.set_xlim(0, 15.4)
    ax.set_ylim(0, 5.6)

    raw = box(ax, 0.35, 3.65, "Biological data\ncounts, variants,\nimages", w=2.2, h=1.0, kind="input")
    clean = box(ax, 3.15, 3.65, "Preprocess\nfilter, normalize,\ntransform", w=2.3, h=1.0)
    split = box(ax, 6.15, 3.65, "Split samples\ntrain and test", w=2.15, h=1.0, kind="process")
    train = box(ax, 9.0, 4.25, "Training set\nfit scaler,\nselect features,\ntrain model", w=2.55, h=1.0, kind="model")
    test = box(ax, 9.0, 2.55, "Test set\ntransform only\nnever fit", w=2.55, h=1.0, kind="warning")
    model = box(ax, 12.15, 4.25, "Trained model", w=2.05, h=0.8, kind="model")
    metrics = box(ax, 12.15, 2.55, "Evaluation\nAUC, F1, RMSE", w=2.05, h=0.8, kind="output")

    for a, b in [(raw, clean), (clean, split)]:
        arrow(ax, a, b)
    arrow_xy(ax, right_at(split, 0.75), left(train), curved=-0.05)
    arrow_xy(ax, right_at(split, 0.25), left(test), curved=0.05)
    arrow(ax, train, model)
    arrow_xy(ax, bottom(model), top(metrics), curved=0.0)
    arrow(ax, test, metrics)

    ax.text(
        0.45,
        0.8,
        "Rule of thumb: anything that learns from data must learn from training samples only.",
        fontsize=12,
        color=COLORS["text"],
    )
    save(fig, "16-ml-workflow.png")


def leakage():
    fig, ax = setup((15.8, 6.6))
    ax.set_xlim(0, 15.8)
    ax.set_ylim(0, 6.6)

    ax.text(2.15, 6.05, "Leaky workflow", ha="center", fontsize=14, weight="bold", color=COLORS["rose"])
    ax.text(12.0, 6.05, "Correct workflow", ha="center", fontsize=14, weight="bold", color=COLORS["green"])

    all_data = box(ax, 0.5, 4.45, "All samples", w=2.2, kind="input")
    global_scale = box(ax, 3.2, 4.45, "Fit scaler and\nselect features\non all samples", w=2.5, h=1.0, kind="warning")
    split_bad = box(ax, 6.25, 4.45, "Then split", w=1.7, kind="process")
    inflated = box(ax, 3.35, 2.95, "Test fold has already\ninfluenced preprocessing", w=3.0, h=0.95, kind="warning")

    split_good = box(ax, 8.25, 4.45, "Split first", w=1.8, kind="process")
    train = box(ax, 10.65, 4.85, "Training samples", w=2.15, kind="input")
    test = box(ax, 10.65, 3.25, "Held-out samples", w=2.15, kind="input")
    fit_train = box(ax, 13.25, 4.85, "Fit preprocessing\non training only", w=2.15, h=0.8, kind="good")
    transform_test = box(ax, 13.25, 3.25, "Apply saved steps\nto held-out data", w=2.15, h=0.95, kind="good")

    arrow(ax, all_data, global_scale)
    arrow(ax, global_scale, split_bad)
    arrow_xy(ax, bottom(global_scale), top(inflated), curved=0.0, color=COLORS["rose"])
    arrow_xy(ax, right(split_bad), (8.2, 4.85), curved=0.0, color=COLORS["muted"])

    arrow_xy(ax, right(split_good), left(train), curved=-0.04)
    arrow_xy(ax, right(split_good), left(test), curved=0.04)
    arrow(ax, train, fit_train)
    arrow_xy(ax, bottom(fit_train), top(transform_test), curved=0.0)
    arrow(ax, test, transform_test)

    ax.text(
        0.7,
        0.85,
        "Data leakage often makes published performance look better than it will be on new patients.",
        fontsize=12,
        color=COLORS["text"],
    )
    save(fig, "16-ml-data-leakage.png")


def cross_validation():
    fig, ax = setup((13, 6.0))
    ax.set_xlim(0, 13)
    ax.set_ylim(0, 6.0)

    n_folds = 5
    x0 = 2.4
    y0 = 4.25
    width = 1.45
    height = 0.5
    gap = 0.08

    ax.text(0.55, 5.55, "5-fold cross-validation", fontsize=15, weight="bold", color=COLORS["text"])
    ax.text(0.55, 5.15, "Each row uses a different fold as the held-out test set.", fontsize=11, color=COLORS["muted"])

    for row in range(n_folds):
        y = y0 - row * 0.75
        ax.text(0.65, y + height / 2, f"Fold {row + 1}", va="center", fontsize=11, color=COLORS["text"])
        for col in range(n_folds):
            color = COLORS["warning"] if col == row else COLORS["good"]
            label = "Test" if col == row else "Train"
            rect = Rectangle((x0 + col * (width + gap), y), width, height, facecolor=color, edgecolor=COLORS["accent"], linewidth=1.2)
            ax.add_patch(rect)
            ax.text(x0 + col * (width + gap) + width / 2, y + height / 2, label, ha="center", va="center", fontsize=10, color=COLORS["text"])

    scores = box(ax, 10.45, 2.9, "Average the\n5 test scores", w=1.95, h=1.0, kind="output")
    arrow_xy(ax, (x0 + 5 * (width + gap) + 0.25, 2.75), left(scores), curved=0.0)

    ax.text(
        0.65,
        0.45,
        "For classification, stratified folds keep class proportions similar in every test fold.",
        fontsize=11,
        color=COLORS["text"],
    )
    save(fig, "16-ml-cross-validation.png")


def metrics():
    fig = plt.figure(figsize=(15, 5.4), constrained_layout=True)
    fig.patch.set_facecolor("white")
    gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.15, 1.15])
    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[0, 2])

    matrix = np.array([[42, 8], [5, 15]])
    im = ax0.imshow(matrix, cmap="YlGn", vmin=0, vmax=45)
    for i in range(2):
        for j in range(2):
            ax0.text(j, i, str(matrix[i, j]), ha="center", va="center", fontsize=17, weight="bold", color=COLORS["text"])
    ax0.set_xticks([0, 1], labels=["Predicted\nnegative", "Predicted\npositive"], fontsize=10)
    ax0.set_yticks([0, 1], labels=["Actual\nnegative", "Actual\npositive"], fontsize=10)
    ax0.set_title("Confusion Matrix", fontsize=14, weight="bold")
    ax0.tick_params(length=0)
    for spine in ax0.spines.values():
        spine.set_visible(False)
    fig.colorbar(im, ax=ax0, fraction=0.046, pad=0.04)

    fpr = np.array([0, 0.03, 0.08, 0.18, 0.32, 0.55, 1.0])
    tpr = np.array([0, 0.35, 0.6, 0.76, 0.86, 0.94, 1.0])
    ax1.plot(fpr, tpr, color=COLORS["green"], linewidth=2.4, label="Model")
    ax1.plot([0, 1], [0, 1], color=COLORS["muted"], linestyle="--", linewidth=1.4, label="Random")
    ax1.fill_between(fpr, tpr, fpr, color=COLORS["green"], alpha=0.12)
    ax1.set_xlabel("False positive rate", fontsize=11)
    ax1.set_ylabel("True positive rate / recall", fontsize=11)
    ax1.set_title("ROC Curve", fontsize=14, weight="bold")
    ax1.legend(frameon=False, fontsize=10, loc="lower right")

    recall = np.array([0.0, 0.15, 0.35, 0.55, 0.75, 0.9, 1.0])
    precision = np.array([1.0, 0.93, 0.82, 0.67, 0.5, 0.32, 0.2])
    baseline = 0.2
    ax2.plot(recall, precision, color=COLORS["orange"], linewidth=2.4, label="Model")
    ax2.axhline(baseline, color=COLORS["muted"], linestyle="--", linewidth=1.4, label="Positive rate")
    ax2.fill_between(recall, precision, baseline, color=COLORS["orange"], alpha=0.12)
    ax2.set_xlabel("Recall", fontsize=11)
    ax2.set_ylabel("Precision", fontsize=11)
    ax2.set_title("Precision-Recall Curve", fontsize=14, weight="bold")
    ax2.legend(frameon=False, fontsize=10, loc="upper right")

    for ax in [ax1, ax2]:
        ax.grid(True, color=COLORS["grid"], linewidth=0.8)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1.02)
        ax.tick_params(labelsize=10, colors=COLORS["text"])
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    save(fig, "16-ml-evaluation-metrics.png")


def bias_variance():
    fig, ax = plt.subplots(figsize=(12, 5.5), constrained_layout=True)
    fig.patch.set_facecolor("white")
    x = np.linspace(0, 10, 300)
    train = 0.65 * np.exp(-x / 2.8) + 0.08
    test = 0.045 * (x - 4.6) ** 2 + 0.28 + 0.06 * np.exp(-x / 1.5)
    ax.plot(x, train, color=COLORS["green"], linewidth=2.6, label="Training error")
    ax.plot(x, test, color=COLORS["rose"], linewidth=2.6, label="Test error")
    best_x = x[np.argmin(test)]
    best_y = test.min()
    ax.axvline(best_x, color=COLORS["accent"], linestyle="--", linewidth=1.5)
    ax.scatter([best_x], [best_y], s=80, color=COLORS["accent"], zorder=5)
    ax.text(best_x + 0.2, best_y + 0.08, "Best balance", fontsize=11, color=COLORS["text"])
    ax.text(0.55, 1.05, "Underfitting\nhigh bias", ha="center", fontsize=11, color=COLORS["text"])
    ax.text(9.1, 1.05, "Overfitting\nhigh variance", ha="center", fontsize=11, color=COLORS["text"])
    ax.set_xlabel("Model complexity", fontsize=12)
    ax.set_ylabel("Prediction error", fontsize=12)
    ax.set_title("Bias-Variance Tradeoff", fontsize=15, weight="bold")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 1.25)
    ax.grid(True, color=COLORS["grid"], linewidth=0.8)
    ax.legend(frameon=False, fontsize=11, loc="upper center")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=10, colors=COLORS["text"])
    save(fig, "16-ml-bias-variance.png")


def feature_selection_cv():
    fig, ax = setup((15, 6.3))
    ax.set_xlim(0, 15)
    ax.set_ylim(0, 6.3)

    ax.text(2.7, 5.95, "Wrong: select features before CV", ha="center", fontsize=14, weight="bold", color=COLORS["rose"])
    ax.text(11.0, 5.95, "Right: put feature selection inside the pipeline", ha="center", fontsize=14, weight="bold", color=COLORS["green"])

    data = box(ax, 0.45, 4.15, "All genes\nall samples", w=2.0, h=0.9, kind="input")
    select_all = box(ax, 3.05, 4.15, "Select top genes\nusing all labels", w=2.25, h=0.9, kind="warning")
    cv_bad = box(ax, 5.9, 4.15, "Cross-validate\nclassifier", w=2.05, h=0.9, kind="process")
    bad_note = box(ax, 2.25, 2.35, "Test folds helped\nchoose the genes", w=2.8, h=0.9, kind="warning")

    fold = box(ax, 8.45, 4.15, "For each fold", w=1.8, h=0.9, kind="input")
    scaler = box(ax, 10.75, 4.6, "Fit scaler\non train fold", w=2.0, h=0.85, kind="good")
    selector = box(ax, 10.75, 3.45, "Select genes\non train fold", w=2.0, h=0.85, kind="good")
    model = box(ax, 10.75, 2.3, "Train model\non train fold", w=2.0, h=0.85, kind="model")
    score = box(ax, 10.75, 1.15, "Score untouched\ntest fold", w=2.0, h=0.85, kind="output")

    arrow(ax, data, select_all)
    arrow(ax, select_all, cv_bad)
    arrow_xy(ax, bottom(select_all), top(bad_note), color=COLORS["rose"])

    arrow(ax, fold, scaler)
    arrow_xy(ax, bottom(scaler), top(selector))
    arrow_xy(ax, bottom(selector), top(model))
    arrow_xy(ax, bottom(model), top(score))

    ax.text(
        0.55,
        0.65,
        "In scikit-learn, use Pipeline so every preprocessing step is refit inside each cross-validation fold.",
        fontsize=11,
        color=COLORS["text"],
    )
    save(fig, "16-ml-feature-selection-cv.png")


def class_imbalance():
    fig, ax = plt.subplots(figsize=(12.5, 5.3), constrained_layout=True)
    fig.patch.set_facecolor("white")
    categories = ["Benign\nvariants", "Pathogenic\nvariants"]
    counts = [950, 50]
    colors = [COLORS["output"], COLORS["warning"]]
    bars = ax.bar(categories, counts, color=colors, edgecolor=COLORS["accent"], linewidth=1.4, width=0.55)
    for bar, value in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 25, f"{value}", ha="center", fontsize=13, weight="bold", color=COLORS["text"])
    ax.axhline(950, color=COLORS["muted"], linewidth=0.8, alpha=0.3)
    ax.set_ylabel("Number of samples", fontsize=12)
    ax.set_title("Class Imbalance Can Hide a Useless Classifier", fontsize=15, weight="bold")
    ax.text(
        0.42,
        0.62,
        "A model that predicts every variant as benign\nhas 95% accuracy but 0% recall for pathogenic variants.",
        transform=ax.transAxes,
        fontsize=12,
        color=COLORS["text"],
        bbox=dict(boxstyle="round,pad=0.35", facecolor=COLORS["warning"], edgecolor=COLORS["accent"], linewidth=1.2),
    )
    ax.set_ylim(0, 1100)
    ax.grid(axis="y", color=COLORS["grid"], linewidth=0.8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=11, colors=COLORS["text"])
    save(fig, "16-ml-class-imbalance.png")


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    ml_workflow()
    leakage()
    cross_validation()
    metrics()
    bias_variance()
    feature_selection_cv()
    class_imbalance()
