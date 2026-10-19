"""Figures for the RECOMB submission.

Design rules followed here:
  * no titles or explanatory text inside the axes; that goes in the caption,
  * large type throughout (base 11 pt at the final printed size),
  * each figure is created at its final width in inches so that LaTeX includes
    it at scale 1 and the type comes out at the size set here,
  * minimal in-figure text: axis labels, tick labels and a short legend.
"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import series as S

OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)

W = 6.5  # text width in inches (letter paper, 1 inch margins)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "xtick.labelsize": 10.5,
    "ytick.labelsize": 10.5,
    "legend.fontsize": 10.5,
    "axes.linewidth": 0.9,
    "lines.linewidth": 2.0,
    "lines.markersize": 6,
    "xtick.major.width": 0.9,
    "ytick.major.width": 0.9,
    "legend.frameon": False,
    "figure.dpi": 200,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
})

C = {
    "oneil":  "#D1660F",
    "dcomb":  "#1B6CA8",
    "b":      "#1B6CA8",
    "eta":    "#1E8A4C",
    "rank":   "#8E44AD",
    "frac":   "#D1660F",
    "grey":   "#666666",
    "band":   "#C44E52",
    "star":   "#1B6CA8",
    "both":   "#6C3483",
}

METH = {
    "null0":    ("none",      "#999999", "-",  "none"),
    "blind":    ("unif",      "#333333", "-",  "none"),
    "passive":  ("struct",    "#1E8A4C", "-",  "none"),
    "band_hat": ("band",      "#1B6CA8", "-",  "none"),
    "lse":      ("adapt",     "#D1660F", "-",  "none"),
    "ig":       ("infogain",  "#8E44AD", "--", "none"),
    "band_inf": ("oracle",    "#C44E52", "-",  "none"),
}


def _logy(ax, lo=0.8, hi=120):
    ax.set_yscale("log")
    ax.set_ylim(lo, hi)
    ax.set_yticks([1, 2, 4, 8, 16, 32, 64])
    ax.set_yticklabels(["1", "2", "4", "8", "16", "32", "64"])
    ax.axhline(1.0, color="black", lw=1.0, zorder=1)


def _tidy(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


# ---------------------------------------------------------------- figure 1 --
def fig_bands():
    """Why the two boundary sets differ: a scatter of prior against truth."""
    rng = np.random.default_rng(0)
    P, m, r, eta, b, q = 300, 60, 8, 0.50, 0.30, 0.90
    U = rng.standard_normal((P, r))
    V = rng.standard_normal((m, r))
    F = (U @ V.T) / np.sqrt(r)
    F = ((F - F.mean()) / F.std()).ravel()
    L = F + eta * rng.standard_normal(F.size)
    tau = np.quantile(F, q)

    star = np.abs(F - tau) <= b
    band = np.abs(L - tau) <= b
    both = star & band

    fig, axes = plt.subplots(1, 2, figsize=(W, 2.65))

    ax = axes[0]
    sel = rng.choice(F.size, 6000, replace=False)
    f, l = F[sel], L[sel]
    s_, b_, o_ = star[sel], band[sel], both[sel]
    rest = ~(s_ | b_)
    ax.scatter(f[rest], l[rest], s=1.5, c="#D5D5D5", linewidths=0, rasterized=True)
    ax.scatter(f[s_ & ~o_], l[s_ & ~o_], s=2.5, c=C["star"], linewidths=0, rasterized=True)
    ax.scatter(f[b_ & ~o_], l[b_ & ~o_], s=2.5, c=C["band"], linewidths=0, rasterized=True)
    ax.scatter(f[o_], l[o_], s=3.0, c=C["both"], linewidths=0, rasterized=True)
    ax.axvspan(tau - b, tau + b, color=C["star"], alpha=0.10, lw=0)
    ax.axhspan(tau - b, tau + b, color=C["band"], alpha=0.10, lw=0)
    ax.set_xlim(tau - 2.0, tau + 2.0)
    ax.set_ylim(tau - 2.0, tau + 2.0)
    ax.set_xlabel(r"withheld value $F_e$  ($\sigma_F$ from $\tau$)")
    ax.set_ylabel(r"reference prior $L_e$")
    ax.set_xticks([tau - 2, tau - 1, tau, tau + 1, tau + 2])
    ax.set_xticklabels(["-2", "-1", r"$\tau$", "1", "2"])
    ax.set_yticks([tau - 2, tau - 1, tau, tau + 1, tau + 2])
    ax.set_yticklabels(["-2", "-1", r"$\tau$", "1", "2"])
    ax.legend(handles=[
        Line2D([], [], marker="o", ls="", color=C["star"], label=r"$\mathcal{B}^{*}$ only"),
        Line2D([], [], marker="o", ls="", color=C["band"], label=r"$\mathcal{B}$ only"),
        Line2D([], [], marker="o", ls="", color=C["both"], label="both"),
    ], loc="upper left", handletextpad=0.2, borderpad=0.2, labelspacing=0.25)
    _tidy(ax)

    # Both knobs are measured in units of sigma_F, so one axis carries both.
    ax = axes[1]
    pts = S.SYN_SWEEP["b"]
    ax.plot([p[0] for p in pts], [p[1] for p in pts], "o-", color=C["b"],
            label=r"band half-width $b$")
    pts = S.SYN_SWEEP["eta"]
    ax.plot([p[0] for p in pts], [p[1] for p in pts], "^-", color=C["eta"],
            label=r"prior error $\eta$")
    ax.set_xlabel(r"$b$ or $\eta$  (units of $\sigma_F$)")
    ax.set_ylabel(r"overlap $|\mathcal{B}\cap\mathcal{B}^{*}|/M^{*}$")
    ax.set_ylim(0, 1)
    ax.set_xlim(0, 1.05)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.legend(loc="upper right", handletextpad=0.4, borderpad=0.2,
              labelspacing=0.25)
    _tidy(ax)

    fig.tight_layout(w_pad=2.0)
    fig.savefig(os.path.join(OUT, "bands.pdf"))
    plt.close(fig)


# ---------------------------------------------------------------- figure 2 --
def fig_effect():
    """The inflation factor, on the screens and in the controlled study."""
    fig, axes = plt.subplots(1, 3, figsize=(W, 2.5))

    ax = axes[0]
    for ds, lab in (("oneil", "O'Neil"), ("dcomb", "DrugComb")):
        pts = S.SCREEN_BAND[ds]
        ax.plot([p[0] for p in pts], [p[3] for p in pts], "o-", color=C[ds], label=lab)
    _logy(ax)
    ax.set_xlabel(r"band half-width $b$  ($\sigma_F$)")
    ax.set_ylabel("inflation factor")
    ax.set_xticks([0.05, 0.1, 0.2, 0.3, 0.5])
    ax.set_xticklabels(["0.05", "0.1", "0.2", "0.3", "0.5"])
    ax.legend(loc="upper right", handletextpad=0.4, borderpad=0.2, labelspacing=0.25)
    _tidy(ax)

    ax = axes[1]
    ax.plot(S.ALPHAS, S.SYN_BASE, "o-", color=C["frac"])
    _logy(ax)
    ax.set_xlabel(r"targeting strength $\alpha$")
    ax.set_xticks([0.25, 0.5, 0.75, 1.0])
    ax.set_xlim(0.15, 1.1)
    _tidy(ax)

    ax = axes[2]
    # rank is drawn as an open marker so the budget points beneath it, which sit
    # at the same overlap, stay visible.
    style = [("b", "s", r"band $b$", True), ("eta", "^", r"prior $\eta$", True),
             ("frac", "v", r"budget $T$", True), ("rank", "D", r"rank $r$", False)]
    for k, mk, lab, fill in style:
        pts = S.SYN_SWEEP[k]
        kw = dict(color=C[k], ls="none", label=lab)
        if not fill:
            kw.update(markerfacecolor="none", markeredgewidth=1.6, markersize=9,
                      zorder=5)
        ax.plot([p[1] for p in pts], [p[2] for p in pts], mk, **kw)
    _logy(ax)
    ax.set_xlabel(r"overlap $|\mathcal{B}\cap\mathcal{B}^{*}|/M^{*}$")
    ax.set_xlim(0, 0.95)
    ax.legend(loc="upper right", ncol=1, handletextpad=0.2,
              borderpad=0.2, labelspacing=0.2)
    _tidy(ax)

    for ax in axes[1:]:
        ax.set_yticklabels([])

    fig.tight_layout(w_pad=1.1)
    fig.savefig(os.path.join(OUT, "effect.pdf"))
    plt.close(fig)


# ---------------------------------------------------------------- figure 3 --
def fig_budget():
    """Boundary error against budget on the model-free set."""
    fig, axes = plt.subplots(1, 2, figsize=(W, 2.45))
    T = np.array(S.BUDGETS)
    order = ["null0", "blind", "passive", "band_hat", "lse", "ig", "band_inf"]
    for ax, ds, lab in zip(axes, ("oneil", "dcomb"), ("O'Neil", "DrugComb")):
        for m in order:
            y = S.CURVES[ds][m]
            xs = [t for t, v in zip(T, y) if v is not None]
            ys = [v for v in y if v is not None]
            name, col, ls, _ = METH[m]
            if m == "null0":
                ax.plot(xs, ys, ls="--", color=col, lw=1.6, zorder=1)
            else:
                ax.plot(xs, ys, ls=ls, color=col, marker="o", ms=3.5, zorder=2)
        ax.set_xscale("log")
        ax.set_xlabel(r"budget $T/N$")
        ax.set_xticks([0.01, 0.1, 1.0])
        ax.set_xticklabels(["0.01", "0.1", "1"])
        _tidy(ax)
    axes[0].set_ylabel(r"boundary error on $\mathcal{B}^{*}$")
    axes[0].set_ylim(0.07, 0.42)
    axes[1].set_ylim(0.12, 0.37)

    handles = [Line2D([], [], color=METH[m][1], lw=2.0,
                      ls="--" if m in ("null0", "ig") else "-",
                      label=METH[m][0]) for m in order]
    fig.legend(handles=handles, loc="lower center", ncol=7, frameon=False,
               handletextpad=0.4, columnspacing=1.1, bbox_to_anchor=(0.5, -0.03))
    fig.tight_layout(w_pad=1.8, rect=(0, 0.09, 1, 1))
    fig.savefig(os.path.join(OUT, "budget.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    fig_bands()
    fig_effect()
    fig_budget()
    print("wrote", sorted(os.listdir(OUT)))
