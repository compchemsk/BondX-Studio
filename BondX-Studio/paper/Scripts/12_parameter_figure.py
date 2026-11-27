#!/usr/bin/env python3
"""Figure: every geometric parameter of BondX tier T1 in one picture.

Panel a: atom acceptor  A-X...Y-R  (d, theta1, theta2, r_X, r_Y, R, sigma-hole)
Panel b: ring acceptor  A-X...centroid  (d, theta1, alpha, ring normal, pseudo-radius)
Panel c: the three score terms g_d, g_theta, g_2 and the score G (default values)

Usage:  python 12_parameter_figure.py [output_dir]
All default values are the BondX 2.0 defaults (engine.js, DEFAULTS.geom and DEFAULTS.score).
"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge, Polygon, FancyArrowPatch, Arc, Ellipse

OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "..", "MS_Figures_V3")

# ---- BondX defaults (Br...O example) -------------------------------------------------------
TH1_MIN, TH1_MAX = 140.0, 180.0
TH2_OPT, TH2_W = 120.0, 40.0
TH2_WIN = (90.0, 150.0)          # optional theta2 filter (off by default)
ALPHA_MAX = 60.0
R_BEST = 0.85
R_BR, R_O, R_RING = 1.85, 1.52, 1.77
DCUT_BR_O = 3.68
DCUT_BR_AR = 3.87
RCUT = DCUT_BR_O / (R_BR + R_O)

# ---- palette ----------------------------------------------------------------------------------
BLUE, ORANGE, GREEN, AMBER, PINK, DGREEN, VIOLET, RED = (
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
INK, INK2, GRID = "#1f2328", "#57606a", "#d0d7de"
ATOM = {"C": "#8a8f96", "Br": "#a0522d", "O": "#e34948", "H": "#f4f4f4", "X": "#a0522d"}

from matplotlib import font_manager as _fm
_F = 'Arial' if any(f.name == 'Arial' for f in _fm.fontManager.ttflist) else 'Liberation Sans'
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": [_F, "Arial", "Liberation Sans"], "font.size": 10, "axes.linewidth": 0.8,
    "mathtext.fontset": "custom", "mathtext.rm": _F, "mathtext.it": _F + ":italic", "mathtext.bf": _F + ":bold", "axes.edgecolor": INK2, "legend.frameon": False, "pdf.fonttype": 42,
})


def unit(deg):
    a = np.radians(deg)
    return np.array([np.cos(a), np.sin(a)])


def atom(ax, p, el, r=0.30, label=None, z=5, lc=INK, fs=10.3):
    ax.add_patch(Circle(p, r, fc=ATOM.get(el, "#999"), ec=INK, lw=0.8, zorder=z))
    # small highlight for a 3-D look
    ax.add_patch(Circle(p + np.array([-0.30, 0.32]) * r, 0.33 * r, fc="white", ec="none", alpha=0.45, zorder=z + 0.1))
    if label:
        ax.text(*p, label, ha="center", va="center", fontsize=fs, color=lc, weight="bold", zorder=z + 0.2)


def bond(ax, p, q, r1=0.30, r2=0.30, lw=3.0, color=INK2, dashed=False, z=3, double=False):
    v = q - p; L = np.linalg.norm(v); u = v / L
    a, b = p + u * r1 * 0.8, q - u * r2 * 0.8
    if dashed:
        ax.plot(*zip(a, b), ls=(0, (1.5, 2.2)), lw=1.6, color=color, zorder=z, solid_capstyle="round")
    elif double:
        n = np.array([-u[1], u[0]]) * 0.085
        for s in (1, -1):
            ax.plot(*zip(a + s * n, b + s * n), lw=lw * 0.55, color=color, zorder=z, solid_capstyle="round")
    else:
        ax.plot(*zip(a, b), lw=lw, color=color, zorder=z, solid_capstyle="round")


def angle_arc(ax, c, a1, a2, r, color, label, lr=None, lw=1.4, fs=11.5, ha="center"):
    lo, hi = min(a1, a2), max(a1, a2)
    if hi - lo > 180:
        lo, hi = hi, lo + 360
    ax.add_patch(Arc(c, 2 * r, 2 * r, theta1=lo, theta2=hi, color=color, lw=lw, zorder=6))
    mid = 0.5 * (lo + hi)
    lr = lr or r + 0.28
    ax.text(*(c + unit(mid) * lr), label, color=color, fontsize=fs, ha=ha, va="center", zorder=7, weight="bold")


def dim_line(ax, p, q, off, color, label, fs=10.9, lpos=0.5, loff=0.22):
    v = q - p; u = v / np.linalg.norm(v); n = np.array([-u[1], u[0]])
    a, b = p + n * off, q + n * off
    for e in (p, q):
        ax.plot(*zip(e + n * off * 0.25, e + n * off * 1.12), color=color, lw=0.7, zorder=4)
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="<|-|>", mutation_scale=8, color=color, lw=1.0, zorder=4))
    m = a + (b - a) * lpos + n * np.sign(off) * loff
    rot = np.degrees(np.arctan2(u[1], u[0]))
    ax.text(*m, label, color=color, fontsize=fs, ha="center", va="center", rotation=rot, zorder=7)


def sigma_hole(ax, X, direction, rX=0.42):
    """Blue positive cap along the C-X axis and a negative belt around X."""
    ax.add_patch(Wedge(X, rX + 0.20, direction - 32, direction + 32, width=0.20, fc=BLUE, ec="none", alpha=0.85, zorder=4))
    for c in (direction + 90, direction - 90):
        ax.add_patch(Wedge(X, rX + 0.16, c - 42, c + 42, width=0.14, fc=RED, ec="none", alpha=0.55, zorder=4))


def panel_label(ax, s):
    ax.text(0.0, 1.0, s, transform=ax.transAxes, fontsize=13.8, weight="bold", va="top", ha="left", color=INK)


# ==============================================================================================
fig = plt.figure(figsize=(7.2, 7.0))
gs = fig.add_gridspec(2, 2, height_ratios=[1.45, 1.0], width_ratios=[1.38, 1.0], hspace=0.20, wspace=0.04,
                      left=0.01, right=0.99, top=0.985, bottom=0.075)

# ---------------------------------------------------------------- (a) atom acceptor
ax = fig.add_subplot(gs[0, 0]); ax.set_aspect("equal"); ax.axis("off")
panel_label(ax, "a")
ax.text(0.06, 0.995, "Atom acceptor   A–X···Y–R", transform=ax.transAxes, fontsize=10.9, va="top", color=INK, weight="bold")

th1, th2, d = 166.0, 120.0, 3.05
X = np.array([0.0, 0.0])
dirXY = 0.0                                   # direction X -> Y
dirXA = dirXY + th1                            # theta1 = angle(A, X, Y)
A = X + unit(dirXA) * 1.95
Y = X + unit(dirXY) * d
dirYX = dirXY + 180.0
dirYR = dirYX - th2                            # theta2 = angle(X, Y, R); R above-right of Y
R = Y + unit(dirYR) * 1.25

# substituents
H1 = A + unit(dirXA + 70) * 1.0; H2 = A + unit(dirXA - 70) * 1.0; H3 = A + unit(dirXA + 180 + 180) * 0  # placeholder
ring_dir = dirXA
Cs = [A + unit(dirXA + 60) * 1.25, A + unit(dirXA - 60) * 1.25]
R2 = R + unit(dirYR + 55) * 1.05; R3 = R + unit(dirYR - 72) * 1.0

# van der Waals spheres (dashed) and overlap
ax.add_patch(Circle(X, R_BR, fc=ORANGE, alpha=0.06, ec=ORANGE, ls=(0, (3, 2)), lw=0.9, zorder=1))
ax.add_patch(Circle(Y, R_O, fc=RED, alpha=0.06, ec=RED, ls=(0, (3, 2)), lw=0.9, zorder=1))

for c in Cs:
    bond(ax, A, c); atom(ax, c, "C", r=0.22)
bond(ax, R, R2); bond(ax, R, R3); atom(ax, R2, "C", r=0.22); atom(ax, R3, "C", r=0.22)
bond(ax, A, X, r1=0.30, r2=0.42)
bond(ax, Y, R, r1=0.30, r2=0.30, double=True, lw=3.4)
bond(ax, X, Y, r1=0.42, r2=0.30, dashed=True, color=INK)
sigma_hole(ax, X, dirXY)
atom(ax, A, "C", label="A")
atom(ax, X, "X", r=0.42, label="X", lc="white", fs=12.6)
atom(ax, Y, "O", r=0.32, label="Y", lc="white")
atom(ax, R, "C", label="R")

# lone pair on Y pointing at X (sp2, 120 deg)
for s in (+1,):
    lp = Y + unit(dirYX) * 0.0
ax.text(*(X + unit(dirXY + 40) * 0.80), "σ-hole (+)", color=BLUE, fontsize=9.2, ha="left", va="center", weight="bold")
ax.text(*(X + unit(dirXY - 90) * 0.70), "belt (−)", color=RED, fontsize=9.2, ha="center", va="top")

# angles
angle_arc(ax, X, dirXA, dirXY, 0.80, BLUE, "θ₁", lr=1.12)
angle_arc(ax, Y, dirYR, dirYX, 0.62, GREEN, "θ₂", lr=0.95)

# distance
dim_line(ax, X, Y, -1.25, INK, "d (X···Y)", lpos=0.55)

# radii
ax.annotate("", xy=X + unit(-162) * R_BR, xytext=X, arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=1.0, mutation_scale=8), zorder=6)
ax.text(*(X + unit(-162) * R_BR * 0.62 + np.array([0.0, -0.30])), "r$_X$", color=ORANGE, fontsize=10.9, ha="center", va="center")
ax.annotate("", xy=Y + unit(-35) * R_O, xytext=Y, arrowprops=dict(arrowstyle="-|>", color=RED, lw=1.0, mutation_scale=8), zorder=6)
ax.text(*(Y + unit(-35) * R_O * 0.55 + np.array([0.12, -0.30])), "r$_Y$", color=RED, fontsize=10.9, ha="center", va="center")

# rule box
txt = (r"$d \leq d_\mathrm{cut}$(X,Y)" "\n"
       r"$140° \leq \theta_1 \leq 180°$" "\n"
       r"$\theta_2$: score optimum $120°$" "\n"
       r"   (filter $90$–$150°$, optional)" "\n"
       r"$R = d\,/\,(r_X + r_Y)$")
ax.text(0.985, 0.012, txt, transform=ax.transAxes, fontsize=9.2, ha="right", va="bottom", color=INK, linespacing=1.45,
        bbox=dict(boxstyle="round,pad=0.45", fc="#f6f8fa", ec=GRID, lw=0.8))
ax.set_xlim(-3.15, 4.95); ax.set_ylim(-3.75, 2.95)

# ---------------------------------------------------------------- (b) ring acceptor
ax = fig.add_subplot(gs[0, 1]); ax.set_aspect("equal"); ax.axis("off")
panel_label(ax, "b")
ax.text(0.09, 0.995, "π acceptor\nA–X···Ct", transform=ax.transAxes, fontsize=10.9, va="top", color=INK, weight="bold", linespacing=1.3)

Ct = np.array([0.0, -1.35]); sq = 0.36; rr = 1.40
ring = [Ct + np.array([rr * np.cos(np.radians(a)), sq * rr * np.sin(np.radians(a))]) for a in range(0, 360, 60)]
ax.add_patch(Ellipse(Ct, 2 * R_RING * 1.0, 2 * R_RING * sq, fc=VIOLET, alpha=0.06, ec=VIOLET, ls=(0, (3, 2)), lw=0.9, zorder=1))
ax.add_patch(Polygon(ring, closed=True, fc="#e9e6fb", ec=INK2, lw=1.6, zorder=2, alpha=0.95))
ax.add_patch(Ellipse(Ct, 1.25 * rr, 1.25 * rr * sq, fc="none", ec=VIOLET, lw=0.9, zorder=2.5))  # aromatic circle
for p in ring:
    ax.add_patch(Circle(p, 0.15, fc=ATOM["C"], ec=INK, lw=0.6, zorder=3))
ax.add_patch(Circle(Ct, 0.08, fc=INK, ec="none", zorder=5))
ax.text(*(Ct + np.array([0.18, 0.17])), "Ct", fontsize=9.8, color=INK, ha="left", va="center")

# normal
ax.add_patch(FancyArrowPatch(Ct, Ct + np.array([0, 3.55]), arrowstyle="-|>", mutation_scale=9, color=INK2, lw=1.0, ls="--", zorder=3))
ax.text(*(Ct + np.array([-0.12, 3.58])), "n", fontsize=11.5, color=INK2, ha="right", va="bottom", style="italic", weight="bold")

alpha = 22.0; dpi = 3.15
dirCX = 90 - alpha
Xp = Ct + unit(dirCX) * dpi
dirXA2 = (dirCX + 180 + (180 - 168)) % 360  # theta1 = 168 deg
dirXA2 = dirCX + 12.0                       # A roughly continues the Ct->X line
Ap = Xp + unit(dirXA2) * 1.70
bond(ax, Ap, Xp, r1=0.28, r2=0.40)
bond(ax, Xp, Ct, r1=0.40, r2=0.05, dashed=True, color=INK)
sigma_hole(ax, Xp, dirCX + 180)
atom(ax, Ap, "C", r=0.28, label="A")
atom(ax, Xp, "X", r=0.40, label="X", lc="white", fs=12.6)
angle_arc(ax, Ct, 90, dirCX, 1.05, VIOLET, "α", lr=1.32)
angle_arc(ax, Xp, dirXA2, dirCX + 180, 0.82, BLUE, "θ₁", lr=1.10)
dim_line(ax, Ct, Xp, -0.62, INK, "d (X···Ct)", lpos=0.42, loff=0.25)
tip = Ct + np.array([R_RING * np.cos(np.radians(235)), R_RING * sq * np.sin(np.radians(235))])
ax.annotate("", xy=tip, xytext=Ct, arrowprops=dict(arrowstyle="-|>", color=VIOLET, lw=1.0, mutation_scale=8), zorder=6)
ax.text(*(Ct + np.array([-1.05, -0.95])), "r$_\\mathrm{ring}$ = 1.77 Å", color=VIOLET, fontsize=9.2, ha="center")

txt = (r"$d \leq d_\mathrm{cut}$(X, ring)" "\n"
       r"$140° \leq \theta_1 \leq 180°$" "\n"
       r"$\alpha \leq 60°$ (to ring normal)")
ax.text(0.03, 0.012, txt, transform=ax.transAxes, fontsize=9.2, ha="left", va="bottom", color=INK, linespacing=1.45,
        bbox=dict(boxstyle="round,pad=0.45", fc="#f6f8fa", ec=GRID, lw=0.8))
ax.set_xlim(-2.45, 2.75); ax.set_ylim(-4.6, 3.85)

# ---------------------------------------------------------------- (c) score terms
sub = gs[1, :].subgridspec(1, 4, width_ratios=[1, 1, 1, 1.45], wspace=0.40)
axc = [fig.add_subplot(sub[0, i]) for i in range(3)]
for a in axc:                       # leave room for the y label of the first plot
    b = a.get_position(); a.set_position([b.x0 + 0.075 - 0.018 * axc.index(a), b.y0, b.width * 0.90, b.height])
axc[0].text(-0.40, 1.16, "c", transform=axc[0].transAxes, fontsize=13.8, weight="bold", va="top", color=INK)


def style(a, xl, yl=True):
    a.set_ylim(-0.04, 1.1); a.set_xlabel(xl, fontsize=9.8, color=INK, labelpad=2)
    a.tick_params(labelsize=7.5, colors=INK2, length=2.5)
    a.grid(color=GRID, lw=0.5); a.set_axisbelow(True)
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
    if not yl:
        a.set_yticklabels([])


# g_d(R)
Rv = np.linspace(0.7, 1.2, 400)
gd = np.clip((RCUT - Rv) / (RCUT - R_BEST), 0, 1)
a = axc[0]; a.plot(Rv, gd, color=ORANGE, lw=2)
a.axvline(R_BEST, color=INK2, lw=0.7, ls=":"); a.axvline(RCUT, color=INK2, lw=0.7, ls=":")
a.text(R_BEST - 0.012, 0.42, "0.85", rotation=90, fontsize=8, color=INK2, ha="right", va="center")
a.text(RCUT + 0.012, 0.42, f"$R_\\mathrm{{cut}}$ = {RCUT:.2f}", rotation=90, fontsize=8, color=INK2, ha="left", va="center")
a.set_title("distance term  $g_d$", fontsize=9.8, color=ORANGE, pad=3, weight="bold")
style(a, "R = d / (r$_X$ + r$_Y$)   [Br···O]"); a.set_ylabel("term value", fontsize=9.8, color=INK)

# g_theta(theta1)
t = np.linspace(120, 180, 300)
gt = np.clip((t - TH1_MIN) / (180 - TH1_MIN), 0, 1)
a = axc[1]; a.plot(t, gt, color=BLUE, lw=2)
a.axvspan(TH1_MIN, TH1_MAX, color=BLUE, alpha=0.07, lw=0)
a.text(160, 1.04, "XB window", fontsize=8, color=BLUE, ha="center", va="bottom")
a.set_title("σ-hole angle term  $g_\\theta$", fontsize=9.8, color=BLUE, pad=3, weight="bold")
a.set_xticks([120, 140, 160, 180]); style(a, "θ₁ (°)", yl=False)

# g_2(theta2) and ring form
t2 = np.linspace(60, 180, 300)
g2 = np.exp(-((t2 - TH2_OPT) / TH2_W) ** 2)
a = axc[2]; a.plot(t2, g2, color=GREEN, lw=2, label="atom: θ₂")
al = np.linspace(0, 90, 200)
a2 = a.twiny()
a2.plot(al, np.clip(1 - al / ALPHA_MAX, 0, 1), color=VIOLET, lw=1.8, ls="--", label="ring: α")
a2.set_xlim(0, 90); a2.tick_params(labelsize=7.5, colors=VIOLET, length=2.5, pad=1)
a2.set_xticks([0, 30, 60, 90]); a2.set_xlabel("α (°)", fontsize=9.8, color=VIOLET, labelpad=2)
for s in ("right",):
    a2.spines[s].set_visible(False)
a2.spines["top"].set_color(VIOLET)
a.set_xticks([60, 90, 120, 150, 180]); style(a, "θ₂ (°)", yl=False)
a.xaxis.label.set_color(GREEN); a.tick_params(axis="x", colors=GREEN)
a.text(152, 0.86, "atom", fontsize=8.6, color=GREEN, ha="left", weight="bold")
a.text(50, 0.18, "1 − α/60°", fontsize=8, color=VIOLET, ha="left", transform=a.transData) if False else \
    a2.text(2, 0.30, "ring", fontsize=8.6, color=VIOLET, ha="left", weight="bold")
a.set_title("acceptor term  $g_2$", fontsize=9.8, color=GREEN, pad=15, weight="bold")

# formula / summary box
ab = fig.add_subplot(sub[0, 3]); ab.axis("off")
bb = ab.get_position(); ab.set_position([bb.x0 + 0.01, bb.y0, bb.width * 0.92, bb.height])
ab.add_patch(plt.Rectangle((-0.05, -0.04), 1.04, 1.14, transform=ab.transAxes, fc="#f6f8fa", ec=GRID, lw=0.8, clip_on=False))
formulas = (
    r"$g_d = \mathrm{clip}[(R_\mathrm{cut}-R)/(R_\mathrm{cut}-0.85)]$" "\n"
    r"$g_\theta = \mathrm{clip}[(\theta_1-140°)/40°]$" "\n"
    r"$g_2 = e^{-[(\theta_2-120°)/40°]^2}$  (atom)" "\n"
    r"$g_2 = 1-\alpha/60°$  (ring)" "\n"
    r"$G = [\,g_d\, g_\theta\, (\frac{1}{2}+\frac{1}{2}g_2)\,]^{1/3}$")
ab.text(-0.02, 1.07, formulas, transform=ab.transAxes, fontsize=8.1, va="top", ha="left", color=INK, linespacing=1.55)
table = ("Cut-off d$_\\mathrm{cut}$ (Å)\n"
         "    N     O     S     π\n"
         "Cl 3.77  3.61  3.95  3.80\n"
         "Br 3.84  3.68  4.02  3.87\n"
         "I  3.97  3.82  4.15  4.00")
ab.text(0.0, 0.47, table.split("\n")[0], transform=ab.transAxes, fontsize=8.5, va="bottom", ha="left", color=INK)
# the cut-off table in Arial: one text item per cell, right-aligned in fixed columns
_rows = [["", "N", "O", "S", "π"], ["Cl", "3.77", "3.61", "3.95", "3.80"], ["Br", "3.84", "3.68", "4.02", "3.87"], ["I", "3.97", "3.82", "4.15", "4.00"]]
_xs = [0.0, 0.22, 0.40, 0.58, 0.76]
for _i, _r in enumerate(_rows):
    for _j, _c in enumerate(_r):
        ab.text(_xs[_j] + (0.0 if _j == 0 else 0.12), 0.45 - _i * 0.085, _c, transform=ab.transAxes, fontsize=7.9, va="top",
                ha="left" if _j == 0 else "right", color=INK, fontweight="bold" if (_i == 0 or _j == 0) else "normal")
ab.text(0.96, -0.01, "all values editable", transform=ab.transAxes, fontsize=7.8, va="bottom", ha="right", color=INK2, style="italic")

os.makedirs(OUT, exist_ok=True)
for ext in ("pdf", "png"):
    fig.savefig(os.path.join(OUT, f"Fig3_T1_model.{ext}"), dpi=400)
print("wrote", os.path.join(OUT, "Fig3_T1_model.pdf/png"))
