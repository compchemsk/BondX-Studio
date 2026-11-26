"""bxstyle3.py - figure style of the V3 paper: Arial (Liberation Sans, the metric-identical Arial clone, where Arial is not
installed), large labels, legends without frames (frameon=False), validated categorical palette, PDF (fonts embedded) + PNG."""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

FONT = 'Arial' if any(f.name == 'Arial' for f in font_manager.fontManager.ttflist) else 'Liberation Sans'
C = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948']
HALC = {'CL': C[2], 'BR': C[1], 'I': C[0], 'F': '#8a8a86'}
VCOL = {'supported': '#2a78d6', 'inconclusive': '#eda100', 'not supported': '#e34948'}
GREY = '#f3e3b5'   # grey-zone band 0.3-0.8 kcal/mol
INK, INK2, GRID, MUTE = '#111111', '#3d3c38', '#e6e5e1', '#9a9993'
FS = 11
RC = {'font.family': 'sans-serif', 'font.sans-serif': [FONT, 'Arial', 'Liberation Sans'], 'mathtext.fontset': 'custom',
      'mathtext.rm': FONT, 'mathtext.it': FONT + ':italic', 'mathtext.bf': FONT + ':bold', 'mathtext.sf': FONT,
      'font.size': FS, 'axes.labelsize': FS + 0.5, 'axes.titlesize': FS + 1, 'xtick.labelsize': FS - 0.5, 'ytick.labelsize': FS - 0.5,
      'legend.fontsize': FS - 1, 'legend.frameon': False, 'axes.titleweight': 'bold', 'axes.titlelocation': 'left',
      'axes.edgecolor': '#8a8983', 'axes.linewidth': 0.9, 'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2,
      'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.7, 'axes.axisbelow': True, 'axes.spines.top': False,
      'axes.spines.right': False, 'lines.linewidth': 1.8, 'lines.markersize': 6, 'savefig.dpi': 300, 'pdf.fonttype': 42, 'ps.fonttype': 42}
plt.rcParams.update(RC)


def save(fig, folder, name):
    os.makedirs(folder, exist_ok=True)
    for ax in fig.axes:
        leg = ax.get_legend()
        if leg is not None: leg.set_frame_on(False)
    for leg in fig.legends: leg.set_frame_on(False)
    fig.savefig(os.path.join(folder, name + '.pdf'), bbox_inches='tight')
    fig.savefig(os.path.join(folder, name + '.png'), bbox_inches='tight', dpi=300)
    plt.close(fig); print('wrote', name)


def grey(ax, axis='y', lo=0.3, hi=0.8, label=True):
    (ax.axhspan if axis == 'y' else ax.axvspan)(lo, hi, color=GREY, zorder=0, lw=0, label='inconclusive (0.3–0.8)' if label else None)
