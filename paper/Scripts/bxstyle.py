"""bxstyle.py - one figure style for every figure of the BondX Studio paper.
Large labels, no legend frames, a validated categorical palette (blue, orange, aqua, yellow, magenta, green, violet, red),
PDF (vector) and PNG (300 dpi) output."""
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

C = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948']
HALC = {'CL': C[2], 'BR': C[1], 'I': C[0], 'F': '#8a8a86'}
INK, INK2, GRID, MUTE = '#0b0b0b', '#45443f', '#e4e3df', '#9a9993'
FS = 10.5          # base font size (pt) at the printed width
RC = {
    'font.family': 'DejaVu Sans', 'font.size': FS, 'axes.labelsize': FS + 0.5, 'axes.titlesize': FS + 1.0,
    'xtick.labelsize': FS - 0.5, 'ytick.labelsize': FS - 0.5, 'legend.fontsize': FS - 1.0, 'legend.frameon': False,
    'axes.titleweight': 'bold', 'axes.titlelocation': 'left', 'axes.edgecolor': '#8a8983', 'axes.linewidth': 0.9,
    'axes.labelcolor': INK, 'xtick.color': INK2, 'ytick.color': INK2, 'xtick.major.width': 0.9, 'ytick.major.width': 0.9,
    'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.7, 'axes.axisbelow': True,
    'axes.spines.top': False, 'axes.spines.right': False, 'lines.linewidth': 1.8, 'lines.markersize': 5.5,
    'savefig.dpi': 300, 'pdf.fonttype': 42, 'ps.fonttype': 42, 'figure.dpi': 100, 'mathtext.fontset': 'dejavusans'}


def apply():
    plt.rcParams.update(RC)


apply()


def save(fig, folder, name):
    os.makedirs(folder, exist_ok=True)
    for ax in fig.axes:
        leg = ax.get_legend()
        if leg is not None: leg.set_frame_on(False)
    for leg in fig.legends: leg.set_frame_on(False)
    fig.savefig(os.path.join(folder, name + '.pdf'), bbox_inches='tight')
    fig.savefig(os.path.join(folder, name + '.png'), bbox_inches='tight', dpi=300)
    plt.close(fig)
    print('wrote', name)
