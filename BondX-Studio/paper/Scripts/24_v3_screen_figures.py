#!/usr/bin/env python3
"""24_v3_screen_figures.py - composes Fig. 2 (interface) and Fig. S11 (screens) from the screenshots of script 23.
Usage: python3 24_v3_screen_figures.py .."""
import os, sys
import numpy as np
from PIL import Image
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bxstyle3 import plt, save
ROOT = sys.argv[1] if len(sys.argv) > 1 else '..'
SH = os.path.join(ROOT, 'MS_Figures_V3', 'shots'); MS = os.path.join(ROOT, 'MS_Figures_V3'); ESI = os.path.join(ROOT, 'ESI_Figures_V3')
def main_area(f, top=80, h_css=760):
    im = Image.open(os.path.join(SH, f)).convert('RGB'); sc = im.width / 1600
    return im.crop((int(262 * sc), int(top * sc), im.width - int(16 * sc), int(h_css * sc)))
A = Image.open(os.path.join(SH, 'a_detect_clip.png')).convert('RGB'); A = A.crop((0, 0, A.width, 1200)); Bm = Image.open(os.path.join(SH, 'b_qm_clip.png')).convert('RGB'); Bm = Bm.crop((0, 0, Bm.width, 1500)); Cm = Image.open(os.path.join(SH, 'c_nbo_card.png')).convert('RGB')
# compact layout: three bands stacked, each cropped to the part that the caption describes
A = A.crop((0, 0, A.width, 700)); Bm = Bm.crop((0, 0, Bm.width, 640))
fig = plt.figure(figsize=(7.2, 7.4)); gs = fig.add_gridspec(3, 1, height_ratios=[A.height / A.width, Bm.height / Bm.width, Cm.height / Cm.width], hspace=0.16)
for k, (im, t) in enumerate([(A, '(a) T1: halogen contacts of 1US0'), (Bm, '(b) T3a: candidates and QM/MM set-up'), (Cm, '(c) T3b: NBO screen of a real 1US0 output')]):
    ax = fig.add_subplot(gs[k]); ax.imshow(im); ax.axis('off'); ax.set_title(t, fontsize=9.5, loc='left', fontweight='bold')
save(fig, MS, 'Fig2_interface')
titles = {'s01': 'T0: load structure', 's02': 'T1: halogen contacts', 's03': 'T3a: QM/MM model and charges (1US0)', 's04': 'T3b: NBO screen',
          's05': 'All parameters', 's06': 'Method page', 's07': 'Tutorial, page 1', 's08': 'Tutorial, page 5'}
files = sorted(f for f in os.listdir(SH) if f[:3] in titles)
fig, axs = plt.subplots(4, 2, figsize=(7.2, 11.0))
for ax in np.ravel(axs): ax.axis('off')
for ax, f in zip(np.ravel(axs), files):
    im = Image.open(os.path.join(SH, f)).convert('RGB'); w, h = im.size
    ax.imshow(im.crop((0, 0, w, min(h, int(w * 0.66))))); ax.set_title(titles[f[:3]], fontsize=10, fontweight='bold', loc='left')
plt.subplots_adjust(wspace=0.04, hspace=0.16); save(fig, ESI, 'FigS11_screens')
