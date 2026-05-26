"""
pTau SASA Publication Figure Script  v3
=========================================
Changes from v2:
  - Gradient bars (light tint → full color, left → right) for selected variants only
  - Dual-site hatch rendered via transparent Rectangle overlay on gradient bars
  - Colormap: plasma → YlOrRd  (pale yellow = borderline, deep red = strongly selected)
  - Terminology: 'phosphoform' removed throughout; 'Dual-site pTau variant' used
  - Panel A: 15×17", y-tick 13pt, value labels 12pt (not bold), title 17pt
  - Panel B: 11×9.5"
  - Cutoff annotation box moved above top bar (arrow from line to box), no overlap
  - PDF output removed; PNG only
  - adjustText arrows: shrinkA=6, shrinkB=4 (suppresses arrow-through-text warning)

Requirements:
    pip install pandas matplotlib openpyxl numpy adjustText
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.colors as mcolors
from matplotlib import rcParams
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Rectangle
from adjustText import adjust_text

# ── Typography ─────────────────────────────────────────────────────────────────
rcParams['font.family']    = 'Liberation Sans'
rcParams['pdf.fonttype']   = 42
rcParams['ps.fonttype']    = 42
rcParams['axes.linewidth'] = 0.9

# ── Config ─────────────────────────────────────────────────────────────────────
INPUT_FILE  = 'pTau_SASA.xlsx'
CUTOFF      = 50
TOP_N_LABEL = 8
BG          = 'white'
COL_REJ     = '#BBBBBB'
COL_CUT     = '#2CA02C'
CMAP_NAME   = 'YlOrRd'

# ── Load & sort numerically ────────────────────────────────────────────────────
df = pd.read_excel(INPUT_FILE)
df['_g'] = df['pTau_Variant'].str.extract(r'pTau_(\d+)\.').astype(int)
df['_v'] = df['pTau_Variant'].str.extract(r'pTau_\d+\.(\d+)').astype(int)
df = df.sort_values(['_g', '_v']).reset_index(drop=True)

# ── Label dual-site variants (i) / (ii) ───────────────────────────────────────
counts = {}
for name in df['pTau_Variant']:
    counts[name] = counts.get(name, 0) + 1

seen, labels, is_dual_list = {}, [], []
for name in df['pTau_Variant']:
    if counts[name] > 1:
        seen[name] = seen.get(name, 0) + 1
        labels.append(name + (' (i)' if seen[name] == 1 else ' (ii)'))
        is_dual_list.append(True)
    else:
        labels.append(name)
        is_dual_list.append(False)

df['label']   = labels
df['is_dual'] = is_dual_list

variants = df['label'].tolist()
ptau     = df['pTau_SASA'].tolist()
htau     = df['hTau_SASA'].tolist()
delta    = df['Δ_SASA'].tolist()
is_dual  = df['is_dual'].tolist()

n_rows_sel   = sum(d >= CUTOFF for d in delta)
n_unique_sel = df[df['Δ_SASA'] >= CUTOFF]['pTau_Variant'].nunique()

# ── Colormap ───────────────────────────────────────────────────────────────────
cmap_sel = plt.get_cmap(CMAP_NAME)
d_sel    = [d for d in delta if d >= CUTOFF]
norm     = mcolors.Normalize(vmin=min(d_sel), vmax=max(d_sel))


# ── Gradient bar helper ────────────────────────────────────────────────────────
def draw_gradient_bar(ax, x_end, y_center, bar_h, color_rgba,
                      light_blend=0.45, n_pts=512, zorder=3):
    """
    Draw a horizontal bar that fades from a lightened tint of `color_rgba`
    at x=0 (left) to the full color at x=x_end (right).

    light_blend : fraction of white mixed into the left-edge tint (0=no tint, 1=white)
    """
    white = np.ones(3)
    color = np.array(color_rgba[:3])
    t     = np.linspace(0, 1, n_pts)                       # 0=left, 1=right
    # Left edge: white*light_blend + color*(1-light_blend)
    # Right edge: color (full saturation)
    left_col = white * light_blend + color * (1.0 - light_blend)
    cols = left_col[np.newaxis, :] * (1 - t)[:, np.newaxis] \
         + color[np.newaxis, :] * t[:, np.newaxis]          # (n_pts, 3)
    grad_rgba = np.clip(
        np.concatenate([cols, np.ones((n_pts, 1))], axis=1),
        0, 1
    ).reshape(1, n_pts, 4)

    ax.imshow(grad_rgba, aspect='auto',
              extent=[0, x_end, y_center - bar_h / 2, y_center + bar_h / 2],
              zorder=zorder, interpolation='bilinear')


# =============================================================================
#  PANEL A — Horizontal bar chart
# =============================================================================
y     = np.arange(len(variants))
BAR_H = 0.70

figA, axA = plt.subplots(figsize=(15, 17))
figA.patch.set_facecolor(BG)
axA.set_facecolor(BG)

# Alternating row background
for i in range(len(variants)):
    axA.axhspan(i - 0.5, i + 0.5,
                color='#EFEFEF' if i % 2 == 0 else BG, zorder=0)

# ── Draw bars ──────────────────────────────────────────────────────────────────
for i, d in enumerate(delta):
    if d >= CUTOFF:
        # Gradient bar for selected variants
        draw_gradient_bar(axA, d, i, BAR_H,
                          cmap_sel(norm(d)), light_blend=0.45, zorder=3)
        # Hatch overlay for dual-site variants (transparent fill, hatch only)
        if is_dual[i]:
            hatch_rect = Rectangle(
                (0, i - BAR_H / 2), d, BAR_H,
                fill=False, hatch='//',
                edgecolor='#444444', linewidth=0, zorder=4
            )
            axA.add_patch(hatch_rect)
    else:
        # Solid grey for excluded; hatch if dual-site
        axA.barh(i, d, height=BAR_H,
                 color=COL_REJ, edgecolor='white', linewidth=0.5,
                 hatch='//' if is_dual[i] else None,
                 zorder=3)

# ── Value labels ───────────────────────────────────────────────────────────────
for i, d in enumerate(delta):
    axA.text(d + 1.2, i, f'{d:.1f}',
             va='center', ha='left',
             fontsize=12, color='#1A1A1A', zorder=5)

# ── Cutoff line ────────────────────────────────────────────────────────────────
axA.axvline(CUTOFF, color=COL_CUT, linewidth=2.2, linestyle='--', zorder=4)

# ── Cutoff annotation — above the topmost bar, connected by arrow ──────────────
axA.annotate(
    f'Δ SASA cutoff = {CUTOFF} Å²',
    xy=(CUTOFF, len(variants) - 0.5),          # arrow tip: top of dashed line
    xytext=(CUTOFF + 6, len(variants) + 0.7),  # box position: above top bar
    fontsize=13, fontweight='bold', color=COL_CUT,
    ha='left', va='bottom',
    bbox=dict(boxstyle='round,pad=0.45', facecolor='white',
              edgecolor=COL_CUT, linewidth=1.2, alpha=0.97),
    arrowprops=dict(arrowstyle='->', color=COL_CUT, lw=1.4)
)

# ── Bracket for selected region ────────────────────────────────────────────────
sel_idx = [i for i, d in enumerate(delta) if d >= CUTOFF]
bx = max(delta) * 1.12
axA.annotate('', xy=(bx, max(sel_idx)), xytext=(bx, min(sel_idx)),
             arrowprops=dict(arrowstyle='-[, widthB=0.55, lengthB=0.4',
                             color='#444499', lw=1.5))
axA.text(bx * 1.02, (max(sel_idx) + min(sel_idx)) / 2,
         f'Selected\nfor docking\n({n_unique_sel} pTau variants\n{n_rows_sel} sites)',
         color='#444499', fontsize=13, fontweight='bold',
         va='center', ha='left', linespacing=1.5)

# ── Colorbar ───────────────────────────────────────────────────────────────────
sm = ScalarMappable(cmap=cmap_sel, norm=norm)
sm.set_array([])
cbar = figA.colorbar(sm, ax=axA, pad=0.015, fraction=0.018, aspect=30)
cbar.set_label('Δ SASA (Å²)', fontsize=14, labelpad=9)
cbar.ax.tick_params(labelsize=12)
cbar.outline.set_linewidth(0.5)

# ── Axis formatting ────────────────────────────────────────────────────────────
axA.set_yticks(y)
axA.set_yticklabels(variants, fontsize=13, color='#1A1A1A')
axA.set_xlabel('Δ SASA  (pTau − hTau)  [Å²]',
               fontsize=15, color='#111111', labelpad=10)
axA.set_xlim(0, max(delta) * 1.50)
axA.set_ylim(-0.7, len(variants) + 1.8)   # top headroom for annotation box
axA.tick_params(axis='x', labelsize=13, length=4, colors='#333333')
axA.tick_params(axis='y', length=0)
axA.xaxis.grid(True, linestyle=':', linewidth=0.6, color='#DDDDDD', zorder=0)
axA.set_axisbelow(True)
for sp in ['top', 'right']:
    axA.spines[sp].set_visible(False)
for sp in ['left', 'bottom']:
    axA.spines[sp].set_color('#AAAAAA')
    axA.spines[sp].set_linewidth(0.9)

axA.set_title(
    'SASA Exposure Difference Across pTau Variants\n'
    '(Δ SASA = pTau − hTau, sorted by variant name)',
    fontsize=17, fontweight='bold', color='#111111',
    pad=16, loc='left', linespacing=1.5
)

# ── Legend ─────────────────────────────────────────────────────────────────────
sel_p = mpatches.Patch(color=cmap_sel(0.75),
                        label=f'Selected  (Δ SASA ≥ {CUTOFF} Å²)')
rej_p = mpatches.Patch(color=COL_REJ,
                        label=f'Excluded  (Δ SASA < {CUTOFF} Å²)')
dup_p = mpatches.Patch(facecolor='#CCCCCC', hatch='//', edgecolor='#444444',
                        label='Dual-site pTau variant (i / ii = individual sites)')
cut_l = plt.Line2D([0], [0], color=COL_CUT, lw=2.2, linestyle='--',
                   label=f'Cutoff = {CUTOFF} Å²')
axA.legend(handles=[sel_p, rej_p, dup_p, cut_l],
           fontsize=12, frameon=True, edgecolor='#BBBBBB',
           facecolor='white', framealpha=0.97,
           loc='upper right',
           borderpad=0.9, labelspacing=0.65, handlelength=1.7)

# ── Caption ────────────────────────────────────────────────────────────────────
figA.text(
    0.5, -0.004,
    f'n = {len(df)} total pTau variants  |  '
    f'{n_unique_sel} variants selected ({n_rows_sel} phosphorylation sites)  |  '
    f'Hatched bars = dual-site pTau variants; (i) and (ii) denote each phosphorylation site separately  |  '
    f'Color: {CMAP_NAME}',
    ha='center', fontsize=10, color='#888888', style='italic'
)

plt.tight_layout(pad=1.8)
figA.savefig('panelA_delta_sasa.png', dpi=300, bbox_inches='tight', facecolor=BG)
plt.close(figA)
print('Panel A saved → panelA_delta_sasa.png')


# =============================================================================
#  PANEL B — Scatter: pTau SASA vs hTau SASA
# =============================================================================
x_max = max(ptau) * 1.10
y_max = max(htau) * 1.18

figB, axB = plt.subplots(figsize=(11, 9.5))
figB.patch.set_facecolor(BG)
axB.set_facecolor('#FAFAFA')

x_fill = np.linspace(0, x_max, 600)

# Excluded zone: full region above the cutoff line (hTau > pTau − 50)
axB.fill_between(x_fill,
                 np.maximum(x_fill - CUTOFF, 0), y_max,
                 alpha=0.08, color='#AAAAAA', zorder=0)

# Identity line
axB.plot([0, x_max], [0, x_max],
         color='#CCCCCC', lw=1.3, linestyle='--', zorder=1)

# Δ SASA = CUTOFF boundary
axB.plot([CUTOFF, x_max], [0, x_max - CUTOFF],
         color=COL_CUT, lw=2.0, linestyle='--', zorder=2)

# Excluded points (below threshold)
for i in range(len(variants)):
    if delta[i] < CUTOFF:
        marker = 's' if is_dual[i] else 'o'
        axB.scatter(ptau[i], htau[i], color=COL_REJ, s=80, zorder=3,
                    edgecolors='white', linewidths=0.8, alpha=0.55, marker=marker)

# Selected points (above threshold)
for i in range(len(variants)):
    if delta[i] >= CUTOFF:
        marker = 's' if is_dual[i] else 'o'
        axB.scatter(ptau[i], htau[i],
                    color=cmap_sel(norm(delta[i])),
                    s=140, zorder=5,
                    edgecolors='white', linewidths=1.0, marker=marker)

# Labels for top N by Δ SASA
top_idx = sorted(
    [i for i, d in enumerate(delta) if d >= CUTOFF],
    key=lambda i: delta[i], reverse=True
)[:TOP_N_LABEL]

texts = []
for i in top_idx:
    t = axB.text(ptau[i], htau[i], variants[i],
                 fontsize=10, color='#111111', fontweight='bold', zorder=8)
    texts.append(t)

adjust_text(
    texts, ax=axB,
    arrowprops=dict(arrowstyle='-', color='#AAAAAA', lw=0.9,
                    shrinkA=6, shrinkB=4),
    expand=(1.5, 1.8),
    force_text=(0.9, 1.1),
    force_points=(0.4, 0.6),
    lim=700
)

# Colorbar
sm2 = ScalarMappable(cmap=cmap_sel, norm=norm)
sm2.set_array([])
cbar2 = figB.colorbar(sm2, ax=axB, pad=0.06, fraction=0.038, aspect=26)
cbar2.set_label('Δ SASA (Å²) — selected variants', fontsize=12.5, labelpad=10)
cbar2.ax.tick_params(labelsize=11)
cbar2.outline.set_linewidth(0.5)

# Axis formatting
axB.set_xlabel('pTau SASA (Å²)', fontsize=14, color='#111111', labelpad=10)
axB.set_ylabel('hTau SASA (Å²)', fontsize=14, color='#111111', labelpad=10)
axB.set_xlim(0, x_max)
axB.set_ylim(0, y_max)
axB.tick_params(axis='both', labelsize=12, length=4, colors='#333333')
axB.xaxis.grid(True, linestyle=':', linewidth=0.6, color='#DDDDDD', zorder=0)
axB.yaxis.grid(True, linestyle=':', linewidth=0.6, color='#DDDDDD', zorder=0)
axB.set_axisbelow(True)
for sp in ['top', 'right']:
    axB.spines[sp].set_visible(False)
for sp in ['left', 'bottom']:
    axB.spines[sp].set_color('#AAAAAA')
    axB.spines[sp].set_linewidth(0.9)

axB.set_title(
    'pTau vs hTau Structural Exposure Landscape\n'
    'Colored by Δ SASA magnitude',
    fontsize=16, fontweight='bold', color='#111111',
    pad=16, loc='left', linespacing=1.5
)

# Legend
sel_dot  = mlines.Line2D([], [], color=cmap_sel(0.75), marker='o',
                          linestyle='None', markersize=11,
                          label=f'Selected  (Δ SASA ≥ {CUTOFF} Å²)',
                          markeredgecolor='white', markeredgewidth=0.9)
rej_dot  = mlines.Line2D([], [], color=COL_REJ, marker='o',
                          linestyle='None', markersize=11,
                          label='Excluded', alpha=0.6,
                          markeredgecolor='white', markeredgewidth=0.9)
dup_dot  = mlines.Line2D([], [], color='#888888', marker='s',
                          linestyle='None', markersize=10,
                          label='Dual-site pTau variant',
                          markeredgecolor='white', markeredgewidth=0.9)
id_line  = mlines.Line2D([], [], color='#CCCCCC', lw=1.5, linestyle='--',
                          label='y = x  (identity)')
cut_line = mlines.Line2D([], [], color=COL_CUT, lw=2.0, linestyle='--',
                          label=f'Δ SASA = {CUTOFF} Å²  boundary')
zone_p   = mpatches.Patch(color='#AAAAAA', alpha=0.2, label='Excluded zone')

axB.legend(handles=[sel_dot, rej_dot, dup_dot, id_line, cut_line, zone_p],
           fontsize=11.5, frameon=True, edgecolor='#CCCCCC',
           facecolor='white', framealpha=0.97,
           loc='upper left',
           borderpad=0.85, labelspacing=0.6, handlelength=1.8)

figB.text(
    0.5, -0.01,
    f'n = {len(df)} pTau variants  |  '
    f'{n_unique_sel} variants selected (Δ SASA ≥ {CUTOFF} Å²)  |  '
    f'Square markers = dual-site pTau variant  |  '
    f'Top {TOP_N_LABEL} labelled  |  Color: {CMAP_NAME}',
    ha='center', fontsize=9.5, color='#888888', style='italic'
)

plt.tight_layout(pad=1.8)
figB.savefig('panelB_scatter.png', dpi=300, bbox_inches='tight', facecolor=BG)
plt.close(figB)
print('Panel B saved → panelB_scatter.png')

print(f'\nSummary:')
print(f'  Total pTau variants  : {len(df)}')
print(f'  Unique variant names : {df["pTau_Variant"].nunique()}')
print(f'  Selected sites       : {n_rows_sel}')
print(f'  Selected variants    : {n_unique_sel}')
print(f'  Cutoff               : {CUTOFF} Å²')
print(f'  x_max (Panel B)      : {x_max:.1f}')
print(f'  y_max (Panel B)      : {y_max:.1f}')