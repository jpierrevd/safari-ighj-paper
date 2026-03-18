#!/usr/bin/env python3
"""Figure 2: Mapped vs Spliced reads for RNA-seq validated species."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

species = [
    'O. virginianus\n(deer, rpLN)',
    'B. mutus\n(yak, spleen)',
    'C. elaphus\n(red deer)',
    'B. bison\n(bison, direct)',
    'B. bonasus\n(wisent, cross-sp.)'
]

mapped =  [6953000, 1031357, 1627661, 469000, 327000]
spliced = [966000,  211019,  209678,  80000,  22000]
n_samples = [8, 3, 6, 5, 1]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={'width_ratios': [2, 1]})

# Panel A: Grouped bar chart
x = np.arange(len(species))
width = 0.35

bars1 = ax1.bar(x - width/2, [m/1e6 for m in mapped], width, label='Mapped reads',
                color='#2e86c1', edgecolor='white', linewidth=0.5)
bars2 = ax1.bar(x + width/2, [s/1e6 for s in spliced], width, label='Spliced reads (JH→CH1)',
                color='#e74c3c', edgecolor='white', linewidth=0.5)

# Sample count annotations
for i, n in enumerate(n_samples):
    total_h = max(mapped[i], spliced[i]) / 1e6
    ax1.text(i, total_h + 0.15, f'n={n}', ha='center', fontsize=9, fontstyle='italic')

ax1.set_ylabel('Reads (millions)', fontsize=12)
ax1.set_title('A. RNA-seq Validation: Mapped and Spliced Reads', fontsize=12, fontweight='bold')
ax1.set_xticks(x)
ax1.set_xticklabels(species, fontsize=9)
ax1.legend(fontsize=10)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

# Panel B: Splice ratio
splice_ratio = [s/m * 100 for s, m in zip(spliced, mapped)]
colors_ratio = ['#1a5276' if r > 10 else '#5dade2' for r in splice_ratio]

short_names = ['O. virg.', 'B. mutus', 'C. elaph.', 'B. bison', 'B. bonas.']
bars3 = ax2.barh(short_names, splice_ratio, color=colors_ratio, edgecolor='white', linewidth=0.5)

for i, (v, s, m) in enumerate(zip(splice_ratio, spliced, mapped)):
    ax2.text(v + 0.5, i, f'{v:.1f}%\n({s//1000}K/{m//1000}K)',
             va='center', fontsize=8)

ax2.set_xlabel('Splice ratio (%)', fontsize=12)
ax2.set_title('B. JH-CH1 Splice Ratio', fontsize=12, fontweight='bold')
ax2.set_xlim(0, 25)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.invert_yaxis()

plt.tight_layout()
plt.savefig('../figures/fig2_splicing_reads.pdf', dpi=300, bbox_inches='tight')
plt.savefig('../figures/fig2_splicing_reads.png', dpi=300, bbox_inches='tight')
print("Figure 2 saved.")
