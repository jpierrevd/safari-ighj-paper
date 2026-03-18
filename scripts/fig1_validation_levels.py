#!/usr/bin/env python3
"""Figure 1: Evidence levels by tribe — stacked bar chart."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

tribes = ['Bovini', 'Hippotragini', 'Cervidae', 'Antilopini', 'Reduncini', 'Alcelaphini', 'Tragelaphini']

# Evidence levels per tribe (number of species at each level)
level4 = [3, 0, 2, 0, 0, 0, 0]   # Bison×2+Yak | Cervus+Odocoileus
level3p = [1, 0, 0, 1, 1, 0, 0]  # Banteng | Saiga | Waterbuck
level3  = [0, 1, 0, 0, 0, 0, 0]  # Oryx (10 ind.)
level2p = [0, 2, 0, 0, 0, 1, 0]  # H.niger+Addax | Wildebeest
level2  = [0, 0, 1, 0, 0, 0, 0]  # Elaphurus
# Tragelaphini validated via SRR5647660 (Level 4 equivalent)
level4[6] = 1  # T. oryx expression confirmed

# Individuals per tribe (for annotation)
n_individuals = [18, 15, 15, 4, 3, 5, 1]
n_species = [4, 3, 3, 1, 1, 1, 1]

x = np.arange(len(tribes))
width = 0.6

fig, ax = plt.subplots(figsize=(12, 6))

colors = {'Level 4': '#1a5276', 'Level 3+': '#2e86c1', 'Level 3': '#5dade2',
          'Level 2+': '#85c1e9', 'Level 2': '#d4e6f1'}

bars4 = ax.bar(x, level4, width, label='Level 4 (Splice + Expression)', color=colors['Level 4'])
bars3p = ax.bar(x, level3p, width, bottom=level4, label='Level 3+ (Multi-individual junction)', color=colors['Level 3+'])
bottom3 = [a+b for a,b in zip(level4, level3p)]
bars3 = ax.bar(x, level3, width, bottom=bottom3, label='Level 3 (V(D)J junction)', color=colors['Level 3'])
bottom2p = [a+b for a,b in zip(bottom3, level3)]
bars2p = ax.bar(x, level2p, width, bottom=bottom2p, label='Level 2+ (Few-individual junction)', color=colors['Level 2+'])
bottom2 = [a+b for a,b in zip(bottom2p, level2p)]
bars2 = ax.bar(x, level2, width, bottom=bottom2, label='Level 2 (Mapping only)', color=colors['Level 2'])

# Annotations
for i, (ns, ni) in enumerate(zip(n_species, n_individuals)):
    total = level4[i] + level3p[i] + level3[i] + level2p[i] + level2[i]
    ax.text(i, total + 0.1, f'{ns} spp.\n{ni} ind.', ha='center', va='bottom', fontsize=8, fontweight='bold')

ax.set_ylabel('Number of species validated', fontsize=12)
ax.set_xlabel('Bovidae Tribe / Cervidae', fontsize=12)
ax.set_title('SAFARI-IGHJ Functional Validation: Evidence Levels Across Seven Tribes', fontsize=13, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(tribes, rotation=30, ha='right', fontsize=10)
ax.legend(loc='upper right', fontsize=9, framealpha=0.9)
ax.set_ylim(0, 6)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig('../figures/fig1_validation_levels.pdf', dpi=300, bbox_inches='tight')
plt.savefig('../figures/fig1_validation_levels.png', dpi=300, bbox_inches='tight')
print("Figure 1 saved.")
