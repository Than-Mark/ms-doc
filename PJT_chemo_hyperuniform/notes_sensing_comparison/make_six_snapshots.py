"""Assemble the six saved final configurations used in the structure-factor scan."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

here = Path(__file__).parent
data = here.parent/'PJT_quorum_sensing'/'results_wu'
configs = [(1.6, 41, 61), (2., 42, 62), (3., 43, 63)]
fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.65), constrained_layout=True)
for column, (u, seed40, seed80) in enumerate(configs):
    for row, (L, seed) in enumerate([(40, seed40), (80, seed80)]):
        z = np.load(data/f'L{L}_R5_dt0.05_seed{seed}_u{u:g}.npz')
        ax = axes[row, column]
        ax.scatter(z['positions'][:, 0], z['positions'][:, 1], s=.55 if L == 40 else .14,
                   c='black', linewidths=0, rasterized=True)
        ax.set(xlim=(0, L), ylim=(0, L), aspect='equal',
               title=fr'$\rho_0={u:g}$, $L={L}$, $N={len(z["positions"])}$')
        ax.set_xticks([])
        ax.set_yticks([])
fig.savefig(here/'six_final_snapshots.pdf')
fig.savefig(here/'six_final_snapshots.png', dpi=220)
