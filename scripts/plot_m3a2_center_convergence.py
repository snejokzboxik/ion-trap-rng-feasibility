"""Lightweight figures from the saved M3a.2 center-force comparison."""
import json
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
d=json.loads((ROOT/'results/m3a2_center_field_comparison.json').read_text())
colors={'hc018':'#4169a1','hc009':'#e07a24','hc006':'#20865a'}
labels={'hc018':'hcenter=0.18 r0','hc009':'hcenter=0.09 r0','hc006':'hcenter=0.06 r0'}
fig,axs=plt.subplots(1,3,figsize=(12,3.8),constrained_layout=True)
for ax,axis,component,title in zip(axs,'xyz',range(3),('Radial x','Radial y','Axial z')):
    for mesh in colors:
        rows=[]
        for row in d['force_point_comparisons']:
            name=row['point']
            if name.startswith(f'axis_{axis}_'):
                coord=float(name.rsplit('_',1)[1]); force=row['by_mesh'][mesh]['dc_force_N'][component]
                rows.append((coord,force))
        rows.sort()
        ax.plot([x[0] for x in rows],[x[1]*1e19 for x in rows],'-o',ms=3,label=labels[mesh],color=colors[mesh])
    ax.axhline(0,color='black',lw=.5);ax.set_xlabel(f'{axis}/r0');ax.set_ylabel('DC force component (10⁻¹⁹ N)');ax.set_title(title);ax.grid(alpha=.25)
axs[0].legend(fontsize=7)
out=ROOT/'results/figures/m3a2_center_force_convergence.png';fig.savefig(out,dpi=180);plt.close(fig)

fig,ax=plt.subplots(figsize=(7.5,4.2),constrained_layout=True)
for mesh in colors:
    rows=[]
    for row in d['force_point_comparisons']:
        if row['point'].startswith('axis_z_'):
            z=float(row['point'].rsplit('_',1)[1]);fz=row['by_mesh'][mesh]['dc_force_N'][2]
            rows.append((z,fz))
    rows.sort();ax.plot([a[0] for a in rows],[a[1]*1e19 for a in rows],'-o',label=labels[mesh],color=colors[mesh])
ax.set_xlabel('z/r0');ax.set_ylabel('Fz (10⁻¹⁹ N)');ax.set_title('Axial DC restoring force at fixed VDC');ax.grid(alpha=.25);ax.legend()
fig.savefig(ROOT/'results/figures/m3a2_axial_symmetry.png',dpi=180);plt.close(fig)
print(out)
