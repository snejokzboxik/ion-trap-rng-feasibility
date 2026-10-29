from pathlib import Path
import json
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
tags=['plateau030','plateau018','plateau014']; labels=['0.030','0.018','0.014']
med=[];p90=[];dof=[]
for tag in tags:
 d=json.loads((ROOT/'results'/f'm3a3_octant_mesh_{tag}.json').read_text())
 st=d['edge_lengths_r0']['tetrahedra_by_centroid_radius']['radius_lt_0.15']
 med.append(st['median']);p90.append(st['p90']);dof.append(d['estimated_P2_dofs'])
rec=json.loads((ROOT/'results/m3a3_octant_recovery.json').read_text())
mesh=rec['mesh_to_mesh']['dc']; generic=[r for r in mesh['pointwise'] if sum(abs(x)>1e-12 for x in r['point_r0'])>=2]
err=max(r['vector_change_pct'] for r in generic)
fig,ax=plt.subplots(1,2,figsize=(10,4.1),constrained_layout=True)
ax[0].plot(labels,med,'o-',label='median edge, tetrahedra centroid r<0.15r0')
ax[0].plot(labels,p90,'s--',label='90th percentile edge')
ax[0].set_xlabel('Requested hcenter / r0');ax[0].set_ylabel('Measured local edge length / r0');ax[0].grid(alpha=.3);ax[0].legend(fontsize=8)
for i,v in enumerate(dof):ax[0].annotate(f'{v:,} P2 DOF',(labels[i],med[i]),xytext=(4,7),textcoords='offset points',fontsize=8)
ax[1].bar(['P2 meshes\n0.030 → 0.018'],[err],color='#b44')
ax[1].axhline(1,color='black',ls=':',label='1% target');ax[1].set_ylim(0,max(3,err*1.35));ax[1].set_ylabel('Maximum DC vector-force change (%)');ax[1].grid(axis='y',alpha=.3);ax[1].legend(fontsize=8)
fig.suptitle('M3a.3: achieved central resolution and remaining DC force change')
fig.savefig(ROOT/'results/m3a3_convergence.png',dpi=180);plt.close(fig)
