"""Compute local safe-radius pseudopotential barrier from the saved P2 field."""
import json,sys
from pathlib import Path
import numpy as np
from scipy.constants import elementary_charge,atomic_mass
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from iontrap_rng.fem_fields import P2ScalarField
p=ROOT/'results'/'m3a_deterministic_results.json';data=json.loads(p.read_text(encoding='utf-8'))
par=data['physical_parameters'];r0=par['r0_m'];omega=par['omega_rad_s'];Vrf=par['Vrf_differential_peak_V'];Vdc=par['Vdc_each_ring_V']
Q=elementary_charge;M=39.96259098*atomic_mass
field=ROOT/data['field_identity']['field_npz'];rf=P2ScalarField(field,'rf_potential');dc=P2ScalarField(field,'dc_potential')
def U(xi):
    g=rf.gradient(xi)/r0
    return Q**2*Vrf**2*np.sum(g*g,axis=0)/(4*M*omega**2)+Q*Vdc*dc.potential(xi)
U0=float(U(np.zeros((3,1)))[0]);radius=.2
mu=np.linspace(-1,1,31);phi=np.arange(64)*2*np.pi/64;dirs=[]
for z in mu:
    s=np.sqrt(max(0.,1-z*z))
    for a in phi:dirs.append([s*np.cos(a),s*np.sin(a),z])
dirs=np.asarray(dirs).T;vals=U(radius*dirs)-U0;i=int(np.argmin(vals))
cuts={}
coordinate=np.linspace(-radius,radius,161)
for axis,label in ((0,'radial_x'),(2,'axial_z')):
    x=np.zeros((3,len(coordinate)));x[axis]=coordinate
    cuts[label]={'xi':coordinate.tolist(),'delta_U_eV':((U(x)-U0)/Q).tolist()}
result={'formula':'U_eff=Q^2 |E_RF,peak|^2/(4m Omega^2)+Q Phi_DC',
 'center_U_J':U0,'sphere_radius_r0':radius,'min_delta_U_J':float(vals[i]),
 'min_delta_U_eV':float(vals[i]/Q),'minimum_direction':dirs[:,i].tolist(),
 'cuts':cuts,'caveat':'local safe-radius acceptance only; mesh and surface fields are not validated sufficiently to infer a physical escape saddle/depth'}
data['pseudopotential_local_final']=result;p.write_text(json.dumps(data,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in result.items() if k!='cuts'},indent=2))
