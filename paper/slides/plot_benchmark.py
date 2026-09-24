from pathlib import Path
import pickle,sys,json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path('src').resolve()))
base=Path('application/data')
scans={}
provenance=[]
hoppings=[]
date_tag = "20260921T144330717373Z"
for kind in ('amoeba','x-strip','y-strip','11-strip'):
    p=base/f'HN2D-{date_tag}-{kind}.pkl'
    with p.open('rb') as f: d=pickle.load(f)
    rs=d['results'][kind]
    assert all(r.success for r in rs)
    E=(d['real_axis'][None,:]+1j*d['imag_axis'][:,None]).ravel(order=d['flatten_order'])
    mask=np.array([r.is_gbz for r in rs])
    scans[kind]=(E,mask)
    hoppings.append(d['hoppings'])
    provenance.append({'file':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'inside':int(mask.sum()),'total':len(mask)})
# The outlines below are analytic, so every panel must share the same model.
assert all(h==hoppings[0] for h in hoppings)
for kind in ('x-strip','y-strip'):
    assert np.array_equal(scans[kind][0],scans['amoeba'][0])
    assert np.array_equal(scans[kind][1],scans['amoeba'][1])
assert np.array_equal(scans['11-strip'][0],scans['amoeba'][0])
assert not np.any(scans['11-strip'][1]&~scans['amoeba'][1])

#### Closed-form GBZ outlines ####
# The amoebic, x-strip, and y-strip GBZs share the Cartesian spectrum
# {u*s+v*t : s,t in [-1,1]}, a parallelogram. The [11] strip instead has
# sigma = {E : E**2 in conv{0,(u+v)**2,(u-v)**2}}, so its outline is the
# square root of that triangle in the squared-energy plane. Both follow from
# Eqs. (S3.16) and (S3.30)-(S3.31) of arXiv:2506.22743v3.
def factorize_coefficients(J1,J2):
    """Return the benchmark's (gamma,delta,J) for one hopping pair."""
    gamma=0.5*np.log(abs(J1)/abs(J2))
    delta=0.5*np.angle(J1*J2)
    return gamma,delta,J1/np.exp(gamma+1j*delta)
def cartesian_halfwidths(model):
    """Return u=2*exp(1j*delta_x)*|Jx| and v=2*exp(1j*delta_y)*|Jy|."""
    _,delta_x,Jx=factorize_coefficients(model['Jx1'],model['Jx2'])
    _,delta_y,Jy=factorize_coefficients(model['Jy1'],model['Jy2'])
    return 2*np.exp(1j*delta_x)*abs(Jx),2*np.exp(1j*delta_y)*abs(Jy)
def parallelogram_border(u,v):
    """Closed outline of the shared Cartesian (amoeba/x-strip/y-strip) spectrum."""
    verts=np.array([u+v,-u+v,-u-v,u-v])
    return np.append(verts,verts[0])
def strip11_border(u,v,n=401):
    """Return one lobe of the [11]-strip outline; the other lobe is its negative.

    The outline of the squared-energy triangle conv{0,(u+v)**2,(u-v)**2} is
    traversed from its vertex at w=0, so both ends of the lifted curve are
    exactly E=0 and the lobe is closed without any branch bookkeeping at the
    branch point. The phase still has to be unwrapped: that outline crosses
    the negative real axis (here near w=-5.51), where the principal square
    root would otherwise jump by a sign. Each lobe is the square root of one
    traversal of the triangle's boundary, the two lobes meeting only at E=0.
    """
    w1,w2=(u+v)**2,(u-v)**2
    t=np.linspace(0.0,1.0,n)
    path=np.concatenate([w1*t,w1+(w2-w1)*t[1:],w2*(1-t)[1:]])
    return np.sqrt(np.abs(path))*np.exp(0.5j*np.unwrap(np.angle(path)))

u,v=cartesian_halfwidths(hoppings[0])
plt.rcParams.update({'font.family':'Arial','font.size':14,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.right':False,'axes.spines.top':False})
fig,axs=plt.subplots(1,2,figsize=(11.5,3.5),sharex=True,sharey=True,layout='constrained')
for ax,kind,title in zip(axs,('amoeba','11-strip'),('Amoeba, x strip, y strip','[11] strip')):
    E,mask=scans[kind]
    ax.scatter(E.real,E.imag,s=2,c='#e5e5e5',rasterized=False)
    ax.scatter(E[mask].real,E[mask].imag,s=3,c='#740c89' if kind=='amoeba' else '#167ba6',rasterized=False)
    if kind=='amoeba':
        border=parallelogram_border(u,v)
        ax.plot(border.real,border.imag, 'r--',lw=1.2)
    else:
        lobe=strip11_border(u,v)
        ax.plot(lobe.real,lobe.imag,'r--',lw=1.2)
        ax.plot(-lobe.real,-lobe.imag,'r--',lw=1.2)
    ax.set(title=title,xlabel='Re E',ylabel='Im E',aspect='equal')
    ax.tick_params(length=3)
for ext in ('pdf','svg'):
    fig.savefig(f'paper/slides/Figures/hn-benchmark-spectra.{ext}',bbox_inches='tight')
plt.close(fig)
print('Saved grid figure with closed-form outlines. Identical Cartesian masks and diagonal-strip inclusion verified.')
