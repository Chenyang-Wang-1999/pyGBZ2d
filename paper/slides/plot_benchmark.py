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
date_tag = "20260921T144330717373Z"
for kind in ('amoeba','x-strip','y-strip','11-strip'):
    p=base/f'HN2D-{date_tag}-{kind}.pkl'
    with p.open('rb') as f: d=pickle.load(f)
    rs=d['results'][kind]
    assert all(r.success for r in rs)
    E=(d['real_axis'][None,:]+1j*d['imag_axis'][:,None]).ravel(order=d['flatten_order'])
    mask=np.array([r.is_gbz for r in rs])
    scans[kind]=(E,mask)
    provenance.append({'file':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'inside':int(mask.sum()),'total':len(mask)})
for kind in ('x-strip','y-strip'):
    assert np.array_equal(scans[kind][0],scans['amoeba'][0])
    assert np.array_equal(scans[kind][1],scans['amoeba'][1])
assert np.array_equal(scans['11-strip'][0],scans['amoeba'][0])
assert not np.any(scans['11-strip'][1]&~scans['amoeba'][1])
plt.rcParams.update({'font.family':'Arial','font.size':14,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.right':False,'axes.spines.top':False})
fig,axs=plt.subplots(1,2,figsize=(11.5,3.5),sharex=True,sharey=True,layout='constrained')
for ax,kind,title in zip(axs,('amoeba','11-strip'),('Amoeba, x strip, y strip','[11] strip')):
    E,mask=scans[kind]
    ax.scatter(E.real,E.imag,s=2,c='#e5e5e5',rasterized=False)
    ax.scatter(E[mask].real,E[mask].imag,s=3,c='#740c89' if kind=='amoeba' else '#167ba6',rasterized=False)
    ax.set(title=title,xlabel='Re E',ylabel='Im E',aspect='equal')
    ax.tick_params(length=3)
for ext in ('pdf','svg'):
    fig.savefig(f'paper/slides/Figures/hn-benchmark-spectra.{ext}',bbox_inches='tight')
plt.close(fig)
Path('paper/slides/.build/figure-provenance.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
print('Saved grid figure. Identical Cartesian masks and diagonal-strip inclusion verified.')
