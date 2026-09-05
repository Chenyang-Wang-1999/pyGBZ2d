'''
Read-only diagnostics for the failed Hermitian-Haldane clustering.

Questions:
  1. exact core-core distance (the two band cores) + WHERE the closest
     approach happens (E, theta1, theta2);
  2. how far the fragments sit from the cores (the eps needed to absorb
     them -- must exceed the core-core bridge for a window to exist);
  3. exact-duplicate count (kNN p10 = 0 suggested coincident points);
  4. LineSubset arc-width distribution (tiny arcs -> single points after
     decimation -> suspected chaining gaps).
'''

import sys
import pickle
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
from scipy.spatial import cKDTree

from demo_band_clustering import (flatten_results, embed, grid_steps,
                                  radius_graph_labels)
from pygbz2d.core import LineSubset

with open("data/Hermitian-Haldane-amoeba.pkl", "rb") as fp:
    data = pickle.load(fp)
results = data["results"]

bp = flatten_results(results)
d_re, d_im = grid_steps(bp.E)
X = embed(bp, alpha_E=0.25, d_re=d_re, d_im=d_im)
print(f"points: {len(X)}")

# ---- 1. band cores at eps just below the merge ----------------------------
lab = radius_graph_labels(X, 0.3748)
sizes = np.bincount(lab)
order = np.argsort(-sizes)
coreA, coreB = int(order[0]), int(order[1])
mA, mB = lab == coreA, lab == coreB
print(f"\ncoreA size={sizes[coreA]}  coreB size={sizes[coreB]}  "
      f"fragment points={len(lab) - sizes[coreA] - sizes[coreB]}")
print(f"coreA ReE [{bp.E[mA].real.min():+.3f},{bp.E[mA].real.max():+.3f}]")
print(f"coreB ReE [{bp.E[mB].real.min():+.3f},{bp.E[mB].real.max():+.3f}]")

idxA, idxB = np.where(mA)[0], np.where(mB)[0]
dA, iA = cKDTree(X[idxA]).query(X[idxB], k=1)
jb = int(np.argmin(dA))
pa, pb = idxA[iA[jb]], idxB[jb]
print(f"\ncore-core min distance = {dA[jb]:.4f} at:")
for tag, p in (("A", pa), ("B", pb)):
    print(f"  {tag}: E={bp.E[p].real:+.4f}  th1={bp.theta1[p]:.4f}  "
          f"th2={bp.theta2[p]:.4f}  mu1={bp.mu1[p]:+.2e}  mu2={bp.mu2[p]:+.2e} "
          f" (slice {bp.slice_idx[p]}, subset {bp.subset_idx[p]})")

# ---- 2. fragments: distance to the NEAREST core ---------------------------
mF = ~(mA | mB)
if mF.any():
    tree_core = cKDTree(np.vstack([X[mA], X[mB]]))
    dF, _ = tree_core.query(X[mF], k=1)
    print(f"\nfragment->core distance percentiles: "
          f"p50={np.percentile(dF, 50):.3f} p90={np.percentile(dF, 90):.3f} "
          f"max={dF.max():.3f}")
    worst = np.argmax(dF)
    pw = np.where(mF)[0][worst]
    print(f"  farthest fragment point: E={bp.E[pw].real:+.4f} "
          f"th1={bp.theta1[pw]:.4f} th2={bp.theta2[pw]:.4f} "
          f"dist={dF[worst]:.3f}")

# ---- 3. exact duplicates ---------------------------------------------------
pairs0 = cKDTree(X).query_pairs(1e-12)
print(f"\nexact duplicate pairs: {len(pairs0)}")

# ---- 4. arc widths ---------------------------------------------------------
widths, Es, npts = [], [], []
for r in results:
    if not (r.success and r.subsets):
        continue
    for s in r.subsets:
        if isinstance(s, LineSubset):
            widths.append(s.theta1_width)
            Es.append(float(np.real(s.E)))
            npts.append(len(s.theta1_arr))
widths, Es = np.array(widths), np.array(Es)
print(f"\narcs: {len(widths)}, width percentiles "
      f"(5/25/50/75/95) = {np.percentile(widths, [5, 25, 50, 75, 95])}")
tiny = widths < 0.15
print(f"tiny arcs (<0.15 rad): {tiny.sum()}")
if tiny.any():
    print(f"  their E range: [{Es[tiny].min():+.3f}, {Es[tiny].max():+.3f}]")
    # per-E count of tiny arcs vs total
    for Ec in [1.0, 0.5, 0.0, -0.5, -3.0]:
        j = np.argmin(abs(Es - Ec))
        sel = abs(Es - Es[j]) < 1e-6
        print(f"  near E={Es[j]:+.3f}: {sel.sum()} arcs, "
              f"{tiny[sel].sum()} tiny, widths={np.sort(widths[sel])[:8]}")
