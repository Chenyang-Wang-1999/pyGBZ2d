"""Standard non-Hermitian QWZ benchmark using the existing pyGBZ2d pipeline.

Default data directory: application/data/QWZ-benchmark.
The computation stages are coarse (20 x 20), fine (201 x 201), and mesh/Chern.
GBZ solving, mesh construction/refinement, and Chern integration use the installed
pygbz2d package. The chern stage reintegrates saved meshes using
reciprocal biorthogonal links, without rerunning sweeps or mesh refinement.
This script can be copied outside the checkout; no sibling Python files are
required. Install pyGBZ2d and the example's third-party dependencies first.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
from dataclasses import asdict
from datetime import datetime, timezone
import json
import multiprocessing as mp
import os
from pathlib import Path
import pickle
import time

# Worker-level root solves should not each start another BLAS thread pool.
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
import numpy as np

APPLICATION_DIR = Path(__file__).resolve().parent
DATA = APPLICATION_DIR / "data/QWZ-benchmark"
CASES = {"nontrivial": 1.2, "trivial": 2.6}
GAMMA = .4
SX = np.array([[0, 1], [1, 0]], complex)
SY = np.array([[0, -1j], [1j, 0]], complex)
SZ = np.diag([1., -1.])


def hamiltonian(beta1, beta2, mass=1.2, gamma=GAMMA):
    sx, sy = (beta1-1/beta1)/(2j), (beta2-1/beta2)/(2j)
    cx, cy = (beta1+1/beta1)/2, (beta2+1/beta2)/2
    return (sx+1j*gamma)*SX + (sy+1j*gamma)*SY + (mass-cx-cy)*SZ


def make_model():
    """Default QWZ Hamiltonian factory (m=1.2)."""
    return lambda b1, b2: hamiltonian(b1, b2, 1.2)


def make_trivial_model():
    return lambda b1, b2: hamiltonian(b1, b2, 2.6)


def build_model(mass=1.2, gamma=GAMMA):
    from BerryPy import TightBinding as tb
    onsite = mass*SZ + 1j*gamma*(SX+SY)
    intra = [[source, dest, onsite[dest, source]] for source in range(2) for dest in range(2)
             if onsite[dest, source] != 0]
    inter = []
    for direction, sigma in [((1, 0), SX), ((0, 1), SY)]:
        # BerryPy multiplies a shift-R bond by beta**(-R).
        positive_power = -.5j*sigma - .5*SZ
        for shift, hopping in [(tuple(-np.array(direction)), positive_power),
                               (direction, positive_power.conj().T)]:
            inter.extend([[source, dest, hopping[dest, source], shift]
                          for source in range(2) for dest in range(2) if hopping[dest, source] != 0])
    return tb.TightBindingModel(2, 2, np.eye(2), intra, inter)


def polynomial(mass, gamma=GAMMA):
    coeffs, degs = build_model(mass, gamma).get_characteristic_polynomial_data()
    coeffs, degs = np.asarray(coeffs, complex), np.asarray(degs, int)
    rng = np.random.default_rng(20260925)
    for _ in range(8):
        e = rng.normal()+1j*rng.normal()
        b = np.exp(rng.normal(0, .2, 2)+1j*rng.uniform(-np.pi, np.pi, 2))
        expected = np.linalg.det(e*np.eye(2)-hamiltonian(*b, mass, gamma))
        actual = np.sum(coeffs*np.prod(np.array([e, *b])**degs, axis=1))
        np.testing.assert_allclose(actual, expected, rtol=1e-11, atol=1e-11)
    return coeffs, degs


def save_pickle(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix+".tmp")
    with temp.open("wb") as stream:
        pickle.dump(payload, stream, protocol=4)
    os.replace(temp, path)


def read_pickle(path):
    with Path(path).open("rb") as stream:
        return pickle.load(stream)


_COEFFS = _DEGS = None


def init_worker(coeffs, degs):
    global _COEFFS, _DEGS
    _COEFFS, _DEGS = coeffs, degs


def solve_one(task):
    from pygbz2d.amoeba import collect_GBZ_subsets
    index, energy = task
    started = time.monotonic()
    result = collect_GBZ_subsets(_COEFFS, _DEGS, complex(energy))
    return index, result, time.monotonic()-started


def counts(data):
    results = data["results"]
    return dict(total=len(results), done=sum(r is not None for r in results),
                in_spectrum=sum(r is not None and r.is_gbz for r in results),
                failed=sum(r is not None and not r.success for r in results))


def boundary_count(data):
    grid=np.array(data["results"],object).reshape(len(data["E_imag"]),len(data["E_real"]))
    rim=np.concatenate([grid[0],grid[-1],grid[1:-1,0],grid[1:-1,-1]])
    return sum(r is not None and r.is_gbz for r in rim)


def scan(mass, real, imag, path, workers=20):
    coeffs, degs = polynomial(mass)
    energy = (real[None, :]+1j*imag[:, None]).ravel()
    if path.exists():
        data = read_pickle(path)
        if not (np.array_equal(data["E_real"], real) and np.array_equal(data["E_imag"], imag)
                and data["mass"] == mass and data["gamma"] == GAMMA
                and np.array_equal(data["degs"], degs) and np.allclose(data["coeffs"], coeffs)):
            raise ValueError(f"Incompatible checkpoint: {path}")
    else:
        data = dict(mass=mass, gamma=GAMMA, coeffs=coeffs, degs=degs,
                    E_real=real, E_imag=imag, results=[None]*len(energy),
                    seconds=np.zeros(len(energy)), method="amoeba", workers=workers,
                    created_utc=datetime.now(timezone.utc).isoformat(), elapsed_seconds=0.)
    pending = [(i, e) for i, e in enumerate(energy) if data["results"][i] is None]
    print(f"{path.name}: {len(pending)}/{len(energy)} pending; {workers} workers", flush=True)
    started = last_save = time.monotonic()
    saved_elapsed = data["elapsed_seconds"]
    save_pickle(path, data)
    try:
        with mp.Pool(workers, initializer=init_worker, initargs=(coeffs, degs)) as pool:
            for done, (index, result, seconds) in enumerate(pool.imap_unordered(solve_one, pending, chunksize=1), 1):
                data["results"][index], data["seconds"][index] = result, seconds
                now = time.monotonic()
                if not result.success:
                    print(f"FAILED {path.name} E={result.E_ref}: {result.error}", flush=True)
                if done % 100 == 0 or now-last_save >= 30 or done == len(pending):
                    data["elapsed_seconds"] = saved_elapsed+now-started
                    save_pickle(path, data)
                    summary = counts(data)
                    rate = done/max(now-started, 1e-9)
                    print(f"{path.name}: {summary}; elapsed={(now-started)/60:.1f} min; "
                          f"ETA={(len(pending)-done)/rate/60:.1f} min", flush=True)
                    last_save = now
    finally:
        data["elapsed_seconds"] = saved_elapsed+time.monotonic()-started
        save_pickle(path, data)
    return data


def coarse(case, directory, workers):
    mass = CASES[case]
    real = np.linspace(-mass-2.2, mass+2.2, 20)
    # An even grid is offset to include the real axis, where thin spectra and
    # continuum subsets would otherwise be missed. Both physical Im bounds
    # +/-sqrt(2)*gamma are still strictly inside this window.
    im_bound = np.sqrt(2)*GAMMA+.15
    imag = np.linspace(-im_bound, im_bound, 20, endpoint=False)
    return scan(mass, real, imag, directory/f"QWZ-{case}-coarse.pkl", workers)


def fine_window(data):
    good = np.array([r.E_ref for r in data["results"] if r is not None and r.is_gbz])
    if not len(good):
        raise RuntimeError("The coarse sweep found no spectrum; inspect it before a fine sweep")
    dr, di = np.diff(data["E_real"])[0], np.diff(data["E_imag"])[0]
    # E -> -E and conjugation are symmetries of the characteristic polynomial.
    # Two coarse cells of margin protect thin spectral tips between grid nodes.
    re_bound = max(abs(good.real))+2*dr
    # The coarse mesh can miss thin complex-energy lobes. A preliminary fine
    # run exposed such a lobe at the window edge. Here H-H^dagger consists only
    # of the onsite gamma*(sigma_x+sigma_y), so +/-sqrt(2)*gamma bounds Im(E)
    # for every finite OBC sample. Keep that whole enclosure plus one coarse
    # cell; the final boundary check below remains a numerical coverage check.
    im_bound = max(max(abs(good.imag))+2*di,np.sqrt(2)*GAMMA+di)
    return np.linspace(-re_bound, re_bound, 201), np.linspace(-im_bound, im_bound, 201)


def fine(case, directory, workers):
    data = read_pickle(directory/f"QWZ-{case}-coarse.pkl")
    real, imag = fine_window(data)
    print(f"{case}: fine bounds Re=[{real[0]},{real[-1]}], Im=[{imag[0]},{imag[-1]}]", flush=True)
    result=scan(CASES[case], real, imag, directory/f"QWZ-{case}-fine.pkl", workers)
    if boundary_count(result):
        raise RuntimeError("Spectrum touches the fine-window boundary; expand the window before constructing a mesh")
    return result


def calculate_chern(mesh_file, Hfun, *, save_flux=False, **options):
    """Example-owned NPZ I/O around the array-based library integrator."""
    from pygbz2d.experimental import integrate_chern
    mesh_file = Path(mesh_file)
    with np.load(mesh_file) as data:
        result, flux = integrate_chern(
            data["verts"], data["triangles"], data["E"], data["mu1"], data["mu2"],
            Hfun, **options)
    if save_flux:
        out = mesh_file.with_name(mesh_file.stem + "_chern_flux.npz")
        np.savez(out, **flux)
        print(f"    flux saved -> {out}")
    return result


def assess_mesh(case, band, level, mesh_path, directory):
    from pygbz2d.experimental.torus_mesh import unwrap_triangle, mesh_topology, edge_lengths
    with np.load(mesh_path) as mesh:
        topology = mesh_topology(mesh["triangles"])
        lengths = edge_lengths(mesh["verts"], mesh["triangles"])[1]
        triangles = mesh["triangles"]
        topology["unused_vertices"] = int(len(mesh["verts"])-len(np.unique(triangles)))
        edges = np.sort(np.vstack([triangles[:,[0,1]],triangles[:,[1,2]],triangles[:,[2,0]]]),axis=1)
        _, incidence = np.unique(edges,axis=0,return_counts=True)
        topology["nonmanifold_edges"] = int(np.sum(incidence != 2))
        p=unwrap_triangle(mesh["verts"],triangles)
        area=.5*abs((p[:,1,0]-p[:,0,0])*(p[:,2,1]-p[:,0,1])
                     -(p[:,1,1]-p[:,0,1])*(p[:,2,0]-p[:,0,0]))
        topology["phase_area_over_4pi2"] = float(area.sum()/(4*np.pi**2))
    rows=[]
    for mode in ["right", "biorthogonal"]:
        result = calculate_chern(mesh_path, lambda b1,b2: hamiltonian(b1,b2,CASES[case]),
                                 mode=mode, save_flux=True)
        flux_path=mesh_path.with_name(mesh_path.stem+"_chern_flux.npz")
        named_flux=mesh_path.with_name(mesh_path.stem+f"_{mode}_flux.npz")
        os.replace(flux_path,named_flux)
        row=dict(case=case,mass=CASES[case],gamma=GAMMA,band=band,level=level,mode=mode,mesh=str(mesh_path),
                 **asdict(result),topology=topology,max_edge=float(lengths.max()),
                 p99_edge=float(np.quantile(lengths,.99)),flux_file=str(named_flux),
                 literature_sign_chern=-result.chern)
        rows.append(row)
        print(json.dumps(row),flush=True)
    return rows


def mesh_and_chern(case, directory, workers, edge_thresh=.2, max_iters=10):
    from pygbz2d.experimental import flatten_results, build_band_mesh, refine_mesh
    path = directory/f"QWZ-{case}-fine.pkl"
    data = read_pickle(path)
    if data["mass"]!=CASES[case] or data["gamma"]!=GAMMA:
        raise ValueError("The sweep parameters do not match the Chern Hamiltonian")
    if np.any(np.asarray(data["degs"])[:,0] % 2):
        raise ValueError("Partner-mesh reuse requires a characteristic polynomial even in E")
    if counts(data)["done"] != len(data["results"]):
        raise RuntimeError("Fine sweep is unfinished")
    if boundary_count(data):
        raise RuntimeError("Spectrum touches the energy-window boundary; inspect coverage first")
    # Two separated half planes identify the two bands without relying on a
    # clustering radius to bridge across their spectral gap.
    summaries = []
    for band, sign in [("lower", -1), ("upper", 1)]:
        selected = [r for r in data["results"] if r.success and r.is_gbz and sign*r.E_ref.real > 0]
        points = flatten_results(selected)
        if points.n_dropped:
            raise RuntimeError(f"{case}/{band}: {points.n_dropped} nonfinite subset samples")
        verts,tri,vdata,_=build_band_mesh(points)
        base_path=directory/f"QWZ-{case}-{band}-base-mesh.npz"
        np.savez_compressed(base_path,verts=verts,triangles=tri,
                            E=vdata["E"],mu1=vdata["mu1"],mu2=vdata["mu2"])
        summaries.extend(assess_mesh(case,band,"before refinement",base_path,directory))
        (directory/f"QWZ-{case}-chern.json").write_text(json.dumps(summaries,indent=2),encoding="utf-8")
        mesh_path = directory/f"QWZ-{case}-{band}-mesh.npz"
        if band=="upper":
            # This trace-zero model has f(E,beta)=f(-E,beta). Its Ronkin
            # minimizers and GBZ are identical for the partner energies, so
            # the geometric refinement is shared. Eigenvectors and Wilson
            # phases of the upper band are still calculated independently.
            lower_path=directory/f"QWZ-{case}-lower-mesh.npz"
            with np.load(lower_path) as lower:
                verts,tri=lower["verts"],lower["triangles"]
                vdata={"E":-lower["E"],"mu1":lower["mu1"],"mu2":lower["mu2"]}
            origin="partner band: exact E -> -E symmetry"
            print(f"{case}/upper: reuse the refined lower-band torus by exact spectral symmetry",flush=True)
        elif mesh_path.exists():
            with np.load(mesh_path) as cached:
                if "edge_target" in cached and float(cached["edge_target"])!=edge_thresh:
                    raise ValueError("Saved mesh has a different refinement target; use a fresh data directory")
                verts,tri=cached["verts"],cached["triangles"]
                vdata={k:cached[k] for k in ["E","mu1","mu2"]}
            origin="cached adaptive refinement"
            print(f"{case}/lower: reuse saved adaptive-refinement output {mesh_path}",flush=True)
        else:
            indices=np.asarray([selected[s].index for s in vdata["slice_idx"]],dtype=int)
            verts, tri, vdata = refine_mesh(
                verts, tri, vdata, data["coeffs"], data["degs"], "amoeba", indices,
                edge_thresh=edge_thresh, max_iters=max_iters, n_procs=workers, verbose=True)
            origin="adaptive refinement from fine-sweep samples"
        np.savez_compressed(mesh_path, verts=verts, triangles=tri,
                            E=vdata["E"], mu1=vdata["mu1"], mu2=vdata["mu2"],
                            mass=CASES[case],gamma=GAMMA,edge_target=edge_thresh,mesh_origin=origin)
        summaries.extend(assess_mesh(case,band,"after refinement",mesh_path,directory))
        (directory/f"QWZ-{case}-chern.json").write_text(json.dumps(summaries,indent=2),encoding="utf-8")
    return summaries


def bz_checks(directory):
    from pygbz2d.experimental import regular_torus_mesh
    rows=[]
    for case,mass in CASES.items():
        for gamma in [0.,GAMMA]:
            for n in [21,41,81]:
                verts,tri=regular_torus_mesh(n)
                spectra=np.array([np.linalg.eigvals(hamiltonian(*np.exp(1j*v),mass,gamma)) for v in verts])
                indices=np.argmin(spectra.real,axis=1)
                energies=spectra[np.arange(len(verts)),indices]
                path=directory/f"QWZ-{case}-BZ-g{gamma:g}-n{n}.npz"
                np.savez_compressed(path,verts=verts,triangles=tri,E=energies,
                                    mu1=np.zeros(len(verts)),mu2=np.zeros(len(verts)))
                for mode in ["right","biorthogonal"]:
                    result=calculate_chern(path,lambda b1,b2:hamiltonian(b1,b2,mass,gamma),mode=mode)
                    row=dict(case=case,mass=mass,gamma=gamma,n=n,mode=mode,mesh=str(path),
                             literature_sign_chern=-result.chern,min_abs_real_energy=float(np.min(abs(energies.real))),
                             **asdict(result))
                    rows.append(row)
        print(f"{case}: BZ and Hermitian checks finished",flush=True)
    (directory/"QWZ-BZ-checks.json").write_text(json.dumps(rows,indent=2),encoding="utf-8")
    return rows


def chern_from_saved_meshes(case, directory):
    """Reintegrate saved inputs without rewriting or refining any mesh."""
    rows = []
    for band in ["lower", "upper"]:
        for suffix, level in [("base-mesh", "before refinement"), ("mesh", "after refinement")]:
            path = directory/f"QWZ-{case}-{band}-{suffix}.npz"
            with np.load(path) as data:
                for key, value in [("mass", CASES[case]), ("gamma", GAMMA)]:
                    if key in data and float(data[key]) != value:
                        raise ValueError(f"{path}: saved {key} does not match the Hamiltonian")
            rows.extend(assess_mesh(case, band, level, path, directory))
    (directory/f"QWZ-{case}-chern.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")
    return rows


def diagnose_meshes(directory):
    """Diagnose saved meshes and links without running any GBZ root solves."""
    from pygbz2d.experimental.chern import (left_right_eigenvectors_on_mesh,
                              triangle_flux_biorthogonal_complex,
                              biorthogonal_edge_links)
    from pygbz2d.experimental.torus_mesh import edge_lengths, orient_triangles
    for case,mass in CASES.items():
        path=directory/f"QWZ-{case}-lower-mesh.npz"
        if not path.exists(): continue
        records=[]
        h=lambda b1,b2:hamiltonian(b1,b2,mass)
        for mode in ["right","biorthogonal"]:
            for orientation in ["positive","negative"]:
                result=calculate_chern(path,h,mode=mode,orientation=orientation)
                records.append(dict(mode=mode,orientation=orientation,raw_chern=result.chern,
                                    link_rule=result.link_rule,total_flux_imag=result.total_flux_imag))
        with np.load(path) as d:
            b1=np.exp(d["mu1"]+1j*d["verts"][:,0])
            b2=np.exp(d["mu2"]+1j*d["verts"][:,1])
            right,left,_=left_right_eigenvectors_on_mesh(h,d["E"],b1,b2,1e-6)
            repeated_edges,lengths=edge_lengths(d["verts"],d["triangles"])
            edges,first=np.unique(np.sort(repeated_edges,axis=1),axis=0,return_index=True)
            lengths=lengths[first]
            forward=np.einsum('ij,ij->i',left[edges[:,0]].conj(),right[edges[:,1]])
            backward=np.einsum('ij,ij->i',left[edges[:,1]].conj(),right[edges[:,0]])
            bias=float(-np.angle(forward*backward).sum()/(2*np.pi))
            tri,_,_=orient_triangles(d["verts"],d["triangles"])
            flux=triangle_flux_biorthogonal_complex(right,left,tri)
            reversed_flux=triangle_flux_biorthogonal_complex(right,left,tri[:,::-1])
            rng=np.random.default_rng(20260926)
            gauge=np.exp(rng.uniform(-12,12,len(right))+1j*rng.uniform(-np.pi,np.pi,len(right)))
            changed=triangle_flux_biorthogonal_complex(
                right*gauge[:,None],left/gauge.conj()[:,None],tri)
            links=biorthogonal_edge_links(right,left,edges)
            reverse_links=biorthogonal_edge_links(right,left,edges[:,::-1])
            link_checks=dict(
                max_local_orientation_error=float(abs(flux+reversed_flux).max()),
                max_local_complex_gauge_error=float(abs(flux-changed).max()),
                max_edge_reciprocity_error=float(abs(links*reverse_links-1).max()),
                min_overlap_product_real=float((forward*backward).real.min()),
                complex_chern_imag=float(flux.imag.sum()/(2*np.pi)))
            long_edges=[]
            for index in np.flatnonzero(lengths>.2):
                a,b=edges[index]
                long_edges.append(dict(length=float(lengths[index]),theta_a=d["verts"][a].tolist(),
                    theta_b=d["verts"][b].tolist(),E_a=[float(d['E'][a].real),float(d['E'][a].imag)],
                    E_b=[float(d['E'][b].real),float(d['E'][b].imag)]))
        result=dict(case=case,orientation_checks=records,
                    legacy_LR_shared_edge_phase_bias_raw=bias,reciprocal_link_checks=link_checks,
                    unique_long_edge_count=len(long_edges),long_edges=long_edges)
        (directory/f"QWZ-{case}-diagnostics.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
        print(f"{case}: legacy raw-overlap bias={bias:.12g}; corrected checks={link_checks}; "
              f"{len(long_edges)} unique edges above 0.2 rad",flush=True)


def seam_display_mask(verts, triangles):
    """Hide seam-crossing faces only in the flat figure, never in integrals."""
    corners = verts[triangles]
    return np.max(abs(corners[:, [1, 2, 0]]-corners[:, [0, 1, 2]]), axis=(1, 2)) <= np.pi


def source_provenance():
    """Record the running script/package sources without assuming a checkout."""
    sources = [("application", Path(__file__))]
    for name in ["core", "amoeba.amoeba", "amoeba.bisect", "experimental.band_clustering",
                 "experimental.torus_mesh", "experimental.mesh_refinement", "experimental.chern"]:
        module = importlib.import_module("pygbz2d."+name)
        location = getattr(module, "__file__", None)
        if location is not None:
            sources.append((module.__name__, Path(location)))
    # A source-less installed package can still run the example; hashes are
    # provenance metadata, not a requirement for calculation or reporting.
    return [dict(module=name, path=str(path.resolve()),
                 sha256=hashlib.sha256(path.read_bytes()).hexdigest())
            for name, path in sources if path.is_file()]


def report(directory):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from pygbz2d.experimental import flatten_results
    from pygbz2d.experimental.torus_mesh import edge_lengths
    from matplotlib.collections import LineCollection
    plt.rcParams.update({"font.size":11,"svg.fonttype":"none","savefig.dpi":180})
    out=APPLICATION_DIR/"Figures/Chern-number-calculation"
    out.mkdir(parents=True,exist_ok=True)
    summary={}
    fig,axes=plt.subplots(2,2,figsize=(9,6.8),layout="constrained")
    for col,(case,mass) in enumerate(CASES.items()):
        summary[case]={"mass":mass,"gamma":GAMMA}
        for row,stage in enumerate(["coarse","fine"]):
            path=directory/f"QWZ-{case}-{stage}.pkl"
            if not path.exists(): continue
            data=read_pickle(path)
            result=data["results"]
            inside=np.array([r.E_ref for r in result if r is not None and r.is_gbz],complex)
            failed=np.array([r.E_ref for r in result if r is not None and not r.success],complex)
            outside=np.array([r.E_ref for r in result if r is not None and r.success and not r.is_gbz],complex)
            summary[case][stage]=dict(**counts(data),elapsed_seconds=data["elapsed_seconds"],
                                     boundary_in_spectrum=boundary_count(data),
                                     real_range=[float(data["E_real"][0]),float(data["E_real"][-1])],
                                     imag_range=[float(data["E_imag"][0]),float(data["E_imag"][-1])])
            if stage=="fine" and all(r is not None for r in result):
                bp=flatten_results(result)
                variables=np.stack([bp.E,np.exp(bp.mu1+1j*bp.theta1),np.exp(bp.mu2+1j*bp.theta2)],axis=1)
                terms=data["coeffs"][None,:]*np.prod(variables[:,None,:]**data["degs"][None,:,:],axis=2)
                relative_residual=abs(terms.sum(axis=1))/np.maximum(1,abs(terms).sum(axis=1))
                summary[case][stage].update(samples=len(bp.E),nonfinite_samples=bp.n_dropped,
                    max_polynomial_residual=float(relative_residual.max()) if len(bp.E) else None,
                    max_mu_symmetry_deviation=float(np.max(abs(bp.mu1-bp.mu2))) if len(bp.E) else None)
            ax=axes[row,col]
            ax.scatter(outside.real,outside.imag,s=.6 if stage=="fine" else 4,color="#c4ccd2",
                       alpha=.35 if stage=="fine" else .6,rasterized=True,label="Outside")
            ax.scatter(inside.real,inside.imag,s=3,color="#315c99" if case=="nontrivial" else "#be5861",
                       rasterized=True,label="Amoeba spectrum")
            if len(failed): ax.scatter(failed.real,failed.imag,s=15,marker="x",color="black",label="Failed")
            ax.set(title=f"m={mass}, {stage}: {counts(data)['done']}/{len(result)}",xlabel=r"Re $E$",ylabel=r"Im $E$")
            ax.set_xlim(data["E_real"][0],data["E_real"][-1])
            ax.set_ylim(data["E_imag"][0],data["E_imag"][-1])
            if stage=="fine" and len(inside) and all(r is not None for r in result) and not len(failed):
                # Show the thin complex spectrum clearly; the full computed
                # window remains in the numerical table and boundary checks.
                im_zoom=max(1.25*float(np.max(abs(inside.imag))),3*float(np.diff(data["E_imag"])[0]))
                ax.set_ylim(-im_zoom,im_zoom)
                ax.set_title(f"m={mass}, fine: {counts(data)['done']}/{len(result)} (zoom)")
            ax.legend(loc="upper right",fontsize=9,frameon=False,markerscale=2)
    for extension in ["svg","png"]:fig.savefig(out/f"sweeps.{extension}")
    plt.close(fig)
    rows=[]
    for case in CASES:
        path=directory/f"QWZ-{case}-chern.json"
        if path.exists(): rows.extend(json.loads(path.read_text(encoding="utf-8")))
    summary["chern"]=rows
    refined=[r for r in rows if r["level"]=="after refinement"]
    unique_meshes={(r["case"],r["band"]) for r in refined}
    summary["refined_band_meshes"]=len(unique_meshes)
    summary["expected_band_meshes"]=2*len(CASES)
    summary["all_refined_meshes_closed"]=(len(unique_meshes)==2*len(CASES) and all(
        r["topology"]["chi"]==0 and r["topology"]["n_components"]==1
        and r["topology"]["nonmanifold_edges"]==0 and r["topology"].get("unused_vertices",0)==0
        and abs(r["topology"].get("phase_area_over_4pi2",0)-1)<1e-8 for r in refined))
    summary["all_refined_edges_below_target"]=(len(unique_meshes)==2*len(CASES) and all(r["max_edge"]<=.2 for r in refined))
    summary["diagnostics"]={}
    for case in CASES:
        path=directory/f"QWZ-{case}-diagnostics.json"
        if path.exists():summary["diagnostics"][case]=json.loads(path.read_text(encoding="utf-8"))
    bz=directory/"QWZ-BZ-checks.json"
    if bz.exists():summary["bz_checks"]=json.loads(bz.read_text(encoding="utf-8"))
    summary["report_source_provenance"]=source_provenance()
    baseline=directory/"biorthogonal-raw-overlap-baseline"
    if baseline.exists():
        summary["legacy_biorthogonal_baseline"]=str(baseline)
        manifest=baseline/"unchanged-input-mesh-hashes.json"
        if manifest.exists():
            original=json.loads(manifest.read_text(encoding="utf-8-sig"))
            summary["integration_input_meshes_unchanged"]=all(
                hashlib.sha256((directory/name).read_bytes()).hexdigest()==digest
                for name,digest in original.items())
    (out/"benchmark-summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    if not rows:return summary
    fig,axes=plt.subplots(2,2,figsize=(9,7),layout="constrained")
    flux_limit=max([r["flux_max_abs"] for r in rows if r["level"]=="after refinement"
                    and r["band"]=="lower" and r["mode"]=="right"]+[1e-12])
    for col,(case,mass) in enumerate(CASES.items()):
        path=directory/f"QWZ-{case}-lower-mesh.npz"
        if not path.exists():continue
        with np.load(path) as d:
            verts,tri=d["verts"],d["triangles"]
            display=seam_display_mask(verts,tri)
            reduced=tri[display]
            axes[0,col].triplot(verts[:,0],verts[:,1],reduced,color="#9ca6b3",lw=.2,rasterized=True)
            sc=axes[0,col].scatter(verts[:,0],verts[:,1],c=d["mu1"],s=1,cmap="viridis",rasterized=True)
            edges,lengths=edge_lengths(verts,tri)
            visible=(lengths>.2)&(np.max(abs(verts[edges[:,0]]-verts[edges[:,1]]),axis=1)<=np.pi)
            axes[0,col].add_collection(LineCollection(verts[edges[visible]],colors="#da3548",linewidths=1.0,zorder=4))
            fig.colorbar(sc,ax=axes[0,col],label=r"$\mu_x$")
            flux_path=path.with_name(path.stem+"_right_flux.npz")
            with np.load(flux_path) as f:
                phase=-f["flux"]
            sc=axes[1,col].tripcolor(verts[:,0],verts[:,1],reduced,
                                    facecolors=phase[display],cmap="coolwarm",vmin=-flux_limit,vmax=flux_limit,rasterized=True)
            fig.colorbar(sc,ax=axes[1,col],label="Triangle Berry flux")
        for row in range(2):axes[row,col].set(xlabel=r"$\theta_x$",ylabel=r"$\theta_y$",title=f"m={mass}, lower band",
                                             xlim=(0,2*np.pi),ylim=(0,2*np.pi))
    for extension in ["svg","png"]:fig.savefig(out/f"gbz-and-flux.{extension}")
    plt.close(fig)
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["coarse", "fine", "mesh", "chern", "bz-check", "diagnose", "report", "all"])
    parser.add_argument("--cases", nargs="+", choices=list(CASES), default=list(CASES))
    parser.add_argument("--data-dir", type=Path, default=DATA)
    parser.add_argument("--workers", type=int, default=20)
    parser.add_argument("--edge-thresh", type=float, default=.2)
    parser.add_argument("--max-iters", type=int, default=10)
    args=parser.parse_args()
    args.data_dir.mkdir(parents=True,exist_ok=True)
    if args.stage=="bz-check":
        bz_checks(args.data_dir)
        return
    if args.stage=="report":
        report(args.data_dir)
        return
    if args.stage=="diagnose":
        diagnose_meshes(args.data_dir)
        return
    for stage in (["coarse", "fine", "mesh"] if args.stage=="all" else [args.stage]):
        for case in args.cases:
            if stage=="coarse": coarse(case,args.data_dir,args.workers)
            elif stage=="fine": fine(case,args.data_dir,args.workers)
            elif stage=="chern": chern_from_saved_meshes(case,args.data_dir)
            else: mesh_and_chern(case,args.data_dir,args.workers,args.edge_thresh,args.max_iters)


if __name__=="__main__":
    mp.freeze_support()
    main()
