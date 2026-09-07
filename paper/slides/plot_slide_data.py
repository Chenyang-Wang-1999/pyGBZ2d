"""Plot computed JSON from collect_slide_data.py as vector SVG and PDF.

Never interpolates missing energies or manufactures failed solver values.
Requires numpy and matplotlib. No model solver is invoked by this script.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def complex_array(values):
    a = np.array(values, dtype=float)
    return a[..., 0] + 1j*a[..., 1]


def torus_segments(x, y):
    cuts = np.flatnonzero((np.abs(np.diff(x)) > np.pi) | (np.abs(np.diff(y)) > np.pi)) + 1
    return zip(np.split(x, cuts), np.split(y, cuts))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input", type=Path)
    p.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent / "Figures")
    p.add_argument("--energy-index", type=int, default=0, help="index in distinct sweep energies for the subset plot")
    args = p.parse_args()
    data = json.loads(args.input.read_text(encoding="utf-8"))
    if data.get("schema_version") != 1:
        p.error("unsupported data schema")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 13,
                         "svg.fonttype": "none", "pdf.fonttype": 42,
                         "axes.spines.top": False, "axes.spines.right": False})
    emitted = []

    def save(fig, name):
        fig.tight_layout()
        for suffix in ("svg", "pdf"):
            fig.savefig(args.output_dir / f"{name}.{suffix}", bbox_inches="tight")
        plt.close(fig)
        emitted.append(name)

    tracks = data["tracks"]
    if tracks["segments"]:
        fig, ax = plt.subplots(figsize=(11.2, 4.6))
        for seg in tracks["segments"]:
            roots = complex_array(seg["beta2"])
            with np.errstate(divide="ignore", invalid="ignore"):
                logabs = np.log(np.abs(roots))
            logabs[~np.isfinite(logabs)] = np.nan
            for j in range(roots.shape[1]):
                ax.plot(seg["theta1"], logabs[:, j], lw=1.3, color=f"C{j%10}")
        for theta in tracks.get("multiple_root_angles", []):
            ax.axvline(theta, color="0.6", ls=":", lw=.8)
        ax.set(xlabel=r"$\theta_1$", ylabel=r"$\ln|\beta_2|$",
               title=f"Computed root tracks: E={complex(*tracks['energy'])}, mu1={tracks['mu1']:g}")
        ax.grid(alpha=.15)
        save(fig, "root-tracks")
    else:
        print("Root-track figure omitted:", tracks.get("error", "no tracks"))

    entries = data["sweep"]
    energies = list(dict.fromkeys(tuple(e["energy"]) for e in entries))
    if not 0 <= args.energy_index < len(energies):
        p.error("energy-index outside the computed sweep")
    chosen = energies[args.energy_index]
    fig, ax = plt.subplots(figsize=(8, 5.4))
    colors = {"sgbz": "#3d8dff", "amoeba": "#00978a"}
    for e in entries:
        if tuple(e["energy"]) != chosen:
            continue
        labeled = False
        for sub in e["subsets"]:
            x = np.angle(complex_array(sub["beta1"])) % (2*np.pi)
            y = np.angle(complex_array(sub["beta2"])) % (2*np.pi)
            for xx, yy in torus_segments(x, y):
                ax.plot(xx, yy, color=colors[e["method"]],
                        marker="o" if sub["kind"] == "point" or len(xx)==1 else None,
                        ls="none" if sub["kind"] == "point" else "-",
                        ms=5, lw=1.4, label=e["method"] if not labeled else None)
                labeled = True
        if not labeled:
            status = "failed" if not e["success"] else "empty"
            ax.plot([], [], color=colors[e["method"]], label=f"{e['method']} ({status})")
    ax.set(xlim=(0,2*np.pi), ylim=(0,2*np.pi), xlabel=r"$\theta_1$", ylabel=r"$\theta_2$",
           title=f"Computed GBZ subsets: E={complex(*chosen)}")
    ax.legend(); ax.grid(alpha=.15)
    save(fig, "subsets")

    fig, ax = plt.subplots(figsize=(10.5, 4.9))
    for method in dict.fromkeys(e["method"] for e in entries):
        for status in ("inside", "outside", "failed"):
            picked = [e for e in entries if e["method"]==method and
                      ("failed" if not e["success"] else "inside" if e["is_gbz"] else "outside")==status]
            if not picked:
                continue
            xy = np.array([e["energy"] for e in picked])
            marker = "x" if status=="failed" else ("o" if method=="sgbz" else "s")
            color = "#c52a37" if status=="failed" else colors[method] if status=="inside" else "0.65"
            kw = {"facecolors": "none"} if status=="inside" else {}
            ax.scatter(xy[:,0], xy[:,1], marker=marker, s=38 if status=="inside" else 18,
                       color=color, label=f"{method}: {status} ({len(picked)})", **kw)
    ax.set(xlabel="Re E", ylabel="Im E", title="Computed energy samples (no interpolation)")
    ax.legend(fontsize=10, loc="best"); ax.grid(alpha=.15)
    save(fig, "spectrum")
    provenance = {"input": str(args.input.resolve()), "sha256": hashlib.sha256(args.input.read_bytes()).hexdigest(),
                  "model": data["model"], "revision": data["revision"], "created_utc": data["created_utc"],
                  "backend": data["backend"], "solver_options": data["solver_options"],
                  "subset_energy": list(chosen), "Figures": emitted,
                  "failed_runs": sum(not e["success"] for e in entries)}
    (args.output_dir / "plot-provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    print("Saved vector Figures:", ", ".join(emitted))


if __name__ == "__main__":
    main()
