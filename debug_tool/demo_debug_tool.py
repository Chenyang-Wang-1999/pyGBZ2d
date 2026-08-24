'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-08-24
Copyright © Department of Physics, Tsinghua University. All rights reserved

Demo / CLI for debug_tool: fixed-(E_ref, mu1) GBZ subsets + winding loops.

Examples
--------
Haldane gain-loss supercell at the E=1.212 debug point (the case from
log/2026-08-16 and simple-plot.py)::

    python debug_tool/demo_debug_tool.py                         # haldane default
    python debug_tool/demo_debug_tool.py --theta2 0.3 1.0 2.0 4.0

2D Hatano-Nelson chain pair (fast, analytic answer mu1 = gamma1)::

    python debug_tool/demo_debug_tool.py --model hn2d --E 1.0 --mu1 0.2
    python debug_tool/demo_debug_tool.py --model hn2d --E 1.0 --mu1 0.25

Use ``--methods sgbz`` / ``--methods amoeba`` to run a single method, and
``--save PATH`` / ``--show`` to control the figure output.
'''

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np


def build_haldane_polynomial():
    """Haldane gain-loss 2-cell supercell characteristic polynomial."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "haldane", ROOT / "demos" / "Haldane-model-gainloss.py")
    h = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(h)
    model = h.Haldane_non_Hermitian_phase(*h.ALL_PARAMS)
    model = model.get_supercell(
        [(0, 0), (1, 0)],
        np.array([[1, 1], [-1, 1]], dtype=int))
    return model.get_characteristic_polynomial_data()


def build_hn2d_polynomial():
    """2D Hatano-Nelson chain pair (conftest's shared test model)."""
    from conftest import build_HN2D_polynomial
    return build_HN2D_polynomial(1.0, 1.0, 0.2, 0.3, 0.0, 0.0)


MODELS = {
    # (builder, default E, default mu1) — the Haldane default is the debug
    # point of log/2026-08-16 (mu1 taken from simple-plot.py); the HN2D
    # default sits exactly on the analytic GBZ radius mu1 = gamma1 = 0.2.
    "haldane": (build_haldane_polynomial, 1.2120000000000002, 0.135328598265),
    "hn2d": (build_hn2d_polynomial, 1.0, 0.2),
}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    parser.add_argument("--model", choices=sorted(MODELS), default="haldane")
    parser.add_argument("--E", type=float, default=None,
                        help="reference energy (default: model's debug point)")
    parser.add_argument("--mu1", type=float, default=None,
                        help="frozen mu1 = ln|beta1| (default: model's debug point)")
    parser.add_argument("--methods", nargs="+", default=["sgbz", "amoeba"],
                        choices=["sgbz", "amoeba"])
    parser.add_argument("--theta2", nargs="*", type=float, default=None,
                        help="loop theta2 values; default: auto gap midpoints")
    parser.add_argument("--n-per-gap", type=int, default=1,
                        help="loops per subset-theta2 gap when auto-picking")
    parser.add_argument("--save", type=str, default=None,
                        help="save the figure to this path (.png)")
    parser.add_argument("--show", action="store_true",
                        help="plt.show() the figure (needs a display)")
    args = parser.parse_args(argv)

    import matplotlib
    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: F401  (backend fixed above)

    from gbz_types import CharPoly
    from debug_tool import collect_debug_subsets, plot_winding_debug

    builder, E_default, mu1_default = MODELS[args.model]
    E_ref = complex(args.E if args.E is not None else E_default)
    mu1 = float(args.mu1 if args.mu1 is not None else mu1_default)

    coeffs, degs = builder()
    poly = CharPoly(coeffs, degs)

    print("=" * 72)
    print(f"model={args.model}  E_ref={E_ref}  mu1={mu1}")
    print("=" * 72)

    report = collect_debug_subsets(poly, E_ref, mu1, methods=args.methods)
    print(report.summary())

    print("\n" + "-" * 72)
    theta2_list = args.theta2  # None -> auto grid per method
    fig, _ax, _results = plot_winding_debug(
        report, theta2_list, n_per_gap=args.n_per_gap)
    print("-" * 72)

    if args.save:
        out = Path(args.save)
    else:
        tag = f"{args.model}_E{E_ref.real:.6f}_mu1{mu1:.9f}"
        out = ROOT / "debug_tool" / "data" / f"gbz_debug_{tag}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=160, bbox_inches="tight")
    print(f"figure saved to {out}")
    if args.show:
        import matplotlib.pyplot as plt
        plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
