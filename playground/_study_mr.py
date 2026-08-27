"""
Study nexus structure + MR positions. Trace MR evolution across E=0.
"""
import sys, numpy as np, networkx as nx
sys.path.insert(0, '.')
sys.path.insert(0, 'demos')
from gbz_types import CharPoly, LineSubset, TWO_PI
from continuation import ZeroManager
from trivial_model import get_model, DEFAULT_PARAMS
from demo_nexus import build_nexus, nexus_info
from demo_build_mesh import _solve_nexus_at_E

model = get_model(**DEFAULT_PARAMS)
coeffs, degs = model.get_characteristic_polynomial_data()
char_poly = CharPoly(coeffs, degs)

for E in [-0.001, 0.0, 0.001]:
    zm = ZeroManager(char_poly, E+0j, 0.0)
    zm.run(h0=0.05)

    print(f'=== E={E:+.4f} === {zm.n_segments} seg, {zm.n_multiple_roots} MR '
          f'(has_boundary_mr={zm.has_boundary_mr})')

    for j, mr in enumerate(zm.multiple_roots):
        rts = ', '.join(f'{float(r.real):.4f}{float(r.imag):+.4f}j' for r in mr.roots)
        print(f'  MR[{j}]: theta1={mr.theta1:.6f} beta2={rts}')
        print(f'         cluster_indices={mr.cluster_indices}')

    # Solve nexus at this E
    G = _solve_nexus_at_E(coeffs, degs, E)
    info = nexus_info(G)
    print(f'  Nexus: {info["n_nodes"]} nodes, {info["n_edges"]} edges, '
          f'{len(info["components"])} components')

    # Merge nodes (nodes with >1 line_id)
    merge_nodes = [(n, d) for n, d in G.nodes(data=True)
                   if len(d.get('line_ids', [])) > 1]
    print(f'  Merge nodes: {len(merge_nodes)}')
    for n, ndata in merge_nodes:
        b1, b2 = ndata['beta1'], ndata['beta2']
        th1 = float(np.angle(b1)) % (TWO_PI)
        th2 = float(np.angle(b2)) % (TWO_PI)
        print(f'    theta1={th1:.6f} theta2={th2:.6f} beta2={b2:.6f}')

    # For each MR, find the closest merge node
    if merge_nodes:
        print(f'  MR → merge node mapping:')
        for j, mr in enumerate(zm.multiple_roots):
            mr_th1 = mr.theta1 % (TWO_PI)
            mr_b2 = mr.roots[0]
            best = None
            for n, ndata in merge_nodes:
                th1 = float(np.angle(ndata['beta1'])) % (TWO_PI)
                d = min(abs(th1 - mr_th1), TWO_PI - abs(th1 - mr_th1))
                if best is None or d < best_d:
                    best, best_d = (n, ndata), d
            if best is not None:
                n, ndata = best
                print(f'    MR[{j}] theta1={mr_th1:.6f} → node theta1='
                      f'{float(np.angle(ndata["beta1"]))%(TWO_PI):.6f} '
                      f'(dist={best_d:.2e})')

    # For E=-0.001 and E=0.001, also test: removing merge nodes splits cycle
    if len(merge_nodes) >= 2:
        G_test = G.copy()
        for n, _ in merge_nodes:
            G_test.remove_node(n)
        comps_after = list(nx.connected_components(G_test))
        print(f'  After removing merge nodes: {len(comps_after)} components')
        for ci, comp in enumerate(comps_after):
            th1s = [np.angle(G_test.nodes[n]['beta1'])%(TWO_PI) for n in comp]
            th2s = [np.angle(G_test.nodes[n]['beta2'])%(TWO_PI) for n in comp]
            print(f'    comp[{ci}]: {len(comp)} nodes, '
                  f'th1=[{min(th1s):.4f}, {max(th1s):.4f}], '
                  f'th2=[{min(th2s):.4f}, {max(th2s):.4f}]')

    print()
