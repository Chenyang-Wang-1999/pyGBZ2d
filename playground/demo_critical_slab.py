'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2026-08-29
Copyright © Department of Physics, Tsinghua University. All rights reserved

Helper functions for the critical-slab patch used by
``demo_torus_prototype.py``.

Near a Morse saddle P of the E-projection, the GBZ surface is locally a
graph over the torus chart (theta1, theta2).  Each flank loop has two
U-shaped passages near P.  These helpers find those passages in the
unwrapped chart and build the cyclic complement arcs between them.
'''

import numpy as np

TWO_PI = 2.0 * np.pi


# ---------------------------------------------------------------------------
# Passage detection in the torus chart
# ---------------------------------------------------------------------------

def angular_dist(a, b):
    return np.abs(np.angle(np.exp(1j * (np.asarray(a) - np.asarray(b)))))


def cyclic_runs(mask):
    """Contiguous True runs on a cyclic boolean mask.

    Returns a list of integer arrays, each in cyclic increasing order.
    """
    mask = np.asarray(mask, dtype=bool)
    n = len(mask)
    if not mask.any():
        return []
    if mask.all():
        return [np.arange(n)]
    start = int(np.where(~mask)[0][0])
    m2 = np.concatenate([mask[start:], mask[:start]])
    runs = []
    i = 0
    while i < n:
        if m2[i]:
            j = i
            while j < n and m2[j]:
                j += 1
            runs.append((np.arange(i, j) + start) % n)
            i = j
        else:
            i += 1
    return runs


def lift_run(loop, run, P):
    """Lift a run's (theta1, theta2) to the universal cover near P."""
    t1 = np.asarray(loop.th1)[run].copy()
    t2 = np.asarray(loop.th2)[run].copy()
    t1u = np.unwrap(t1)
    t2u = np.unwrap(t2)
    k1 = int(round((P[0] - t1u.mean()) / TWO_PI))
    k2 = int(round((P[1] - t2u.mean()) / TWO_PI))
    return np.column_stack([t1u + k1 * TWO_PI, t2u + k2 * TWO_PI])


def passages_on_loop(loop, P, R):
    """All passages of a closed loop within chart-distance R of P.

    Each passage is a dict:
      idx   : loop indices in cyclic increasing order
      uv    : lifted (theta1, theta2) coordinates near P
      s, e  : first / last loop index of the passage (cyclic order)
      dmin  : minimum chart distance to P
    """
    th1 = np.asarray(loop.th1) % TWO_PI
    th2 = np.asarray(loop.th2) % TWO_PI
    d = np.sqrt(angular_dist(th1, P[0]) ** 2
                + angular_dist(th2, P[1]) ** 2)
    runs = cyclic_runs(d < R)
    out = []
    for run in runs:
        uv = lift_run(loop, run, P)
        out.append({
            "idx": run,
            "uv": uv,
            "s": int(run[0]),
            "e": int(run[-1]),
            "dmin": float(d[run].min()),
        })
    return out


# ---------------------------------------------------------------------------
# Cyclic arcs
# ---------------------------------------------------------------------------

def cyclic_range(n, a, b):
    """Loop indices a -> b inclusive, in increasing cyclic order."""
    if b >= a:
        return np.arange(a, b + 1)
    return np.concatenate([np.arange(a, n), np.arange(0, b + 1)])


def complement_arc_from_endpoint(n, passages, L):
    """Complement arc adjacent to passage endpoint L, oriented away from
    the passage (first index = L)."""
    npass = len(passages)
    for pi, p in enumerate(passages):
        if L == p["s"]:
            prev = passages[(pi - 1) % npass]
            arc = cyclic_range(n, prev["e"], p["s"])
            return arc[::-1]           # starts at L = p['s'], ends at prev['e']
        if L == p["e"]:
            nxt = passages[(pi + 1) % npass]
            return cyclic_range(n, p["e"], nxt["s"])
    raise ValueError(f"index {L} is not a passage endpoint")


def find_complement_arc(n, passages, a, b):
    """Complement arc with endpoints {a,b}, oriented a -> b."""
    for pi, p in enumerate(passages):
        nxt = passages[(pi + 1) % len(passages)]
        arc = cyclic_range(n, p["e"], nxt["s"])
        if {int(arc[0]), int(arc[-1])} == {int(a), int(b)}:
            if int(arc[0]) == int(a):
                return arc
            return arc[::-1]
    raise ValueError(f"no complement arc with endpoints {a}, {b}")
