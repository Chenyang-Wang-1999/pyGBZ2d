"""Minimal repro: does _join_runs_across_mrs cap merges at n_seg?

Setup: 2 segments (1 interior MR at theta=pi + the 0/2pi seam), a 2-fold
continuum (tracks A and B) closed around the circle.  4 pieces total;
joining into 2 closed loops requires 4 merges (2 per track: MR + seam).
The outer loop `for _ in range(n_seg)` allows at most n_seg = 2 merges.
"""
import numpy as np
from types import SimpleNamespace

from brute_force_SGBZ.continuum_lines import _join_runs_across_mrs
from brute_force_amoeba.zm_extract import _LinePiece

# --- mock zm: 2 segments, no MR clusters, K=2 ---
vals = [1.0, 10.0, 100.0]           # 3-fold continuum: tracks A, B, C
seg0 = SimpleNamespace(
    theta1_arr=np.array([0.0, np.pi / 2, np.pi]),
    tracked_roots=np.array([[v + 0j for v in vals],
                            [v + 0.5 for v in vals],
                            [v + 1.0 for v in vals]]),
    left_mr=-1, right_mr=0,           # seam on the left, MR 0 on the right
)
seg1 = SimpleNamespace(
    theta1_arr=np.array([np.pi, 3 * np.pi / 2, 2 * np.pi]),
    tracked_roots=np.array([[v + 1.0 + 0j for v in vals],   # closed at seam
                            [v + 0.5 for v in vals],
                            [v + 0.0 for v in vals]]),
    left_mr=0, right_mr=-1,           # MR 0 on the left, seam on the right
)
zm = SimpleNamespace(segments=[seg0, seg1], multiple_roots=[], K=3,
                     _event_groups=[])

# 6 pieces: per segment per track (what _runs_to_pieces would emit)
pieces = []
for s_idx, seg in enumerate((seg0, seg1)):
    for j in range(len(vals)):
        pieces.append(_LinePiece(
            0j, 0.0, seg.theta1_arr.copy(),
            seg.tracked_roots[:, j].copy(), ml=s_idx, mr=s_idx))

joined = _join_runs_across_mrs(zm, pieces)
print(f"n_seg = 2, pieces in = 6 (3-fold continuum), merges needed = 3")
print(f"pieces out = {len(joined)}  (expected 3 closed loops)")
for p in joined:
    print(f"  ml={p.ml} mr={p.mr} theta=[{p.theta1_arr[0]:.3f},"
          f" {p.theta1_arr[-1]:.3f}] b2=[{p.beta2_arr[0]:.2f},"
          f" {p.beta2_arr[-1]:.2f}] len={len(p.theta1_arr)}")
