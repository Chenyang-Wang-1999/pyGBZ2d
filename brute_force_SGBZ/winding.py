'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2025-11-18
Copyright © Department of Physics, Tsinghua University. All rights reserved
'''

''' Functions for winding number calculations '''

import numpy as np

from scipy import integrate
from math import pi
from scipy import sparse
from scipy.sparse import linalg as spla

from gbz_types import CharPoly


class WindingFun:
    """Winding integrand of the characteristic polynomial along a loop.

    Calling the instance at parameter t returns Im[f'(t)/f(t)], i.e.
    d/dt Im log f(param(t)); its integral over loop_range divided by 2*pi
    is the winding number of f around zero (see get_winding_number).
    """
    char_poly: CharPoly
    loop_fun: callable  # Input: t. Output: param, dparam/dt
    loop_range: tuple[float, float]

    def __init__(self, char_poly: CharPoly, loop_fun: callable, loop_range: tuple[float, float]):
        """
        Parameters:
            char_poly: characteristic Laurent polynomial (CharPoly instance).
            loop_fun: loop in variable space, t -> (param, dparam/dt) where
                param collects the polynomial variables at t.
            loop_range: (t_start, t_end) parameter range of the closed loop.
        """
        self.char_poly = char_poly
        self.loop_fun = loop_fun
        self.loop_range = loop_range

    def __call__(self, t: float) -> complex:
        param, dparam_dt = self.loop_fun(t)
        val = self.char_poly.eval_val(param)
        partials = self.char_poly.eval_partials(param)
        dval_dt = sum(partials[param_ind] * dparam_dt[param_ind] for param_ind in range(len(partials)))
        return (dval_dt / val).imag


class MatWindingFun(WindingFun):
    """Winding integrand of det(E_ref - H(beta)) built from the Bloch matrix.

    Matrix counterpart of WindingFun for models given as H(beta) rather
    than an explicit characteristic polynomial, using
    d/dt log det(E - H) = tr[(E - H)^{-1} (-dH/dt)].
    Only the sparse-matrix path is implemented.
    """
    mat_fun: callable  # foo(beta)
    dmat_fun: callable  # foo(beta, diff_orders)
    param_fun: callable # t -> beta(t), dbeta(t)
    E_ref: complex
    is_sparse: bool

    def __init__(self, mat_fun: callable, dmat_fun: callable, param_fun, loop_range, E_ref: complex, is_sparse: bool = True):
        self.mat_fun = mat_fun
        self.dmat_fun = dmat_fun
        self.param_fun = param_fun
        self.is_sparse = is_sparse
        self.loop_range = loop_range
        self.E_ref = E_ref

    def call_sparse(self, t: float) -> complex:
        beta, dbeta_dt = self.param_fun(t)
        mat = self.mat_fun(beta)
        dmat_dt = sparse.csc_matrix(mat.shape, dtype=complex)
        for dim_ind in range(len(beta)):
            diff_orders = np.zeros(len(beta), dtype=int)
            diff_orders[dim_ind] = 1
            dmat_dt += self.dmat_fun(beta, diff_orders) * dbeta_dt[dim_ind]

        try:
            return spla.spsolve(self.E_ref * sparse.eye(mat.shape[0]) - mat, -dmat_dt).trace().imag
        except RuntimeError:
            return -np.inf
    
    def __call__(self, t: float) -> complex:
        if self.is_sparse:
            return self.call_sparse(t)
        else:
            raise NotImplementedError("Dense matrix method is not implemented yet.")


def get_winding_number(winding_fun: WindingFun, N_seg=1) -> float:
    '''
    Integrate winding_fun over its loop range and divide by 2*pi.

    The loop is split into N_seg equal segments, each handled by
    scipy.integrate.quad; splitting helps convergence when the integrand
    is sharply peaked (loop passing close to a root).

    Returns:
        Real-valued winding number (unrounded; integer only in exact
        arithmetic — callers round as needed).
    '''
    int_tot = 0.0
    intervals = np.linspace(winding_fun.loop_range[0], winding_fun.loop_range[1], N_seg + 1)
    for i in range(N_seg):
        int_tot += integrate.quad(winding_fun, intervals[i], intervals[i+1], epsabs=1e-2, epsrel=1e-2, limit=500)[0]

    return int_tot / (2 * pi)
