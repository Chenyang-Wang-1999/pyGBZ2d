'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2025-11-18
Copyright © Department of Physics, Tsinghua University. All rights reserved
'''

''' Functions for winding number calculations '''

import poly_tools as pt
import numpy as np

from scipy import integrate
from math import pi
from scipy import sparse
from scipy.sparse import linalg as spla


class PolyDiffContext:
    """Bundle polynomial value/derivative evaluation for reuse."""
    char_poly: pt.CLaurent
    dchar_poly: list[pt.CLaurent]

    def __init__(self, char_poly: pt.CLaurent):
        self.char_poly = char_poly
        self.dchar_poly = [char_poly.derivative(var_ind) for var_ind in range(char_poly.dim)]

    def eval_val(self, var: tuple[complex]) -> complex:
        return self.char_poly.eval(pt.CScalarVec(var))

    def eval_partials(self, var: tuple[complex]) -> list[complex]:
        var_ctype = pt.CScalarVec(var)
        return [dchar.eval(var_ctype) for dchar in self.dchar_poly]
    
    def eval_dmu2(self, var: tuple[complex]):
        '''Reture (dmu2/dmu1, dmu2/dtheta1)'''
        partials = self.eval_partials(var)
        complex_diff = var[1] * partials[1] / (var[2] * partials[2])
        return (-complex_diff.real, complex_diff.imag)


class WindingFun:
    """Winding integrand of the characteristic polynomial along a loop.

    Calling the instance at parameter t returns Im[f'(t)/f(t)], i.e.
    d/dt Im log f(param(t)); its integral over loop_range divided by 2*pi
    is the winding number of f around zero (see get_winding_number).
    """
    poly_diff: PolyDiffContext
    char_poly: pt.CLaurent
    loop_fun: callable  # Input: t. Output: param, dparam/dt
    dchar_poly: list[pt.CLaurent]
    loop_range: tuple[float, float]

    def __init__(self, char_poly: pt.CLaurent, loop_fun: callable, loop_range: tuple[float, float]):
        """
        Parameters:
            char_poly: characteristic Laurent polynomial f(vars).
            loop_fun: loop in variable space, t -> (param, dparam/dt) where
                param collects the polynomial variables at t.
            loop_range: (t_start, t_end) parameter range of the closed loop.
        """
        self.poly_diff = PolyDiffContext(char_poly)
        # Keep compatibility with existing code that accesses these attributes.
        self.char_poly = self.poly_diff.char_poly
        self.dchar_poly = self.poly_diff.dchar_poly
        self.loop_fun = loop_fun
        self.loop_range = loop_range

    def __call__(self, t: float) -> complex:
        param, dparam_dt = self.loop_fun(t)
        val = self.poly_diff.eval_val(param)
        partials = self.poly_diff.eval_partials(param)
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
