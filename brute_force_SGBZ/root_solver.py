'''
author:        wangchenyang <cy-wang21@mails.tsinghua.edu.cn>
date:          2025-11-18 00:00:00
Copyright © YourCompanyName All rights reserved
'''

'''
Calcualtions related to root solving of polynomials
'''

import numpy as np
import poly_tools as pt
from scipy import optimize


class ComplexEqConverter:
    def __init__(self, fun, cdim) -> None:
        self.fun = fun
        self.cdim = cdim
    
    def __call__(self, x_re, *args):
        x = x_re[:self.cdim] + 1j * x_re[self.cdim:]
        eq_LHS, eq_jac = self.fun(x, *args)

        # convert to real equation
        eq_LHS_re = np.zeros(2*self.cdim)
        eq_LHS_re[:self.cdim] = eq_LHS.real
        eq_LHS_re[self.cdim:] = eq_LHS.imag

        eq_jac_re = np.zeros((2*self.cdim, 2*self.cdim))
        eq_jac_re[:self.cdim, :self.cdim] = eq_jac.real 
        eq_jac_re[:self.cdim, self.cdim:] = - eq_jac.imag
        eq_jac_re[self.cdim:, :self.cdim] = eq_jac.imag
        eq_jac_re[self.cdim:, self.cdim:] = eq_jac.real
        return eq_LHS_re, eq_jac_re


def complex_root(complex_fun, x0, **options):
    fun = ComplexEqConverter(complex_fun, len(x0))
    x0_re = np.zeros(2*fun.cdim)
    x0_re[:fun.cdim] = x0.real
    x0_re[fun.cdim:] = x0.imag
    res= optimize.root(fun, x0_re, **options)
    x_re = res.x
    # jac_re = res.jac
    res.x = x_re[:fun.cdim] + 1j * x_re[fun.cdim:]
    if('args' in options.keys()):
        _, jac = complex_fun(res.x,*options['args'])
    else:
        _, jac = complex_fun(res.x)
    res.jac = jac
    return res


def poly_to_np_coefficients(coeffs: list[complex], degs: list[int]) -> np.ndarray:
    '''
    Convert a polynomial to a numpy array of coefficients.

    Parameters:
    coeffs (list[complex]): The coefficients of the polynomial.
    degs (list[int]): The degrees of the coefficients.

    Returns:
    np.ndarray: The numpy array of coefficients.
    '''
    max_deg = max(degs)
    np_coeffs = np.zeros(max_deg + 1, dtype=complex)
    for coeff, deg in zip(coeffs, degs):
        np_coeffs[max_deg - deg] = coeff
    return np_coeffs


def calculate_point_roots(
    char_poly: pt.CLaurent,
    param_ind_ctype: pt.CIndexVec,
    param_val: np.ndarray,
    var_ind_ctype: pt.CIndexVec,
    M_max: int,
    N_max: int
):
    poly_1d = char_poly.partial_eval(
        pt.CScalarVec(param_val),
        param_ind_ctype,
        var_ind_ctype
    )
    coeffs = pt.CScalarVec([])
    degs = pt.CIndexVec([])
    poly_1d.num.batch_get_data(coeffs, degs)
    deg_M = poly_1d.denom_orders[0]
    np_coeffs = poly_to_np_coefficients(coeffs, degs)
    curr_roots = list(np.roots(np_coeffs))
    new_root_len = len(curr_roots)
    if new_root_len < M_max + N_max:
        if deg_M < M_max:
            curr_roots.append(0)
        if new_root_len - deg_M < N_max:
            curr_roots.append(np.inf)
    return curr_roots
