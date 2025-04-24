#  -*- coding: utf-8 -*-
"""
Author: Rafael R. L. Benevides
Date: 04/08/2023
"""

# ---------- ---------- ---------- ---------- ---------- ---------- py
import numpy
import scipy.linalg

from scipy.linalg import svd, eigh, pinv
from scipy.optimize import minimize

# ---------- ---------- ---------- ---------- ---------- ---------- cy
cimport numpy as numpyc

from numpy cimport ndarray, complex128_t

from aftk.mathtools.psdm cimport inverse_and_cond

# ---------- ---------- ---------- ---------- ---------- ---------- ty
from numpy.typing import NDArray


# cpdef tuple blue(ndarray[double complex, ndim=2] A,
#     ndarray[double complex, ndim=1] b,
#        ndarray[double, ndim=2] R):
#     """Best Linear Unbiased Estimator
#
#     Parameters
#     ----------
#     A : (s, p)
#     b : (p, 1)
#     R : (s, s)
#
#     Returns
#     -------
#     x : (p, 1)
#     P : (p, p)
#
#     """

# cpdef tuple moore_penrose_inverse(ndarray[double complex, ndim=2] A,
#                                   double cond_max=None):
#     r"""
#
#     Parameters
#     ----------
#     A
#     cond_max
#
#     Returns
#     -------
#
#     """
#     cdef:
#         Py_ssize_t s, p  # samples, parameters
#
#         numpyc.ndarray[double, ndim=1] sigma_diag
#         numpyc.ndarray[double, ndim=1] old_eigenvalues, old_cond
#         numpyc.ndarray[double, ndim=1] new_eigenvalues, new_cond
#
#         numpyc.ndarray[double complex, ndim = 2] U, V, VH, inverse
#
#         double cond
#
#     s, p = A.shape
#
#     if s < p:
#         raise ValueError('Moore-Penrose inverse of A cannot be regularised: '
#                          'A.H@A is singular')
#
#     U, sigma_diag, VH = svd(A)
#
#     if sigma_diag[p-1] == 0:
#         raise ValueError('Moore-Penrose inverse of A cannot be regularised: '
#                          'A has zero-valued singular values')
#
#     old_eigenvalues = sigma_diag*sigma_diag
#
#     old_cond = old_eigenvalues[0] / old_eigenvalues
#
#     if not numpy.any(old_cond > cond_max):
#         inverse_and_cond(A, inverse, &cond)


# def compressed_sensing(A, b):
#
#         """
#
#         Parameters
#         ----------
#         A
#         b
#
#         Returns
#         -------
#
#         """
#
#         # cdef:
#         #     ndarray[double complex, ndim=2] x0 = numpy.linalg.pinv(A) @ b
#         #
#         #     object constraint, result
#
#         x0 = numpy.linalg.pinv(A) @ b
#
#         constraint = {'type': 'eq', 'fun': lambda x : A@x - b}
#
#         return minimize(lambda x : numpy.linalg.norm(x, ord=1),
#                           x0, method='SLSQP', constraints=constraint)


cdef class LinearRegression:
    r"""General Linear Regression Solver

Given a sample array :math:`b\in\mathbb{C}^{s \times 1}`, a design matrix
:math:`A\in\mathbb{C}^{s \times p}`, a positive-definite weight matrix
:math:`W\in\mathbb{C}^{s\times s}` and a Tikhonov matrix
:math:`\Gamma\in\mathbb{C}^{n\times  p}`;

.. math::
    \big\|
        A \widehat{x} - \widetilde{b}
    \big\|_{W}^2
    +
    \big\|
        \Gamma \widehat{x}
    \big\|^2

Parameters
----------
A : complex-dtyped (s, p)-shaped :py:class:`ndarray <numpy.ndarray>`
    The design (or model) matrix.

b : complex-dtyped (s, 1)-shaped :py:class:`ndarray <numpy.ndarray>`
    The samples matrix

W : complex-dtyped (s, s)-shaped :py:class:`ndarray <numpy.ndarray>`
    The weight matrix. Must be positive-definite (no check is performed). If
    `None`, identity matrix is assumed.

Gamma : complex-dtyped (n, p)-shaped :py:class:`ndarray <numpy.ndarray>`
    The Tikhonov-Phillips regularisation matrix. If ``None``, zero matrix is
    assumed (not regularised problem).

unbiased_restriction : :py:class:`bool`
    Only applicable when `Gamma` is set. If ``True``, algorithm will provide
    unbiased estimator.

Notes
-----

The general solution of the problem is given by

.. math::
    \widehat{x}
    =
    \left[
        A^{\mathrm{H}}
        W
        A
        +
        \Gamma^{\mathrm{H}}
        \Gamma
    \right]^{-1}
    \left[
        A^{\mathrm{H}}
        W
        +
        \Lambda
        A^{\mathrm{H}}
    \right]
    \widetilde{b}

where

.. math::
    \Lambda
    =
    \begin{cases}
        0, & \text{no unbiased restriction}\\[3mm]
        \Gamma^{\mathrm{H}}
        \Gamma
        \left[
        A^{\mathrm{H}}
        A
    \right]^{-1}, & \text{unbiased restriction}
    \end{cases}

    """

    # ========== ========== ========== ========== ========== class attributes
    cdef:
        double complex[:, :] _x, _P, _K, _W,_Gamma, _A, _b, _R, _Lambda, _res
        double complex[:, :] _bias_matrix

        readonly double condition_number
        readonly Py_ssize_t s, p

    # ========== ========== ========== ========== ========== special methods
    def __cinit__(self,
                  ndarray[double complex, ndim=2] A,
                  ndarray[double complex, ndim=2] b,
                  ndarray[double complex, ndim=2] W = None,
                  ndarray[double complex, ndim=2] Gamma = None,
                  ndarray[double complex, ndim=2] R = None,
                  bint unbiased_restriction=False) -> None:

        cdef:

            ndarray[double, ndim=1] eigenvalues_W
            ndarray[double, ndim=2] eigenvalues_W_sqrt_diag_inv
            ndarray[double complex, ndim=2] eigenvectors_W, eigenvectors_W_H

            ndarray[double complex, ndim=2] AUX_W, A_H_times_W, Lambda

            ndarray[double complex, ndim=2] moment_matrix

            ndarray[double complex, ndim=2] Gamma_H, Gamma_H_times_Gamma

            double condition_number

            ndarray[double complex, ndim=2] inv_AHA, A_H

        # ---------- ---------- check main matrices dimensions
        self.s = A.shape[0]
        self.p = A.shape[1]

        assert b.shape[0] == self.s, "A and b must have the same number of columns"

        # ---------- ---------- solve
        if Gamma is None:
            self._Gamma = numpy.zeros((1, self.p))

            if W is None:
                self._W = numpy.eye(self.s)

                self._K = pinv(A)
                self._Lambda = numpy.zeros(self.p)

            else:
                self._W = W

                eigenvalues_W, eigenvectors_W = eigh(W)

                eigenvalues_W_sqrt_diag_inv = 1 / numpy.diag(numpy.sqrt(eigenvalues_W))

                eigenvectors_W_H = eigenvectors_W.conjugate().T

                AUX_W = eigenvalues_W_sqrt_diag_inv @ eigenvectors_W_H

                self._K = pinv(AUX_W @ A) @ AUX_W
                self._Lambda = numpy.zeros(self.p)

        else:
            self._Gamma = Gamma

            Gamma_H = Gamma.conjugate().T
            Gamma_H_times_Gamma = Gamma_H @ Gamma

            A_H = A.conjugate().T

            A_H_times_W = A_H @ W

            moment_matrix = A_H_times_W @ A + Gamma_H_times_Gamma

            if unbiased_restriction:

                inv_AHA = numpy.zeros((self.p, self.p), dtype=numpy.complex)

                inverse_and_cond(A, inv_AHA, &condition_number)

                Lambda = Gamma_H_times_Gamma @ inv_AHA

                self._K = scipy.linalg.solve(moment_matrix, A_H_times_W + Lambda @ A_H)

                self._Lambda = Lambda

            else:
                self._K = scipy.linalg.solve(moment_matrix, A_H_times_W)

                self._Lambda = numpy.zeros(self.p)

        self._A = A
        self._b = b
        self._R = numpy.eye(self.s) if R is None else R

        # ---------- ---------- ---------- ---------- ---------- ----------
        self._P = self.K @ self.R @ self.K.conjugate().T
        self._x = self.K @ self.b
        self._res = self.b - self.A @ self.x
        self._bias_matrix = self.K @ self.A - numpy.eye(self.p)

    def __init__(self,
                 ndarray[double complex, ndim=2] A,
                 ndarray[double complex, ndim=2] b,
                 ndarray[double complex, ndim=2] W = None,
                 ndarray[double complex, ndim=2] Gamma = None,
                 ndarray[double complex, ndim=2] R = None,
                 bint unbiased_restriction=False) -> None:
        ...


    # ========== ========== ========== ========== ========== private methods
    ...

    # ========== ========== ========== ========== ========== protected methods
    ...

    # ========== ========== ========== ========== ========== public methods
    def is_unbiased(self):
        return numpy.abs(self._bias_matrix).max() == 0

    # ---------- ---------- ---------- ---------- ---------- properties
    @property
    def A(self) -> NDArray:
        return numpy.asarray(self._A)

    @property
    def b(self) -> NDArray:
        return numpy.asarray(self._b)

    @property
    def x(self) -> NDArray:
        return numpy.asarray(self._x)

    @property
    def res(self) -> NDArray:
        return numpy.asarray(self._res)

    @property
    def R(self) -> NDArray:
        return numpy.asarray(self._R)
    
    @property
    def P(self) -> NDArray:
        return numpy.asarray(self._P)

    @property
    def W(self) -> NDArray:
        return numpy.asarray(self._W)

    @property
    def Gamma(self) -> NDArray:
        return numpy.asarray(self._Gamma)

    @property
    def bias_matrix(self) -> NDArray:
        return numpy.asarray(self._bias_matrix)

    @property
    def condition_number(self) -> tuple:
        raise NotImplementedError()
