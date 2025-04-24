#  -*- coding: utf-8 -*-
"""
Author: Rafael R. L. Benevides
"""

# ---------- ---------- ---------- ---------- ---------- ---------- py
import numpy

from scipy.linalg import svd

# ---------- ---------- ---------- ---------- ---------- ---------- cy
cimport numpy as numpyc

from numpy cimport ndarray, complex128_t

# ---------- ---------- ---------- ---------- ---------- ---------- ty
from numpy.typing import NDArray


cpdef ndarray[double, ndim=1] eigenvalues(ndarray[complex128_t, ndim=2] A):
    """Calculates the eigenvalues of the PSD matrix :math:`M = A^H A`.

    Parameters
    ----------
    A : complex-dtyped (m, n)-shaped :py:class:`ndarray <numpy.ndarray>`
        Matrix used to retrieve :math:`M` from the decomposition :math:`A^H A`.

    Returns
    -------
    eigv : double-dtyped (n,)-shaped :py:class:`ndarray <numpy.ndarray>`
        List of eigenvalues of :math:`A^H A`.
    """

    cdef:
        numpyc.ndarray[double, ndim=1] sigma_diag = svd(A, compute_uv=False)
        numpyc.ndarray[double, ndim=2] sigma = numpy.zeros_like(A, dtype=numpy.double)

    sigma[:sigma_diag.shape[0], :sigma_diag.shape[0]] = numpy.diag(sigma_diag)

    return (sigma.T @ sigma).diagonal()


cpdef double condition_number(ndarray[complex128_t, ndim=2] A):
    """ Calculates the condition number of the PSD matrix :math:`M = A^H A`.
    
    Parameters
    ----------
    A : complex-dtyped (m, n)-shaped :py:class:`ndarray <numpy.ndarray>`
        Matrix used to retrieve :math:`M` from the decomposition :math:`A^H A`.

    Returns
    -------
    cond : :py:class:`double <numpy.float_>`
        Condition number of :math:`A^H A`.

    """

    cdef:
        numpyc.ndarray[double, ndim=1] eig = eigenvalues(A)

        double eig_max = eig.max(), eig_min = eig.min()

    if eig_min > 0.0:
        return eig_max / eig_min

    if eig_min == 0.0:
        return numpy.inf

    raise ValueError(f'Expected all eigenvalues to be non-negative. '
                     f'But found {eig_min}')


cpdef ndarray[double complex, ndim=2] inverse(ndarray[complex128_t, ndim=2] A):
    """Calculates the inverse of the PD matrix :math:`M = A^H A`.

    Parameters
    ----------
    A : complex-dtyped (m, n)-shaped :py:class:`ndarray <numpy.ndarray>`
        Matrix used to retrieve :math:`M` from the decomposition :math:`A^H A`.

    Returns
    -------
    inverse : complex-dtyped (m, m)-shaped :py:class:`ndarray <numpy.ndarray>`
        inverse matrix of :math:`A^H A` if it exists.
    """

    cdef:
        numpyc.ndarray[double, ndim=1] sigma_diag
        numpyc.ndarray[double, ndim=2] inv_SHS
        numpyc.ndarray[numpyc.complex128_t, ndim=2] U, V, VH

    if A.shape[0] < A.shape[1]:
        raise ValueError('AH@A is positive semi-definite')

    U, sigma_diag, VH = svd(A)
    inv_SHS = numpy.diag(1 / sigma_diag / sigma_diag)
    V = VH.conjugate().T

    return V @ inv_SHS @ VH


cdef int inverse_and_cond(ndarray[double complex, ndim=2] A,
                          ndarray[double complex, ndim=2] inv_A,
                          double *cond_number) except -1:

    """
    Calculates the inverse and the condition number of A.H @ A
    """

    cdef:
        numpyc.ndarray[double, ndim=1] sing_value_diag, eigenvalue_diag
        numpyc.ndarray[double, ndim=2] inv_SHS
        numpyc.ndarray[numpyc.complex128_t, ndim=2] U, V, VH

        double eigenvalue_max, eigenvalue_min

    if A.shape[0] <= A.shape[1]:

        inv_A = numpy.full_like(A, numpy.nan + 1j*numpy.nan)
        cond_number[0] = numpy.inf

        raise ValueError('A.H@A is singular')

    U, sing_value_diag, VH = svd(A)
    eigenvalue_diag = sing_value_diag * sing_value_diag

    eigenvalue_max = eigenvalue_diag.max()
    eigenvalue_min = eigenvalue_diag.min()

    if eigenvalue_min == 0.0:

        inv_A = numpy.full_like(A, numpy.nan + 1j * numpy.nan)
        cond_number[0] = numpy.inf

        raise ValueError('A.H@A is singular')

    inv_SHS = numpy.diag(1 / eigenvalue_diag)
    V = VH.conjugate().T

    inv_A = V @ inv_SHS @ VH
    cond_number[0] = eigenvalue_max / eigenvalue_min

    return 0

