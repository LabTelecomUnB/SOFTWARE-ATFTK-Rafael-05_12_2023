#  -*- coding: utf-8 -*-
"""

Author: Rafael R. L. Benevides
Date: 21/07/2023

"""

import numpy

from matplotlib import pyplot

from scipy.sparse.linalg import eigsh

from warnings import warn

import pyvista

import scipy

from scipy.linalg import svd

# ---------- ---------- ---------- ---------- ---------- ----------
cimport numpy as numpyc

from libc.math cimport sqrt

# ---------- ---------- ---------- ---------- ---------- ----------
from numpy.typing import ArrayLike, NDArray


def plot_matrix(A: ArrayLike, cmap_mag='jet'):
    """this function is old and must be replaced by a new one in RULE.py"""
    if numpy.issubdtype(A.dtype, numpy.complexfloating):
        figure, axes = pyplot.subplots(2, 1, sharex=True, sharey=True)
        figure.set_size_inches(8, 12)

        A_mag = numpy.abs(A)
        A_pha = numpy.angle(A, deg=True)

        mappable = axes[0].imshow(A_mag, cmap=cmap_mag, origin='lower')

    else:
        figure, axes = pyplot.subplots()


def condition_number_of_PDM(numpyc.ndarray[numpyc.complex128_t, ndim=2] A):
    """Calculates the condition number of a positive-definite matrix.

    Parameters
    ----------
    A : (N, N) shaped NDArray
        Positive definite matrix

    Returns
    -------
    out : float
        condition number (ratio between the greatest and smallest eigenvalue)

    """
    cdef double lambda_max, lambda_min

    lambda_max = eigsh(A, k=1, which='LA', return_eigenvectors=False)
    lambda_min = eigsh(A, k=1, which='SA', return_eigenvectors=False)

    if lambda_min <= 0:
        warn(f'Given Matrix is not positive-definite. '
             f'Smallest eigenvalue: {lambda_min}. '
             f'Using numpy.linalg.cond instead')

        return numpy.linalg.cond(A)

    return lambda_max / lambda_min


# ========== ========== ========== ========== ========== ========== planning directions
def create_spherical_random_points(int N) -> NDArray:

    cdef numpyc.ndarray[double, ndim=2] r = numpy.random.randn(N, 3)

    return r / numpy.linalg.norm(r, axis=1).reshape(-1, 1)


def create_grid_points(int n_theta, int n_phi) -> NDArray:

    cdef:
        numpyc.ndarray[double, ndim=2] r = numpy.zeros((n_theta*n_phi, 3), dtype=numpy.double)
        numpyc.ndarray[double, ndim=1] theta = numpy.linspace(0, numpy.pi, n_theta+2)[1:-1]
        numpyc.ndarray[double, ndim=1] phi = numpy.linspace(0, 2*numpy.pi, n_phi+1)[:-1]

        numpyc.ndarray[double, ndim=2] THETA, PHI, X, Y, Z

    THETA, PHI = numpy.meshgrid(theta, phi)

    X = numpy.cos(PHI)*numpy.sin(THETA)
    Y = numpy.sin(PHI)*numpy.sin(THETA)
    Z = numpy.cos(THETA)

    X = X.reshape(-1, 1)
    Y = Y.reshape(-1, 1)
    Z = Z.reshape(-1, 1)

    return numpy.hstack([X, Y, Z])


def create_directions(numpyc.ndarray[double, ndim=2] r0, double max_error=1e-2, int maxiter = 50) -> NDArray:

    cdef:
        numpyc.ndarray[double, ndim=2] r = r0
        numpyc.ndarray[double, ndim=2] r_next = numpy.zeros_like(r)

        double aux, energy, error = 1.0, energy_plus, energy_minus

        N = r0.shape[0]

        Py_ssize_t i, j
        int count = 0

    while error > max_error and count < maxiter:
        count += 1

        energy = 0.0

        for i in range(N):

            r_next[i, 0] = 0.0
            r_next[i, 1] = 0.0
            r_next[i, 2] = 0.0

            for j in range(N):
                if i != j:

                    aux = r[i, 0]*r[j, 0] + r[i, 1]*r[j, 1] + r[i, 2]*r[j, 2]
                    if aux >= 1.0:
                        aux = 0.999999

                    aux = sqrt(1.0 - aux)

                    energy+= 1/aux

                    aux = aux * aux * aux

                    r_next[i, 0] += r[j, 0] / aux
                    r_next[i, 1] += r[j, 1] / aux
                    r_next[i, 2] += r[j, 2] / aux

            aux = sqrt(r_next[i, 0]*r_next[i, 0] + r_next[i, 1]*r_next[i, 1] + r_next[i, 2]*r_next[i, 2])

            r_next[i, 0] = r_next[i, 0] / aux
            r_next[i, 1] = r_next[i, 1] / aux
            r_next[i, 2] = r_next[i, 2] / aux

            energy_plus = 0.0
            energy_minus = 0.0

            for j in range(N):
                if i != j:

                    aux = r_next[i, 0] * r[j, 0] + r_next[i, 1] * r[j, 1] + r_next[i, 2] * r[j, 2]

                    if aux >= 1.0:
                        aux = 0.999999

                    energy_plus += 1 / sqrt(1.0 - aux)
                    energy_minus += 1 / sqrt(1.0 + aux)

            if energy_minus < energy_plus:
                r_next[i, 0] *= -1
                r_next[i, 1] *= -1
                r_next[i, 2] *= -1

        error = numpy.linalg.norm(r_next - r, axis=1).max()
        print(f'iter: {count:04d}\terror: {error:.5e}\tenergy: {energy:.5e}')

        r = r_next.copy()

    return r_next


def plot_points_3D(numpyc.ndarray[double, ndim=2] r) -> None:

    cdef object plotter

    plotter = pyvista.Plotter()
    plotter.background_color = 'gray'
    plotter.add_points(r, render_points_as_spheres=True)
    plotter.add_axes_at_origin(xlabel='x', ylabel='y', zlabel='z', line_width=10)
    plotter.show()


# ========== ========== ========== ========== ========== ========== evaluate rank


# ========== ========== ========== ========== Positive Definite Matrices
def psdm_eigenvalues(numpyc.ndarray[numpyc.complex128_t, ndim=2] A):
    """
    Calculate the eigenvalues of the matrix $A^H A$.

    Parameters
    ----------
    A

    Returns
    -------

    """

    cdef:
        numpyc.ndarray[double, ndim=1] sigma_diag = svd(A, compute_uv=False)
        numpyc.ndarray[double, ndim=2] sigma = numpy.zeros_like(A, dtype=numpy.double)

    sigma[:sigma_diag.shape[0], :sigma_diag.shape[0]] = numpy.diag(sigma_diag)

    return (sigma.T @ sigma).diagonal()


def psdm_condition_number(numpyc.ndarray[numpyc.complex128_t, ndim=2] A,
                          bint return_nan=False):
    """
    Calculate the condition number of the matrix $A^H A$.

    Parameters
    ----------
    A
    return_nan

    Returns
    -------

    """

    cdef:
        numpyc.ndarray[double, ndim=1] eig = psdm_eigenvalues(A)

        double eig_max = eig.max(), eig_min = eig.min()

    if eig_min > 0.0:
        return eig_max / eig_min

    if eig_min == 0.0:
        return numpy.nan if return_nan else numpy.inf

    raise ValueError(f'Expected all eigenvalues to be non-negative. But found'
                     f' {eig_min}')


def pdm_inverse(numpyc.ndarray[numpyc.complex128_t, ndim=2] A):
    """
    
    Parameters
    ----------
    A

    Returns
    -------

    """
    cdef:
        numpyc.ndarray[double, ndim=1] sigma_diag
        numpyc.ndarray[double, ndim=2] inv_SHS
        numpyc.ndarray[numpyc.complex128_t, ndim=1] U, V, VH

    if A.shape[0] < A.shape[1]:
        raise ValueError('AH@A is positive semi-definite')

    U, sigma_diag, VH = svd(A)
    inv_SHS = numpy.diag(1/sigma_diag/sigma_diag)
    V = VH.conjugate().T

    return V @ inv_SHS @ VH



