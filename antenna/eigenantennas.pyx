#  -*- coding: utf-8 -*-
"""

Author: Rafael R. L. Benevides
Date: 8/16/23

"""


# ---------- ---------- ---------- ---------- ---------- ---------- py
import numpy
from scipy.linalg import svd

# ---------- ---------- ---------- ---------- ---------- ---------- cy
from numpy cimport ndarray
from libc.math cimport sqrt, pi

from aftk.mathtools.spharm cimport T, eta
# from aftk.model.eigenfunctions import Tk


# ---------- ---------- ---------- ---------- ---------- ---------- ty
cpdef ndarray[double complex, ndim=2] get_eigenantennas(int lmax,
                                                        double theta,
                                                        double phi):
    """
    
    Parameters
    ----------
    lmax
    theta
    phi

    Returns
    -------

    """

    cdef:
        ndarray[double complex, ndim=2] _T = T(lmax, theta, phi)
        ndarray[double complex, ndim = 2] U, QH
        ndarray[double, ndim = 1] sigma_diag

    U, sigma_diag, QH = svd(_T * sqrt(4*pi/eta))

    return QH.conjugate().T


cpdef ndarray[double, ndim=1] directivity(ndarray[double complex, ndim=2] q,
                                          ndarray[double, ndim=1] theta,
                                          ndarray[double, ndim=1] phi):
    """
    
    Parameters
    ----------
    q
    theta
    phi

    Returns
    -------

    """
    cdef:
        int lmax = <int> sqrt(q.shape[0] / 2 + 1) - 1
        Py_ssize_t index

        double qHq = (q.conjugate().T @ q)[0, 0].real

        ndarray[double, ndim=1] D = numpy.zeros_like(theta)
        ndarray[double complex, ndim=2] E

    for index in range(theta.shape[0]):

        E = T(lmax, theta[index], phi[index]) @ q
        # E = Tk(lmax, theta[index], phi[index], eta) @ q

        D[index] = (E.conjugate().T @ E)[0, 0].real

    return 4*pi/eta * D / qHq





