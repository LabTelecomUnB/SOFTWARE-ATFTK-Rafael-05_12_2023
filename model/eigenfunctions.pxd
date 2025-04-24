#  -*- coding: utf-8 -*-
"""
Author: Rafael R. L. Benevides
"""

cimport numpy as numpyc


cdef inline double complex Ylm(int l, int m, double theta, double phi)

cpdef numpyc.ndarray[numpyc.complex128_t, ndim=2] Tk(int k,
                                                     double theta,
                                                     double phi,
                                                     double eta)

cpdef numpyc.ndarray[numpyc.complex128_t, ndim=2] Upsilonk(int k,
                                                           double[:] theta,
                                                           double[:] phi,
                                                           double eta)