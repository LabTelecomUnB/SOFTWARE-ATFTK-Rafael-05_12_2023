#  -*- coding: utf-8 -*-
"""
Author: Rafael R. L. Benevides
"""

from numpy cimport ndarray


cpdef ndarray[double complex, ndim=2] get_eigenantennas(int lmax,
                                                        double theta,
                                                        double phi)

cpdef ndarray[double, ndim=1] directivity(ndarray[double complex, ndim=2] q,
                                          ndarray[double, ndim=1] theta,
                                          ndarray[double, ndim=1] phi)