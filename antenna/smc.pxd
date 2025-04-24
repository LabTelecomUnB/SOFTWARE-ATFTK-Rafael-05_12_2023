#  -*- coding: utf-8 -*-
"""
Author: Rafael R. L. Benevides
"""

from numpy cimport ndarray

cpdef ndarray[double, ndim=1] radiation_intensity(ndarray[double complex, ndim=2] q,
                                                  ndarray[double, ndim=1] theta,
                                                  ndarray[double, ndim=1] phi)

cpdef ndarray[double, ndim=1] directivity(ndarray[double complex, ndim=2] q,
                                          ndarray[double, ndim=1] theta,
                                          ndarray[double, ndim=1] phi)