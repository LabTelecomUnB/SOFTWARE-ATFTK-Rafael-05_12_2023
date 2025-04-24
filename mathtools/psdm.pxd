#  -*- coding: utf-8 -*-
"""
Author: Rafael R. L. Benevides
"""

from numpy cimport ndarray


cdef int inverse_and_cond(ndarray[double complex, ndim=2] A,
                          ndarray[double complex, ndim=2] inv_A,
                          double *cond_number) except -1