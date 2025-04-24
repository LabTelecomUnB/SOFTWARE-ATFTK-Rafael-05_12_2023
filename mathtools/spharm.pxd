#  -*- coding: utf-8 -*-
"""

Author: Rafael Rodrigues Luz Benevides
Date: 8/16/23

"""

from numpy cimport ndarray


cdef double eta

cpdef ndarray[double complex, ndim=2] T(int lmax, double theta, double phi)