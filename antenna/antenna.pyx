#  -*- coding: utf-8 -*-
"""
Author: Rafael R. L. Benevides
"""

from numpy.typing import ArrayLike


cdef class Antenna:

    cdef readonly double[3] position
    cdef readonly double roll, pitch, yaw

    def radiation_field(self, theta: ArrayLike, phi: ArrayLike) -> ArrayLike:
        """
        :param theta:
        :param phi:
        :return:
        """
