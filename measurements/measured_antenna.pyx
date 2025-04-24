#  -*- coding: utf-8 -*-
"""
Author: Rafael R. L. Benevides
"""


cimport numpy as numpyc

from numpy.typing import ArrayLike


cdef class ModeCoefficients:
    r"""

    """

    def plot_power_distribution(self):
        ...

    def plot_radiation_pattern(self):
        ...

    def as_ndarray(self):
        ...

    @property
    def radiated_power(self) -> ArrayLike:
        ...

    @property
    def total_radiated_power(self) -> double:
        ...


cdef class RadiationVectorData:
    r"""
    Measured samples of an antenna radiation vector



    .. math::
        \mathbf{E}\left(r, \theta, \phi \right)
        =
        \frac{
            e^{-jkr}
        }{
            r
        }
        \mathbf{\mathtt{E}}\left( \theta, \phi \right)
    """

    def __cinit__(self,
                  numpyc.float64_t[:] theta,
                  numpyc.float64_t[:] phi,
                  numpyc.complex128_t[:] E_theta,
                  numpyc.complex128_t[:] E_phi):

        ...

    def estimate_mode_coeff(self, highest_degree: int, method: str = 'batch')





class MeasuredAntenna:

    def __init__(self) -> None:
        pass


    def add_radiation_vector_data(self,
                                  freq: float,
                                  theta: ArrayLike,
                                  phi: ArrayLike,
                                  E_theta: ArrayLike,
                                  E_phi: ArrayLike):
        pass