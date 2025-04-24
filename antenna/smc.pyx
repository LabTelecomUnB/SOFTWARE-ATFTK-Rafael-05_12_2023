#  -*- coding: utf-8 -*-
"""
Author: Rafael R. L. Benevides
"""

# ---------- ---------- ---------- ---------- ---------- ---------- py
import numpy

import h5py

from pathlib import Path
from scipy.linalg import svd

from abc import ABCMeta, abstractmethod

# ---------- ---------- ---------- ---------- ---------- ---------- cy
from numpy cimport ndarray
from libc.math cimport sqrt, pi

from aftk.mathtools.spharm cimport T, eta

# ---------- ---------- ---------- ---------- ---------- ---------- ty
from numpy.typing import NDArray


cpdef ndarray[double, ndim=1] radiation_intensity(ndarray[double complex, ndim=2] q,
                                                  ndarray[double, ndim=1] theta,
                                                  ndarray[double, ndim=1] phi):
    r"""Calculates the radiation intensity of an antenna.
    
    Parameters
    ----------
    q : complex-dtyped (N,1)-shaped :py:class:`ndarray <numpy.ndarray>`
        Spherical mode coefficients (SMC) in :math:`W^{\frac{1}{2}} 
        (square root of watt)`.
    
    theta : double-dtyped (N,)-shaped :py:class:`ndarray <numpy.ndarray>`
        polar angle in radians.

    phi : double-dtyped (N,)-shaped :py:class:`ndarray <numpy.ndarray>`
        azimuthal angle in radians.

    Returns
    -------
    out : double-dtyped (N,)-shaped :py:class:`ndarray <numpy.ndarray>`
        Radiation intensity at the specified directions in units of 
        :math:`\mathrm{W}/\mathrm{sr}` (watt/steradian).

    """
    cdef:
        int lmax = <int> sqrt(q.shape[0] / 2 + 1) - 1
        Py_ssize_t index

        ndarray[double, ndim = 1] U = numpy.zeros_like(theta)
        ndarray[double complex, ndim = 2] E

    for index in range(theta.shape[0]):

        E = T(lmax, theta[index], phi[index]) @ q

        U[index] = (E.conjugate().T @ E)[0, 0].real / 2 / eta

    return U


cpdef ndarray[double, ndim=1] directivity(ndarray[double complex, ndim=2] q,
                                          ndarray[double, ndim=1] theta,
                                          ndarray[double, ndim=1] phi):
    r"""Calculates the directivity of an antenna.
    
    Parameters
    ----------
    q : complex-dtyped (N,)-shaped :py:class:`ndarray <numpy.ndarray>`
        Spherical mode coefficients (SMC) in :math:`W^{\frac{1}{2}} 
        (square root of watt)`.
    
    theta : double-dtyped (N,)-shaped :py:class:`ndarray <numpy.ndarray>`
        polar angle in radians.

    phi : double-dtyped (N,)-shaped :py:class:`ndarray <numpy.ndarray>`
        azimuthal angle in radians.
    
    Returns
    -------
    out : double-dtyped (N,)-shaped :py:class:`ndarray <numpy.ndarray>`
        Directivity at the specified directions (dimensionless).

    """
    cdef:
        int lmax = <int> sqrt(q.shape[0] / 2 + 1) - 1
        Py_ssize_t index

        double qHq = (q.conjugate().T @ q)[0, 0].real

        ndarray[double, ndim=1] D = numpy.zeros_like(theta)
        ndarray[double complex, ndim=2] E

    for index in range(theta.shape[0]):
        E = T(lmax, theta[index], phi[index]) @ q

        D[index] = (E.conjugate().T @ E)[0, 0].real

    return 4 * pi / eta * D / qHq


class SMC:
    """
    Spherical Mode Coefficients

    Parameters
    ----------
    q : complex-dtyped (p,) or (p, 1)-shaped :py:class:`ndarray <numpy.ndarray>`
        Spherical mode coefficients values in :math:`W^{\frac{1}{2}}
        (square root of watt)`.

    freq : :py:class:`double <float>`
        Frequency in GHz
    """
    # ========== ========== ========== ========== ========== class attributes
    ...

    # ========== ========== ========== ========== ========== special methods
    def __init__(self, ndarray[double complex, ndim=1] q, double freq):
        ...

    # ========== ========== ========== ========== ========== private methods
    ...

    # ========== ========== ========== ========== ========== protected methods
    ...

    # ========== ========== ========== ========== ========== public methods
    def radiation_field(self,
                        ndarray[double, ndim=1] theta,
                        ndarray[double, ndim=1] phi) -> NDArray:
        ...

    # ---------- ---------- ---------- ---------- ---------- properties
    ...




# ========== ========== ========== ========== ========== ==========
class SMCEstimator(metaclass=ABCMeta):


    
    # ========== ========== ========== ========== ========== class attributes
    ...

    # ========== ========== ========== ========== ========== special methods
    def __new__(cls, *args, **kwargs):

        cdef str method

        if cls is SMCEstimator:

            method = kwargs.pop('method', 'blue').lower()

            if method == 'blue':
                return super().__new__(SMCBLUE, *args, **kwargs)

            elif method == 'ridge':
                return super().__new__(SMCRidge, *args, **kwargs)

            else:
                raise NotImplementedError()

        return super().__new__(cls, *args, **kwargs)

    # ========== ========== ========== ========== ========== private methods
    ...

    # ========== ========== ========== ========== ========== protected methods
    ...

    # ========== ========== ========== ========== ========== public methods
    def show_dashboard(self, *args, **kwargs) -> None:
        raise NotImplementedError()

    def save_report(self, *args, **kwargs) -> None:
        raise NotImplementedError()

    def save(self, path: Path | str, overwrite: bool = False) -> None:

        cdef:
            object path = Path(path)
            object file

        if not path.suffix:
            path = path.with_suffix('.smce')

        if path.is_file() and not overwrite:
            raise FileExistsError()

        file = h5py.File(path, 'x')

        file.attrs['description'] = self.description
        file.attrs['residual_energy'] = self.residual_energy

    # ---------- ---------- ---------- properties (regression related)
    @property
    def description(self) -> str | None:
        try:
            return self.__description

        except AttributeError:
            self.__description = None

        return self.__description

    @description.setter
    def description(self, value: str | None) -> None:
        self.__description = value

    # ---------- ---------- ---------- properties (regression related)
    @property
    @abstractmethod
    def residual(self) -> NDArray:
        raise NotImplementedError()

    @property
    @abstractmethod
    def residual_energy(self) -> float:
        raise NotImplementedError()

    @property
    @abstractmethod
    def param_cov_matrix(self) -> NDArray:
        raise NotImplementedError()

    @property
    def P(self) -> NDArray:
        return self.param_cov_matrix

    @property
    @abstractmethod
    def estimation_matrix(self) -> NDArray:
        raise NotImplementedError()

    @property
    def K(self) -> NDArray:
        return self.estimation_matrix

    @property
    @abstractmethod
    def regressor_matrix(self) -> NDArray:
        raise NotImplementedError()

    @property
    def Upsilon(self) -> NDArray:
        return self.regressor_matrix

    @property
    @abstractmethod
    def smc_matrix(self) -> NDArray:
        raise NotImplementedError()

    @property
    def q(self) -> NDArray:
        return self.smc_matrix

    @property
    @abstractmethod
    def bias_matrix(self) -> NDArray:
        raise NotImplementedError()

    @property
    def b(self):
        return self.bias_matrix

    @property
    @abstractmethod
    def condition_number(self) -> NDArray:
        raise NotImplementedError()

    @property
    @abstractmethod
    def lagrange_multiplier_matrix(self) -> NDArray:
        raise NotImplementedError()

    @property
    def Lambda(self) -> NDArray:
        return self.lagrange_multiplier_matrix






class SMCBLUE(SMCEstimator):

    def __init__(self, *args, **kwargs):
        print('Running Blue')


class SMCRidge(SMCEstimator):

    def __init__(self, *args, **kwargs):
        print('Running Ridge')