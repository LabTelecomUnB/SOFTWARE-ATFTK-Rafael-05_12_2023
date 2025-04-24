#  -*- coding: utf-8 -*-
"""

Author: Rafael R. L. Benevides
Date: 4/24/23

"""


cimport numpy as numpyc

import numpy

from libc.math cimport sqrt


from numpy.typing import ArrayLike


cdef class ModeCoeff:
    """

    """

    cdef:
        readonly int k, Nk

        double complex[:, :] _q

    # ========== ========== ========== ========== ========== class attributes
    ...

    # ========== ========== ========== ========== ========== special methods
    def __cinit__(self, double complex[:, :] q) -> None:

        self._q = q

        assert self._q.shape[0] % 2 == 0
        assert self._q.shape[1] == 1

        self.Nk = <int>(self._q.shape[0] // 2)
        self.k = <int>(sqrt(self.Nk + 1) - 1)

    # ========== ========== ========== ========== ========== private methods
    ...

    # ========== ========== ========== ========== ========== protected methods
    ...

    # ========== ========== ========== ========== ========== public methods
    def plot_power_distribution(self) -> None:
        pass

    # ---------- ---------- ---------- ---------- ---------- properties
    @property
    def q(self) -> ArrayLike:
        return numpy.asarray(self._q)

    @property
    def qTE(self) -> ArrayLike:
        return self.q.reshape(-1, 2)[:, 0].reshape(-1, 1)

    @property
    def qTM(self) -> ArrayLike:
        return self.q.reshape(-1, 2)[:, 1].reshape(-1, 1)

    @property
    def n_modes(self) -> int:
        return self.Nk

    @property
    def highest_degree(self) -> int:
        return self.k