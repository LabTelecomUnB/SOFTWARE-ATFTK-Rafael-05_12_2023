#  -*- coding: utf-8 -*-
"""
Author: Rafael R. L. Benevides
"""

import pytest

import numpy

from aftk.mathtools import psdm, spharm


class TestPSDM:

    @pytest.fixture
    def A(self):
        shape = (5, 3)

        random_generator = numpy.random.default_rng(seed=0)

        return random_generator.random(shape) + 1j * random_generator.random(shape)

    def test_eigenvalues(self, A):
        eigvals_numpy = numpy.sort(numpy.linalg.eigvalsh(A.conjugate().T @ A))
        eigvals_aftk = numpy.sort(psdm.eigenvalues(A))

        assert numpy.allclose(eigvals_aftk, eigvals_numpy)

    def test_condition_number(self, A):

        cond_numpy = numpy.linalg.cond(A.conjugate().T @ A)
        cond_aftk = psdm.condition_number(A)

        assert numpy.isclose(cond_aftk, cond_numpy)

    def test_inverse(self, A):

        inv_numpy = numpy.linalg.inv(A.conjugate().T @ A)
        inv_aftk = psdm.inverse(A)

        assert numpy.allclose(inv_aftk, inv_numpy)


class TestSPHARM:

    def test_continuity_of_eigenfunctions_at_the_poles(self):

        phis = numpy.linspace(0, 2*numpy.pi, 10)[:-1]

        modes = [(l, m) for l in range(1, 6) for m in range(-l, l+1)]

        theta = 1.0e-7

        for mode in modes:
            for phi in phis:

                # ---------- ---------- ---------- ---------- upper pole
                value_at_pole = spharm.m_times_Ylm_over_sin_theta(*mode, 0.0, phi)
                value_vicinity = spharm.m_times_Ylm_over_sin_theta(*mode, theta, phi)

                assert numpy.isclose(value_vicinity, value_at_pole, atol=1e-3), \
                    f'Problem at {mode=} and {phi=}:\n{value_at_pole=}\n{value_vicinity=}'

                # ---------- ---------- ---------- ---------- lower pole
                value_at_pole = spharm.m_times_Ylm_over_sin_theta(*mode, numpy.pi, phi)
                value_vicinity = spharm.m_times_Ylm_over_sin_theta(*mode, numpy.pi-theta, phi)

                assert numpy.isclose(value_vicinity, value_at_pole, atol=1e-3), \
                    f'Problem at {mode=} and {phi=}:\n{value_at_pole=}\n{value_vicinity=}'

                # ---------- ---------- ---------- ---------- upper pole
                value_at_pole = spharm.dYlm_dtheta(*mode, 0.0, phi)
                value_vicinity = spharm.dYlm_dtheta(*mode, theta, phi)

                assert numpy.isclose(value_vicinity, value_at_pole, atol=1e-3), \
                    f'Problem at {mode=} and {phi=}:\n{value_at_pole=}\n{value_vicinity=}'

                # ---------- ---------- ---------- ---------- lower pole
                value_at_pole = spharm.dYlm_dtheta(*mode, numpy.pi, phi)
                value_vicinity = spharm.dYlm_dtheta(*mode, numpy.pi-theta, phi)

                assert numpy.isclose(value_vicinity, value_at_pole, atol=1e-3), \
                    f'Problem at {mode=} and {phi=}:\n{value_at_pole=}\n{value_vicinity=}'
