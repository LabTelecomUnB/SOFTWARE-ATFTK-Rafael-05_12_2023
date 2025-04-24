#  -*- coding: utf-8 -*-
"""
Author: Rafael R. L. Benevides
"""
import numpy

from aftk.utils.misc import plot_matrix, condition_number_of_PDM
from aftk.model.eigenfunctions import Upsilonk
from aftk.measurements.radiation_field import plot_regressors

from matplotlib import rc


rc('text', usetex=True)
rc('font', family='serif')


# ========== ========== ========== ========== ========== ========== misc
def test_plot_matrix():

    theta = numpy.deg2rad(numpy.arange(0, 180, 30))[1:]
    phi = numpy.deg2rad(numpy.arange(0, 360, 30))

    theta, phi = numpy.meshgrid(theta, phi)
    theta = theta.flatten()
    phi = phi.flatten()

    U = Upsilonk(73, theta, phi, 1)

    plot_regressors(U, show=True)


def test_condition_number_of_PDM():
    random_generator = numpy.random.default_rng(seed=0)

    A = random_generator.random((3, 3)) + 1j*random_generator.random((3, 3))
    A += A.conjugate().T

    cond_A_est = condition_number_of_PDM(A)
    cond_A = numpy.linalg.cond(A)

    assert numpy.isclose(cond_A_est, cond_A)