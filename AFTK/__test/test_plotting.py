#  -*- coding: utf-8 -*-
"""

Author: Rafael R. L. Benevides
Date: 8/25/23

"""

import pytest

import numpy

from aftk.antenna.eigenantennas import get_eigenantennas
from aftk.plotting.smc import plot_radiation_pattern, plot_radiation_pattern_sections


@pytest.fixture
def eigenantennas():
    return get_eigenantennas(20, 0.0, 0.0)


@pytest.fixture
def laurent_lmax_60():
    return numpy.load('bias_analysis_lmax60.npz')['q']


class Test_plot_radiation_pattern:

    def test_plotting_in_dBi(self, eigenantennas):

        q = eigenantennas[:, 1].reshape(-1, 1)

        plot_radiation_pattern(q,
                               sm_resolution=15,
                               sm_refinement_levels=4,
                               quantity='dir_dBi',
                               dB_floor=-10)


class Test_plot_radiation_patter_section:

    def test_plotting(self, eigenantennas):

        q = eigenantennas[:, 1].reshape(-1, 1)
        #
        # q = numpy.load('bias_analysis_lmax60.npz')['q']

        plot_radiation_pattern_sections(q, n_xz=360, dB_floor=-10)



