#  -*- coding: utf-8 -*-
"""

Author: Rafael R. L. Benevides
Date: 8/16/23

"""

# ---------- ---------- ---------- ---------- ---------- ----------
import numpy
import pyvista

from matplotlib import pyplot

# ---------- ---------- ---------- ---------- ---------- cy
cimport numpy as numpyc
from libc.math cimport sqrt, sin, cos, acos, pi, atan2

# from aftk.mathtools.spharm cimport T as Tk
# from aftk.model.eigenfunctions import Tk

from aftk.antenna cimport directivity


def plot_radiation_pattern_3D(numpyc.ndarray[numpyc.complex128_t, ndim=2] q,
                              int n_theta=60, int n_phi=60,
                              tuple theta_range=(0.0, numpy.pi),
                              tuple phi_range=(0.0, 2 * numpy.pi)):
    cdef:
        numpyc.ndarray[numpyc.float64_t, ndim=1] _theta, _phi
        numpyc.ndarray[numpyc.float64_t, ndim=2] theta, phi
        numpyc.ndarray[numpyc.float64_t, ndim=2] _directivity, directivity_aux
        numpyc.ndarray[numpyc.float64_t, ndim=2] directiviy_dB, directivity_dB_aux

        numpyc.ndarray[numpyc.float64_t, ndim=2] x, y, z

        object plotter, radiation_pattern

        Py_ssize_t i, j, shape_0, shape_1

    _theta = numpy.linspace(theta_range[0], theta_range[1], n_theta)
    _phi = numpy.linspace(phi_range[0], phi_range[1], n_phi)

    theta, phi = numpy.meshgrid(_theta, _phi)

    shape_0 = theta.shape[0]
    shape_1 = theta.shape[1]

    # ---------- ---------- ---------- ---------- ---------- ----------
    _directivity = directivity(q, theta.flatten(), phi.flatten()).reshape(shape_0, shape_1)

    print(f'max directivity = {10 * numpy.log10(_directivity.max())} dBi')

    # for i in range(theta.shape[0]):
    #     for j in range(theta.shape[1]):
    #         _directivity[i, j] = directivity(q, theta[i, j], phi[i, j])

    directivity_aux = numpy.copy(_directivity)

    directivity_aux[directivity_aux == 0.0] = 0.1 * directivity_aux[directivity_aux != 0.0].min()

    directivity_dB = 10.0 * numpy.log10(directivity_aux)

    directivity_dB_aux = directivity_dB - directivity_dB.min()

    x = directivity_dB_aux * numpy.cos(phi) * numpy.sin(theta)
    y = directivity_dB_aux * numpy.sin(phi) * numpy.sin(theta)
    z = directivity_dB_aux * numpy.cos(theta)

    plotter = pyvista.Plotter()
    plotter.background_color = 'gray'

    radiation_pattern = pyvista.StructuredGrid(x, y, z)

    # directivity_dB[directivity_dB <= directivity_dB.max() - 30] = directivity_dB.max() - 30

    plotter.add_mesh(radiation_pattern,
                     cmap='jet',
                     smooth_shading=True,
                     # scalars=numpy.linspace(directivity_dB.max() - 30, directivity_dB.max(), directivity_dB.T.flatten().shape[0]),
                     scalars=directivity_dB.T.flatten(),
                     # scalars="Directivity Gain (dBi)",
                     show_scalar_bar=True,
                     name='radiation_pattern',
                     # opacity=0.5
                     )
    # plotter.show_grid()
    plotter.show_axes()
    plotter.show()

    return _directivity


def plot_directive_gain(numpyc.ndarray[numpyc.complex128_t, ndim=2] q):

    cdef:
        int k = <int> sqrt(q.shape[0] / 2 + 1) - 1
        double eta = 376.730_313_668  # ohm

    # ---------- ---------- ---------- ---------- ----------

    # ---------- ---------- ---------- ---------- ---------- plotting
    figure, axes = pyplot.subplots(1, 2, subplot_kw={'projection': 'polar'})




