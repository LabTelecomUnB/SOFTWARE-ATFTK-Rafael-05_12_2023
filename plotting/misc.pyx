#  -*- coding: utf-8 -*-
"""
Author: Rafael R. L. Benevides
"""

# ---------- ---------- ---------- ---------- ---------- ---------- py
import numpy
import pyvista

# ---------- ---------- ---------- ---------- ---------- ---------- cy
from numpy cimport ndarray


def plot_spherical_data(ndarray[double, ndim=1] theta,
                        ndarray[double, ndim=1] phi,
                        ndarray[double, ndim=1] data, **kwargs) -> None:
    """
    Plot data in spherical coordinates.

    Parameters
    ----------
    theta : double-dtyped (n,)-shaped :py:class:`ndarray <numpy.ndarray>`
        Polar angle in radians

    phi : double-dtyped (n,)-shaped :py:class:`ndarray <numpy.ndarray>`
        Azimuthal angle in radians

    data : double-dtyped (n,)-shaped :py:class:`ndarray <numpy.ndarray>`
        The data to be plotted. Must be non-negative. Negative values will be
        considered zero.

    kwargs : :py:class:`dict`
        Additional keyword arguments to be passed to `plotter.add_mesh`.
    """
    cdef:
        ndarray[double, ndim=1] sin_theta = numpy.sin(theta)
        ndarray[double, ndim=1] _data = data.copy()
        ndarray[double, ndim=1] x, y, z

        object plotter

    _data[data < 0.0] = 0.0

    x = _data * numpy.cos(phi) * sin_theta
    y = _data * numpy.sin(phi) * sin_theta
    z = _data * numpy.cos(theta)

    plotter = pyvista.Plotter()
    plotter.background_color = 'gray'

    plotter.add_mesh(pyvista.StructuredGrid(x, y, z),
                     cmap='jet',
                     smooth_shading=True,
                     scalars=_data,
                     show_scalar_bar=True,
                     **kwargs)

    plotter.show_axes()
    plotter.show()
