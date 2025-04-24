#  -*- coding: utf-8 -*-
"""

Author: Rafael R. L. Benevides
Date: 8/25/23

"""

# ---------- ---------- ---------- ---------- ---------- ---------- py
import numpy
import pyvista
from matplotlib import pyplot, rc_context

from aftk.mathtools.sptri import SphericalMesh

# ---------- ---------- ---------- ---------- ---------- ---------- cy
from numpy cimport ndarray

from libc.math cimport INFINITY, NAN

from aftk.antenna cimport directivity, radiation_intensity

# ---------- ---------- ---------- ---------- ---------- ---------- ty



# ========== ========== ========== ========== ========== Radiation Pattern
cdef ndarray[double] to_dB(ndarray[double] value, double dB_floor=-INFINITY):
    """
    
    Parameters
    ----------
    value
    floor

    Returns
    -------

    """
    cdef:
        ndarray[double] value_dB = numpy.zeros_like(value)
        ndarray[double] value_dB_aux = 10.0 * numpy.log10(value[value > 0.0])

    value_dB_aux[value_dB_aux < dB_floor] = dB_floor

    value_dB[value < 0.0] = NAN
    value_dB[value == 0.0] = dB_floor
    value_dB[value > 0.0] = value_dB_aux

    return value_dB


def plot_radiation_pattern(ndarray[double complex, ndim=2] q,
                           int sm_resolution=15,
                           int sm_refinement_levels=3,
                           str quantity='dir_dBi', # or dB, W/sr, dBW, dimensionless,
                           double dB_floor = -60.0,
                           double dB_range = 30.0):
    """

    Parameters
    ----------
    q : complex-dtyped (N,)-shaped :py:class:`ndarray <numpy.ndarray>`
        Spherical mode coefficients (SMC) in :math:`W^{\frac{1}{2}}
        (square root of watt)`.

    sm_resolution : :py:class:`int`
        TODO

    sm_refinement_levels : :py:class:`int`
        TODO

    quantity : :py:class:`str`
        The plotted quantity. Acceptable values: `'dir_dBi'`, `'dBW'`, `'dB'`

    dB_floor : :py:class:`double <float>`
        TODO

    dB_range : :py:class:`double <float>`
        TODO

    Returns
    -------
    values :
        TODO

    radiation_pattern :
        TODO

    """
    cdef:

        object sph_mesh = SphericalMesh(resolution=sm_resolution,
                                        refinement_levels=sm_refinement_levels)

        ndarray[double, ndim=2] r = sph_mesh.points, _value
        ndarray[double, ndim=1] value

        str title

    print(len(sph_mesh))

    # ---------- ---------- ---------- ---------- creating value
    if quantity == 'dir':
        value = directivity(q, sph_mesh.theta, sph_mesh.phi)
        title = 'Directivity (dimensionless)'

    elif quantity == 'dir_dBi':
        value = directivity(q, sph_mesh.theta, sph_mesh.phi)
        value = to_dB(value, dB_floor=dB_floor)
        title = 'Directivity (dBi)'

    elif quantity == 'rad_int':
        value = radiation_intensity(q, sph_mesh.theta, sph_mesh.phi)
        title = 'Radiation Intensity (W/sr)'

    elif quantity == 'rad_int_norm':
        value = radiation_intensity(q, sph_mesh.theta, sph_mesh.phi)
        value = value / value.max()
        title = 'Normalised Radiation Intensity (dimensionless)'

    elif quantity == 'rad_int_dBW':
        value = radiation_intensity(q, sph_mesh.theta, sph_mesh.phi)
        value = to_dB(value, dB_floor=dB_floor)
        title = 'Radiation Intensity (dBW)'

    elif quantity == 'rad_int_dBm':
        value = radiation_intensity(q, sph_mesh.theta, sph_mesh.phi)
        value = to_dB(value*1000, dB_floor=dB_floor)
        title = 'Radiation Intensity (dBm)'

    elif quantity == 'rad_int_norm_dB':
        value = radiation_intensity(q, sph_mesh.theta, sph_mesh.phi)
        value = to_dB(value / value.max(), dB_floor=dB_floor)
        title = 'Normalised Radiation Intensity (dB)'

    else:
        raise ValueError(f'{quantity=} is not valid.')

    # ---------- ---------- ---------- ---------- plotting
    _value = r * (value - value.min()).reshape(-1, 1)

    plotter = pyvista.Plotter()
    plotter.background_color = 'gray'

    radiation_pattern = sph_mesh.unstructured_grid(_value)

    if 'dB' in quantity:
        value[value <= value.max() - dB_range] = value.max() - dB_range

    plotter.add_mesh(radiation_pattern,
                     cmap='jet',
                     smooth_shading=True,
                     scalars=value,
                     scalar_bar_args = {'title': title},
                     show_scalar_bar=True,
                     name='radiation_pattern',
                     # opacity=0.5
                     )
    # plotter.show_grid()
    plotter.show_axes()
    plotter.show()

    return value, radiation_pattern


def plot_radiation_pattern_sections(ndarray[double complex, ndim=2] q,
                                    int n_xy=360,
                                    int n_xz=360,
                                    int n_yz=360,
                                    double dB_floor = -60.0,):

    cdef:

        object figure, axes

        ndarray[double, ndim=1] theta, phi, _theta

        ndarray[double, ndim=1] value_xy, value_xz, value_yz


    figure, axes = pyplot.subplots(1, 3, subplot_kw={'projection': 'polar'})

    # ---------- ---------- ---------- ---------- ---------- xy plane
    phi = numpy.linspace(0, 2*numpy.pi, n_xy)
    theta = numpy.full_like(phi, numpy.pi/2)

    value_xy = to_dB(directivity(q, theta, phi), dB_floor)

    axes[0].set_title('XY - Plane')
    axes[0].set_rmax(1.05*value_xy.max())
    axes[0].plot(phi, value_xy)

    # ---------- ---------- ---------- ---------- ---------- xz plane
    theta = numpy.linspace(0, numpy.pi, n_xz//2)
    phi = numpy.zeros_like(theta)
    _theta = numpy.hstack([theta, numpy.pi + theta[1:]])

    theta = numpy.hstack([theta, theta[-2::-1]])
    phi = numpy.hstack([phi, numpy.full_like(phi, numpy.pi)[1:]])

    value_xz = to_dB(directivity(q, theta, phi), dB_floor)
    value_xz -= value_xz.min()

    axes[1].set_title('XZ - Plane')
    axes[1].set_theta_zero_location("N")
    # axes[1].set_rmin(dB_floor)
    axes[1].set_rmax(1.05 * value_xz.max())
    axes[1].plot(_theta, value_xz)

    # ---------- ---------- ---------- ---------- ---------- yz plane
    theta = numpy.linspace(0, numpy.pi, n_yz//2)
    phi = numpy.zeros_like(theta) + numpy.pi/2
    _theta = numpy.hstack([theta, numpy.pi + theta[1:]])

    theta = numpy.hstack([theta, theta[-2::-1]])
    phi = numpy.hstack([phi, numpy.full_like(phi, 1.5*numpy.pi)[1:]])

    value_yz = to_dB(directivity(q, theta, phi), dB_floor)
    value_yz -= value_yz.min()

    axes[2].set_title('YZ - Plane')
    axes[2].set_theta_zero_location("N")
    # axes[2].set_rmin(dB_floor) ;
    axes[2].set_rmax(1.05 * value_yz.max())
    axes[2].plot(_theta, value_yz)



    pyplot.show()







