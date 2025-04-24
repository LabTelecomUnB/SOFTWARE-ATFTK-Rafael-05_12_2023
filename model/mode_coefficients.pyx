#  -*- coding: utf-8 -*-
"""
Author: Rafael R. L. Benevides
"""

import h5py

cimport numpy as numpyc
import numpy
import pandas
import pyvista

from matplotlib import pyplot, ticker, rc

from numpy.typing import ArrayLike


cdef class RadiationFieldData:
    """

    Parameters
    ----------
    theta : ndarray
        Spherical polar angle (radians)

    phi : ndarray
        Spherical azimuthal angle (radians)

    E_theta : ndarray
        Radiation field polar component (volts)

    E_phi : ndarray
        Radiation field azimuthal component (volts)

    **kwargs : ...


    """

    # ========== ========== ========== ========== ========== class attributes
    cdef:
        double[:] _theta, _phi
        double complex[:] _E_theta, _E_phi
        readonly int M
        readonly double freq
        readonly str name, description

    # ========== ========== ========== ========== ========== special methods
    def __cinit__(self, double[:] theta,
                        double[:] phi,
                        double complex[:] E_theta,
                        double complex[:] E_phi,
                        double freq = None,
                        str name = None,
                        str description = None):

        # ---------- ---------- ---------- parse lengths
        assert theta.shape[0] == phi.shape[0]
        assert theta.shape[0] == E_theta.shape[0]
        assert theta.shape[0] == E_phi.shape[0]

        self.M = theta.shape[0]

        # ---------- ---------- ----------
        self._theta = theta
        self._phi = phi
        self._E_theta = E_theta
        self._E_phi = E_phi

        self.freq = freq
        self.name = name
        self.description = description

    # ========== ========== ========== ========== ========== private methods
    ...

    # ========== ========== ========== ========== ========== protected methods
    ...

    # ========== ========== ========== ========== ========== public methods
    def to_dataframe(self) -> pandas.DataFrame:

        return pandas.DataFrame({
            'theta': self.theta,
            'phi': self.phi,
            'E_theta': self.E_theta,
            'E_phi': self.E_phi
        })

    def plot_data_directions(self) -> None:

        cdef numpyc.ndarray[double, ndim=2] theta, phi, x, y, z, r
        cdef object plotter

        theta = self.theta.reshape(-1, 1)
        phi = self.phi.reshape(-1, 1)
        x = numpy.sin(theta) * numpy.cos(phi)
        y = numpy.sin(theta) * numpy.sin(phi)
        z = numpy.cos(theta)
        r = numpy.hstack([x, y, z])

        plotter = pyvista.Plotter()
        plotter.add_points(r, render_points_as_spheres=True)
        plotter.add_axes_at_origin(xlabel='x', ylabel='y', zlabel='z', line_width=10)
        plotter.show()

    def plot_data_magnitude(self,
                            double dB_floor=None,
                            int levels=50,
                            str cmap='jet',
                            bint show=False) -> None:

        cdef:
            numpyc.ndarray[numpyc.complex128_t, ndim=1] E_theta, E_phi
            numpyc.ndarray[numpyc.float64_t, ndim=1] E_theta_mag_sq, E_phi_mag_sq
            numpyc.ndarray[numpyc.float64_t, ndim=1] E_theta_mag_dB, E_phi_mag_dB
            numpyc.float64_t E_max_sq, _dB_floor

            object figure, axes, contour, colour_axes, colour_bar

        E_theta_mag_sq = E_theta.conjugate() * E_theta
        E_phi_mag_sq = E_phi.conjugate() * E_phi

        E_max_sq = numpy.max(E_theta_mag_sq + E_phi_mag_sq)

        E_theta_mag_dB = 10 * numpy.log10(E_theta_mag_sq / E_max_sq)
        E_phi_mag_dB = 10 * numpy.log10(E_phi_mag_sq / E_max_sq)

        if dB_floor is None:
            _dB_floor = min(numpy.min(E_theta_mag_dB), numpy.min(E_phi_mag_dB))
        else:
            _dB_floor = dB_floor

        # ---------- ---------- ---------- ---------- ----------
        figure, axes = pyplot.subplots(2, sharex=True, sharey=True)
        figure.set_size_inches(12, 5)
        # figure.suptitle(f'', y=0.93, fontsize=20)

        contour = axes[0].tricontourf(
            numpy.rad2deg(self.phi), numpy.rad2deg(self.theta), E_theta_mag_dB,
            levels=levels, cmap=cmap, vmax=0, vmin=_dB_floor)

        axes[0].set_xlabel(r'$\phi$ (deg)', fontsize=16)
        axes[0].set_ylabel(r'$\theta$ (deg)', fontsize=16)
        axes[0].set_title(r'$\|E_\theta\|$ (dB)', fontsize=18)
        axes[0].set_aspect('equal')

        axes[1].tricontourf(
            numpy.rad2deg(self.phi), numpy.rad2deg(self.theta), E_phi_mag_dB,
            levels=levels, cmap=cmap, vmax=0, vmin=_dB_floor)

        axes[1].set_xlabel(r'$\phi$ (deg)', fontsize=16)
        axes[1].set_ylabel(r'$\theta$ (deg)', fontsize=16)
        axes[1].set_title(r'$\|E_\phi\|$ (dB)', fontsize=18)
        axes[1].set_aspect('equal')

        colour_axes = figure.add_axes([])
        colour_bar.ax.set_ylabel('Radiation Field (dB)', fontsize=16)
        colour_bar.ax.tick_params(labelsize=12)
        colour_bar.ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

        if show:
            pyplot.show()

    # ========== ========== ========== ========== ========== properties
    @property
    def theta(self) -> ArrayLike:
        return numpy.asarray(self._theta)

    @property
    def phi(self) -> ArrayLike:
        return numpy.asarray(self._phi)

    @property
    def E_theta(self) -> ArrayLike:
        return numpy.asarray(self._E_theta)

    @property
    def E_phi(self) -> ArrayLike:
        return numpy.asarray(self._E_phi)

    @property
    def E(self) -> ArrayLike:
        return numpy.hstack([
            self.E_theta.reshape(-1, 1),
            self.E_phi.reshape(-1, 1)
        ]).reshape(-1, 1)