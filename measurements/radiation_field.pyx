#  -*- coding: utf-8 -*-
"""
Author: Rafael R. L. Benevides
"""

# ---------- ---------- ---------- ---------- ---------- py
import h5py

from pathlib import Path

import numpy
import pandas
import pyvista
import scipy.linalg
from matplotlib import pyplot, ticker, rc, patches

from datetime import datetime

# ---------- ---------- ---------- ---------- ---------- cy
cimport numpy as numpyc
from libc.math cimport sqrt, sin, cos, acos, pi, atan2

from numpy cimport ndarray, complex128_t, float64_t

from aftk.model cimport Upsilonk
from aftk.model.eigenfunctions import Tk
from aftk.utils.misc import pdm_inverse

# ---------- ---------- ---------- ---------- ---------- typing
from numpy.typing import ArrayLike


cpdef ndarray[float64_t, ndim=2] rot_x(float64_t ang):
    """
    
    Parameters
    ----------
    ang

    Returns
    -------

    """
    cdef:
        double sin_ang = sin(ang), cos_ang = cos(ang)
        double[3][3] rot

    rot[0][0] = 1.0
    rot[0][1] = 0.0
    rot[0][2] = 0.0

    rot[1][0] = 0.0
    rot[1][1] = cos_ang
    rot[1][2] = sin_ang

    rot[2][0] = 0.0
    rot[2][1] = -sin_ang
    rot[2][2] = cos_ang

    return numpy.asarray(rot)


cpdef ndarray[float64_t, ndim=2] rot_y(float64_t ang):
    """

    Parameters
    ----------
    ang

    Returns
    -------

    """
    cdef:
        double sin_ang = sin(ang), cos_ang = cos(ang)
        double[3][3] rot

    rot[0][0] = cos_ang
    rot[0][1] = 0.0
    rot[0][2] = -sin_ang

    rot[1][0] = 0.0
    rot[1][1] = 1.0
    rot[1][2] = 0.0

    rot[2][0] = sin_ang
    rot[2][1] = 0.0
    rot[2][2] = cos_ang

    return numpy.asarray(rot)


cpdef ndarray[float64_t, ndim=2] rot_z(float64_t ang):
    """

    Parameters
    ----------
    ang

    Returns
    -------

    """
    cdef:
        double sin_ang = sin(ang), cos_ang = cos(ang)
        double[3][3] rot

    rot[0][0] = cos_ang
    rot[0][1] = sin_ang
    rot[0][2] = 0.0

    rot[1][0] = -sin_ang
    rot[1][1] = cos_ang
    rot[1][2] = 0.0

    rot[2][0] = 0.0
    rot[2][1] = 0.0
    rot[2][2] = 1.0

    return numpy.asarray(rot)


def inv_hermitian_matrix(numpyc.ndarray[numpyc.complex128_t, ndim = 2] A):

    cdef:
        numpyc.ndarray[numpyc.complex128_t, ndim = 2] L
        numpyc.ndarray[numpyc.complex128_t, ndim = 2] eye = numpy.eye(A.shape[0], dtype=numpy.complex128)
        numpyc.ndarray[numpyc.complex128_t, ndim = 2] y
        # numpyc.ndarray[numpyc.complex128_t, ndim = 2] x = numpy.zeros_like(A, dtype=numpy.complex128)

    try:
        L = numpy.linalg.cholesky(A)
        y = numpy.zeros_like(A, dtype=numpy.complex128)
        y = scipy.linalg.solve_triangular(L, eye, lower=True)

        return scipy.linalg.solve_triangular(L.conj().T, y, lower=False)

    except numpy.linalg.LinAlgError:
        return scipy.linalg.solve(A, eye, check_finite=False, assume_a='her')


cdef ndarray[float64_t, ndim=2] rot_matrix_sph_from_cart(float64_t theta, float64_t phi):
    """
    
    Parameters
    ----------
    theta
    phi

    Returns
    -------

    """

    cdef:
        double sin_theta = sin(theta), cos_theta = cos(theta)
        double sin_phi = sin(phi), cos_phi = cos(phi)

        double[3][3] rot

    rot[0][0] = sin_theta * cos_phi
    rot[0][1] = sin_theta * sin_phi
    rot[0][2] = cos_theta

    rot[1][0] = cos_theta * cos_phi
    rot[1][1] = cos_theta * sin_phi
    rot[1][2] = -sin_theta

    rot[2][0] = -sin_phi
    rot[2][1] = cos_phi
    rot[2][2] = 0.0

    return numpy.asarray(rot)


def plot_regressors(U: ArrayLike, show=False):

    figure, axes = pyplot.subplots(2, 1, sharex=True, sharey=True)
    # figure.set_size_inches(8, 12)

    A_mag = numpy.abs(U)
    A_pha = numpy.angle(U, deg=True)

    mappable_mag = axes[0].imshow(A_mag, cmap='jet', origin='upper', interpolation='nearest')
    axes[0].set_xlabel(r'Regressors degree ($\ell$) and order ($m$)', fontsize=24)
    axes[0].set_ylabel(r'Samples', fontsize=24)
    axes[0].set_title(r'Regressors Magnitude', fontsize=26)
    axes[0].set_aspect('auto')

    colourbar = figure.colorbar(mappable_mag, ax=axes[0])
    colourbar.set_label(f'Regressors Magnitude', fontsize=20)
    colourbar.ax.tick_params(labelsize=18)

    mappable_pha = axes[1].imshow(A_pha, cmap='jet', origin='upper', interpolation='nearest')
    axes[1].set_xlabel(r'Regressors degree ($\ell$) and order ($m$)', fontsize=24)
    axes[1].set_ylabel(r'Samples', fontsize=24)
    axes[1].set_title(r'Regressors Phase', fontsize=26)
    axes[1].set_aspect('auto')

    colourbar = figure.colorbar(mappable_pha, ax=axes[1])
    colourbar.set_label(f'Regressors Phase', fontsize=20)
    colourbar.ax.tick_params(labelsize=18)

    if show:
        pyplot.show()


def plot_covariance(U: ArrayLike, show=False):

    cdef:
        Py_ssize_t M = U.shape[0] // 2
        Py_ssize_t Nk = U.shape[1] // 2, N_k
        Py_ssize_t l_max = <int>(numpy.sqrt(Nk + 1) - 1)

        ndarray[Py_ssize_t, ndim=1] k = numpy.arange(20, l_max+1)
        ndarray[numpyc.float64_t, ndim=1] eigenvalues
        ndarray[complex128_t, ndim=2] Uk
        double trace

        object trace_P

    if l_max <= 20:
        raise ValueError()

    trace_P = []

    for _k in k:

        N_k = _k * (_k + 2)
        Uk = U[:, :2*N_k]
        eigenvalues = numpy.linalg.eigvalsh(Uk.conjugate().T @ Uk)
        trace = (1 / eigenvalues).sum()
        trace_P.append(trace)

        print(f'{datetime.now()}\t{_k}\t{trace}')

    # ---------- ---------- ---------- ---------- ---------- ----------
    figure, axes = pyplot.subplots(2, 1, sharex=True, sharey=True)

    axes[0].plot(k, numpy.array(trace_P), '-o')

    mappable_mag = axes[1].imshow(numpy.abs(U), cmap='jet', origin='upper', interpolation='nearest')
    axes[1].set_xlabel(r'Regressors degree ($\ell$) and order ($m$)', fontsize=16)
    axes[1].set_ylabel(r'Samples', fontsize=16)
    axes[1].set_title(r'Regressors Magnitude', fontsize=18)
    axes[1].set_aspect('auto')

    colourbar = figure.colorbar(mappable_mag, ax=axes[1])
    colourbar.set_label(f'Regressors Magnitude', fontsize=12)
    colourbar.ax.tick_params(labelsize=10)

    if show:
        pyplot.show()


def get_gershgorin_discs(numpyc.ndarray[numpyc.complex128_t, ndim=2] A: ArrayLike):

    cdef:
        Py_ssize_t i, j
        double _radius

        double complex[:] centre = numpy.zeros(A.shape[0], dtype=numpy.complex_)
        double[:] radius = numpy.zeros(A.shape[0], dtype=numpy.double)

    for i in range(A.shape[0]):

        centre[i] = A[i, i]

        _radius = 0.0

        for j in range(A.shape[0]):
            if j != i:
                _radius += numpy.abs(A[i, j])

        radius[i] = _radius

    return numpy.asarray(centre), numpy.asarray(radius)


def plot_discs(numpyc.ndarray[numpyc.complex128_t, ndim=1] centre,
               numpyc.ndarray[numpyc.float64_t, ndim=1] radius,
               **kwargs):

    cdef:
        Py_ssize_t index

        object figure, axes, circle

        numpyc.ndarray[numpyc.float64_t, ndim=1] x, y

    x = centre.real
    y = centre.imag

    figure, axes = pyplot.subplots()

    for index in range(centre.shape[0]):

        circle = patches.Circle(
            (x[index], y[index]),
            radius=radius[index],
            facecolor='C0',
            edgecolor='C0',
            alpha=0.3,
            **kwargs)

        axes.add_patch(circle)

    axes.plot(x, y, marker='o', linestyle='', color='C3')

    axes.set_aspect('equal')
    axes.grid()

    return figure, axes





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
        double eta
        readonly numpyc.npy_bool[:] no_pole_cond
        double[:] _theta, _phi
        double complex[:] _E_theta, _E_phi
        double complex[:, :] _regressors, EEH_inv

        readonly int M
        readonly double frequency
        readonly str name, description


    # ========== ========== ========== ========== ========== special methods
    def __cinit__(self, double[:] theta,
                        double[:] phi,
                        double complex[:] E_theta,
                        double complex[:] E_phi,
                        frequency = None,
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

        if frequency is not None:
            self.frequency = frequency

        self.name = name
        self.description = description

        self.eta = 376.730_313_668  # ohm

        self.no_pole_cond = (self.theta > 0) & (self.theta < numpy.pi)

        self.EEH_inv = None

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

    def plot_measured_directions(self) -> None:

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

    def plot_2D_measured_radiation_pattern(self,
                                           dB_floor=None,
                                           int levels=50,
                                           str cmap='jet',
                                           bint show=False) -> None:

        cdef:
            numpyc.ndarray[numpyc.complex128_t, ndim=1] E_theta, E_phi
            numpyc.ndarray[numpyc.float64_t, ndim=1] E_theta_mag_sq, E_phi_mag_sq
            numpyc.ndarray[numpyc.float64_t, ndim=1] E_theta_mag_dB, E_phi_mag_dB
            numpyc.float64_t E_max_sq, _dB_floor

            object figure, axes, contour, colour_axes, colour_bar

        E_theta = self.E_theta
        E_phi = self.E_phi

        E_theta_mag_sq = numpy.real(E_theta.conjugate() * E_theta)
        E_phi_mag_sq = numpy.real(E_phi.conjugate() * E_phi)

        E_max_sq = numpy.max(E_theta_mag_sq + E_phi_mag_sq)

        E_theta_mag_dB = 10 * numpy.log10(E_theta_mag_sq / E_max_sq)
        E_phi_mag_dB = 10 * numpy.log10(E_phi_mag_sq / E_max_sq)

        if dB_floor is None:
            _dB_floor = min(numpy.min(E_theta_mag_dB), numpy.min(E_phi_mag_dB))
        else:
            _dB_floor = dB_floor

        # ---------- ---------- ---------- ---------- ----------
        rc('text', usetex=True)
        rc('font', family='serif')

        figure, axes = pyplot.subplots(2, sharex=True, sharey=True)
        figure.set_size_inches(8.5, 9.5)
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

        colour_axes = figure.add_axes([
            axes[1].get_position().x1 + 0.015,
            axes[1].get_position().y0,
            0.02,
            axes[0].get_position().y1 - axes[1].get_position().y0
        ])
        colour_bar = figure.colorbar(contour, cax=colour_axes)
        colour_bar.ax.set_ylabel('Radiation Field (dB)', fontsize=16)
        colour_bar.ax.tick_params(labelsize=12)
        colour_bar.ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

        if show:
            pyplot.show()

    def plot_3D_measured_radiation_pattern(self) -> None:

        cdef:
            numpyc.ndarray[double, ndim=1] theta, phi, x, y, z
            numpyc.ndarray[double, ndim=2] r
            numpyc.ndarray[numpyc.complex128_t, ndim=1] E_theta, E_phi
            numpyc.ndarray[numpyc.float64_t, ndim=1] E_mag_sq, E_mag_dB
            numpyc.float64_t E_mag_sq_max, _dB_floor

            object plotter, radiation_pattern

        E_theta = self.E_theta
        E_phi = self.E_phi

        E_mag_sq = numpy.real(E_theta.conjugate() * E_theta + E_phi.conjugate() * E_phi)
        E_mag_sq_max = numpy.max(E_mag_sq)
        E_mag_dB =  10 * numpy.log10(E_mag_sq / E_mag_sq_max)
        E_mag_dB -= numpy.min(E_mag_dB)

        theta = self.theta
        phi = self.phi
        x = E_mag_dB * numpy.sin(theta) * numpy.cos(phi)
        y = E_mag_dB * numpy.sin(theta) * numpy.sin(phi)
        z = E_mag_dB * numpy.cos(theta)

        r = numpy.hstack([
            x.reshape(-1, 1),
            y.reshape(-1, 1),
            z.reshape(-1, 1)
        ])

        radiation_pattern = pyvista.PolyData(r)

        plotter = pyvista.Plotter()
        plotter.background_color = 'gray'
        plotter.add_mesh(radiation_pattern,
                         cmap='jet',
                         smooth_shading=True,
                         scalars= E_mag_dB - E_mag_dB.max(),
                         show_scalar_bar=True,
                         name='Radiation Pattern')
        # plotter.add_points(r, render_points_as_spheres=True)
        # plotter.add_axes_at_origin(xlabel='x', ylabel='y', zlabel='z', line_width=10)
        plotter.show_grid()
        plotter.show_axes()
        plotter.show()

    def create_regressors_matrix(self, int k) -> None:
        """

        Parameters
        ----------
        k

        Returns
        -------

        """

        cdef numpyc.ndarray[double, ndim=1] theta, phi

        # ---------- ---------- ---------- ---------- ---------- remove poles
        theta = self.theta[self.no_pole_cond]
        phi = self.phi[self.no_pole_cond]

        print(f'start creating regressors at {datetime.utcnow()}')
        self._regressors = Upsilonk(k, theta, self.phi, self.eta)
        print(f'finish creating regressors at {datetime.utcnow()}')

    def save(self, path: Path|str, overwrite: bool = False) -> None:

        path = Path(path)

        if not path.suffix:
            path = path.with_suffix('.rfd')

        if path.is_file() and not overwrite:
            raise FileExistsError()

        file = h5py.File(path, 'x')
        file.attrs['frequency'] = self.frequency
        file.attrs['name'] = self.name
        file.attrs['description'] = self.description

        theta_dataset = file.create_dataset('theta', data=self.theta)
        theta_dataset.attrs['unit'] = 'rad'
        theta_dataset.attrs['description'] = 'polar angle'

        phi_dataset = file.create_dataset('phi', data=self.phi)
        phi_dataset.attrs['unit'] = 'rad'
        phi_dataset.attrs['description'] = 'azimuthal angle'

        E_theta_dataset = file.create_dataset('E_theta', data=self.E_theta)
        E_theta_dataset.attrs['unit'] = 'volt'
        E_theta_dataset.attrs['description'] = 'radiation field polar component'

        E_phi_dataset = file.create_dataset('E_phi', data=self.E_phi)
        E_phi_dataset.attrs['unit'] = 'volt'
        E_phi_dataset.attrs['description'] = 'radiation field azimuthal component'

        regressor_dataset = file.create_dataset('regressors', data=self.regressors)
        regressor_dataset.attrs['unit'] = 'ohm**0.5'
        regressor_dataset.attrs['description'] = 'radiation field regressors model'

        # ---------- ---------- ---------- ---------- ---------- ----------
        regressions_dataset = file.create_group('regressions')

        file.close()

    def estimate_mode_coeff(self, 
                            int k,
                            str method = 'batch',
                            double tkp = 0.0) -> tuple:
        """

        Parameters
        ----------
        k
        method
        tkp

        Returns
        -------

        """

        cdef:
            numpyc.ndarray[numpyc.float64_t, ndim=1] theta, phi
            numpyc.ndarray[numpyc.complex128_t, ndim=1] E_theta, E_phi
            numpyc.ndarray[numpyc.complex128_t, ndim=2] Uk, Uk_H, P, E, q, res
            numpyc.ndarray[numpyc.complex128_t, ndim=2] E_est
            numpyc.ndarray[numpyc.complex128_t, ndim=1] E_est_theta, E_est_phi

            double E_max, max_rel_error

        theta = self.theta[self.no_pole_cond]
        phi = self.phi[self.no_pole_cond]
        E_theta = self.E_theta[self.no_pole_cond]
        E_phi = self.E_phi[self.no_pole_cond]

        E = numpy.hstack([
            E_theta.reshape(-1, 1),
            E_phi.reshape(-1, 1)
        ]).reshape(-1, 1)

        E_max = numpy.max(numpy.abs(E).real)

        if method == 'batch':
            Uk = Upsilonk(k, theta, phi, self.eta)
            Uk_H = Uk.conjugate().T

            if tkp == 0.0:
                P = inv_hermitian_matrix(Uk_H @ Uk)
                q = P @ Uk_H @ E
                E_est = Uk @ q
                res = E - E_est

                E_est_theta = E_est.reshape(-1, 2)[:, 0]
                E_est_phi = E_est.reshape(-1, 2)[:, 1]

                max_rel_error = numpy.max(numpy.abs(res).real) / E_max

            else:
                raise NotImplementedError()

        elif method == 'recursive':
            raise NotImplementedError()
        
        else:
            raise NotImplementedError()

        return q, P, res, max_rel_error, theta, phi, E_theta, E_phi, E_est_theta, E_est_phi

    def estimate_mode_coeff_unbiased_and_regularised(self,
                                                     int k,
                                                     double tkp = 0.0) -> tuple:

        cdef:
            numpyc.ndarray[numpyc.float64_t, ndim=1] theta, phi
            numpyc.ndarray[numpyc.complex128_t, ndim=1] E_theta, E_phi
            numpyc.ndarray[numpyc.complex128_t, ndim=2] Uk, Uk_H, P, E, q, res
            numpyc.ndarray[numpyc.complex128_t, ndim=2] Lambda, EEH_inv, K
            numpyc.ndarray[numpyc.complex128_t, ndim=2] E_est
            numpyc.ndarray[numpyc.complex128_t, ndim=1] E_est_theta, E_est_phi

            double E_max, max_rel_error

        theta = self.theta[self.no_pole_cond]
        phi = self.phi[self.no_pole_cond]
        E_theta = self.E_theta[self.no_pole_cond]
        E_phi = self.E_phi[self.no_pole_cond]

        E = numpy.hstack([
            E_theta.reshape(-1, 1),
            E_phi.reshape(-1, 1)
        ]).reshape(-1, 1)

        E_max = numpy.max(numpy.abs(E).real)

        Uk = Upsilonk(k, theta, phi, self.eta)
        Uk_H = Uk.conjugate().T

        if self.EEH_inv is None:
            self.EEH_inv = inv_hermitian_matrix(E@E.conjugate().T)

        # Lagrange multiplier:
        Lambda = tkp * inv_hermitian_matrix(Uk_H @ self.EEH_inv @ Uk)

        # Gain
        K = inv_hermitian_matrix(Uk_H @ Uk + tkp*numpy.eye(Uk.shape[1]))@\
            (Uk_H + Lambda @ Uk_H @ self.EEH_inv)

        q = K @ E
        P = K @ K.conjugate().T

        E_est = Uk @ q
        res = E - E_est

        E_est_theta = E_est.reshape(-1, 2)[:, 0]
        E_est_phi = E_est.reshape(-1, 2)[:, 1]

        max_rel_error = numpy.max(numpy.abs(res).real) / E_max

        return q, P, res, max_rel_error, theta, phi, E_theta, \
            E_phi, E_est_theta, E_est_phi

    def estimate_eigenantenna_weights(self, int lmax) -> tuple:

        cdef:
            numpyc.ndarray[numpyc.float64_t, ndim=1] theta, phi, S
            numpyc.ndarray[numpyc.complex128_t, ndim=1] E_theta, E_phi
            numpyc.ndarray[numpyc.complex128_t, ndim=2] T, Uk, Uk_H, P, E, q, res, w
            numpyc.ndarray[numpyc.complex128_t, ndim=2] Lambda, EEH_inv, K
            numpyc.ndarray[numpyc.complex128_t, ndim=2] A, AH
            numpyc.ndarray[numpyc.complex128_t, ndim=2] E_est, U, Vh, V
            numpyc.ndarray[numpyc.complex128_t, ndim=1] E_est_theta, E_est_phi

            double E_max, max_rel_error

        theta = self.theta[self.no_pole_cond]
        phi = self.phi[self.no_pole_cond]
        E_theta = self.E_theta[self.no_pole_cond]
        E_phi = self.E_phi[self.no_pole_cond]

        T = Tk(lmax, 0, 0, self.eta) * numpy.sqrt(4*numpy.pi/self.eta)

        U, S, Vh = scipy.linalg.svd(T)
        V = Vh.conjugate().T

        E = numpy.hstack([
            E_theta.reshape(-1, 1),
            E_phi.reshape(-1, 1)
        ]).reshape(-1, 1)

        E_max = numpy.max(numpy.abs(E).real)

        Uk = Upsilonk(lmax, theta, phi, self.eta)
        A = Uk @ V
        AH = A.conjugate().T

        Lambda = pdm_inverse(A)

        K = Lambda @ AH

        w = K @ E
        P = Lambda

        E_est = A @ w
        res = E - E_est

        E_est_theta = E_est.reshape(-1, 2)[:, 0]
        E_est_phi = E_est.reshape(-1, 2)[:, 1]

        max_rel_error = numpy.max(numpy.abs(res).real) / E_max

        return w, P, res, max_rel_error, theta, phi, E_theta, \
            E_phi, E_est_theta, E_est_phi


    def estimate_mode_coeff_recursive(self,
                                      int k0,
                                      int kmax,
                                      double alpha,
                                      object path,
                                      bint overwrite = False):

        # ---------- ---------- ---------- ---------- declarations
        cdef:
            numpyc.ndarray[numpyc.float64_t, ndim=1] theta, phi
            numpyc.ndarray[numpyc.float64_t, ndim=2] R, I, F, G
            numpyc.ndarray[numpyc.complex128_t, ndim=1] E_theta, E_phi

            numpyc.ndarray[numpyc.complex128_t, ndim=2] E_sam
            numpyc.ndarray[numpyc.complex128_t, ndim=2] q_plus, q_minus
            numpyc.ndarray[numpyc.complex128_t, ndim=2] P_plus, P_minus
            numpyc.ndarray[numpyc.complex128_t, ndim=2] r_plus, r_minus
            numpyc.ndarray[numpyc.complex128_t, ndim=2] E_plus, E_minus
            numpyc.ndarray[numpyc.complex128_t, ndim=2] C_plus, C_minus
            numpyc.ndarray[numpyc.complex128_t, ndim=2] U, UH
            numpyc.ndarray[numpyc.complex128_t, ndim=2] dq, dP
            numpyc.ndarray[numpyc.complex128_t, ndim=2] K
            numpyc.ndarray[numpyc.complex128_t, ndim=2] AUX, AUX2

            int N, dN

            double rre, EHE

        theta = self.theta[self.no_pole_cond]
        phi = self.phi[self.no_pole_cond]
        E_theta = self.E_theta[self.no_pole_cond]
        E_phi = self.E_phi[self.no_pole_cond]

        E_sam = numpy.hstack([
            E_theta.reshape(-1, 1),
            E_phi.reshape(-1, 1)
        ]).reshape(-1, 1)

        EHE = (E_sam.conjugate().T @ E_sam).real[0, 0]

        R = numpy.eye(E_sam.shape[0])

        # ---------- ---------- ---------- ---------- initialisation
        U = Upsilonk(k0, theta, phi, self.eta)
        UH = U.conjugate().T

        q_plus, P_plus, r_plus, *rest = self.estimate_mode_coeff(k0)
        E_plus = E_sam - r_plus
        C_plus = P_plus @ UH

        path = Path(path)

        if not path.suffix:
            path = path.with_suffix('.amc')

        if path.is_file():
            if overwrite:
                path.unlink()
            else:
                raise FileExistsError()

        file = h5py.File(path, 'w')
        file.attrs['description'] = "recursive regression fo antenna mode coeff"
        file.attrs['frequency'] = "to be implemented"

        group = file.create_group(f'{k0}')
        group.create_dataset('q', data=q_plus)
        group.create_dataset('P', data=P_plus.diagonal().real)
        group.create_dataset('r', data=r_plus)
        group.create_dataset('E_est', data=E_plus)

        rre = (r_plus.conjugate().T @ r_plus).real[0, 0] / EHE
        print(f'{k0}\t{2*P_plus.diagonal().real.mean()}\t{rre}')

        # ---------- ---------- ---------- ---------- loop
        for k in range(k0 + 1, kmax + 1):

            N = (k - 1) * (k + 1)
            dN = 2*k + 1

            F = numpy.eye(2*(N + dN), 2*N)
            G = numpy.roll(numpy.eye(2*(N + dN), 2*dN), 2*N, axis=0)

            dq = numpy.vstack([[[0j], [0j]], q_plus[2*(k-2)*k:2*N, :], [[0j], [0j]]])
            dP = alpha*numpy.max(P_plus.diagonal().real) * numpy.eye(2*dN, dtype=numpy.complex128)

            U = Upsilonk(k, theta, phi, self.eta)
            UH = U.conjugate().T

            # Mode coeff propagation
            q_minus = F @ q_plus + G @ dq

            # Covariance propagation
            P_minus = F @ P_plus @ F.T + G @ dP @ G.T

            # Auxiliary cov propagation
            C_minus = F @ C_plus

            # A priori estimate
            E_minus = U @ q_minus

            # Innovation
            r_minus = E_sam - E_minus

            # Kalman Gain
            AUX = U @ C_minus
            AUX = U @ P_minus @ UH + R - AUX - AUX.conjugate().T
            AUX2 = (P_minus @ UH - C_minus).conjugate().T
            K = scipy.linalg.solve(AUX, AUX2, assume_a='her').conjugate().T

            # Covariance Correction
            P_plus = P_minus - K @ AUX2

            # Mode coeff correction
            q_plus = q_minus + K @ r_minus

            # Auxiliary cov correction
            I = numpy.eye(K.shape[0])
            C_plus = (I - K @ U) @ C_minus + K @ R

            # A posteriori estimate
            E_plus = U @ q_plus

            # residual
            r_plus = E_sam - E_plus

            group = file.create_group(f'{k}')
            group.create_dataset('q', data=q_plus)
            group.create_dataset('P', data=P_plus.diagonal().real)
            group.create_dataset('r', data=r_plus)
            group.create_dataset('E_est', data=E_plus)

            rre = (r_plus.conjugate().T @ r_plus).real[0, 0] / EHE
            print(f'{k}\t{2 * P_plus.diagonal().real.mean()}\t{rre}')

            self.plot_mode_power_distribution(q_plus, show=True)

        file.close()

    def estimate_mode_coeff_recursive_2(self,
                                        int k0,
                                        int kmax,
                                        double sigma2,
                                        object path,
                                        bint overwrite=False):

        # ---------- ---------- ---------- ---------- ---------- declarations
        cdef:
            ndarray[float64_t, ndim=1] theta, phi
            ndarray[complex128_t, ndim=1] E_theta, E_phi

            int N, dN
            double w_cov, w_res, rre, EHE

            ndarray[complex128_t, ndim=2] E_sam, U, UH, dq, dP, K
            ndarray[complex128_t, ndim=2] q_plus, q_minus
            ndarray[complex128_t, ndim=2] P_plus, P_minus
            ndarray[complex128_t, ndim=2] r_plus, r_minus
            ndarray[complex128_t, ndim=2] E_plus, E_minus
            ndarray[complex128_t, ndim=2] C_plus, C_minus

        theta = self.theta[self.no_pole_cond]
        phi = self.phi[self.no_pole_cond]
        E_theta = self.E_theta[self.no_pole_cond]
        E_phi = self.E_phi[self.no_pole_cond]

        E_sam = numpy.hstack([
            E_theta.reshape(-1, 1),
            E_phi.reshape(-1, 1)
        ]).reshape(-1, 1)

        EHE = (E_sam.conjugate().T @ E_sam).real[0, 0]

        R = numpy.eye(E_sam.shape[0]) * sigma2

        # ---------- ---------- ---------- ---------- ---------- initialisation
        U = Upsilonk(k0, theta, phi, self.eta)
        UH = U.conjugate().T

        q_plus, P_plus, r_plus, *rest = self.estimate_mode_coeff(k0)
        E_plus = E_sam - r_plus
        C_plus = P_plus @ UH

        path = Path(path)

        if not path.suffix:
            path = path.with_suffix('.amc')

        if path.is_file():
            if overwrite:
                path.unlink()
            else:
                raise FileExistsError()

        file = h5py.File(path, 'w')
        file.attrs['description'] = "recursive regression fo antenna mode coeff"
        file.attrs['frequency'] = "to be implemented"

        group = file.create_group(f'{k0}')
        group.create_dataset('q', data=q_plus)
        group.create_dataset('P', data=P_plus.diagonal().real)
        group.create_dataset('r', data=r_plus)
        group.create_dataset('E_est', data=E_plus)

        rre = (r_plus.conjugate().T @ r_plus).real[0, 0] / EHE
        print(f'{k0}\t{2*P_plus.diagonal().real.mean()}\t{rre}')



    def plot_mode_power_distribution(self,
                numpyc.ndarray[numpyc.complex128_t, ndim=2] q, bint show=True):

        cdef:
            double total_rad_power = (q.conjugate().T @ q)[0, 0].real
            numpyc.ndarray[numpyc.float64_t, ndim=1] qTE, qTM
            int k = <int>sqrt(q.shape[0] / 2 + 1) - 1
            int l_max = k + 1
            numpyc.ndarray[numpyc.float64_t, ndim = 2] QTE, QTM

        qTE = (q.reshape(-1, 2)[:, 0].conjugate() * q.reshape(-1, 2)[:, 0]).real
        qTM = (q.reshape(-1, 2)[:, 1].conjugate() * q.reshape(-1, 2)[:, 1]).real

        QTE = numpy.zeros(shape=(l_max, 2 * l_max + 1), dtype=numpy.float64)
        QTM = numpy.zeros(shape=(l_max, 2 * l_max + 1), dtype=numpy.float64)

        l = []
        m = []

        for _l in range(1, k + 1):
            for _m in range(-_l, _l + 1):
                l.append(_l)
                m.append(_m)

        for _l, _m, _qTE, _qTM in zip(l, m, qTE, qTM):
            QTE[_l, _m] = _qTE
            QTM[_l, _m] = _qTM

        QTE = numpy.roll(QTE, l_max, axis=1).T
        QTE = QTE / total_rad_power

        QTM = numpy.roll(QTM, l_max, axis=1).T
        QTM = QTM / total_rad_power

        Q_min = min(QTE.min(), QTM.min())

        QTE[QTE < Q_min] = Q_min
        QTM[QTM < Q_min] = Q_min

        QTE = 10 * numpy.log10(QTE)
        QTM = 10 * numpy.log10(QTM)

        # ---------- ---------- ---------- ---------- ---------- ----------
        figure, axes = pyplot.subplots(1, 2, sharex=True, sharey=True)

        figure.set_size_inches(10, 20 / 3)

        # figure.suptitle(f'Mode Power Distribution at {freq} GHZ', fontsize=16)
        # figure.suptitle(f'Batch regression with $\\left|\\lambda \\right|^2 = {TKP:.1f}$ and $k = {k}$', va='bottom')

        extent = (0, k, -k, k)

        # ---------- ---------- ---------- ---------- ---------- left_axes
        cmap = 'turbo'
        cmap = 'nipy_spectral'
        cmap = 'gist_ncar'
        cmap = 'rainbow'
        cmap = 'jet'

        mappable = axes[0].imshow(QTE,
                                  cmap=cmap,
                                  origin='lower',
                                  extent=extent,
                                  interpolation='nearest')

        axes[0].set_xlabel(r'$\ell$', fontsize=24)
        axes[0].set_ylabel(r'$m$', fontsize=24)
        axes[0].set_title(rf'Normalised $\left| q_{{\ell m}}^{{\mathrm{{TE}}}}\right|^2$',
                          fontsize=26)
        axes[0].tick_params(axis='x', labelsize=22)
        axes[0].tick_params(axis='y', labelsize=22)

        axes[1].imshow(QTM,
                       cmap=cmap,
                       origin='lower',
                       extent=extent,
                       interpolation='nearest')

        axes[1].set_xlabel(r'$\ell$', fontsize=24)
        axes[1].set_title(rf'Normalised $\left| q_{{\ell m}}^{{\mathrm{{TM}}}}\right|^2$',
                          fontsize=26)
        axes[1].tick_params(axis='x', labelsize=22)

        colourbar = figure.colorbar(mappable, ax=axes)
        colourbar.set_label(f'Mode Coefficients Power - Normalised (dB)', fontsize=20)
        colourbar.ax.tick_params(labelsize=18)

        # figure.savefig(there / f'_fig_mode_batch_'
        #                        f'{freq}_k_{k}_TKP_'
        #                        f'{TKP}_NEW.png',
        #                dpi=300, transparent=True, bbox_inches='tight')
        # figure.savefig(there / f'_fig_mode_batch_'
        #                        f'{freq}_k_{k}_TKP_'
        #                        f'{TKP}_NEW_eps.eps',
        #                dpi=300, transparent=True, bbox_inches='tight')

        if show:
            pyplot.show()


    def plot_radiation_pattern_3D(self,
                numpyc.ndarray[numpyc.complex128_t, ndim=2] q,
                int n_theta=60, int n_phi=60,
                tuple theta_range=(0.0, numpy.pi),
                tuple phi_range=(0.0, 2*numpy.pi)):

        cdef:
            int k = <int> sqrt(q.shape[0] / 2 + 1) - 1
            double eta = 376.730_313_668  # ohm
            numpyc.ndarray[numpyc.float64_t, ndim=1] _theta, _phi
            numpyc.ndarray[numpyc.float64_t, ndim=2] theta, phi
            numpyc.ndarray[numpyc.float64_t, ndim=2] radiation, radiation_dB
            numpyc.ndarray[numpyc.float64_t, ndim=2] directivity, directiviy_dB
            double radiation_dB_max
            double qHq = (q.conj().T @ q).real[0, 0]
            numpyc.ndarray[numpyc.complex128_t, ndim=2] E, T
            numpyc.ndarray[numpyc.float64_t, ndim=2] x, y, z

            object plotter, radiation_pattern

            Py_ssize_t i, j

        _theta = numpy.linspace(theta_range[0], theta_range[1], n_theta)
        _phi = numpy.linspace(phi_range[0], phi_range[1], n_phi)
        theta, phi = numpy.meshgrid(_theta, _phi)
        radiation = numpy.zeros_like(theta)
        directivity = numpy.zeros_like(theta)


        for i in range(theta.shape[0]):
            for j in range(theta.shape[1]):

                T = Tk(k, theta[i, j], phi[i, j], eta)
                E = T @ q
                radiation[i, j] = ((E.conjugate().T) @ E).real[0, 0]
                directivity[i, j] = radiation[i, j]

        # radiation = radiation / 2 / eta
        directivity = directivity * 4 * numpy.pi  / eta / qHq

        # radiation[radiation == 0.0] = 0.1*radiation[radiation != 0.0].min()
        directivity[directivity == 0.0] = 0.1*directivity[directivity != 0.0].min()

        # radiation = radiation / radiation.max()

        # radiation_dB = 10.0 * numpy.log10(radiation / radiation.max())
        directivity_dB = 10.0 * numpy.log10(directivity)

        radiation_dB = directivity_dB - directivity_dB.min()

        # radiation_dB[radiation_dB == numpy.nan] = -80

        # radiation_dB -= numpy.min(radiation_dB)

        x = radiation_dB * numpy.cos(phi) * numpy.sin(theta)
        y = radiation_dB * numpy.sin(phi) * numpy.sin(theta)
        z = radiation_dB * numpy.cos(theta)

        plotter = pyvista.Plotter()
        plotter.background_color = 'gray'

        radiation_pattern = pyvista.StructuredGrid(x, y, z)

        directivity_dB[directivity_dB <= directivity_dB.max() - 30] = directivity_dB.max()-30

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

        return radiation

    # def plot_directivity




    @staticmethod
    def rotate_data(ndarray[float64_t, ndim=2] MBA,
                    ndarray[float64_t, ndim=1] theta_A,
                    ndarray[float64_t, ndim=1] phi_A,
                    ndarray[complex128_t, ndim=1] E_theta_A,
                    ndarray[complex128_t, ndim=1] E_phi_A):
        """
        
        Parameters
        ----------
        MBA : ndarray shape = (3, 3)
            Rotation matrix from system A to system B. r_B = MBA @ r_A, where
            r_A and r_B are cartesian components  of the same vector in systems
            A and B respectively.
            
        theta_A : ndarray shape = (N,)
            polar angle in system A
            
        phi_A : ndarray shape (N,)
            azimuthal angle in system A
            
        E_theta_A : ndarray shape (N,)
            radiation field theta component in system A
            
        E_phi_A : ndarray shape (N,)
            radiation field phi component in system A

        Returns
        -------
        theta_B : ndarray shape = (N,)
            polar angle in system B
            
        phi_B : ndarray shape (N,)
            azimuthal angle in system B
            
        E_theta_B : ndarray shape (N,)
            radiation field theta component in system B
            
        E_phi_B : ndarray shape (N,)
            radiation field phi component in system B

        """

        cdef:
            Py_ssize_t index

            # ---------- ---------- ---------- internal
            float64

            ndarray[float64_t, ndim=2] SA = numpy.zeros((3, 3), dtype=numpy.float64)
            ndarray[float64_t, ndim=2] SB = numpy.zeros((3, 3), dtype=numpy.float64)
            ndarray[float64_t, ndim=2] MAB_SAT = numpy.zeros((3, 3), dtype=numpy.float64)
            ndarray[float64_t, ndim=2] AUX = numpy.zeros((3, 1), dtype=numpy.float64)

            ndarray[complex128_t, ndim=2] E_sph_A = numpy.zeros((3, 1), dtype=numpy.complex128)
            ndarray[complex128_t, ndim=2] E_sph_B = numpy.zeros((3, 1), dtype=numpy.complex128)

            # ---------- ---------- ---------- output
            ndarray[float64_t, ndim = 1] theta_B = numpy.zeros_like(theta_A)
            ndarray[float64_t, ndim = 1] phi_B = numpy.zeros_like(phi_A)
            ndarray[complex128_t, ndim = 1] E_theta_B = numpy.zeros_like(E_theta_A)
            ndarray[complex128_t, ndim = 1] E_phi_B = numpy.zeros_like(E_theta_B)

        for index in range(theta_B.shape[0]):

            SA = rot_matrix_sph_from_cart(theta_A[index], phi_A[index])

            E_sph_A[0, 0] = 0.0
            E_sph_A[1, 0] = E_theta_A[index]
            E_sph_A[2, 0] = E_phi_A[index]

            MAB_SAT = MBA @ SA.T

            AUX = MAB_SAT @ numpy.array([[1, 0, 0]], dtype=numpy.float64).T

            if AUX[2] >= 1:
                theta_B[index] = 0.0

            elif AUX[2] <= -1:
                theta_B[index] = pi

            else:
                theta_B[index] = acos(AUX[2])

            phi_B[index] = atan2(AUX[1], AUX[0])

            SB = rot_matrix_sph_from_cart(theta_B[index], phi_B[index])

            E_sph_B = SB @ MBA @ SA.T @ E_sph_A

            # ---------- ---------- ----------
            E_theta_B[index] = E_sph_B[1]
            E_phi_B[index] = E_sph_B[2]

        return numpy.asarray(theta_B), \
            numpy.asarray(phi_B), \
            numpy.asarray(E_theta_B), \
            numpy.asarray(E_phi_B)

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

    @property
    def regressors(self) -> ArrayLike:
        if self._regressors is None:
            return

        return numpy.asarray(self._regressors)