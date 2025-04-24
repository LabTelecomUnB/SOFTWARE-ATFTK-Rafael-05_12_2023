#  -*- coding: utf-8 -*-
"""
Author: Rafael R. L. Benevides
"""

# ---------- ---------- ---------- ---------- ---------- ---------- py
import numpy

from scipy.special.cython_special import sph_harm

# ---------- ---------- ---------- ---------- ---------- ---------- cy
cimport numpy as numpyc

from numpy cimport ndarray

from libc.math cimport sin, cos, sqrt, pi
from libc.stdlib cimport abs


# ---------- ---------- ---------- ---------- ---------- ---------- ty


cdef double eta = 376.730_313_668  # ohm


cpdef double complex Ylm(int l, int m, double theta, double phi):
    r"""
    Calculates the spherical harmonic :math:`Y_{\ell}^m\left(\theta,
    \phi\right)`.
    
    .. math::
        Y_{\ell}^{m}\left( \theta, \phi \right)
        =
        \sqrt{
            \frac{2\ell + 1}{4\pi}
            \frac{
                \left( \ell - m \right)!
            }{
                \left( \ell + m \right)!
            }
        }
        P_{\ell}^{m}\left( \cos\theta \right)
        e^{jm\phi}
    
    Parameters
    ----------
    l : int
        degree. Must be greater than zero.

    m: int
        order.

    theta: double
        polar angle in radians

    phi: double
        azimuthal angle in radians

    Returns
    -------
    out : complex
        The :math:`\ell`-degree :math:`m`-order spherical harmonic 
        :math:`Y_{\ell}^m\left(\theta, \phi\right)`
    
    See also
    --------
    :func:`m_times_Ylm_over_sin_theta`
    :func:`dYlm_dtheta`
    
    """

    if abs(m) > l:
        return 0.0 + 0.0j

    return sph_harm(m, l, phi, theta)


cpdef double complex m_times_Ylm_over_sin_theta(int l, int m, double theta,  double phi):
    r"""
    Calculates :math:`\dfrac{mY_{\ell}^m\left(\theta,
    \phi\right)}{\sin\theta}`.
    
    Parameters
    ----------
    l : int
        degree. Must be greater than zero.

    m: int
        order.

    theta: double
        polar angle in radians

    phi: double
        azimuthal angle in radians

    Returns
    -------
    out : complex
        The :math:`\ell`-degree :math:`m`-order first spherical eigenfunction
        :math:`\dfrac{mY_{\ell}^m\left(\theta, \phi\right)}{\sin\theta}`

    See also
    --------
    :func:`dYlm_dtheta`
    :func:`Ylm`
    :func:`Tlm`

    """
    cdef:
        double complex coeff


    if theta == 0.0:

        if abs(m) != 1:
            return 0j

        coeff = -1.0 / 4.0 * sqrt((2.0*l + 1.0) / pi * l * (l + 1.0))

        return coeff * numpy.exp(1j*m*phi)

    if theta == pi:

        if abs(m) != 1:
            return 0j

        coeff = 1.0 / 4.0 * sqrt((2.0*l + 1.0) / pi * l * (l + 1.0)) * (-1)**l

        return coeff * numpy.exp(1j*m*phi)

    return Ylm(l, m, theta, phi) / sin(theta) * m


cpdef double complex dYlm_dtheta(int l, int m, double theta,  double phi):
    r"""
    Calculates
    :math:`\dfrac{\partial Y_{\ell}^m}{\partial\theta}\left(\theta,\phi\right)`.
    
    Parameters
    ----------
    l : int
        degree. Must be greater than zero.

    m: int
        order.

    theta: double
        polar angle in radians

    phi: double
        azimuthal angle in radians

    Returns
    -------
    out : complex
        The :math:`\ell`-degree :math:`m`-order second spherical eigenfunction
        :math:`\dfrac{\partial Y_{\ell}^m}{\partial\theta}\left(\theta,\phi\right)`

    See also
    --------
    :func:`m_times_Ylm_over_sin_theta`
    :func:`Ylm`
    :func:`Tlm`

    """
    cdef:
        double sin_theta, cos_theta

        double complex coeff, numerator


    if theta == 0.0:

        if abs(m) != 1:
            return 0j

        coeff = -m / 4.0 * sqrt((2.0 * l + 1.0) / pi * l * (l + 1.0))

        return coeff * numpy.exp(1j * m * phi)

    if theta == pi:

        if abs(m) != 1:
            return 0j

        coeff = -m / 4.0 * sqrt((2.0 * l + 1.0) / pi * l * (l + 1.0)) * (
            -1) ** l

        return coeff * numpy.exp(1j * m * phi)

    sin_theta = sin(theta)
    cos_theta = cos(theta)

    coeff = sqrt((2.0 * l + 1.0) / (2.0 * l - 1.0) * (l - m) * (l + m))

    numerator = l * cos_theta * Ylm(l, m, theta, phi) - coeff * Ylm(l - 1, m, theta, phi)

    return numerator / sin_theta


cpdef ndarray[double complex, ndim=2] Tlm(int l, int m, double theta,  double phi):
    r"""
    Computes the matrix-eigenfunction 
    :math:`\mathbf{T}_{\ell m}\left(\theta, \phi\right)`.
    
    The matrix-eigenfunction 
    :math:`\mathbf{T}_{\ell m}\left(\theta, \phi\right)\in\mathbb{C}^{2\times2}`
    is define as
     
    .. math::
        \mathbf{T}_{\ell m}\left(\theta, \phi\right)
        =
        \frac{
            j^{\ell}
            \sqrt{\eta}
        }{
            \sqrt{\ell \left(\ell + 1\right)}
        }
        \begin{bmatrix}
            \dfrac{
                -m Y_{l}^{m}\left(\theta, \phi\right)
            }{
                \sin\theta
            }
            &
            \dfrac{
                \partial Y_{l}^{m}
            }{
                \partial \theta
            }
            \left(\theta, \phi\right)
            \\[2mm]
            -j
            \dfrac{
                \partial Y_{l}^{m}
            }{
                \partial \theta
            }
            \left(\theta, \phi\right)
            &
            \dfrac{
                j m Y_{l}^{m}\left(\theta, \phi\right)
            }{
                \sin\theta
            }
        \end{bmatrix}
    
    and it is used to reconstruct the antenna radiation field as
    
    .. math::
        \mathbf{\mathtt{E}}\left(\theta, \phi\right)
        =
        \sum_{\ell=1}^{\infty}
        \sum_{m=-\ell}^{\ell}
            \mathbf{T}_{\ell m}\left(\theta, \phi\right)
            \mathbf{q}_{\ell m}
    
    where :math:`\mathbf{q}_{\ell m}\in\mathbb{C}^{2\times 1}` are the 
    *spherical mode coefficients* (SMC).
            
    Parameters
    ----------
    l : int
        degree. Must be greater than zero.

    m: int
        order.

    theta: double
        polar angle in radians

    phi: double
        azimuthal angle in radians

    Returns
    -------
    out : complex
        The :math:`\ell`-degree :math:`m`-order matrix-eigenfunction
        :math:`\mathbf{T}_{\ell m}\left(\theta, \phi\right)` in units of 
        :math:`\Omega^{\frac{1}{2}}` (square root of ohm)
    
    See also
    --------
    :func:`m_times_Ylm_over_sin_theta`
    :func:`dYlm_dtheta`
    :func:`T`

    """
    cdef:
        double complex AUX1 = m_times_Ylm_over_sin_theta(l, m, theta, phi)
        double complex AUX2 = dYlm_dtheta(l, m, theta, phi)
        double complex coeff

        double complex[:, :] _Tlm = numpy.zeros((2, 2), dtype=numpy.complex128)

    if l % 4 == 0:
        coeff = 1.0 + 0.0j

    elif l % 4 == 1:
        coeff = 0.0 + 1.0j

    elif l % 4 == 2:
        coeff = -1.0 + 0.0j

    else:
        coeff = 0.0 - 1.0j

    coeff *= sqrt(eta / l / (l + 1.0))

    _Tlm[0, 0] = -AUX1
    _Tlm[0, 1] = +AUX2
    _Tlm[1, 0] = -1j * AUX2
    _Tlm[1, 1] = +1j * AUX1

    return coeff * numpy.asarray(_Tlm)


cpdef ndarray[double complex, ndim=2] T(int lmax, double theta, double phi):
    r"""Computes the model matrix :math:`\mathbf{\mathsf{T}}\left(\theta, \phi\right)`.

    :math:`\mathbf{\mathsf{T}}\left(\theta, \phi\right) \in \mathbb{C}^{2\times p}` is a 
    matrix of spherical eigenfunctions defined by the blocks

    .. math::
        \mathbf{\mathsf{T}}\left(\theta, \phi\right) = \begin{bmatrix}
            \mathbf{T}_{(1, -1)} & \mathbf{T}_{(1, 0)} & \mathbf{T}_{(1, 
            1)} & \cdots & \mathbf{T}_{(\ell_{\max}, \ell_{\max})}
        \end{bmatrix}

    where

    .. math::
        \mathbf{T}_{\ell m}\left(\theta, \phi\right)
        =
        \frac{
            j^{\ell}
            \sqrt{\eta}
        }{
            \sqrt{\ell \left(\ell + 1\right)}
        }
        \begin{bmatrix}
            \dfrac{
                -m Y_{l}^{m}\left(\theta, \phi\right)
            }{
                \sin\theta
            }
            &
            \dfrac{
                \partial Y_{l}^{m}
            }{
                \partial \theta
            }
            \left(\theta, \phi\right)
            \\[2mm]
            -j
            \dfrac{
                \partial Y_{l}^{m}
            }{
                \partial \theta
            }
            \left(\theta, \phi\right)
            &
            \dfrac{
                j m Y_{l}^{m}\left(\theta, \phi\right)
            }{
                \sin\theta
            }
        \end{bmatrix},

    while :math:`p = 2\ell_{\max}\cdot\left(\ell_{\max} + 2\right)` is the number of
    spherical modes, where :math:`\ell_{\max}` is the degree of the highest 
    mode, used to project the radiation field as the linear form
    :math:`\mathbf{\mathtt{
    E}}\left(\theta, \phi\right) = \mathbf{
    \mathsf{T}}\left(\theta, \phi\right) \mathbf{\mathsf{q}}`, 
    :math:`\mathbf{\mathsf{q}} \in \mathbb{C}^{p \times 1}`.

    Parameters
    ----------
    lmax : int
        Highest degree (or truncation degree). Must be greater than zero.

    theta : double
        Polar angle in radians

    phi : double
        Azimuthal angle in radians

    Returns
    -------
    out : complex-dtyped (2, p)-shaped :py:class:`ndarray <numpy.ndarray>`
        The matrix :math:`\mathbf{\mathsf{T}}\left(\theta, \phi\right)` in units of 
        :math:`\Omega^{\frac{1}{2}}` (square root of ohm).
        
    See also
    --------
    :func:`Tlm`
    :func:`Upsilon`

    """
    cdef:
        Py_ssize_t p = 2 * lmax * (lmax + 2)
        Py_ssize_t l, m, col = 0
        ndarray[double complex, ndim=2] _T = numpy.zeros(shape=(2, p), dtype=numpy.complex128)

    for l in range(1, lmax+1):
        for m in range(-l, l+1):

            _T[:, col:col+2] = Tlm(l, m, theta, phi)

            col += 2

    return _T


cpdef ndarray[double complex, ndim=2] Upsilon(int lmax, double[:] theta, double[:] phi):
    r"""
    Computes the regressor matrix :math:`\mathbf{\Upsilon}\left(\theta, \phi\right)`.
    
    Given a sequence of :math:`\frac{s}{2}` directions
    :math:`\left(\theta_1, \phi_1\right)`,
    :math:`\left(\theta_2, \phi_2\right)`,
    :math:`\cdots`,
    :math:`\left(\theta_{\frac{s}{2}}, \phi_{\frac{s}{2}}\right)`
    (possibly at which samples of the radiation field have been collected) and
    :math:`p = 2\ell_{\max}\cdot\left(\ell_{\max} + 2\right)`, the regressor 
    matrix :math:`\mathbf{\Upsilon} \in \mathbb{C}^{s\times p}` is defined by
    the blocks
    
    .. math::
        \mathbf{\Upsilon}
        =
        \begin{bmatrix}
            \mathbf{\mathsf{T}}\left(\theta_1, \phi_1\right) \\ 
            \mathbf{\mathsf{T}}\left(\theta_2, \phi_2\right) \\
            \vdots \\
            \mathbf{\mathsf{T}}\left(\theta_{\frac{s}{2}}, \phi_M\right) \\
        \end{bmatrix}.
    
    Parameters
    ----------
    lmax : :py:class:`int`
        Highest degree (or truncation degree). Must be greater than zero.

    theta : double-dtyped (s/2,)-shaped :py:class:`ndarray <numpy.ndarray>`
        Polar angle in radians

    phi : double-dtyped (s/2,)-shaped :py:class:`ndarray <numpy.ndarray>`
        Azimuthal angle in radians

    Returns
    -------
    out : complex-dtyped (s, p)-shaped :py:class:`ndarray <numpy.ndarray>`
        The matrix :math:`\mathbf{\Upsilon}\left(\theta, \phi\right)` in units of 
        :math:`\Omega^{\frac{1}{2}}` (square root of ohm).

    See also
    --------
    :func:`T`

    """
    cdef:
        Py_ssize_t p = 2 * lmax * (lmax + 2)
        Py_ssize_t s = 2 * theta.shape[0]

        Py_ssize_t index

        ndarray[double complex, ndim = 2] _Treg = numpy.zeros(shape=(s, p), dtype=numpy.complex128)

    for index in range(theta.shape[0]):
        _Treg[2*index:2*index+2, :] = T(lmax, theta[index], phi[index])

    return _Treg


