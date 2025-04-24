#  -*- coding: utf-8 -*-
"""
Author: Rafael R. L. Benevides
"""

# ---------- ---------- ---------- ---------- ---------- py
import numpy
from scipy.linalg import solve_triangular

# ---------- ---------- ---------- ---------- ---------- cy
cimport numpy as numpyc

from libc.math cimport sin, cos, sqrt, pi
from libc.stdlib cimport abs

from scipy.special.cython_special import sph_harm

# ---------- ---------- ---------- ---------- ---------- typing


# ========== ========== ========== ========== ========== ==========
# Spherical Harmonic Y_l^m
# ========== ========== ========== ========== ========== ==========
cdef double complex Ylm(int l, int m, double theta, double phi):
    r"""Calculates the spherical harmonic :math:`Y_{\ell}^m\left(\theta,
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
    e^{im\phi}
        

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

    """

    if abs(m) > l:
        return 0.0 + 1j * 0.0

    return sph_harm(m, l, phi, theta)


# ========== ========== ========== ========== ========== ==========
# T_k matrix
# ========== ========== ========== ========== ========== ==========
cpdef numpyc.ndarray[numpyc.complex128_t, ndim=2] Tk(int k,
                                                     double theta,
                                                     double phi,
                                                     double eta):
    r"""
    Computes the matrix :math:`\mathbf{\mathsf{T}}_k\left(\theta, \phi\right) 
    \in \mathbb{C}^{2\times2N_k}`.

    :math:`\mathbf{\mathsf{T}}_k\left(\theta, \phi\right)` is a 
    matrix of spherical eigenfunctions defined by
    the blocks

    .. math::
        \mathbf{\mathsf{T}}_k\left(\theta, \phi\right) = \begin{bmatrix}
            \mathbf{T}_{(1, -1)} & \mathbf{T}_{(1, 0)} & \mathbf{T}_{(1, 
            1)} & \cdots & \mathbf{T}_{(k, k)}
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

    while :math:`N_k = k^2 + 2k` is the number of spherical modes, where 
    :math:`k` is the degree of the highest mode, used to project the
    radiation field as the linear form  :math:`\mathbf{\mathtt{
    E}}_k\left(\theta, \phi\right) = \mathbf{
    \mathsf{T}}_k\left(\theta, \phi\right) \mathbf{\mathsf{q}}_k`, 
    :math:`\mathbf{\mathsf{q}}_k \in \mathbb{C}^{2N_k \times 1}`.


    Parameters
    ----------
    k: int
        highest degree. Must be greater than zero.

    theta: double
        polar angle in radians

    phi: double
        azimuthal angle in radians

    eta: double
        medium impedance :math:`\sqrt{\dfrac{\mu}{\varepsilon}}`


    Returns
    -------
    out : numpy.ndarray[numpy.complex128_t, ndim=2] with shape (2, 2Nk)
        The matrix :math:`\mathbf{\mathsf{T}}_k\left(\theta, \phi\right)` as 
        explained above.
    """
    cdef:
        unsigned int Nk = k * (k + 2)
        int l, m, col = 0
        double sin_theta, cos_theta, Clm
        double delta = 0.0
        double complex coeff, Ylm_curr, A, B, aux, coeff_exp_jmphi
        double complex[:, :] T = numpy.zeros(shape=(2, 2*Nk), dtype=numpy.complex128)

    # if theta == 0.0 or theta == pi:

    #     coeff = -1j / 4 * sqrt(3 * eta / pi) * (cos(phi) + 1j * sin(phi))

    #     T[0, 0] = coeff
    #     T[1, 1] = -1j * coeff

    #     T[0, 4] = -1 * coeff
    #     T[1, 5] = 1j * coeff

    #     return numpy.asarray(T)

    if theta <= delta:
        # theta = 0.001

        for l in range(1, k+1):

            if l % 4 == 0:
                aux = 1.0 + 0.0j
            elif l % 4 == 1:
                aux = 1.0j
            elif l % 4 == 2:
                aux = -1.0 + 0.0j
            else:
                aux = -1.0j

            coeff = aux * sqrt(2.0*l + 1.0)

            for m in range(-l, l+1):

                if abs(m) == 1:

                    # coeff_exp_jmphi = coeff * (cos(m*phi) + 1j*sin(m*phi))
                    if phi == pi:
                        coeff_exp_jmphi = coeff * -1
                    else:
                        coeff_exp_jmphi = coeff * numpy.exp(1j*m*phi)

                    T[0, col] = coeff_exp_jmphi
                    T[0, col + 1] = -m*coeff_exp_jmphi
                    T[1, col] = 1j*m*coeff_exp_jmphi
                    T[1, col + 1] = -1j*coeff_exp_jmphi

                col += 2

        return numpy.asarray(T) * sqrt(eta/4/pi) / 2

    if theta >= pi - delta:
        #theta = pi - 0.001

        for l in range(1, k + 1):

            if l % 4 == 0:
                aux = 1.0 + 0.0j
            elif l % 4 == 1:
                aux = 1.0j
            elif l % 4 == 2:
                aux = -1.0 + 0.0j
            else:
                aux = -1.0j

            coeff = sqrt(2.0 * l + 1.0) * aux * (-1)**l

            for m in range(-l, l + 1):

                if abs(m) == 1:

                    # coeff_exp_jmphi = coeff * (cos(m * phi) + 1j * sin(m * phi))
                    # coeff_exp_jmphi = coeff * numpy.exp(1j*m*phi)

                    if phi == pi:
                        coeff_exp_jmphi = coeff * -1
                    else:
                        coeff_exp_jmphi = coeff * numpy.exp(1j*m*phi)

                    T[0, col] = -coeff_exp_jmphi
                    T[0, col + 1] = -m * coeff_exp_jmphi
                    T[1, col] = 1j * m * coeff_exp_jmphi
                    T[1, col + 1] = 1j * coeff_exp_jmphi

                col += 2

        return numpy.asarray(T) * sqrt(eta / 4.0 / pi) / 2.0

    sin_theta = sin(theta)
    cos_theta = cos(theta)

    for l in range(1, k+1):

        if l % 4 == 0:
            aux = 1.0 + 0.0j
        elif l % 4 == 1:
            aux = 1.0j
        elif l % 4 == 2:
            aux = -1.0 + 0.0j
        else:
            aux = -1.0j

        # coeff = 1j ** (l % 4) / sqrt(l * (l + 1)) / sin_theta
        coeff = aux / sqrt(l * (l + 1.0)) / sin_theta

        for m in range(-l, l+1):
            Clm = sqrt((2.0*l + 1.0) / (2.0*l - 1.0) * (l - m) * (l + m))
            Ylm_curr = Ylm(l, m, theta, phi)

            # ---------- ---------- ---------- ---------- ---------- ----------
            A = Clm * Ylm(l-1, m, theta, phi) - l * cos_theta * Ylm_curr
            B = m * Ylm_curr

            # ---------- ---------- ---------- ---------- ---------- ----------
            T[0, col] = -coeff * B
            T[0, col + 1] = -coeff * A
            T[1, col] = 1j * coeff * A
            T[1, col + 1] = 1j * coeff * B

            # ---------- ---------- ---------- ---------- ---------- ----------
            col += 2

    return numpy.asarray(T) * sqrt(eta)


# ========== ========== ========== ========== ========== ==========
# Upsilon_k matrix
# ========== ========== ========== ========== ========== ==========
cpdef numpyc.ndarray[numpyc.complex128_t, ndim=2] Upsilonk(int k,
                                                           double[:] theta,
                                                           double[:] phi,
                                                           double eta):
    r"""
    Computes the regressor matrix :math:`\mathbf{\Upsilon}_k \in \mathbb{C}^{
    2M\times2N_k}`.

    Given a sequence of :math:`M` directions
    :math:`\left(\theta_1, \phi_1\right)`,
    :math:`\left(\theta_2, \phi_2\right)`,
    :math:`\cdots`,
    :math:`\left(\theta_M, \phi_M\right)`
    (possibly at which samples of the radiation field have been collected), 
    the regressor matrix :math:`\mathbf{\Upsilon}_k` is defined by the blocks

    .. math::
        \mathbf{\Upsilon}_k
        =
        \begin{bmatrix}
            \mathbf{\mathsf{T}}_k\left(\theta_1, \phi_1\right) \\ 
            \mathbf{\mathsf{T}}_k\left(\theta_2, \phi_2\right) \\
            \vdots \\
            \mathbf{\mathsf{T}}_k\left(\theta_M, \phi_M\right) \\
        \end{bmatrix}.

    Parameters
    ----------
    k: int
        highest degree. Must be greater than zero.

    theta: numpy.ndarray[numpy.float, ndim=1] with shape (M,)
        polar angle in radians

    phi: numpy.ndarray[numpy.float, ndim=1]  with shape (M,)
        azimuthal angle in radians

    eta: double
        medium impedance :math:`\sqrt{\dfrac{\mu}{\varepsilon}}`

    Returns
    -------
    out : numpy.ndarray[numpy.complex128_t, ndim=2] with shape (2M, 2Nk)
        The matrix :math:`\mathbf{\Upsilon}_k` as 
        explained above.

    See also
    --------
    :func:`Tk`

    """
    cdef:
        unsigned int Nk = k * (k + 2), M = theta.shape[0]
        numpyc.ndarray[numpyc.complex128_t, ndim=2] gamma = numpy.zeros(shape=(2*M, 2*Nk), dtype=numpy.complex128)
        numpyc.ndarray[numpyc.complex128_t, ndim=2] T = numpy.zeros(shape=(2, 2*Nk), dtype=numpy.complex128)
        Py_ssize_t index

    for index in range(M):
        gamma[2*index:2*index+2, :] = Tk(k, theta[index], phi[index], eta)

    return gamma

