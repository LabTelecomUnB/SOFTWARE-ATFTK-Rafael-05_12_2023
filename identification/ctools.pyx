
#  -*- coding: utf-8 -*-
"""

Author: Rafael R. L. Benevides
Date: 12/1/22

"""


import numpy
cimport numpy
import scipy.linalg
from libc.math cimport sin, cos, sqrt, pi
from libc.stdlib cimport abs

from scipy.linalg import solve, solve_triangular
from scipy.special.cython_special cimport sph_harm


numpy.import_array()


# ========== ========== ========== ========== ========== ==========
# Spherical Harmonic
# ========== ========== ========== ========== ========== ==========
cpdef inline double complex _ylm(int l, int m, double theta, double phi):
    r"""
    Calculates the spherical harmonic :math:`Y_{\ell}^m\left(\theta, 
    \phi\right)`.
    
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
        return 0.0 + 1j*0.0

    return sph_harm(m, l, phi, theta)


# ========== ========== ========== ========== ========== ==========
# Tk matrix
# ========== ========== ========== ========== ========== ==========
cpdef numpy.ndarray[numpy.complex128_t, ndim=2] Tk(int k, double theta,
                                                   double phi, double eta):
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
        double complex coeff, Ylm_curr, A, B, aux
        double complex[:, :] T = numpy.zeros(shape=(2, 2*Nk), dtype=numpy.complex128)

    if theta == 0.0 or theta == pi:

        coeff = -1j / 4 * sqrt(3 * eta / pi) * (cos(phi) + 1j * sin(phi))

        T[0, 0] = coeff
        T[1, 1] = -1j * coeff

        T[0, 4] = -1 * coeff
        T[1, 5] = 1j * coeff

        return numpy.asarray(T)

    sin_theta = sin(theta)
    cos_theta = cos(theta)

    for l in range(1, k+1):

        if l % 4 == 0:
            aux = 1 + 0j
        elif l % 4 == 1:
            aux = 1j
        elif l % 4 == 2:
            aux = -1 + 0j
        else:
            aux = -1j

        # coeff = 1j ** (l % 4) / sqrt(l * (l + 1)) / sin_theta
        coeff = aux / sqrt(l * (l + 1)) / sin_theta

        for m in range(-l, l+1):
            Clm = sqrt((2*l + 1) / (2*l - 1) * (l - m) * (l + m))
            Ylm_curr = _ylm(l, m, theta, phi)

            # ---------- ---------- ---------- ---------- ---------- ----------
            A = Clm * _ylm(l-1, m, theta, phi) - l * cos_theta * Ylm_curr
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
# Gamma_k matrix
# ========== ========== ========== ========== ========== ==========
cpdef numpy.ndarray[numpy.complex128_t, ndim=2] Upsilonk(int k, double[:] theta, double[:] phi, double eta):
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
        numpy.ndarray[numpy.complex128_t, ndim=2] gamma = numpy.zeros(shape=(2*M, 2*Nk), dtype=numpy.complex128)
        numpy.ndarray[numpy.complex128_t, ndim=2] T = numpy.zeros(shape=(2, 2*Nk), dtype=numpy.complex128)
        Py_ssize_t index

    for index in range(M):
        gamma[2*index:2*index+2, :] = Tk(k, theta[index], phi[index], eta)

    return gamma


# ========== ========== ========== ========== ========== ==========
# Batch Estimator
# ========== ========== ========== ========== ========== ==========
def q_estimate_batch(
        numpy.ndarray[numpy.complex128_t, ndim=2] Upsilon,
        numpy.ndarray[numpy.complex128_t, ndim=1] E,
        double lambda2):
    r"""
    Computes :math:`\hat{\mathbf{\mathsf{q}}}_k\in\mathbb{C}^{2N_k\times1}`
    using the Tikhonov-Philips Regularisation method.

    Given a regressor matrix :math:`\mathbf{\Upsilon}_k \in \mathbb{C}^{
    2M\times2N_k}` (``Upsilon``) and radiation field samples :math:`\widetilde{\mathbf{
    \mathsf{E}}} \in \mathbb{C}^{2M\times1}` (``E``) related by the linear
    model :math:`\widetilde{\mathbf{\mathsf{E}}} = \mathbf{\Upsilon}_k \mathbf{\mathsf{
    q}}_k + \mathbf{\epsilon}_k` (where the noise :math:`\mathbf{
    \epsilon}_k` is the sum of the measurement noise with the truncation
    error), the estimator :math:`\hat{\mathbf{\mathsf{q}}}_k` is defined as
    the solution of the minimisation of

    .. math::
        \left\|
            \mathbf{\mathsf{E}}
            -
            \mathbf{\Upsilon}_k
            \mathbf{\mathsf{q}}_k
        \right\|^2
        +
        \left|\lambda\right|^2
        \left\|
            \mathbf{\mathsf{q}}_k
        \right\|^2,

    which happens to be given by
    :math:`\hat{\mathbf{\mathsf{q}}}_k = \left[\mathbf{\Upsilon}_k^{
    \mathrm{H}}\mathbf{\Upsilon}_k + \left|\lambda\right|^2 \mathbf{
    I}_{2N_k}\right]^{-1} \mathbf{\Upsilon}_k^{
    \mathrm{H}} \widetilde{\mathbf{\mathsf{E}}}`.

    Parameters
    ----------
    Upsilon: numpy.ndarray[numpy.complex128_t, ndim=2] with shape (2M, 2Nk)
        Regressor matrix.

    E: numpy.ndarray[numpy.complex128_t, ndim=1] with shape (2M,)
        Radiation field samples matrix.

    Returns
    -------
    q: numpy.ndarray[numpy.complex128_t, ndim=1] with shape (2Nk, 1)
        Estimated mode coefficients :math:`\hat{\mathbf{\mathsf{q}}}_k`.

    P: numpy.ndarray[numpy.complex128_t, ndim=2] with shape (2Nk, 2Nk)
        Normalised Covariance Matrix :math:`\mathbf{\mathsf{P}}_k = \mathrm{
        cov}\left(\hat{\mathbf{\mathsf{q}}}_k\right) \in \mathbb{C}^{2N_k\times 2N_k}`.

    """
    cdef:
        numpy.ndarray[numpy.complex128_t, ndim=2] Upsilon_H = Upsilon.conj().T
        numpy.ndarray[numpy.complex128_t, ndim=2] I = numpy.eye(
            Upsilon_H.shape[0], dtype=numpy.complex128)
        numpy.ndarray[numpy.complex128_t, ndim=2] P
        numpy.ndarray[numpy.complex128_t, ndim=1] q

    if lambda2 == 0.0:
        P = solve(Upsilon_H @ Upsilon, I, assume_a='her')
        q = P @ Upsilon_H @ E

    else:
        P = solve(Upsilon_H @ Upsilon + lambda2 * I, I, assume_a='her')
        q = P @ Upsilon_H @ E
        P = P - lambda2 * P @ P

    return q, P


cpdef numpy.ndarray[numpy.complex128_t, ndim=1] solve_jacobi(
        numpy.ndarray[numpy.complex128_t, ndim=2] A,
        numpy.ndarray[numpy.complex128_t, ndim=2] b,
        double rtol):

    cdef:
        Py_ssize_t index
        numpy.ndarray[numpy.complex128_t, ndim = 2] L = -A

        numpy.ndarray[numpy.complex128_t, ndim = 2] x_prev = numpy.zeros_like(b)
        numpy.ndarray[numpy.complex128_t, ndim = 2] x_curr = b.copy()
        double eps = 1.0


    for index in range(A.shape[0]):
        L[index, index] += 1

    while eps >= rtol:
        x_curr = L @ x_prev + b

        eps = numpy.max(numpy.abs(x_curr - x_prev)) / numpy.max(numpy.abs(x_curr))

        print(eps)

        x_prev = x_curr

    return x_curr


cpdef numpy.ndarray[numpy.complex128_t, ndim = 2] inv_hermitian_matrix(
        numpy.ndarray[numpy.complex128_t, ndim = 2] A):

    cdef:
        numpy.ndarray[numpy.complex128_t, ndim = 2] L = numpy.linalg.cholesky(A)
        numpy.ndarray[numpy.complex128_t, ndim = 2] eye = numpy.eye(A.shape[0], dtype=numpy.complex128)
        numpy.ndarray[numpy.complex128_t, ndim = 2] y = numpy.zeros_like(A, dtype=numpy.complex128)
        numpy.ndarray[numpy.complex128_t, ndim = 2] x = numpy.zeros_like(A, dtype=numpy.complex128)

    y = scipy.linalg.solve_triangular(L, eye, lower=True)
    x = scipy.linalg.solve_triangular(L.conj().T, y, lower=False)

    return x

def q_estimate_batch_v2(
        numpy.ndarray[numpy.complex128_t, ndim=2] Upsilon,
        numpy.ndarray[numpy.complex128_t, ndim=1] E,
        double lambda2):
    r"""
    Computes :math:`\hat{\mathbf{\mathsf{q}}}_k\in\mathbb{C}^{2N_k\times1}`
    using the Tikhonov-Philips Regularisation method.

    Given a regressor matrix :math:`\mathbf{\Upsilon}_k \in \mathbb{C}^{
    2M\times2N_k}` (``Upsilon``) and radiation field samples :math:`\widetilde{\mathbf{
    \mathsf{E}}} \in \mathbb{C}^{2M\times1}` (``E``) related by the linear
    model :math:`\widetilde{\mathbf{\mathsf{E}}} = \mathbf{\Upsilon}_k \mathbf{\mathsf{
    q}}_k + \mathbf{\epsilon}_k` (where the noise :math:`\mathbf{
    \epsilon}_k` is the sum of the measurement noise with the truncation
    error), the estimator :math:`\hat{\mathbf{\mathsf{q}}}_k` is defined as
    the solution of the minimisation of

    .. math::
        \left\|
            \mathbf{\mathsf{E}}
            -
            \mathbf{\Upsilon}_k
            \mathbf{\mathsf{q}}_k
        \right\|^2
        +
        \left|\lambda\right|^2
        \left\|
            \mathbf{\mathsf{q}}_k
        \right\|^2,

    which happens to be given by
    :math:`\hat{\mathbf{\mathsf{q}}}_k = \left[\mathbf{\Upsilon}_k^{
    \mathrm{H}}\mathbf{\Upsilon}_k + \left|\lambda\right|^2 \mathbf{
    I}_{2N_k}\right]^{-1} \mathbf{\Upsilon}_k^{
    \mathrm{H}} \widetilde{\mathbf{\mathsf{E}}}`.

    Parameters
    ----------
    Upsilon: numpy.ndarray[numpy.complex128_t, ndim=2] with shape (2M, 2Nk)
        Regressor matrix.

    E: numpy.ndarray[numpy.complex128_t, ndim=1] with shape (2M,)
        Radiation field samples matrix.

    Returns
    -------
    q: numpy.ndarray[numpy.complex128_t, ndim=1] with shape (2Nk, 1)
        Estimated mode coefficients :math:`\hat{\mathbf{\mathsf{q}}}_k`.

    P: numpy.ndarray[numpy.complex128_t, ndim=2] with shape (2Nk, 2Nk)
        Normalised Covariance Matrix :math:`\mathbf{\mathsf{P}}_k = \mathrm{
        cov}\left(\hat{\mathbf{\mathsf{q}}}_k\right) \in \mathbb{C}^{2N_k\times 2N_k}`.

    """
    cdef:
        numpy.ndarray[numpy.complex128_t, ndim=2] Upsilon_H = Upsilon.conj().T
        numpy.ndarray[numpy.complex128_t, ndim=2] I = numpy.eye(
            Upsilon_H.shape[0], dtype=numpy.complex128)
        numpy.ndarray[numpy.complex128_t, ndim=2] P
        numpy.ndarray[numpy.complex128_t, ndim=1] q

    if lambda2 == 0.0:
        # P = solve(Upsilon_H @ Upsilon, I, assume_a='her')
        P = inv_hermitian_matrix(Upsilon_H @ Upsilon)
        q = P @ Upsilon_H @ E

    else:
        P = solve(Upsilon_H @ Upsilon + lambda2 * I, I, assume_a='her')
        q = P @ Upsilon_H @ E
        P = P - lambda2 * P @ P

    return q, P

















