#  -*- coding: utf-8 -*-
"""

Author: Rafael R. L. Benevides
Date: 27/04/2023

"""
import numpy

from libc.math cimport sin, cos, sqrt

# ========== ========== ========== ========== ========== Ludwig I (cart)
cpdef tuple sph_from_cart(double[:] theta,
                          double[:] phi,
                          double complex[:] Ex,
                          double complex[:] Ey,
                          double complex[:] Ez):
    """
    Ludwig II (spherical) from Ludwig I
    
    Parameters
    ----------
    theta
    phi
    Ex
    Ey
    Ez

    Returns
    -------

    """

    cdef:
        Py_ssize_t index
        double complex[:] Eth = numpy.zeros_like(Ex)
        double complex[:] Eph = numpy.zeros_like(Ex)

        double cos_th, sin_th
        double cos_ph, sin_ph

    for index in range(theta.shape[0]):

        cos_th = cos(theta[index])
        sin_th = sin(theta[index])

        cos_ph = cos(phi[index])
        sin_ph = sin(phi[index])

        Eth[index] = cos_th*cos_ph*Ex[index] + cos_th*sin_ph*Ey[index] - sin_th*Ez[index]
        Eph[index] = -sin_th*Ex[index] + cos_ph*Ey[index]

    return numpy.asarray(Eth), numpy.asarray(Eph)


# ========== ========== ========== ========== ========== Ludwig II (azel)
cpdef tuple sph_from_azel(double[:] theta,
                          double[:] phi,
                          double complex[:] Eel,
                          double complex[:] Eaz):
    """
    Ludwig II (spherical) from Ludwig II (Az/El)
    
    In this definition poles are at y-axis, elevation is the polar angle and azimuth is 
    the longitudinal angle. 
    
    Parameters
    ----------
    theta
    phi
    Eel
    Eaz

    Returns
    -------

    """

    cdef:
        Py_ssize_t index
        double complex[:] Eth = numpy.zeros_like(Eaz)
        double complex[:] Eph = numpy.zeros_like(Eaz)

        double cos_th, sin_th
        double cos_ph, sin_ph

    for index in range(theta.shape[0]):

        cos_th = cos(theta[index])
        sin_th = sin(theta[index])

        cos_ph = cos(phi[index])
        sin_ph = sin(phi[index])

        Eth[index] = cos_ph*Eaz[index] + cos_th*sin_ph*Eel[index]
        Eph[index] = -cos_th*sin_ph*Eaz[index] + cos_ph*Eel[index]

    return numpy.asarray(Eth), numpy.asarray(Eph)


cpdef tuple azel_from_sph(double[:] theta,
                          double[:] phi,
                          double complex[:] Eth,
                          double complex[:] Eph):
    """
    Ludwig II (Az/El) from Ludwig II (spherical)
    
    In this definition poles are at y-axis, elevation is the polar angle and azimuth is 
    the longitudinal angle. 
    
    Parameters
    ----------
    theta
    phi
    Eth
    Eph

    Returns
    -------

    """

    cdef:
        Py_ssize_t index
        double complex[:] Eaz = numpy.zeros_like(Eth)
        double complex[:] Eel = numpy.zeros_like(Eth)

        double cos_th, sin_th
        double cos_ph, sin_ph
        double cos_el

    for index in range(theta.shape[0]):

        cos_th = cos(theta[index])
        sin_th = sin(theta[index])

        cos_ph = cos(phi[index])
        sin_ph = sin(phi[index])
        cos_el = sqrt(1 - sin_th*sin_th*sin_ph*sin_ph)

        Eaz[index] = (cos_ph*Eth[index] - cos_th*sin_ph*Eph[index]) / cos_el
        Eel[index] = (cos_th*sin_ph*Eth[index] + cos_ph*Eph[index]) / cos_el

    return numpy.asarray(Eel), numpy.asarray(Eaz)


# ========== ========== ========== ========== ========== Ludwig II (elza)
cpdef tuple sph_from_elaz(double[:] theta,
                          double[:] phi,
                          double complex[:] Eel,
                          double complex[:] Eaz):
    """
    Ludwig II (spherical) from Ludwig II (El/Az)

    In this definition poles are at x-axis, azimuth is the polar angle and elevation is 
    the longitudinal angle. 

    Parameters
    ----------
    theta
    phi
    Eel
    Eaz

    Returns
    -------

    """

    cdef:
        Py_ssize_t index
        double complex[:] Eth = numpy.zeros_like(Eaz)
        double complex[:] Eph = numpy.zeros_like(Eaz)

        double cos_th, sin_th
        double cos_ph, sin_ph

    for index in range(theta.shape[0]):

        cos_th = cos(theta[index])
        sin_th = sin(theta[index])

        cos_ph = cos(phi[index])
        sin_ph = sin(phi[index])

        Eth[index] = cos_ph*cos_th*Eaz[index] + sin_ph*Eel[index]
        Eph[index] = -sin_ph*Eaz[index] + cos_ph*cos_th*Eel[index]

    return numpy.asarray(Eth), numpy.asarray(Eph)


cpdef tuple elaz_from_sph(double[:] theta,
                          double[:] phi,
                          double complex[:] Eth,
                          double complex[:] Eph):
    """
    Ludwig II (El/Az) from Ludwig II (spherical)

    In this definition poles are at x-axis, azimuth is the polar angle and elevation is 
    the longitudinal angle. 
    
    Parameters
    ----------
    theta
    phi
    Eth
    Eph

    Returns
    -------

    """

    cdef:
        Py_ssize_t index
        double complex[:] Eaz = numpy.zeros_like(Eth)
        double complex[:] Eel = numpy.zeros_like(Eth)

        double cos_th, sin_th
        double cos_ph, sin_ph
        double cos_az

    for index in range(theta.shape[0]):

        cos_th = cos(theta[index])
        sin_th = sin(theta[index])

        cos_ph = cos(phi[index])
        sin_ph = sin(phi[index])

        cos_az = sqrt(1 - sin_th*sin_th*cos_ph*cos_ph)

        Eaz[index] = (cos_th*cos_ph*Eth[index] - sin_ph*Eph[index]) / cos_az
        Eel[index] = (sin_ph*Eth[index] + cos_th*cos_ph*Eph[index]) / cos_az

    return numpy.asarray(Eaz), numpy.asarray(Eel)


# ========== ========== ========== ========== ========== Ludwig III (xpol)
cpdef tuple sph_from_xpol(double[:] phi,
                          double complex[:] Eco,
                          double complex[:] Exp):
    """
    Ludwig II (spherical) from Ludwig III
    
    Parameters
    ----------
    phi : rad, longitudinal angle
    Eco : co-polar component
    Exp : cross-polar component

    Returns
    -------

    """

    cdef:
        Py_ssize_t index
        double complex[:] Eth = numpy.zeros_like(Eco)
        double complex[:] Eph = numpy.zeros_like(Eco)

        double cos_ph, sin_ph

    for index in range(phi.shape[0]):

        cos_ph = cos(phi[index])
        sin_ph = sin(phi[index])

        Eth[index] = +cos_ph*Eco[index] + sin_ph*Exp[index]
        Eph[index] = -sin_ph*Eco[index] + cos_ph*Exp[index]

    return numpy.asarray(Eth), numpy.asarray(Eph)


cpdef tuple xpol_from_sph(double[:] phi,
                          double complex[:] Eth,
                          double complex[:] Eph):
    """
    Ludwig III (spherical) from Ludwig II
    
    Parameters
    ----------
    phi
    Eth
    Eph

    Returns
    -------

    """

    cdef:
        Py_ssize_t index
        double complex[:] Eco = numpy.zeros_like(Eth)
        double complex[:] Exp = numpy.zeros_like(Eth)

        double cos_ph, sin_ph

    for index in range(phi.shape[0]):

        cos_ph = cos(phi[index])
        sin_ph = sin(phi[index])

        Eco[index] = cos_ph*Eth[index] - sin_ph*Eph[index]
        Exp[index] = sin_ph*Eth[index] + cos_ph*Eph[index]

    return numpy.asarray(Eco), numpy.asarray(Exp)
