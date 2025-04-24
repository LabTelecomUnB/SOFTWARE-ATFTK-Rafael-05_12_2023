###############################
Math tools (``aftk.mathtools``)
###############################


This package is a collection of useful routines and high-level services involving mathematical calculations and data visualisation that supports AFTK ecossystem.


.. currentmodule:: aftk.mathtools

.. toctree:: 
    :hidden:
    :maxdepth: 1
    :caption: Mathtools

    aftk.mathtools.psdm.eigenvalues
    aftk.mathtools.psdm.condition_number
    aftk.mathtools.psdm.inverse
    aftk.mathtools.linreg.LinearRegression
    aftk.mathtools.spharm.Ylm
    aftk.mathtools.spharm.m_times_Ylm_over_sin_theta
    aftk.mathtools.spharm.dYlm_dtheta
    aftk.mathtools.spharm.Tlm
    aftk.mathtools.spharm.T
    aftk.mathtools.spharm.Upsilon


**************************************
PSD Matrices (``aftk.mathtools.psdm``)
**************************************

This module provides performance-optimised services regarding complex-valued positive-semidefinite matrices :math:`A^H A` when :math:`A\in\mathbb{C}^{m\times n}` is given.

In several situations one needs to analyse some mathematical properties of square matrices in the form :math:`A^H A \in\mathbb{C}^{n\times n}` where :math:`A\in\mathbb{C}^{m\times n}` is already known. In scenarios like those, the previous knowledge of :math:`A` can be used to increase the performance of the required calculations.

.. Recall that a positive-semidefinite matrix :math:`M` can always be  (refer to `Decomposition of PSDM <https://en.wikipedia.org/wiki/Definite_matrix#Decomposition>`_)

Recommended import
==================

.. code-block:: python

    from aftk.mathtools import psdm


Implemented Services
====================

.. autosummary:: 
    :nosignatures:

    psdm.eigenvalues
    psdm.condition_number
    psdm.inverse


*********************************************
Linear Regression (``aftk.mathtools.linreg``)
*********************************************

Recommended import
==================

.. code-block:: python

    from aftk.mathtools import linreg


Implemented Services
====================

.. autosummary:: 
    :nosignatures:

    linreg.LinearRegression     


***********************************************
Spherical Harmonics (``aftk.mathtools.spharm``)
***********************************************

Antenna fields are most commonly described in `spherical coordinates <https://en.wikipedia.org/wiki/Spherical_coordinate_system>`_ where they are usually projected in orthogonal set of functions involving `spherical harmonics <https://en.wikipedia.org/wiki/Spherical_harmonics>`_. In this context, this module provides services to quickly calculate those functions.


Convention: :math:`\theta \in \left[0, \pi\right]` and :math:`\phi \in \left[0, 2\pi\right)` are the polar and azimuthal angle, respectively.  I


Recommended import
==================

.. code-block:: python

    from aftk.mathtools import spharm


Implemented Services
====================

.. autosummary:: 
    :nosignatures:

    spharm.Ylm
    spharm.m_times_Ylm_over_sin_theta
    spharm.dYlm_dtheta
    spharm.Tlm
    spharm.T
    spharm.Upsilon
    

**************************************************
Spherical Triangulation (``aftk.mathtools.sptri``)
**************************************************

Recommended import
==================

.. code-block:: python

    from aftk.mathtools import sptri


Implemented Services
====================

.. autosummary:: 
    :nosignatures:

