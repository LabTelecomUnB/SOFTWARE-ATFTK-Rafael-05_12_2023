.. _api-typing:

###########################
Typehints (``smet.typing``)
###########################


.. currentmodule:: smet.typing

This module provides the most important typehints defined and used in SMET.

***************
Time Like Types
***************

.. todo:: sorry!


*******************
Quantity Like Types
*******************


The following aliases represent types that can be coerced into (can be understood as) :py:class:`quantities <astropy.units.Quantity>`. In short, those types are numbers (scalars), numerical data (arrays, list, etc) or even quantities itself.


.. note::

   The aliases ``_ArrayLike<X>_co`` used in the following definitions are imported from ``numpy._typing._array_like`` and denote array-like objects that can be coerced into ndarrays with <X> dtype. 


Generic Physical Type
=====================


.. autodata:: QuantityLike
.. autodata:: RealQuantityLike

.. .. data:: smet.typing.QuantityLike
    :type: TypeAlias
    :value: _ArrayLikeNumber_co | Quantity

    Alias representing types that can be coerced into generic quantities.


.. .. data:: smet.typing.RealQuantityLike
    :type: TypeAlias
    :value: _ArrayLikeFloat_co | Quantity

    Alias representing types that can be coerced into real-valued generic quantities.


Length
======


.. data:: smet.typing.LengthLike
    :type: TypeAlias
    :value: _ArrayLikeNumber_co | Quantity['length']

    Alias representing types that can be coerced into length quantities.


.. data:: smet.typing.RealLengthLike
    :type: TypeAlias
    :value: _ArrayLikeFloat_co | Quantity['length']

    Alias representing types that can be coerced into real-valued length quantities.


Angle
=====


.. data:: smet.typing.AngleLike
    :type: TypeAlias
    :value: _ArrayLikeNumber_co | Quantity['angle']

    Alias representing types that can be coerced into angle quantities.


.. data:: smet.typing.RealAngleLike
    :type: TypeAlias
    :value: _ArrayLikeFloat_co | Quantity['angle']

    Alias representing types that can be coerced into real-valued angle quantities.






















******
Legacy
******

.. autodata:: NumericType
