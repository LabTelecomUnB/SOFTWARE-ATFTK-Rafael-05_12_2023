##########################
Utilities (``smet.utils``)
##########################


.. toctree::
    :hidden:
    :maxdepth: 1
    :caption: Exceptions

    smet.utils.Time
    smet.utils.AlgorithmError
    smet.utils.DTypeError
    smet.utils.ReadOnlyError
    smet.utils.SignatureError


.. currentmodule:: smet.utils


This package comprehends several tools on top of which the project relies.


****
Time
****


.. autosummary::
    :nosignatures:
    
    Time

*****************
Custom Exceptions
*****************

.. autosummary::
    :nosignatures:
    
    AlgorithmError
    SignatureError
    DTypeError
    ReadOnlyError


***********
Convenience
***********

.. toctree::
    :hidden:
    :maxdepth: 1
    :caption: Convenience

    smet.utils.check_real_dtype
    smet.utils.combine
    smet.utils.stream_combinations
    smet.utils.compare
    smet.utils.check_types

.. autosummary::
    :nosignatures:
    
    check_real_dtype
    check_types
    combine
    stream_combinations
    compare



****************
Quantity Related
****************

.. toctree::
    :hidden:
    :maxdepth: 1
    :caption: Quantity Related
    
    smet.utils.reduce_angle
    smet.utils.parse_angle
    smet.utils.parse_angular_speed
    smet.utils.parse_eccentricity
    smet.utils.parse_inclination
    smet.utils.parse_dimensionless
    smet.utils.parse_length
    smet.utils.parse_speed
    smet.utils.parse_timespan
    smet.utils.parse_distance
    smet.utils.parse_latitude
    smet.utils.parse_longitude
    smet.utils.parse_time
    smet.utils.parse_types


Parsing by physical type
========================

.. autosummary::
    :nosignatures:

    parse_dimensionless
    parse_timespan
    parse_angle
    parse_angular_speed
    parse_length
    parse_speed


Parsing by value
================

.. autosummary::
    :nosignatures:

    reduce_angle
    parse_latitude
    parse_longitude
    parse_distance
    parse_eccentricity
    parse_inclination

Auxiliary Parsers
=================

.. autosummary::
    :nosignatures:

    parse_types

************
Base Classes
************

.. toctree::
    :hidden:
    :maxdepth: 1
    :caption: Base Classes

    smet.utils.Identifiable
    smet.utils.Taggable
    smet.utils.Nameable
    smet.utils.Describable
    smet.utils.Representable
    smet.utils.Activatable
    smet.utils.Serialisable
    smet.utils.Storable
    smet.utils.Propagable

.. autosummary::
    :nosignatures:

    Identifiable
    Taggable
    Nameable
    Describable
    Representable
    Serialisable
    Storable


