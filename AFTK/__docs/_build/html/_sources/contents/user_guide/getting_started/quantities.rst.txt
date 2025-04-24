#########################
Quantities and Timestamps
#########################

Needless to say, SMET needs to deal with several physical quantities as positions, velocities, angle, angular speeds, etc. Therefore, most of its numerical inputs (usually provided by the user) and outputs are required to come with physical units in tow to preserve their meaning. Although providing such units is not actually a problem, specifying them every single time will certainly make the user interface not so user-friendly.

To avoid a scenario like that, one might argue that the best development strategy is to standardise the tools to work only with a specific unit system, like `SI <https://docs.astropy.org/en/stable/units/index.html#module-astropy.units.si>`_ for example. Even though this sounds like a very good idea, it surely has some drawbacks regarding the fact that some of SI units are not so common in the context of space missions. In fact, try to imagine how inconvenient it would be always having to write the semi-major axis of a satellite orbit in metres instead of kilometres; its orbital period in seconds instead of minutes or hours; its mean motion in hertz instead revolutions per day; or its speed in metres per second instead of kilometres per hour or per second. 

Most certainly, a nice interface for dealing with physical quantities would not only consider their numerical values, but also their units and transforms. Fortunately, the `unit package <https://docs.astropy.org/en/stable/units/index.html>`_ of Astropy already provides such interface through `Quantity <https://docs.astropy.org/en/stable/api/astropy.units.Quantity.html>`_ objects, which hence makes them perfect for SMET purposes. This tutorial presents how SMET tools handle physical quantities.

.. warning:: 
    For a tutorial on how to use Astropy's units and quantities, please refer to `Units and Quantities <https://docs.astropy.org/en/stable/units/index.html>`_.


******************
Parsing Quantities
******************

Internally, all physical quantities (dimensionless included) in SMET are instances of `Quantity <https://docs.astropy.org/en/stable/api/astropy.units.Quantity.html>`_, hence it can be said that no unit system is actually preferred. In short, it is worth remembering the following rules:

1. 


.. code-block:: python
    :linenos:
    
    >>> from astropy import units
    
    >>> x = 10.0 * units.km
    
    

.. code-block:: python
    :linenos:
    :caption: How SMET services internally parses quantities
    :emphasize-lines: 5
    
    from smet.utils import parse_length

    def func(..., length_quantity, ...):
        ...
        length_quantity = parse_length(length_quantity)
        ...


**********
Timestamps
**********

******************
Parsing Quantities
******************


