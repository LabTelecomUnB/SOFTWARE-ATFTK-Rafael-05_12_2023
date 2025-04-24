##########
Satellites
##########



Satellites are certainly the central objects in :mod:`smet` workflow. The most interesting services SMET provides depend on such objects at some level. The goal of this tutorial is to present the usage and main aspects of Satellite objects.  


***************************
Creating Satellites Objects
***************************

All satellites are instances of the :class:`Satellite <smet.satellite.Satellite>` class, which can be directly imported from :mod:`smet` package as follows.

.. code-block:: python
    :linenos:
    :name: import_satellite_class

    >>> from smet import Satellite

    >>> sat = Satellite()
    >>> sat
    =======================================
                   SATELLITE               
    ---------------------------------------
    NAME:           None
    DESCRIPTION:    None
    TAG:            None
    ---------------------------------------
    None
    ---------------------------------------
    ID:    22a281625d894784a5133517ca60d7d8
    =======================================




******************
Propagating Orbits
******************




***************
Adding Payloads
***************



