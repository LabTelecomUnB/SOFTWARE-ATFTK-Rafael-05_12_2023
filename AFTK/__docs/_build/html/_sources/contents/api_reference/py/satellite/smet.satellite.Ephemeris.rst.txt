#########
Ephemeris
#########

.. currentmodule:: smet.satellite

.. autoclass:: Ephemeris
    :show-inheritance:



*******
Summary
*******

==========
Attributes
==========


Epoch
-----

    .. autosummary::

        ~Ephemeris.epoch


Anomalies
---------

    .. autosummary::

        ~Ephemeris.mean_anomaly
        ~Ephemeris.eccentric_anomaly
        ~Ephemeris.true_anomaly


Orbital Parameters
------------------

    .. autosummary::

        ~Ephemeris.semi_major_axis
        ~Ephemeris.eccentricity
        ~Ephemeris.inclination
        ~Ephemeris.right_ascension
        ~Ephemeris.argument_of_perigee


Auxiliary Parameters
--------------------

    .. autosummary::

        ~Ephemeris.perigee
        ~Ephemeris.apogee
        ~Ephemeris.orbital_period
        ~Ephemeris.mean_motion


Aliases
-------

    .. autosummary::

        ~Ephemeris.sma
        ~Ephemeris.ecc
        ~Ephemeris.inc
        ~Ephemeris.raa
        ~Ephemeris.arp
        ~Ephemeris.mea
        ~Ephemeris.eca
        ~Ephemeris.tra



=======
Methods
=======


Parsers
-------

    .. autosummary::

        ~Ephemeris.parse_eccentricity



Conversions
-----------

Involving true, eccentric and mean anomalies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    .. autosummary::

        ~Ephemeris.eccentric_anomaly_from_mean_anomaly
        ~Ephemeris.eccentric_anomaly_from_true_anomaly
        ~Ephemeris.mean_anomaly_from_eccentric_anomaly
        ~Ephemeris.mean_anomaly_from_true_anomaly
        ~Ephemeris.true_anomaly_from_eccentric_anomaly
        ~Ephemeris.true_anomaly_from_mean_anomaly


Involving perigee, apogee, eccentricity and semi-major axis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    .. autosummary::

        ~Ephemeris.sma_from_perigee_and_apogee
        ~Ephemeris.ecc_from_perigee_and_apogee
        ~Ephemeris.perigee_from_sma_and_ecc
        ~Ephemeris.apogee_from_sma_and_ecc


Involving semi-major axis, orbital period and mean angular motion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    .. autosummary::

        ~Ephemeris.mean_motion_from_orbital_period
        ~Ephemeris.mean_motion_from_sma
        ~Ephemeris.orbital_period_from_mean_motion
        ~Ephemeris.orbital_period_from_sma
        ~Ephemeris.sma_from_mean_motion
        ~Ephemeris.sma_from_orbital_period



*************
Documentation
*************

    .. autoattribute:: smet.satellite.Ephemeris.epoch

    .. autoattribute:: smet.satellite.Ephemeris.semi_major_axis
    .. autoattribute:: smet.satellite.Ephemeris.sma

    .. autoattribute:: smet.satellite.Ephemeris.eccentricity
    .. autoattribute:: smet.satellite.Ephemeris.ecc

    .. autoattribute:: smet.satellite.Ephemeris.inclination
    .. autoattribute:: smet.satellite.Ephemeris.inc

    .. autoattribute:: smet.satellite.Ephemeris.right_ascension
    .. autoattribute:: smet.satellite.Ephemeris.raa

    .. autoattribute:: smet.satellite.Ephemeris.argument_of_perigee
    .. autoattribute:: smet.satellite.Ephemeris.arp

    .. autoattribute:: smet.satellite.Ephemeris.mean_anomaly
    .. autoattribute:: smet.satellite.Ephemeris.mea

    .. autoattribute:: smet.satellite.Ephemeris.eccentric_anomaly
    .. autoattribute:: smet.satellite.Ephemeris.eca

    .. autoattribute:: smet.satellite.Ephemeris.true_anomaly
    .. autoattribute:: smet.satellite.Ephemeris.tra

    .. autoattribute:: smet.satellite.Ephemeris.perigee
    .. autoattribute:: smet.satellite.Ephemeris.apogee
    .. autoattribute:: smet.satellite.Ephemeris.orbital_period
    .. autoattribute:: smet.satellite.Ephemeris.mean_motion

    .. automethod:: smet.satellite.Ephemeris.parse_eccentricity

    .. automethod:: smet.satellite.Ephemeris.eccentric_anomaly_from_mean_anomaly
    .. automethod:: smet.satellite.Ephemeris.eccentric_anomaly_from_true_anomaly
    .. automethod:: smet.satellite.Ephemeris.mean_anomaly_from_eccentric_anomaly
    .. automethod:: smet.satellite.Ephemeris.mean_anomaly_from_true_anomaly
    .. automethod:: smet.satellite.Ephemeris.true_anomaly_from_eccentric_anomaly
    .. automethod:: smet.satellite.Ephemeris.true_anomaly_from_mean_anomaly

    .. automethod:: smet.satellite.Ephemeris.perigee_from_sma_and_ecc
    .. automethod:: smet.satellite.Ephemeris.apogee_from_sma_and_ecc
    .. automethod:: smet.satellite.Ephemeris.sma_from_perigee_and_apogee
    .. automethod:: smet.satellite.Ephemeris.ecc_from_perigee_and_apogee

    .. automethod:: smet.satellite.Ephemeris.mean_motion_from_orbital_period
    .. automethod:: smet.satellite.Ephemeris.mean_motion_from_sma
    .. automethod:: smet.satellite.Ephemeris.orbital_period_from_mean_motion
    .. automethod:: smet.satellite.Ephemeris.orbital_period_from_sma
    .. automethod:: smet.satellite.Ephemeris.sma_from_mean_motion
    .. automethod:: smet.satellite.Ephemeris.sma_from_orbital_period
