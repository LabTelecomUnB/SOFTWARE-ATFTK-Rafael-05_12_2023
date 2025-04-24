##########
Quickstart
##########


*****************
The SMC Estimator
*****************

.. code-block:: python3
    :caption: Estimation SMC

    >>> from aftk.solvers import SMCEstimator

    >>> smce = SMCEstimator(lmax: int,
                            theta: float[:], 
                            phi: float[:], 
                            E_theta: complex[:], 
                            E_phi: complex[:], 
                            R: complex[:, :] | None = None,
                            W: complex[:, :] | None = None, 
                            method: str = 'blue',
                            tikhonov_matrix : complex[:, :] = zero,
                            unbiased_restriction: bool = False)

    # regression related:
    >>> smce.residual                     # complex[:]
    >>> smce.residual_energy              # float
    >>> smce.param_cov_matrix             # complex[:, :] 
    >>> smce.P                            # same as param_cov_matrix
    >>> smce.estimation_matrix            # complex[:, :]
    >>> smce.K                            # same as estimation_matrix
    >>> smce.regressor_matrix             # complex[:, :]
    >>> smce.U                            # same as regressor_matrix
    >>> smce.smc_matrix                   # q: complex[:, 1]
    >>> smce.q                            # same as smc_matrix
    >>> smce.bias_matrix                  # KU - I : complex[:, :]
    >>> smce.b                            # same as bias_matrix
    >>> smce.condition_number             # float, condition number of (U.H@W@U + Gamma.H@Gamma)
    >>> smce.lagrange_multiplier_matrix   # Lambda: complex[:, :]

    # identification:
    >>> smce.description
    
    # persisting
    >>> smce.save('')
    
    # reporting
    >>> smce.show_result_dashboard(**kwargs)            # shows a dashboard containing results of the regression
    >>> smce.save_result_report('path.pdf', **kwargs)   # 

    >>> from aftk import SMC 
    >>> smc = SMC(smce.q, freq: float|None)

    >>> smcs = [SMC(smce.q, f) for smce in zip(smce_list, freqs)]

    >>> ant = Antenna(smcs)
    >>> ant(freq).plot_directivity()



******************
The Antenna Object
******************

As the most central tool, AFTK provides :class:`Antenna` objects.

.. code-block:: python3
    :caption: Antenna attributes and services

    # Operational

    Antenna.position          # [x, y, z]
    Antenna.band_tx           # [fmin, fmax]
    Antenna.band_rx           # [fmin, fmax]
    Antenna.smc(freq: float)  # q: complex[:] 

    Antenna.radiation_intensity(freq: float, theta: float[:], phi: float[:]): float[:]  # U: float[:]
    Antenna.directivity(freq: float, theta: float[:], phi: float[:]): float[:]          # D: float[:]
    Antenna.radiation_field(float freq, float theta[:], float phi[:])                   # [E_theta[:], E_phi[:]]
    Antenna.electric_field(input_signal, r, theta, phi)


.. code-block:: python3
    :linenos:
    :caption: Antenna Object
    :name: some_name
    :emphasize-lines: 0
    
    >>> from aftk import Antenna


*******************
Creating an Antenna
*******************






========
From TLE
========


.. code-block::
    :linenos:
    :caption: Creating an antenna from radiation field data
    :name:
    :emphasize-lines: 0

    >>> from sopt import Satellite
    >>>
    >>> line1 = '1 08808U 76035A   22147.63373556  .00000101  00000+0  00000+0 0  9999'
    >>> line2 = '2 08808   8.3547 308.1412 0019716 286.4296  84.8969  0.99845237732488'
    >>> sat = Satellite(tle=[line1, line2])


===============
Real Satellites
===============

  >>> sat = Satellite.from_last_tle(name='cosmo-skymed 1')


==================
Virtual Satellites
==================

  >>> parameters = {
  ...     'sma': 7000,         # Semi-major axis (in km)
  ...     'ecc': 0.2,          # eccentricity
  ...     'raa_deg': 30.0,     # Right ascension of the ascending node (in degrees)
  ...     'inc_deg': 15.0,     # Orbit inclination (in degrees)
  ...     'arp_deg': 45.0,     # Argument of perigee (in degrees)
  ...     'tra_deg': 45.0,     # True anomaly (in degrees)
  ...     'epoch': '2021-03-18 14:20:00.000'
  ... }
  >>> sat = Satellite(**parameters)


*********************
Propagating the Orbit
*********************

  >>> from pandas import date_range
  >>>
  >>> time_range = date_range(start='2021-03-18', end='2021-03-21', freq='min')
  >>> sat.propagate(time_range, propagator='two-body')
  >>> sat.ground_track


*************************
Plotting the Ground track
*************************

  >>> from sopt import Map
  >>>
  >>> map = Map()
  >>> map.add(sat.ground_track)
  >>> map.show()
