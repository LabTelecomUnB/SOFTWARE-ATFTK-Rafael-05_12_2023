.. todo::
    We are sorry. This page is under construction!

###########
Development
###########


.. code-block:: python
    :linenos:
    :caption: Sensor data
    :name:
    :emphasize-lines: 0

    {
        "__class__": "FieldOfRegard",
        "id": "36a84sd36asd6as",
        "tag": "sensorA",
        "name": "Sensor A",
        "description": "Range of possibilities for sensor A",
        "frame": "ITRS",
        "ellipsoid": "WGS84",
        "aperture": 35*u.deg,
        "time": Time(...),  # (N,) like instance of astropy.time.Time
        "footprint": Quantity(..., u.deg)  # (M, 2, N) angle Quantity
        "alt": Quantity(..., u.km),  # (N,) length Quantity
        "active": True,
        "visible": True,
        "facecolor": "#012345",
        "edgecolor": "#012345",
        "linestyle": "solid",
        "linewidth": 1.0,
        "alpha": 1.0
    }



.. code-block:: python
    :linenos:
    :caption: Ephemeris data
    :name:
    :emphasize-lines: 0

    {
        "__class__": "TLE",
        "frame": "TEME"
        "line1": "...",
        "line2": "..."
    }

    or

    {
        "__class__": "KeplerianParameters",
        "epoch": "2022-03-18 14:20:35",
        "frame": "GCRS",
        "sma": 7000*u.km,
        "ecc": 0.05,
        "inc": 40.0*u.deg,
        "raa": 15.0*u.deg,
        "arp": 23.0*u.deg,
        "tra": 71.0*u.deg
    }


.. code-block:: python
    :linenos:
    :caption: Satellite data
    :name:
    :emphasize-lines: 0

    {
        "__class__": "Satellite",
        "id": "36a84sd36asd6as",
        "tag": "sgdc1",
        "name": "SGDC-1",
        "description": "Geostationary Defence and Strategic Communications Satellite",
        "ephemerides": [...], # or {} with epoch as key
        "propagator": "SGP4", # "TWOBODY", "J2", ...
        "time": Time(...),  # (N,) like instance of astropy.time.Time
        "coords_ECI": {
            "__class__": "CartesianCoords",
            "time": Time(...),  # (N,) like instance of astropy.time.Time
            "frame": "TEME",
            "r": Quantity(..., u.km),  # (N, 3) length Quantity
            "v": Quantity(..., u.km/u.s),  # (N, 3) length Quantity
        },
        "coords_ECEF": {
            "__class__": "CartesianCoords",
            "time": Time(...),  # (N,) like instance of astropy.time.Time
            "frame": "ITRS",
            "r": Quantity(..., u.km),  # (N, 3) length Quantity
            "v": Quantity(..., u.km/u.s),  # (N, 3) length Quantity
        },
        "coords_geo": {
            "__class__": "GeodesicCoords",
            "time": Time(...),  # (N,) like instance of astropy.time.Time
            "frame": "ITRS",
            "ellipsoid": "WGS84",
            "lon": Quantity(..., u.deg),  # (N,) length Quantity
            "lat": Quantity(..., u.deg),  # (N,) length Quantity
            "alt": Quantity(..., u.km),  # (N,) length Quantity
            "v": Quantity(..., u.km/u.s),  # (N, 3) length Quantity (east, north, up),
        }
        "groundtrack": {
            "__class__": "GroundTrack",
            "future": 30*u.min, # or None
            "past": 30*u.min, # or None
            "visible": True,
            "facecolor": "#012345",
            "edgecolor": "#012345",
            "linestyle": "solid",
            "linewidth": 1.0,
            "alpha": 1.0
        }
        "active": True,
        "visible": True,
        "facecolor": "#012345",
        "edgecolor": "#012345",
        "linestyle": "solid",
        "linewidth": 1.0,
        "alpha": 1.0,
        "markersize": 10,
        "children": {...} # sensors
    }
