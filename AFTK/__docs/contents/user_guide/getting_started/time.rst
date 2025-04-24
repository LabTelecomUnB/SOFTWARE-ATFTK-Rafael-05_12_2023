.. _tutorial-time:

############
Time Objects
############


In order to deal with different time scales and formats, SMET provides its own time type implemented by the :class:`Time <smet.utils.Time>` class.


****************************************************
Constructing :class:`Time <smet.utils.Time>` objects
****************************************************


.. list-table::
    :widths: 10 90
    :header-rows: 1

    *   -   Name
        -   Description

    *   -   ``jd1``     
        -   First part of two-part julian date (See [1]_). If it is an
            array, must be 0 or 1-dimensional and if it is a
            :py:class:`quantity <astropy.units.Quantity>` must be of
            time :py:class:`physical type
            <astropy.units.PhysicalType>` and unit must be days.

            Acceptable types:
            
            *   :data:`NumericType <smet.typing.NumericType>`

    *   -   ``jd2``
        -   Second part of two-part julian date. Always matching the
            type and shape of ``jd1``. :data:`NumericType <smet.typing.NumericType>`

    *   -   ``jd12``
        -   Concatenation of ``jd1`` and ``jd2``.

            If ``jd1`` and ``jd2`` are scalars, then:

            .. code-block:: python

                jd12 = numpy.hstack([jd1, jd2])

            with shape ``(2,)``.

            Else if ``jd1`` and ``jd2`` are arrays, then:

            .. code-block:: python

                jd12 = numpy.hstack([
                    jd1.reshape(-1, 1),
                    jd2.reshape(-1, 1)
                ])

            with shape ``(N, 2)``.

            :data:`NumericType <smet.typing.NumericType>`

    *   -   ``jd``
        -   :data:`NumericType <smet.typing.NumericType>`
            Julian date: :python:`jd = jd1 + jd2`

    *   -   ``mjd``
        -   :data:`NumericType <smet.typing.NumericType>`
            Modified julian date: :python:`mjd = jd - 2_400_000.5`

    *   -   ``byear``
        -   :data:`NumericType <smet.typing.NumericType>`
            Fractional years from Besselian Epoch

    *   -   ``jyear``
        -   :data:`NumericType <smet.typing.NumericType>`
            Fractional years from Julian Epoch

    *   -   ``ptp``
        -   `array_like <https://numpy.org/doc/stable/glossary.html#term-array_like>`_ of 
            :py:data:`int64 <numpy.int64>`
            `Precision Time Protocol <https://en.wikipedia.org/wiki/Precision_Time_Protocol>`_ 
            in nanoseconds. Always associated to the TAI scale.

From two-part julian date
=========================


1.  If ``jd1`` and ``jd2`` are available, then ``t`` shall be constructed by one of the following signatures:

    >>> t = Time(jd1, jd2)  
 
    >>> t = Time(jd1=jd1, jd2=jd2)

    where the first is the only case where 2 positional arguments are accepted.

2.  In the other hand, if ``jd12`` is the one available, then the following are equivalent

    >>> t = Time(jd12)
    
    >>> t = Time(jd12=jd12)


From (modified) julian date
===========================


1.  If ``jd`` is available, then:
    
    >>> t = Time(jd=jd)

1.  If ``mjd`` is available, then:
    
    >>> t = Time(mjd=mjd)


(Modified) julian date, julian and besselian epoch 
==================================================

Whenever ``t`` must be created from ``jd``, ``mjd``, ``jyear``, ``byear`` or ``ptp``, then a keyword must be provided as follows:

    >>> t = Time(jd=jd)

    >>> t = Time(mjd=mjd)

    >>> t = Time(jyear=jyear)

    >>> t = Time(byear=byear)

    >>> t = Time(ptp=ptp)  # if scale is set, then it is ignored and TAI is assumed.


***************************************************************
Differences from :py:class:`Astropy's Time <astropy.time.Time>`
***************************************************************


.. **Difficulty with Astropy** :py:class:`Time <astropy.time.Time>`

This class may be seen as an operational revisit of the great :py:class:`Astropy's Time <astropy.time.Time>` class, by which it was inspired (even though it does not extend it). In fact, the original idea behind SMET design was to use :py:class:`Astropy's Time <astropy.time.Time>` as the main time objects used by SMET functions (quite much in the same way it is still done with :py:class:`Quantity <astropy.units.Quantity>`).

.. math:: X(e^{j\omega }) = x(n)e^{ - j\omega n}

