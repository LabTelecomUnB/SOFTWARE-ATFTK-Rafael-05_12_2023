.. AFTK documentation master file, created by
   sphinx-quickstart on Fri Apr 21 08:07:25 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

##################
AFTK Documentation
##################

.. toctree::
   :hidden:
   :maxdepth: 1

   Services Overview <contents/introduction/__main>
   contents/user_guide/__main
   contents/api_reference/__main
   contents/gallery/__main
   contents/development/__main


AFTK, which stands for `Antenna Fields Tool Kit`, is an open-source MIT-licensed `Python <https://www.python.org>`_ package providing a whole ecosystem of services regarding the electromagnetic fields of antennas.


************
Installation
************

.. code-block::

    pip install aftk

*****************
Where do I start?
*****************

.. .. image:: _static/logo/logo_light.png
   :scale: 30%
   :align: center
   :class: only-light

.. .. image:: _static/logo/logo_dark.png
   :scale: 30%
   :align: center
   :class: only-dark

.. grid:: 2

    .. grid-item-card:: 
        :img-top: _static/logo/index_getting_started.svg


        Getting Started
        ^^^^^^^^^^^^^^^


        .. div:: sd-text-justify

            New to AFTK? Take a quick look at the main functionalities of the package and learn what kind of problems it solves.


        +++

        .. button-ref:: contents/introduction/__main
            :expand:
            :color: primary
            :click-parent:

            To the services overview

    .. grid-item-card::
        :img-top: _static/logo/index_user_guide.svg

        User Guide
        ^^^^^^^^^^

        .. div:: sd-text-justify

            Learn how to use AFTK! The user guide provides instructions on how to install AFTK and basic tutorials for all AFTK services.

        +++

        .. button-ref:: contents/user_guide/__main
            :expand:
            :color: primary
            :click-parent:

            To the user guide

    .. grid-item-card::
        :img-top: _static/logo/index_api.svg

        API Reference
        ^^^^^^^^^^^^^
        
        .. div:: sd-text-justify:

            The reference guide contains a detailed description of the entities (classes and functions) implemented by AFTK.

        +++

        .. button-ref:: contents/api_reference/__main
            :expand:
            :color: primary
            :click-parent:

            To the API reference

    .. grid-item-card::
        :img-top: _static/logo/index_contribute.svg

        Gallery 
        ^^^^^^^
        
        .. div:: sd-text-justify

            Need to see some inspirational examples? Access the source code and results of some analysis performed using AFTK.

        +++

        .. button-ref:: contents/gallery/__main
            :expand:
            :color: primary
            :click-parent:

            To the Gallery
