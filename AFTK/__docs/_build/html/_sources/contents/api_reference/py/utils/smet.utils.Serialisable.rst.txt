``Serialisable``
================

.. currentmodule:: smet.utils

.. autoclass:: Serialisable


Typehints
---------

.. data:: Serialisable.SubType
    :value: TypeVar('SubType', bound=Serialisable)

    Typehint representing instances of :class:`Serialisable` subclasses.


.. data:: Serialisable.PrimitiveSType
    :value: 
        Union[None, str, Number, ndarray, Quantity, Timestamp, DatetimeIndex]
    :annotation: TypeAlias

    Alias representing the primitive types whose instances are allowed as
    Attributes in :class:`Serialisable` subclasses.

.. data:: Serialisable.SType
    :value: PrimitiveStype | SubType | list[SType] | dict[str, SType]
    :annotation: TypeAlias

    **Serialisable Type**: alias representing all types whose instances are allowed as serialisable attributes in :class:`Serialisable` subclasses.

.. data:: Serialisable.RType
    :value: PrimitiveStype | SubType | list[SType] | dict[str, SType]
    :annotation: TypeAlias

    **Retrievable Type**:Alias representing all types whose instances can be unserialised into in :class:`Serialisable` subclasses.

Methods
-------
      
.. automethod:: Serialisable._get_serialisable_property_names
.. automethod:: Serialisable._get_copiable_property_names
.. automethod:: Serialisable.serialise
.. automethod:: Serialisable.unserialise
.. automethod:: Serialisable._copy


Examples
--------

Point and Triangle
^^^^^^^^^^^^^^^^^^

The first and simplest example regards

.. code-block::
    :linenos:
    :name: serialisable-code-subclassing
    :caption: Position: a subclass of Serialisable

    >>> class Point(Serialisable):
    ...
    ...    @classmethod
    ...    def _get_serialisable_property_names(cls):
    ...        return ['x', 'y']
    ...
    ...    @classmethod
    ...    def _get_copiable_property_names(cls):
    ...        return ['x', 'y']
    ...
    >>> p0 = Point()
    >>> Serialisable.serialise(p0)
    {'x': None, 'y': None}
    >>> p1 = Point(x=1)
    >>> p2 = Point(x=4, y=3)


