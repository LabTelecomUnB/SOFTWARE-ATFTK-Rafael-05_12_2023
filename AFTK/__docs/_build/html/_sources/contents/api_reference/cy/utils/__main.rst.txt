##################################
C-level Utilities (``smet.utils``)
##################################



.. c:function:: numpy.complex128_t smet.utils.foo(double a[:, :], str scale)

   :param type: description of the first parameter.
   :param nitems: description of the second parameter.
   :returns: a result.
   :retval NULL: under some conditions.
   :retval NULL: under some other conditions as well.


.. c:struct:: Data

    .. c:union:: @data

        .. c:var:: int a

        .. c:var:: double b

Explicit :c:var:`Data.@data.a`. Short-hand :c:var:`Data.a`.
    

.. c:function:: PyObject *PyType_GenericAlloc(PyTypeObject *type, Py_ssize_t nitems)

   :param type: description of the first parameter.
   :param nitems: description of the second parameter.
   :returns: a result.
   :retval NULL: under some conditions.
   :retval NULL: under some other conditions as well.
