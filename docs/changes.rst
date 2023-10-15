=======
Changes
=======

The parameter processing has undergone significant changes.

.. index::
   single: parameters ; dTheta
   single: parameters ; dTheta2
   single: parameters ; dTheta20

:ref:`dTheta <dTheta>`, dTheta2, dTheta20
  Only dTheta remains. To retain the old behaviour see the :ref:`accuracy extension <accuracy-module>`.


.. index::
   single: parameters ; dxPeriod
   single: parameters ; dyPeriod
   single: parameters ; dzPeriod

dxPeriod, dyPeriod, dzPeriod
  These parameters have been removed.
  Instead, use :ref:`dPeriod <dPeriod>` and set it to a list of three values.

.. index::
   single: parameters ; hxLCP
   single: parameters ; hyLCP
   single: parameters ; hzLCP

hxLCP,hyLCP,hzLCP
  These parameters have been replaced with :ref:`hLCP <hLCP>`, a vector of three values.
