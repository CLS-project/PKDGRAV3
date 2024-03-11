=======
Changes
=======

**The format of the checkpoint files has been updated.**

Before, the restart parameters write written to checkpoint.chk.par.
The old behaviour of ommitting the ".par" for the checkpoint
has been reintroduced. In addition, all of the state variables
are now saved in Python Dill (Pickle) format in a ".chk.pkl" file.
To restart, use the ".chk" file as input, for example::

  srun pkdgrav3 checkpoint.00010.chk

To override parameters you can add them as arguments to the ``restore`` call.
For example if you wanted to turn on parallel reading and update
the number of readers to 100 you would change::

  msr.restore(__file__)

to::

  msr.restore(__file__,bParaRead=True,nParaRead=100)

The new values of those parameters will persist on subsequent restarts.

**The parameter processing has undergone significant changes.**

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
