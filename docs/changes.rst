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

achOutTimes
  This parameter has been removed. The functionality can be duplicated by using the list
  format for nSteps and dRedTo. For example, if you have a list of redshift to output,
  you can set the number of steps to be one for each interval::

    dRedFrom = 49               # Start at z=49
    dRedTo = [10,2,1,0.5,0]     # Step to z=10, 2, 1, 0.5 and 0
    nSteps = [1] * len(dRedTo)  # Taking one step for each interval
    iOutInterval = 1            # Output every interval

  If you have expansion factors you can just convert them to redshift::

    dRedTo = [1 / a - 1 for a in [0.1,0.5,1] ]

nStepsSync,dRedSync
  These parameter have been removed. The same functionality is now more generally
  available with nSteps and dRedTo. For example, before you might have had::

    dRedFrom = 49
    dRedTo = 0
    nSteps = 160
    nStepsSync = 60
    dRedSync = 10

  This would take 60 steps to redshift 10, then 100 steps to redshift 0.
  The new way of specifying this is::

    dRedFrom = 40
    dRedTo = [10,0]
    nSteps = [60,100]