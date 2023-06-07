========
Analysis
========

--------------------------
Power Spectrum Measurement
--------------------------

The code can measure the power spectrum P(k) of a periodic simulation.
It uses grid interlacing and high order mass assignment as described
in Sefusatti et. al [#sefusatti]_

Parameters
==========

The following parameters control P(k) measurement:

.. index:: single: parameters ; nGridPk

nGridPk (default 0 or off)
    This is the width of the 3D mass assignment grid used to measure P(k).
    A good value for this is the same as the simulation IC grid size
    give or take a factor of two::

        nGridPk = nGrid // 2

.. index:: single: parameters ; iPkInterval

iPkInterval (default 1)
    Measure and output P(k) every step. You can change this to have P(k)
    measurements done less frequently. For example, setting iPkInterval to 10
    will cause P(k) to be measured every 10 steps (10, 20, 30, etc.).

.. index:: single: parameters ; iPkOrder

iPkOrder (default 4)
    This selects the mass assignment scheme. Possible values are:

    1. Nearest Grid Point (NGP)
    2. Cloud in Cell (CIC)
    3. Triangular Shaped Cloud (TSC)
    4. Piecewise Cubic Spline (PCS)

    It is recommended to keep the default (PCS).

.. index:: single: parameters ; bPkInterlace

bPkInterlace (default True)
    This controls the use of interlaced grids to reduce aliasing. If enabled
    (the default), it creates two mass assignment grids and deposits the mass
    on the second grid offset by half a grid spacing. This greatly reduces the
    aliasing effect at the cost of additional memory. We recommend that this be left on.

Output
======

The code outputs a ``pk`` file for each selected step with four columns:

k
    The :math:`k` bin in :math:`h\text{Mpc}^{-1}`. If ``dBoxSize`` has been specified
    (often the case when using the build-in IC generator), :math:`k` will be in
    physical units (scale by ``dBoxSize``). Otherwise you need to divide by your box length.

P(k)
    The measured P(k) in :math:`h^{-1}\text{Mpc}^{3}`. Again, if ``dBoxSize`` was not
    specified then this will need to be corrected by multiplying by
    :math:`\texttt{dBoxSize}^3`.

Number of Samples
    This is the number of ``Delta(k)`` cells that went into the P(k) measurement for each bin.

Linear P(k)
    This is the same as the P(k) column, except it includes the linear part (if present),
    for example if you are using the linear neutrino treatment.



.. rubric:: Footnotes

.. [#sefusatti] https://arxiv.org/abs/1512.07295

