========
Analysis
========

.. load_parameters:: ../parameters.toml

--------------------------
Power Spectrum Measurement
--------------------------

The code can measure the power spectrum P(k) of a periodic simulation.
It uses grid interlacing (see :cite:t:`HockneyEastwood1988`, section 7-8 "Interlacing")
and high order mass assignment (see :cite:t:`2016MNRAS.460.3624S`).

Parameters
==========

The following parameters control P(k) measurement:

.. show_parameters:: nGridPk iPkInterval iPkOrder bPkInterlace

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

----------
Light Cone
----------

The code can output light cone data as either healpix maps, or raw particle data (or both).

Parameters
==========

.. show_parameters:: bLightCone dRedshiftLCP nSideHealpix bLightConeParticles bBowtie sqdegLCP hLCP


Output
======

Healpix
=======

The healpix output is a binary format consisting of 32-bit integer counts of grouped and ungrouped particles
(grouped is only valid if group finding is enabled) as well as the potential. The following script
(found in tools/hpb2fits.py) will convert the file to fits format.

.. literalinclude:: ../tools/hpb2fits.py

Particles
---------

The particle output is also a binary format. Note that files may be empty if the assigned processor does not
output any particles during the step. This is normal. Each output particle is 40 bytes as follows.

.. tabularcolumns:: |r|l|l|
.. csv-table::    
   :header: Offset, Type, Field

   0,  64-bit integer, particle ID
   8,  3 x float,      position (x y z)
   20, 3 x float,      velocity (x y z)
   32, float,          potential
   36, 4 bytes,        padding
