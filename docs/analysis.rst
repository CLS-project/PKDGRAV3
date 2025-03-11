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

------------------------------
Friend of Friends Group Finder
------------------------------

The code can output FoF groups at each step.
The following parameters control the FoF Group Finder:

.. show_parameters:: bFindGroups dTau nMinMembers dEnvironment0 dEnvironment1

The output is a binary format with the following structure.

.. tabularcolumns:: |r|l|l|
.. csv-table::
   :header: Offset, Type, Field

    0, float[3], Position (x y z) of deepest potential
    12, float, Value of deepest potential
    16, float[3], Shrinking Sphere Center (x y z) (see :cite:t:`2003MNRAS.338...14P`)
    28, float[3], Center of Mass Position (x y z)
    40, float[3], Center of Mass Velocity
    52, float[3], Angular momentum
    64, float[6], Moment of Inertia
    88, float, Velocity Dispersion :math:`\sigma`
    92, float, :math:`r_{\text{max}}` (maximum distance from center of mass)
    96, float, Mass
    100, float, Mass enclosed in dEnvironment0
    104, float, Mass enclosed in dEnvironment1
    108, float, Half mass radius
    112, int, Number of black holes
    116, int, Number of stars
    120, int, Number of gas particles
    124, int, Number of dark matter particles
    128, uint64_t, Global group ID
    136

----------------------
General Analysis Hooks
----------------------

Aside from the legacy analysis functions, it is possible to add custom analysis hooks.
You can use :func:`PKDGRAV.add_analysis` to add any object instance that is callable by Python.
The callable object will be passed the following named parameters.

.. tabularcolumns:: |l|l|
.. csv-table::
   :header: Name, Description

   msr, the PKDGRAV module
   step, current integer step number
   time, simulation time
   a, the expansion factor (for cosmological simulations)
   theta, currently calculated theta


For example::

    from PKDGRAV import add_analysis, ASSIGNMENT_ORDER

    class MassGrid:
        grid = 0
        order = ASSIGNMENT_ORDER.PCS
        def __init__(self,grid,order=ASSIGNMENT_ORDER.PCS):
            self.grid = grid
            self.order = order
        def __call__(self,msr,step,time,**kwargs):
            print('calculating density grid')
            msr.grid_create(self.grid)
            msr.assign_mass(order=self.order)
            msr.grid_write('output.{:05d}'.format(step))
            msr.grid_delete()
        def ephemeral(self,msr,**kwargs):
            return msr.grid_ephemeral(self.grid)


This class can be used to output a density grid at each step.
To enable it, use::

    add_analysis(MassGrid(nGrid))



