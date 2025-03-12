============
File Formats
============

------------
Tipsy Format
------------

The Tipsy format is a binary format for storing particle data. It consists of a
fixed length header followed by particle data. The data can either be stored
in big-endian ("standard") or "native" (usually little-endian) format. The header is as follows:

.. code-block:: none

    struct TipsyHeader {
        double   dTime;     // time of snapshot
        uint32_t nBodies;   // number of particles 
        uint32_t nDim;      // normally 3
        uint32_t nSph;      // number of gas particles
        uint32_t nDark;     // number of dark matter particles
        uint32_t nStar;     // number of star particles
        uint32_t nPad;      // padding (normally 0)
    };

The header is followed by the particle data. The number of particles is given by
``nBodies``. Immediately following the header is the particle data. The data is
ordered with all of the gas (sph) particles first, followed by the dark matter
particles, and finally the star particles. The number of particles of each type
is given by ``nSph``, ``nDark``, and ``nStar``, the sum of which should equal
``nBodies``.

This differs from the original Tipsy format in that the original Tipsy format
used a signed integer for the number of particles which limited the number of
particles to :math:`2^{31} - 1`.

The pkdgrav3 code uses an other extension to the Tipsy format that uses the ``nPad``
field to extend the particles counts by eight bits each. This allows for
particle counts up to :math:`2^{40} - 1`. The ``nPad`` field is used as follows:

.. code-block:: none

    header.nPad  = ((N>>32)&0xff)
                 + ((nSph>>24)&0xff00)
                 + ((nDark>>16)&0xff0000)
                 + ((nStar>>8)&0xff000000);

If this extension is not use then the ``nPad`` field must be set to zero.

The format of the gas (sph) particles is as follows:

.. code-block:: none

    struct tipsySph {
        float mass;         // mass of particle
        float pos[3];       // position of particle
        float vel[3];       // velocity of particle
        float rho;          // density of particle
        float temp;         // temperature of particle
        float hsmooth;      // smoothing length of particle
        float metals;       // metallicity of particle
        float phi;          // gravitational potential of particle
    };

The format of the dark particles is as follows:

.. code-block:: none

    struct tipsyDark {
        float mass;         // mass of particle
        float pos[3];       // position of particle
        float vel[3];       // velocity of particle
        float eps;          // softening length of particle
        float phi;          // gravitational potential of particle
    };

The format of the star particles is as follows:

.. code-block:: none

    struct tipsyStar {
        float mass;         // mass of particle
        float pos[3];       // position of particle
        float vel[3];       // velocity of particle
        float metals;       // metallicity of particle
        float tform;        // formation time of particle
        float eps;          // softening length of particle
        float phi;          // gravitational potential of particle
    };

All fields use code units (see the section on :ref:`units <units>`).


--------------
Gadget2 Format
--------------

The pkdgrav3 code can also read the Gadget2 format.
Information on this format can be found in the Gadget2 documentation.

See: https://wwwmpa.mpa-garching.mpg.de/gadget/users-guide.pdf


--------------
HDF5 Format
--------------

The HDF5 format has a structure similar to what is found in codes such as
Gadget2 and its descendants. Following is a description of each group that
is stored and their contents.

The ``Header`` group contains general information about the simulation and
the snapshot.

.. code-block:: none

    Header {
        BoxSize                // Size of the box in Mpc/h
        Cosmological run       // Whether or not this is a cosmological run
        Flag_Cooling           // Whether or not the code was compiled with the COOLING flag
        Flag_DoublePrecision   // Whether or not position and velocities are stored in double precision
        Flag_Feedback          // Whether or not the code was compiled with the FEEDBACK flag
        Flag_Metals            // Whether or not gas and star particles contain a metallicity dataset
        Flag_Sfr               // Whether or not the code was compiled with the STAR_FORMATION flag
        Flag_StellarAge        // Same as Flag_Sfr above
        Hubble0                // Hubble constant in code units
        HubbleParam            // Dimensionless Hubble constant, h, in units of 100 km/(s Mpc)
        MassTable              // UNUSED. It is just an array of zeros
        NumFilesPerSnapshot    // The number of files the snapshot has been split into
        NumPart_ThisFile       // Number of particles of each PartType in the present file
        NumPart_Total          // Total number of particles of each PartType
        Omega0                 // Present-day total-matter parameter
        OmegaB                 // Present-day baryonic-matter parameter
        OmegaLambda            // Present-day dark-energy parameter
        PKDGRAV version        // Code's version
        Redshift               // The snapshot's redshift
        Softening              // Comoving softening length in code units
        Time                   // The snapshot's cosmic time in code units
    }

The ``Cosmology`` group contains the cosmological parameters. They are the same as
those included in the ``Header`` group. This is done for compatibility reasons
with existing simulation analysis libraries.

.. code-block:: none

    Cosmology {
        Cosmological run       // Whether or not this is a cosmological run
        Hubble0                // Hubble constant in code units
        HubbleParam            // Dimensionless Hubble constant, h, in units of 100 km/(s Mpc)
        Omega_b                // Present-day baryonic-matter parameter
        Omega_lambda           // Present-day dark-energy parameter
        Omega_m                // Present-day total-matter parameter
    }

The ``Units`` group contains unit conversion factors from code units to any of the
units below. These factors should be multiplied to values read from the snapshot to
obtain the desired units.

.. code-block:: none

    Units {
        ErgPerGmUnit           // Energy per unit mass to erg/g
        ErgUnit                // Energy to erg
        GasConst               // Boltzmann constant over hydrogen mass in code units
        GmPerCcUnit            // Mass density to g/cm^3
        KmPerSecUnit           // Velocity to km/s
        KpcUnit                // Length to kpc
        MsolUnit               // Mass to solar mass
        SecUnit                // Time to seconds
    }

Following the Gadget convention, six particle types are allowed in the snapshots.
They map to pkdgrav3 particle species as follows.

* ``PartType0``: Gas
* ``PartType1``: Dark matter
* ``PartType2``: UNUSED
* ``PartType3``: UNUSED
* ``PartType4``: Star
* ``PartType5``: Black hole

A group is created for every particle type present in the simulation. In general,
the datasets contained in each of them will depend on the compile- and run-time
parameters that were used in running the simulation. Following is the minimum
contents each particle group has.

.. code-block:: none

    PartType0 {
        Coordinates            // Coordinates (comoving in cosmological runs)
        ParticleIDs            // Particle ID
        Potential              // Gravitational potential
        Velocities             // Velocity (comoving in cosmological runs; i.e., dx/dt)
    }

