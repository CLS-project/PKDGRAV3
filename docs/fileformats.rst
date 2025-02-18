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

