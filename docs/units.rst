============
Units
============

---------------
Code Convention
---------------

We adopt the convention that:

.. math::

    G &= 1\\
    H_0 &= \sqrt{\frac{8}{3}\pi}

where :math:`G` is the gravitational constant and :math:`H_0` is the Hubble parameter. Given the definition of
:math:`\rho_c` (the critical density of the universe), it follows that:

.. math::

    \rho_{c} = \frac{3H_0^2}{8\pi G} = 1

------------------
Cosmological Units
------------------

For cosmological runs we adopt the usual comoving coordinates and periodic boundary conditions.
We also adopt the convention:

.. math::

    L = 1

where :math:`L` is the length of the periodic box centered on :math:`0,0,0`.
The coordinates x, y and z have the range :math:`\left[-0.5,0.5\right)`.
This also means that the volume of the box is :math:`1`,
and thus the total matter (dark matter and baryons) in the box must be:

.. math::

    M_{box} = \Omega_m\times\rho_c\times V = \Omega_m

where :math:`\Omega_m` is the matter density of the Universe.

----------------------------
Conversion to Physical Units
----------------------------

Summary
=======

+----------+-----------------------------------------------------------+--------------------------------+
| Unit     | Multiply by                                               | Units                          |
+==========+===========================================================+================================+
| Length   | L                                                         | comoving :math:`h^{-1}` Mpc    |
+----------+-----------------------------------------------------------+--------------------------------+
| Mass     | :math:`L^3\times 2.77536627208\times 10^{11}`             | :math:`h^{-1}M_\odot`          |
+----------+-----------------------------------------------------------+--------------------------------+
| Velocity | :math:`\sqrt{\frac{8\pi}{3}}\div({L\times a\times 100})`  | km/s                           |
+----------+-----------------------------------------------------------+--------------------------------+

Where:

:math:`L`
  The physical length of the box in comoving :math:`h^{-1}` Mpc.

:math:`a`
  The expansion factor.

Derivation
==========

It is helpful to introduce the value of the gravitational constant.

.. math:: 
    G = 4.30091727063\times 10^{-9}~\text{Mpc}~M_\odot^{-1}\left(\text{km}/\text{s}\right)^2

Recall that :math:`1` Mpc is :math:`\pi^{-1}\times 9.69394202136\times 10^{19}` km, so we can rewrite this as the product of a density and time squared.

.. math:: 
    G &= 4.30091727063\times 10^{-9}\frac{Mpc}{M_\odot}\left(\frac{\text{km}}{\text{s}}\right)^2 \times
        \left(\frac{\pi}{9.69394202136\times 10^{19}}\right)^2\left(\frac{\text{Mpc}}{\text{km}}\right)^2\\
      &= \frac{4.30091727063\times 10^{-9}\times\pi^2}{9.69394202136^2\times 10^{38}}
      \frac{\text{Mpc}^3}{M_\odot}\frac{1}{s^2}\\
      &= 4.51710305052\times 10^{-48}\frac{\text{Mpc}^3}{M_\odot}\frac{1}{s^2}

Also recall that:

.. math:: 
    H_0 = 100h \left[\text{km}/\text{s}/\text{Mpc}\right]

Length
------

The length unit is simply the width of the cosmological box in :math:`L` in :math:`h^{-1}` Mpc.

Mass
----

The density unit is given by:

.. math::
    \rho_{unit} &= \rho_c\left[\frac{M_\odot}{\text{Mpc}^3}\right]\\
            &= \frac{3H_0^2}{8\pi G}\left[\frac{M_\odot}{\text{Mpc}^3}\right]\\
            &= \frac{3\times 100^2\times h^2}{8\pi G}\left[\frac{M_\odot}{\text{Mpc}^3}\right]\\
            &= \frac{3\times 100^2\times h^2}{8\pi \times 4.30091727063\times 10^{-9}}\left[\frac{M_\odot}{\text{Mpc}^3}\right]\\
            &= 2.77536627208\times 10^{11}~\left[\frac{h^2 M_\odot}{\text{Mpc}^3}\right]


The mass unit is given by:

.. math::
    M_{box} &= L^3\left[h^{-3}\text{Mpc}^3\right]\times \rho_c\left[\frac{h^2 M_\odot}{\text{Mpc}^3}\right]\\
            &= L^3\times 2.77536627208\times 10^{11}~h^{-1}M_\odot

Velocity
--------

The :math:`G=1` criteria allows us to derive the velocity unit by factoring out the mass unit and the length unit.

.. math:: 
    \left(\frac{\text{km}}{\text{s}}\right)^2 &=
        G\times \frac{M_{box}}{L}\\
        &= G\times \frac{M_{box}}{L}\times\rho_c\frac{L^3}{M_{box}}\\
        &= G\times L^2\times\rho_c\\
        &= G\times L^2\times \frac{3H_0^2}{8\pi G}\\
        &= \frac{3}{8\pi}\times H_0^2\times L^2\\
    \frac{\text{km}}{\text{s}} &= \sqrt{\frac{3}{8\pi}}\times H_0 \times L\\
        &= \sqrt{\frac{3}{8\pi}}\times 100h \times L\\

Thus when :math:`L` is expressed in :math:`h^{-1}` Mpc:

.. math:: 
    V_{unit} = \frac{\sqrt{\frac{8}{3}\pi}}{L\times 100}

Or in proper (non-comoving) units:

.. math:: 
    V_{unit} = \frac{\sqrt{\frac{8}{3}\pi}}{L\times a\times 100}


Time
----

The :math:`G=1` criteria also allows us to derive the time unit by factoring out the density unit.

.. math:: 
    G &= 4.51710305052\times 10^{-48}\frac{\text{Mpc}^3}{M_\odot}\frac{1}{s^2}\\

Multiplying by our density unit gives:

.. math:: 
    \frac{1}{t_{unit}^2} &= G\times\rho_{unit}\\
        &= 4.51710305052\times 10^{-48}\frac{\text{Mpc}^3}{M_\odot}\frac{1}{s^2} \times 2.77536627208\times 10^{11}\frac{M_\odot}{\text{Mpc}^3}\\
      &= 1.25366154539\times 10^{-36}\frac{1}{s^2}\\
    \frac{1}{t_{unit}} &= \sqrt{1.25366154539\times 10^{-36}\frac{1}{s^2}}\\
      &= 1.1196702842\times 10^{-18}\frac{1}{s}\\
    t_{unit} &= 8.9312006765\times 10^{17} s

------------------
Other Unit Systems
------------------

kpc & km/s
==========

It is often convenient to fix the length unit to be kpc, and the velocity unit to be km/s. With :math:`G=1` as before,
we can calculate the mass unit:

.. math:: 

    G &= 4.30091727063\times 10^{-9}~\text{Mpc}~M_\odot^{-1}(\text{km}/\text{s})^2\\
      &= 4.30091727063\times 10^{-6}~\text{kpc}~M_\odot^{-1}(\text{km}/\text{s})^2\\

Thus, for :math:`G=1`, the mass unit must be:

.. math:: 
    M_{unit} = \frac{1}{4.30091727063\times 10^{-6}} = 2.32508541103\times 10^5~M_\odot

The length unit is kpc, so the density unit is:

.. math:: 
    \rho_{unit} = 2.32508541103\times 10^5\frac{M_\odot}{\text{kpc}^3}


Similarily the time unit must be:

.. math:: 
    \frac{1}{t_{unit}^2} &= G\times\rho_{unit}\\
        &= 4.51710305052\times 10^{-48}\frac{\text{Mpc}^3}{M_\odot}\frac{1}{\text{s}^2} \times 2.32508541103\times 10^5\frac{M_\odot}{\text{kpc}^3}\\
        &= 4.51710305052\times 10^{-48}\frac{\text{Mpc}^3}{M_\odot}\frac{1}{\text{s}^2} \times 2.32508541103\times 10^{14}\frac{M_\odot}{\text{Mpc}^3}\\
        &= 1.05026504029\times 10^{-33}\frac{1}{\text{s}^2}\\
    \frac{1}{t_{unit}} &= \sqrt{1.05026504029\times 10^{-33}\frac{1}{\text{s}^2}}\\
      &= 3.24077928944\times 10^{-17}\frac{1}{\text{s}}\\
    t_{unit} &= 3.08567758149\times 10^{16}~\text{s}\\
             &\approx 0.978461942~\text{Gyrs}

kpc & Gyrs
==========

One can also fix the length unit to kpc, the time unit to Gyrs and thus the velocity unit to kpc/Gyr. The mass unit is then:


.. math:: 
    G &= 4.51710305052\times 10^{-48}\frac{Mpc^3}{M_\odot}\frac{1}{s^2}\\
      &= 4.51710305052\times 10^{-39}\frac{kpc^3}{M_\odot}\frac{1}{s^2}\\

With :math:`3.15581184\times 10^7` seconds in a year and :math:`3.15581184\times 10^{16}` seconds in a Gyr, we get:

.. math:: 
    G &= 4.51710305052\times 10^{-39}\frac{kpc^3}{M_\odot}\frac{1}{s^2} \times \left(3.15581184\times 10^{16}\frac{s}{Gyr}\right)^2\\
      &= 4.49864994804\times 10^{-6}\frac{kpc^3}{M_\odot Gyr}

This means that the mass unit is:

.. math:: 
    M_{unit} = \frac{1}{4.49864994804\times 10^{-6}}~M_\odot = 2.22288911462\times 10^5~M_\odot

Velocity
--------

It should be noted that kpc/Gyr is very nearly km/s. This should be obvious because the time unit in the kpc and km/s unit system is very nearly 1 Gyr (but not exactly).

.. math:: 
    V = \frac{kpc}{Gyr}
      = \frac{\pi^{-1}\times 9.69394202136\times 10^{16}\frac{km}{kpc}}{3.15581184\times 10^{16}\frac{\text{s}}{\text{Gyr}}}
      = 0.97777615965~\text{km}/\text{s}

