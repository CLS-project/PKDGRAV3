============
Units
============

We adopt the convention that:

.. math::

    G &= 1\\
    H_0 &= \sqrt{\frac{8}{3}\pi}

where :math:`G` is the gravitational constant and :math:`H_0` is the Hubble parameter. Given the definition of
:math:`\rho_c` (the critical density of the universe):

.. math::

    \rho_{c} = \frac{3H_0^2}{8\pi G}

it follows that :math:`\rho_c` is :math:`1` is code units. We further adopt the convention that:

.. math::

    L = 1

where :math:`L` is the length of the periodic box centered on :math:`0,0,0`.
The coordinates x, y and z have the range :math:`\left[-0.5,0.5\right)`.
This also means that the volume of the box is :math:`1`,
and thus the total mass in the box must be :math:`\Omega_m` (matter density of the Universe).

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

It is helpfull to introduce the value of the gravitational constant.

.. math:: 
    G = 4.30091727063\times 10^{-9}~Mpc~M_\odot^{-1}(km/s)^2

Recall that 1 Mpc is :math:`\pi^{-1}\times 9.69394202136\times 10^{19}` km, so we can rewrite this as the product of a density and time squared.

.. math:: 
    G &= 4.30091727063\times 10^{-9}\frac{Mpc}{M_\odot}\left(\frac{km}{s}\right)^2 \times
        \left(\frac{\pi}{9.69394202136\times 10^{19}}\right)^2\left(\frac{Mpc}{km}\right)^2\\
      &= \frac{4.30091727063\times 10^{-9}\times\pi^2}{9.69394202136^2\times 10^{38}}
      \frac{Mpc^3}{M_\odot}\frac{1}{s^2}\\
      &= 4.51710305052\times 10^{-48}\frac{Mpc^3}{M_\odot}\frac{1}{s^2}

Also recall that:

.. math:: 
    H_0 = 100h \left[km/s/Mpc\right]

Length
------

The length unit is simply the width of the cosmological box in :math:`L` in :math:`h^{-1}` Mpc.

Mass
----

The mass unit is given by:

.. math::
    M_{box} &= L^3\left[h^{-3}Mpc^3\right]\times \rho_c\left[\frac{M_\odot}{Mpc^3}\right]\\
            &= \frac{L^3}{h^3}\times\frac{3H_0^2}{8\pi G}~M_\odot\\
            &= \frac{L^3}{h^3}\times\frac{3\times 100^2\times h^2}{8\pi G}~M_\odot\\
            &= L^3\times\frac{3\times 100^2\times h^2}{8\pi\times 4.30091727063\times 10^{-9}}~h^{-1}M_\odot\\
            &= L^3\times 2.77536627208\times 10^{11}~h^{-1}M_\odot

Since the volume unit is 1, the density unit is the same as the mass unit.

Velocity
--------

The :math:`G=1` criteria allows us to derive the velocity unit by factoring out the mass unit and the length unit.

.. math:: 
    \left(\frac{km}{s}\right)^2 &= \frac{G\times M_{box}}{L}\\
                 &= \frac{G\times L^3\times \rho_{c}}{L}\\
                 &= \frac{G\times L^3\times 3H_0^2 }{L\times 8\pi \times G}\\
                 &= \frac{3}{8\pi}\times L^2 \times H_0^2\\
    \frac{km}{s} &= \sqrt{\frac{3}{8\pi}}\times L \times H_0\\
                 &= \sqrt{\frac{3}{8\pi}}\times L \times 100h

Thus when :math:`L` is expressed in :math:`h^{-1}` Mpc:

.. math:: 
    V_{unit} = \frac{\sqrt{\frac{8\pi}{3}}}{L\times 100}

Or in proper (non-comoving) units:

.. math:: 
    V_{unit} = \frac{\sqrt{\frac{8\pi}{3}}}{L\times a\times 100}


Time
----

The :math:`G=1` criteria also allows us to derive the time unit by factoring out the density unit.

.. math:: 
    G &= 4.51710305052\times 10^{-48}\frac{Mpc^3}{M_\odot}\frac{1}{s^2}\\

Multiplying by our density unit gives:

.. math:: 
    \frac{1}{t_{unit}^2} &= G\times\rho_{unit}\\
        &= 4.51710305052\times 10^{-48}\frac{Mpc^3}{M_\odot}\frac{1}{s^2} \times 2.77536627208\times 10^{11}\frac{M_\odot}{Mpc^3}\\
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

    G &= 4.30091727063\times 10^{-9}~Mpc~M_\odot^{-1}(km/s)^2\\
      &= 4.30091727063\times 10^{-6}~kpc~M_\odot^{-1}(km/s)^2\\

Thus, for :math:`G=1`, the mass unit must be:

.. math:: 
    M_{unit} = \frac{1}{4.30091727063\times 10^{-6}} = 2.32508541103\times 10^5~M_\odot

The length unit is kpc, so the density unit is:

.. math:: 
    p_{unit} = 2.32508541103\times 10^5\frac{M_\odot}{kpc^3}


Similarily the time unit must be:

.. math:: 
    \frac{1}{t_{unit}^2} &= G\times\rho_{unit}\\
        &= 4.51710305052\times 10^{-48}\frac{Mpc^3}{M_\odot}\frac{1}{s^2} \times 2.32508541103\times 10^5\frac{M_\odot}{kpc^3}\\
        &= 4.51710305052\times 10^{-48}\frac{Mpc^3}{M_\odot}\frac{1}{s^2} \times 2.32508541103\times 10^{14}\frac{M_\odot}{Mpc^3}\\
        &= 1.05026504029\times 10^{-33}\frac{1}{s^2}\\
    \frac{1}{t_{unit}} &= \sqrt{1.05026504029\times 10^{-33}\frac{1}{s^2}}\\
      &= 3.24077928944\times 10^{-17}\frac{1}{s}\\
    t_{unit} &= 3.08567758149\times 10^{16} s\\
             &\approx 0.978461942~Gyrs

