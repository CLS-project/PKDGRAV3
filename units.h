/*  This file is part of PKDGRAV3 (http://www.pkdgrav.org/).
 *  Copyright (c) 2001-2018 Joachim Stadel & Douglas Potter
 *
 *  PKDGRAV3 is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  PKDGRAV3 is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with PKDGRAV3.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef UNITS_HINCLUDED
#define UNITS_HINCLUDED

#define KBOLTZ	1.38e-16     /* bolzman constant in cgs */
#define MHYDR 1.67e-24       /* mass of hydrogen atom in grams */
#define MSOLG 1.99e33        /* solar mass in grams */
#define GCGS 6.67e-8         /* G in cgs */
#define KPCCM 3.085678e21    /* kiloparsec in centimeters */
#define SIGMAT 6.6524e-25    /* Thompson cross-section (cm^2) */
#define LIGHTSPEED 2.9979e10 /* Speed of Light cm/s */
#define SECONDSPERYEAR   31557600.

typedef struct units{
   double dMsolUnit;
   double dKpcUnit;
   double dKBoltzUnit;
   double dGmPerCcUnit;
   double dComovingGmPerCcUnit;
   double dErgPerGmUnit;
   double dErgUnit;
   double dSecUnit;
   double dKmPerSecUnit;
   double dGasConst;
} UNITS;

#endif
