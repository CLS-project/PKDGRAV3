/* -----------------------------------------------------------------------------
 * SUBSET OF THE ORIGINAL FUNCTIONS. ONLY THOSE REQUIRED BY PKDGRAV.
 * -----------------------------------------------------------------------------
 *
 *  Copyright (C) 1997-2012 Krzysztof M. Gorski, Eric Hivon, Martin Reinecke,
 *                          Benjamin D. Wandelt, Anthony J. Banday,
 *                          Matthias Bartelmann,
 *                          Reza Ansari & Kenneth M. Ganga
 *
 *
 *  This file is part of HEALPix.
 *
 *  HEALPix is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  HEALPix is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with HEALPix; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  For more information about HEALPix see http://healpix.sourceforge.net
 *
 *---------------------------------------------------------------------------*/

#ifndef HEALPIX_H
#define HEALPIX_H
#include <stdint.h>

typedef int64_t hpint64;

hpint64 nside2npix64(hpint64 nside);
hpint64 ang2pix_ring_z_phi64(hpint64 nside_,double z,double s,double phi);
int64_t vec2pix_ring64(hpint64 nside, const double *vec);

#endif
