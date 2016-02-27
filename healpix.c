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

#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include "healpix.h"

static const double twothird=2.0/3.0;
static const double twopi=6.283185307179586476925286766559005768394;
static const double inv_halfpi=0.6366197723675813430755350534900574;

/*! Returns the remainder of the division \a v1/v2.
    The result is non-negative.
    \a v1 can be positive or negative; \a v2 must be positive. */
static double fmodulo (double v1, double v2)
  {
  if (v1>=0)
    return (v1<v2) ? v1 : fmod(v1,v2);
  double tmp=fmod(v1,v2)+v2;
  return (tmp==v2) ? 0. : tmp;
/*  return (v1>=0) ? ((v1<v2) ? v1 : fmod(v1,v2)) : (fmod(v1,v2)+v2); */
  }

static hpint64 imodulo64 (hpint64 v1, hpint64 v2)
  { hpint64 v=v1%v2; return (v>=0) ? v : v+v2; }

hpint64 nside2npix64(hpint64 nside)
  { return 12*nside*nside; }

static hpint64 ang2pix_ring_z_phi64 (hpint64 nside_, double z, double s,
  double phi)
  {
  double za = fabs(z);
  double tt = fmodulo(phi,twopi) * inv_halfpi; /* in [0,4) */

  if (za<=twothird) /* Equatorial region */
    {
    double temp1 = nside_*(0.5+tt);
    double temp2 = nside_*z*0.75;
    hpint64 jp = (hpint64)(temp1-temp2); /* index of  ascending edge line */
    hpint64 jm = (hpint64)(temp1+temp2); /* index of descending edge line */

    /* ring number counted from z=2/3 */
    hpint64 ir = nside_ + 1 + jp - jm; /* in {1,2n+1} */
    int kshift = 1-(ir&1); /* kshift=1 if ir even, 0 otherwise */

    hpint64 ip = (jp+jm-nside_+kshift+1)/2; /* in {0,4n-1} */
    ip = imodulo64(ip,4*nside_);

    return nside_*(nside_-1)*2 + (ir-1)*4*nside_ + ip;
    }
  else  /* North & South polar caps */
    {
    double tp = tt-(int)(tt);
    double tmp = (s>-2.) ? nside_*s/sqrt((1.+za)/3.) : nside_*sqrt(3*(1-za));

    hpint64 jp = (hpint64)(tp*tmp); /* increasing edge line index */
    hpint64 jm = (hpint64)((1.0-tp)*tmp); /* decreasing edge line index */

    hpint64 ir = jp+jm+1; /* ring number counted from the closest pole */
    hpint64 ip = (hpint64)(tt*ir); /* in {0,4*ir-1} */
    ip = imodulo64(ip,4*ir);

    if (z>0)
      return 2*ir*(ir-1) + ip;
    else
      return 12*nside_*nside_ - 2*ir*(ir+1) + ip;
    }
  }

int64_t vec2pix_ring64(hpint64 nside, const double *vec)
  {
  double vlen=sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
  double cth = vec[2]/vlen;
  double sth=(fabs(cth)>0.99) ? sqrt(vec[0]*vec[0]+vec[1]*vec[1])/vlen : -5;
  return ang2pix_ring_z_phi64 (nside,cth,sth,atan2(vec[1],vec[0]));
  }
