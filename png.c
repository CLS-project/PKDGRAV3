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

#ifdef HAVE_CONFIG_H
#include "config.h"
#else
#include "pkd_config.h"
#endif

#include <stdint.h>
#include <assert.h>
#include <malloc.h>
#include <math.h>

#include "png.h"

/* Famous white, red, blue, black colour map */
static int wrbb( float v ) {
    if ( v < 0.0 ) v = 0.0;
    else if ( v > 1.0 ) v = 1.0;

    v *= 6;

    if      ( v < 1 ) return gdTrueColor(0,0,(int)(255*v));
    else if ( v < 2 ) return gdTrueColor((int)(255*(v-1.0)),0,255);
    else if ( v < 3 ) return gdTrueColor(255,0,(int)(255*(3.0-v)));
    else if ( v < 5 ) return gdTrueColor(255,(int)(255*(v-3.0)*0.5),0);
    else if ( v < 6 ) return gdTrueColor(255,255,(int)(255*(v-5.0)));
    else              return gdTrueColor(255,255,255);
    }

void pngWrite( PNG png, FILE *fp, float *Density ) {
    float v, color_scale;
    int x, y;

    color_scale = 1.0 / (png->maxValue - png->minValue);
    for ( x=0; x<png->iResolution; x++ ) {
	for ( y=0; y<png->iResolution; y++ ) {
	    v = Density[x+png->iResolution*y];
	    if ( v < 1e-10 ) v = -10;
	    else v = log10f(v);
	    v = color_scale * ( v - png->minValue );
	    gdImageSetPixel(png->ic,x,y,wrbb(v));
	    }
	}
    gdImagePng( png->ic, fp );
    }

PNG pngInitialize( int iResolution, float minValue, float maxValue ) {
    PNG png;

    png = (PNG)malloc( sizeof(struct png) );
    assert( png !=NULL );

    png->iResolution = iResolution;
    png->minValue = minValue;
    png->maxValue = maxValue;
    png->ic = gdImageCreateTrueColor(png->iResolution,png->iResolution);

    return png;
    }

void pngFinish( PNG png ) {
    gdImageDestroy(png->ic);
    free(png);
    }
