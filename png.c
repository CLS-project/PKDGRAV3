#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdint.h>
#include <assert.h>
#include <malloc.h>
#include <math.h>

#include "png.h"

/* Famous white, red, blue, black colour map */
static int wrbb( float v )
{
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

void pngWrite( PNG png, FILE *fp, float *Density )
{
    float v, color_scale;
    uint_fast32_t i;
    int x, y;

    color_scale = 1.0 / (png->maxValue - png->minValue);
    for( x=0; x<png->iResolution; x++ ) {
	for( y=0; y<png->iResolution; y++ ) {
	    v = Density[x+png->iResolution*y];
	    if ( v < 1e-10 ) v = -10;
	    else v = log10f(v);
	    v = color_scale * ( v - png->minValue );
	    gdImageSetPixel(png->ic,x,y,wrbb(v));
	}
    }
    gdImagePng( png->ic, fp );
}

PNG pngInitialize( int iResolution, float minValue, float maxValue )
{
    PNG png;

    png = (PNG)malloc( sizeof(struct png) );
    assert( png !=NULL );

    png->iResolution = iResolution;
    png->minValue = minValue;
    png->maxValue = maxValue;
    png->ic = gdImageCreateTrueColor(png->iResolution,png->iResolution);

    return png;
}

void pngFinish( PNG png )
{
    gdImageDestroy(png->ic);
    free(png);
}
