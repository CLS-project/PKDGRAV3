#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdint.h>
#include <assert.h>
#include <malloc.h>
#include <math.h>

#include "png.h"

void wrbb( int *cmap, gdImagePtr im ) {
    double slope, offset;
    int i;

    // Shamefully stolen from Tipsy

    slope = 255./20. ;
    for(i = 0 ;i < 21 ;i++){
	cmap[i] = gdImageColorAllocate(im,0,0,(int)(slope * (double)(i) + .5)) ;
    }
    slope = 191./20. ;
    offset = 64. - slope * 21. ;
    for(i = 21 ;i < 42 ;i++){
	cmap[i] = gdImageColorAllocate(im,(int)(slope*(double)(i) + offset + .5),0,255);
    }
    slope = -205./21. ;
    offset = 255. - slope * 41. ;
    for(i = 42 ;i < 63 ;i++){
	cmap[i] = gdImageColorAllocate(im, 255,0,(int)(slope*(double)(i) + offset + .5));
    }
    slope = 205./40. ;
    offset = 50. - slope * 63. ;
    for(i = 63 ;i < 104 ;i++){
	cmap[i] = gdImageColorAllocate(im,255,(int)(slope*(double)(i) + offset + .5),0);
    }
    slope = 255./21. ;
    offset = -slope * 103. ;
    for(i = 104 ;i < 125 ;i++){
	cmap[i] = gdImageColorAllocate(im,255,255,(int)(slope*(double)(i) + offset +.5));
    }

    // The MARK color
    cmap[125] = gdImageColorAllocate(im,255,255,255);
    cmap[126] = gdBrushed;
}

static const int BWIDTH = 2; //!< Width of the brush

void pngWrite( PNG png, FILE *fp, float *Density )
{
    float v, color_slope;
    uint_fast32_t i;
    int x, y;

    color_slope = 124.0 / (png->maxValue - png->minValue);
    for( x=0; x<png->iResolution; x++ ) {
	for( y=0; y<png->iResolution; y++ ) {
	    int clr;
	    v = Density[x+png->iResolution*y];
	    if ( v < 1e-10 ) v = -10;
	    else v = log10f(v);
	    v = color_slope * ( v - png->minValue ) + 0.5;
	    clr = (int)(v);
	    if ( clr < 0 ) clr = 0;
	    if ( clr > 124 ) clr = 124;
	    assert( clr >= 0 && clr < 125 );
	    gdImageSetPixel(png->ic,x,y,png->cmap[clr]);
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

    png->brush = gdImageCreate(BWIDTH,BWIDTH);
    gdImageFilledRectangle(png->brush,0,0,BWIDTH-1,BWIDTH-1,
			   gdImageColorAllocate(png->brush,0,255,0));
    png->ic = gdImageCreate(png->iResolution,png->iResolution);
    wrbb(png->cmap,png->ic);

    return png;
}

void pngFinish( PNG png )
{
    gdImageDestroy(png->ic);
    gdImageDestroy(png->brush);
    free(png);
}
