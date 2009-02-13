#ifndef PNG_H
#define PNG_H
#define PNG_H_MODULE_ID "$Id$"

#include <gd.h>
#include <gdfontg.h>

typedef struct png {
    int iResolution;
    float minValue;
    float maxValue;
    gdImagePtr ic;
    } *PNG;

PNG pngInitialize( int iResolution, float minValue, float maxValue );
void pngWrite( PNG png, FILE *fp, float *Density );
void pngFinish( PNG png );
#endif
