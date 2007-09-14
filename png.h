#include <gd.h>
#include <gdfontg.h>

typedef struct png
{
    int iResolution;
    float minValue;
    float maxValue;

    gdImagePtr ic;
    gdImagePtr brush;

    int cmap[128];
} *PNG;

PNG pngInitialize( int iResolution, float minValue, float maxValue );
