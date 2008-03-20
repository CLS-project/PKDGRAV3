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

