#ifndef OUTPUT_H
#define OUTPUT_H
#include "pkd.h"

typedef enum {
    OUT_TINY_GROUP,
    } outType;

struct inOutputSend {
    int iPartner;     /* Who to send the data to */
    outType eOutputType;  /* What kind of output */
    };

void pkdOutputSend(PKD pkd, outType eOutputType, int iPartner);
void pkdOutput(PKD pkd, outType eOutputType, int iProcessor,int nProcessor,
    int iPartner,int nPartner, const char *fname );

#endif
