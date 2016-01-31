#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "output.h"
#include "pkd.h"


/* Generic context: all things/stuff from iIndex (starts at zero) */
struct packCtx {
    PKD pkd;
    int iIndex;
    };





void pkdOutput(PKD pkd,int eOutputType, int iPartner,int nPartner) {
    struct packCtx ctx = {pkd,0};


    /* I do all of the writing */
    if (iPartner == pkd->idSelf) {
	}

    /* We just send all of our data onward */
    else {
	}




    }

