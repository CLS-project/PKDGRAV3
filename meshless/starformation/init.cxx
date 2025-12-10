#include "init.h"
#include "master.h"

void MSR::StarFormInit(double dTime) {
    struct inStarFormInit in;
    struct outStarForm out;
    in.dTime = dTime;
    in.bFeedback = DoStellarFeedback();
    in.bCCSNFeedback = parameters.get_bCCSNFeedback();
    in.bSNIaFeedback = parameters.get_bSNIaFeedback();
    in.dCCSNFBDelay = calc.dCCSNFBDelay;
    in.dSNIaFBDelay = calc.dSNIaFBDelay;

    mdl->RunService<ServiceStarFormInit>(sizeof(in), &in, &out, sizeof(out));

    printf("%d star particles are about to explode in this IC\n", out.nFormed);
}


