#include "starformation/starformation.h"
#include "master.h"

void MSR::StarFormInit(double dTime) {
#ifdef FEEDBACK
    struct inStarFormInit in;
    struct outStarForm out;
    in.dTime = dTime;
    in.bCCSNFeedback = parameters.get_bCCSNFeedback();
    in.bSNIaFeedback = parameters.get_bSNIaFeedback();
    in.dCCSNFBDelay = calc.dCCSNFBDelay;
    in.dSNIaFBDelay = calc.dSNIaFBDelay;

    pstStarFormInit(pst, &in, sizeof(in), &out, sizeof(out));

    printf("%d star particles are about to explode in this IC\n", out.nFormed);
#endif
}

int pstStarFormInit(PST pst,void *vin,int nIn,void *vout,int nOut) {
#ifdef FEEDBACK
    struct inStarFormInit *in = (struct inStarFormInit *) vin;
    struct outStarForm *out = (struct outStarForm *) vout;
    int rID;

    mdlassert(pst->mdl,nIn == sizeof(struct inStarFormInit));
    if (pst->nLeaves > 1) {
        struct outStarForm fsStats;

        rID = mdlReqService(pst->mdl,pst->idUpper,PST_STARFORMINIT,in,nIn);
        pstStarFormInit(pst->pstLower,in,nIn,vout,nOut);
        pst->mdl->GetReply(rID,fsStats);
        out->nFormed += fsStats.nFormed;
    }
    else {
        pkdStarFormInit(pst->plcl->pkd, *in, &out->nFormed);
    }
#endif
    return sizeof(struct outStarForm);
}

void pkdStarFormInit(PKD pkd, struct inStarFormInit in, int *nFormed) {
#ifdef FEEDBACK
    *nFormed = 0;
    for (auto &p : pkd->particles) {
        if (p.is_star()) {
            auto &star = p.star();
            if (star.fTimer >= 0) { // fTimer < 0 can be used for stars
                // in the IC that are not supossed to explode
                if ((in.dTime - star.fTimer) < in.dCCSNFBDelay)
                    star.bCCSNFBDone = in.bCCSNFeedback ? 0 : 1;

                if ((in.dTime - star.fTimer) < in.dSNIaFBDelay)
                    star.bSNIaFBDone = in.bSNIaFeedback ? 0 : 1;

                if (!(star.bCCSNFBDone && star.bSNIaFBDone))
                    *nFormed += 1;
            }
        }
    }
#endif
}
