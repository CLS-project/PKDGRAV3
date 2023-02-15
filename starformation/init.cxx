#include "starformation/starformation.h"
#include "master.h"

void MSR::StarFormInit(double dTime) {
    struct inStarFormInit in;
    struct outStarForm out;
    in.dTime = dTime;
#ifdef FEEDBACK
    in.dSNFBDelay = param.dSNFBDelay;
#else
    in.dSNFBDelay = -1.0;
#endif

    pstStarFormInit(pst, &in, sizeof(in), &out, sizeof(out));

    printf("%d star particles are about to explode in this IC\n", out.nFormed);
}

#ifdef __cplusplus
extern "C" {
#endif

int pstStarFormInit(PST pst,void *vin,int nIn,void *vout,int nOut) {
    struct inStarFormInit *in = (struct inStarFormInit *) vin;
    struct outStarForm *out = (struct outStarForm *) vout;
    int rID;

    mdlassert(pst->mdl,nIn == sizeof(struct inStarFormInit));
    if (pst->nLeaves > 1) {
        struct outStarForm fsStats;

        rID = mdlReqService(pst->mdl,pst->idUpper,PST_STARFORMINIT,in,nIn);
        pstStarFormInit(pst->pstLower,in,nIn,vout,nOut);
        mdlGetReply(pst->mdl,rID,&fsStats,NULL);
        out->nFormed += fsStats.nFormed;
    }
    else {
        pkdStarFormInit(pst->plcl->pkd, in->dTime, in->dSNFBDelay, &out->nFormed);
    }
    return sizeof(struct outStarForm);
}

void pkdStarFormInit(PKD pkd, double dTime, double dSNFBDelay, int *nFormed) {
    *nFormed = 0;
    for (auto &p : pkd->particles) {
        if (p.is_star()) {
            auto &star = p.star();
            if (star.fTimer >= 0) { // fTimer < 0 can be used for stars
                // in the IC that are not supossed to explode
#ifdef FEEDBACK
                if ( (dTime-star.fTimer) < dSNFBDelay)
#endif
                {
                    star.hasExploded = 0; // This particle has not exploded before
                    *nFormed += 1;
                }
            }
        }
    }
}
#ifdef __cplusplus
}
#endif
