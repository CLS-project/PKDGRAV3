#include "starformation/starformation.h"
#include "master.h"

void MSR::StarFormInit(double dTime) {
#ifdef FEEDBACK
    struct inStarFormInit in;
    struct outStarForm out;
    in.dTime = dTime;
    in.bCCSNFeedback = param.bCCSNFeedback;
    in.bSNIaFeedback = param.bSNIaFeedback;
    in.dCCSNFBDelay = param.dCCSNFBDelay;
    in.dSNIaFBDelay = param.dSNIaFBDelay;

    pstStarFormInit(pst, &in, sizeof(in), &out, sizeof(out));

    printf("%d star particles are about to explode in this IC\n", out.nFormed);
#endif
}

#ifdef __cplusplus
extern "C" {
#endif

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
        mdlGetReply(pst->mdl,rID,&fsStats,NULL);
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
    for (int i = 0; i < pkdLocal(pkd); ++i) {
        PARTICLE *p = pkdParticle(pkd,i);
        if (pkdIsStar(pkd,p)) {
            STARFIELDS *pStar = pkdStar(pkd,p);
            if (pStar->fTimer >= 0) { // fTimer < 0 can be used for stars
                // in the IC that are not supossed to explode
                if ((in.dTime - pStar->fTimer) < in.dCCSNFBDelay)
                    pStar->bCCSNFBDone = in.bCCSNFeedback ? 0 : 1;

                if ((in.dTime - pStar->fTimer) < in.dSNIaFBDelay)
                    pStar->bSNIaFBDone = in.bSNIaFeedback ? 0 : 1;

                if (!(pStar->bCCSNFBDone && pStar->bSNIaFBDone))
                    *nFormed += 1;
            }
        }
    }
#endif
}
#ifdef __cplusplus
}
#endif
