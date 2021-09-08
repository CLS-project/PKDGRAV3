#include "starformation/starformation.h"

void msrStarFormInit(MSR msr, double dTime){
   struct inStarForm in;
   struct outStarForm out;
   in.dTime = dTime;

   pstStarFormInit(msr->pst, &in, sizeof(in), &out, sizeof(out));

   printf("%d star particles are about to explode in this IC\n", out.nFormed);
}


int pstStarFormInit(PST pst,void *vin,int nIn,void *vout,int nOut) {
    struct inStarForm *in = vin;
    struct outStarForm *out = vout;
    int rID;

    mdlassert(pst->mdl,nIn == sizeof(struct inStarForm));
    if (pst->nLeaves > 1) {
       struct outStarForm fsStats;

       rID = mdlReqService(pst->mdl,pst->idUpper,PST_STARFORMINIT,in,nIn);
       pstStarFormInit(pst->pstLower,in,nIn,vout,nOut);
       mdlGetReply(pst->mdl,rID,&fsStats,NULL);
       out->nFormed += fsStats.nFormed;
       }
    else {
       pkdStarFormInit(pst->plcl->pkd, in->dTime, &out->nFormed);
       }
    return sizeof(struct outStarForm);
    }

void pkdStarFormInit(PKD pkd, double dTime, int *nFormed){
    *nFormed = 0;
    for (int i=0;i<pkdLocal(pkd);++i) {
       PARTICLE *p = pkdParticle(pkd,i);
      if (pkdIsStar(pkd,p)){
          STARFIELDS* pStar = pkdStar(pkd,p);
          if (pStar->fTimer >= 0){// fTimer < 0 can be used for stars
                                  // in the IC that are not supossed to explode
#ifdef FEEDBACK
            if ( (dTime-pStar->fTimer) < pkd->param.dFeedbackDelay)
#endif
            {
               pStar->hasExploded = 0; // This particle has not exploded before
               *nFormed += 1;
            }
          }
      }
    }
}
