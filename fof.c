#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <math.h>
#include "hop.h"

/*
** We have a super fancy mindist in smooth.c which we should be using...
*/
#define MINDIST(bnd,pos,min2) {\
    double BND_dMin;\
    int BND_j;\
    (min2) = 0;					\
    for (BND_j=0;BND_j<3;++BND_j) {\
	BND_dMin = fabs((bnd)->fCenter[BND_j] - (pos)[BND_j]) - (bnd)->fMax[BND_j]; \
	if (BND_dMin > 0) (min2) += BND_dMin*BND_dMin;			\
	}\
    }

void smFofGatherLocal(SMX smx,FLOAT fBall2,FLOAT r[3],int32_t iGroup) {
    PKD pkd = smx->pkd;
    KDN *kdn;
    PARTICLE *p;
    double p_r[3];
    FLOAT min2,dx,dy,dz,fDist2;
    int *S = smx->S;
    int sp = 0;
    int iCell,pj,pEnd;
    const BND *bnd;
    int32_t *pPartGroup;

    kdn = pkdTreeNode(pkd,iCell = ROOT);
    while (1) {
        bnd = pkdNodeBnd(pkd, kdn);
	MINDIST(bnd,r,min2);
	if (min2 > fBall2) {
	    goto NoIntersect;
	}
	/*
	** We have an intersection to test.
	*/
	if (kdn->iLower) {
	    kdn = pkdTreeNode(pkd,iCell = kdn->iLower);
	    S[sp++] = iCell+1;
	    continue;
	}
	else {
	    pEnd = kdn->pUpper;
	    for (pj=kdn->pLower;pj<=pEnd;++pj) {
		p = pkdParticle(pkd,pj);
		pPartGroup = pkdInt32(p,pkd->oGroup);
		if (*pPartGroup) continue;		    
		pkdGetPos1(pkd,p,p_r);
		dx = r[0] - p_r[0];
		dy = r[1] - p_r[1];
		dz = r[2] - p_r[2];
		fDist2 = dx*dx + dy*dy + dz*dz;
		if (fDist2 <= fBall2) {
		    /*
		    **  Mark particle and add it to the do-fifo
		    */
		    *pPartGroup = iGroup;
		    smx->Fifo[smx->iTail++] = pj;
		}
	    }
	}
    NoIntersect:
	if (sp) kdn = pkdTreeNode(pkd,iCell = S[--sp]);
	else break;
    }
}


int smNewFof(SMX smx,SMF *smf) {
    PKD pkd = smx->pkd;
    MDL mdl = pkd->mdl;
    PARTICLE *p;
    double R[3];
    double dTau2 = smf->dTau2;
    struct smGroupArray *ga = smx->ga;
    int32_t iGroup,*pGroup;
    int pn,pi,i;

    assert(pkd->oGroup); /* Validate memory model */
    smx->Fifo = (int32_t *)(pkd->pLite); /* Used only for SPH */
    pkd->nGroups = 0;    
    iGroup = 0;
    for (pn=0;pi<pkd->nLocal;++pn) {
	p = pkdParticle(pkd,pn);
	pGroup = pkdInt32(p,pkd->oGroup);
	if (*pGroup ) continue;
	iGroup++;
	/*
	** Mark particle and add it to the do-fifo
	*/
	*pGroup = iGroup;
	smx->iHead = smx->iTail = 0;
	smx->Fifo[smx->iTail++] = pn;
	while (smx->iHead != smx->iTail) {
	    pi = smx->Fifo[smx->iHead++];
	    p = pkdParticle(pkd,pi);
	    pkdGetPos1(pkd,p,R);
	    smFofGatherLocal(smx,dTau2,R,iGroup);
	    }
	assert(smx->iTail < pkd->nLocal);
	}
    smx->Fifo = NULL;  /* done with the Fifo, can use the storage for other stuff */
    pkd->nLocalGroups = iGroup;
    pkd->nGroups = iGroup + 1;
    /*
    ** Now we need to find links to remote groups and make global groups.
    */



    /*
    ** Create final group table.
    */

    pkd->hopGroups = mdlMalloc(mdl, pkd->nGroups * sizeof(HopGroupTable));
    assert(pkd->hopGroups!=NULL);
    for(i=0; i<pkd->nGroups; ++i) {
	pkd->hopGroups[i].id.iPid      = pkd->idSelf;
	pkd->hopGroups[i].id.iIndex    = i;
	pkd->hopGroups[i].bNeedGrav    = 1;
	pkd->hopGroups[i].bComplete    = 0;
	}
//    purgeSmallGroups(pkd,nMinGroupSize,bPeriodic,dPeriod);
    pkd->hopSavedRoots = 0;

    return pkd->nLocalGroups;
    }
