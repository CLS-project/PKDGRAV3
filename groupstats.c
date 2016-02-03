#include <math.h>
#include "groupstats.h"
#include "group.h"

/*
** Stuff to calculate group properties...
*/
static void initMinPot(void *vpkd,void *v) {
    }

static void combMinPot(void *vpkd, void *v1, void *v2) {
    TinyGroupTable *g1 = (TinyGroupTable *)v1;
    TinyGroupTable *g2 = (TinyGroupTable *)v2;
    int j;
   
    if (g2->minPot < g1->minPot) {
	g1->minPot = g2->minPot;
	for (j=0;j<3;++j) {
	    g1->rPot[j] = g2->rPot[j];
	    }
	}
    }

static void initTinyGroup(void *vpkd, void *v) {
    TinyGroupTable * g = (TinyGroupTable *)v;
    int j;

    g->fMass = 0.0;
    for (j=0;j<3;j++) {
	g->rcom[j] = 0.0;
	g->vcom[j] = 0.0;
	}
    g->rMax =  0.0;
    }

static void combTinyGroup(void *vpkd, void *v1, void *v2) {
    TinyGroupTable *g1 = (TinyGroupTable *)v1;
    TinyGroupTable *g2 = (TinyGroupTable *)v2;
    int j;

    g1->fMass += g2->fMass;
    for (j=0;j<3;j++) {
	g1->rcom[j] += g2->rcom[j];
	g1->vcom[j] += g2->vcom[j];
	}
    if (g2->rMax > g1->rMax) g1->rMax = g2->rMax;
    }

struct packHopCtx {
    PKD pkd;
    int iIndex;
    };

static int packHop(void *vctx, int *id, size_t nSize, void *vBuff) {
    struct packHopCtx *ctx = (struct packHopCtx *)vctx;
    int nLeft = ctx->pkd->nLocalGroups - ctx->iIndex;
    int n = nSize / sizeof(TinyGroupTable);
    if ( n > nLeft ) n = nLeft;
    memcpy(vBuff,ctx->pkd->tinyGroupTable + 1 + ctx->iIndex, n*sizeof(TinyGroupTable) );
    ctx->iIndex += n;
    return n*sizeof(TinyGroupTable);
    }

/* Send the group information to processor 0 */
void pkdHopSendStats(PKD pkd) {
    struct packHopCtx ctx;
    ctx.pkd = pkd;
    ctx.iIndex = 0;
    mdlSend(pkd->mdl,0,packHop, &ctx);
    }

void pkdCalculateGroupStats(PKD pkd, int bPeriodic, double *dPeriod) {
    MDL mdl = pkd->mdl;
    PARTICLE *p;
    float fPot;
    double fMass;
    double r[3],r2,rMax;
    double dHalf[3];
    int i,j,gid;
    vel_t *v;
    TinyGroupTable *g;
    int nLocalGroups;

    pkd->tinyGroupTable = (TinyGroupTable *)(&pkd->ga[pkd->nGroups]);
    /*
    ** Initialize the table.
    */
    nLocalGroups = 0;
    for (gid=0;gid<pkd->nGroups;++gid) {
	if (pkd->ga[gid].id.iPid == pkd->idSelf && gid) ++nLocalGroups;
	pkd->tinyGroupTable[gid].n = pkd->ga[gid].nTotal;
	pkd->tinyGroupTable[gid].minPot = FLOAT_MAXVAL;
	pkd->tinyGroupTable[gid].rMax = 0;
	pkd->tinyGroupTable[gid].fMass = 0;
	for (j=0;j<3;++j) {
	    pkd->tinyGroupTable[gid].rcom[j] = 0;
	    pkd->tinyGroupTable[gid].vcom[j] = 0;
	    }
	}
    /*
    ** First determine the minimum potential particle for each group.
    ** This will be the reference position for the group as well.
    */
    for (i=0;i<pkd->nLocal;++i) {
	p = pkdParticle(pkd,i);
	gid = pkdGetGroup(pkd,p);
	if (!gid) continue;
	fPot = *pkdPot(pkd,p);
	if (fPot < pkd->tinyGroupTable[gid].minPot) {
	    pkd->tinyGroupTable[gid].minPot = fPot;
	    for (j=0;j<3;++j) {
		pkd->tinyGroupTable[gid].rPot[j] = pkdPosRaw(pkd,p,j);
		}
	    }
	}
    /*
    ** Now look at remote group particle to see if we have a lower potential.
    ** Note that we only expose the local groups to the cache! This allows us 
    ** to make any desired update to the remote group entries of the table.
    */
    mdlCOcache(mdl,CID_GROUP,NULL,pkd->tinyGroupTable,sizeof(TinyGroupTable),nLocalGroups+1,
	NULL,initMinPot,combMinPot);
    for(gid=1+nLocalGroups;gid<pkd->nGroups;++gid) {
	assert(pkd->ga[gid].id.iPid != pkd->idSelf);
	g = mdlFetch(mdl,CID_GROUP,pkd->ga[gid].id.iIndex,pkd->ga[gid].id.iPid);
	if (g->minPot <= pkd->tinyGroupTable[gid].minPot) {
	    pkd->tinyGroupTable[gid].minPot = g->minPot;
	    for (j=0;j<3;++j) {
		pkd->tinyGroupTable[gid].rPot[j] = g->rPot[j];
		}
	    }
	else {
	    /*
	    ** We update the remote center (the inverse of the above).
	    */
	    g->minPot = pkd->tinyGroupTable[gid].minPot;
	    for (j=0;j<3;++j) {
		g->rPot[j] = pkd->tinyGroupTable[gid].rPot[j];
		}	   
	    }
	}
    mdlFinishCache(mdl,CID_GROUP);
    /*
    ** Now based on the newly defined reference point for the groups we can
    ** determine the other group properties.
    */
    for (j=0;j<3;++j) dHalf[j] = bPeriodic ? 0.5 * dPeriod[j] : FLOAT_MAXVAL;
    for (i=0;i<pkd->nLocal;++i) {
	p = pkdParticle(pkd,i);
	gid = pkdGetGroup(pkd,p);
	fMass = pkdMass(pkd,p);
	v = pkdVel(pkd,p);
	pkd->tinyGroupTable[gid].fMass += fMass;
	if (gid > 0) {
	    r2 = 0.0;
	    for (j=0;j<3;++j) {
		r[j] =  pkdPosToDbl(pkd,pkdPosRaw(pkd,p,j) - pkd->tinyGroupTable[gid].rPot[j]);
		if      (r[j] < -dHalf[j]) r[j] += dPeriod[j];
		else if (r[j] > +dHalf[j]) r[j] -= dPeriod[j];
		r2 += r[j]*r[j];
		pkd->tinyGroupTable[gid].rcom[j] += fMass*r[j];
		pkd->tinyGroupTable[gid].vcom[j] += fMass*v[j];
		}
	    rMax = sqrt(r2);
	    if (rMax > pkd->tinyGroupTable[gid].rMax) pkd->tinyGroupTable[gid].rMax = rMax;
	    }
	}
    /* 
    ** Now accumulate totals globally
    */
    mdlCOcache(mdl,CID_GROUP,NULL,pkd->tinyGroupTable,sizeof(TinyGroupTable),nLocalGroups+1,
	NULL,initTinyGroup,combTinyGroup);
    for(gid=1+nLocalGroups;gid<pkd->nGroups;++gid) {
	TinyGroupTable *g;
	g = mdlVirtualFetch(mdl,CID_GROUP,pkd->ga[gid].id.iIndex,pkd->ga[gid].id.iPid);
	g->fMass += pkd->tinyGroupTable[gid].fMass;
	for (j=0;j<3;j++) {
	    g->rcom[j] += pkd->tinyGroupTable[gid].rcom[j];
	    g->vcom[j] += pkd->tinyGroupTable[gid].vcom[j];
	    }
	if (pkd->tinyGroupTable[gid].rMax > g->rMax) g->rMax = pkd->tinyGroupTable[gid].rMax;
	}
    mdlFinishCache(mdl,CID_GROUP);
    for(gid=1;gid<=nLocalGroups;++gid) {
	for (j=0;j<3;j++) {
	    pkd->tinyGroupTable[gid].rcom[j] /= pkd->tinyGroupTable[gid].fMass;
	    pkd->tinyGroupTable[gid].rcom[j] += pkdPosToDbl(pkd,pkd->tinyGroupTable[gid].rPot[j]);
	    pkd->tinyGroupTable[gid].vcom[j] /= pkd->tinyGroupTable[gid].fMass;
	    }
	}
    /*
    ** Fetch remote group properties (we probably only need rMax to be fetched here).
    */
    mdlROcache(mdl,CID_GROUP,NULL,pkd->tinyGroupTable,sizeof(TinyGroupTable),nLocalGroups+1);
    for(gid=1+nLocalGroups;gid<pkd->nGroups;++gid) {
	g = mdlFetch(mdl,CID_GROUP,pkd->ga[gid].id.iIndex,pkd->ga[gid].id.iPid);
	pkd->tinyGroupTable[gid].fMass = g->fMass;
	for (j=0;j<3;++j) {
	    pkd->tinyGroupTable[gid].rcom[j] = g->rcom[j];
	    pkd->tinyGroupTable[gid].vcom[j] = g->vcom[j];
	    }
	pkd->tinyGroupTable[gid].rMax = g->rMax;
	}
    mdlFinishCache(mdl,CID_GROUP);
    /*
    ** Important to really make sure pkd->nLocalGroups is set correctly before output!
    */
    pkd->nLocalGroups = nLocalGroups;
    }
