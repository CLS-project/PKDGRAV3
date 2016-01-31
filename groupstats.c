#include <math.h>
#include "groupstats.h"

/*
** Stuff to calculate group properties...
*/

static void initHopGroupProperties(void *vpkd, void *v) {
    HopGroupTable * g = (HopGroupTable *)v;
    int j;

    g->nLocal = 0;
    g->nTotal = 0;
    g->nRemote = 0;
    for (j=0;j<3;j++) {
	g->ravg[j] = 0.0;
	g->rcom[j] = 0.0;
	g->vcom[j] = 0.0;
	}
    g->fRMSRadius =  0.0;
    g->fMass = 0.0;
    }

static void combHopGroupProperties(void *vpkd, void *v1, void *v2) {
    HopGroupTable * g1 = (HopGroupTable *)v1;
    HopGroupTable * g2 = (HopGroupTable *)v2;
    int j;

    g1->nTotal += g2->nTotal;
    g1->nRemote += g2->nRemote;
    for (j=0;j<3;j++) {
	if (g2->rmin[j] < g1->rmin[j]) g1->rmin[j] = g2->rmin[j];
	if (g2->rmax[j] > g1->rmax[j]) g1->rmax[j] = g2->rmax[j];
	g1->ravg[j] += g2->ravg[j];
	g1->rcom[j] += g2->rcom[j];
	g1->vcom[j] += g2->vcom[j];
	}
    g1->fRMSRadius += g2->fRMSRadius;
    g1->fMass += g2->fMass;
    g1->bComplete = g1->bComplete && g2->bComplete;
    }

struct packHopCtx {
    PKD pkd;
    int iIndex;
    };

static int packHop(void *vctx, int *id, size_t nSize, void *vBuff) {
    struct packHopCtx *ctx = (struct packHopCtx *)vctx;
    int nLeft = ctx->pkd->nLocalGroups - ctx->iIndex;
    int n = nSize / sizeof(HopGroupTable);
    if ( n > nLeft ) n = nLeft;
    memcpy(vBuff,ctx->pkd->hopGroups + 1 + ctx->iIndex, n*sizeof(HopGroupTable) );
    ctx->iIndex += n;
    return n*sizeof(HopGroupTable);
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
    HopGroupTable * g;
    int i,j,gid;
    PARTICLE *p;
    double dHalf[3];
    float fMass;
    double r[3];
    vel_t *v;

    pkd->hopGroups = (HopGroupTable *)&pkd->ga[pkd->nGroups];
    /*
    ** Copy the name of the group to the group table structure for now, 
    ** but actually we don't need to duplicate this information.
    */
    for (i=0;i<pkd->nGroups;++i) {
	pkd->hopGroups[i].id.iPid = pkd->ga[i].id.iPid;
	pkd->hopGroups[i].id.iIndex = pkd->ga[i].id.iIndex;
	pkd->hopGroups[i].bNeedGrav = 1;
	pkd->hopGroups[i].bComplete = 0;	
	}

#ifdef TEST_SINGLE_GRAVITY
    for(gid=2; gid<=pkd->nLocalGroups; ++gid)
	pkd->hopGroups[gid].bNeedGrav = 0;
#endif

    for(i=0; i<pkd->nGroups; ++i) {
	pkd->hopGroups[i].iGlobalId = 0;
	pkd->hopGroups[i].nLocal = 0;
	pkd->hopGroups[i].nTotal = 0;
	pkd->hopGroups[i].nRemote = 0;
	pkd->hopGroups[i].fMass = 0.0;
	for (j=0;j<3;j++) {
	    pkd->hopGroups[i].ravg[j] = 0.0;
	    pkd->hopGroups[i].rcom[j] = 0.0;
	    pkd->hopGroups[i].vcom[j] = 0.0;
	    }
	pkd->hopGroups[i].fRMSRadius = 0.0;
	}

    /* Get a reference point for each group; any particle will do */
    for (i=0;i<pkd->nLocal;++i) {
	p = pkdParticle(pkd,i);
	gid = *pkdGroup(pkd,p);
	pkdGetPos1(pkd,p,pkd->hopGroups[gid].rref);
	}
    for(gid=1; gid<=pkd->nLocalGroups; ++gid) {
	for (j=0;j<3;j++) {
	    pkd->hopGroups[gid].rmin[j] = pkd->hopGroups[gid].rmax[j] = pkd->hopGroups[gid].rref[j];
	    }
	}
    mdlROcache(mdl,CID_GROUP,NULL,pkd->hopGroups,sizeof(HopGroupTable), pkd->nGroups);
    for(i=1+pkd->nLocalGroups; i<pkd->nGroups; ++i) {
	assert(pkd->hopGroups[i].id.iPid != pkd->idSelf);
	g = mdlAcquire(mdl,CID_GROUP,pkd->hopGroups[i].id.iIndex,pkd->hopGroups[i].id.iPid);
	for(j=0;j<3;++j) pkd->hopGroups[i].rref[j] = g->rref[j];
	mdlRelease(mdl,CID_GROUP,g);
	}
    mdlFinishCache(mdl,CID_GROUP);

    for(j=0; j<3; ++j) dHalf[j] = bPeriodic ? 0.5 * dPeriod[j] : FLOAT_MAXVAL;
    for (i=0;i<pkd->nLocal;++i) {
	p = pkdParticle(pkd,i);
	gid = *pkdGroup(pkd,p);
	pkdGetPos1(pkd,p,r);
	v = pkdVel(pkd,p);
	fMass = pkdMass(pkd,p);

	/*
	** Groups cannot be larger than half a box width -- this is a requirement.
	** We correct the particle position by wrapping if this constraint is
	** violated. This cannot happen with non-periodic boxes of course.
	*/
	for (j=0;j<3;j++) {
	    double Gr = pkd->hopGroups[gid].rref[j];
	    if      (r[j] < Gr - dHalf[j]) r[j] += dPeriod[j];
	    else if (r[j] > Gr + dHalf[j]) r[j] -= dPeriod[j];
	    assert(r[j] > Gr - dHalf[j] && r[j] < Gr + dHalf[j] );
	    }

	for (j=0;j<3;j++) {
	    if (r[j]<pkd->hopGroups[gid].rmin[j]) pkd->hopGroups[gid].rmin[j] = r[j];
	    else if (r[j]>pkd->hopGroups[gid].rmax[j]) pkd->hopGroups[gid].rmax[j] = r[j];
	    pkd->hopGroups[gid].ravg[j] += r[j];
	    pkd->hopGroups[gid].rcom[j] += r[j]*fMass;
	    pkd->hopGroups[gid].fRMSRadius +=  r[j]*r[j]*fMass;
	    pkd->hopGroups[gid].vcom[j] += v[j]*fMass;
	    }
	pkd->hopGroups[gid].fMass += fMass;
	++pkd->hopGroups[gid].nLocal;
	++pkd->hopGroups[gid].nTotal;
	}

    /* Now accumulate totals globally */
    mdlCOcache(mdl,CID_GROUP,NULL,pkd->hopGroups,sizeof(HopGroupTable), pkd->nGroups,
	pkd, initHopGroupProperties, combHopGroupProperties );

    for(i=1+pkd->nLocalGroups; i<pkd->nGroups; ++i) {
	HopGroupTable *g;
	if (pkd->hopGroups[i].id.iPid == pkd->idSelf) continue;
	g = mdlVirtualFetch(mdl,CID_GROUP,pkd->hopGroups[i].id.iIndex,pkd->hopGroups[i].id.iPid);
	g->nTotal += pkd->hopGroups[i].nLocal;
	g->nRemote = 1;
	for (j=0;j<3;j++) {
	    g->rmin[j] = pkd->hopGroups[i].rmin[j];
	    g->rmax[j] = pkd->hopGroups[i].rmax[j];
	    g->ravg[j] += pkd->hopGroups[i].ravg[j];
	    g->rcom[j] += pkd->hopGroups[i].rcom[j];
	    g->vcom[j] += pkd->hopGroups[i].vcom[j];
	    }
	g->fRMSRadius += pkd->hopGroups[i].fRMSRadius;
	g->fMass += pkd->hopGroups[i].fMass;
	g->bNeedGrav = pkd->hopGroups[i].bNeedGrav;
	g->bComplete = pkd->hopGroups[i].bComplete;
	}
    mdlFinishCache(mdl,CID_GROUP);

    int nActive = 0;
    for(i=1; i<=pkd->nLocalGroups; ++i) {
	for (j=0;j<3;j++) {
	    pkd->hopGroups[i].ravg[j] /= pkd->hopGroups[i].nTotal;
	    pkd->hopGroups[i].rcom[j] /= pkd->hopGroups[i].fMass;
	    pkd->hopGroups[i].vcom[j] /= pkd->hopGroups[i].fMass;
	    }
	/* 
	** Do not calculate fDeltaR2 with the corrected positions!
	*/
	pkd->hopGroups[i].fRMSRadius /= pkd->hopGroups[i].fMass;
	pkd->hopGroups[i].fRMSRadius -= pkd->hopGroups[i].rcom[0]*pkd->hopGroups[i].rcom[0]
	    + pkd->hopGroups[i].rcom[1]*pkd->hopGroups[i].rcom[1]
	    + pkd->hopGroups[i].rcom[2]*pkd->hopGroups[i].rcom[2];
	pkd->hopGroups[i].fRMSRadius = sqrt(pkd->hopGroups[i].fRMSRadius);
	if (pkd->hopGroups[i].bComplete) pkd->hopGroups[i].bNeedGrav = 0;
	if (!pkd->hopGroups[i].bComplete) nActive++;
	}

   /* Fetch remote group properties */
    mdlROcache(mdl,CID_GROUP,NULL,pkd->hopGroups,sizeof(HopGroupTable), pkd->nGroups);
    for(i=1+pkd->nLocalGroups; i<pkd->nGroups; ++i) {
	assert(pkd->hopGroups[i].id.iPid != pkd->idSelf);
	g = mdlAcquire(mdl,CID_GROUP,pkd->hopGroups[i].id.iIndex,pkd->hopGroups[i].id.iPid);
	assert(pkd->hopGroups[i].id.iIndex==g->id.iIndex && pkd->hopGroups[i].id.iPid==g->id.iPid);
	pkd->hopGroups[i].nRemote = g->nRemote;
	pkd->hopGroups[i].nTotal = g->nTotal;
	pkd->hopGroups[i].fMass = g->fMass;
	pkd->hopGroups[i].fRMSRadius = g->fRMSRadius;
	pkd->hopGroups[i].bNeedGrav = g->bNeedGrav;
	pkd->hopGroups[i].bComplete = g->bComplete;
	for(j=0;j<3;++j) {
	    pkd->hopGroups[i].rmin[j] = g->rmin[j];
	    pkd->hopGroups[i].rmax[j] = g->rmax[j];
	    pkd->hopGroups[i].ravg[j] = g->ravg[j];
	    pkd->hopGroups[i].rcom[j] = g->rcom[j];
	    pkd->hopGroups[i].vcom[j] = g->vcom[j];
	    }
	mdlRelease(mdl,CID_GROUP,g);
	}
    mdlFinishCache(mdl,CID_GROUP);
    }


