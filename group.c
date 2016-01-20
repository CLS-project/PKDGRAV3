#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <math.h>
#include "group.h"

/*
** We are about to renumber the groups locally. Each group has the new ID in iNewGid.
*/
static void updateGroupIds(PKD pkd, int nGroups, struct smGroupArray *ga, int bIndexIsGID) {
    MDL mdl = pkd->mdl;
    PARTICLE *p;
    int pi, gid;
    int *g, *newGid = NULL;

    /* Update the group for all local particles */
    for (pi=0;pi<pkd->nLocal;++pi) {
	p = pkdParticle(pkd,pi);
	gid = *pkdGroup(pkd,p);
	if (gid<=0) continue;
	*pkdGroup(pkd,p) = ga[gid].iNewGid;
	}
    /* Now gid has the new position -- a reorder is necessary */
    if (bIndexIsGID) {
	newGid = mdlMalloc(mdl, nGroups * sizeof(int));
	assert(newGid!=NULL);
	for(gid=1; gid<nGroups; ++gid) {
	    newGid[gid] = ga[gid].iNewGid;
	    }
	mdlROcache(mdl,CID_GROUP,NULL,newGid,sizeof(int), nGroups);
	for(gid=1; gid<nGroups; ++gid) {
	    if (ga[gid].id.iPid==pkd->idSelf) {
		ga[gid].id.iIndex = newGid[ga[gid].id.iIndex];
		}
	    else {
		g = mdlFetch(pkd->mdl,CID_GROUP,ga[gid].id.iIndex,ga[gid].id.iPid);
		ga[gid].id.iIndex = *g;
		}
	    }
	mdlFinishCache(mdl,CID_GROUP);
	mdlFree(mdl,newGid);
	}
    }

static int smga_cmp(const void *a0, const void *b0) {
    struct smGroupArray *a = (struct smGroupArray *)a0;
    struct smGroupArray *b = (struct smGroupArray *)b0;
    if (a->id.iPid < b->id.iPid)     return -1;
    if (a->id.iPid > b->id.iPid)     return +1;
    if (a->id.iIndex < b->id.iIndex) return -1;
    if (a->id.iIndex > b->id.iIndex) return +1;
    return 0;
    }

/* Put the groups back in "group" order  */
static int reorderGroups(PKD pkd, int nGroups, struct smGroupArray *ga) {
    struct smGroupArray tmp;
    int pi,nNew;

    /* Move the dead groups (if any) to the end by partitioning */
    nNew = 1;
    pi = nGroups-1;
    PARTITION(nNew<pi,nNew<=pi,++nNew,--pi,
    {tmp = ga[pi]; ga[pi]=ga[nNew]; ga[nNew]=tmp;},
	ga[nNew].iGid>0,ga[pi].iGid<=0);
    for(pi=1; pi<nGroups; ++pi) {
	if (pi<nNew) {
	    assert(ga[pi].iGid>0);
	    assert(ga[pi].iGid<nNew);
	    }
	else assert(ga[pi].iGid==0);
	}
    /* Now just do a simple reorder */
    for(pi=1; pi<nNew; ++pi) {
	while(ga[pi].iGid != pi) {
	    /* If the swap destination is already correct, then we are a duplicate - drop it */
	    if (ga[ga[pi].iGid].iGid == ga[pi].iGid) {
		if (--nNew == pi) break;
		tmp = ga[nNew];
		}
	    else {
		tmp = ga[ga[pi].iGid];
		ga[ga[pi].iGid] = ga[pi];
		}
	    ga[pi] = tmp;
	    }
	}
    for(pi=1; pi<nNew; ++pi) {
	assert(ga[pi].iGid == pi);
	assert(ga[pi].iNewGid >= 0);
	}
    return nNew;
    }

/*
** Sort groups by Processor and ID, with local groups first
** The new group ID is in iNewGid and groups are in still in the old order.
*/
static int renumberGroups(PKD pkd, int nGroups,struct smGroupArray *ga) {
    int i, gid, nNew;

    /* We want local groups first */
    for(i=1; i<nGroups; ++i) {
	assert(ga[i].iGid == i);
	assert(ga[i].id.iPid>=0);
	if (ga[i].id.iPid == pkd->idSelf) ga[i].id.iPid = -1;
	}
    qsort(ga+1,nGroups-1,sizeof(struct smGroupArray),smga_cmp);
    for(i=1; i<nGroups; ++i) if (ga[i].id.iPid == -1) ga[i].id.iPid = pkd->idSelf;
    gid = 0; /* Sentinal node always as iPid=idSelf, iIndex == -1 */
    ga[0].id.iPid = pkd->idSelf;
    ga[0].id.iIndex = -1;
    /* Count and note *unique* groups */
    for(i=nNew=1; i<nGroups; ++i) {
	if (ga[i].id.iPid!=ga[gid].id.iPid || ga[i].id.iIndex!=ga[gid].id.iIndex) {
	    ++nNew;
	    gid = i;
	    }
	ga[i].iNewGid = nNew-1;
	}
    reorderGroups(pkd,nGroups,ga);
    return nNew;
    }

/*
** Merge duplicate groups (iPid/iIndex is the same), and update particle pointers.
*/
int pkdCombineDuplicateGroupIds(PKD pkd, int nGroups, struct smGroupArray *ga,int bIndexIsGID) {
    int gid, nNew;

    nNew = renumberGroups(pkd,nGroups,ga);
    updateGroupIds(pkd,nGroups,ga,bIndexIsGID);
    /* We will now put the groups in the new order */
    for(gid=1; gid<nGroups; ++gid) ga[gid].iGid = ga[gid].iNewGid;
    nGroups = reorderGroups(pkd,nGroups,ga);
    return nGroups;
    }


static void initMaxHopGroup(void *vpkd, void *v) {}
static void combMaxHopGroup(void *vctx, void *v1, void *v2) {
    HopGroupTable * g1 = (HopGroupTable *)v1;
    HopGroupTable * g2 = (HopGroupTable *)v2;
    /* Remember: nTotal is really MAX(nLocal) here. */
    if (g2->nTotal > g1->nTotal) {
	g1->nTotal = g2->nTotal;
	g1->id.iPid = g2->id.iPid;
	g1->id.iIndex = g2->id.iIndex;
	}
    }

/*
** For each group, find the processor that has the most particles
** and make that processor that master for that group. If any
** group has zero particles, then remove it.
*/
static void hopRelocateGroups(PKD pkd) {
    MDL mdl = pkd->mdl;
    struct smGroupArray *ga = (struct smGroupArray *)(pkd->pLite);
    HopGroupTable *g;
    PARTICLE *p;
    int i, gid, n;

    /* Count local members of all groups */
    for(i=0; i<pkd->nGroups; ++i) pkd->hopGroups[i].nLocal = 0;
    for (i=0;i<pkd->nLocal;++i) {
	p = pkdParticle(pkd,i);
	gid = *pkdGroup(pkd,p);
	++pkd->hopGroups[gid].nLocal;
	}
    for(i=0; i<pkd->nGroups; ++i) /* Reuse nTotal */
	pkd->hopGroups[i].nTotal = pkd->hopGroups[i].nLocal;

    /* Now find the real processor with the most particles for each group */
    mdlCOcache(mdl,CID_GROUP,NULL,pkd->hopGroups,sizeof(HopGroupTable), pkd->nGroups,
	pkd, initMaxHopGroup, combMaxHopGroup );
    for(i=1+pkd->nLocalGroups; i<pkd->nGroups; ++i) {
	assert(pkd->hopGroups[i].id.iPid != pkd->idSelf);
	g = mdlAcquire(mdl,CID_GROUP,pkd->hopGroups[i].id.iIndex,pkd->hopGroups[i].id.iPid);
	assert(g->id.iPid!=pkd->idSelf);
	g->id.iPid = pkd->idSelf;
	g->id.iIndex = i;
	g->nTotal = pkd->hopGroups[i].nLocal;
	mdlRelease(mdl,CID_GROUP,g);
	}
    mdlFinishCache(mdl,CID_GROUP);

    /* Now update the new group location */
    mdlROcache(mdl,CID_GROUP,NULL,pkd->hopGroups,sizeof(HopGroupTable), pkd->nGroups);
    for(i=1+pkd->nLocalGroups; i<pkd->nGroups; ++i) {
	g = mdlAcquire(mdl,CID_GROUP,pkd->hopGroups[i].id.iIndex,pkd->hopGroups[i].id.iPid);
	pkd->hopGroups[i].id.iPid = g->id.iPid;
	pkd->hopGroups[i].id.iIndex = g->id.iIndex;
	mdlRelease(mdl,CID_GROUP,g);
	}
    mdlFinishCache(mdl,CID_GROUP);

    for(i=1; i<pkd->nGroups; ++i) {
	ga[i].iGid = i;
	if (pkd->hopGroups[i].nLocal) {
	    ga[i].id.iPid = pkd->hopGroups[i].id.iPid;
	    ga[i].id.iIndex = pkd->hopGroups[i].id.iIndex;
	    assert(ga[i].id.iPid>=0);
	    }
	else { /* Remove this group */
	    ga[i].id.iPid = pkd->idSelf;
	    ga[i].id.iIndex = -1;
	    }
	}

    renumberGroups(pkd,pkd->nGroups,ga);
    updateGroupIds(pkd,pkd->nGroups,ga,1);

    n = pkd->nGroups;
    for(gid=1; gid<n; ) {
	struct smGroupArray tmp1;
	HopGroupTable tmp2;
	if (ga[gid].iNewGid == 0) {
	    if (--n == gid) break;
	    ga[gid] = ga[n];
	    pkd->hopGroups[gid] = pkd->hopGroups[n];
	    }
	else if (ga[gid].iNewGid != gid) {
	    assert(ga[gid].iNewGid > gid);
	    tmp1 = ga[gid];
	    tmp2 = pkd->hopGroups[gid];
	    ga[gid] = ga[tmp1.iNewGid];
	    pkd->hopGroups[gid] = pkd->hopGroups[tmp1.iNewGid];
	    ga[tmp1.iNewGid] = tmp1;
	    pkd->hopGroups[tmp1.iNewGid] = tmp2;
	    }
	else ga[gid].iGid = gid++; /*FIXME: correct? */
	}
    pkd->nGroups = n;

    pkd->nLocalGroups = 0;
    for(i=1; i<pkd->nGroups; ++i) {
	assert(ga[i].iGid == i);
	pkd->hopGroups[i].id.iPid = ga[i].id.iPid;
	pkd->hopGroups[i].id.iIndex = ga[i].id.iIndex;
	if (pkd->hopGroups[i].id.iPid==pkd->idSelf) ++pkd->nLocalGroups;
	}
    }

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

static void hopCalculateGroupStats(PKD pkd, int bPeriodic, double *dPeriod) {
    MDL mdl = pkd->mdl;
    HopGroupTable * g;
    int i,j,gid;
    PARTICLE *p;
    double dHalf[3];
    float fMass;
    double r[3];
    vel_t *v;

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
	g = mdlAcquire(mdl,CID_GROUP,pkd->hopGroups[i].id.iIndex,pkd->hopGroups[i].id.iPid);
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
	mdlRelease(mdl,CID_GROUP,g);
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


void pkdPurgeSmallGroups(PKD pkd,int nMinGroupSize, int bPeriodic, double *dPeriod) {
    struct smGroupArray *ga = (struct smGroupArray *)(pkd->pLite);
    int i, j, gid;

    hopRelocateGroups(pkd);
    hopCalculateGroupStats(pkd,bPeriodic,dPeriod);
    /* Purge groups with too few particles */
    for(i=j=1; i<pkd->nGroups; ++i) {
	assert(pkd->hopGroups[i].nTotal==0 || pkd->hopGroups[i].id.iPid!=pkd->idSelf || pkd->hopGroups[i].nLocal>0);
	if (pkd->hopGroups[i].nLocal==0 || pkd->hopGroups[i].nTotal<nMinGroupSize) gid=0;
	else gid = j++;
	ga[i].id.iPid = pkd->hopGroups[i].id.iPid;
	ga[i].id.iIndex = pkd->hopGroups[i].id.iIndex;
	ga[i].iNewGid = gid;
	}
    updateGroupIds(pkd,pkd->nGroups,ga,1);
    for(i=gid=j=1; i<pkd->nGroups; ++i) {
	if (ga[i].iNewGid && ga[i].id.iIndex) {
	    if (ga[i].id.iPid==pkd->idSelf) ++j;
	    if (i != gid) pkd->hopGroups[gid] = pkd->hopGroups[i];
	    pkd->hopGroups[gid].id.iIndex = ga[i].id.iIndex;
	    ++gid;
	    }
	}
    pkd->nGroups = gid;
    pkd->nLocalGroups = j-1;
    }


int pkdHopCountGID(PKD pkd) {
    MDL mdl = pkd->mdl;
    int nLocal, i;
    for(i=1,nLocal=0; i<pkd->nGroups; ++i) {
	if (pkd->hopGroups[i].id.iPid==mdlSelf(mdl)) ++nLocal;
	}
    return nLocal;
    }


void pkdHopAssignGID(PKD pkd,uint64_t iStartGID) {
    MDL mdl = pkd->mdl;
    int i, nLocal;
    PARTICLE *p;
    HopGroupTable *g;

    if (pkd->hopRootIndex) { free(pkd->hopRootIndex); pkd->hopRootIndex = NULL; }
    if (pkd->hopRoots)     { free(pkd->hopRoots);     pkd->hopRoots = NULL;     }

    for(i=1,nLocal=0; i<pkd->nGroups; ++i) {
	if (pkd->hopGroups[i].id.iPid==mdlSelf(mdl)) ++nLocal;
	}
    for(i=1; i<=nLocal; ++i) {
	assert(pkd->hopGroups[i].id.iPid==mdlSelf(mdl));
	pkd->hopGroups[i].iGlobalId = iStartGID + i;
	}
    mdlROcache(mdl,CID_GROUP,NULL,pkd->hopGroups,sizeof(HopGroupTable), pkd->nGroups);
    for(; i<pkd->nGroups ; ++i) {
	assert(pkd->hopGroups[i].id.iPid!=mdlSelf(mdl));
	g = mdlAcquire(mdl,CID_GROUP,pkd->hopGroups[i].id.iIndex,pkd->hopGroups[i].id.iPid);
	assert(g->id.iPid==pkd->hopGroups[i].id.iPid);
	pkd->hopGroups[i].iGlobalId = g->iGlobalId;
	mdlRelease(mdl,CID_GROUP,g);
	}
    mdlFinishCache(mdl,CID_GROUP);

    for (i=0;i<pkd->nLocal;++i) {
	p = pkdParticle(pkd,i);
	assert( *pkdGroup(pkd,p) < pkd->nGroups );
	*pkdGroup(pkd,p) = pkd->hopGroups[*pkdGroup(pkd,p)].iGlobalId;
	}
    }
