#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "hop.h"
#include <math.h>

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

static void initJoinLoops(void *vctx, void *v) {}
static void combJoinLoops(void *vctx, void *v1, void *v2) {
    SMF *smf = (SMF *)vctx;
    GHtmpGroupTable * g1 = (GHtmpGroupTable *)v1;
    GHtmpGroupTable * g2 = (GHtmpGroupTable *)v2;
    if ( g1->iPid>g2->iPid || (g1->iPid==g2->iPid && g1->iIndex>g2->iIndex) ) {
	g1->iPid = g2->iPid;
	g1->iIndex = g2->iIndex;
	smf->bDone = 0;
	}
    }

static void initSetArc(void *vpkd, void *v) {}
static void combSetArc(void *vpkd, void *v1, void *v2) {
    PKD pkd = (PKD) vpkd;
    PARTICLE * p1 = (PARTICLE *)v1;
    PARTICLE * p2 = (PARTICLE *)v2;
    if (p2->bMarked) p1->bMarked = 1;
    assert( *pkdGroup(pkd,p1) == *pkdGroup(pkd,p2) );
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

/*
** Follow the link for the given particle
*/
static int traverseLink(PKD pkd,struct smGroupArray *ga,int pi,
    int *iPid2,int *iIndex2,int *iMinPartPid,int *iMinPartIndex) {
    MDL mdl = pkd->mdl;
    PARTICLE *p;
    int gid2;
    GHtmpGroupTable *g2;
    int bDone = 0;
    int bRemote = *iPid2 != pkd->idSelf;

    if (!bRemote) p = pkdParticle(pkd,*iIndex2);
    else p = mdlAcquire(pkd->mdl,CID_PARTICLE,*iIndex2,*iPid2);
    gid2 = *pkdGroup(pkd,p);
    g2 = mdlAcquire(pkd->mdl,CID_GROUP,gid2,*iPid2);

    /* Update the minimum group id in case we need it below */
    if (iMinPartPid && (*iPid2 < *iMinPartPid || (*iPid2 == *iMinPartPid && *iIndex2 < *iMinPartIndex))) {
	*iMinPartPid = *iPid2;
	*iMinPartIndex = *iIndex2;
	}

    /*
    ** This group points to the same node, so it is a terminal group.
    ** The particle ids may not match, but the group will point to
    ** a particle that points back to this group.
    */
    if (g2->iPid == *iPid2) {
	ga[pi].id.iPid   = g2->iPid;
	ga[pi].id.iIndex = g2->iIndex;
	bDone = 1;
	}
    /*
    ** If we have linked back to our group on our node, then we are done.
    ** We choice the lowest gid we have seen so far. We also have to mark
    ** the remote particle that started this, because it starts an arc!
    */
    else if (*iPid2 == mdlSelf(mdl) && gid2==pi) { /*LOOP*/
	PARTICLE *p2 = mdlAcquire(pkd->mdl,CID_PARTICLE,ga[pi].id.iIndex,ga[pi].id.iPid);
	p2->bMarked = 1;
	mdlRelease(pkd->mdl,CID_PARTICLE,p2);
	assert(ga[pi].id.iPid != pkd->idSelf);
	ga[pi].id.iPid = *iMinPartPid;
	ga[pi].id.iIndex = *iMinPartIndex;
	bDone = 1;
	}

    /* Follow the link - ignored if bDone is set. */
    *iPid2 = g2->iPid;
    *iIndex2 = g2->iIndex;

    mdlRelease(pkd->mdl,CID_GROUP,g2);
    if (bRemote) mdlRelease(pkd->mdl,CID_PARTICLE,p);

    return bDone;
    }

/*
** Link particles based on density gradients
** After we finish, all "chains" will be complete globally across all domains.
** Chains still need to be joined if their terminating loops are close together.
** Prerequisites: Density and fBall from smooth
*/
int smHopLink(SMX smx,SMF *smf) {
    PKD pkd = smx->pkd;
    MDL mdl = pkd->mdl;
    PARTICLE *p;
    int pi, j, gid, nLoop, nSpur;
    int nGroups, nLocal, nRemote;
    int iIndex1, iIndex2, iMinPartIndex, iPid1, iPid2, iMinPartPid, iParticle;
    struct smGroupArray *ga = smx->ga;
    uint32_t *pl = smx->pl;

    ga[0].iGid = 0;
    ga[0].id.iPid = mdlSelf(mdl);
    ga[0].id.iIndex = -1; /* Sentinel */
    for (pi=0;pi<pkd->nLocal;++pi) {
	p = pkdParticle(pkd,pi);
	*pkdGroup(pkd,p) = 0; /* Ungrouped */
	p->bMarked = p->bSrcActive; /* Used by smooth to determine active particles */
	}
    smSmoothInitialize(smx);

    pkd->nGroups = 0;
    pkd->nLocalGroups = 0;
    pkd->tmpHopGroups = NULL;

    /*
    ** Run through all of the particles, and find a particle to link to
    ** by calculating the gradient and looking in that direction. We make
    ** a link to this particle, then follow the link. Once we have reached
    ** the end of the chain, we consider the next particle (which may have
    ** been already done in which case we skip it).
    */
    nGroups = 1;
    nSpur = nLoop = nRemote = 0;
    for (pi=0;pi<pkd->nLocal;++pi) {
	p = pkdParticle(pkd,iParticle=pi);
	if ( *pkdGroup(pkd,p) > 0 ) continue; /* Already done (below) */
	if ( !pkdIsDstActive(p,0,MAX_RUNG) ) continue;
	ga[nGroups].iGid = nGroups;
	for(;;) {
	    *pkdGroup(pkd,p) = nGroups;
	    /* Need: fDensity and fBall but pl[iParticle] is not yet used */
	    /* fDensity MUST persist but fBall could be stashed into pl[] */
	    smReSmoothSingle(smx,smf,p,pkdBall(pkd,p));
	    /* Done: fBall */
	    ga[nGroups].id = smf->hopParticleLink;
	    /*
	    ** Call mdlCacheCheck to make sure we are making progress!
	    */
	    mdlCacheCheck(mdl);
	    /*
	    ** If we link to a remote node, then this group chain is done for now.
	    ** Later we will merge remote groups.
	    **/
	    if (smf->hopParticleLink.iPid != pkd->idSelf ) { ++nGroups; ++nRemote; break; }
	    else {
		iParticle=pl[iParticle]=smf->hopParticleLink.iIndex;
		p = pkdParticle(pkd,iParticle);
		gid = *pkdGroup(pkd,p);
		/* We link to ourselves: this forms a new "loop" */
		if ( gid == nGroups ) { ++nGroups; ++nLoop; break; }
		/* Ok, some other group: merge this "spur" on to it */
		else if ( gid > 0 ) {
		    for(j=pi; j!=iParticle; j=pl[j]) {
			p = pkdParticle(pkd,j);
			*pkdGroup(pkd,p) = gid;
			}
		    ++nSpur;
		    break;
		    }
		}
	    /* None of the above? Just keep following the links. */
	    }
	}
    smSmoothFinish(smx);
    /*
    ** All chains now terminate at a specific *particle*, either local or remote.
    ** Spurs were automatically joined, and point to their corresponding "loop".
    ** This "loop" particle is guaranteed to be part of the terminating loop.
    ** Remote particles are still pending, so we combine duplicates so there are
    ** fewer remote searches to perform. This results in substantially fewer chains.
    */
    mdlFinishCache(mdl,CID_PARTICLE);
    nGroups = pkdCombineDuplicateGroupIds(pkd,nGroups,ga,0);
    mdlROcache(mdl,CID_PARTICLE,NULL,pkdParticleBase(pkd),pkdParticleSize(pkd),pkdLocal(pkd));
    /*
    ** Now lets go looking for remote loops and see if we are part of them.
    ** The group table still contains PARTICLE links, so we use it to
    ** accelerate the search. We will find three cases:
    **  1. We find a terminating loop. Update the link to this loop.
    **  2. We are part of a terminating loop. Update the link to the lowest
    **     loop id and we are done. All members will do this.
    **  3. We encounter a loop that we are not a part of. Save any member
    **     of the loop. Later, we will fetch the correct id from it.
    */

    /* We will update bMarked (an arc) remotely, so we need a combiner cache */
    mdlFinishCache(mdl,CID_PARTICLE);
    for (pi=0;pi<pkd->nLocal;++pi)
	pkdParticle(pkd,pi)->bMarked = 0; /* Not an arc (yet) */
    mdlCOcache(pkd->mdl,CID_PARTICLE,NULL,
	pkdParticleBase(pkd),pkdParticleSize(pkd),
	pkdLocal(pkd),pkd,initSetArc,combSetArc);
    pkd->tmpHopGroups = mdlMalloc(mdl, nGroups * sizeof(GHtmpGroupTable));
    assert(pkd->tmpHopGroups!=NULL);
    for(pi=1; pi<nGroups; ++pi) {
	pkd->tmpHopGroups[pi].iPid   = ga[pi].id.iPid;
	pkd->tmpHopGroups[pi].iIndex = ga[pi].id.iIndex;
	}
    mdlROcache(mdl,CID_GROUP,NULL,pkd->tmpHopGroups,sizeof(GHtmpGroupTable), nGroups);
    for(pi=1; pi<nGroups; ++pi) {
	/* Completely local loop: any particle can start the "arc" */
	if (ga[pi].id.iPid == pkd->idSelf) {
	    p = pkdParticle(pkd,ga[pi].id.iIndex);
	    p->bMarked = 1;
	    }
	/* Remote: may form a loop, or we may simply be a spur */
	else {
	    iMinPartPid = iPid1 = iPid2 = ga[pi].id.iPid;
	    iMinPartIndex = iIndex1 = iIndex2 = ga[pi].id.iIndex;
	    for(;;) {
		/* We use Floyd's cycle-finding algorithm here */
		if (traverseLink(pkd,ga,pi,&iPid2,&iIndex2,&iMinPartPid,&iMinPartIndex)) break;
		if (traverseLink(pkd,ga,pi,&iPid2,&iIndex2,&iMinPartPid,&iMinPartIndex)) break;
		traverseLink(pkd,ga,pi,&iPid1,&iIndex1,NULL,NULL);
		/*
		** If we see a loop here, then we are not a part of it (we are a spur).
		** We point ourselves at one of the loop members so we can later
		** retrieve the proper terminating particle. We can be pointing
		** to our node of course, but to a different group.
		*/
		if (iPid1==iPid2 && iIndex1 == iIndex2) {
		    ga[pi].id.iPid   = iPid1;
		    ga[pi].id.iIndex = iIndex1;
		    break;
		    }
		}
	    }
	}
    mdlFinishCache(mdl,CID_GROUP);
    /*
    ** By now, we have marked the particle that starts each arc on the node,
    ** and exactly one particle if the group is a local loop.
    ** We now want to mark the rest of the particles in the arcs and loops.
    */
    mdlFinishCache(mdl,CID_PARTICLE);
    for (pi=0;pi<pkd->nLocal;++pi) {
	p = pkdParticle(pkd,pi);
	if (!p->bMarked) continue;
	iIndex1 = pl[pi];
	while(iIndex1!=-1 && !(p=pkdParticle(pkd,iIndex1))->bMarked) {
	    p->bMarked = 1;
	    iIndex1 = pl[iIndex1];
	    }
	}
    /* Done: pl[] is no longer needed */
    pl = NULL;

    nGroups = pkdCombineDuplicateGroupIds(pkd,nGroups,ga,0);
    /*
    ** Now handle deferred groups (ones pointing to, but not part of, loops).
    */
    for(pi=1; pi<nGroups; ++pi) {
	pkd->tmpHopGroups[pi].iPid   = ga[pi].id.iPid;
	pkd->tmpHopGroups[pi].iIndex = ga[pi].id.iIndex;
	}
    mdlROcache(mdl,CID_PARTICLE,NULL,pkdParticleBase(pkd),pkdParticleSize(pkd),pkdLocal(pkd));
    mdlROcache(mdl,CID_GROUP,NULL,pkd->tmpHopGroups,sizeof(GHtmpGroupTable), nGroups);
    nSpur = 0;
    for(pi=1; pi<nGroups; ++pi) {
	PARTICLE *p1 = mdlAcquire(mdl,CID_PARTICLE,ga[pi].id.iIndex,ga[pi].id.iPid);
	int gid1 = *pkdGroup(pkd,p1);
	GHtmpGroupTable *g1 = mdlAcquire(mdl,CID_GROUP,gid1,ga[pi].id.iPid);
	PARTICLE *p2 = mdlAcquire(mdl,CID_PARTICLE,g1->iIndex,g1->iPid);
	int gid2 = *pkdGroup(pkd,p2);
	GHtmpGroupTable *g2 = mdlAcquire(mdl,CID_GROUP,gid2,g1->iPid);
	assert (g2->iPid == g1->iPid);
	ga[pi].id.iPid   = g2->iPid;
	ga[pi].id.iIndex = g2->iIndex;
	mdlRelease(mdl,CID_GROUP,g2);
	mdlRelease(mdl,CID_PARTICLE,p2);
	mdlRelease(mdl,CID_GROUP,g1);
	mdlRelease(mdl,CID_PARTICLE,p1);
	}
    mdlFinishCache(mdl,CID_GROUP);

    /* If we end up point on-node, then we might have to combine the groups */
    for(pi=1; pi<nGroups; ++pi) {
	if (ga[pi].id.iPid==pkd->idSelf) {
	    gid = *pkdGroup(pkd,pkdParticle(pkd,ga[pi].id.iIndex));
	    if (gid!=pi) {
		assert(ga[gid].id.iPid==pkd->idSelf);
		ga[pi].iGid = pi;
		ga[pi].id.iPid = ga[gid].id.iPid;
		ga[pi].id.iIndex = ga[gid].id.iIndex;
		}
	    }
	}
    mdlFinishCache(mdl,CID_PARTICLE);
    nGroups = pkdCombineDuplicateGroupIds(pkd,nGroups,ga,0);
    mdlROcache(mdl,CID_PARTICLE,NULL,pkdParticleBase(pkd),pkdParticleSize(pkd),pkdLocal(pkd));

    /*
    ** Whew. The situation for groups is now one of the following.
    **  1. Groups that point locally form a local terminal loop or arc,
    **     so it must point to itself.
    **  2. Groups that point to a remote group must point to a remote
    **     group that points to itself.
    ** All particles that form a loop (either local or remote) have
    ** been marked (bMarked=1).
    */

    /* Turn particle IDs into group IDs */
    for(pi=1; pi<nGroups; ++pi) {
	pkd->tmpHopGroups[pi].iPid   = ga[pi].id.iPid;
	pkd->tmpHopGroups[pi].iIndex = ga[pi].id.iIndex;
	}
    mdlROcache(mdl,CID_GROUP,NULL,pkd->tmpHopGroups,sizeof(GHtmpGroupTable), nGroups);
    for(pi=1; pi<nGroups; ++pi) {
	if (ga[pi].id.iPid == pkd->idSelf) {
	    p = pkdParticle(pkd,ga[pi].id.iIndex);
	    gid = *pkdGroup(pkd,p);
	    assert(gid==pi); /* I point to myself! */
	    }
	else {
	    p = mdlAcquire(mdl,CID_PARTICLE,ga[pi].id.iIndex,ga[pi].id.iPid);
	    gid = *pkdGroup(pkd,p);
	    mdlRelease(mdl,CID_PARTICLE,p);
	    GHtmpGroupTable *g = mdlAcquire(mdl,CID_GROUP,gid,ga[pi].id.iPid);
	    /* The remote particle must point to itself */
	    assert(g->iPid==ga[pi].id.iPid && g->iIndex==ga[pi].id.iIndex);
	    mdlRelease(mdl,CID_GROUP,g);
	    }
	ga[pi].id.iIndex = gid;
	assert(ga[pi].id.iPid!=pkd->idSelf || ga[pi].id.iIndex ==pi);
	}
    mdlFinishCache(mdl,CID_GROUP);

    pkd->nGroups = nGroups;
    nLocal = nRemote = 0;
    for(pi=1; pi<pkd->nGroups; ++pi) {
	if (ga[pi].id.iPid==mdlSelf(mdl)) nLocal++;
	else nRemote++;
	pkd->tmpHopGroups[pi].iPid   = ga[pi].id.iPid;
	pkd->tmpHopGroups[pi].iIndex = ga[pi].id.iIndex;
	}
    pkd->nLocalGroups = nLocal;

    /* tmpHopGroups is still allocated: freed in pkdHopFinishUp */
    return nLocal;
    }

/*
** Join together chains if their terminating loops are withing a linking length.
**
** OPTIMIZATION: we should mark individual groups as done/not done so we can skip
**               subsequent iterations for those groups.
*/
int smHopJoin(SMX smx,SMF *smf, double dHopTau, int *nLocal) {
    PKD pkd = smx->pkd;
    MDL mdl = pkd->mdl;
    KDN *pRoot = pkdTreeNode(pkd,ROOT);
    struct smGroupArray *ga = smx->ga;
    PARTICLE *p;
    int pi;

    smf->bDone = 1;
    mdlCOcache(mdl,CID_GROUP,NULL,pkd->tmpHopGroups,sizeof(GHtmpGroupTable), pkd->nGroups,
	smf, initJoinLoops, combJoinLoops );
    /*
    ** Check all particles that are part of a loop or an arc (they are marked).
    ** We have contructed a tree with only marked particles, so just check those.
    */
    for (pi=pRoot->pLower;pi<=pRoot->pUpper;++pi) {
	float fBall;
	p = pkdParticle(pkd,pi);
	assert(p->bMarked);
	if (dHopTau<0.0) fBall = -dHopTau * pkdSoft(pkd,p);
	else fBall = dHopTau;
	fBall = fmaxf(fBall,pkdBall(pkd,p)*0.5f);
	smReSmoothSingle(smx,smf,p,fBall);
	}
    mdlFinishCache(mdl,CID_GROUP);

    /* Follow the chains to the end */
    mdlROcache(mdl,CID_GROUP,NULL,pkd->tmpHopGroups,sizeof(GHtmpGroupTable), pkd->nGroups);
    for(pi=1; pi<pkd->nGroups; ++pi) {
	int iPid = pkd->idSelf;
	int iIndex = pi;
	int iNextPid = pkd->tmpHopGroups[iIndex].iPid;
	int iNextIndex = pkd->tmpHopGroups[iIndex].iIndex;
	while(iPid!=iNextPid || iIndex!=iNextIndex) {
	    GHtmpGroupTable *g = mdlAcquire(mdl,CID_GROUP,iNextIndex,iNextPid);
	    iPid = iNextPid;
	    iIndex = iNextIndex;
	    iNextPid = g->iPid;
	    iNextIndex = g->iIndex;
	    mdlRelease(mdl,CID_GROUP,g);
	    }
	pkd->tmpHopGroups[pi].iPid = iPid;
	pkd->tmpHopGroups[pi].iIndex = iIndex;
	}
    mdlFinishCache(mdl,CID_GROUP);

    /* Merge duplicates for the next round */
    for(pi=1; pi<pkd->nGroups; ++pi) {
	ga[pi].iGid   = pi;
	ga[pi].id.iPid   = pkd->tmpHopGroups[pi].iPid;
	ga[pi].id.iIndex = pkd->tmpHopGroups[pi].iIndex;
	}
    mdlFinishCache(mdl,CID_PARTICLE);
    pkd->nGroups = pkdCombineDuplicateGroupIds(pkd,pkd->nGroups,ga,1);
    mdlROcache(mdl,CID_PARTICLE,NULL,pkdParticleBase(pkd),pkdParticleSize(pkd),pkdLocal(pkd));
    int nRemote = 0;
    *nLocal = 0;
    for(pi=1; pi<pkd->nGroups; ++pi) {
	if (ga[pi].id.iPid==mdlSelf(mdl)) (*nLocal)++;
	else nRemote++;
	pkd->tmpHopGroups[pi].iPid   = ga[pi].id.iPid;
	pkd->tmpHopGroups[pi].iIndex = ga[pi].id.iIndex;
	}
    pkd->nLocalGroups = *nLocal;
    return smf->bDone;
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

    printf("%1d:relocate\n",pkd->idSelf);
    hopRelocateGroups(pkd);
    printf("%1d:group stats\n",pkd->idSelf);
    hopCalculateGroupStats(pkd,bPeriodic,dPeriod);
    /* Purge groups with too few particles */
    printf("%1d:purge\n",pkd->idSelf);
   for(i=j=1; i<pkd->nGroups; ++i) {
	assert(pkd->hopGroups[i].nTotal==0 || pkd->hopGroups[i].id.iPid!=pkd->idSelf || pkd->hopGroups[i].nLocal>0);
	if (pkd->hopGroups[i].nLocal==0 || pkd->hopGroups[i].nTotal<nMinGroupSize) gid=0;
	else gid = j++;
	ga[i].id.iPid = pkd->hopGroups[i].id.iPid;
	ga[i].id.iIndex = pkd->hopGroups[i].id.iIndex;
	ga[i].iNewGid = gid;
	}
    printf("%1d:update group ids\n",pkd->idSelf);
    updateGroupIds(pkd,pkd->nGroups,ga,1);
    printf("%1d:final loop\n",pkd->idSelf);
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


int pkdHopFinishUp(PKD pkd,int nMinGroupSize, int bPeriodic, double *dPeriod) {
    MDL mdl = pkd->mdl;
    struct smGroupArray *ga = (struct smGroupArray *)(pkd->pLite);
    int i;

    for(i=1; i<pkd->nGroups; ++i) {
	ga[i].iGid = i;
	ga[i].id.iPid = pkd->tmpHopGroups[i].iPid;
	ga[i].id.iIndex = pkd->tmpHopGroups[i].iIndex;
	}
    mdlFree(mdl,pkd->tmpHopGroups); pkd->tmpHopGroups = NULL;

    /* Create the final group table */
    pkd->hopGroups = mdlMalloc(mdl, pkd->nGroups * sizeof(HopGroupTable));
    assert(pkd->hopGroups!=NULL);
    for(i=0; i<pkd->nGroups; ++i) {
	pkd->hopGroups[i].id.iPid      = ga[i].id.iPid;
	pkd->hopGroups[i].id.iIndex    = ga[i].id.iIndex;
	pkd->hopGroups[i].bNeedGrav    = 1;
	pkd->hopGroups[i].bComplete    = 0;
	}
    pkdPurgeSmallGroups(pkd,nMinGroupSize,bPeriodic,dPeriod);
    pkd->hopSavedRoots = 0;

    /* free: We have allocated pkd->hopGroups */
    return pkd->nLocalGroups;
    }

static void initHopGetRoots(void *vpkd, void *v) {
    HopGroupTable * g = (HopGroupTable *)v;
    g->rmt.iPid = -1;
    g->rmt.iIndex = 0;
    g->iGlobalId = 0;
    }

static void combHopGetRoots(void *vctx, void *v1, void *v2) {
    PKD pkd = (PKD)vctx;
    HopGroupTable * g1 = (HopGroupTable *)v1;
    HopGroupTable * g2 = (HopGroupTable *)v2;

    if (g2->iGlobalId) {
	int gid = g1 - pkd->hopGroups;
	int iRoot = pkd->hopRootIndex[gid]++;
	assert(gid<pkd->nLocalGroups); /* We should only get data for local groups */
	assert(pkd->hopRoots[iRoot].iPid==-1 && pkd->hopRoots[iRoot].iIndex==-1);
	pkd->hopRoots[iRoot].iPid = g2->rmt.iPid;
	pkd->hopRoots[iRoot].iIndex = g2->rmt.iIndex;
	}
    }

void pkdHopTreeBuild(PKD pkd, int nBucket) {
    MDL mdl = pkd->mdl;
    int nDomains = mdlThreads(mdl);
    HopGroupTable * g;
    int i, gid, iRoot;
    int iPid, iLastPid;
    int nRootsTotal;

    /*
    ** The remote groups should be sorted by processor and remote group id,
    ** but we will rely on this behaviour later so we verify it now.
    */
    for(i=2+pkd->nLocalGroups; i<pkd->nGroups; ++i) {
	assert( pkd->hopGroups[i-1].id.iPid<pkd->hopGroups[i].id.iPid
	    || ( pkd->hopGroups[i-1].id.iPid==pkd->hopGroups[i].id.iPid
		&& pkd->hopGroups[i-1].id.iIndex<pkd->hopGroups[i].id.iIndex) );
	}

    /* Setup the buffer for tree roots */
    pkd->hopRootIndex= malloc((pkd->nGroups+1) * sizeof(*pkd->hopRootIndex) );
    assert(pkd->hopRootIndex!=NULL);

    /* Calculate the index into the tree roots table */
    nRootsTotal = 0;
    for(i=1; i<pkd->nGroups; ++i) {
	pkd->hopRootIndex[i] = nRootsTotal;
	pkd->hopGroups[i].iAllRoots = nRootsTotal;
	nRootsTotal += 1 + pkd->hopGroups[i].nRemote;
	}

    pkd->hopRoots= malloc( (nRootsTotal+1) * sizeof(*pkd->hopRoots) );
    assert(pkd->hopRoots!=NULL);

    /* Invalidate all entries (for debug checks) */
    for(iRoot=0; iRoot<=nRootsTotal; ++iRoot)
	pkd->hopRoots[iRoot].iPid = pkd->hopRoots[iRoot].iIndex = -1;

    /* We build a tree of each group of particles - hopGroups[].iTreeRoot is set */
    pkdTreeBuildByGroup(pkd,nBucket);

    /* Add ourselves for local groups */
    for(gid=1; gid<pkd->nLocalGroups; ++gid) {
	iRoot = pkd->hopRootIndex[gid]++;
	assert(iRoot < nRootsTotal);
	assert(pkd->hopRoots[iRoot].iPid==-1 && pkd->hopRoots[iRoot].iIndex==-1);
	pkd->hopRoots[iRoot].iPid = mdlSelf(mdl);
	pkd->hopRoots[iRoot].iIndex = i;
	}

    /*
    ** Let the group master know about the remote tree roots. We can use the
    ** combiner cache because there is at most 1 root per group per node.
    */
    mdlCOcache(mdl,CID_GROUP,NULL,pkd->hopGroups,sizeof(HopGroupTable), pkd->nGroups,
	pkd, initHopGetRoots, combHopGetRoots );
    for(i=1+pkd->nLocalGroups; i<pkd->nGroups; ++i) {
	assert(pkd->hopGroups[i].id.iPid != pkd->idSelf);
//	if (pkd->hopNumRoots[i] == 0) continue; 
	g = mdlAcquire(mdl,CID_GROUP,pkd->hopGroups[i].id.iIndex,pkd->hopGroups[i].id.iPid);
	g->rmt.iPid = pkd->idSelf;
	g->rmt.iIndex = pkd->hopGroups[i].iTreeRoot;
	g->iGlobalId = i;
	pkd->hopGroups[i].iAllRoots = g->iAllRoots; /* We will fix this later! */
	mdlRelease(mdl,CID_GROUP,g);
	}
    mdlFinishCache(mdl,CID_GROUP);

    /*
    ** Now we just need to fetch all of the roots for remote groups
    */
    mdlROcache(mdl,CID_GROUP,NULL,pkd->hopRoots,sizeof(*pkd->hopRoots), nRootsTotal);
    for(gid=1+pkd->nLocalGroups; gid<pkd->nGroups; ++gid) {
	int iIndex = pkd->hopGroups[gid].iAllRoots;
	for(i=0; i<pkd->hopGroups[gid].nRemote+1; ++i) {
	    remoteID *pRoot = mdlFetch(mdl,CID_GROUP,iIndex+i,pkd->hopGroups[gid].id.iPid);
	    int iRoot = pkd->hopRootIndex[gid]++;
	    assert(pkd->hopRoots[iRoot].iPid==-1 && pkd->hopRoots[iRoot].iIndex==-1);
	    pkd->hopRoots[iRoot].iPid = pRoot->iPid;
	    pkd->hopRoots[iRoot].iIndex = pRoot->iIndex;
	    }
	}
    mdlFinishCache(mdl,CID_GROUP);

    /* Fix the indexes */
    nRootsTotal = 0;
    for(gid=1; gid<pkd->nGroups; ++gid) {
	pkd->hopGroups[gid].iAllRoots = nRootsTotal;
	nRootsTotal += 1 + pkd->hopGroups[gid].nRemote;
	}

    for(iRoot=1; iRoot<nRootsTotal; ++iRoot)
	assert(pkd->hopRoots[iRoot].iPid>=0 && pkd->hopRoots[iRoot].iIndex>0);

    for(gid=1; gid<pkd->nGroups; ++gid) {
	int iRootUs = pkd->hopRootIndex[gid];
	for (iRoot=pkd->hopRootIndex[gid]; iRoot<pkd->hopRootIndex[gid+1]; ++iRoot) {
	    assert(iRoot<nRootsTotal);
	    assert(pkd->hopRoots[iRoot].iPid>=0 && pkd->hopRoots[iRoot].iIndex>0);
	    /* Move our local root to the start */
	    if (pkd->hopRoots[iRoot].iPid==pkd->idSelf && iRoot>iRootUs) {
		remoteID t;
		assert(pkd->hopRoots[iRootUs].iPid!=pkd->idSelf);
		t = pkd->hopRoots[iRootUs];
		pkd->hopRoots[iRootUs] = pkd->hopRoots[iRoot];
		pkd->hopRoots[iRoot] = t;
		}
	    }
	}
    free(pkd->hopRootIndex); pkd->hopRootIndex = NULL;
    /* free: We have allocated pkd->hopRoots -- needed for gravity */
    }

typedef struct EnergyElement {
    double dPot;
    double dSqrtKin;
    double dTot;
    int i;
    } EE;


static int cmpEE(const void *p1,const void *p2) {
    EE *a = (EE *)p1;
    EE *b = (EE *)p2;
    if (a->dTot > b->dTot) return 1;
    if (a->dTot < b->dTot) return -1;
    return 0;
    }

static void initMaxHopEnergy(void *vpkd, void *v) {}
static void combMaxHopEnergy(void *vctx, void *v1, void *v2) {
    HopGroupTable * g1 = (HopGroupTable *)v1;
    HopGroupTable * g2 = (HopGroupTable *)v2;
    if (g2->dEnergy > g1->dEnergy) g1->dEnergy = g2->dEnergy;
    }

int pkdHopUnbind(PKD pkd, double dTime, int nMinGroupSize, int bPeriodic, double *dPeriod) {
    MDL mdl = pkd->mdl;
    EE *ee = (EE *)(pkd->pLite);
    int gid, iRoot, i, j, n;
    int nEvaporated;
    PARTICLE *p;
    KDN *pNode;
    HopGroupTable *g;
    double a = csmTime2Exp(pkd->param.csm,dTime);
    double ia = 1.0 / a;
    double a2 = a * a;
    double dv, dv2;
    vel_t *v;
    double dEnergy;

    mdlCOcache(mdl,CID_GROUP,NULL,pkd->hopGroups,sizeof(HopGroupTable), pkd->nGroups,
	pkd, initMaxHopEnergy, combMaxHopEnergy);
    for(gid=1; gid<pkd->nGroups; ++gid) {
	if (pkd->hopGroups[gid].bComplete) continue;
	iRoot = pkd->hopGroups[gid].iTreeRoot;
	pNode = pkdTreeNode(pkd,iRoot);
	n = pNode->pUpper - pNode->pLower + 1;
	/* Calculate kinetic energy & total energy */
	for(i=pNode->pLower; i<=pNode->pUpper; ++i) {
	    p = pkdParticle(pkd,i);
	    assert(*pkdGroup(pkd,p)==gid);
	    v = pkdVel(pkd,p);
	    dv2 = 0.0;
	    for (j=0;j<3;++j) {
		dv = v[j] - pkd->hopGroups[gid].vcom[j];
		dv2 += dv*dv;
		}
	    ee[i].i = i;
	    ee[i].dPot = *pkdPot(pkd,p) * ia;
	    if (pkdIsGas(pkd,p) && pkd->oSph) {  /* TODO: is this correct? */
		SPHFIELDS *pSph = pkdSph(pkd,p);
		ee[i].dPot += pSph->u;
		}
	    ee[i].dSqrtKin = a*sqrt(0.5*dv2);
	    ee[i].dTot = 0.5*a2*dv2 + ee[i].dPot;
	    }
	qsort(ee+pNode->pLower,n,sizeof(EE),cmpEE);
	i = n>100 ? pNode->pUpper-100 : pNode->pLower;
	while( i<=pNode->pUpper && ee[i].dTot <= 0.0 ) ++i;
	if (i>pNode->pUpper) dEnergy = 0.0;
	else  dEnergy = ee[i].dTot;
	pkd->hopGroups[gid].dEnergy = dEnergy;

	if (pkd->hopGroups[gid].id.iPid != pkd->idSelf) {
	    g = mdlAcquire(mdl,CID_GROUP,pkd->hopGroups[gid].id.iIndex,pkd->hopGroups[gid].id.iPid);
	    g->dEnergy = dEnergy;
	    mdlRelease(mdl,CID_GROUP,g);
	    }
	}
    mdlFinishCache(mdl,CID_GROUP);

    mdlROcache(mdl,CID_GROUP,NULL,pkd->hopGroups,sizeof(HopGroupTable), pkd->nGroups);
    for(gid=1+pkd->nLocalGroups; gid<pkd->nGroups; ++gid) {
	if (pkd->hopGroups[gid].bComplete) continue;
	assert(pkd->hopGroups[gid].id.iPid != pkd->idSelf);
	g = mdlAcquire(mdl,CID_GROUP,pkd->hopGroups[gid].id.iIndex,pkd->hopGroups[gid].id.iPid);
	assert(pkd->hopGroups[gid].id.iIndex==g->id.iIndex && pkd->hopGroups[gid].id.iPid==g->id.iPid);
	pkd->hopGroups[gid].dEnergy = g->dEnergy;
	mdlRelease(mdl,CID_GROUP,g);
	}
    mdlFinishCache(mdl,CID_GROUP);
    nEvaporated = 0;
    for(gid=1; gid<pkd->nGroups; ++gid) {
	if (pkd->hopGroups[gid].bComplete) continue;
	iRoot = pkd->hopGroups[gid].iTreeRoot;
	dEnergy = pkd->hopGroups[gid].dEnergy;
	pNode = pkdTreeNode(pkd,iRoot);
	n = pNode->pUpper - pNode->pLower + 1;
	for( i=pNode->pUpper; i>=pNode->pLower && ee[i].dTot >= dEnergy; --i) {
	    p = pkdParticle(pkd,ee[i].i);
	    *pkdGroup(pkd,p) = 0;
	    ++nEvaporated;
	    }
	if (i==pNode->pUpper) pkd->hopGroups[gid].bComplete = 1;
	/* Move evaporated particles to the end */
	else for(i=pNode->pLower; i<=pNode->pUpper; ) {
	    p = pkdParticle(pkd,i);
	    if (*pkdGroup(pkd,p)) ++i;
	    else pkdSwapParticle(pkd,p,pkdParticle(pkd,pNode->pUpper--));
	    }
	}
    pkdPurgeSmallGroups(pkd,nMinGroupSize,bPeriodic,dPeriod);

    return nEvaporated;
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
