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
    struct smGroupArray *g;

    /* Update the group for all local particles */
    for (pi=0;pi<pkd->nLocal;++pi) {
	p = pkdParticle(pkd,pi);
	gid = *pkdGroup(pkd,p);
	if (gid<=0) continue;
	*pkdGroup(pkd,p) = ga[gid].iNewGid;
	}

    /* Now gid has the new position -- a reorder is necessary */
    if (bIndexIsGID) {
	mdlROcache(mdl,CID_GROUP,NULL,ga,sizeof(struct smGroupArray), nGroups);
	for(gid=1; gid<nGroups; ++gid) {
	    if (ga[gid].id.iPid==pkd->idSelf) {
		ga[gid].id.iIndex = ga[gid].iNewGid;
		}
	    else {
		g = mdlFetch(pkd->mdl,CID_GROUP,ga[gid].id.iIndex,ga[gid].id.iPid);
		ga[gid].id.iIndex = g->iNewGid;
		}
	    }
	mdlFinishCache(mdl,CID_GROUP);
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
    pkd->nLocalGroups = 0;
    for(pi=1; pi<nGroups; ++pi) {
	if (ga[pi].id.iPid == pkd->idSelf && pi<nNew) ++pkd->nLocalGroups;
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
    pkd->nLocalGroups = 0;
    for(i=1; i<nGroups; ++i) {
	assert(ga[i].iGid == i);
	assert(ga[i].id.iPid>=0);
	if (ga[i].id.iPid == pkd->idSelf) {
	    ga[i].id.iPid = -1;
	    ++pkd->nLocalGroups;
	    }
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
    nGroups = reorderGroups(pkd,nGroups,ga);
    return(nGroups);
    }

/*
** Merge duplicate groups (iPid/iIndex is the same), and update particle pointers.
** When the remote links are links to groups (basically the names of groups) then
** we want bIndexIsGID set to 1 always so they will be updated. HOP has in early 
** phases remote links to actual particles, and sets this to 0 in those cases.
*/
int pkdGroupCombineDuplicateIds(PKD pkd, int nGroups, struct smGroupArray *ga,int bIndexIsGID) {
    int gid;

    nGroups = renumberGroups(pkd,nGroups,ga);
    updateGroupIds(pkd,nGroups,ga,bIndexIsGID);
    /* We will now put the groups in the new order */
    for(gid=1; gid<nGroups; ++gid) ga[gid].iGid = ga[gid].iNewGid;
    nGroups = reorderGroups(pkd,nGroups,ga);
    return nGroups;
    }


static void initMaxnGroup(void *vpkd, void *v) {}
static void combMaxnGroup(void *vctx, void *v1, void *v2) {
    struct smGroupArray * g1 = (struct smGroupArray *)v1;
    struct smGroupArray * g2 = (struct smGroupArray *)v2;
    /* Remember: nTotal is really MAX(nLocal) here. */
    if (g2->nTotal > g1->nTotal) {
	g1->nTotal = g2->nTotal;
	g1->id.iPid = g2->id.iPid;
	g1->id.iIndex = g2->id.iIndex;
	}
    }

/*
** For each group, find the processor that has the most particles
** and make that processor that master for that group.
*/
int pkdGroupRelocate(PKD pkd,int nGroups,struct smGroupArray *ga) {
    MDL mdl = pkd->mdl;
    PARTICLE *p;
    struct smGroupArray *g;
    int i, gid, n,nLocalGroups;

    /* Count local members of all groups */
    nLocalGroups = 0;
    for(i=0; i<pkd->nGroups; ++i) {
	if (ga[i].id.iPid == pkd->idSelf && i) ++nLocalGroups;
	ga[i].nTotal = 0;
	}
    for (i=0;i<pkd->nLocal;++i) {
	p = pkdParticle(pkd,i);
	gid = *pkdGroup(pkd,p);
	++ga[gid].nTotal;
	}
    /* Now find the real processor with the most particles for each group */
    mdlCOcache(mdl,CID_GROUP,NULL,ga,sizeof(struct smGroupArray),nGroups,
	NULL,initMaxnGroup,combMaxnGroup);
    for(i=1+nLocalGroups; i<nGroups; ++i) {
	assert(ga[i].id.iPid != pkd->idSelf);
	g = mdlVirtualFetch(mdl,CID_GROUP,ga[i].id.iIndex,ga[i].id.iPid);
	g->id.iPid = pkd->idSelf;
	g->id.iIndex = i;
	g->nTotal = ga[i].nTotal;
	}
    mdlFinishCache(mdl,CID_GROUP);
    /* Now update the new group location */
    mdlROcache(mdl,CID_GROUP,NULL,ga,sizeof(struct smGroupArray), pkd->nGroups);
    for(i=1+nLocalGroups; i<nGroups; ++i) {
	g = mdlFetch(mdl,CID_GROUP,ga[i].id.iIndex,ga[i].id.iPid);
	ga[i].id.iPid = g->id.iPid;
	ga[i].id.iIndex = g->id.iIndex;
	}
    mdlFinishCache(mdl,CID_GROUP);
    nGroups = pkdGroupCombineDuplicateIds(pkd,nGroups,ga,1);
    return(nGroups);
    }


static void initTotalnGroup(void *vpkd, void *v) {
    struct smGroupArray * g = (struct smGroupArray *)v;
    g->nTotal = 0;
    }
static void combTotalnGroup(void *vctx, void *v1, void *v2) {
    struct smGroupArray * g1 = (struct smGroupArray *)v1;
    struct smGroupArray * g2 = (struct smGroupArray *)v2;
    g1->nTotal += g2->nTotal;
    }

int pkdGroupCounts(PKD pkd,int nGroups,struct smGroupArray *ga) {
    MDL mdl = pkd->mdl;
    PARTICLE *p;
    struct smGroupArray *g;
    int i, gid;
    int nLocalGroups;

    /*
    ** First count the number of particles in each group.
    */
    nLocalGroups = 0;
    for(i=0; i<nGroups; ++i) {
	if (ga[i].id.iPid == pkd->idSelf && i) ++nLocalGroups;
	ga[i].nTotal = 0;
	}
    for (i=0;i<pkd->nLocal;++i) {
	p = pkdParticle(pkd,i);
	gid = *pkdGroup(pkd,p);
	++ga[gid].nTotal;
	}
    /*
    ** Then add to counts for all groups owned by a remote processor.
    */
    mdlCOcache(mdl,CID_GROUP,NULL,ga,sizeof(struct smGroupArray),nGroups,
	NULL,initTotalnGroup,combTotalnGroup);
    for(i=1+nLocalGroups; i<nGroups; ++i) {
	assert(ga[i].id.iPid != pkd->idSelf);   /* assumes local groups come first, then remotes */
	g = mdlVirtualFetch(mdl,CID_GROUP,ga[i].id.iIndex,ga[i].id.iPid);
	g->nTotal += ga[i].nTotal;
	}
    mdlFinishCache(mdl,CID_GROUP);
    /*
    ** Finally update remote group counts.
    */
    mdlROcache(mdl,CID_GROUP,NULL,ga,sizeof(struct smGroupArray),nGroups);
    for(i=1+nLocalGroups; i<nGroups; ++i) {
	g = mdlFetch(mdl,CID_GROUP,ga[i].id.iIndex,ga[i].id.iPid);
	ga[i].nTotal = g->nTotal;
	}
    mdlFinishCache(mdl,CID_GROUP);

    return(nLocalGroups);
    }


int pkdPurgeSmallGroups(PKD pkd,int nGroups,struct smGroupArray *ga,int nMinGroupSize) {
    int i,j,gid;
    int nLocalGroups;

    nLocalGroups = pkdGroupCounts(pkd,nGroups,ga);
    /* Purge groups with too few particles */
    for(i=j=1;i<nGroups;++i) {
	if (ga[i].nTotal < nMinGroupSize) gid=0;
	else gid = j++;
	ga[i].iNewGid = gid;
	}
   updateGroupIds(pkd,nGroups,ga,1);
   /* We will now put the groups in the new order */
   for(gid=1; gid<nGroups; ++gid) ga[gid].iGid = ga[gid].iNewGid;
   nGroups = reorderGroups(pkd,nGroups,ga);
   return(nGroups);
   }


/*
** pkd->nLocalGroups is the count of groups owned by this processor.
*/
int pkdGroupCountGID(PKD pkd) {
    return pkd->nLocalGroups;
    }


void pkdGroupAssignGID(PKD pkd,uint64_t iStartGID) {
    pkd->iStartGID = iStartGID;
    }


