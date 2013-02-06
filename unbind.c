#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>
#include <assert.h>
#include "pkd.h"
#include "unbind.h"

static void ubBuildLocalGroupTrees(PKD pkd, int nBucket);
static void ubBroadcastTreeRoots(PKD pkd);
static void ubGravity(PKD pkd);
static void ubFindGroupCenters(PKD pkd);
static void ubDoUnbind(PKD pkd);

void ubInitializePLiteParticles(PKD pkd) {
    PLITE *pLite = pkd->pLite;
    PLITE t;
    PARTICLE *p;
    KDN *pNode;
    pBND bnd;
    int i,j;
    int iRoot;

    int grp_compar(const void *a0, const void *b0)
    {
	PLITE *a = (PLITE *)a0;
	PLITE *b = (PLITE *)b0;
	int ga = *pkdGroup(pkd, pkdParticle(pkd, a->i));
	int gb = *pkdGroup(pkd, pkdParticle(pkd, b->i));
	if (ga < gb) return -1;
	if (ga > gb) return +1;
	return 0;
    }

    //if (pkd->nNodes > 0) mdlFinishCache(pkd->mdl,CID_CELL);

    /*
    ** Initialize the temporary particles.
    */
    for (i=0;i<pkd->nLocal;++i) {
	p = pkdParticle(pkd,i);
	for (j=0;j<3;++j) pLite[i].r[j] = p->r[j];
	pLite[i].i = i;
	}

    qsort(pLite, pkd->nLocal, sizeof(*pLite), grp_compar);

    /*
    **It is only forseen that there are 4 reserved nodes at present 0-NULL, 1-ROOT, 2-UNUSED, 3-VAROOT.
    */
    pkd->nNodes = NRESERVED_NODES;
    pkd->nVeryActive = 0;
    /*
    ** Set up the root node for each group.
    */
    struct psGroup *gd = pkd->psGroupTable.pGroup;
    int pLower = gd[0].nLocal;
    for (i=1; i < pkd->psGroupTable.nGroups; i++)
    {
	//assert(gd[i].nLocal>0);
	if (gd[i].nLocal == 0) continue;

	gd[i].nTreeRoots = 1;
	gd[i].treeRoots = malloc(gd[i].nTreeRoots * sizeof(*(gd[i].treeRoots)));

	pkdTreeAllocRootNode(pkd,&iRoot);
	gd[i].treeRoots[0].iPid = pkd->idSelf;
	gd[i].treeRoots[0].iLocalRootId = iRoot;
	gd[i].treeRoots[0].iLocalGroupId = i;
	pNode = pkdTreeNode(pkd,gd[i].treeRoots[0].iLocalRootId);
	pNode->iLower = 0;
	pNode->iParent = 0;
	pNode->pLower = pLower;
	pNode->pUpper = pNode->pLower + gd[i].nLocal - 1;
	pLower += gd[i].nLocal;
	pkdNodeBnd(pkd, pNode, &bnd);
	/*
	** TODO: Make tight bounds. I think this might happen in Create()...
	*/
	for (j=0;j<3;++j) {
	    bnd.fCenter[j] = pkd->bnd.fCenter[j];
	    bnd.fMax[j] = pkd->bnd.fMax[j];
	}
    }

    pkd->nNodesFull = pkd->nNodes;
    //mdlROcache(pkd->mdl,CID_CELL,pkdTreeNodeGetElement,pkd, pkd->iTreeNodeSize,pkd->nNodes);
    //pkdCopyNode(pkd,pkdn,pkdTreeNode(pkd,ROOT));
}

static void ubBuildLocalGroupTrees(PKD pkd, int nBucket)
{
    int i;
    struct psGroup *gd = pkd->psGroupTable.pGroup;

    ubInitializePLiteParticles(pkd);

    int S = gd[0].nLocal;
    for (i=1; i < pkd->psGroupTable.nGroups; i++)
    {
	if (gd[i].nLocal == 0) continue;
	S += gd[i].nLocal;
	BuildTemp(pkd,gd[i].treeRoots[0].iLocalRootId,nBucket);
    }
    assert(S == pkd->nLocal);
    ShuffleParticles(pkd,0);

    for (i=1; i < pkd->psGroupTable.nGroups; i++)
    {
	if (gd[i].nLocal == 0) continue;
	Create(pkd,gd[i].treeRoots[0].iLocalRootId);
    }
}

static void initTreeRoot(void *vpkd, void *a)
{
    struct tree_root_cache_msg *msg = (struct tree_root_cache_msg *)a;
    msg->iPid = 0;
    msg->iLocalId = 0;
    msg->tr.iPid = 0;
    msg->tr.iLocalRootId = 0;
    msg->tr.iLocalGroupId = 0;
}

static void combTreeRoot(void *vpkd, void *a, void *b)
{
    PKD pkd = (PKD)vpkd;
    struct tree_root_cache_msg *msg = (struct tree_root_cache_msg *)b;
    assert(msg->iPid == pkd->idSelf);
    assert(msg->tr.iPid != pkd->idSelf);

    struct psGroup *gd = pkd->psGroupTable.pGroup + msg->iLocalId;

    gd->nTreeRoots++;
    gd->treeRoots = realloc(gd->treeRoots, gd->nTreeRoots * sizeof(*gd->treeRoots));
    assert(gd->treeRoots != NULL);

    gd->treeRoots[gd->nTreeRoots-1] = msg->tr;
}

struct AllBroadcastData
{
    int nDomains;
    int *scounts, *rcounts, *sdispls, *rdispls, *ioffset;
    MDL_Datatype dataType;
};

static void AllBroadcastAlloc(struct AllBroadcastData *d, int nDomains, MDL_Datatype dataType)
{
    d->nDomains = nDomains;
    d->dataType = dataType;
    d->scounts = malloc(sizeof(*d->scounts) * nDomains); assert(d->scounts != NULL);
    d->rcounts = malloc(sizeof(*d->rcounts) * nDomains); assert(d->rcounts != NULL);
    d->sdispls = malloc(sizeof(*d->sdispls) * nDomains); assert(d->sdispls != NULL);
    d->rdispls = malloc(sizeof(*d->rdispls) * nDomains); assert(d->rdispls != NULL);
    d->ioffset = malloc(sizeof(*d->ioffset) * nDomains); assert(d->ioffset != NULL);
}

static void AllBroadcastFree(struct AllBroadcastData *d)
{
    free(d->scounts);
    free(d->rcounts);
    free(d->sdispls);
    free(d->rdispls);
    free(d->ioffset);
}

static int AllBroadcastCount(PKD pkd, struct AllBroadcastData *d)
{
    int nDomains = d->nDomains;
    //int nLocal   = pkd->nLocal;
    int *scounts, *rcounts, *sdispls, *rdispls, *ioffset;
    int nSelf, iSelf, iTarget;
    int i;

    scounts = d->scounts;
    rcounts = d->rcounts;
    sdispls = d->sdispls;
    rdispls = d->rdispls;
    ioffset = d->ioffset;

    mdlAlltoall( pkd->mdl, scounts, 1, MDL_INT, rcounts, 1, MDL_INT );

    ioffset[0] = sdispls[0] = rdispls[0] = 0;
    for(i=1; i<nDomains; i++) {
	ioffset[i] = sdispls[i] = sdispls[i-1] + scounts[i-1];
	rdispls[i] = rdispls[i-1] + rcounts[i-1];
	}

    return rdispls[nDomains-1] + rcounts[nDomains-1]; // + nSelf;
}

static void AllBroadcast(PKD pkd, struct AllBroadcastData *d, void *src, void *dst)
{
    /* Send the particles to their correct processors and update our local count */
    mdlAlltoallv(pkd->mdl,src, d->scounts, d->sdispls, d->dataType, dst, d->rcounts, d->rdispls, d->dataType);
}

static void ubBroadcastTreeRoots(PKD pkd)
{
    int i,j,k,sTotal;

    int nID      = mdlSelf(pkd->mdl);
    MDL_Datatype dataType;
    struct AllBroadcastData bd;
    struct psGroup *gd = pkd->psGroupTable.pGroup;

    struct btr
    {
	struct remote_root_id treeRootId;
	int iLocalDestId;
    };

    mdlTypeContiguous(pkd->mdl, sizeof(struct btr), MDL_BYTE, &dataType);
    mdlTypeCommit(pkd->mdl,&dataType);

    struct tree_root_cache_msg tree_root_msg;
    struct tree_root_cache_msg *remote_msg;

    mdlCOcache(pkd->mdl,CID_TREE_ROOT,NULL,&tree_root_msg,sizeof(tree_root_msg),1,pkd,initTreeRoot,combTreeRoot);
    for (i=1; i < pkd->psGroupTable.nGroups; i++)
    {
	if (gd[i].nLocal == 0) continue;
	if (gd[i].iPid != pkd->idSelf)
	{
	    remote_msg = mdlAquire(pkd->mdl,CID_TREE_ROOT,0, gd[i].iPid);
	    remote_msg->iPid     = gd[i].iPid;
	    remote_msg->iLocalId = gd[i].iLocalId;
	    remote_msg->tr       = gd[i].treeRoots[0];
	    mdlRelease(pkd->mdl,CID_TREE_ROOT,remote_msg);
	    mdlFlushCache(pkd->mdl, CID_TREE_ROOT);
	}
    }
    mdlFinishCache(pkd->mdl,CID_TREE_ROOT);


    AllBroadcastAlloc(&bd, mdlThreads(pkd->mdl), dataType);

    for (i=0; i < bd.nDomains; i++) bd.scounts[i] = 0;
    for (i=1; i < pkd->psGroupTable.nGroups; i++)
    {
	if (gd[i].nLocal == 0) continue;
	if (gd[i].iPid == pkd->idSelf)
	{
	    for (j=1; j < gd[i].nTreeRoots; j++)
	    {
		assert(gd[i].treeRoots[j].iPid != nID);
		int nItems = gd[i].nTreeRoots-1;
		bd.scounts[gd[i].treeRoots[j].iPid] += nItems;
		sTotal += nItems;
	    }
	}
    }
    assert(bd.scounts[nID] == 0);

    int rTotal = AllBroadcastCount(pkd, &bd);

    struct btr *src = malloc(sTotal * sizeof(*src)); assert(src != NULL);
    struct btr *dst = malloc(rTotal * sizeof(*dst)); assert(dst != NULL);

    for (i=1; i < pkd->psGroupTable.nGroups; i++)
    {
	if (gd[i].nLocal == 0) continue;
	if (gd[i].iPid == pkd->idSelf)
	{
	    for (j=1; j < gd[i].nTreeRoots; j++)
	    {
		int pid = gd[i].treeRoots[j].iPid;
		int lid = gd[i].treeRoots[j].iLocalGroupId;
		assert(pid != pkd->idSelf);

		for (k=0; k < gd[i].nTreeRoots; k++)
		{
		    if (gd[i].treeRoots[k].iPid == pid) continue;
		    assert(bd.ioffset[pid] < sTotal);
		    src[bd.ioffset[pid]].treeRootId = gd[i].treeRoots[k];
		    src[bd.ioffset[pid]].iLocalDestId = lid;
		    bd.ioffset[pid]++;
		}
	    }
	}
    }

    AllBroadcast(pkd, &bd,src,dst);
    AllBroadcastFree(&bd);

    int *nTreeRoots = malloc(pkd->psGroupTable.nGroups * sizeof(*nTreeRoots));

    for (i=1; i < pkd->psGroupTable.nGroups; i++) nTreeRoots[i] = gd[i].nLocal==0 ? 0 : 1;
    for (i=0; i < rTotal; i++)
    {
	//fprintf(stderr, "%i\n", dst[i].iLocalDestId);
	assert(dst[i].iLocalDestId != 0);
	assert(dst[i].iLocalDestId < pkd->psGroupTable.nGroups);
	assert(dst[i].treeRootId.iPid != pkd->idSelf);
	nTreeRoots[dst[i].iLocalDestId]++;
    }

    for (i=1; i < pkd->psGroupTable.nGroups; i++)
    {
	if (gd[i].nLocal == 0) continue;
	if (gd[i].iPid != pkd->idSelf)
	{
	    gd[i].treeRoots = realloc(gd[i].treeRoots, nTreeRoots[i] * sizeof(*gd[i].treeRoots));
	    gd[i].nTreeRoots = 1;
	}
    }

    for (i=0; i < rTotal; i++)
    {
	int lid = dst[i].iLocalDestId;
	int offs = gd[lid].nTreeRoots;
	assert(offs != 0);
	gd[lid].treeRoots[offs] = dst[i].treeRootId;
	gd[lid].nTreeRoots++;
    }
    for (i=1; i < pkd->psGroupTable.nGroups; i++)
    {
	if (gd[i].nLocal == 0) continue;
	if (gd[i].iPid != pkd->idSelf)
	{
	    if (gd[i].nTreeRoots != nTreeRoots[i])
		fprintf(stderr, "gd[%i].nTreeRoots=%i  !=  nTreeRoots[%i]=%i\n", i, gd[i].nTreeRoots,i, nTreeRoots[i]);
	    assert(gd[i].nTreeRoots == nTreeRoots[i]);
	}
    }

    //for (i=1; i < pkd->psGroupTable.nGroups; i++)
	//fprintf(stderr, "iPid %i  Group %i  nRoots %i  nLocal %i %s\n", pkd->idSelf, gd[i].iGlobalId, gd[i].nTreeRoots, gd[i].nLocal, gd[i].iPid==pkd->idSelf?"(Owner)":"");

    free(nTreeRoots);
    free(src);
    free(dst);
}

static void ubGravity(PKD pkd)
{
    double dFlop = 0.0;
    double dPartSum = 0.0;
    double dCellSum = 0.0;
    double dTime = 1;
    double dThetaMin = 0.8;
    double dThetaMax = 1.0;
    int nActive;
    int nGroup = 64;

    int i;

    for (i=0; i < pkd->nLocal; i++)
    {
	*pkdPot(pkd, pkdParticle(pkd, i)) = 0;
    }

    mdlROcache(pkd->mdl,CID_PARTICLE,NULL,pkdParticleBase(pkd),pkdParticleSize(pkd), pkdLocal(pkd));
    mdlROcache(pkd->mdl,CID_CELL,pkdTreeNodeGetElement,pkd, pkd->iTreeNodeSize,pkd->nNodes);
    pkdGravWalkGroups(pkd,dTime,nGroup,dThetaMin,dThetaMax,&dFlop,&dPartSum,&dCellSum);
    mdlFinishCache(pkd->mdl,CID_PARTICLE);
    mdlFinishCache(pkd->mdl,CID_CELL);
}

#if 0
static void initGroups(void *vpkd, void *a)
{
    struct psGroup *g = (struct psGroup *)a;
    g->fRMSRadius = 0;
    g->fMass = 0;
    g->nTotal = 0;
    g->nLocal = 0;
    g->fMass_com = 0;
    g->vcom[0] = 0;
    g->vcom[1] = 0;
    g->vcom[2] = 0;
}

static void combGroups(void *vpkd, void *a, void *b)
{
    struct psGroup * g1 = (struct psGroup *)a;
    struct psGroup * g2 = (struct psGroup *)b;
    g1->nTotal += g2->nTotal;
    g1->fMass  += g2->fMass;
    if (g2->fRMSRadius > g1->fRMSRadius)
	g1->fRMSRadius = g2->fRMSRadius;
    g1->fMass_com  += g2->fMass_com;
    g1->vcom[0] += g2->vcom[0];
    g1->vcom[1] += g2->vcom[1];
    g1->vcom[2] += g2->vcom[2];
}
#endif

static void ubDoUnbind(PKD pkd)
{
    int i;
    struct psGroup *g = pkd->psGroupTable.pGroup;

    int nUnbound = 0;

    for (i=0; i < pkd->nLocal; i++)
    {
	PARTICLE *p = pkdParticle(pkd, i);
	int gid = *pkdGroup(pkd, p);

	if (gid == 0) continue;

	double *v = pkdVel(pkd, p);
	double v2;

	v2 = pow(v[0] - g[gid].vcom[0], 2) 
	   + pow(v[1] - g[gid].vcom[1], 2) 
	   + pow(v[2] - g[gid].vcom[2], 2);

	FLOAT T = 0.5 * v2;

	if (T > -*pkdPot(pkd,p))
	{
	    g[gid].nLocal--;
	    nUnbound++;
	    *pkdGroup(pkd,p) = 0;
	    g[0].nLocal++;
	}
    }
    fprintf(stderr, "nUnbound=%i\n", nUnbound);
}

void ubUnbind(PKD pkd, int nBucket)
{
    int i;
    if (pkd->psGroupTable.nGroups > 0)
    {
	int i;
	for (i=0; i < 10; i++)
	{
	    fprintf(stderr, "Pass %i *************************\n", i+1);
	    psdUpdateGroupProperties(pkd);
	    ubBuildLocalGroupTrees(pkd, nBucket);
	    ubBroadcastTreeRoots(pkd);
	    ubGravity(pkd);
	    ubDoUnbind(pkd);
	}
        psdUpdateGroupProperties(pkd);
    }
    mdlROcache(pkd->mdl,CID_CELL,pkdTreeNodeGetElement,pkd, pkd->iTreeNodeSize,pkd->nNodes);
}


