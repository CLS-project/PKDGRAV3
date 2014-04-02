#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <math.h>
#include <limits.h>
#include <assert.h>
#include "smooth.h"
#include "psd.h"
#include "pkd.h"
#include "knn6d.h"
#include "rbtree.h"
#include <sys/stat.h>
#include "qsort.h"

#define USE_POTENTIAL_GRADIENT 0

#if USE_POTENTIAL_GRADIENT
#define PROP(pkd,p) (-*pkdPot((pkd), (p)))
#define PROPGRAD  psdPotentialGrad
#else
#define PROP(pkd,p) (p->fDensity)
#define PROPGRAD  psdDensityGrad
#endif

#define PSM(i) (psx->psm[i])


/*
** The Epanechnikov kernel
*/
static inline double ekernel(const double u2) {
    return (0 <= u2 && u2 <= 1) * (1-u2);
}

static inline double grad_ekernel(const double u) {
    return (fabs(u) <= 1) * (-2*u);
}

/* 
** A spline kernel (unused)
*/
static inline float spline_kernel(float u) {
    if (0 < u && u <= 0.5)
	return 1 - (6*u*u) * (1 + u);
    else if (0.5 < u && u <= 1)
	return 2 * pow(1-u, 3);
    return 0;
}

static inline float grad_spline_kernel(float x, float u) {
    if (0 < u && u <= 0.5)
	return (6*x) * (2 - 3*u);
    else if (0.5 < u && u <= 1)
	return -3*x/u * (1 - 2*u + u*u);
    return 0;
}

/*
** The smoothed phase-space density.
*/
static inline float psdDensity(PKD pkd, PARTICLE *p,int nSmooth, PQ6 *nnList, FLOAT *rscale, FLOAT *vscale) {
    double ih2,r2,fDensity,fMass;
    int i;

    ih2 = BALL2(p);
    assert(ih2 > 0);
    fDensity = 0.0;

    fDensity = -ekernel(0) * pkdMass(pkd,p);
    fDensity += pkdMass(pkd, p) * ekernel(pow(1 * 6 / (1+6.), 2));
    for (i=0;i<nSmooth;++i) {
	fMass = pkdMass(pkd,nnList[i].pPart);
	r2 = nnList[i].fDist2/ih2;
	fDensity += ekernel(r2) * fMass;
	}

    float V0=1, V1=1;
    for (i=0; i < 3; i++) {
	if (rscale[i] > 0)
	{
	    V0 *= p->fBall;
	    V1 *= rscale[i];
	}
    }

    for (i=0; i < 3; i++) {
	if (vscale[i] > 0)
	{
	    V0 *= p->fBall;
	    V1 *= vscale[i];
	}
    }

    return 0.77403670 * fDensity * V1 / V0;

}

/*
** Compute the normalized density gradient. Also correct the E0 error
** to reduce the noise.
*/
static inline void psdDensityGrad(PKD pkd, PSX psx, KNN6D knn, int pid, FLOAT *fDensityGrad, int normalize) {
    int i,j;
    FLOAT *rscale, *vscale;
    PARTICLE *p = pkdParticle(pkd, pid);
    const FLOAT fBall = p->fBall;
    const FLOAT fDensity = p->fDensity;

    rscale = PSM(pid).rscale;
    vscale = PSM(pid).vscale;

    for (j=0; j < 6; j++) fDensityGrad[j] = 0;

    for (i=0;i < psx->nSmooth;++i)
    {
	if (knn->pq[i].fDist2 > 0)
	{
	    const double r = sqrt(knn->pq[i].fDist2) / fBall;
	    const double fMass = pkdMass(pkd, knn->pq[i].pPart);
	    const double c = fMass * (fDensity - knn->pq[i].pPart->fDensity) / knn->pq[i].pPart->fDensity;
	    for (j=0; j < 6; j++)
	    {
		const double dx = knn->pq[i].dr[j]/fBall;
		fDensityGrad[j] += c * grad_ekernel(dx);
	    }
	}
    }

    if (normalize)
    {
	double L = 0;
	for (i=0; i < 6; i++) L += pow(fDensityGrad[i], 2);
	L = sqrt(L);
	for (i=0; i < 6; i++) fDensityGrad[i] /= L;
    }
}

int psdInitialize(PKD pkd, PSX psx, int nSmooth,int bPeriodic)
{
    psx->nSmooth = nSmooth;
    psx->bPeriodic = bPeriodic;

    psx->psm = malloc(pkd->nLocal * sizeof(*psx->psm)); assert(psx->psm != NULL);
    psx->knn = malloc(sizeof(*psx->knn)); assert(psx->knn != NULL);

    knn6dInitialize(pkd, psx->knn, nSmooth, bPeriodic);
    psx->knn->psm = psx->psm;
    return 0;
}

void psdFinish(PKD pkd, PSX psx)
{
    knn6dFree(pkd, psx->knn);
    free(psx->psm); psx->psm = NULL;
}


/*
** Compute the normalized potential gradient. Also correct the E0 error
** to reduce the noise.
*/
#if 0
void psdPotentialGrad(PKD pkd, PSX psx, int pid, FLOAT *fPotentialGrad) {
    int i,j;
    FLOAT *rscale, *vscale;
    PARTICLE *p = pkdParticle(pkd, pid);
    const FLOAT fBall = p->fBall;
    const FLOAT fPot = *pkdPot(pkd, p);
    const FLOAT fDensity = p->fDensity;

    rscale = PSM(pid).rscale;
    vscale = PSM(pid).vscale;

    for (j=0; j < 6; j++) fPotentialGrad[j] = 0;

    for (i=0;i < psx->nSmooth;++i)
    {
	double r = sqrt(knn->pq[i].fDist2) / fBall;
	if (r > 0)
	{
	    for (j=0; j < 6; j++)
	    {
		double dx = knn->pq[i].dr[j]/fBall;
		double c = (fPot - *pkdPot(pkd, knn->pq[i].pPart)) / knn->pq[i].pPart->fDensity;
		fPotentialGrad[j] += (pkdMass(pkd, knn->pq[i].pPart) * c) * grad_ekernel(dx);
	    }
	}
    }

    double L = 0;
    for (i=0; i < 6; i++) L += pow(fPotentialGrad[i], 2);
    L = sqrt(L);
    for (i=0; i < 6; i++) fPotentialGrad[i] /= L;
}
#endif

#if 0
void psdSmooth(PSX smx) {
    PKD pkd = smx->pkd;
    PARTICLE *p;
    int pi;
    int N;
    printf("psdSmooth (file)\n");

    //FILE *fp = fopen("/home/itp/jonathan/zbox3b/MAD/MockHalos/Hernquist/main.den", "r"); assert(fp != NULL);
    //FILE *fp = fopen("subsubhalo_nfw.den6", "r"); assert(fp != NULL);
    FILE *fp = fopen("A-5.den6", "r"); assert(fp != NULL);
    //FILE *fp = fopen("test.den", "r"); assert(fp != NULL);
    fscanf(fp, "%i\n", &N);
    assert(N == pkd->nLocal);
    float *den = malloc(N * sizeof(*den)); assert(den != NULL);
    for (pi=0;pi<pkd->nLocal;++pi) {
	fscanf(fp, "%f\n", &den[pi]);
    }
    fclose(fp);

    for (pi=0;pi<pkd->nLocal;++pi) {
	p = pkdParticle(pkd,pi);
	p->fDensity = den[p->iOrder];
    }
    free(den);

#if 1
    //fp = fopen("test.ball", "r"); assert(fp != NULL);
    //fp = fopen("/home/itp/jonathan/zbox3b/MAD/MockHalos/Hernquist/main.ball", "r"); assert(fp != NULL);
    fp = fopen("A-5.ball", "r"); assert(fp != NULL);
    fscanf(fp, "%i\n", &N);
    assert(N == pkd->nLocal);
    for (pi=0;pi<pkd->nLocal;++pi) {
	p = pkdParticle(pkd,pi);
	fscanf(fp, "%f\n", &p->fBall);
    }
    fclose(fp);
#endif
    printf("psdSmooth (file) finished\n");
}

#else


/*
** Compute the smoothed phase-space densities for all particles.
*/
void psdSmooth(PKD pkd, PSX psx) {
    PARTICLE *p;
    FLOAT fBall;
    int pi,i,bDone=0;

    int pqSize = psx->knn->pqSize;

    mdlROcache(pkd->mdl,CID_PARTICLE,NULL, pkdParticleBase(pkd),pkdParticleSize(pkd), pkd->nLocal);
    mdlROcache(pkd->mdl,CID_CELL,pkdTreeNodeGetElement,pkd, pkd->iTreeNodeSize,pkd->nNodes);

    for (pi=0;pi<pkd->nLocal;++pi) 
    {
	p = pkdParticle(pkd,pi);
	/* if ( !pkdIsDstActive(p,0,MAX_RUNG) ) continue; */

	knn6d(pkd, psx->knn, pi, &p->fBall, pi==0);

	/*
	** Apply smooth function to the neighbor list.
	*/
	p->fDensity = psdDensity(pkd, p,pqSize,psx->knn->pq,PSM(pi).rscale,PSM(pi).vscale);
	/*
	** Call mdlCacheCheck to make sure we are making progress!
	*/
	mdlCacheCheck(pkd->mdl);
    }

    knn6dFinish(pkd, psx->knn);

    mdlFinishCache(pkd->mdl,CID_PARTICLE);
    mdlFinishCache(pkd->mdl,CID_CELL);

}
#endif

/*
** For each of the neighbors stored in psx->pq, compute the arc length
** between the vector a and the position of the neighbor.
*/
void calc_arclen(int nSmooth, PQ6 *pq, double *a, double *arclen)
{
    int i,pj;
    double L = 0;
    double b[6];

    for (pj=0; pj < nSmooth; pj++)
    {
	L = 0;
	for (i=0; i < 6; i++)
	{
	    b[i] = pq[pj].dr[i];
	    L += pow(b[i], 2);
	}
	L = sqrt(L);
	for (i=0; i < 6; i++)
	    b[i] /= L;

	arclen[pj] = 0;
	for (i=0; i < 6; i++)
	    arclen[pj] += a[i] * b[i];
	arclen[pj] = acos(arclen[pj]);
	arclen[pj] *= L;
    }
}

/*
** Form local groups by forming chains of particles that end at local density maxima.
**
** Each particle links to its neighbor that lies best along the density gradient and
** is denser. Best is defined as having the smallest arc length between the gradient
** vector and the relative neighbor position vector.
**
** If a neighbor already belongs to a chain then the current chain attaches to it.
** If no neighbor is denser then the chain terminates.
** A chain may also terminate if a neighbor exists but is on another processor. In
** this case, the terminal particle and its neighbor form a bridge which will be 
** joined at a later stage after the local groups are built.
*/
void psdSmoothLink(PKD pkd, PSX psx) {
    PARTICLE *p;
    PQ6 *pq;
    int64_t pi,i;
    int64_t pj;
    int64_t idx;
    int32_t *piGroup;

#define TEMP_S_INCREASE 100
    int *C; NEW_STACK(C, TEMP_S_INCREASE);
    int *G; NEW_STACK(G, TEMP_S_INCREASE);

#ifdef _MSC_VER
    int sorted_nbrs[256];
    double arclen[256];
#else
	int sorted_nbrs[psx->nSmooth];
	double arclen[psx->nSmooth];
#endif
	double fDensityGrad[6];

    struct bridge *B; NEW_STACK(B, TEMP_S_INCREASE);
    struct bridge bi;
#if 0
    float den_max = 0;
    for (pi=0;pi<pkd->nLocal;pi++) 
    {
	if (pkdParticle(pkd, pi)->fDensity > den_max)
	    den_max = pkdParticle(pkd, pi)->fDensity;
    }

#endif
    for (pi=0;pi<pkd->nLocal;pi++) 
    {
	p = pkdParticle(pkd,pi);
	*pkdGroup(pkd,p) = 0;
    }

    pkd->psGroupTable.nGroups = 0;
    pkd->nMaxRm = 0;

    psx->nBridges = 0;

    int nGroups = 1;
    int32_t trial_group;
    int first_time = 1;

#if 0
    char fname[256];
    sprintf(fname, "link.%i", pkd->idSelf);
    FILE *fp = fopen(fname, "w");
#endif

    mdlROcache(pkd->mdl,CID_PARTICLE,NULL, pkdParticleBase(pkd),pkdParticleSize(pkd), pkd->nLocal);
    mdlROcache(pkd->mdl,CID_CELL,pkdTreeNodeGetElement,pkd, pkd->iTreeNodeSize,pkd->nNodes);

    pq = psx->knn->pq;

    int nPeaks = 0;
    int nSingles = 0;
    for (idx=0;idx<pkd->nLocal;++idx) 
    {
	pi = idx;

	PARTICLE *p0 = pkdParticle(pkd,pi);
	if (*pkdGroup(pkd,p0) != 0) continue;

	pkdAccel(pkd, p0)[0] = 0;
	pkdAccel(pkd, p0)[1] = 0;
	pkdAccel(pkd, p0)[2] = 0;

	int chain_len = 0;

	trial_group = -1;

	assert(STACK_EMPTY(C));
	while (trial_group == -1)
	{
	    EXTEND_STACK(C);

	    PUSH(C,pi);
	    chain_len++;

	    p = pkdParticle(pkd,pi);

	    assert(*pkdGroup(pkd,p) != -1);
	    *pkdGroup(pkd,p) = -1;

	    /* Find our neighbors */
	    knn6d(pkd, psx->knn, pi, NULL, first_time);
	    first_time = 0;

	    for (pj=0; pj < psx->nSmooth; pj++)
		sorted_nbrs[pj] = pj;

	    PROPGRAD(pkd, psx, psx->knn, pi, fDensityGrad, 1);
	    calc_arclen(psx->nSmooth, pq, fDensityGrad, arclen);
#define cmp_arclen(a,b) (arclen[*a] < arclen[*b])
	    QSORT(int, sorted_nbrs, psx->nSmooth, cmp_arclen);

	    int bridge = 0;

	    /* Find the first sorted particle with a density greater than ours */
	    float max_den=PROP(pkd, pkdParticle(pkd, pi));
	    PARTICLE *next_p = NULL;
	    int max_den_j=-1;
	    PQ6 *nbr;
	    for (pj=0; pj < psx->nSmooth; pj++)
	    {
		nbr = pq + sorted_nbrs[pj];

		if (PROP(pkd, nbr->pPart) > max_den)
		{
/* Use this to find the most dense particle greater than ours */
#if 0
		    max_den = nbr->pPart->fDensity;
		    max_den_j = pj;
		}
	    }

	    while (max_den_j != -1)
	    {
		pj = max_den_j;
		nbr = psx->pq + sorted_nbrs[pj];
		{

#endif
		    bridge = nbr->iPid != pkd->idSelf;

		    if (bridge)
		    {
			bi.iPid   = nbr->iPid;
			bi.iIndex = nbr->iIndex;
		    }
		    else
		    {
			int32_t nbr_grp = *pkdGroup(pkd, nbr->pPart);

			/* This would be false if we formed a loop. We shouldn't be able to form loops. */
			assert(nbr_grp != -1);

			if (nbr_grp > 0) /* We've found a local chain. Take that group. */
			    trial_group = nbr_grp;

			next_p = nbr->pPart;
			pi = nbr->iIndex;
		    }

		    /* fprintf(fp, "LINK %ld %ld\n", p->iOrder, nbr->pPart->iOrder); */

		    if (pkd->oAcceleration)
		    {
			pkdAccel(pkd, p)[0] = nbr->pPart->r[0] - p->r[0];
			pkdAccel(pkd, p)[1] = nbr->pPart->r[1] - p->r[1];
			pkdAccel(pkd, p)[2] = nbr->pPart->r[2] - p->r[2];
		    }

		    break;
		}
	    }

	    /* We didn't find a new particle to link to. We must be at a peak. */
	    if (next_p == NULL)
	    {
		if (!bridge)
		{
		    if (chain_len == 1)
			nSingles++;

		    /* We still couldn't find something to link to. Just create a new group. */
		    if (trial_group == -1 && !bridge)
		    {
			trial_group = nGroups++;
			EXTEND_STACK(G);
			PUSH(G,pi);
		    }

		    nPeaks++;
		}


		if (bridge)
		{
		    trial_group = nGroups++;
		    EXTEND_STACK(G);
		    PUSH(G,pi);

		    bi.local_gid = trial_group;
		    bi.remote_gid = -1;
		    bi.done = 0;
		    EXTEND_STACK(B);
		    PUSH(B, bi);
		    psx->nBridges++;
		}
	    }

#if 0
	    pkdAccel(pkd, p)[0] = fDensityGrad[0];
	    pkdAccel(pkd, p)[1] = fDensityGrad[1];
	    pkdAccel(pkd, p)[2] = fDensityGrad[2];
#endif

	    /*
	    ** Call mdlCacheCheck to make sure we are making progress!
	    */
	    mdlCacheCheck(pkd->mdl);
	}

	assert(trial_group != -1);
	assert(trial_group != 0);

	/* Assign particles to the (new) group */
	while (!STACK_EMPTY(C))
	{
	    int pid = POP(C);
	    *pkdGroup(pkd, pkdParticle(pkd,pid)) = trial_group;
	}
    }

#if 0
    fclose(fp);
#endif
#if 0
    fprintf(stdout, "%i] nPeaks is %i\n",   pkd->idSelf, nPeaks);
    fprintf(stdout, "%i] nBridges is %i\n",   pkd->idSelf, psx->nBridges);
    fprintf(stdout, "%i] nSingles is %i\n",  pkd->idSelf,  nSingles);
    fprintf(stdout, "%i] nGroups is %i\n",  pkd->idSelf,  nGroups);
#endif

    knn6dFinish(pkd, psx->knn);

    mdlFinishCache(pkd->mdl,CID_PARTICLE);
    mdlFinishCache(pkd->mdl,CID_CELL);

    /*
    ** Get rid of the stack of bridges and store them in the psx context.
    */
    psx->bridges = malloc(psx->nBridges * sizeof(*(psx->bridges)));
    i=0;
    while (!STACK_EMPTY(B)) { psx->bridges[i++] = POP(B); }
    assert(STACK_EMPTY(B));

    /* When will we free this? - in msrDeletePSGroups() */
    pkd->psGroupTable.nGroups = nGroups;
    pkd->psGroupTable.pGroup = mdlMalloc(pkd->mdl, pkd->psGroupTable.nGroups * sizeof(struct psGroup));

    /*
    ** Create the local group table
    */
    memset(pkd->psGroupTable.pGroup, 0, pkd->psGroupTable.nGroups * sizeof(struct psGroup));
    struct psGroup *gd = pkd->psGroupTable.pGroup;

    gd[0].iGlobalId = 0;
    gd[0].iLocalId = 0;
    gd[0].iPid = pkd->idSelf;
    gd[0].fDensity = 0;
    gd[0].nTotal = pkd->nLocal;
    gd[0].bridge = 0;

    for (i=nGroups-1; i > 0; i--)
    {
	int pid = POP(G);
	PARTICLE *p = pkdParticle(pkd,pid);
	assert(*pkdGroup(pkd,p) == i);
	gd[i].iGlobalId = i;
	gd[i].iLocalId = i;
	gd[i].iPid = pkd->idSelf;
	gd[i].fDensity = PROP(pkd, p);
	gd[i].r[0] = p->r[0];
	gd[i].r[1] = p->r[1];
	gd[i].r[2] = p->r[2];
	gd[i].fMass += pkdMass(pkd, p);
	gd[i].bridge = 0;
	gd[i].dup = 0;
    }
    assert(STACK_EMPTY(G));

    for (i=0; i < psx->nBridges; i++)
	gd[psx->bridges[i].local_gid].bridge = 1;

    FREE_STACK(C);
    FREE_STACK(G);
    FREE_STACK(B);
}

/*
** Join groups that span multiple processors. This will be called many times from master.
**
** Group information is propagated down the chains from the processor that owns the group.
** Ownership is defined as the processor where the density peak is located. The group id
** from the owner defines the group id on other processors.
**
** One edge case occurs when a group chain returns to the owners processor (possible many times).
** This creates a duplicate entry in the group table. Such occurences are noted and later
** the particles in duplicated groups will be reassigned to a single group.
**
** Group information within the particles is not considered here except to retrieve the group id
** of the remote group in a bridge.
*/
int psdJoinBridges(PKD pkd, PSX psx) {
    int done = 1;
    int i,j;
    MDL mdl = pkd->mdl;
    PARTICLE *p, *p_remote;
    struct bridge *bi = psx->bridges;
    struct psGroup *gd = pkd->psGroupTable.pGroup;

    struct store
    {
	int i;
	struct psGroup gd;
    };

#define TEMP_S_INCREASE 100
    struct store store;
    struct store *S; NEW_STACK(S, TEMP_S_INCREASE);

    assert(pkd->oGroup);

    for (i=0; i < pkd->nLocal; i++)
    {
	PARTICLE *p = pkdParticle(pkd,i);
	if (*pkdGroup(pkd, p) == 0)
	{
	    fprintf(stderr, "%d] Particle %d has group 0\n", pkd->idSelf, i);
	    assert(0);
	}
    }


    mdlROcache(pkd->mdl,CID_PARTICLE,NULL, pkdParticleBase(pkd),pkdParticleSize(pkd), pkd->nLocal);
    mdlROcache(mdl,CID_GROUP,NULL,pkd->psGroupTable.pGroup,sizeof(struct psGroup), pkd->psGroupTable.nGroups);

    int nBridgesLeft = 0;

    /*
    ** Find once the group id of the remote particle in a bridge.
    */
    for (i=0; i < psx->nBridges; i++)
    {
	if (!bi[i].done && bi[i].remote_gid == -1)
	{
	    assert(bi[i].iPid != pkd->idSelf);
	    p_remote = mdlAquire(mdl,CID_PARTICLE,bi[i].iIndex,bi[i].iPid);
	    bi[i].remote_gid = *pkdGroup(pkd, p_remote);
	    if (bi[i].remote_gid == 0)
	    {
		fprintf(stderr, "%i] %i %i %ld id=%llu\n",
		    pkd->idSelf, bi[i].remote_gid, bi[i].iPid, bi[i].iIndex,
		    p_remote->iOrder); 
	    }
	    assert(bi[i].remote_gid != 0);
	    mdlRelease(mdl,CID_PARTICLE,p_remote);
	}
    }

    /*
    ** Copy the remote group information from higher up the chain. Store this on a stack rather
    ** than overwrite our local data because another processor might be trying to read it.
    */
    for (i=0; i < psx->nBridges; i++)
    {
	if (bi[i].done) continue;

	struct psGroup *new_gd = mdlAquire(mdl, CID_GROUP, bi[i].remote_gid, bi[i].iPid);
	store.i = i;
	store.gd = *new_gd;
	EXTEND_STACK(S);
	PUSH(S, store);
	if (store.gd.iLocalId == 0)
	{
	    fprintf(stderr, "%i] store.iLocalId=%i  %i %i %ld\n", pkd->idSelf, store.gd.iLocalId, bi[i].remote_gid, bi[i].iPid, bi[i].iIndex); 
	}
	assert(store.gd.iLocalId != 0);
	mdlRelease(mdl,CID_GROUP,new_gd);
    }

    mdlFinishCache(mdl,CID_GROUP);
    mdlFinishCache(mdl,CID_PARTICLE);

    /*
    ** Now that the cache is closed, we can safely update our local group table.
    */
    while (!STACK_EMPTY(S))
    {
	store = POP(S);
	i = store.i;

	assert(!bi[i].done);

	gd = pkd->psGroupTable.pGroup + bi[i].local_gid;

	if(gd->bridge == 0 && store.gd.bridge == 1)
	    assert(0);

	*gd = store.gd;
	gd->dup = (gd->iPid == pkd->idSelf) * gd->iLocalId;

	bi[i].done = !gd->bridge;
	done = done && bi[i].done;

	nBridgesLeft += !bi[i].done;
    }


#if 0
    for (i=0; i < psx->nBridges; i++)
    {
	if (bi[i].done) continue;

	gd = pkd->psGroupTable.pGroup + bi[i].local_gid;

	if(gd->bridge == 0 && store[i].bridge == 1)
	    assert(0);

	*gd = store[i];
	gd->dup = (gd->iPid == pkd->idSelf) * gd->iLocalId;

	bi[i].done = !gd->bridge;
	done = done && bi[i].done;

	nBridgesLeft += !bi[i].done;
    }
#endif

    if (nBridgesLeft == 0 && psx->nBridges != 0)
    {
	psx->nBridges = 0;
	free(psx->bridges);
	psx->bridges = NULL;
    }

    FREE_STACK(S);

#if 0
    if (nBridgesLeft > 0)
	fprintf(stdout, "%i] %i bridges left\n", pkd->idSelf, nBridgesLeft);
#endif

    return done;
}

int psdJoinGroupBridges(PKD pkd, PSX psx) {
    int done = 1;
#if 0
    int i,j;
    MDL mdl = pkd->mdl;
    PSGD *gd = pkd->psGroupTable.pGroup;

    struct store
    {
	int i;
	PSGD gd;
    };

#define TEMP_S_INCREASE 100
    struct store store;
    struct store *S; NEW_STACK(S, TEMP_S_INCREASE);

    assert(pkd->oGroup);

    //mdlROcache(pkd->mdl,CID_PARTICLE,NULL, pkdParticleBase(pkd),pkdParticleSize(pkd), pkd->nLocal);
    mdlROcache(mdl,CID_GROUP,NULL,pkd->psGroupTable.pGroup,sizeof(PSGD), pkd->psGroupTable.nGroups);

    int nBridgesLeft = 0;

    /*
    ** Copy the remote group information from higher up the chain. Store this on a stack rather
    ** than overwrite our local data because another processor might be trying to read it.
    */
    for (i=1; i < pkd->psGroupTable.nGroups; i++)
    {
	if (gd[i].iPid != pkd->idSelf) continue;
	if (!gd[i].bridge) continue;

	PSGD *new_gd;

	int is_remote = gd[i].sp.nbr.iPid != pkd->idSelf;
	if (is_remote)
	    new_gd = mdlAquire(mdl, CID_GROUP, gd[i].sp.nbr.iLocalId, gd[i].sp.nbr.iPid);
	else
	    new_gd = gd + gd[i].sp.nbr.iLocalId;

	store.i = i;
	store.gd = *new_gd;
	EXTEND_STACK(S);
	PUSH(S, store);
	if (store.gd.iLocalId == 0)
	{
	    fprintf(stderr, "%i] store.iLocalId=%i  %i %ld\n", pkd->idSelf, store.gd.iLocalId, gd[i].sp.nbr.iLocalId, gd[i].sp.nbr.iPid);
	}
	assert(store.gd.iLocalId != 0);

	if (is_remote)
	    mdlRelease(mdl,CID_GROUP,new_gd);
    }

    mdlFinishCache(mdl,CID_GROUP);
    //mdlFinishCache(mdl,CID_PARTICLE);

    /*
    ** Now that the cache is closed, we can safely update our local group table.
    */
    while (!STACK_EMPTY(S))
    {
	store = POP(S);
	i = store.i;
	assert(gd[i].bridge);

	gd[i] = store.gd;
	gd[i].dup = (gd[i].iPid == pkd->idSelf) * gd[i].iLocalId;

	done = done && !gd[i].bridge;

	nBridgesLeft += gd[i].bridge;
    }

    FREE_STACK(S);

#endif
    return done;
}

/*
** Count unique local groups. Do not count duplicates that will later be removed. 
*/
int psdCountLocalGroups(PKD pkd) {
    int i;
    int nGroups = pkd->psGroupTable.nGroups;
    struct psGroup *gd = pkd->psGroupTable.pGroup;

    int count = 0;
    for (i=1; i < nGroups; i++)
    {
	assert(gd[i].iLocalId != 0);
	if (!gd[i].dup && gd[i].iPid == pkd->idSelf)
	    count++;
    }

    return count;
}

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

    if (g1->iLocalId == 0) return;

    if (g2->fRMSRadius > g1->fRMSRadius)
	g1->fRMSRadius = g2->fRMSRadius;
    g1->fMass_com  += g2->fMass_com;
    g1->vcom[0] += g2->vcom[0];
    g1->vcom[1] += g2->vcom[1];
    g1->vcom[2] += g2->vcom[2];
}

void psdUpdateGroupProperties(PKD pkd)
{
    int i;

    struct psGroup *gd = pkd->psGroupTable.pGroup;

    for (i=0; i < pkd->psGroupTable.nGroups; i++)
    {
	gd[i].nTotal = 0;
	gd[i].nLocal = 0;
	gd[i].fMass = 0;
	gd[i].fRMSRadius = 0;
	gd[i].dup = 0;
	gd[i].fMass_com = 0;
	gd[i].vcom[0] = 0;
	gd[i].vcom[1] = 0;
	gd[i].vcom[2] = 0;
    }

    for (i=0; i < pkd->nLocal; i++)
    {
	PARTICLE *p = pkdParticle(pkd, i);
	int gid = *pkdGroup(pkd, p);

	gd[gid].nLocal++;

	if (gid == 0) continue;

	double d = pow(gd[gid].r[0] - p->r[0],2)
		 + pow(gd[gid].r[1] - p->r[1],2)
		 + pow(gd[gid].r[2] - p->r[2],2);

	d = sqrt(d);
	if (d > gd[gid].fRMSRadius)
	    gd[gid].fRMSRadius = d;

	if (p->fDensity > 0.9*gd[gid].fDensity)
	{
	    FLOAT fMass = pkdMass(pkd,p);
	    double *v   = pkdVel(pkd,p);
	    gd[gid].fMass_com += fMass;
	    gd[gid].vcom[0] += fMass * v[0];
	    gd[gid].vcom[1] += fMass * v[1];
	    gd[gid].vcom[2] += fMass * v[2];
	}
    }

    /*
    ** Combine local data across domains.
    */
    mdlCOcache(pkd->mdl,CID_GROUP,NULL,pkd->psGroupTable.pGroup,sizeof(struct psGroup), pkd->psGroupTable.nGroups,pkd,initGroups,combGroups);
    for (i=1; i < pkd->psGroupTable.nGroups; i++)
    {
	if (gd[i].iPid != pkd->idSelf)
	{
	    struct psGroup *remote_gd = mdlAquire(pkd->mdl,CID_GROUP, gd[i].iLocalId, gd[i].iPid);
	    *remote_gd = gd[i];
	    mdlRelease(pkd->mdl,CID_GROUP,remote_gd);
	}
    }
    mdlFinishCache(pkd->mdl,CID_GROUP);

    mdlROcache(pkd->mdl,CID_GROUP,NULL,pkd->psGroupTable.pGroup,sizeof(struct psGroup), pkd->psGroupTable.nGroups);
    for (i=1; i < pkd->psGroupTable.nGroups; i++)
    {
	if (gd[i].iPid != pkd->idSelf)
	{
	    struct psGroup *remote_gd = mdlAquire(pkd->mdl,CID_GROUP, gd[i].iLocalId, gd[i].iPid);
	    assert(remote_gd->iPid == gd[i].iPid);
	    assert(remote_gd->iLocalId == gd[i].iLocalId);
	    gd[i].nTotal = remote_gd->nTotal;
	    gd[i].fDensity = remote_gd->fDensity;
	    gd[i].fMass = remote_gd->fMass;
	    gd[i].fRMSRadius = remote_gd->fRMSRadius;
	    gd[i].v[0] = remote_gd->v[0];
	    gd[i].v[1] = remote_gd->v[1];
	    gd[i].v[2] = remote_gd->v[2];
	    gd[i].vcom[0] = remote_gd->vcom[0];
	    gd[i].vcom[1] = remote_gd->vcom[1];
	    gd[i].vcom[2] = remote_gd->vcom[2];
	    gd[i].r[0] = remote_gd->r[0];
	    gd[i].r[1] = remote_gd->r[1];
	    gd[i].r[2] = remote_gd->r[2];
	    gd[i].rcom[0] = remote_gd->rcom[0];
	    gd[i].rcom[1] = remote_gd->rcom[1];
	    gd[i].rcom[2] = remote_gd->rcom[2];
	    gd[i].fMass_com = remote_gd->fMass_com;
	    mdlRelease(pkd->mdl,CID_GROUP,remote_gd);
	}
    }
    mdlFinishCache(pkd->mdl,CID_GROUP);

    for (i=1; i < pkd->psGroupTable.nGroups; i++)
    {
	gd[i].vcom[0] /= gd[i].fMass_com;
	gd[i].vcom[1] /= gd[i].fMass_com;
	gd[i].vcom[2] /= gd[i].fMass_com;

	gd[i].rcom[0] /= gd[i].fMass_com;
	gd[i].rcom[1] /= gd[i].fMass_com;
	gd[i].rcom[2] /= gd[i].fMass_com;
    }
}

static void _CompactGroupTable(PKD pkd)
{
#if 0
    int i;
    struct psGroup *gd = pkd->psGroupTable.pGroup;

    /*
    ** Now we need to compact the group table. Many entries will point to the same
    ** remote group and we will later need to have one unique entry for each remote
    ** group. To get this, we sort a copy of the table by the global group id,
    ** being careful to partition the table into local and remote groups at the
    ** beginning and end of the table, respectively. Once sorted the table is
    ** compacted and all the particles must update their group pointer to the new
    ** location in the table.
    */
    nested function are not portable! Fix me!

    int groupid_cmp(const void *a0, const void *b0)
    {
	struct psGroup *a = (struct psGroup *)a0;
	struct psGroup *b = (struct psGroup *)b0;
	if (a->iPid == pkd->idSelf && b->iPid != pkd->idSelf) return -1;
	if (b->iPid == pkd->idSelf && a->iPid != pkd->idSelf) return +1;

	if (a->iGlobalId < b->iGlobalId) return -1;
	if (a->iGlobalId > b->iGlobalId) return +1;
	return 0;
    }

    /* Copy only none duplicate groups. We are now free of this hassle. */
    int nGroups=0;
    struct psGroup *psGroupData = mdlMalloc(pkd->mdl, pkd->psGroupTable.nGroups * sizeof(struct psGroup));
    for (i=0; i < pkd->psGroupTable.nGroups; i++)
    {
	if (!gd[i].dup)
	{
	    psGroupData[nGroups] = gd[i];
	    psGroupData[nGroups].dup = i;
	    nGroups++;
	}
    }
    qsort(psGroupData, nGroups, sizeof(struct psGroup), groupid_cmp);

    /* 
    ** Compact the table. Use the iLocalId of the original table to store
    ** the new location.
    */
    int gid = psGroupData[0].iGlobalId;
    int k=0;
    for (i=1; i < nGroups; i++)
    {
	assert(i > k);
	if (psGroupData[i].iGlobalId != gid)
	{
	    k++;
	    psGroupData[k] = psGroupData[i];
	    gid = psGroupData[i].iGlobalId; 
	}

	gd[psGroupData[i].dup].iLocalId = k;
    }
    nGroups = k+1;

    /* Now update the particles with the new table entry locations. */
    for (i=0; i < pkd->nLocal; i++)
    {
	PARTICLE *p = pkdParticle(pkd, i);
	*pkdGroup(pkd, p) = gd[*pkdGroup(pkd, p)].iLocalId;
    }

    mdlFree(pkd->mdl, pkd->psGroupTable.pGroup);
    pkd->psGroupTable.nGroups = nGroups;
    pkd->psGroupTable.pGroup = mdlMalloc(pkd->mdl, pkd->psGroupTable.nGroups * sizeof(struct psGroup));
    memcpy(pkd->psGroupTable.pGroup, psGroupData, pkd->psGroupTable.nGroups * sizeof(struct psGroup));
#else
    assert(0);
#endif
}


/*
** Update the local group table with global ids. The offset into the range of globally
** unique ids comes from master.
*/
void psdAssignGlobalIds(PKD pkd, int offs, int count)
{
    int64_t i;
    struct psGroup *gd = pkd->psGroupTable.pGroup;

    int new_local_id=1;

    struct gdstore
    {
	int i;
	struct psGroup gd;
    };

#define TEMP_S_INCREASE 100
    struct gdstore *G; NEW_STACK(G, TEMP_S_INCREASE);
    struct gdstore gdstore;

    /*
    ** Update local groups with a unique global id
    */
    for (i=1; i < pkd->psGroupTable.nGroups; i++)
    {
	if (!gd[i].dup && gd[i].iPid == pkd->idSelf)
	{
	    gd[i].iGlobalId = new_local_id + offs;
	    gd[i].iLocalId = new_local_id;
	    new_local_id++;
	}
    }
    assert(new_local_id == count+1);

    /*
    ** To avoid hassles later, we fix up particles that belong to 
    ** duplicate group entries.
    */
    for (i=0; i < pkd->nLocal; i++)
    {
	PARTICLE *p = pkdParticle(pkd, i);
	int gid = *pkdGroup(pkd, p);

	if (gd[gid].dup != 0)
	{
	    *pkdGroup(pkd, p) = gd[gid].dup;
	}
    }

    /*************/
    /*************/

    /*
    ** Now bring over the global ids from remote groups. 
    */

    mdlROcache(pkd->mdl,CID_GROUP,NULL,pkd->psGroupTable.pGroup,sizeof(struct psGroup), pkd->psGroupTable.nGroups);
    for (i=1; i < pkd->psGroupTable.nGroups; i++)
    {
	if (gd[i].iPid != pkd->idSelf)
	{
	    assert(gd[i].dup == 0);
	    struct psGroup *new_gd = mdlAquire(pkd->mdl,CID_GROUP, gd[i].iLocalId, gd[i].iPid);
	    assert(new_gd->iPid == gd[i].iPid);
	    assert(new_gd->dup == 0);
	    gd[i] = *new_gd;
	    mdlRelease(pkd->mdl,CID_GROUP,new_gd);
	}
    }
    mdlFinishCache(pkd->mdl,CID_GROUP);

    _CompactGroupTable(pkd);

    psdUpdateGroupProperties(pkd);

    /*************/
    /*************/

#if 0

    /*
    ** Compute some group quantities like total mass and radius.
    */
    for (i=1; i < pkd->psGroupTable.nGroups; i++)
    {
	gd[i].nTotal = 0;
	gd[i].nLocal = 0;
	gd[i].fMass = 0;
	gd[i].fRMSRadius = 0;
	gd[i].dup = 0;
    }

    for (i=0; i < pkd->nLocal; i++)
    {
	PARTICLE *p = pkdParticle(pkd, i);
	int gid = *pkdGroup(pkd, p);

	gd[gid].nLocal++;
	gd[gid].nTotal++;
	gd[gid].fMass += pkdMass(pkd,p);

	double d = pow(gd[gid].r[0] - p->r[0],2)
		 + pow(gd[gid].r[1] - p->r[1],2)
		 + pow(gd[gid].r[2] - p->r[2],2);

	d = sqrt(d);

	//fprintf(stderr, "@@ %i/%i %e %e\n", gid, pkd->psGroupTable.nGroups, d, gd[gid].fRMSRadius);
	if (d > gd[gid].fRMSRadius)
	    gd[gid].fRMSRadius = d;
    }

    /*
    ** Combine local data across domains.
    */
    mdlCOcache(pkd->mdl,CID_GROUP,NULL,pkd->psGroupTable.pGroup,sizeof(struct psGroup), pkd->psGroupTable.nGroups,pkd,initGroups,combGroups);
    //fprintf(stderr, "%i] nGroups %i\n", pkd->idSelf, pkd->psGroupTable.nGroups);

    for (i=1; i < pkd->psGroupTable.nGroups; i++)
    {
	if (gd[i].iPid != pkd->idSelf)
	{
	    assert(gd[i].dup == 0);
	    struct psGroup *remote_gd = mdlAquire(pkd->mdl,CID_GROUP, gd[i].iLocalId, gd[i].iPid);
	    remote_gd->nTotal = gd[i].nTotal;
	    remote_gd->fMass  = gd[i].fMass;
	    if (gd[i].fRMSRadius > remote_gd->fRMSRadius)
		remote_gd->fRMSRadius = gd[i].fRMSRadius;
	    mdlRelease(pkd->mdl,CID_GROUP,remote_gd);
	}
    }
    mdlFinishCache(pkd->mdl,CID_GROUP);

#endif

#if 0
    mdlROcache(pkd->mdl,CID_GROUP,NULL,pkd->psGroupTable.pGroup,sizeof(PSGD), pkd->psGroupTable.nGroups);

    for (i=0; i < pkd->nLocal; i++)
    {
	PARTICLE *p = pkdParticle(pkd, i);
	int gid = *pkdGroup(pkd, p);
	if (gd[gid].dup != 0)
	    gid = gd[gid].dup;

	double d = pow(gd[gid].r[0] - p->r[0],2) 
		 + pow(gd[gid].r[1] - p->r[1],2)
		 + pow(gd[gid].r[2] - p->r[2],2);

	double fRMSRadius;
	if (gd[gid].iPid == pkd->idSelf)
	{
	    fRMSRadius = gd[gid].fRMSRadius;
	}
	else
	{
	    PSGD *remote_gd = mdlAquire(pkd->mdl,CID_GROUP, gd[gid].iLocalId, gd[gid].iPid);
	    fRMSRadius = remote_gd->fRMSRadius;
	    mdlRelease(pkd->mdl,CID_GROUP,remote_gd);
	}

	p->fDensity = 1 - pow(d / fRMSRadius,2);
    }
    mdlFinishCache(pkd->mdl,CID_GROUP);

#endif
    FREE_STACK(G);
}

void psdSetGlobalId(PKD pkd)
{
    int i;
    struct psGroup *gd = pkd->psGroupTable.pGroup;

    /*
    ** Update particles with their global group id.
    */
    for (i=0; i < pkd->nLocal; i++)
    {
	PARTICLE *p = pkdParticle(pkd, i);
	*pkdGroup(pkd, p) = gd[*pkdGroup(pkd, p)].iGlobalId;
    }
}

static void initSaddle(void *vpkd, void *a)
{
    struct saddle_point_buffer * spb = (struct saddle_point_buffer *)a;
    spb->nSaddlePoints = 0;
}

static void combSaddle(void *vpkd, void *a, void *b)
{
    struct saddle_point_buffer * spb = (struct saddle_point_buffer *)b;
    int i,j;

    struct psGroup *g = ((PKD)vpkd)->psGroupTable.pGroup + spb->iLocalId;

    struct saddle_point_list *spl = &((PKD)vpkd)->saddle_points;

    for (i=0; i < spb->nSaddlePoints; i++)
    {
	int better_saddle = 0;

	struct saddle_point *grp_sp;
	for (j=0; j < g->nSaddlePoints; j++)
	{
	    grp_sp = spl->sp + g->sp[j];
	    if (spb->sp[i].nbr.iGlobalId == grp_sp->nbr.iGlobalId)
	    {
		better_saddle = spb->sp[i].fDensity > grp_sp->fDensity;
		break;
	    }
	}

	/* If the saddle point doesn't exist in the group then we need to
	** extend the list of saddle points for that group and the global
	** list of saddle points.
	*/
	if (j == g->nSaddlePoints)
	{
	    g->nSaddlePoints++;
	    g->sp = realloc(g->sp, g->nSaddlePoints * sizeof(*g->sp));
	    assert(g->sp != NULL);

	    if (spl->n+1 > spl->size)
	    {
		spl->size += 32;
		spl->sp = realloc(spl->sp, spl->size * sizeof(*spl->sp));
		assert(spl->sp != NULL);
	    }

	    g->sp[j] = spl->n;
	    grp_sp = spl->sp + g->sp[j];
	    spl->n++;

	    better_saddle = 1;
	}

	if (better_saddle)
	    *grp_sp = spb->sp[i];
    }

#if 0
    if (g2->sp.fDensity > g1->sp.fDensity)
    {
	assert(g2->bridge);
	assert(g2->sp.nbr_grp_fDensity > 0);

	g1->bridge = g2->bridge;
	g1->sp     = g2->sp;
#if 0
	g1->sp.fDensity = g2->sp.fDensity;
	g1->sp.iLocalId = g2->sp.iLocalId;
	g1->sp.iPid     = g2->sp.iPid;
	//if (g2->sp.nbr_fDensity > g1->sp.nbr_fDensity)
	{
	    g1->sp.nbr_fDensity     = g2->sp.nbr_fDensity;
	    g1->sp.nbr_grp_iLocalId = g2->sp.nbr_grp_iLocalId;
	    g1->sp.nbr_grp_iPid     = g2->sp.nbr_grp_iPid;
	    g1->sp.nbr_grp_fDensity = g2->fDensity;
	}
#endif
    }
#endif
}

static void initMarkSaddlePoint(void *vpkd, void *p0)
{
    PKD pkd = (PKD)vpkd;
    PARTICLE *p = (PARTICLE *)p0;
    float *a = pkdAccel(pkd, p);
    a[0] = a[1] = a[2] = 0;
}

static void combMarkSaddlePoint(void *vpkd, void *a, void *b)
{
    PKD pkd = (PKD)vpkd;
    PARTICLE *p1 = (PARTICLE *)a;
    PARTICLE *p2 = (PARTICLE *)b;
    float *a1 = pkdAccel(pkd, p1);
    float *a2 = pkdAccel(pkd, p2);
    if (a2[0] != 0 || a2[1] != 0 || a2[2] != 0)
    {
	a1[0] = a2[0];
	a1[1] = a2[1];
	a1[2] = a2[2];
	*pkdGroup(pkd, p1) = *pkdGroup(pkd, p2);
    }
}

void psdMergeNoisyGroups(PKD pkd, PSX psx)
{
#if 0
    PQ6 *pq;
    int64_t i,j;
    int pi, pj;
    struct psGroup *gd = pkd->psGroupTable.pGroup;
    int c;
    double fDensityGrad[6];

#define TEMP_S_INCREASE 100
    int *G; NEW_STACK(G, TEMP_S_INCREASE);
    int *S; NEW_STACK(S, TEMP_S_INCREASE);

    for (i=1; i < pkd->psGroupTable.nGroups; i++)
    {
	if (!gd[i].dup)
	{
	    gd[i].nSaddlePoints = 0;
	    gd[i].sp = NULL;
	    //gd[i].foreign_fraction = 0;
	}

	if (gd[i].iPid == pkd->idSelf && !gd[i].dup)
	{

	    //gd[i].sp.fDensity = 0;
	    //gd[i].sp.iPid = -1;
	    //gd[i].sp.iLocalId = -1;
	    //gd[i].sp.nbr_fDensity = 0;
	    //gd[i].sp.nbr_grp_iLocalId = -1;
	    //gd[i].sp.nbr_grp_iPid = -1;
	    gd[i].bridge = 0;
	}
    }

    pkd->saddle_points.n = 0;
    pkd->saddle_points.size = 0;
    pkd->saddle_points.sp = NULL;

    mdlROcache(pkd->mdl,CID_CELL,pkdTreeNodeGetElement,pkd, pkd->iTreeNodeSize,pkd->nNodes);
    mdlROcache(pkd->mdl,CID_PARTICLE,NULL, pkdParticleBase(pkd),pkdParticleSize(pkd), pkd->nLocal);
    mdlROcache(pkd->mdl,CID_GROUP,NULL,pkd->psGroupTable.pGroup,sizeof(struct psGroup), pkd->psGroupTable.nGroups);

    for (pi=0;pi<pkd->nLocal;++pi) 
    {
	PARTICLE *p = pkdParticle(pkd, pi);

	knn6d(pkd, psx->knn, pi, NULL, pi==0);

	struct psGroup *gdentry = gd + *pkdGroup(pkd, p);

	//int remote_group = gdentry->iPid != pkd->idSelf;
	
	//if (remote_group)
	    //gdentry = mdlAquire(pkd->mdl,CID_GROUP, gdentry->iLocalId, gdentry->iPid);

	int gid = gdentry->iGlobalId;

	int is_saddle = 0;

	//gdentry->sp.nbr_fDensity = 0;

	//PROPGRAD(pkd, psx, pi, fDensityGrad, 0);

	c=0;
	for (pj=0; pj < psx->nSmooth; pj++)
	{
	    int nbr_in_other_group = 0;
	    struct psGroup *nbr_gd;
	    int nbr_gid   = *pkdGroup(pkd, psx->knn->pq[pj].pPart);
	    int remote_nbr = pq[pj].iPid != pkd->idSelf;

	    /* 
	    ** Get the neighbor's group data. If it's a remote particle then lookup the data in the table
	    ** of the machine with the particle. It may not be the master entry for that group but all
	    ** we need is the global id of the group, which will be correct.
	    */
	    if (remote_nbr)
		nbr_gd = mdlAquire(pkd->mdl,CID_GROUP, nbr_gid, psx->knn->pq[pj].iPid);
	    else
		nbr_gd = gd + nbr_gid;

	    assert(nbr_gd->dup == 0);

	    /* Is the neighbor in another group? */
	    nbr_in_other_group = nbr_gd->iGlobalId != gid;

	    /* If so, and it is close enough, we could be a saddle point. */
	    if (nbr_in_other_group)
	    {

#if 0
		c=0;
		for (i=0; i < 6; i++)
		{
		    fprintf(stderr, "%e ", fDensityGrad[i]);
		    if (fabs(fDensityGrad[i]) < 1e-1)
			c++;
		}
		fprintf(stderr, "\n");
		if (c == 6)
		    is_saddle = 1;
#endif

#if 1
		if (sqrt(psx->knn->pq[pj].fDist2) < 0.01)
		{
		    /* Look through the list of saddle points for a group matching the neighbor group.
		     * Update or create the saddle point with the current particle if it is better.
		     */
		    c = 1;

		    struct saddle_point *sp;
		    for (i=0; i < gdentry->nSaddlePoints; i++)
		    {
			sp = pkd->saddle_points.sp + gdentry->sp[i];
			if (sp->nbr.iGlobalId == nbr_gd->iGlobalId)
			{
			    is_saddle = p->fDensity > sp->fDensity;
			    break;
			}
		    }

		    if (i == gdentry->nSaddlePoints)
		    {
			gdentry->nSaddlePoints++;
			gdentry->sp = realloc(gdentry->sp, gdentry->nSaddlePoints * sizeof(*gdentry->sp));
			assert(gdentry->sp != NULL);

			if (pkd->saddle_points.n+1 > pkd->saddle_points.size)
			{
			    pkd->saddle_points.size += 32;
			    pkd->saddle_points.sp = realloc(pkd->saddle_points.sp, pkd->saddle_points.size * sizeof(*pkd->saddle_points.sp));
			    assert(pkd->saddle_points.sp != NULL);
			}

			gdentry->sp[i] = pkd->saddle_points.n;
			sp = pkd->saddle_points.sp + gdentry->sp[i];
			pkd->saddle_points.n++;

			is_saddle = 1;
		    }

		    if (is_saddle)
		    {
			sp->iLocalId = pi;
			sp->iPid     = pkd->idSelf;
			sp->fDensity = p->fDensity;

			sp->nbr.iGlobalId = nbr_gd->iGlobalId;
			sp->nbr.iLocalId  = nbr_gd->iLocalId;
			sp->nbr.iPid      = nbr_gd->iPid;
			sp->nbr.fDensity  = nbr_gd->fDensity;

			sp->owner.iGlobalId = gdentry->iGlobalId;
			sp->owner.iLocalId  = gdentry->iLocalId;
			sp->owner.iPid      = gdentry->iPid;
			sp->owner.fDensity  = gdentry->fDensity;
		    }

#if 0
		    is_saddle = p->fDensity > gdentry->sp.fDensity;

		    if (is_saddle && psx->pq[pj].pPart->fDensity > gdentry->sp.nbr_fDensity)
		    {
			gdentry->sp.nbr_fDensity     = psx->pq[pj].pPart->fDensity;
			gdentry->sp.nbr_grp_iLocalId = nbr_gd->iLocalId;
			gdentry->sp.nbr_grp_iPid     = nbr_gd->iPid;
			gdentry->sp.nbr_grp_fDensity = nbr_gd->fDensity;
		    }
#endif
		}
#endif
	    }

	    //if (remote_nbr)
		//mdlRelease(pkd->mdl,CID_GROUP,nbr_gd);
	}

#if 0
	if (is_saddle)
	{
	    EXTEND_STACK(G);
	    PUSH(G,pi);
	}
#endif

#if 0
	if (0)//is_saddle)
	{
	    //fprintf(stderr, "?\n");
	    gdentry->sp.fDensity = p->fDensity;
	    gdentry->sp.iLocalId = pi;
	    gdentry->sp.iPid = pkd->idSelf;
	    gdentry->bridge = 1;
	    assert(gdentry->sp.nbr_fDensity > 0);

	    assert(!(gdentry->sp.nbr_grp_iPid == gdentry->iPid && gdentry->sp.nbr_grp_iLocalId == gdentry->iLocalId));
	}
#endif

#if 0
	if (c==0)
	{
	    EXTEND_STACK(G);
	    PUSH(G,pi);
	}
#endif

	//if (remote_group)
	 //   mdlRelease(pkd->mdl,CID_GROUP,gdentry);

	mdlCacheCheck(pkd->mdl);
    }

    knn6dFinish(pkd, psx->knn);
    mdlFinishCache(pkd->mdl,CID_PARTICLE);
    mdlFinishCache(pkd->mdl,CID_GROUP);
    mdlFinishCache(pkd->mdl,CID_CELL);

    /*
    ** Send back the saddle points we found locally to the owner of the group.
    */

    pkd->saddle_points.buf = mdlMalloc(pkd->mdl, sizeof(*pkd->saddle_points.buf));
    mdlCOcache(pkd->mdl,CID_SADDLE_BUF,NULL,pkd->saddle_points.buf,sizeof(*pkd->saddle_points.buf),1,pkd,initSaddle,combSaddle);
    for (i=1; i < pkd->psGroupTable.nGroups; i++)
    {
	struct saddle_point_buffer *sp_buffer;

	if (!gd[i].dup)
	{
	    int remote_group = gd[i].iPid != pkd->idSelf;
	    
	    if (remote_group)
	    {
		assert(gd[i].nSaddlePoints < 32);

		int offs = 0;
		while (offs < gd[i].nSaddlePoints)
		{
		    sp_buffer = mdlAquire(pkd->mdl,CID_SADDLE_BUF,0, gd[i].iPid);

		    sp_buffer->iLocalId  = gd[i].iLocalId;
		    sp_buffer->iGlobalId = gd[i].iGlobalId;
		    sp_buffer->nSaddlePoints = gd[i].nSaddlePoints - offs;
		    if (sp_buffer->nSaddlePoints > 32)
			sp_buffer->nSaddlePoints = 32;

		    for (j=0; j < sp_buffer->nSaddlePoints; j++, offs++)
			sp_buffer->sp[j] = pkd->saddle_points.sp[gd[i].sp[offs]];
		    mdlRelease(pkd->mdl,CID_SADDLE_BUF,sp_buffer);
#ifdef MPI_VERSION
		    mdlFlushCache(pkd->mdl, CID_SADDLE_BUF);
#else
		    assert(0);
#endif
		}
	    }
	}
    }
    mdlFinishCache(pkd->mdl,CID_SADDLE_BUF);

    /*
    ** 
    */
    mdlCOcache(pkd->mdl,CID_SADDLE_BUF,NULL,pkd->saddle_points.buf,sizeof(*pkd->saddle_points.buf),1,pkd,initSaddle,combSaddle);
    for (i=1; i < pkd->psGroupTable.nGroups; i++)
    {
	struct saddle_point_buffer *sp_buffer;

	if (!gd[i].dup)
	{
	    int remote_group = gd[i].iPid != pkd->idSelf;
	    
	    if (remote_group)
	    {
		assert(gd[i].nSaddlePoints < 32);

		int offs = 0;
		while (offs < gd[i].nSaddlePoints)
		{
		    sp_buffer = mdlAquire(pkd->mdl,CID_SADDLE_BUF,0, gd[i].iPid);

		    sp_buffer->iLocalId  = gd[i].iLocalId;
		    sp_buffer->iGlobalId = gd[i].iGlobalId;
		    sp_buffer->nSaddlePoints = gd[i].nSaddlePoints - offs;
		    if (sp_buffer->nSaddlePoints > 32)
			sp_buffer->nSaddlePoints = 32;

		    for (j=0; j < sp_buffer->nSaddlePoints; j++, offs++)
			sp_buffer->sp[j] = pkd->saddle_points.sp[gd[i].sp[offs]];
		    mdlRelease(pkd->mdl,CID_SADDLE_BUF,sp_buffer);
#ifdef MPI_VERSION
		    mdlFlushCache(pkd->mdl, CID_SADDLE_BUF);
#else
		    assert(0);
#endif
		}
	    }
	}
    }
    mdlFinishCache(pkd->mdl,CID_SADDLE_BUF);


    mdlFree(pkd->mdl, pkd->saddle_points.buf);

#if 0
    for (i=1; i < pkd->psGroupTable.nGroups; i++)
    {
	if (gd[i].iPid == pkd->idSelf && !gd[i].dup)
	{
	    if (gd[i].bridge)
	    {
		fprintf(stderr, "%f %f %f\n", gd[i].fDensity, gd[i].sp.nbr_grp_fDensity, gd[i].sp.fDensity / gd[i].fDensity);
		gd[i].bridge = gd[i].fDensity < gd[i].sp.nbr_grp_fDensity;

		//gd[i].bridge = gd[i].bridge && (gd[i].sp.fDensity / gd[i].fDensity) > 0.90;
	    }
	}
    }
#endif

    //psdSetGlobalId(psx);

#if 0
    for (pi=0;pi<pkd->nLocal;++pi) 
    {
	*pkdGroup(pkd, pkdParticle(pkd,pi)) = 0;
    }
#endif

#if 0
    while (!STACK_EMPTY(G))
    {
	*pkdGroup(pkd, pkdParticle(pkd,POP(G))) = 1;
    }
#endif

#if 0
    for (i=1; i < pkd->psGroupTable.nGroups; i++)
    {
	if (gd[i].sp.iPid == pkd->idSelf)
	{
	    *pkdGroup(pkd, pkdParticle(pkd, gd[i].sp.iLocalId)) = 500;
	}
    }
#endif

    for (pi=0;pi<pkd->nLocal;++pi) 
    {
	PARTICLE *p = pkdParticle(pkd,pi);
	float *a = pkdAccel(pkd, p);
	a[0] = a[1] = a[2] = 0;
    }

#if 1
    mdlROcache(pkd->mdl,CID_GROUP,NULL,pkd->psGroupTable.pGroup,sizeof(struct psGroup), pkd->psGroupTable.nGroups);
    mdlCOcache(pkd->mdl,CID_PARTICLE,NULL,pkdParticleBase(pkd),pkdParticleSize(pkd), pkd->nLocal, pkd, initMarkSaddlePoint, combMarkSaddlePoint);

    for (i=1; i < pkd->psGroupTable.nGroups; i++)
    {
	if (gd[i].iPid == pkd->idSelf && !gd[i].dup)
	{

#if 0
	    if (!gd[i].bridge)
	    {
		fprintf(stderr, "%i] %i not in a bridge\n", pkd->idSelf, gd[i].iGlobalId);
	    }
#endif

	    PARTICLE *p;
	    struct psGroup *gdentry;
	    //if (!gd[i].bridge) continue;

	    for (j=0; j < gd[i].nSaddlePoints; j++)
	    {
		struct saddle_point *sp = pkd->saddle_points.sp + gd[i].sp[j];

		int remote_p   = sp->iPid != pkd->idSelf;
		int remote_nbr = sp->nbr.iPid != pkd->idSelf;

		if (remote_p)
		    p = mdlAquire(pkd->mdl,CID_PARTICLE,sp->iLocalId,sp->iPid);
		else
		    p = pkdParticle(pkd, sp->iLocalId);

		assert(sp->nbr.iLocalId > 0);
		if (remote_nbr)
		    gdentry = mdlAquire(pkd->mdl,CID_GROUP, sp->nbr.iLocalId, sp->nbr.iPid);
		else
		    gdentry = gd + sp->nbr.iLocalId;

		float *a = pkdAccel(pkd, p);
		a[0] = gdentry->r[0] - p->r[0];
		a[1] = gdentry->r[1] - p->r[1];
		a[2] = gdentry->r[2] - p->r[2];
		//fprintf(stderr, "%f %f %f\n", a[0], a[1], a[2]);

		*pkdGroup(pkd, p) = 1;

		if (remote_p)   mdlRelease(pkd->mdl,CID_PARTICLE,p);
		if (remote_nbr) mdlRelease(pkd->mdl,CID_GROUP,gdentry);
		mdlCacheCheck(pkd->mdl);
	    }
	}
    }

    mdlFinishCache(pkd->mdl,CID_PARTICLE);
    mdlFinishCache(pkd->mdl,CID_GROUP);
#endif

    FREE_STACK(G);
#endif
}
