#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <math.h>
#include "groupstats.h"
#include "group.h"
#include "vqsort.h"


#if 0
double illinois(double (*func)(double),double r,double s,double xacc,int *pnIter) {
    const int maxIter = 100;
    double t,fr,fs,ft,phis,phir,gamma;
    int i;

    fr = func(r);
    fs = func(s);
    assert(fr*fs < 0);
    t = (s*fr - r*fs)/(fr - fs);
    for (i=0;i < maxIter && fabs(t-s) > xacc;++i) {

	ft = func(t);

	if (ft*fs < 0) {
	    /*
	    ** Unmodified step.
	    */
	    r = s;
	    s = t;
	    fr = fs;
	    fs = ft;
	}
	else {
	    /*
	    ** Modified step to make sure we do not retain the 
	    ** endpoint r indefinitely.
	    */
#if 1
	    phis = ft/fs;
	    phir = ft/fr;
	    gamma = 1 - (phis/(1-phir));  /* method 3 */
	    if (gamma < 0) gamma = 0.5;
#else
	    gamma = 0.5;    /* illinois */
#endif
	    fr *= gamma;
	    s = t;
	    fs = ft;
	}
	t = (s*fr - r*fs)/(fr - fs);
/*
	printf("%d %d r:%.14f fr:%.14f s:%.14f fs:%.14f t:%.14f\n",
	       i,r,fr,s,fs,t);
*/
    }
    *pnIter = i;
    return(t);
}
#endif



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
    g->sigma = 0.0;
    g->rMax =  0.0;
    }

static void combTinyGroup(void *vpkd, void *v1, void *v2) {
    TinyGroupTable *g1 = (TinyGroupTable *)v1;
    TinyGroupTable *g2 = (TinyGroupTable *)v2;
    float x;
    int j;

    g1->fMass += g2->fMass;
    x = g2->fMass/g1->fMass;
    for (j=0;j<3;j++) {
	g1->rcom[j] = (1-x)*g1->rcom[j] + x*g2->rcom[j];
	g1->vcom[j] = (1-x)*g1->vcom[j] + x*g2->vcom[j];
	}
    g1->sigma = (1-x)*g1->sigma + x*g2->sigma;
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

typedef struct {
    float fMass;
    float dr;
    }  MassRadius;

#define mrLessThan(a,b) (((MassRadius *)a)->dr < ((MassRadius *)b)->dr)

static void initRoot(void *vpkd, void *v) {
    MassRadius * g = (MassRadius *)v;

    g->fMass = 0.0;
    }

static void combRoot(void *vpkd, void *v1, void *v2) {
    MassRadius *g1 = (MassRadius *)v1;
    MassRadius *g2 = (MassRadius *)v2;

    g1->fMass += g2->fMass;
    }

static KDN *getCell(PKD pkd, int iCell, int id) {
    if (id==pkd->idSelf) return pkdTreeNode(pkd,iCell);
    return mdlFetch(pkd->mdl,CID_CELL,iCell,id);
    }

static double gatherMass(PKD pkd,remoteID *S,double fBall2,double ri2,double ro2,double r[3]) {
    KDN *kdn;
    MDL mdl = pkd->mdl;
    double min2,max2;
    int iCell,id;
    int sp = 0;
    const BND *bnd;
    PARTICLE *p;
    double p_r[3];
    double dx,dy,dz,fDist2;
    int pj,pEnd;
    double fMass = 0;

    kdn = getCell(pkd,iCell=pkd->iTopTree[ROOT],id = pkd->idSelf);
    while (1) {
        bnd = pkdNodeBnd(pkd,kdn);
	MINDIST(bnd,r,min2);
	if (min2 > ri2) goto NoIntersect;
	MAXDIST(bnd,r,max2);
	if (max2 <= ro2) {
	    fMass += pkdNodeMom(pkd,kdn)->m;
	    goto NoIntersect;
	    }
	/*
	** We have an intersection to test.
	*/
	if (kdn->iLower) {
	    int idUpper,iUpper;
	    pkdGetChildCells(kdn,id,id,iCell,idUpper,iUpper);
	    kdn = getCell(pkd,iCell,id);
	    S[sp].iPid = idUpper;
	    S[sp].iIndex = iUpper;
	    ++sp;
	    continue;
	    }
	else {
	    if (id == pkd->idSelf) {
		pEnd = kdn->pUpper;
		for (pj=kdn->pLower;pj<=pEnd;++pj) {
		    p = pkdParticle(pkd,pj);
		    pkdGetPos1(pkd,p,p_r);
		    dx = r[0] - p_r[0];
		    dy = r[1] - p_r[1];
		    dz = r[2] - p_r[2];
		    fDist2 = dx*dx + dy*dy + dz*dz;
		    if (fDist2 <= fBall2) {
			fMass += pkdMass(pkd,p);
			}
		    }
		}
	    else {
		pEnd = kdn->pUpper;
		for (pj=kdn->pLower;pj<=pEnd;++pj) {
		    p = mdlFetch(mdl,CID_PARTICLE,pj,id);
		    pkdGetPos1(pkd,p,p_r);
		    dx = r[0] - p_r[0];
		    dy = r[1] - p_r[1];
		    dz = r[2] - p_r[2];
		    fDist2 = dx*dx + dy*dy + dz*dz;
		    if (fDist2 <= fBall2) {
			fMass += pkdMass(pkd,p);
			}
		    }
		}
	    }
    NoIntersect:
	if (sp) {
	    --sp;
	    id = S[sp].iPid;
	    iCell = S[sp].iIndex;
	    kdn = getCell(pkd,iCell,id);
	    }
	else break;
	}
    return(fMass);
    }


static double gatherLocalMass(PKD pkd,remoteID *S,double fBall2,double ri2,double ro2,double r[3]) {
    KDN *kdn;
    MDL mdl = pkd->mdl;
    double min2,max2;
    int iCell;
    int sp = 0;
    const BND *bnd;
    PARTICLE *p;
    double p_r[3];
    double dx,dy,dz,fDist2;
    int pj,pEnd;
    double fMass = 0;

    kdn = pkdTreeNode(pkd,iCell=ROOT);
    while (1) {
        bnd = pkdNodeBnd(pkd,kdn);
	MINDIST(bnd,r,min2);
	if (min2 > ri2) goto NoIntersect;
	MAXDIST(bnd,r,max2);
	if (max2 <= ro2) {
	    fMass += pkdNodeMom(pkd,kdn)->m;
	    goto NoIntersect;
	    }
	/*
	** We have an intersection to test.
	*/
	if (kdn->iLower) {
	    kdn = pkdTreeNode(pkd,iCell = kdn->iLower);
	    S[sp++].iIndex = iCell+1;
	    continue;
	    }
	else {
	    pEnd = kdn->pUpper;
	    for (pj=kdn->pLower;pj<=pEnd;++pj) {
		p = pkdParticle(pkd,pj);
		pkdGetPos1(pkd,p,p_r);
		dx = r[0] - p_r[0];
		dy = r[1] - p_r[1];
		dz = r[2] - p_r[2];
		fDist2 = dx*dx + dy*dy + dz*dz;
		if (fDist2 <= fBall2) {
		    fMass += pkdMass(pkd,p);
		    }
		}
	    }
    NoIntersect:
	if (sp) kdn = pkdTreeNode(pkd,iCell = S[--sp].iIndex);
	else break;
	}
    return(fMass);
    }


double pkdGatherMass(PKD pkd,remoteID *S,double fBall,double r[3],int bPeriodic,double dPeriod[3]) {
    double rp[3];
    int ix,iy,iz,nReps;
    double fMass=0;
    double fBall2 = fBall*fBall;
    double ri,ri2,ro2;

    assert(pkd->oNodeMom);    
    ri = 0.9*fBall;        /* we use an approximate radius but make sure it is unbiased by volume */
    ri2 = ri*ri;
    ro2 = pow(2*fBall2*fBall - ri2*ri,2.0/3.0);
    if (bPeriodic) nReps = 1;
    else nReps = 0;
    for (ix=-nReps;ix<=nReps;++ix) {
	rp[0] = r[0] + ix*dPeriod[0];
	for (iy=-nReps;iy<=nReps;++iy) {
	    rp[1] = r[1] + iy*dPeriod[1];
	    for (iz=-nReps;iz<=nReps;++iz) {
		rp[2] = r[2] + iz*dPeriod[2];
		fMass += gatherMass(pkd,S,fBall2,ri2,ro2,rp);
		}
	    }
	}
    return(fMass);
    }


typedef struct {
    float a;
    float b;
    int bConverged;
    int iter;            /* this is only used for a diagnostic */
    int gid;
    } RootFindingTable;


void pkdCalculateGroupStats(PKD pkd,int bPeriodic,double *dPeriod,double rEnvironment[2]) {
    MDL mdl = pkd->mdl;
    PARTICLE *p;
    float fPot;
    double fMass;
    double r[3],v2;
    double dHalf[3];
    double diVol0,diVol1;
    float r2,rMax;
    int i,j,gid,n;
    vel_t *v;
    TinyGroupTable *g;
    int nLocalGroups;
    double *dAccumulate;
    MassRadius *mr,*mrFree,*rootFunction,*pmr;
    remoteID *S;
    uint32_t *iGrpOffset,*iGrpEnd;
    int nRootFind,bIncomplete,nMaxIter,iter;
    int *mrIndex,iRoot,*bRemoteDone;
    RootFindingTable *rootFindingTable;

    assert(pkd->nGroups*(sizeof(*pkd->ga)+sizeof(*pkd->tinyGroupTable)+4*sizeof(double)) < EPHEMERAL_BYTES*pkd->nStore);
    pkd->tinyGroupTable = (TinyGroupTable *)(&pkd->ga[pkd->nGroups]);
    dAccumulate = (double *)(&pkd->tinyGroupTable[pkd->nGroups]);
    /*
    ** Initialize the table.
    */
    nLocalGroups = 0;
    for (gid=0;gid<pkd->nGroups;++gid) {
	if (pkd->ga[gid].id.iPid == pkd->idSelf && gid) ++nLocalGroups;
	pkd->tinyGroupTable[gid].minPot = FLOAT_MAXVAL;
	pkd->tinyGroupTable[gid].rMax = 0;
	dAccumulate[4*gid] = 0;
	dAccumulate[4*gid+1] = 0;
	dAccumulate[4*gid+2] = 0;
	dAccumulate[4*gid+3] = 0;
	}
    /*
    ** First determine the minimum potential particle for each group.
    ** This will be the reference position for the group as well.
    ** In this first pass we also sum up the v2 and vcom locally.
    */
    for (i=0;i<pkd->nLocal;++i) {
	p = pkdParticle(pkd,i);
	gid = pkdGetGroup(pkd,p);
	fMass = pkdMass(pkd,p);
	if (!gid) continue;
	fPot = *pkdPot(pkd,p);
	v = pkdVel(pkd,p);
	v2 = 0;
	for (j=0;j<3;++j) {
	    dAccumulate[4*gid+j] += fMass*v[j];
	    v2 += v[j]*v[j];
	    }
	dAccumulate[4*gid+3] += fMass*v2;
	if (fPot < pkd->tinyGroupTable[gid].minPot) {
	    pkd->tinyGroupTable[gid].minPot = fPot;
	    for (j=0;j<3;++j) {
		pkd->tinyGroupTable[gid].rPot[j] = pkdPosRaw(pkd,p,j);
		}
	    }
	}
    for (gid=1;gid<pkd->nGroups;++gid) {
	for (j=0;j<3;++j) {
	    pkd->tinyGroupTable[gid].vcom[j] = dAccumulate[4*gid+j];
	    dAccumulate[4*gid+j] = 0;
	    }
	pkd->tinyGroupTable[gid].sigma = dAccumulate[4*gid+3];
	dAccumulate[4*gid+3] = 0;
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
	g = mdlAcquire(mdl,CID_GROUP,pkd->ga[gid].id.iIndex,pkd->ga[gid].id.iPid);
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
	mdlRelease(mdl,CID_GROUP,g);
	}
    mdlFinishCache(mdl,CID_GROUP);
    /*
    ** Now find rcom and rMax and total mass locally.
    */
    for (j=0;j<3;++j) dHalf[j] = bPeriodic ? 0.5 * dPeriod[j] : FLOAT_MAXVAL;
    for (i=0;i<pkd->nLocal;++i) {
	p = pkdParticle(pkd,i);
	gid = pkdGetGroup(pkd,p);
	fMass = pkdMass(pkd,p);
	dAccumulate[4*gid+3] += fMass;
	if (gid > 0) {
	    r2 = 0.0;
	    for (j=0;j<3;++j) {
		r[j] =  pkdPosToDbl(pkd,pkdPosRaw(pkd,p,j) - pkd->tinyGroupTable[gid].rPot[j]);
		if      (r[j] < -dHalf[j]) r[j] += dPeriod[j];
		else if (r[j] > +dHalf[j]) r[j] -= dPeriod[j];
		dAccumulate[4*gid+j] += fMass*r[j];
		r2 += r[j]*r[j];
		}
	    rMax = sqrtf(r2);
	    if (rMax > pkd->tinyGroupTable[gid].rMax) pkd->tinyGroupTable[gid].rMax = rMax;
	    }
	}
    for (gid=1;gid<pkd->nGroups;++gid) {
	for (j=0;j<3;++j) {
	    pkd->tinyGroupTable[gid].rcom[j] = dAccumulate[4*gid+j]/dAccumulate[4*gid+3];
	    pkd->tinyGroupTable[gid].vcom[j] /= dAccumulate[4*gid+3];
	    dAccumulate[4*gid+j] = 0;
	    }
	pkd->tinyGroupTable[gid].sigma /= dAccumulate[4*gid+3];  /* "sigma" is actually still v^2 */
	pkd->tinyGroupTable[gid].fMass = dAccumulate[4*gid+3];
	dAccumulate[4*gid+3] = 0;  /* if we still need masses later we should not clear this field */
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
	    g->sigma += pkd->tinyGroupTable[gid].sigma;
	    }
	if (pkd->tinyGroupTable[gid].rMax > g->rMax) g->rMax = pkd->tinyGroupTable[gid].rMax;
	}
    mdlFinishCache(mdl,CID_GROUP);
    for(gid=1;gid<=nLocalGroups;++gid) {
	v2 = 0;
	for (j=0;j<3;j++) {
            /*
	    ** If we want absolute rcom instead of relative to the minimum potential then we
	    ** need to uncomment the line below.
	    ** pkd->tinyGroupTable[gid].rcom[j] += pkdPosToDbl(pkd,pkd->tinyGroupTable[gid].rPot[j]);
	    */
	    v2 += pkd->tinyGroupTable[gid].vcom[j]*pkd->tinyGroupTable[gid].vcom[j];
	    }
	/*
	** Now convert from mass weighted average v2 to sigma.
	*/
	pkd->tinyGroupTable[gid].sigma -= v2;
	pkd->tinyGroupTable[gid].sigma = sqrtf(pkd->tinyGroupTable[gid].sigma);
	}
    /*
    ** Scoop environment mass (note the full tree must be present).
    */
    S = (remoteID *)(&pkd->tinyGroupTable[pkd->nGroups]);
    mdlROcache(mdl,CID_PARTICLE,NULL,pkdParticleBase(pkd),pkdParticleSize(pkd),pkdLocal(pkd));
    diVol0 = 3.0/4.0*M_1_PI*pow(rEnvironment[0],-3);
    diVol1 = 3.0/4.0*M_1_PI*pow(rEnvironment[1],-3);
    for (gid=1;gid<=nLocalGroups;++gid) {
	for (j=0;j<3;++j) r[j] = pkdPosToDbl(pkd,pkd->tinyGroupTable[gid].rPot[j]);
	if (rEnvironment[0] > 0.0) {
	    pkd->tinyGroupTable[gid].fEnvironDensity0 = 
		pkdGatherMass(pkd,S,rEnvironment[0],r,bPeriodic,dPeriod)*diVol0;
	    }
	else pkd->tinyGroupTable[gid].fEnvironDensity0 = 0.0;
	if (rEnvironment[1] > 0.0) {
	    pkd->tinyGroupTable[gid].fEnvironDensity1 = 
		pkdGatherMass(pkd,S,rEnvironment[1],r,bPeriodic,dPeriod)*diVol1;
	    }
	else pkd->tinyGroupTable[gid].fEnvironDensity1 = 0.0; 
	}
    mdlFinishCache(mdl,CID_PARTICLE);
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
	pkd->tinyGroupTable[gid].sigma = g->sigma;
	pkd->tinyGroupTable[gid].rMax = g->rMax;
	pkd->tinyGroupTable[gid].fEnvironDensity0 = g->fEnvironDensity0;
	pkd->tinyGroupTable[gid].fEnvironDensity1 = g->fEnvironDensity1;	
	}
    mdlFinishCache(mdl,CID_GROUP);
    /*
    ** Now find the half mass radii of the groups.
    ** Start by reordering particles into group order.
    ** Note that we must not access tree cells until the next tree build! 
    */
    mdlFinishCache(mdl,CID_CELL);
    iGrpOffset = (uint32_t *)(&pkd->tinyGroupTable[pkd->nGroups]);
    mrIndex = (int *)(&iGrpOffset[2*pkd->nGroups]);    /* this stores the ending index of each group in the mr below */
    mr = (MassRadius *)(&mrIndex[pkd->nGroups]);  /* we must be careful how much of this we use! */
    pkdGroupOrder(pkd,iGrpOffset);
    iGrpEnd = &iGrpOffset[pkd->nGroups];
    iGrpEnd[0] = 0;
    /*
    ** First do all purely local groups. We can do these one at a time without 
    ** worrying about remote particles. We should use the bounds determined by 
    ** fof previously.
    */
    mrFree = mr;
    nRootFind = 0;
    mrIndex[0] = 0;
    for (gid=1;gid<=nLocalGroups;++gid) {
	rMax = pkd->tinyGroupTable[gid].rMax;
	for (j=0;j<3;++j) {
	    double rt = pkd->bndInterior.fMax[j] - 
		fabs(pkd->bndInterior.fCenter[j] - pkdPosToDbl(pkd,pkd->tinyGroupTable[gid].rPot[j]));
	    if (rt < rMax) rMax = rt;
	    }
	n = 0;
	for (i=iGrpEnd[gid-1];i<iGrpEnd[gid];++i,++n) {
	    p = pkdParticle(pkd,i);
	    assert(gid == pkdGetGroup(pkd,p)); /* make sure it is really in group order */
	    r2 = 0.0;
	    for (j=0;j<3;++j) {
		r[j] =  pkdPosToDbl(pkd,pkdPosRaw(pkd,p,j) - pkd->tinyGroupTable[gid].rPot[j]);
		if      (r[j] < -dHalf[j]) r[j] += dPeriod[j];
		else if (r[j] > +dHalf[j]) r[j] -= dPeriod[j];
		r2 += r[j]*r[j];
		}
	    mrFree[n].fMass = pkdMass(pkd,p);
	    mrFree[n].dr = sqrtf(r2);
	    }
	QSORT(sizeof(MassRadius),mrFree,n,mrLessThan);
	fMass = mrFree[i].fMass;
	for (i=1;i<n;++i) {
	    fMass += mrFree[i].fMass;
	    mrFree[i].fMass = fMass;
	    }
	/*
	** If the radius in this sum ever exceeds rMax, then we have to resolve this
	** group's rHalf by using particles from the other remote processors.
	*/
	bIncomplete = 0;
	for (i=0;i<n;++i) {
	    if (mrFree[i].dr > rMax) {
		bIncomplete = 1;
		break;
		}
	    else if (mrFree[i].fMass > 0.5*pkd->tinyGroupTable[gid].fMass) break;
	    }
	pkd->tinyGroupTable[gid].rHalf = (i > 0)?mrFree[--i].dr:0;
	if (bIncomplete) {
	    /*
	    ** Compact the mrFree array.
	    */
	    for (j=0;i<n;++i,++j) mrFree[j] = mrFree[i];
	    /*
	    ** We need to do some root finding on this group in order to determine the 
	    ** correct rHalf. The stored rHalf is only a lower bound at present in this
	    ** case. Only the rootTable is exposed to the remote processors via the cache.
	    */
	    mrIndex[gid] = mrFree - mr + j;
	    /*
	    ** Set the new mrFree to unallocated space.
	    */
	    mrFree = mr + mrIndex[gid];
	    ++nRootFind;
	    }
	else {
	    mrIndex[gid] = mrFree - mr;
	    }
	}
    for (gid=nLocalGroups+1;gid<pkd->nGroups;++gid) {
	n = 0;
	for (i=iGrpEnd[gid-1];i<iGrpEnd[gid];++i,++n) {
	    p = pkdParticle(pkd,i);
	    assert(gid == pkdGetGroup(pkd,p)); /* make sure it is really in group order */
	    r2 = 0.0;
	    for (j=0;j<3;++j) {
		r[j] =  pkdPosToDbl(pkd,pkdPosRaw(pkd,p,j) - pkd->tinyGroupTable[gid].rPot[j]);
		if      (r[j] < -dHalf[j]) r[j] += dPeriod[j];
		else if (r[j] > +dHalf[j]) r[j] -= dPeriod[j];
		r2 += r[j]*r[j];
		}
	    mrFree[n].fMass = pkdMass(pkd,p);
	    mrFree[n].dr = sqrtf(r2);
	    }
	QSORT(sizeof(MassRadius),mrFree,n,mrLessThan);
	fMass = mrFree[i].fMass;
	for (i=1;i<n;++i) {
	    fMass += mrFree[i].fMass;
	    mrFree[i].fMass = fMass;
	    }
	mrIndex[gid] = mrFree - mr + n;
	/*
	** Set the new mrFree to unallocated space.
	*/
	mrFree = mr + mrIndex[gid];
	}
    /* 
    ** Set up a table which defines the rootFunction of all enclosed 
    ** group mass at a given dr. Also set up extra table needed to keep 
    ** track of lower and upper bounds for root finding after all the 
    ** saved mass-radius values.
    ** We will only need nRootFind entries in this table.
    */
    iRoot = 0;
    rootFunction = mrFree;
    bRemoteDone = (int *)(&rootFunction[nLocalGroups+1]);
    for (gid=nLocalGroups+1;gid<pkd->nGroups;++gid) bRemoteDone[gid-(nLocalGroups+1)] = 0;
    rootFindingTable = (RootFindingTable *)(&bRemoteDone[pkd->nGroups-(nLocalGroups+1)]);
    for (gid=0;gid<=nLocalGroups;++gid) {
	if (mrIndex[gid] > mrIndex[gid-1]) {
	    rootFindingTable[iRoot].bConverged = 0;
	    rootFindingTable[iRoot].gid = gid;
	    rootFindingTable[iRoot].iter = 65536;
	    rootFindingTable[iRoot].a = pkd->tinyGroupTable[gid].rHalf;
	    rootFindingTable[iRoot].b = pkd->tinyGroupTable[gid].rMax;	    
	    ++iRoot;
	    }
	else {
	    rootFunction[gid].dr = -1.0; /* marks no root finding needed */
	    }
	rootFunction[gid].fMass = 0;
	}
    assert(iRoot == nRootFind);
    /*
    ** Now everything is set up for the root finding phases.
    */
    nMaxIter = 16;
    for (iter=0;iter<nMaxIter;++iter) {
	/*
	** First calculate the local mass contained in a trial dr.
	*/
	for (iRoot=0;iRoot<nRootFind;++iRoot) {
	    gid = rootFindingTable[iRoot].gid;
	    if (!rootFindingTable[iRoot].bConverged) {
		rootFunction[gid].dr = 0.5*(rootFindingTable[iRoot].a + rootFindingTable[iRoot].b);
		/*
		** Binary search for rootFunction[gid].dr
		*/
		int ia = mrIndex[gid-1];
		int ib = mrIndex[gid] - 1;
		while (ia < ib) {
		    int ic = (ia+ib+1)/2;
		    if (mr[ic].dr > rootFunction[gid].dr) ib = ic-1;
		    else ia = ic;
		    }
		assert(ia == ib);
		rootFunction[gid].fMass = (rootFunction[gid].dr > mr[ia].dr)?mr[ia].fMass:0;
		}
	    else {
		rootFunction[gid].dr = -1.0;  /* communicates to the remote that nothing needs to be done for this gid */
		}
	    }
	/*
	** Now scatter the enclosed local mass to remote groups at their trial dr.
	*/
	mdlCOcache(mdl,CID_RM,NULL,rootFunction,sizeof(MassRadius),nLocalGroups+1,
	    NULL,initRoot,combRoot);
  	for (gid=nLocalGroups+1;gid<pkd->nGroups;++gid) {
	    if (bRemoteDone[gid-(nLocalGroups+1)]) continue;
	    pmr = mdlAcquire(mdl,CID_RM,pkd->ga[gid].id.iIndex,pkd->ga[gid].id.iPid);
	    if (pmr->dr > 0) {
		/*
		** Binary search for pmr->dr!
		*/
		int ia = mrIndex[gid-1];
		int ib = mrIndex[gid] - 1;
		while (ia < ib) {
		    int ic = (ia+ib+1)/2;
		    if (mr[ic].dr > pmr->dr) ib = ic-1;
		    else ia = ic;
		    }
		assert(ia == ib);
		pmr->fMass += (pmr->dr > mr[ia].dr)?mr[ia].fMass:0;
		}
	    else {
		bRemoteDone[gid-(nLocalGroups+1)] = 1;  /* This will even bypass the mdlAcquire in the next iteration. */
		}
	    mdlRelease(mdl,CID_RM,pmr);
	    }
	mdlFinishCache(mdl,CID_RM);
	/*
	** Now the root functions are complete.
	*/
	for (iRoot=0;iRoot<nRootFind;++iRoot) {	    
	    if (!rootFindingTable[iRoot].bConverged) {
		gid = rootFindingTable[iRoot].gid;
		if (rootFunction[gid].fMass > (0.5 + 0.6/pkd->ga[gid].nTotal)*pkd->tinyGroupTable[gid].fMass) 
		    rootFindingTable[iRoot].b = rootFunction[gid].dr;
		else if (rootFunction[gid].fMass < (0.5 - 0.6/pkd->ga[gid].nTotal)*pkd->tinyGroupTable[gid].fMass)
		    rootFindingTable[iRoot].a = rootFunction[gid].dr;
		else {
		    rootFindingTable[iRoot].bConverged = 1;
		    pkd->tinyGroupTable[gid].rHalf = rootFunction[gid].dr;
		    rootFindingTable[iRoot].iter = iter;
		    }
		}
	    }
	}
/*
    for (iRoot=0;iRoot<nRootFind;++iRoot) {
	gid = rootFindingTable[iRoot].gid;
	printf("%d: iter:%d gid:%d MassRatio:%g rHalf:%g\n",pkd->idSelf,rootFindingTable[iRoot].iter,gid,
	    rootFunction[gid].fMass/pkd->tinyGroupTable[gid].fMass,
	    pkd->tinyGroupTable[gid].rHalf);
	}
*/
    /*
    ** Important to really make sure pkd->nLocalGroups is set correctly before output!
    */
    pkd->nLocalGroups = nLocalGroups;
    }
