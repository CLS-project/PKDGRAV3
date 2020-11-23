/*  This file is part of PKDGRAV3 (http://www.pkdgrav.org/).
 *  Copyright (c) 2001-2018 Joachim Stadel & Douglas Potter
 *
 *  PKDGRAV3 is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  PKDGRAV3 is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with PKDGRAV3.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#else
#include "pkd_config.h"
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
    TinyGroupTable *g = (TinyGroupTable *)v;
    
    g->minPot = HUGE_VALF;  /* make sure it is always set */
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
    g->nBH = 0;
    g->nDM = 0;
    g->nGas = 0;
    g->nStar = 0;
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
    g1->nBH += g2->nBH;
    g1->nDM += g2->nDM;
    g1->nGas += g2->nGas;
    g1->nStar += g2->nStar;
    //printf("g1 %d \t g2 %d \n", g1->nBH, g2->nBH);
    }

static void initTinyRmax(void *vpkd, void *v) {
    TinyGroupTable * g = (TinyGroupTable *)v;
    int j;

    g->rMax =  0.0;
    }

static void combTinyRmax(void *vpkd, void *v1, void *v2) {
    TinyGroupTable *g1 = (TinyGroupTable *)v1;
    TinyGroupTable *g2 = (TinyGroupTable *)v2;

    if (g2->rMax > g1->rMax) g1->rMax = g2->rMax;
    }

typedef struct {
    double fMass;
    double rcom[3];
    int nEnclosed;
    } ShrinkStruct;

static void initShrink(void *vpkd, void *v) {
    ShrinkStruct * g = (ShrinkStruct *)v;
    int j;

    g->fMass = 0.0;
    for (j=0;j<3;j++) {
	g->rcom[j] = 0.0;
	}
    g->nEnclosed = 0;
    }

static void combShrink(void *vpkd, void *v1, void *v2) {
    ShrinkStruct *g1 = (ShrinkStruct *)v1;
    ShrinkStruct *g2 = (ShrinkStruct *)v2;
    float x;
    int j;

    g1->fMass += g2->fMass;
    x = g2->fMass/g1->fMass;
    for (j=0;j<3;j++) {
	g1->rcom[j] = (1-x)*g1->rcom[j] + x*g2->rcom[j];
	}
    g1->nEnclosed += g2->nEnclosed;
    }


static void initAngular(void *vpkd, void *v) {
    TinyGroupTable * g = (TinyGroupTable *)v;
    int j;

    for (j=0;j<3;j++) {
	g->angular[j] = 0.0;
	}
    for (j=0;j<6;j++) {
	g->inertia[j] = 0.0;
	}
    }

static void combAngular(void *vpkd, void *v1, void *v2) {
    TinyGroupTable *g1 = (TinyGroupTable *)v1;
    TinyGroupTable *g2 = (TinyGroupTable *)v2;
    float x;
    int j;

    for (j=0;j<3;j++) {
	g1->angular[j] += g2->angular[j];
	}
    for (j=0;j<6;j++) {
	g1->inertia[j] += g2->inertia[j];
	}
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
    BND bnd;
    PARTICLE *p;
    double p_r[3];
    double dx,dy,dz,fDist2;
    int pj,pEnd;
    double fMass = 0;

    kdn = getCell(pkd,iCell=pkd->iTopTree[ROOT],id = pkd->idSelf);
    while (1) {
        bnd = pkdNodeGetBnd(pkd,kdn);
	MINDIST(&bnd,r,min2);
	if (min2 > ri2) goto NoIntersect;
	MAXDIST(&bnd,r,max2);
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
    BND bnd;
    PARTICLE *p;
    double p_r[3];
    double dx,dy,dz,fDist2;
    int pj,pEnd;
    double fMass = 0;

    kdn = pkdTreeNode(pkd,iCell=ROOT);
    while (1) {
        bnd = pkdNodeGetBnd(pkd,kdn);
	MINDIST(&bnd,r,min2);
	if (min2 > ri2) goto NoIntersect;
	MAXDIST(&bnd,r,max2);
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
    const int bDoShrinkingSphere = 0;
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
    double vrel[3];
    TinyGroupTable *g;
    int nLocalGroups;
    double *dAccumulate;
    MassRadius *mr,*mrFree,*rootFunction,*pmr;
    remoteID *S;
    uint32_t *iGrpOffset,*iGrpEnd;
    int nRootFind,bIncomplete,nMaxIter,iter;
    int *mrIndex,iRoot,*bRemoteDone;
    RootFindingTable *rootFindingTable;
    ShrinkStruct *shrink;
    int bShrink;
    double f2;

    assert(pkd->nGroups*(sizeof(*pkd->ga)+sizeof(*pkd->tinyGroupTable)+sizeof(ShrinkStruct)) < 1ul*pkd->nEphemeralBytes*pkd->nStore);
    pkd->tinyGroupTable = (TinyGroupTable *)(&pkd->ga[pkd->nGroups]);
    dAccumulate = (double *)(&pkd->tinyGroupTable[pkd->nGroups]);
    /*
    ** Initialize the table.
    */
    nLocalGroups = 0;
    for (gid=0;gid<pkd->nGroups;++gid) {
	if (pkd->ga[gid].id.iPid == pkd->idSelf && gid) ++nLocalGroups;
	/*
	** Here we assume that the ga table was setup prior to calling
	** gravity so that we have filled in the minimum potential particle
	** for each group.
	*/
	pkd->tinyGroupTable[gid].minPot = pkd->ga[gid].minPot;
	if (gid) {
	    assert(pkd->ga[gid].minPot < FLOAT_MAXVAL);
	    assert(pkd->ga[gid].iMinPart < pkd->nLocal);
	    p = pkdParticle(pkd,pkd->ga[gid].iMinPart);
	    for (j=0;j<3;++j) {
		pkd->tinyGroupTable[gid].rPot[j] = pkdPos(pkd,p,j);
		pkd->tinyGroupTable[gid].rcen[j] = 0.0; /* relative to rPot */
		}
	    }
	pkd->tinyGroupTable[gid].rMax = 0;
      pkd->tinyGroupTable[gid].nBH = 0;
      pkd->tinyGroupTable[gid].nDM = 0;
      pkd->tinyGroupTable[gid].nGas = 0;
      pkd->tinyGroupTable[gid].nStar = 0;
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
	v = pkdVel(pkd,p);
	v2 = 0;
	for (j=0;j<3;++j) {
	    dAccumulate[4*gid+j] += fMass*v[j];
	    v2 += v[j]*v[j];
	    }
	dAccumulate[4*gid+3] += fMass*v2;
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
	g = mdlVirtualFetch(mdl,CID_GROUP,pkd->ga[gid].id.iIndex,pkd->ga[gid].id.iPid);
	/*
	** We update the remote center, the combiner makes sure the minimum minPot is set.
	** We don't have to check the virtual fetch value in this case since we fetch each
	** remote group only once (there cannot be more than one reference) otherwise we
	** would have to check the fetched value.
	*/
	g->minPot = pkd->tinyGroupTable[gid].minPot;
	for (j=0;j<3;++j) {
	    g->rPot[j] = pkd->tinyGroupTable[gid].rPot[j];
	    }	   
	}
    mdlFinishCache(mdl,CID_GROUP);
    mdlROcache(mdl,CID_GROUP,NULL,pkd->tinyGroupTable,sizeof(TinyGroupTable),nLocalGroups+1);
    for(gid=1+nLocalGroups;gid<pkd->nGroups;++gid) {
	g = mdlFetch(mdl,CID_GROUP,pkd->ga[gid].id.iIndex,pkd->ga[gid].id.iPid);
	/*
	** We update the local version of the remote center, here we can just copy it.
	*/
	pkd->tinyGroupTable[gid].minPot = g->minPot;
	for (j=0;j<3;++j) {
	    pkd->tinyGroupTable[gid].rPot[j] = g->rPot[j];
	    }	   
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
      int pSpecies = pkdSpecies(pkd,p);
      switch (pSpecies){
         case FIO_SPECIES_BH:
            pkd->tinyGroupTable[gid].nBH++;
            break;
         case FIO_SPECIES_DARK:
            pkd->tinyGroupTable[gid].nDM++;
            break;
         case FIO_SPECIES_SPH:
            pkd->tinyGroupTable[gid].nGas++;
            break;
         case FIO_SPECIES_STAR:
            pkd->tinyGroupTable[gid].nStar++;
            break;
         default: // Maybe a deleted particle, skip it
            continue;
      }
	dAccumulate[4*gid+3] += fMass;
	if (gid > 0) {
	    r2 = 0.0;
	    for (j=0;j<3;++j) {
		r[j] =  pkdPos(pkd,p,j) - pkd->tinyGroupTable[gid].rPot[j];
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
      g->nBH += pkd->tinyGroupTable[gid].nBH;
      g->nDM += pkd->tinyGroupTable[gid].nDM;
      g->nGas += pkd->tinyGroupTable[gid].nGas;
      g->nStar += pkd->tinyGroupTable[gid].nStar;
      //printf("g->nBH %d \t pkd %d \n", g->nBH, pkd->tinyGroupTable[gid].nBH);
	for (j=0;j<3;j++) {
	    g->rcom[j] += pkd->tinyGroupTable[gid].rcom[j];
	    g->vcom[j] += pkd->tinyGroupTable[gid].vcom[j];
	    g->sigma += pkd->tinyGroupTable[gid].sigma;
	    }
	if (pkd->tinyGroupTable[gid].rMax > g->rMax) g->rMax = pkd->tinyGroupTable[gid].rMax;
	}
    mdlFinishCache(mdl,CID_GROUP);
    mdlROcache(mdl,CID_GROUP,NULL,pkd->tinyGroupTable,sizeof(TinyGroupTable),nLocalGroups+1);
    for(gid=1+nLocalGroups;gid<pkd->nGroups;++gid) {
	TinyGroupTable *g;
	g = mdlFetch(mdl,CID_GROUP,pkd->ga[gid].id.iIndex,pkd->ga[gid].id.iPid);
	pkd->tinyGroupTable[gid].fMass = g->fMass;
      pkd->tinyGroupTable[gid].nBH = g->nBH;
      pkd->tinyGroupTable[gid].nDM = g->nDM;
      pkd->tinyGroupTable[gid].nGas = g->nGas;
      pkd->tinyGroupTable[gid].nStar = g->nStar;
	for (j=0;j<3;j++) {
	    pkd->tinyGroupTable[gid].rcom[j] = g->rcom[j];
	    pkd->tinyGroupTable[gid].vcom[j] = g->vcom[j];
	    pkd->tinyGroupTable[gid].sigma = g->sigma;
	    }
	pkd->tinyGroupTable[gid].rMax = g->rMax;
	}
    mdlFinishCache(mdl,CID_GROUP);

    for(gid=1;gid<pkd->nGroups;++gid) {
	v2 = 0;
	for (j=0;j<3;j++) {
            /*
	    ** If we want absolute rcom instead of relative to the minimum potential then we
	    ** need to uncomment the line below.
	    ** pkd->tinyGroupTable[gid].rcom[j] += pkd->tinyGroupTable[gid].rPot[j];
	    */
	    v2 += pkd->tinyGroupTable[gid].vcom[j]*pkd->tinyGroupTable[gid].vcom[j];
	    }
	/*
	** Now convert from mass weighted average v2 to sigma.
	*/
	pkd->tinyGroupTable[gid].sigma -= v2;
	pkd->tinyGroupTable[gid].sigma = sqrtf(pkd->tinyGroupTable[gid].sigma);
	}
    if (bDoShrinkingSphere) {
	/*
	** Set the center of each group to be the center of mass, at least initially...
	*/
	for (gid=1;gid<pkd->nGroups;++gid) {
	    for (j=0;j<3;++j) {
		pkd->tinyGroupTable[gid].rcen[j] = pkd->tinyGroupTable[gid].rcom[j];
		}
	    pkd->tinyGroupTable[gid].rMax = 0; /* reset this here since we will recalculate it about the new center */
	    }
	/*
	** Now find rMax about the center (which is currently the center of mass).
	*/
	for (i=0;i<pkd->nLocal;++i) {
	    p = pkdParticle(pkd,i);
	    gid = pkdGetGroup(pkd,p);
	    if (gid > 0) {
		r2 = 0.0;
		for (j=0;j<3;++j) {
		    r[j] =  pkdPos(pkd,p,j) - pkd->tinyGroupTable[gid].rPot[j];
		    if      (r[j] < -dHalf[j]) r[j] += dPeriod[j];
		    else if (r[j] > +dHalf[j]) r[j] -= dPeriod[j];
		    /*
		    ** Make the position relative to the current best center.
		    ** This starts as being the center of mass of the entire group.
		    */
		    r[j] -= pkd->tinyGroupTable[gid].rcen[j];
		    r2 += r[j]*r[j];
		    }
		rMax = sqrtf(r2);
		if (rMax > pkd->tinyGroupTable[gid].rMax) pkd->tinyGroupTable[gid].rMax = rMax;
		}
	    }
	/* 
	** Now find the largest rmax for the local groups globally.
	*/
	mdlCOcache(mdl,CID_GROUP,NULL,pkd->tinyGroupTable,sizeof(TinyGroupTable),nLocalGroups+1,
	    NULL,initTinyRmax,combTinyRmax);
	for(gid=1+nLocalGroups;gid<pkd->nGroups;++gid) {
	    TinyGroupTable *g;
	    g = mdlVirtualFetch(mdl,CID_GROUP,pkd->ga[gid].id.iIndex,pkd->ga[gid].id.iPid);
	    g->rMax = pkd->tinyGroupTable[gid].rMax;
	    }
	mdlFinishCache(mdl,CID_GROUP);
	mdlROcache(mdl,CID_GROUP,NULL,pkd->tinyGroupTable,sizeof(TinyGroupTable),nLocalGroups+1);
	for(gid=1+nLocalGroups;gid<pkd->nGroups;++gid) {
	    TinyGroupTable *g;
	    g = mdlFetch(mdl,CID_GROUP,pkd->ga[gid].id.iIndex,pkd->ga[gid].id.iPid);
	    pkd->tinyGroupTable[gid].rMax = g->rMax;
	    }
	mdlFinishCache(mdl,CID_GROUP);

	shrink = (ShrinkStruct *)dAccumulate; /* overlay the shrink structure on the accumulators */
	for (gid=0;gid<pkd->nGroups;++gid) shrink[gid].nEnclosed = pkd->ga[gid].nTotal; /* all particles in the group */

	/*
	** We also may want the "shrinking sphere" center for many applications.
	** At this point we have a sphere of maximum radius enclosing all the particles
	** as well as an initial guess of the center. We now shrink the sphere by
	** 2.5% (a la Power, et. al. 2003) and find a new center of mass. We continue this
	** process until we reach a minumum number of particles within the sphere (~1000).
	*/
	f2 = 1.0;
	while (1) {
	    /*
	    ** Start of the main shrinking iteration.
	    */
	    bShrink = 0;  /* by default we will not shrink anything */
	    /*
	    ** Clear the sphere totals
	    */
	    for (gid=0;gid<pkd->nGroups;++gid) { 
		for (j=0;j<3;++j) shrink[gid].rcom[j] = 0.0;
		shrink[gid].fMass = 0.0;
		if (shrink[gid].nEnclosed > 1000) {
		    shrink[gid].nEnclosed = 0;
		    bShrink = 1;
		    }
		else {
		    shrink[gid].nEnclosed = -1;  /* marks not to shrink this group */
		    }
		}
	    shrink[0].nEnclosed = -1;  /* make sure ungrouped particles (gid==0) are not involved */
	    if (!bShrink) break; /* exits the main shrinking iteration...we are done */
	    f2 *= 0.975*0.975; /* make the shrink factor 2.5 % smaller */
	    for (i=0;i<pkd->nLocal;++i) {
		p = pkdParticle(pkd,i);
		gid = pkdGetGroup(pkd,p);
		if (shrink[gid].nEnclosed < 0) continue; /* skip this particle */
		r2 = 0.0;
		for (j=0;j<3;++j) {
		    r[j] =  pkdPos(pkd,p,j) - pkd->tinyGroupTable[gid].rPot[j];
		    if      (r[j] < -dHalf[j]) r[j] += dPeriod[j];
		    else if (r[j] > +dHalf[j]) r[j] -= dPeriod[j];
		    r[j] -= pkd->tinyGroupTable[gid].rcen[j];
		    r2 += r[j]*r[j];
		    }
		/*
		** If the particle lies within the shrinking sphere, then add use it to calculate the center of mass.
		*/
		if (r2 < f2*pkd->tinyGroupTable[gid].rMax*pkd->tinyGroupTable[gid].rMax) {
		    fMass = pkdMass(pkd,p);
		    for (j=0;j<3;++j) shrink[gid].rcom[j] += fMass*(r[j]+pkd->tinyGroupTable[gid].rcen[j]);
		    shrink[gid].fMass += fMass;
		    shrink[gid].nEnclosed += 1; /* count number of particles inside rmax */
		    }
		}
	    for (gid=1;gid<pkd->nGroups;++gid) {
		for (j=0;j<3;++j) shrink[gid].rcom[j] /= shrink[gid].fMass;
		}
	    /* 
	    ** Now accumulate totals globally, we just sum the 5 fields for each group.
	    */
	    mdlCOcache(mdl,CID_GROUP,NULL,shrink,sizeof(ShrinkStruct),nLocalGroups+1,NULL,initShrink,combShrink);
	    for(gid=1+nLocalGroups;gid<pkd->nGroups;++gid) {
		ShrinkStruct *g;
		if (shrink[gid].nEnclosed < 0) continue;
		g = mdlVirtualFetch(mdl,CID_GROUP,pkd->ga[gid].id.iIndex,pkd->ga[gid].id.iPid);
		for (j=0;j<3;j++) {
		    g->rcom[j] = shrink[gid].rcom[j];
		    }
		g->fMass = shrink[gid].fMass;
		g->nEnclosed = shrink[gid].nEnclosed;
		}
	    mdlFinishCache(mdl,CID_GROUP);
	    mdlROcache(mdl,CID_GROUP,NULL,shrink,sizeof(ShrinkStruct),nLocalGroups+1);
	    for(gid=1+nLocalGroups;gid<pkd->nGroups;++gid) {
		ShrinkStruct *g;
		if (shrink[gid].nEnclosed < 0) continue;
		g = mdlFetch(mdl,CID_GROUP,pkd->ga[gid].id.iIndex,pkd->ga[gid].id.iPid);
		for (j=0;j<3;j++) {
		    shrink[gid].rcom[j] = g->rcom[j];
		    }
		shrink[gid].fMass = g->fMass;
		shrink[gid].nEnclosed = g->nEnclosed;
		}    
	    mdlFinishCache(mdl,CID_GROUP);
	    for (gid=1;gid<pkd->nGroups;++gid) {
		if (shrink[gid].nEnclosed > 0) {
		    assert(shrink[gid].fMass > 0.0);
		    for (j=0;j<3;++j) {
			pkd->tinyGroupTable[gid].rcen[j] = shrink[gid].rcom[j];
			}
		    }
		}
	    }
	shrink = NULL; /* make sure we can't use it anymore by accident */
	} /* end of if (doShrinkingSphere) */ 
    /*
    ** Calculate angular momentum about the center of the group.
    */
    for (i=0;i<pkd->nLocal;++i) {
	p = pkdParticle(pkd,i);
	gid = pkdGetGroup(pkd,p);
	if (gid > 0) {
	    fMass = pkdMass(pkd,p);
	    v = pkdVel(pkd,p);
	    for (j=0;j<3;++j) {
		r[j] =  pkdPos(pkd,p,j) - pkd->tinyGroupTable[gid].rPot[j];
		if      (r[j] < -dHalf[j]) r[j] += dPeriod[j];
		else if (r[j] > +dHalf[j]) r[j] -= dPeriod[j];
		r[j] -= pkd->tinyGroupTable[gid].rcen[j];
		}
	    dAccumulate[4*gid+0] += fMass*(r[1]*v[2] - r[2]*v[1]);
	    dAccumulate[4*gid+1] += fMass*(r[2]*v[0] - r[0]*v[2]);
	    dAccumulate[4*gid+2] += fMass*(r[0]*v[1] - r[1]*v[0]);
	    }
	}
    for (gid=1;gid<pkd->nGroups;++gid) {
	for (j=0;j<3;++j) {
	    pkd->tinyGroupTable[gid].angular[j] = dAccumulate[4*gid+j];
	    dAccumulate[4*gid+j] = 0;
	    }
	}
    /*
    ** Calculate the moment of inertia (we do this in 2 stages to save accumulator memory.
    */
    for (i=0;i<pkd->nLocal;++i) {
	p = pkdParticle(pkd,i);
	gid = pkdGetGroup(pkd,p);
	if (gid > 0) {
	    fMass = pkdMass(pkd,p);
	    for (j=0;j<3;++j) {
		r[j] =  pkdPos(pkd,p,j) - pkd->tinyGroupTable[gid].rPot[j];
		if      (r[j] < -dHalf[j]) r[j] += dPeriod[j];
		else if (r[j] > +dHalf[j]) r[j] -= dPeriod[j];
		r[j] -= pkd->tinyGroupTable[gid].rcen[j];
		}
	    dAccumulate[4*gid+0] += fMass*r[0]*r[0];
	    dAccumulate[4*gid+1] += fMass*r[0]*r[1];
	    dAccumulate[4*gid+2] += fMass*r[0]*r[2];
	    }
	}
    for (gid=1;gid<pkd->nGroups;++gid) {
	for (j=0;j<3;++j) {
	    pkd->tinyGroupTable[gid].inertia[j] = dAccumulate[4*gid+j];
	    dAccumulate[4*gid+j] = 0;
	    }
	}
    for (i=0;i<pkd->nLocal;++i) {
	p = pkdParticle(pkd,i);
	gid = pkdGetGroup(pkd,p);
	if (gid > 0) {
	    fMass = pkdMass(pkd,p);
	    for (j=1;j<3;++j) {  /* careful, this loop skips r[0]! */
		r[j] =  pkdPos(pkd,p,j) - pkd->tinyGroupTable[gid].rPot[j];
		if      (r[j] < -dHalf[j]) r[j] += dPeriod[j];
		else if (r[j] > +dHalf[j]) r[j] -= dPeriod[j];
		r[j] -= pkd->tinyGroupTable[gid].rcen[j];
		}
	    dAccumulate[4*gid+0] += fMass*r[1]*r[1];
	    dAccumulate[4*gid+1] += fMass*r[1]*r[2];
	    dAccumulate[4*gid+2] += fMass*r[2]*r[2];
	    }
	}
    for (gid=1;gid<pkd->nGroups;++gid) {
	for (j=0;j<3;++j) {
	    pkd->tinyGroupTable[gid].inertia[j+3] = dAccumulate[4*gid+j];
	    dAccumulate[4*gid+j] = 0;
	    }
	}
    /* 
    ** Now accumulate angular momentum globally
    */
    mdlCOcache(mdl,CID_GROUP,NULL,pkd->tinyGroupTable,sizeof(TinyGroupTable),nLocalGroups+1,
	NULL,initAngular,combAngular);
    for(gid=1+nLocalGroups;gid<pkd->nGroups;++gid) {
	TinyGroupTable *g;
	g = mdlVirtualFetch(mdl,CID_GROUP,pkd->ga[gid].id.iIndex,pkd->ga[gid].id.iPid);
	for (j=0;j<3;j++) {
	    g->angular[j] += pkd->tinyGroupTable[gid].angular[j];
	    }
	for (j=0;j<6;++j) {
	    g->inertia[j] += pkd->tinyGroupTable[gid].inertia[j];
	    }
	}
    mdlFinishCache(mdl,CID_GROUP);
    /*
    ** Scoop environment mass (note the full tree must be present).
    ** At the very least we need to have masses in the cells, either as part of the FMOMR structure
    ** or a seperate mass field!
    */
    S = (remoteID *)(&pkd->tinyGroupTable[pkd->nGroups]);
    mdlROcache(mdl,CID_PARTICLE,NULL,pkdParticleBase(pkd),pkdParticleSize(pkd),pkdLocal(pkd));
    diVol0 = 3.0/4.0*M_1_PI*pow(rEnvironment[0],-3);
    diVol1 = 3.0/4.0*M_1_PI*pow(rEnvironment[1],-3);
    for (gid=1;gid<=nLocalGroups;++gid) {
	for (j=0;j<3;++j) r[j] = pkd->tinyGroupTable[gid].rPot[j];
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
	    pkd->tinyGroupTable[gid].angular[j] = g->angular[j];
	    }
	for (j=0;j<6;++j) {
	    pkd->tinyGroupTable[gid].inertia[j] = g->inertia[j];
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
		fabs(pkd->bndInterior.fCenter[j] - pkd->tinyGroupTable[gid].rPot[j]);
	    if (rt < rMax) rMax = rt;
	    }
	n = 0;
	for (i=iGrpEnd[gid-1];i<iGrpEnd[gid];++i,++n) {
	    p = pkdParticle(pkd,i);
	    assert(gid == pkdGetGroup(pkd,p)); /* make sure it is really in group order */
	    r2 = 0.0;
	    for (j=0;j<3;++j) {
		r[j] =  pkdPos(pkd,p,j) - pkd->tinyGroupTable[gid].rPot[j];
		if      (r[j] < -dHalf[j]) r[j] += dPeriod[j];
		else if (r[j] > +dHalf[j]) r[j] -= dPeriod[j];
		r[j] -= pkd->tinyGroupTable[gid].rcen[j];
		r2 += r[j]*r[j];
		}
	    mrFree[n].fMass = pkdMass(pkd,p);
	    mrFree[n].dr = sqrtf(r2);
	    }
	QSORT(sizeof(MassRadius),mrFree,n,mrLessThan);
	fMass = mrFree[0].fMass;
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
		r[j] =  pkdPos(pkd,p,j) - pkd->tinyGroupTable[gid].rPot[j];
		if      (r[j] < -dHalf[j]) r[j] += dPeriod[j];
		else if (r[j] > +dHalf[j]) r[j] -= dPeriod[j];
		r[j] -= pkd->tinyGroupTable[gid].rcen[j];
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


    /*
    ** IA: We copy part of the stats to a even smaller structure that will survive later in memory
    **  such that it can be used for planting the BH seeds
    */
    if (pkd->veryTinyGroupTable!=NULL) free(pkd->veryTinyGroupTable);
    pkd->veryTinyGroupTable = (VeryTinyGroupTable *) malloc( (1+nLocalGroups) * sizeof(VeryTinyGroupTable)  );
    for (gid=1; gid<=nLocalGroups; gid++){
       pkd->veryTinyGroupTable[gid].rPot[0] = pkd->tinyGroupTable[gid].rPot[0];
       pkd->veryTinyGroupTable[gid].rPot[1] = pkd->tinyGroupTable[gid].rPot[1];
       pkd->veryTinyGroupTable[gid].rPot[2] = pkd->tinyGroupTable[gid].rPot[2];
       pkd->veryTinyGroupTable[gid].fMass = pkd->tinyGroupTable[gid].fMass;
       pkd->veryTinyGroupTable[gid].nBH = pkd->tinyGroupTable[gid].nBH;
    }


    /*
    ** Clear pkd->ga so that later gravity calculations don't set the minimum potential again.
    ** This means we cannot directly output ga and could use this space for something else.
    */
    pkd->ga = NULL;
    }
