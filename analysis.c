#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <assert.h>
#include <math.h>
#include "pkd.h"

/*
** Combiner cache functions
*/
static void combProfileBins1(void *vpkd, void *b1, void *b2) {
    PROFILEBIN * pBin1 = (PROFILEBIN *)b1;
    PROFILEBIN * pBin2 = (PROFILEBIN *)b2;

    pBin1->dMassInBin += pBin2->dMassInBin;
    pBin1->nParticles += pBin2->nParticles;
    pBin1->L[0] += pBin2->L[0];
    pBin1->L[1] += pBin2->L[1];
    pBin1->L[2] += pBin2->L[2];
    pBin1->vel_radial += pBin2->vel_radial;
    pBin1->vel_radial_sigma += pBin2->vel_radial_sigma;
    }

static void initProfileBins1(void *vpkd, void *b) {
    PROFILEBIN * pBin = (PROFILEBIN *)b;
    pBin->dMassInBin = 0.0;
    pBin->L[0] = pBin->L[1] = pBin->L[2] = 0.0;
    pBin->vel_radial = pBin->vel_radial_sigma = 0.0;
    pBin->nParticles = 0;
    }

static void combProfileBins2(void *vpkd, void *b1, void *b2) {
    PROFILEBIN * pBin1 = (PROFILEBIN *)b1;
    PROFILEBIN * pBin2 = (PROFILEBIN *)b2;
    pBin1->vel_tang_sigma += pBin2->vel_tang_sigma;
    }

static void initProfileBins2(void *vpkd, void *b) {
    PROFILEBIN * pBin = (PROFILEBIN *)b;
    pBin->vel_tang_sigma = 0.0;
    }

/*
** This function will calculate the distance between a particle and a
** reference point.  If a periodic boundary is in effect then the smallest
** possible distance is returned.
*/
double pkdGetDistance2(PKD pkd,PARTICLE *p, const double *dCenter ) {
    double d2;
    double dx,dx2;
    int j;

    d2 = 0.0;
    for( j=0; j<3; j++ ) {
	dx = p->r[j] - dCenter[j];
	/*
	** If a periodic wrap results in a smaller distance, then use that.
	*/
	if ( pkd->param.bPeriodic ) {
	    if ( dx<0.0 ) dx2 = dx + pkd->fPeriod[j];
	    else dx2 = dx - pkd->fPeriod[j];
	    if ( dx2*dx2 < dx*dx ) dx = dx2;
	    }
	d2 += dx*dx;
	}
    return d2;
    }

/*
** Count the number of particles in a given shell from
** [dMinRadius,dMaxRadius).
*/
int pkdShellCount(PKD pkd, uint8_t uRungLo, uint8_t uRungHi,
		  double *dCenter, double dMinRadius, double dMaxRadius ) {
    PARTICLE *p;
    int n, i, iCount;
    double d2,r2min, r2max;

    r2min = dMinRadius*dMinRadius;
    r2max = dMaxRadius*dMaxRadius;
    n = pkdLocal(pkd);
    iCount = 0;
    for (i=0;i<n;++i) {
	p = pkdParticle(pkd,i);
	if (pkdIsSrcActive(p,uRungLo,uRungHi)) {
	    d2 = pkdGetDistance2(pkd,p,dCenter);
	    if ( d2>=r2min || d2 < r2max )
		iCount ++;
	    }
	}

    return iCount;
    }

static int cmpRadiusLite(const void *pva,const void *pvb) {
    PLITE *pa = (PLITE *)pva;
    PLITE *pb = (PLITE *)pvb;
    double d = pa->r[0] - pb->r[0];
    if ( d > 0 ) return 1;
    else if ( d < 0 ) return -1;
    return 0;
    }

/*
** Use the pLite structure to calculate the distance to each particle
** Sort by distance when finished.
*/
void pkdCalcDistance(PKD pkd, double *dCenter) {
    PLITE *pl = pkd->pLite;
    int i;

    /*
    ** Initialize the temporary particles.
    */
    for (i=0;i<pkd->nLocal;++i) {
	PARTICLE *p = pkdParticle(pkd,i);
	double m = pkdMass(pkd,p);
	double *v = pkdVel(pkd,p);
	pl[i].r[0] = pkdGetDistance2(pkd,p,dCenter);
	pl[i].r[1] = m;
	pl[i].r[2] = 0.0;
	pl[i].i = i;
	}
    qsort(pkd->pLite,pkdLocal(pkd),sizeof(PLITE),cmpRadiusLite);
    }

/*
** Return the mass weighted center of mass and velocity
*/
void pkdCalcCOM(PKD pkd, double *dCenter, double dRadius,
		double *com, double *vcm, double *L,
		double *M, uint64_t *N) {
    double d2, dRadius2, T[3];
    int i;

    for( i=0; i<3; i++ ) com[i] = vcm[i] = L[i] = 0.0;
    *M = 0.0;
    *N = 0;
    dRadius2 = dRadius * dRadius;
    for (i=0;i<pkd->nLocal;++i) {
	PARTICLE *p = pkdParticle(pkd,i);
	double m = pkdMass(pkd,p);
	double *v = pkdVel(pkd,p);
	d2 = pkdGetDistance2(pkd,p,dCenter );
	if ( d2 < dRadius2 ) {
	    *M += m;
	    vec_add_const_mult(com, com, m, p->r);
	    vec_add_const_mult(vcm, vcm, m, v);
	    cross_product(T, p->r, v);
	    vec_add_const_mult(L, L, m, T);
	    (*N)++;
	    }
	}
    }

/*
** Count the number of elements that are interior to r2
*/
uint_fast32_t pkdCountDistance(PKD pkd, double r2i, double r2o ) {
    PLITE *pl = pkd->pLite;
    uint64_t lo,hi,i,upper;

    lo = 0;
    hi = pkd->nLocal;
    while( lo<hi ) {
	i = (lo+hi) / 2;
	if ( pl[i].r[0] >= r2o ) hi = i;
	else lo = i+1;
	}
    upper = hi;

    lo = 0;
    while( lo<hi ) {
	i = (lo+hi) / 2;
	if ( pl[i].r[0] >= r2i ) hi = i;
	else lo = i+1;
	}

    return upper-hi;
    }

/*
** Density Profile
*/
void pkdProfile(PKD pkd, uint8_t uRungLo, uint8_t uRungHi,
		const double *dCenter, const double *dRadii, int nBins,
		const double *com, const double *vcm, const double *L) {
    PLITE *pl = pkd->pLite;
    local_t n = pkdLocal(pkd);
    double r0, r, r2;
    int i,iBin;
    PROFILEBIN *pBin;

    if (pkd->idSelf == 0) {
	if ( pkd->profileBins != NULL ) mdlFree(pkd->mdl,pkd->profileBins);
	pkd->profileBins = mdlMalloc(pkd->mdl,nBins*sizeof(PROFILEBIN));
	assert( pkd->profileBins != NULL );
	r0 = 0.0;
	for(iBin=0;iBin<nBins;iBin++) {
	    pBin = &pkd->profileBins[iBin];
	    r = dRadii[iBin];
	    r2 = r*r;
	    assert( r > r0 );
	    pBin->nParticles = 0;
	    pBin->dMassInBin = 0.0;
	    pBin->dRadius = r;
	    pBin->dVolume = (4.0/3.0) * M_PI * (r*r2 - r0*r0*r0);
	    pBin->L[0] = pBin->L[1] = pBin->L[2] = 0.0;
	    pBin->vel_radial = pBin->vel_radial_sigma = 0.0;
	    pBin->vel_tang_sigma = 0.0;
	    r0 = r;
	    }
	}
    mdlCOcache(pkd->mdl,CID_BIN,NULL,pkd->profileBins,sizeof(PROFILEBIN),pkd->idSelf==0?nBins:0,pkd,initProfileBins1,combProfileBins1);

    /*
    ** Now we add all of the particles to the appropriate bin.  NOTE than both
    ** the bins, and the particles are already sorted by distance from the center.
    */
    r0 = 0.0;
    i = 0;
    for(iBin=0;iBin<nBins;iBin++) {
	pBin = mdlAquire(pkd->mdl,CID_BIN,iBin,0);
	r = pBin->dRadius;
	r = dRadii[iBin];
	r2 = r*r;
	assert( r > r0 );

	while( pl[i].r[0] <= r2 && i<n) {
	    PARTICLE *p = pkdParticle(pkd,pl[i].i);
	    double m = pkdMass(pkd,p);
	    double *v = pkdVel(pkd,p);
	    double delta_x[3], delta_v[3], ang_mom[3], dx2, vel;
	    double vel_tang[3], vel_shell[3], vel_tang_pec[3];

	    pBin->dMassInBin += pl[i].r[1];
	    pBin->nParticles++;

	    vec_sub(delta_x,p->r,com);
	    vec_sub(delta_v,v,vcm);
	    cross_product(ang_mom,delta_x,delta_v);
	    vec_add_const_mult(pBin->L,pBin->L,m,ang_mom);
	    dx2 = dot_product(delta_x,delta_x);
	    if(dx2 != 0.0)
		vel = dot_product(delta_x,delta_v)/sqrt(dx2);
	    else
		vel = sqrt(dot_product(delta_v,delta_v));
	    pBin->vel_radial += m * vel;
	    pBin->vel_radial_sigma += m * vel * vel;
	    i++;
	    }
	r0 = r;
	mdlRelease(pkd->mdl,CID_BIN,pBin);
	}
    mdlFinishCache(pkd->mdl,CID_BIN);

    /* We need angular momentum to calculate tangental velocity sigma, so we reopen the cache */
    mdlCOcache(pkd->mdl,CID_BIN,NULL,pkd->profileBins,sizeof(PROFILEBIN),pkd->idSelf?0:nBins,pkd,initProfileBins2,combProfileBins2);
    r0 = 0.0;
    i = 0;
    for(iBin=0;iBin<nBins;iBin++) {
	pBin = mdlAquire(pkd->mdl,CID_BIN,iBin,0);
	r = dRadii[iBin];
	r2 = r*r;
	assert( r > r0 );
	while( pl[i].r[0] <= r2 && i<n) {
	    PARTICLE *p = pkdParticle(pkd,pl[i].i);
	    double m = pkdMass(pkd,p);
	    double *v = pkdVel(pkd,p);
	    double delta_x[3], delta_v[3], dx2;
	    double vel_tang[3], vel_shell[3], vel_tang_pec[3];

	    vec_sub(delta_x,p->r,com);
	    vec_sub(delta_v,v,vcm);
	    dx2 = dot_product(delta_x,delta_x);

	    if(dx2 != 0.0) {
		vec_add_const_mult(vel_tang,delta_v,-dot_product(delta_v,delta_x) / dx2, delta_x);
		cross_product(vel_shell,pBin->L,delta_x);
		vec_add_const_mult(vel_tang_pec,vel_tang,-1.0 / dx2,vel_shell);
		pBin->vel_tang_sigma += m * dot_product(vel_tang_pec,vel_tang_pec);
		}
	    i++;
	    }
	r0 = r;
	mdlRelease(pkd->mdl,CID_BIN,pBin);
	}
    mdlFinishCache(pkd->mdl,CID_BIN);

    /* Only the main processor needs the result */
#if 0
    if (pkd->idSelf != 0) {
	mdlFree(pkd->mdl,pkd->profileBins);
	pkd->profileBins = NULL;
	}
#endif
    }

/*
** Given a position, find the processor on which it would be found.
*/
int pkdFindProcessor(const PKD pkd, const FLOAT *R) {
    const KDN *Top = pkd->kdTop;
    int iCell, iLower;

    iCell = ROOT;
    assert ( IN_BND(R,&Top[iCell].bnd) );

    /* Descend the tree until we find the processor for point R */
    while( (iLower=Top[iCell].iLower) ) {
	if ( !IN_BND(R,&Top[++iCell].bnd) ) iCell = iLower;
	}

    return Top[iCell].pLower;
    }

/*
** Perform a transformation on all particles.  The source position is taken
** from the oSource field (normally r[]) and the result is put into the
** oResult field (which can also be r[]).
*/
void pkdTransform(PKD pkd, int oSource, int oResult, const double *dRotCenter,
		  const double *dRotate, const double *dRecenter) {
    PARTICLE *p;
    double *ps, *pr;
    double r[3];
    int i,j;

    assert(0);
    for (i=0;i<pkd->nLocal;++i) {
	p = pkdParticle(pkd,i);
	ps = pkdField(p,oSource);
	pr = pkdField(p,oResult);

	for(j=0; j<3; j++) r[j] = ps[j] - dRotCenter[j];



	for(j=0; j<3; j++) pr[j] = r[j] + dRotCenter[j] - dRecenter[j];
	}
    }
