#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
const char *analysis_module_id = "$Id$";
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "pkd.h"

#define SHAPES

static void transpose(double mat[3][3],double trans_mat[3][3]) {
    int i,j ;

    for (i=0; i<3; i++) {
        for (j=0; j<3; j++) {
            trans_mat[i][j] = mat[j][i];
	    }
	}
    }


#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
        a[k][l]=h+s*(g-h*tau);

static void jacobi(double a[4][4],int n,double d[4],double v[4][4], int *nrot) {
    int j,iq,ip,i;
    double tresh,theta,tau,t,sm,s,h,g,c ;
    double b[4],z[4] ;

    for (ip=1;ip<=n;ip++) {
	for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
	v[ip][ip]=1.0;
        }
    for (ip=1;ip<=n;ip++) {
	b[ip]=d[ip]=a[ip][ip];
	z[ip]=0.0;
        }
    *nrot=0;
    for (i=1;i<=100;i++) {
	sm=0.0;
	for (ip=1;ip<=n-1;ip++) {
	    for (iq=ip+1;iq<=n;iq++)
		sm += fabs(a[ip][iq]);
	    }
	if (sm <= 1.0e-12) {
	    return;
	    }
	if (i < 4)
	    tresh=0.2*sm/(n*n);
	else
	    tresh=0.0;
	for (ip=1;ip<=n-1;ip++) {
	    for (iq=ip+1;iq<=n;iq++) {
		g=100.0*fabs(a[ip][iq]);
		if (i > 4 && fabs(d[ip])+g == fabs(d[ip])
		    && fabs(d[iq])+g == fabs(d[iq]))
		    a[ip][iq]=0.0;
		else if (fabs(a[ip][iq]) > tresh) {
		    h=d[iq]-d[ip];
		    if (fabs(h)+g == fabs(h))
			t=(a[ip][iq])/h;
		    else {
			theta=0.5*h/(a[ip][iq]);
			t=1.0/(fabs(theta)+
			       sqrt(1.0+theta*theta));
			if (theta < 0.0) t = -t;
			}
		    c=1.0/sqrt(1+t*t);
		    s=t*c;
		    tau=s/(1.0+c);
		    h=t*a[ip][iq];
		    z[ip] -= h;
		    z[iq] += h;
		    d[ip] -= h;
		    d[iq] += h;
		    a[ip][iq]=0.0;
		    for (j=1;j<=ip-1;j++) {
			ROTATE(a,j,ip,j,iq)
			    }
		    for (j=ip+1;j<=iq-1;j++) {
			ROTATE(a,ip,j,j,iq)
			    }
		    for (j=iq+1;j<=n;j++) {
			ROTATE(a,ip,j,iq,j)
			    }
		    for (j=1;j<=n;j++) {
			ROTATE(v,j,ip,j,iq)
			    }
		    ++(*nrot);
		    }
		}
	    }
	for (ip=1;ip<=n;ip++) {
	    b[ip] += z[ip];
	    d[ip]=b[ip];
	    z[ip]=0.0;
	    }
        }
    printf("<error in jacobi>\n") ;
    }

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

#ifdef SHAPES
static void combShapesBins1(void *vpkd, void *b1, void *b2) {
    SHAPESBIN * pBin1 = (SHAPESBIN *)b1;
    SHAPESBIN * pBin2 = (SHAPESBIN *)b2;
    int j,k;
    pBin1->dMassEnclosed += pBin2->dMassEnclosed;
    for (j=0; j<3; j++) {
        pBin1->com[j] += pBin2->com[j];
        for (k=0; k<3; k++) {
            pBin1->dInertia[j][k] += pBin2->dInertia[j][k];
	    }
	}
    }
static void initShapesBins1(void *vpkd, void *b) {
    SHAPESBIN * pBin = (SHAPESBIN *)b;
    int j,k;
    pBin->dMassEnclosed = 0.0;
    for (j=0; j<3; j++) {
        pBin->com[j] = 0.0;
        for (k=0; k<3; k++) {
            pBin->dInertia[j][k] = 0.0;
	    }
	}
    }
#endif

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

#ifdef SHAPES

double ell_distance2(const double *r,SHAPESBIN *pShape, double ba, double ca) {
    double dx[3], dx_rot[3];
    int i;

    for ( i=0; i<3; i++ ) dx[i] = r[i] - pShape->ell_center[i];
    matrix_vector_mult(dx_rot,pShape->ell_matrix,dx);
    return dx_rot[0] * dx_rot[0]
	+ dx_rot[1] * dx_rot[1] / ba / ba
	+ dx_rot[2] * dx_rot[2] / ca / ca;
    }


/*
** NB: We assume that shells do NOT cross.  In other words, if a particle is part
**     of bin n, then it is always part of bin n+1, n+2, ... etc.  In principle
**     this is not necessarily true, but it would be extremely odd for the centers
**     to be off by that much.  This assumption greatly improves the performance,
**     taking the algorithm from Nbin * Npart down to Npart.
*/
static void CalculateInertia(PKD pkd,int nBins, const double *dRadii, SHAPESBIN *shapesBins) {
    PLITE *pl = pkd->pLite;
    local_t n = pkdLocal(pkd);
    SHAPESBIN *pShape;
    double r, r2;
    int i, j, k;
    int iBin;
    double ell_matrix[3][3], ell_matrix_inv[3][3];


    mdlCOcache(pkd->mdl,CID_SHAPES,NULL,shapesBins,sizeof(SHAPESBIN),pkd->idSelf==0?nBins:0,pkd,initShapesBins1,combShapesBins1);

    /*
    ** The most efficient way to handle this is to do the calculations for all bins
    */
    for(i=iBin=0;i<n;i++) {

	PARTICLE *p = pkdParticle(pkd,pl[i].i);
	double m = pkdMass(pkd,p);

	r = dRadii[iBin];
	r2 = r*r;

	/* Find the bin: Assume that the last particle was close to the correct bin */
	while( pl[i].r[0]<dRadii[iBin]*dRadii[iBin] && iBin < nBins ) ++iBin;
	while( pl[i].r[0]>dRadii[iBin]*dRadii[iBin] ) --iBin;

	pShape = mdlAquire(pkd->mdl,CID_SHAPES,iBin,0);
	pShape->dMassEnclosed += m;
	for (j=0; j<3; j++) {
	    pShape->com[j] += m * p->r[j] ;
	    for (k=0; k<=j; k++) {
		pShape->dInertia[j][k] += m * p->r[j] * p->r[k];
		}
	    }
	mdlRelease(pkd->mdl,CID_SHAPES,pShape);

	}
    mdlFinishCache(pkd->mdl,CID_SHAPES);

    if (pkd->idSelf == 0) {
	double inertia_cm[4][4];
	double evectors[4][4];
	double evalues[4];
	double VecProd[4], ScalProd;
	double ba, ca, theta, phi, psi;
	int nrot;
	int ia, ib, ic;
	for(i=iBin=0;iBin<nBins;iBin++) {
	    pShape = &shapesBins[iBin];

	    for (j=0; j<3; j++) {
		for (k=0; k<3; k++) {
		    if (k<=j) inertia_cm[j+1][k+1] = ((pShape->dInertia[j][k] - pShape->com[j] * pShape->com[k] / pShape->dMassEnclosed) / pShape->dMassEnclosed);
		    else      inertia_cm[j+1][k+1] = ((pShape->dInertia[k][j] - pShape->com[k] * pShape->com[j] / pShape->dMassEnclosed) / pShape->dMassEnclosed);
		    }
		evalues[j+1] = inertia_cm[j+1][j+1] ;
		}
	    jacobi(inertia_cm,3,evalues,evectors,&nrot) ;
	    if(evalues[1] >= evalues[2] && evalues[1] >= evalues[3]){
		ia = 1 ;
		if(evalues[2] >= evalues[3]){
		    ib = 2 ;
		    ic = 3 ;
		    }
		else{
		    ib = 3 ;
		    ic = 2 ;
		    }
		}
	    else if(evalues[2] > evalues[1] && evalues[2] >= evalues[3]){
		ia = 2 ;
		if(evalues[1] >= evalues[3]){
		    ib = 1 ;
		    ic = 3 ;
		    }
		else{
		    ib = 3 ;
		    ic = 1 ;
		    }
		}
	    else{
		ia = 3 ;
		if(evalues[1] >= evalues[2]){
		    ib = 1 ;
		    ic = 2 ;
		    }
		else{
		    ib = 2 ;
		    ic = 1 ;
		    }
		}

	    /* Check if Eigenvectors are righthanded in 3D :
	       ev[ib] x ev[ic] = ev[ia] */

	    VecProd[1] =  evectors[2][ib]*evectors[3][ic]
		- evectors[3][ib]*evectors[2][ic];
	    VecProd[2] = -evectors[1][ib]*evectors[3][ic]
		+ evectors[3][ib]*evectors[1][ic];
	    VecProd[3] =  evectors[1][ib]*evectors[2][ic]
		- evectors[2][ib]*evectors[1][ic];
	    ScalProd   =  evectors[1][ia]*VecProd[1] + evectors[2][ia]*VecProd[2]
		+ evectors[3][ia]*VecProd[3];
	    if (ScalProd < 0.0) {
		for(i=0; i<3; i++) evectors[i+1][ia] = -evectors[i+1][ia];
		}

	    ba = sqrt((double)(evalues[ib]/evalues[ia])) ;
	    ca = sqrt((double)(evalues[ic]/evalues[ia])) ;

	    /* euler angles for a zyz rotation */
	    theta = 180. / M_PI * acos((double) evectors[3][ic]);
	    phi =   180. / M_PI * acos((double) evectors[1][ic]/sqrt(evectors[1][ic]*evectors[1][ic] + evectors[2][ic]*evectors[2][ic]));
	    psi =   180. / M_PI * acos((double) (-evectors[2][ic]*evectors[1][ib] + evectors[1][ic]*evectors[2][ib])/
				     sqrt(evectors[1][ic]*evectors[1][ic] + evectors[2][ic]*evectors[2][ic]));

	    /* inverse acos is only defined between 0 and pi therefore we must
	       deal with pi to 2*pi */
	    if(evectors[2][ic] < 0.0) phi = 360. - phi; /* phi always positive */
	    if(evectors[3][ib] < 0.0) psi = 360. - psi; /* psi always positive */ 

/*
	    for(i = 0; i<3; i++){
		ell_center[i] = pShape->com[i] / pShape->dMassEnclosed;
		}
*/
	    for(i=0;i<3;i++){
		ell_matrix_inv[i][0] = evectors[i+1][ia];
		ell_matrix_inv[i][1] = evectors[i+1][ib];
		ell_matrix_inv[i][2] = evectors[i+1][ic];
		}
	    transpose(ell_matrix_inv,ell_matrix);
	    }
	}


    }

/*
** Calculate shapes for existing profile bins
*/
void pkdShapes(PKD pkd, int nBins, const double *dCenter, const double *dRadii) {
    SHAPESBIN *shapesBins;
    int i, j, k;

    if (pkd->idSelf == 0) {
	shapesBins = mdlMalloc(pkd->mdl,nBins*sizeof(SHAPESBIN));
	assert( shapesBins != NULL );
	/* Start with the given center for every bin */
	for( i=0; i<nBins; i++ ) {
	    for (j=0; j<3; j++) {
		shapesBins[i].ell_center[j] = dCenter[j];
		for (k=0; k<=j; k++) {
		    shapesBins[i].ell_matrix[j][k] = (j==k) ? 1.0 : 0.0;
		    }
		}
	    }
	}
    else shapesBins = NULL;

    mdlROcache(pkd->mdl,CID_BIN,NULL,pkd->profileBins,sizeof(PROFILEBIN),pkd->idSelf?0:nBins);
    CalculateInertia(pkd,nBins,dRadii,shapesBins);
    mdlFinishCache(pkd->mdl,CID_BIN);

    if (pkd->idSelf == 0) {
	mdlFree(pkd->mdl,shapesBins);
	}


    }
#endif

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
	    /*double vel_tang[3], vel_shell[3], vel_tang_pec[3];*/

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

#ifdef MDL_FFTW

static void initPk(void *vpkd, void *g) {
    fftw_real * r = (fftw_real *)g;
    *r = 0.0;
    }
static void combPk(void *vpkd, void *g1, void *g2) {
    fftw_real * r1 = (fftw_real *)g1;
    fftw_real * r2 = (fftw_real *)g2;

    *r1 += *r2;
    }

static int wrap(int i,int w) {
    if (i>=w) i -= w;
    else if (i<0) i += w;
    return i;
    }

static inline float pow2(float x) {
    return x*x;
    }

double grid_rms(MDL mdl, int nGrid, float *fGrid) {
    int i;
    double sum2, total;

    sum2 = 0.0;
    for( i=0; i<nGrid; i++ )
	sum2 += fGrid[i] * fGrid[i];
    mdlReduce(mdl,&sum2,&total,1,MDL_DOUBLE,MDL_SUM,0);
    return sqrt(total);
    }

double grid_mass(MDL mdl, int nGrid, float *fGrid) {
    int i,j,k;
    double sum, total;

    sum = 0.0;
    for( i=0; i<nGrid; i++ )
	sum += fGrid[i];
    mdlAllreduce(mdl,&sum,&total,1,MDL_DOUBLE,MDL_SUM);
    return total;
    }

/*
** Given grid coordinates x,y,z (integer), find out on which processor
** this cell can be found and accumulate mass into it.
*/
static void cell_accumulate(PKD pkd, MDLFFT fft,int x,int y,int z, float m) {
    fftw_real *p;
    int id, idx;

    /* Map coordinate to processor/index */
    id = mdlFFTrId(fft,x,y,z);
    idx = mdlFFTrIdx(fft,x,y,z);

    //if (id == pkd->idSelf) {
    p = mdlAquire(fft->mdl,CID_PK,idx,id);
    *p += m;
    mdlRelease(fft->mdl,CID_PK,p);
    }

static void tsc_assign(PKD pkd, MDLFFT fft, double x, double y, double z, double mass) {
    int           ix, iy, iz, ixp1, iyp1, izp1, ixm1, iym1, izm1;
    float         rrx, rry, rrz;
    float         hx, hy, hz;
    float         hx0, hy0, hz0, hxp1, hyp1, hzp1, hxm1, hym1, hzm1;
    float         pw;
   
    /* coordinates in subcube units [0,NGRID] */
    rrx = x * (float)(fft->n1);
    rry = y * (float)(fft->n2);
    rrz = z * (float)(fft->n3);
               
    /* index of nearest grid point [0,NGRID] */
    ix  = (int)(rrx+0.5);
    iy  = (int)(rry+0.5);
    iz  = (int)(rrz+0.5);
               
    /* distance to nearest grid point */
    hx  = rrx - (float)ix;
    hy  = rry - (float)iy;
    hz  = rrz - (float)iz;
               
    /* keep track of peridoc boundaries -> [0,NGRID-1] ; NGRID=0  */
    ix = wrap(ix,fft->n1);
    iy = wrap(iy,fft->n2);
    iz = wrap(iz,fft->n3);
               
    /* particle mass */
    pw = mass;
               
    /* calculate TSC weights */
    hx0=0.75 - hx*hx;
    hxp1=0.5* pow2(0.5 + hx);
    hxm1=0.5* pow2(0.5 - hx);
    hy0=0.75 - hy*hy;
    hyp1=0.5* pow2(0.5 + hy);
    hym1=0.5* pow2(0.5 - hy);
    hz0= 0.75 - hz*hz;
    hzp1=0.5* pow2(0.5 + hz);
    hzm1=0.5* pow2(0.5 - hz);
               
    /* keep track of peridoc boundaries */
    ixp1=wrap(ix+1,fft->n1);
    iyp1=wrap(iy+1,fft->n2);
    izp1=wrap(iz+1,fft->n3);
    ixm1=wrap(ix-1,fft->n1);
    iym1=wrap(iy-1,fft->n2);
    izm1=wrap(iz-1,fft->n3);

    /* assign particle according to weights to 27 neighboring nodes */
    cell_accumulate(pkd,fft,ixm1,iym1,izm1,hxm1*hym1 *hzm1 * pw);
    cell_accumulate(pkd,fft,ix,  iym1,izm1,hx0 *hym1 *hzm1 * pw);
    cell_accumulate(pkd,fft,ixp1,iym1,izm1,hxp1*hym1 *hzm1 * pw);
    cell_accumulate(pkd,fft,ixm1,  iy,izm1,hxm1*hy0  *hzm1 * pw);
    cell_accumulate(pkd,fft,  ix,  iy,izm1,hx0 *hy0  *hzm1 * pw);
    cell_accumulate(pkd,fft,ixp1,  iy,izm1,hxp1*hy0  *hzm1 * pw);
    cell_accumulate(pkd,fft,ixm1,iyp1,izm1,hxm1*hyp1 *hzm1 * pw);
    cell_accumulate(pkd,fft,  ix,iyp1,izm1,hx0 *hyp1 *hzm1 * pw);
    cell_accumulate(pkd,fft,ixp1,iyp1,izm1,hxp1*hyp1 *hzm1 * pw);
    cell_accumulate(pkd,fft,ixm1,iym1,  iz,hxm1*hym1 *hz0 * pw);
    cell_accumulate(pkd,fft,  ix,iym1,  iz,hx0 *hym1 *hz0 * pw);
    cell_accumulate(pkd,fft,ixp1,iym1,  iz,hxp1*hym1 *hz0 * pw);
    cell_accumulate(pkd,fft,ixm1,  iy,  iz,hxm1*hy0  *hz0 * pw);
    cell_accumulate(pkd,fft,  ix,  iy,  iz,hx0 *hy0  *hz0 * pw);
    cell_accumulate(pkd,fft,ixp1,  iy,  iz,hxp1*hy0  *hz0 * pw);
    cell_accumulate(pkd,fft,ixm1,iyp1,  iz,hxm1*hyp1 *hz0 * pw);
    cell_accumulate(pkd,fft,  ix,iyp1,  iz,hx0 *hyp1 *hz0 * pw);
    cell_accumulate(pkd,fft,ixp1,iyp1,  iz,hxp1*hyp1 *hz0 * pw);
    cell_accumulate(pkd,fft,ixm1,iym1,izp1,hxm1*hym1 *hzp1 * pw);
    cell_accumulate(pkd,fft,  ix,iym1,izp1,hx0 *hym1 *hzp1 * pw);
    cell_accumulate(pkd,fft,ixp1,iym1,izp1,hxp1*hym1 *hzp1 * pw);
    cell_accumulate(pkd,fft,ixm1,  iy,izp1,hxm1*hy0  *hzp1 * pw);
    cell_accumulate(pkd,fft,  ix,  iy,izp1,hx0 *hy0  *hzp1 * pw);
    cell_accumulate(pkd,fft,ixp1,  iy,izp1,hxp1*hy0  *hzp1 * pw);
    cell_accumulate(pkd,fft,ixm1,iyp1,izp1,hxm1*hyp1 *hzp1 * pw);
    cell_accumulate(pkd,fft,  ix,iyp1,izp1,hx0 *hyp1 *hzp1 * pw);
    cell_accumulate(pkd,fft,ixp1,iyp1,izp1,hxp1*hyp1 *hzp1 * pw);
}

void pkdMeasurePk(PKD pkd, double dCenter[3], double dRadius,
		  int nGrid, float *fPower, int *nPower) {
    PARTICLE *p;
    MDLFFT fft;
    fftw_real *fftData;
    fftw_complex *fftDataK;
    double dTotalMass, dRhoMean, diRhoMean;
    double fftNormalize;
    double rms;
    double r[3];
    double dScale;
    int i,j,k, idx, ks;
    int iNyquist = nGrid / 2;

    /* Box scaling factor */
    dScale = 0.5 / dRadius;

    fftNormalize = 1.0 / (1.0*nGrid*nGrid*nGrid);

    mdlFFTInitialize(pkd->mdl,&fft,nGrid,nGrid,nGrid,0);
    fftData = mdlFFTMAlloc( fft );
    fftDataK = (fftw_complex *)fftData;

    for( i=0; i<fft->nlocal; i++ ) fftData[i] = 0.0;

    mdlCOcache(pkd->mdl,CID_PK,NULL,fftData,sizeof(fftw_real), fft->nlocal,pkd,initPk,combPk);
    for (i=0;i<pkd->nLocal;++i) {
	p = pkdParticle(pkd,i);
	if ( !pkdIsSrcActive(p,0,MAX_RUNG) ) continue;
	/* Recenter, apply periodic boundary and scale to the correct size */
	for(j=0;j<3;j++) {
	    r[j] = p->r[j] - dCenter[j] + dRadius;
	    if ( pkd->param.bPeriodic ) {
		if ( r[j] >= pkd->fPeriod[j] ) r[j] -= pkd->fPeriod[j];
		else if ( r[j] < 0.0 ) r[j] += pkd->fPeriod[j];
		}
	    r[j] *= dScale;
	    }
	/* 
	** The position has been rescaled to [0,1).  If it is not in that range,
	** then this particle is outside of the given box and should be ignored.
	*/
	if( r[0]>=0.0 && r[0]<1.0 && r[1]>=0.0 && r[1]<1.0 && r[2]>=0.0 && r[2]<1.0 )
	    tsc_assign(pkd, fft, r[0], r[1], r[2], pkdMass(pkd,p));
	}
    mdlFinishCache(pkd->mdl,CID_PK);

    dTotalMass = grid_mass(pkd->mdl,fft->nlocal, fftData);
    dRhoMean = dTotalMass * fftNormalize;
    diRhoMean = 1.0 / dRhoMean;

    /*rms = grid_rms(pkd->mdl,fft->nlocal, fftData);
      if (pkd->idSelf==0)	printf( "RMS after TSC: %g\n",rms);*/

    /*printf( "Calculating density contrast\n" );*/
    for( i=0; i<fft->nlocal; i++ ) {
	fftData[i] = (fftData[i] - dRhoMean) * diRhoMean;
	}
    /*rms = grid_rms(pkd->mdl,fft->nlocal, fftData);
      if (pkd->idSelf==0) printf( "RMS after Density: %g\n",rms);*/

    mdlFFT(fft,fftData,0);
    // Remember, the grid is now transposed to x,z,y (from x,y,z)

    /*rms = grid_rms(pkd->mdl,fft->nlocal, fftData);
      if (pkd->idSelf==0)	printf( "RMS after FFT: %g\n",rms);*/


    /*printf("Calculating delta^2\n");*/
    for( i=0; i<fft->nlocal/2; i++ ) {
	c_re(fftDataK[i]) = pow2(c_re(fftDataK[i])) + pow2(c_im(fftDataK[i]));
	c_im(fftDataK[i]) = 0.0;
	}

    /*rms = grid_rms(pkd->mdl,fft->nlocal, fftData);
      if (pkd->idSelf==0) printf( "RMS after Delta: %g\n",rms);*/

    /*printf("Calculating P(k)\n");*/
    for( i=0; i<=iNyquist; i++ ) {
	fPower[i] = 0.0;
	nPower[i] = 0;
	}

    for(j=fft->sy;j<fft->sy+fft->ny;j++) {
	int jj = j>iNyquist ? fft->n2 - j : j;
	for(k=0;k<fft->n3;k++) {
	    int kk = k>iNyquist ? fft->n3 - k : k;
	    for(i=0;i<fft->a1k;i++) {
		ks = sqrtl(i*i + jj*jj + kk*kk);
		idx = mdlFFTkIdx(fft,i,j,k);
		if ( ks >= 1 && ks <= iNyquist ) {
		    fPower[ks] += c_re(fftDataK[idx]);
		    nPower[ks] += 1;
		    }
		}
	    }
	}

    mdlFFTFree(fft,fftData);
    mdlFFTFinish(fft);
    }
#endif
