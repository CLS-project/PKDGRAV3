#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <assert.h>
#include <math.h>
#include "pkd.h"

double illinois(double (*func)(double,void *),void *ctx,double r,double s,double xacc,int *pnIter) {
    const int maxIter = 100;
    double t,fr,fs,ft,phis,phir,gamma;
    int i;

    fr = func(r,ctx);
    fs = func(s,ctx);
    assert(fr*fs < 0);
    t = (s*fr - r*fs)/(fr - fs);
    for (i=0;i < maxIter && fabs(t-s) > xacc;++i) {
	ft = func(t,ctx);
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
    }
    *pnIter = i;
    return(t);
}


static void combProfileBins(void *b1, void *b2) {
    PROFILEBIN * pBin1 = (PROFILEBIN *)b1;
    PROFILEBIN * pBin2 = (PROFILEBIN *)b2;

    pBin1->dMassInBin += pBin2->dMassInBin;
    pBin1->nParticles += pBin2->nParticles;
    }

static void initProfileBins(void *b) {
    PROFILEBIN * pBin = (PROFILEBIN *)b;
    pBin->dRadius = 0.0;
    pBin->dMassInBin = 0.0;
    pBin->dVolume = 0.0;
    pBin->nParticles = 0;
    }




/*
** This function will calculate the distance between a particle and a
** reference point.  If a periodic boundary is in effect then the smallest
** possible distance is returned.
*/
double pkdGetDistance(PKD pkd,PARTICLE *p, double *dCenter ) {
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
	p = &pkd->pStore[i];
	if (pkdIsSrcActive(p,uRungLo,uRungHi)) {
	    d2 = pkdGetDistance(pkd,p,dCenter);
	    if ( d2>=r2min || d2 < r2max )
		iCount ++;
	    }
	}

    return iCount;
    }

static int cmpRadiusLite(const void *pva,const void *pvb) {
    PLITE *pa = (PLITE *)pva;
    PLITE *pb = (PLITE *)pvb;
    return(pa->r[0] - pb->r[0]);
    }

/*
** Use the pLite structure to calculate the distance to each particle
** Sort by Radius when finished.
*/
void pkdCalcRadius(PKD pkd, double *dCenter) {
    PLITE *pl = pkd->pLite;
    int i;

    /*
    ** Initialize the temporary particles.
    */
    for (i=0;i<pkd->nLocal;++i) {
	PARTICLE *p = &pkd->pStore[i];
	pl[i].r[0] = pkdGetDistance(pkd,p,dCenter);
	pl[i].r[1] = pkdMass(pkd,p);
	pl[i].r[2] = 0.0;
	pl[i].i = i;
	}
    qsort(pkd->pLite,pkdLocal(pkd),sizeof(PLITE),cmpRadiusLite);
    }

/*
** Density Profile
*/
void pkdProfile(PKD pkd, uint8_t uRungLo, uint8_t uRungHi,
		double *dCenter, double dMinRadius, double dMaxRadius,
		int nBins) {
    PARTICLE *p;
    double d2, r2, r, r1;
    double deltaR;
    double dLogMinRadius, dLogMaxRadius, dLogRadius;
    int i,n,iBin;
    PROFILEBIN *pBin;

    assert(dMinRadius>0.0);
    assert(dMaxRadius>dMinRadius);

    assert( pkd->profileBins == NULL );
    pkd->profileBins = mdlMalloc(pkd->mdl,nBins*sizeof(PROFILEBIN));

    r2 = dMaxRadius*dMaxRadius;
    dLogMinRadius = log10(dMinRadius);
    dLogMaxRadius = log10(dMaxRadius);
    dLogRadius = dLogMaxRadius - dLogMinRadius;
    deltaR = dLogRadius / nBins;
    r = 0.0;
    for( iBin=0; iBin<nBins; iBin++ ) {
	pkd->profileBins[iBin].nParticles = 0;
	pkd->profileBins[iBin].dMassInBin = 0.0;
	r1 = pow(10,dLogMinRadius+deltaR*(iBin+1));
	pkd->profileBins[iBin].dRadius = r1;
	pkd->profileBins[iBin].dVolume = (4.0/3.0) * M_PI * (r1*r1*r1 - r*r*r);
	r = r1;
	}

    n = pkdLocal(pkd);
    for (i=0;i<n;++i) {
	p = &pkd->pStore[i];
	if (pkdIsSrcActive(p,uRungLo,uRungHi)) {
	    d2 = pkdGetDistance(pkd,p,dCenter);
	    if ( d2 < dMinRadius*dMinRadius ) d2 = dMinRadius*dMinRadius;
	    if ( d2 < r2 ) {
		iBin = (log10(sqrt(d2))-dLogMinRadius)/dLogRadius * nBins;
		if ( iBin >= nBins ) iBin = nBins-1;
		assert(iBin>=0);
		pkd->profileBins[iBin].dMassInBin += pkdMass(pkd,p);
		pkd->profileBins[iBin].nParticles++;
		}
	    }
	}

    /* Combine the work from all processors */
    mdlCOcache(pkd->mdl,CID_BIN,pkd->profileBins,sizeof(PROFILEBIN),nBins,initProfileBins,combProfileBins);
    if (pkd->idSelf != 0) {
	for (i=0; i<nBins; i++) {
	    if (pkd->profileBins[i].dMassInBin > 0.0) {
		pBin = mdlAquire(pkd->mdl,CID_BIN,i,0);
		*pBin = pkd->profileBins[i];
		mdlRelease(pkd->mdl,CID_BIN,pBin);
		}
	    }
	}
    mdlFinishCache(pkd->mdl,CID_BIN);

    /* Only the main processor needs the result */
    if (pkd->idSelf != 0) {
	mdlFree(pkd->mdl,pkd->profileBins);
	pkd->profileBins = NULL;
	}
    }
