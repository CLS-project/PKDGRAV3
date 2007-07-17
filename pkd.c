#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#define _LARGEFILE_SOURCE 
#define _FILE_OFFSET_BITS 64 
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <math.h>
#include <assert.h>
#include <errno.h>
#include <sys/time.h>
#include <time.h> /* added MZ */

#include <rpc/types.h>
#include <rpc/xdr.h>

#include "pkd.h"
#include "ewald.h"
#include "walk.h"
#include "grav.h"
#include "mdl.h"
#include "tipsydefs.h"
#include "ssio.h"

#include "parameters.h"
#include "cosmo.h"

#ifdef USE_HDF5
#include "iohdf5.h"
#endif

#ifdef BSC
#include "mpitrace_user_events.h"
#endif

double Zeit() { /* added MZ */
    struct timeval tv;

    gettimeofday(&tv,NULL);
    return (tv.tv_sec+(tv.tv_usec*1e-6));
    }

double pkdGetTimer(PKD pkd,int iTimer)
    {
    return(pkd->ti[iTimer].sec);
    }

double pkdGetSystemTimer(PKD pkd,int iTimer)
    {
    return(pkd->ti[iTimer].system_sec);
    }

double pkdGetWallClockTimer(PKD pkd,int iTimer)
    {
    return(pkd->ti[iTimer].wallclock_sec);
    }


void pkdClearTimer(PKD pkd,int iTimer)
    {
    int i;

    if (iTimer >= 0) {
	pkd->ti[iTimer].sec = 0.0;
	pkd->ti[iTimer].system_sec = 0.0;
	pkd->ti[iTimer].wallclock_sec = 0.0;
	pkd->ti[iTimer].iActive = 0;
	}
    else {
	for (i=0;i<MAX_TIMERS;++i) {
	    pkd->ti[i].sec = 0.0;
	    pkd->ti[i].system_sec = 0.0;
	    pkd->ti[i].wallclock_sec = 0.0;
	    pkd->ti[i].iActive = 0;
	    }
	}
    }


void pkdStartTimer(PKD pkd,int iTimer)
    {
    struct timeval tv;

    pkd->ti[iTimer].iActive++;

    if (pkd->ti[iTimer].iActive == 1) {
	pkd->ti[iTimer].stamp = mdlCpuTimer(pkd->mdl);
	gettimeofday(&tv,NULL);
	pkd->ti[iTimer].wallclock_stamp = tv.tv_sec + 1e-6*(double) tv.tv_usec;
	    {
	    struct rusage ru;
	    
	    getrusage(0,&ru);
	    pkd->ti[iTimer].system_stamp = (double)ru.ru_stime.tv_sec + 1e-6*(double)ru.ru_stime.tv_usec;
	    }
	}
    }


void pkdStopTimer(PKD pkd,int iTimer)
    {
    double sec;
    struct timeval tv;

    sec = -pkd->ti[iTimer].stamp;
    pkd->ti[iTimer].stamp = mdlCpuTimer(pkd->mdl);
    sec += pkd->ti[iTimer].stamp;
    if (sec < 0.0) sec = 0.0;
    pkd->ti[iTimer].sec += sec;

    sec = -pkd->ti[iTimer].wallclock_stamp;
    gettimeofday( &tv, NULL );
    pkd->ti[iTimer].wallclock_stamp = tv.tv_sec + 1e-6*(double)tv.tv_usec;
    sec += pkd->ti[iTimer].wallclock_stamp;
    if (sec < 0.0) sec = 0.0;
    pkd->ti[iTimer].wallclock_sec += sec;

#ifndef _CRAYMPP
	{
	struct rusage ru;

	sec = -pkd->ti[iTimer].system_stamp;
	getrusage(0,&ru);
	pkd->ti[iTimer].system_stamp = ((double)ru.ru_stime.tv_sec + 1e-6*(double)ru.ru_stime.tv_usec);
	sec += pkd->ti[iTimer].system_stamp;
	if (sec < 0.0) sec = 0.0;
	pkd->ti[iTimer].system_sec += sec;
	}
#endif
	pkd->ti[iTimer].iActive--;
    }


void pkdInitialize(PKD *ppkd,MDL mdl,int nStore,FLOAT *fPeriod,
		   int nDark,int nGas,int nStar,double dSunMass)
    {
    PKD pkd;
    int j;
	
    pkd = (PKD)malloc(sizeof(struct pkdContext));
    mdlassert(mdl,pkd != NULL);
    pkd->mdl = mdl;
    pkd->idSelf = mdlSelf(mdl);
    pkd->nThreads = mdlThreads(mdl);
    pkd->nStore = nStore;
    pkd->nLocal = 0;
    pkd->nDark = nDark;
    pkd->nGas = nGas;
    pkd->nStar = nStar;
    pkd->nMaxOrderGas = nGas - 1;
    pkd->nMaxOrderDark = nGas + nDark - 1;
    pkd->nRejects = 0;
#ifdef PLANETS
    pkd->dSunMass = dSunMass;
#endif
    for (j=0;j<3;++j) {
	pkd->fPeriod[j] = fPeriod[j];
	}

    pkd->uMinRungActive  = 0;
    pkd->uMaxRungActive  = 255;
    pkd->uRungVeryActive = 255;

    /*
    ** Allocate the main particle store.
    ** Need to use mdlMalloc() since the particles will need to be
    ** visible to all other processors thru mdlAquire() later on.
    **
    ** We need one EXTRA storage location at the very end to use for 
    ** calculating acceleration on arbitrary positions in space, for example
    ** determining the force on the sun. The easiest way to do this is to
    ** allocate one hidden particle, which won't interfere with the rest of
    ** the code (hopefully). pkd->pStore[pkd->nStore] is this particle.
    */
    pkd->pStore = mdlMalloc(pkd->mdl,(nStore+1)*sizeof(PARTICLE));
    mdlassert(mdl,pkd->pStore != NULL);
    pkd->pLite = malloc((nStore+1)*sizeof(PLITE));
    mdlassert(mdl,pkd->pLite != NULL);
    pkd->nNodes = 0;
    pkd->kdNodes = NULL;
    pkd->kdTop = NULL;
    /*
    ** Allocate initial temporary node storage.
    */
    pkd->nMaxNodes = 10000;
    pkd->kdTemp = malloc(pkd->nMaxNodes*sizeof(KDT));
    mdlassert(mdl,pkd->kdTemp != NULL);
    /*
    ** Ewald stuff!
    */
    pkd->nMaxEwhLoop = 100;
    pkd->ewt = malloc(pkd->nMaxEwhLoop*sizeof(EWT));
    mdlassert(mdl,pkd->ewt != NULL);
    *ppkd = pkd;
    /*
    ** Allocate initial particle pointer arrays for active/inactive particles.
    */
    pkd->nMaxBucketActive = 1000;
    pkd->piActive = malloc(pkd->nMaxBucketActive*sizeof(PARTICLE *));
    mdlassert(mdl,pkd->piActive != NULL);
    pkd->piInactive = malloc(pkd->nMaxBucketActive*sizeof(PARTICLE *));
    mdlassert(mdl,pkd->piInactive != NULL);
    }


void pkdFinish(PKD pkd)
    {
    if (pkd->kdNodes) {
	/*
	** Close caching space and free up nodes.
	*/
	mdlFinishCache(pkd->mdl,CID_CELL);
	mdlFree(pkd->mdl,pkd->kdNodes);
	}
    if (pkd->kdTop) free(pkd->kdTop);
    free(pkd->ewt);
    mdlFree(pkd->mdl,pkd->pStore);
    free(pkd->pLite);
    csmFinish(pkd->param.csm);
    free(pkd);
    }


void pkdSeek(PKD pkd,FILE *fp,int nStart,int bStandard,int bDoublePos) {
    off_t MAX_OFFSET = 2147483640;
    long long nStart64 = nStart;
    long long nGas64 = pkd->nGas;
    long long nDark64 = pkd->nDark;
    off_t lStart;
    int iErr;
    /*
    ** Seek according to true XDR size structures when bStandard is true.
    ** This may be a bit dicey, but it should work as long
    ** as no one changes the tipsy binary format!
    */
    if (bStandard) lStart = 32;
    else lStart = sizeof(struct dump);
    if (nStart64 > nGas64) {
	if (bStandard) lStart += nGas64*(bDoublePos?60:48);
	else lStart += nGas64*sizeof(struct gas_particle);
	nStart64 -= nGas64;
	if (nStart64> nDark64) {
	    if (bStandard) lStart += nDark64*(bDoublePos?48:36);
	    else lStart += nDark64*sizeof(struct dark_particle);
	    nStart64 -= nDark64;
	    if (bStandard) lStart += nStart64*(bDoublePos?56:44);
	    else lStart += nStart64*sizeof(struct star_particle);
	    }
	else {
	    if (bStandard) lStart += nStart64*(bDoublePos?48:36);
	    else lStart += nStart64*sizeof(struct dark_particle);
	    }
	} 
    else {
	if (bStandard) lStart += nStart64*(bDoublePos?60:48);
	else lStart += nStart64*sizeof(struct gas_particle);
	}
    
    /*fseek fails for offsets >= 2**31; this is an ugly workaround;*/
    if(lStart > MAX_OFFSET){
	iErr = fseek(fp,0,SEEK_SET);
	if (iErr) {
	    perror("pkdSeek failed");
	    exit(errno);
	    }
	while(lStart > MAX_OFFSET){
	    fseek(fp,MAX_OFFSET,SEEK_CUR);
	    lStart -= MAX_OFFSET;
	    }
	iErr = fseek(fp,lStart,SEEK_CUR);
	if (iErr) {
	    perror("pkdSeek failed");
	    exit(errno);
	    }
	} 
    else { 
	iErr = fseek(fp,lStart,SEEK_SET);
	if (iErr) {
	    perror("pkdSeek failed");
	    exit(errno);
	    }
	}
    }


void IOCheck(int nout) {
    if (nout != 1) {
	perror("IOCheck failed");
	exit(errno);
	}
    }


#ifdef USE_HDF5
void pkdReadHDF5(PKD pkd, IOHDF5 io, double dvFac,double duTFac, int nStart, int nLocal ) {
    PARTICLE *p;
    FLOAT dT1, dT2;
    int i, j;

    pkd->nLocal = nLocal;
    pkd->nActive = nLocal;

    /*
    ** General initialization.
    */
    for (i=0;i<nLocal;++i) {
	p = &pkd->pStore[i];
	TYPEClear(p);
	p->iRung = 0;
	p->fWeight = 1.0;
	p->fDensity = 0.0;
	p->fBall = 0.0;
	}

    /*TODO: add tracker file */

    for (i=0;i<nLocal;++i) {
	p = &pkd->pStore[i];
	p->iOrder = nStart + i;

	if (pkdIsDark(pkd,p)) {
	    ioHDF5GetDark( io, &p->iOrder, p->r, p->v,
			   &p->fMass, &p->fSoft, &p->fPot );
	    for (j=0;j<3;++j) p->v[j] *= dvFac;
#ifdef CHANGESOFT
	    p->fSoft0 = p->fSoft;
#endif
	}
	else if (pkdIsGas(pkd,p)) {
	    ioHDF5GetGas( io, &p->iOrder, p->r, p->v,
			  &p->fMass, &p->fSoft, &p->fPot,
			  &dT1, &dT2 );
	    for (j=0;j<3;++j) p->v[j] *= dvFac;
#ifdef CHANGESOFT
	    p->fSoft0 = p->fSoft;
#endif
	}
	else if (pkdIsStar(pkd,p)) {
	    ioHDF5GetStar( io, &p->iOrder, p->r, p->v,
			   &p->fMass, &p->fSoft, &p->fPot,
			   &dT1, &dT2 );
	    for (j=0;j<3;++j) p->v[j] *= dvFac;
#ifdef CHANGESOFT
	    p->fSoft0 = p->fSoft;
#endif
	}
	else mdlassert(pkd->mdl,0);
	if(p->fSoft < sqrt(2.0e-38)) { /* set minimum softening */
	    p->fSoft = sqrt(2.0e-38);
	}
    }
}
#endif


void pkdReadTipsy(PKD pkd,char *pszFileName, char *achOutName, int nStart,int nLocal,
		  int bStandard,double dvFac,double dTuFac,int bDoublePos)
    {
    FILE *fp;
    int i,j;
    PARTICLE *p;
    struct dark_particle dp;
    struct gas_particle gp;
    struct star_particle sp;
    float fTmp;
    double dTmp;
    float mass=0.0;
    /*
    ** Variables for particle tracking
    */
    FILE *fpTrack = NULL;
    char aTrack[256];
    int nBodies,nGas,nStar,iRet;
    int iTracker;	

    pkd->nLocal = nLocal;
    pkd->nActive = nLocal;
    /*
    ** General initialization.
    */
    for (i=0;i<nLocal;++i) {
	p = &pkd->pStore[i];
	TYPEClear(p);
	p->iRung = 0;
	p->fWeight = 1.0;
	p->fDensity = 0.0;
	p->fBall = 0.0;
	}
    /*
    ** Seek past the header and up to nStart.
    */
    fp = fopen(pszFileName,"r");
    mdlassert(pkd->mdl,fp != NULL);
    /*
    ** Seek to right place in file.
    */
    pkdSeek(pkd,fp,nStart,bStandard,bDoublePos);
    /*
    ** See if the user has specified a .track file.
    */
    iTracker = nStart-1;
    sprintf(aTrack,"%s.track",achOutName);
    fpTrack = fopen(aTrack,"r");
    if (fpTrack) {
	/*
	** Get to the right place in the file.
	*/
	iRet = fscanf(fpTrack,"%d %d %d",&nBodies,&nGas,&nStar);
	if (!iRet || iRet == EOF) {
	    fclose(fpTrack);
	    goto SkipCheck;
	    }
	while (1) {
	    iRet = fscanf(fpTrack,"%d",&iTracker);
	    iTracker = iTracker - 1;
	    if (!iRet || iRet == EOF) {
		fclose(fpTrack);
		iTracker = nStart-1;
		break;
		}
	    if (iTracker >= nStart) break;
	    }
	}
SkipCheck:
    /*
    ** Read Stuff!
    */
    if (bStandard) {
	FLOAT vTemp;
	XDR xdrs;
	xdrstdio_create(&xdrs,fp,XDR_DECODE);
	for (i=0;i<nLocal;++i) {
	    p = &pkd->pStore[i];
	    p->iOrder = nStart + i;
	    if (p->iOrder == iTracker) {
		TYPESet(p,TYPE_TRACKER);
		iRet = fscanf(fpTrack,"%d",&iTracker);
		iTracker = iTracker-1;
		if (!iRet || iRet == EOF) {
		    fclose(fpTrack);
		    iTracker = nStart-1;
		    }
		}
	    if (pkdIsDark(pkd,p)) {
		xdr_float(&xdrs,&fTmp);
		p->fMass = fTmp;
		mass += fTmp;
		if (bDoublePos) {
		    for (j=0;j<3;++j) {
			xdr_double(&xdrs,&dTmp);
			p->r[j] = dTmp;
			}
		    }
		else {
		    for (j=0;j<3;++j) {
			xdr_float(&xdrs,&fTmp);
			p->r[j] = fTmp;
			}
		    }
		for (j=0;j<3;++j) {
		    xdr_float(&xdrs,&fTmp);
		    vTemp = fTmp;
		    p->v[j] = dvFac*vTemp;			
		    }
		xdr_float(&xdrs,&fTmp);
		p->fSoft = fTmp;
#ifdef CHANGESOFT				
		p->fSoft0 = fTmp;
#endif
		xdr_float(&xdrs,&fTmp);
		p->fPot = fTmp;
		}
	    else if (pkdIsGas(pkd,p)) {
		xdr_float(&xdrs,&fTmp);
		p->fMass = fTmp;
		if (bDoublePos) {
		    for (j=0;j<3;++j) {
			xdr_double(&xdrs,&dTmp);
			p->r[j] = dTmp;
			}
		    }
		else {
		    for (j=0;j<3;++j) {
			xdr_float(&xdrs,&fTmp);
			p->r[j] = fTmp;
			}
		    }
		for (j=0;j<3;++j) {
		    xdr_float(&xdrs,&fTmp);
		    vTemp = fTmp;
		    p->v[j] = dvFac*vTemp;			
		    }

		xdr_float(&xdrs,&fTmp);
		xdr_float(&xdrs,&fTmp);
		xdr_float(&xdrs,&fTmp);
		p->fSoft = fTmp;
#ifdef CHANGESOFT
		p->fSoft0 = fTmp;
#endif
		xdr_float(&xdrs,&fTmp);
		xdr_float(&xdrs,&fTmp);
		p->fPot = fTmp;
		}
	    else if (pkdIsStar(pkd,p)) {
		xdr_float(&xdrs,&fTmp);
		p->fMass = fTmp;
		if (bDoublePos) {
		    for (j=0;j<3;++j) {
			xdr_double(&xdrs,&dTmp);
			p->r[j] = dTmp;
			}
		    }
		else {
		    for (j=0;j<3;++j) {
			xdr_float(&xdrs,&fTmp);
			p->r[j] = fTmp;
			}
		    }
		for (j=0;j<3;++j) {
		    xdr_float(&xdrs,&fTmp);
		    vTemp = fTmp;
		    p->v[j] = dvFac*vTemp;			
		    }
		xdr_float(&xdrs,&fTmp);
		xdr_float(&xdrs,&fTmp);
		xdr_float(&xdrs,&fTmp);
		p->fSoft = fTmp;
#ifdef CHANGESOFT
		p->fSoft0 = fTmp;
#endif
		xdr_float(&xdrs,&fTmp);
		p->fPot = fTmp;
		}
	    else mdlassert(pkd->mdl,0);
	    if(p->fSoft < sqrt(2.0e-38)) { /* set minimum softening */
		p->fSoft = sqrt(2.0e-38);
		}
	    }
	xdr_destroy(&xdrs);
	}
    else {
	for (i=0;i<nLocal;++i) {
	    p = &pkd->pStore[i];
	    p->iOrder = nStart + i;
	    if (p->iOrder == iTracker) {
		TYPESet(p,TYPE_TRACKER);
		iRet = fscanf(fpTrack,"%d",&iTracker);
		iTracker = iTracker -1;
		if (!iRet || iRet == EOF) {
		    fclose(fpTrack);
		    iTracker = nStart-1;
		    }
		}
	    if (pkdIsDark(pkd,p)) {
		fread(&dp,sizeof(struct dark_particle),1,fp);
		for (j=0;j<3;++j) {
		    p->r[j] = dp.pos[j];
		    p->v[j] = dvFac*dp.vel[j];
		    }
		p->fMass = dp.mass;
		mass += dp.mass;
		p->fSoft = dp.eps;
#ifdef CHANGESOFT
		p->fSoft0 = dp.eps;
#endif
		p->fPot = dp.phi;
		}
	    else if (pkdIsGas(pkd,p)) {
		fread(&gp,sizeof(struct gas_particle),1,fp);
		for (j=0;j<3;++j) {
		    p->r[j] = gp.pos[j];
		    p->v[j] = dvFac*gp.vel[j];
		    }
		p->fMass = gp.mass;
		p->fSoft = gp.hsmooth;
#ifdef CHANGESOFT
		p->fSoft0 = gp.hsmooth;
#endif
		p->fPot = gp.phi;
		}
	    else if (pkdIsStar(pkd,p)) {
		fread(&sp,sizeof(struct star_particle),1,fp);
		for (j=0;j<3;++j) {
		    p->r[j] = sp.pos[j];
		    p->v[j] = dvFac*sp.vel[j];
		    }
		p->fMass = sp.mass;
		p->fSoft = sp.eps;
#ifdef CHANGESOFT
		p->fSoft0 = sp.eps;
#endif
		p->fPot = sp.phi;
		}
	    else mdlassert(pkd->mdl,0);
	    if(p->fSoft < sqrt(2.0e-38)) { /* set minimum softening */
		p->fSoft = sqrt(2.0e-38);
		}
	    }
	}
    fclose(fp);
    }


void pkdCalcBound(PKD pkd,BND *pbnd)
    {
    PARTICLE *p = pkd->pStore;
    FLOAT fMin[3],fMax[3];
    int i = 0;
    int j;

    mdlassert(pkd->mdl,pkd->nLocal > 0);
    for (j=0;j<3;++j) {
	fMin[j] = p[i].r[j];
	fMax[j] = p[i].r[j];
	}
    for (++i;i<pkd->nLocal;++i) {
	for (j=0;j<3;++j) {
	    if (p[i].r[j] < fMin[j]) fMin[j] = p[i].r[j];
	    else if (p[i].r[j] > fMax[j]) fMax[j] = p[i].r[j];
	    }
	}
    for (j=0;j<3;++j) {
	pbnd->fCenter[j] = 0.5*(fMax[j] + fMin[j]);
	pbnd->fMax[j] = 0.5*(fMax[j] - fMin[j]);
	}
    }


void pkdRungDDWeight(PKD pkd, int iMaxRung, double dWeight)
    {
    PARTICLE *p;
    int i;
    float fRungWeight[50],sum;

    mdlassert(pkd->mdl,iMaxRung<50);
    fRungWeight[0]=1.0;
    sum=1.0;
    for (i=1;i<=iMaxRung;i++) {
	sum*=2.0;
	fRungWeight[i] = dWeight* sum + (1-dWeight);
	}
  
    for(i=0;i<pkdLocal(pkd);++i) {
	p = &pkd->pStore[i];
	p->fWeight *= fRungWeight[p->iRung];
	}
    }

/*
** Partition particles between iFrom and iTo into those < fSplit and
** those >= to fSplit.  Find number and weight in each partition.
*/
int pkdWeight(PKD pkd,int d,FLOAT fSplit,int iSplitSide,int iFrom,int iTo,
	      int *pnLow,int *pnHigh,FLOAT *pfLow,FLOAT *pfHigh)
    {
    int i,iPart;
    FLOAT fLower,fUpper;

    /*
    ** First partition the memory about fSplit for particles iFrom to iTo.
    */
    if (iSplitSide) {
	iPart = pkdLowerPart(pkd,d,fSplit,iFrom,iTo);
	*pnLow = pkdLocal(pkd)-iPart;
	*pnHigh = iPart;
	}
    else {
	iPart = pkdUpperPart(pkd,d,fSplit,iFrom,iTo);
	*pnLow = iPart;
	*pnHigh = pkdLocal(pkd)-iPart;
	}
    /*
    ** Calculate the lower weight and upper weight BETWEEN the particles
    ** iFrom to iTo!
    */
    fLower = 0.0;
    for (i=iFrom;i<iPart;++i) {
	fLower += pkd->pStore[i].fWeight;
	}
    fUpper = 0.0;
    for (i=iPart;i<=iTo;++i) {
	fUpper += pkd->pStore[i].fWeight;
	}
    if (iSplitSide) {
	*pfLow = fUpper;
	*pfHigh = fLower;
	}
    else {
	*pfLow = fLower;
	*pfHigh = fUpper;
	}
    return(iPart);
    }


void pkdCountVA(PKD pkd,int d,FLOAT fSplit,int *pnLow,int *pnHigh) {
    int i;

    *pnLow = 0;
    *pnHigh = 0;
    for (i=0;i<pkd->nLocal;++i) {
	if (pkdIsVeryActive(pkd,pkd->pStore+i)) {
	    if (pkd->pStore[i].r[d] < fSplit) *pnLow += 1;
	    else *pnHigh += 1;
	    }
	}
    }


/*
** Partition particles between iFrom and iTo into those < fSplit and
** those >= to fSplit.  Find number and weight in each partition.
*/
int pkdWeightWrap(PKD pkd,int d,FLOAT fSplit,FLOAT fSplit2,int iSplitSide,int iVASplitSide,
		  int iFrom,int iTo,int *pnLow,int *pnHigh) {
    int iPart;

    /*
    ** First partition the memory about fSplit for particles iFrom to iTo.
    */
    if (!iSplitSide) {
	iPart = pkdLowerPartWrap(pkd,d,fSplit,fSplit2,iVASplitSide,iFrom,iTo);
	*pnLow = iPart;
	*pnHigh = pkdLocal(pkd)-iPart;
	}
    else {
	iPart = pkdUpperPartWrap(pkd,d,fSplit,fSplit2,iVASplitSide,iFrom,iTo);
	*pnHigh = iPart;
	*pnLow = pkdLocal(pkd)-iPart;
	}
    return(iPart);
    }


int pkdOrdWeight(PKD pkd,int iOrdSplit,int iSplitSide,int iFrom,int iTo,
		 int *pnLow,int *pnHigh)
    {
    int iPart;
	
    /*
    ** First partition the memory about fSplit for particles iFrom to iTo.
    */
    if (iSplitSide) {
	iPart = pkdLowerOrdPart(pkd,iOrdSplit,iFrom,iTo);
	*pnLow = pkdLocal(pkd)-iPart;
	*pnHigh = iPart;
	}
    else {
	iPart = pkdUpperOrdPart(pkd,iOrdSplit,iFrom,iTo);
	*pnLow = iPart;
	*pnHigh = pkdLocal(pkd)-iPart;
	}
    return(iPart);
    }


int pkdLowerPart(PKD pkd,int d,FLOAT fSplit,int i,int j)
    {
    PARTICLE pTemp;

    PARTITION(pkd->pStore,pTemp,i,j,
	      pkd->pStore[i].r[d] >= fSplit,
	      pkd->pStore[j].r[d] < fSplit);
    return(i);
    }


int pkdUpperPart(PKD pkd,int d,FLOAT fSplit,int i,int j)
    {
    PARTICLE pTemp;

    PARTITION(pkd->pStore,pTemp,i,j,
	      pkd->pStore[i].r[d] < fSplit,
	      pkd->pStore[j].r[d] >= fSplit);
    return(i);
    }


int pkdLowerPartWrap(PKD pkd,int d,FLOAT fSplit1,FLOAT fSplit2,int iVASplitSide,int i,int j) {
    PARTICLE pTemp;

    if (fSplit1 > fSplit2) {
	if (iVASplitSide < 0) {
	    PARTITION(pkd->pStore,pTemp,i,j,
		      (pkd->pStore[i].r[d] < fSplit2 || pkd->pStore[i].r[d] >= fSplit1) &&
		      !pkdIsVeryActive(pkd,pkd->pStore+i),
		      (pkd->pStore[j].r[d] >= fSplit2 && pkd->pStore[j].r[d] < fSplit1) ||
		      pkdIsVeryActive(pkd,pkd->pStore+j));
	    }
	else if (iVASplitSide > 0) {
	    PARTITION(pkd->pStore,pTemp,i,j,
		      (pkd->pStore[i].r[d] < fSplit2 || pkd->pStore[i].r[d] >= fSplit1) ||
		      pkdIsVeryActive(pkd,pkd->pStore+i),
		      (pkd->pStore[j].r[d] >= fSplit2 && pkd->pStore[j].r[d] < fSplit1) &&
		      !pkdIsVeryActive(pkd,pkd->pStore+j));
	    }
	else {
	    PARTITION(pkd->pStore,pTemp,i,j,
		      (pkd->pStore[i].r[d] < fSplit2 || pkd->pStore[i].r[d] >= fSplit1),
		      (pkd->pStore[j].r[d] >= fSplit2 && pkd->pStore[j].r[d] < fSplit1));
	    }
	}
    else {
	if (iVASplitSide < 0) {
	    PARTITION(pkd->pStore,pTemp,i,j,
		      (pkd->pStore[i].r[d] < fSplit2 && pkd->pStore[i].r[d] >= fSplit1) &&
		      !pkdIsVeryActive(pkd,pkd->pStore+i),
		      (pkd->pStore[j].r[d] >= fSplit2 || pkd->pStore[j].r[d] < fSplit1) ||
		      pkdIsVeryActive(pkd,pkd->pStore+j));
	    }
	else if (iVASplitSide > 0) {
	    PARTITION(pkd->pStore,pTemp,i,j,
		      (pkd->pStore[i].r[d] < fSplit2 && pkd->pStore[i].r[d] >= fSplit1) ||
		      pkdIsVeryActive(pkd,pkd->pStore+i),
		      (pkd->pStore[j].r[d] >= fSplit2 || pkd->pStore[j].r[d] < fSplit1) &&
		      !pkdIsVeryActive(pkd,pkd->pStore+j));
	    }
	else {
	    PARTITION(pkd->pStore,pTemp,i,j,
		      (pkd->pStore[i].r[d] < fSplit2 && pkd->pStore[i].r[d] >= fSplit1),
		      (pkd->pStore[j].r[d] >= fSplit2 || pkd->pStore[j].r[d] < fSplit1));
	    }
	}
    return(i);
    }


int pkdUpperPartWrap(PKD pkd,int d,FLOAT fSplit1,FLOAT fSplit2,int iVASplitSide,int i,int j) {
    PARTICLE pTemp;

    if (fSplit1 > fSplit2) {
	if (iVASplitSide < 0) {
	    PARTITION(pkd->pStore,pTemp,i,j,
		      (pkd->pStore[i].r[d] >= fSplit2 && pkd->pStore[i].r[d] < fSplit1) ||
		      pkdIsVeryActive(pkd,pkd->pStore+i),
		      (pkd->pStore[j].r[d] < fSplit2 || pkd->pStore[j].r[d] >= fSplit1) &&
		      !pkdIsVeryActive(pkd,pkd->pStore+j));
	    }
	else if (iVASplitSide > 0) {
	    PARTITION(pkd->pStore,pTemp,i,j,
		      (pkd->pStore[i].r[d] >= fSplit2 && pkd->pStore[i].r[d] < fSplit1) &&
		      !pkdIsVeryActive(pkd,pkd->pStore+i),
		      (pkd->pStore[j].r[d] < fSplit2 || pkd->pStore[j].r[d] >= fSplit1) ||
		      pkdIsVeryActive(pkd,pkd->pStore+j));
	    }
	else {
	    PARTITION(pkd->pStore,pTemp,i,j,
		      (pkd->pStore[i].r[d] >= fSplit2 && pkd->pStore[i].r[d] < fSplit1),
		      (pkd->pStore[j].r[d] < fSplit2 || pkd->pStore[j].r[d] >= fSplit1));
	    }
	}
    else {
	if (iVASplitSide < 0) {
	    PARTITION(pkd->pStore,pTemp,i,j,
		      (pkd->pStore[i].r[d] >= fSplit2 || pkd->pStore[i].r[d] < fSplit1) ||
		      pkdIsVeryActive(pkd,pkd->pStore+i),
		      (pkd->pStore[j].r[d] < fSplit2 && pkd->pStore[j].r[d] >= fSplit1) &&
		      !pkdIsVeryActive(pkd,pkd->pStore+j));
	    }
	else if (iVASplitSide > 0) {
	    PARTITION(pkd->pStore,pTemp,i,j,
		      (pkd->pStore[i].r[d] >= fSplit2 || pkd->pStore[i].r[d] < fSplit1) &&
		      !pkdIsVeryActive(pkd,pkd->pStore+i),
		      (pkd->pStore[j].r[d] < fSplit2 && pkd->pStore[j].r[d] >= fSplit1) ||
		      pkdIsVeryActive(pkd,pkd->pStore+j));
	    }
	else {
	    PARTITION(pkd->pStore,pTemp,i,j,
		      (pkd->pStore[i].r[d] >= fSplit2 || pkd->pStore[i].r[d] < fSplit1),
		      (pkd->pStore[j].r[d] < fSplit2 && pkd->pStore[j].r[d] >= fSplit1));
	    }
	}
    return(i);
    }


int pkdLowerOrdPart(PKD pkd,int nOrdSplit,int i,int j) {
    PARTICLE pTemp;

    PARTITION(pkd->pStore,pTemp,i,j,
	      pkd->pStore[i].iOrder >= nOrdSplit,
	      pkd->pStore[j].iOrder < nOrdSplit);
    return(i);
    }


int pkdUpperOrdPart(PKD pkd,int nOrdSplit,int i,int j) {
    PARTICLE pTemp;

    PARTITION(pkd->pStore,pTemp,i,j,
	      pkd->pStore[i].iOrder < nOrdSplit,
	      pkd->pStore[j].iOrder >= nOrdSplit);
    return(i);
    }


int pkdActiveOrder(PKD pkd) {
    PARTICLE pTemp;
    int i=0;
    int j=pkdLocal(pkd)-1;

    PARTITION(pkd->pStore,pTemp,i,j,
	      pkdIsActive(pkd,&(pkd->pStore[i])),
	      !pkdIsActive(pkd,&(pkd->pStore[j])));
    return (pkd->nActive = i);
    }


int pkdColRejects_Active_Inactive(PKD pkd,int d,FLOAT fSplit,FLOAT fSplitInactive,
				  int iSplitSide)
    {
    PARTICLE pTemp;
    int nSplit,nSplitInactive,iRejects,i,j;

    mdlassert(pkd->mdl,pkd->nRejects == 0);
    if (iSplitSide) {
	nSplit = pkdLowerPart(pkd,d,fSplit,0,pkdActive(pkd)-1);
	}
    else {
	nSplit = pkdUpperPart(pkd,d,fSplit,0,pkdActive(pkd)-1);
	}
    if (iSplitSide) {
	nSplitInactive = pkdLowerPart(pkd,d,fSplitInactive,
				      pkdActive(pkd),pkdLocal(pkd)-1);
	}
    else {
	nSplitInactive = pkdUpperPart(pkd,d,fSplitInactive,
				      pkdActive(pkd),pkdLocal(pkd)-1);
	}
    /*
      for(i = 0; i < nSplit; ++i)
      mdlassert(pkd->mdl,pkdIsActive(pkd,&(pkd->pStore[i])));
      for(i = pkdActive(pkd); i < nSplitInactive; ++i)
      mdlassert(pkd->mdl,!pkdIsActive(pkd,&(pkd->pStore[i])));
    */

    nSplitInactive -= pkdActive(pkd);
    /*
    ** Now do some fancy rearrangement.
    */
    i = nSplit;
    j = nSplit+nSplitInactive;
    while (j < pkdActive(pkd) + nSplitInactive) {
	pTemp = pkd->pStore[i];
	pkd->pStore[i] = pkd->pStore[j];
	pkd->pStore[j] = pTemp;
	++i; ++j;
	}
    pkd->nRejects = pkdLocal(pkd) - nSplit - nSplitInactive;
    iRejects = pkdFreeStore(pkd) - pkd->nRejects;
    pkd->nActive = nSplit;
    pkd->nLocal = nSplit + nSplitInactive;
    /*
    ** Move rejects to High memory.
    */
    for (i=pkd->nRejects-1;i>=0;--i)
	pkd->pStore[iRejects+i] = pkd->pStore[pkd->nLocal+i];
    return(pkd->nRejects);
    }


int pkdColRejects(PKD pkd,int nSplit)
    {
    int iRejects,i;

    mdlassert(pkd->mdl,pkd->nRejects == 0);

    pkd->nRejects = pkdLocal(pkd) - nSplit;
    iRejects = pkdFreeStore(pkd) - pkd->nRejects;
    pkd->nLocal = nSplit;
    /*
    ** Move rejects to High memory.
    */
    for (i=pkd->nRejects-1;i>=0;--i)
	pkd->pStore[iRejects+i] = pkd->pStore[pkd->nLocal+i];
    return(pkd->nRejects);
    }


int pkdSwapRejects(PKD pkd,int idSwap)
    {
    size_t nBuf;
    size_t nOutBytes,nSndBytes,nRcvBytes;

    if (idSwap != -1) {
	nBuf = (pkdSwapSpace(pkd))*sizeof(PARTICLE);
	nOutBytes = pkd->nRejects*sizeof(PARTICLE);
	mdlassert(pkd->mdl,pkdLocal(pkd) + pkd->nRejects <= pkdFreeStore(pkd));
	mdlSwap(pkd->mdl,idSwap,nBuf,&pkd->pStore[pkdLocal(pkd)],
		nOutBytes,&nSndBytes,&nRcvBytes);
	pkd->nLocal += nRcvBytes/sizeof(PARTICLE);
	pkd->nRejects -= nSndBytes/sizeof(PARTICLE);
	}
    return(pkd->nRejects);
    }

void pkdSwapAll(PKD pkd, int idSwap)
    {
    size_t nBuf;
    size_t nOutBytes,nSndBytes,nRcvBytes;
    int i;
    int iBuf;
    
    /*
    ** Move particles to High memory.
    */
    iBuf = pkdSwapSpace(pkd);
    for (i=pkdLocal(pkd)-1;i>=0;--i)
	pkd->pStore[iBuf+i] = pkd->pStore[i];

    nBuf = pkdFreeStore(pkd)*sizeof(PARTICLE);
    nOutBytes = pkdLocal(pkd)*sizeof(PARTICLE);
    mdlSwap(pkd->mdl,idSwap,nBuf,&pkd->pStore[0], nOutBytes,
	    &nSndBytes, &nRcvBytes);
    mdlassert(pkd->mdl,nSndBytes/sizeof(PARTICLE) == pkdLocal(pkd));
    pkd->nLocal = nRcvBytes/sizeof(PARTICLE);
    }

int pkdSwapSpace(PKD pkd)
    {
    return(pkdFreeStore(pkd) - pkdLocal(pkd));
    }


int pkdFreeStore(PKD pkd)
    {
    return(pkd->nStore);
    }

int pkdActive(PKD pkd)
    {
    return(pkd->nActive);
    }

int pkdInactive(PKD pkd)
    {
    return(pkd->nLocal - pkd->nActive);
    }

int pkdLocal(PKD pkd)
    {
    return(pkd->nLocal);
    }

int pkdNodes(PKD pkd)
    {
    return(pkd->nNodes);
    }


int pkdColOrdRejects(PKD pkd,int nOrdSplit,int iSplitSide)
    {
    int nSplit,iRejects,i;

    if (iSplitSide) nSplit = pkdLowerOrdPart(pkd,nOrdSplit,0,pkdLocal(pkd)-1);
    else nSplit = pkdUpperOrdPart(pkd,nOrdSplit,0,pkdLocal(pkd)-1);
    pkd->nRejects = pkdLocal(pkd) - nSplit;
    iRejects = pkdFreeStore(pkd) - pkd->nRejects;
    pkd->nLocal = nSplit;
    /*
    ** Move rejects to High memory.
    */
    for (i=pkd->nRejects-1;i>=0;--i)
	pkd->pStore[iRejects+i] = pkd->pStore[pkd->nLocal+i];
    return(pkd->nRejects);
    }


int cmpParticles(const void *pva,const void *pvb)
    {
    PARTICLE *pa = (PARTICLE *)pva;
    PARTICLE *pb = (PARTICLE *)pvb;

    return(pa->iOrder - pb->iOrder);
    }


void pkdLocalOrder(PKD pkd)
    {
    qsort(pkd->pStore,pkdLocal(pkd),sizeof(PARTICLE),cmpParticles);
    }

int pkdPackIO(PKD pkd,PIO *io,int nStart,int nMax)
{
    int nCopy, i, d;

    mdlassert(pkd->mdl,nStart<=pkd->nLocal);

    /* Calculate the number of particles to copy to the output buffer */
    nCopy = pkd->nLocal - nStart;
    if ( nCopy > nMax ) nCopy = nMax;

    for( i=0; i<nCopy; i++ ) {
	if ( pkdIsDark(pkd,pkd->pStore+nStart+i) ) {
	}
	for( d=0; d<3; d++ ) {
	    io[i].r[d] = pkd->pStore[nStart+i].r[d];
	    io[i].v[d] = pkd->pStore[nStart+i].v[d];
	}
	io[i].fMass = pkd->pStore[nStart+i].fMass;
    }

    return nCopy;
}

#ifdef USE_HDF5
void pkdWriteHDF5(PKD pkd, IOHDF5 io, double dvFac,double duTFac)
{
    PARTICLE *p;
    FLOAT v[3], fSoft;
    int i;

    for (i=0;i<pkdLocal(pkd);++i) {
	p = &pkd->pStore[i];

	v[0] = p->v[0] * dvFac;
	v[1] = p->v[1] * dvFac;
	v[2] = p->v[2] * dvFac;
#ifdef CHANGESOFT
	fSoft = p->fSoft0;
#else
	fSoft = p->fSoft;
#endif


	if (pkdIsDark(pkd,p)) {
	    ioHDF5AddDark(io,p->iOrder,p->r,v,
			  p->fMass,fSoft,p->fPot );
	}
	else if (pkdIsGas(pkd,p)) {
	    /* Why are temp and metals always set to zero? */
	    ioHDF5AddGas( io,p->iOrder,p->r,v,
			  p->fMass,fSoft,p->fPot,0.0,0.0);
	}
	else if (pkdIsStar(pkd,p)) {
	    /* Why are metals and tform always set to zero? */
	    ioHDF5AddStar(io, p->iOrder, p->r, v,
			  p->fMass,fSoft,p->fPot,0.0,0.0);
	}
	else mdlassert(pkd->mdl,0);
    }
}


#endif

void pkdWriteTipsy(PKD pkd,char *pszFileName,int nStart,
		   int bStandard,double dvFac,double duTFac,int bDoublePos) {
    PARTICLE *p;
    FILE *fp;
    int i,j;
    struct dark_particle dp;
    struct gas_particle gp;
    struct star_particle sp;
    int nout;
    float fTmp;
    double dTmp;
	
    /*
    ** Seek past the header and up to nStart.
    */
    fp = fopen(pszFileName,"r+");
    mdlassert(pkd->mdl,fp != NULL);
    pkdSeek(pkd,fp,nStart,bStandard,bDoublePos);

    if (bStandard) {
	FLOAT vTemp;
	XDR xdrs;
	/* 
	** Write Stuff!
	*/
	xdrstdio_create(&xdrs,fp,XDR_ENCODE);
	for (i=0;i<pkdLocal(pkd);++i) {
	    p = &pkd->pStore[i];
	    if (pkdIsDark(pkd,p)) {
		fTmp = p->fMass;
		assert(fTmp > 0.0);
		IOCheck(xdr_float(&xdrs,&fTmp));
		if (bDoublePos) {
		    for (j=0;j<3;++j) {
			dTmp = p->r[j];
			IOCheck(xdr_double(&xdrs,&dTmp));
			}
		    }
		else {
		    for (j=0;j<3;++j) {
			fTmp = p->r[j];
			IOCheck(xdr_float(&xdrs,&fTmp));
			}
		    }
		for (j=0;j<3;++j) {
		    vTemp = dvFac*p->v[j];			
		    fTmp = vTemp;
		    IOCheck(xdr_float(&xdrs,&fTmp));
		    }
#ifdef CHANGESOFT
		fTmp = p->fSoft0;
#else
		fTmp = p->fSoft;
#endif
		assert(fTmp > 0.0);
		IOCheck(xdr_float(&xdrs,&fTmp));
		fTmp = p->fPot;
		IOCheck(xdr_float(&xdrs,&fTmp));
		}
	    else if (pkdIsGas(pkd,p)) {
		fTmp = p->fMass;
		IOCheck(xdr_float(&xdrs,&fTmp));
		if (bDoublePos) {
		    for (j=0;j<3;++j) {
			dTmp = p->r[j];
			IOCheck(xdr_double(&xdrs,&dTmp));
			}
		    }
		else {
		    for (j=0;j<3;++j) {
			fTmp = p->r[j];
			IOCheck(xdr_float(&xdrs,&fTmp));
			}
		    }
		for (j=0;j<3;++j) {
		    vTemp = dvFac*p->v[j];			
		    fTmp = vTemp;
		    IOCheck(xdr_float(&xdrs,&fTmp));
		    }
		fTmp = p->fDensity;
		IOCheck(xdr_float(&xdrs,&fTmp));
		fTmp = 0.0;
		IOCheck(xdr_float(&xdrs,&fTmp));
#ifdef CHANGESOFT
		fTmp = p->fSoft0;
#else
		fTmp = p->fSoft;
#endif
		IOCheck(xdr_float(&xdrs,&fTmp));
		fTmp = 0.0;
		IOCheck(xdr_float(&xdrs,&fTmp));
		fTmp = p->fPot;
		IOCheck(xdr_float(&xdrs,&fTmp));
		}
	    else if (pkdIsStar(pkd,p)) {
		fTmp = p->fMass;
		IOCheck(xdr_float(&xdrs,&fTmp));
		if (bDoublePos) {
		    for (j=0;j<3;++j) {
			dTmp = p->r[j];
			IOCheck(xdr_double(&xdrs,&dTmp));
			}
		    }
		else {
		    for (j=0;j<3;++j) {
			fTmp = p->r[j];
			IOCheck(xdr_float(&xdrs,&fTmp));
			}
		    }
		for (j=0;j<3;++j) {
		    vTemp = dvFac*p->v[j];			
		    fTmp = vTemp;
		    IOCheck(xdr_float(&xdrs,&fTmp));
		    }
		fTmp = 0.0;
		IOCheck(xdr_float(&xdrs,&fTmp));
		IOCheck(xdr_float(&xdrs,&fTmp));
#ifdef CHANGESOFT
		fTmp = p->fSoft0;
#else
		fTmp = p->fSoft;
#endif
		IOCheck(xdr_float(&xdrs,&fTmp));
		fTmp = p->fPot;
		IOCheck(xdr_float(&xdrs,&fTmp));
		}
	    else mdlassert(pkd->mdl,0);
	    }
	xdr_destroy(&xdrs);
	}
    else {
	/* 
	** Write Stuff!
	*/
	for (i=0;i<pkdLocal(pkd);++i) {
	    p = &pkd->pStore[i];
	    if (pkdIsDark(pkd,p)) {
		for (j=0;j<3;++j) {
		    dp.pos[j] = p->r[j];
		    dp.vel[j] = dvFac*p->v[j];
		    }
		dp.mass = p->fMass;
#ifdef CHANGESOFT
		dp.eps = p->fSoft0;
#else
		dp.eps = p->fSoft;
#endif
		dp.phi = p->fPot;
		IOCheck(fwrite(&dp,sizeof(struct dark_particle),1,fp));
		}
	    else if (pkdIsGas(pkd,p)) {
		for (j=0;j<3;++j) {
		    gp.pos[j] = p->r[j];
		    gp.vel[j] = dvFac*p->v[j];
		    }
		gp.mass = p->fMass;
#ifdef CHANGESOFT
		gp.hsmooth = p->fSoft0;
#else
		gp.hsmooth = p->fSoft;
#endif
		gp.phi = p->fPot;
		gp.rho = p->fDensity;
		gp.temp = 0.0;
		gp.metals = 0.0;
		IOCheck(fwrite(&gp,sizeof(struct gas_particle),1,fp));
		}
	    else if (pkdIsStar(pkd,p)) {
		for (j=0;j<3;++j) {
		    sp.pos[j] = p->r[j];
		    sp.vel[j] = dvFac*p->v[j];
		    }
		sp.mass = p->fMass;
#ifdef CHANGESOFT
		sp.eps = p->fSoft0;
#else
		sp.eps = p->fSoft;
#endif
		sp.phi = p->fPot;
		sp.metals = 0.0;
		sp.tform = 0.0;
		IOCheck(fwrite(&sp,sizeof(struct star_particle),1,fp));
		}
	    else mdlassert(pkd->mdl,0);
	    }
	}
    nout = fclose(fp);
    mdlassert(pkd->mdl,nout == 0);
    }


void pkdSetSoft(PKD pkd,double dSoft)
    {
    PARTICLE *p;
    int i,n;

    p = pkd->pStore;
    n = pkdLocal(pkd);
    if(dSoft < sqrt(2.0e-38)) { /* set minimum softening */
	dSoft = sqrt(2.0e-38);
	}
    for (i=0;i<n;++i) {
#ifdef CHANGESOFT
	p[i].fSoft0 = dSoft;
#else
	p[i].fSoft = dSoft;
#endif
	}
    }

#ifdef CHANGESOFT
void pkdPhysicalSoft(PKD pkd,double dSoftMax,double dFac,int bSoftMaxMul)
    {
    PARTICLE *p;
    int i,n;

    p = pkd->pStore;
    n = pkdLocal(pkd);
	
    mdlassert(pkd->mdl,dFac > 0);
    if (bSoftMaxMul) {
	for (i=0;i<n;++i) {
	    mdlassert(pkd->mdl,p[i].fSoft0 > 0);
	    p[i].fSoft = p[i].fSoft0*dFac;
	    mdlassert(pkd->mdl,p[i].fSoft > 0);
	    }
	}
    else {
	mdlassert(pkd->mdl,dSoftMax > 0);
	for (i=0;i<n;++i) {
	    mdlassert(pkd->mdl,p[i].fSoft0 > 0);
	    p[i].fSoft = p[i].fSoft0*dFac;
	    if (p[i].fSoft > dSoftMax) p[i].fSoft = dSoftMax;
	    mdlassert(pkd->mdl,p[i].fSoft > 0);
	    }
	}
    }

void pkdPreVariableSoft(PKD pkd)
    {
    PARTICLE *p;
    int i,n;

    p = pkd->pStore;
    n = pkdLocal(pkd);
	
    for (i=0;i<n;++i) {
	if (pkdIsActive(pkd,&(p[i]))) p[i].fSoft = 0.5*p[i].fBall;
	}
    }

void pkdPostVariableSoft(PKD pkd,double dSoftMax,int bSoftMaxMul)
    {
    PARTICLE *p;
    int i,n;
    double dTmp;

    p = pkd->pStore;
    n = pkdLocal(pkd);
	
    if (bSoftMaxMul) {
	for (i=0;i<n;++i) {
	    if (pkdIsActive(pkd,&(p[i]))) {
		dTmp = 0.5*p[i].fBall;
		p[i].fBall = 2.0*p[i].fSoft;
		p[i].fSoft = (dTmp <= p[i].fSoft0*dSoftMax ? dTmp : p[i].fSoft0*dSoftMax);
		}
	    }
	}
    else {
	for (i=0;i<n;++i) {
	    if (pkdIsActive(pkd,&(p[i]))) {
		dTmp = 0.5*p[i].fBall;
		p[i].fBall = 2.0*p[i].fSoft;
		p[i].fSoft = (dTmp <= dSoftMax ? dTmp : dSoftMax);
		}
	    }
	}	
    }
#endif


void pkdBucketWeight(PKD pkd,int iBucket,FLOAT fWeight)
    {
    KDN *pbuc;
    int pj;
	
    pbuc = &pkd->kdNodes[iBucket];
    for (pj=pbuc->pLower;pj<=pbuc->pUpper;++pj) {
	if (pkdIsActive(pkd,&(pkd->pStore[pj])))
	    pkd->pStore[pj].fWeight = fWeight;
	}
    }


void
pkdGravAll(PKD pkd,double dTime,int nReps,int bPeriodic,int iOrder,int bEwald,
	   int iEwOrder, double fEwCut,double fEwhCut,int *nActive, 
	   double *pdPartSum, double *pdCellSum,CASTAT *pcs, double *pdFlop)
    {
    int bVeryActive = 0;

    pkdClearTimer(pkd,1);
    pkdClearTimer(pkd,2);
    pkdClearTimer(pkd,3);

#ifdef BSC
    MPItrace_event(10000,0);
#endif

    /*
    ** Set up Ewald tables and stuff.
    */
    if (bPeriodic && bEwald) {
	pkdEwaldInit(pkd,fEwhCut,iEwOrder);	/* ignored in Flop count! */
	}
    /*
    ** Start particle caching space (cell cache already active).
    */
    mdlROcache(pkd->mdl,CID_PARTICLE,pkd->pStore,sizeof(PARTICLE),
	       pkdLocal(pkd));
    /*
    ** Calculate newtonian gravity, including replicas if any.
    */
    *pdFlop = 0.0;
    *pdPartSum = 0.0;
    *pdCellSum = 0.0;
    pkdStartTimer(pkd,1);
    *nActive = pkdGravWalk(pkd,dTime,nReps,bPeriodic && bEwald,bVeryActive,fEwCut,pdFlop,pdPartSum,pdCellSum);
    pkdStopTimer(pkd,1);

#ifdef BSC
    /*MPItrace_event(10001, (int)(pkdGetWallClockTimer(pkd,1)*1000000) );*/
    /*MPItrace_event(10001, (int)(pkd->nActive/(pkd->nLocal/100)) );*/
    MPItrace_event(10001, 3 );
#endif

    /*
    ** Get caching statistics.
    */
    pcs->dcNumAccess = mdlNumAccess(pkd->mdl,CID_CELL);
    pcs->dcMissRatio = mdlMissRatio(pkd->mdl,CID_CELL);
    pcs->dcCollRatio = mdlCollRatio(pkd->mdl,CID_CELL);
    pcs->dcMinRatio = mdlMinRatio(pkd->mdl,CID_CELL);
    pcs->dpNumAccess = mdlNumAccess(pkd->mdl,CID_PARTICLE);
    pcs->dpMissRatio = mdlMissRatio(pkd->mdl,CID_PARTICLE);
    pcs->dpCollRatio = mdlCollRatio(pkd->mdl,CID_PARTICLE);
    pcs->dpMinRatio = mdlMinRatio(pkd->mdl,CID_PARTICLE);
    /*
    ** Stop particle caching space.
    */
    mdlFinishCache(pkd->mdl,CID_PARTICLE);
    }


void pkdCalcEandL(PKD pkd,double *T,double *U,double *Eth,double L[])
    {
    /* L is calculated with respect to the origin (0,0,0) */

    PARTICLE *p;
    FLOAT rx,ry,rz,vx,vy,vz;
    int i,n;

    p = pkd->pStore;
    n = pkdLocal(pkd);
    *T = 0.0;
    *U = 0.0;
    *Eth = 0.0;
    L[0] = L[1] = L[2] = 0;
    for (i=0;i<n;++i) {
	rx = p[i].r[0]; ry = p[i].r[1]; rz = p[i].r[2];
	vx = p[i].v[0]; vy = p[i].v[1]; vz = p[i].v[2];
	*T += 0.5*p[i].fMass*(vx*vx + vy*vy + vz*vz);
	*U += 0.5*p[i].fMass*p[i].fPot;
	L[0] += p[i].fMass*(ry*vz - rz*vy);
	L[1] += p[i].fMass*(rz*vx - rx*vz);
	L[2] += p[i].fMass*(rx*vy - ry*vx);
	}
    }


void
pkdDrift(PKD pkd,double dTime,double dDelta,FLOAT fCenter[3],int bPeriodic,int bFandG,
	 FLOAT fCentMass) {

    PARTICLE *p;
    int i,j,n;
    int ipt;
    double px,py,pz,pm;
    double dDeltaM;
    
    mdlDiag(pkd->mdl, "Into pkdDrift\n");
    
    if (dTime >= pkd->param.dGrowStartT && dTime < pkd->param.dGrowEndT) {
	dDeltaM = pkd->param.dGrowDeltaM*dDelta/
	    (pkd->param.dGrowEndT - pkd->param.dGrowStartT);
	}
    else {
	dDeltaM = 0;
	}

#ifdef BSC
    MPItrace_event(10000,0);
#endif
    
    p = pkd->pStore;
    n = pkdLocal(pkd);
    for (i=0;i<n;++i) {
	/*
	** Do the Growmass stuff if needed...
	*/
	if (p[i].iOrder < pkd->param.nGrowMass) {
	    p[i].fMass += dDeltaM;
	    }
	/*
	** Update particle positions
	*/
	for (j=0;j<3;++j) {
	    p[i].r[j] += dDelta*p[i].v[j];
	    if (bPeriodic) {
		if (p[i].r[j] >= fCenter[j] + 0.5*pkd->fPeriod[j]) {
		    p[i].r[j] -= pkd->fPeriod[j];
		    }
		if (p[i].r[j] < fCenter[j] - 0.5*pkd->fPeriod[j]) {
		    p[i].r[j] += pkd->fPeriod[j];
		    }
		}
	    }


	/*
	** Detailed output for particles that are tracked at dTime + dDelta/2
	** Correct positons & mass for dDelta/2 step which is done before the actual update
	*/
	if (TYPETest(&p[i],TYPE_TRACKER) && dDelta != 0) {
      	    px = p[i].r[0] - dDelta/2*p[i].v[0];
	    py = p[i].r[1] - dDelta/2*p[i].v[1];
	    pz = p[i].r[2] - dDelta/2*p[i].v[2];
	    pm = p[i].fMass - dDeltaM/2;
	    ipt = p[i].iOrder + 1; /* to be consistent with Tipsy numbering */
	    printf("PDID %d %g %g %d %g %g %g %g %g %g %g %g %g %g %g %g %g\n",ipt,dTime+dDelta/2,dDelta,
		   p[i].iRung,p[i].dt,px,py,pz,p[i].v[0],p[i].v[1],p[i].v[2],p[i].a[0],p[i].a[1],p[i].a[2],p[i].fPot,pm,p[i].fSoft);
	    }
	}

#ifdef BSC
    MPItrace_event(10001,4);
#endif
    mdlDiag(pkd->mdl, "Out of pkdDrift\n");
    }

void
pkdDriftInactive(PKD pkd,double dTime, double dDelta,FLOAT fCenter[3],int bPeriodic,
		 int bFandG, FLOAT fCentMass) {

    PARTICLE *p;
    int i,j,n;
    int ipt;
    double dDeltaM;
    
    mdlDiag(pkd->mdl, "Into pkdDriftInactive\n");
    
    if (dTime >= pkd->param.dGrowStartT && dTime < pkd->param.dGrowEndT) {
	dDeltaM = pkd->param.dGrowDeltaM*dDelta/
	    (pkd->param.dGrowEndT - pkd->param.dGrowStartT);
	}
    else {
	dDeltaM = 0;
	}
    
    p = pkd->pStore;
    n = pkdLocal(pkd);
    for (i=0;i<n;++i) {
	if (pkdIsActive(pkd,&p[i])) continue;
	/*
	** Do the Growmass stuff if needed...
	*/
	if (p[i].iOrder < pkd->param.nGrowMass) {
	    p[i].fMass += dDeltaM;
	    }
	/*
	** Update particle positions
	*/
	for (j=0;j<3;++j) {
	    p[i].r[j] += dDelta*p[i].v[j];
	    if (bPeriodic) {
		if (p[i].r[j] >= fCenter[j] + 0.5*pkd->fPeriod[j]) {
		    p[i].r[j] -= pkd->fPeriod[j];
		    }
		if (p[i].r[j] < fCenter[j] - 0.5*pkd->fPeriod[j]) {
		    p[i].r[j] += pkd->fPeriod[j];
		    }
		}
	    }
	/*
	** Detailed output for particles that are tracked at dTime + dDelta
	** Inactive Drift => dDelta is already a half step! No correction needed!
	*/
	if (TYPETest(&p[i],TYPE_TRACKER) && dDelta != 0) {
	    ipt = p[i].iOrder + 1; /* to be consistent with Tipsy numbering */
	    printf("PDID %d %g %g %d %g %g %g %g %g %g %g %g %g %g %g %g %g\n",ipt,dTime+dDelta/2,dDelta,
		   p[i].iRung,p[i].dt,p[i].r[0],p[i].r[1],p[i].r[2],p[i].v[0],p[i].v[1],p[i].v[2],p[i].a[0],p[i].a[1],p[i].a[2],p[i].fPot,p[i].fMass,p[i].fSoft);
	    }
	}
    mdlDiag(pkd->mdl, "Out of pkdDriftInactive\n");
    }

void
pkdDriftActive(PKD pkd,double dTime,double dDelta) {

    PARTICLE *p;
    int i,j,n;
    int ipt;
    double px,py,pz,pm;
    double dDeltaM;
    
    mdlDiag(pkd->mdl, "Into pkdDriftActive\n");

    if (dTime >= pkd->param.dGrowStartT && dTime < pkd->param.dGrowEndT) {
	dDeltaM = pkd->param.dGrowDeltaM*dDelta/
	    (pkd->param.dGrowEndT - pkd->param.dGrowStartT);
	}
    else {
	dDeltaM = 0;
	}

    p = pkd->pStore;
    n = pkdLocal(pkd);
    for (i=0;i<n;++i) {
	if (!pkdIsActive(pkd,&p[i])) continue;
	/*
	** Do the Growmass stuff if needed...
	*/
	if (p[i].iOrder < pkd->param.nGrowMass) {
	    p[i].fMass += dDeltaM;
	    }
	/*
	** Update particle positions
	*/
	for (j=0;j<3;++j) {
	    p[i].r[j] += dDelta*p[i].v[j];
	    }
	/*
	** Detailed output for particles that are tracked at dTime + dDelta/2
	** Correct positons & mass for dDelta/2 step
	*/
	if (TYPETest(&p[i],TYPE_TRACKER) && dDelta != 0) {
	    px = p[i].r[0] - dDelta/2*p[i].v[0];
	    py = p[i].r[1] - dDelta/2*p[i].v[1];
	    pz = p[i].r[2] - dDelta/2*p[i].v[2];
	    pm = p[i].fMass - dDeltaM/2;
	    ipt = p[i].iOrder + 1; /* to be consistent with Tipsy numbering */
	    printf("PDID %d %g %g %d %g %g %g %g %g %g %g %g %g %g %g %g %g\n",ipt,dTime+dDelta/2,dDelta,
		   p[i].iRung,p[i].dt,px,py,pz,p[i].v[0],p[i].v[1],p[i].v[2],p[i].a[0],p[i].a[1],p[i].a[2],p[i].fPot,pm,p[i].fSoft);
	    }
	}
    mdlDiag(pkd->mdl, "Out of pkdDriftActive\n");
    }

void pkdGravityVeryActive(PKD pkd,double dTime,int bEwald,int nReps,double fEwCut,double dStep) {
    int nActive;
    int bVeryActive = 1;
    double dFlop,dPartSum,dCellSum;

    /*
    ** Calculate newtonian gravity for the very active particles ONLY, including replicas if any.
    */
    dFlop = 0.0;
    dPartSum = 0.0;
    dCellSum = 0.0;
    nActive = pkdGravWalk(pkd,dTime,nReps,bEwald,bVeryActive,fEwCut,&dFlop,&dPartSum,&dCellSum);
    }


void
pkdStepVeryActiveKDK(PKD pkd, double dStep, double dTime, double dDelta,
		     int iRung, int iKickRung, int iRungVeryActive,int iAdjust, double diCrit2,
		     int *pnMaxRung, double aSunInact[], double adSunInact[], double dSunMass)
    {
    int nRungCount[256];
    double dDriftFac;
#ifdef PLANETS
    double aSun[3], adSun[3];
    int j;
#endif

    double time1,time2; /* added MZ 1.6.2006 */
    
    if(iAdjust && (iRung < pkd->param.iMaxRung-1)) {
	pkdActiveRung(pkd, iRung, 1);
	pkdInitDt(pkd, pkd->param.dDelta);
	if (pkd->param.bGravStep) {
	    double a = csmTime2Exp(pkd->param.csm,dTime);
	    double dRhoFac = 1.0/(a*a*a);
	    pkdGravStep(pkd,pkd->param.dEta,dRhoFac);
	    }
	if (pkd->param.bAccelStep) {
	    double a = csmTime2Exp(pkd->param.csm,dTime);
	    double dVelFac = 1.0/(a*a);
	    double dAccFac = 1.0/(a*a*a);
	    double dhMinOverSoft = 0;
	    pkdAccelStep(pkd,pkd->param.dEta, dVelFac,dAccFac,pkd->param.bDoGravity,
			 pkd->param.bEpsAccStep,pkd->param.bSqrtPhiStep,dhMinOverSoft);
	    }
	*pnMaxRung = pkdDtToRung(pkd,iRung,dDelta,pkd->param.iMaxRung-1, 1, nRungCount);
    
	if (pkd->param.bVDetails) {
	    printf("%*cAdjust at iRung: %d, nMaxRung:%d nRungCount[%d]=%d\n",
		   2*iRung+2,' ',iRung,*pnMaxRung,*pnMaxRung,nRungCount[*pnMaxRung]);
	}
	
	}
    if(iRung > iRungVeryActive) {	/* skip this if we are
					   entering for the first
					   time: Kick is taken care of
					   in master(). 
					*/
	pkdActiveRung(pkd,iRung,0);
	if (pkd->param.bVDetails) {
	    printf("%*cVeryActive pkdKickOpen  at iRung: %d, 0.5*dDelta: %g\n",
		   2*iRung+2,' ',iRung,0.5*dDelta);
	    }
	pkdKickKDKOpen(pkd, dTime, 0.5*dDelta);
	}
    if (*pnMaxRung > iRung) {
	/*
	** Recurse.
	*/
	pkdStepVeryActiveKDK(pkd,dStep,dTime,0.5*dDelta,iRung+1,iRung+1,iRungVeryActive,0,
			     diCrit2,pnMaxRung,aSunInact,adSunInact,dSunMass);
	dStep += 1.0/(2 << iRung);
	dTime += 0.5*dDelta;
	pkdActiveRung(pkd,iRung,0);
	pkdStepVeryActiveKDK(pkd,dStep,dTime,0.5*dDelta,iRung+1,iKickRung,iRungVeryActive,1,
			     diCrit2,pnMaxRung,aSunInact,adSunInact,dSunMass);
	}
    else {
	if (pkd->param.bVDetails) {
	    printf("%*cVeryActive Drift at iRung: %d, drifting %d and higher with dDelta: %g\n",
		   2*iRung+2,' ',iRung,iRungVeryActive+1,dDelta);
	    }
	/*
	** This should drift *all* very actives!
	*/
	pkdActiveRung(pkd,iRungVeryActive+1,1);
	/*
	** We need to account for cosmological drift factor here!
	** Normally this is done at the MASTER level in msrDrift.
	** Note that for kicks we have written new "master-like" functions
	** KickOpen and KickClose which do this same job at PKD level.
	*/
	if (pkd->param.bCannonical) {
	    dDriftFac = csmComoveDriftFac(pkd->param.csm,dTime,dDelta);
	    }
	else {
	    dDriftFac = dDelta;
	    }

	pkdDriftActive(pkd,dTime,dDriftFac);
	dTime += dDelta;
	dStep += 1.0/(1 << iRung);

	if(iKickRung > iRungVeryActive) {	/* skip this if we are
						   entering for the first
						   time: Kick is taken care of
						   in master(). 
						*/

	    if(pkd->param.bVDetails) {
		printf("%*cGravityVA: iRung %d Gravity for rungs %d to %d ... ",
		       2*iRung+2,' ',iRung,iKickRung,*pnMaxRung);
		}

	    time1 = Zeit(); /* added MZ 1.6.2006 */

	    pkdActiveRung(pkd,iKickRung,1);
	    pkdVATreeBuild(pkd,pkd->param.nBucket,diCrit2,0,dTime);
	    pkdGravityVeryActive(pkd,dTime,pkd->param.bEwald && pkd->param.bPeriodic,pkd->param.nReplicas,pkd->param.dEwCut,dStep);

#ifdef PLANETS
	    /* Sun's gravity */
	    if(pkd->param.bHeliocentric){
	      /* Sun's indirect gravity due to very active particles*/
	      pkdActiveRung(pkd,iRungVeryActive+1,1);
	      pkdSunIndirect(pkd,aSun,adSun,1); 
	      for (j=0;j<3;++j)
	       	{                 
                 aSun[j] += aSunInact[j];
		 adSun[j] += adSunInact[j];
	       	}
	      pkdActiveRung(pkd,iKickRung,1);
	      pkdGravSun(pkd,aSun,adSun,dSunMass); 
	    }
#endif

	    time2 = Zeit();
	    if(pkd->param.bVDetails)
		printf("Time: %g\n",time2-time1);
	    }
	/*
	 * move time back to 1/2 step so that KickClose can integrate
	 * from 1/2 through the timestep to the end.
	 */
	dTime -= 0.5*dDelta;
	}
    if(iKickRung > iRungVeryActive) {	/* skip this if we are
						   entering for the first
						   time: Kick is taken care of
						   in master(). 
						*/
	pkdActiveRung(pkd,iRung,0);
	if (pkd->param.bVDetails) {
	    printf("%*cVeryActive pkdKickClose at iRung: %d, 0.5*dDelta: %g\n",
		   2*iRung+2,' ',iRung,0.5*dDelta);
	    }
	pkdKickKDKClose(pkd,dTime,0.5*dDelta);
	}
    }

#ifdef HERMITE
void
pkdStepVeryActiveHermite(PKD pkd, double dStep, double dTime, double dDelta,
		     int iRung, int iKickRung, int iRungVeryActive,int iAdjust, double diCrit2,
		     int *pnMaxRung, double aSunInact[], double adSunInact[], double dSunMass)
    {
    int nRungCount[256];
    double dDriftFac;
    int i;

    double aSun[3], adSun[3];
    int j;

    double time1,time2; /* added MZ 1.6.2006 */
    
    if(iAdjust && (iRung < pkd->param.iMaxRung-1)) {
	pkdActiveRung(pkd, iRung, 1);
	pkdInitDt(pkd, pkd->param.dDelta);
	if (pkd->param.bAarsethStep) {
	  pkdAarsethStep(pkd,pkd->param.dEta);
	}
	if (pkd->param.bGravStep) {
	    double a = csmTime2Exp(pkd->param.csm,dTime);
	    double dRhoFac = 1.0/(a*a*a);
	    pkdGravStep(pkd,pkd->param.dEta,dRhoFac);
	    }
	if (pkd->param.bAccelStep) {
	    double a = csmTime2Exp(pkd->param.csm,dTime);
	    double dVelFac = 1.0/(a*a);
	    double dAccFac = 1.0/(a*a*a);
	    double dhMinOverSoft = 0;
	    pkdAccelStep(pkd,pkd->param.dEta, dVelFac,dAccFac,pkd->param.bDoGravity,
			 pkd->param.bEpsAccStep,pkd->param.bSqrtPhiStep,dhMinOverSoft);
	    }
	*pnMaxRung = pkdDtToRung(pkd,iRung,dDelta,pkd->param.iMaxRung-1, 1, nRungCount);
    
	if (pkd->param.bVDetails) {
	    printf("%*cAdjust at iRung: %d, nMaxRung:%d nRungCount[%d]=%d\n",
		   2*iRung+2,' ',iRung,*pnMaxRung,*pnMaxRung,nRungCount[*pnMaxRung]);
	}
	
	}
   
    if (*pnMaxRung > iRung) {
	/*
	** Recurse.
	*/
	pkdStepVeryActiveHermite(pkd,dStep,dTime,0.5*dDelta,iRung+1,iRung+1,iRungVeryActive,0,
			     diCrit2,pnMaxRung,aSunInact,adSunInact,dSunMass);
	dStep += 1.0/(2 << iRung);
	dTime += 0.5*dDelta;
	pkdActiveRung(pkd,iRung,0);
	pkdStepVeryActiveHermite(pkd,dStep,dTime,0.5*dDelta,iRung+1,iKickRung,iRungVeryActive,1,
			     diCrit2,pnMaxRung,aSunInact,adSunInact,dSunMass);
	}
    else {
	if (pkd->param.bVDetails) {
	    printf("%*cVeryActive Predictor at iRung: %d, drifting %d and higher with dDelta: %g\n",
		   2*iRung+2,' ',iRung,iRungVeryActive+1,dDelta);
	    }
	/*
	** This should predict *all* very actives!
	*/
	pkdActiveRung(pkd,iRungVeryActive+1,1);
	/*
	** We need to account for cosmological drift factor here!
	** Normally this is done at the MASTER level in msrDrift.
	** Note that for kicks we have written new "master-like" functions
	** KickOpen and KickClose which do this same job at PKD level.
	*/
	/*if (pkd->param.bCannonical) {
	    dDriftFac = csmComoveDriftFac(pkd->param.csm,dTime,dDelta);
	    }
	else {
	    dDriftFac = dDelta;
	    }*/

	dTime += dDelta;
	dStep += 1.0/(1 << iRung);

	pkdPredictor(pkd,dTime);

	if(iKickRung > iRungVeryActive) {	/* skip this if we are
						   entering for the first
						   time: Kick is taken care of
						   in master(). 
						*/

	    if(pkd->param.bVDetails) {
		printf("%*cGravityVA: iRung %d Gravity for rungs %d to %d ...\n ",
		       2*iRung+2,' ',iRung,iKickRung,*pnMaxRung);
		}

	    time1 = Zeit(); /* added MZ 1.6.2006 */

	    pkdActiveRung(pkd,iKickRung,1);
	    pkdVATreeBuild(pkd,pkd->param.nBucket,diCrit2,0,dTime);
	    pkdGravityVeryActive(pkd,dTime,pkd->param.bEwald && pkd->param.bPeriodic,pkd->param.nReplicas,pkd->param.dEwCut,dStep);


#ifdef PLANETS
	    /* Sun's gravity */
           if(pkd->param.bHeliocentric){
	       /* Sun's indirect gravity due to very active particles*/
	   pkdActiveRung(pkd,iRungVeryActive+1,1);
	   pkdSunIndirect(pkd,aSun,adSun,1); 
           for (j=0;j<3;++j)
	       	{                 
                 aSun[j] += aSunInact[j];
		 adSun[j] += adSunInact[j];
	       	}
	    pkdActiveRung(pkd,iKickRung,1);
	    pkdGravSun(pkd,aSun,adSun,dSunMass); 
           }
#endif

	   if (pkd->param.bVDetails) {
	     printf("%*cVeryActive pkdCorrector at iRung: %d, 0.5*dDelta: %g\n",
		    2*iRung+2,' ',iRung,0.5*dDelta);
	    }
	   pkdCorrector(pkd,dTime);

#ifdef PLANETS
	   if(pkd->param.bHeliocentric){
	     int nite = 0;
	     int nitemax = 3;                                                    
	     do{    
	       nite += 1;
	       pkdSunCorrector(pkd,dTime,dSunMass);  		
	     }while(nite < nitemax);
           }
	   if(pkd->param.bCollision && pkd->iCollisionflag) pkdDoCollisionVeryActive(pkd,dTime); 
#endif
	   pkdCopy0(pkd,dTime); 

	   time2 = Zeit();
	   if(pkd->param.bVDetails)
	     printf("Time: %g\n",time2-time1);
	    
	}
      }   
    }

void
pkdCopy0(PKD pkd,double dTime) 
{  
    PARTICLE *p;
    int i,j,n;
    
    mdlDiag(pkd->mdl, "Into pkdCopy0\n");
            
    p = pkd->pStore;
    n = pkdLocal(pkd);
    for (i=0;i<n;++i) {	
	if (pkdIsActive(pkd,&p[i])) { 			      
	    for (j=0;j<3;++j) {    
		p[i].r0[j] = p[i].r[j];
		p[i].v0[j] = p[i].v[j];
		p[i].a0[j] = p[i].a[j];	         	       
		p[i].ad0[j] = p[i].ad[j];	 
	    }
	    p[i].dTime0 = dTime;
#ifdef PLANETS
	    if(pkd->param.bCollision){
		p[i].iColflag = 0;  /* just in case */
	    }
#endif
	}
    }
    mdlDiag(pkd->mdl, "Out of pkdCopy0\n");
}

void
pkdPredictor(PKD pkd,double dTime) 
{
    PARTICLE *p;
    int i,j,n;
    double dt; /* time since last evaluation of force */
    
    mdlDiag(pkd->mdl, "Into pkdPredictor\n");

    if (pkd->param.bVDetails) printf("Into pkdPredictor\n");      

    p = pkd->pStore;
    n = pkdLocal(pkd);
    for (i=0;i<n;++i) {		
	if (pkdIsActive(pkd,&p[i])) { 	
	    dt =  dTime - p[i].dTime0; 
      /*if (pkd->param.bVDetails)
	{
         printf("Particle=%d iRung=%d dt=%g \n",p[i].iOrder,p[i].iRung,dt);
	 }*/
	    for (j=0;j<3;++j) {                
		p[i].r[j] = p[i].r0[j] + dt*(p[i].v0[j] + 0.5*dt*(p[i].a0[j]+dt*p[i].ad0[j]/3.0));
		p[i].v[j] = p[i].v0[j] + dt*(p[i].a0[j]+0.5*dt*p[i].ad0[j]);   
                p[i].rp[j] = p[i].r[j];
                p[i].vp[j] = p[i].v[j];  	       
		}
           }	    
	}
    mdlDiag(pkd->mdl, "Out of pkdPredictor\n");
    }

void
pkdCorrector(PKD pkd,double dTime) 
{
    PARTICLE *p;
    int i,j,n;
    double dt, add, addd, am; 
    double a0d2, a1d2, a2d2, a3d2, dti, dt2i;
    /*double alpha = 7.0/6.0;*/    

    mdlDiag(pkd->mdl, "Into pkdCorrector\n");

    if (pkd->param.bVDetails) printf("Into pkdCorrector\n");

    p = pkd->pStore;
    n = pkdLocal(pkd);
    for (i=0;i<n;++i) {	
	if (pkdIsActive(pkd,&p[i])) {	
	      dt =  dTime - p[i].dTime0;  
	      if(pkd->param.bAarsethStep){	    
		a0d2  = 0.0;
		a1d2  = 0.0;
		a2d2  = 0.0;
		a3d2  = 0.0;
		dti = 1.0/dt;
		dt2i = dti*dti;
	      }
	      
	    for (j=0;j<3;++j) { 	  
              am = p[i].a0[j]-p[i].a[j];

	      if(pkd->param.bAarsethStep){
	      add = 2.0*(3.0*am*dt2i+(2.0*p[i].ad0[j]+p[i].ad[j])*dti);
	      addd = 6.0*(2.0*am*dti+(p[i].ad0[j]+p[i].ad[j]))*dt2i;
              add += dt*addd; 

              a0d2 += p[i].a[j]*p[i].a[j];
              a1d2 += p[i].ad[j]*p[i].ad[j];
	      a2d2 += add*add;  
              a3d2 += addd*addd;	      
              }

	      add = 16.0*am+(13.0*p[i].ad0[j]+3.0*p[i].ad[j])*dt;
	      addd = 6.0*am+(5.0*p[i].ad0[j]+p[i].ad[j])*dt;
	      p[i].r[j] = p[i].rp[j] - add/120.0*dt*dt;
	      p[i].v[j] = p[i].vp[j] - addd/12.0*dt; 
	
	    }
	    if(pkd->param.bAarsethStep){
	      p[i].dtGrav = (sqrt(a0d2*a2d2)+a1d2)/(sqrt(a1d2*a3d2)+a2d2); 
	    }
	    /*if(p[i].dtGrav <=0.0) p[i].dtGrav = 0.1*(a0d2/a1d2);*/
	    }	   
    }
    mdlDiag(pkd->mdl, "Out of pkdCorrector\n");
}

void
pkdSunCorrector(PKD pkd,double dTime,double dSunMass) 
{
  /* iteratively correct only sun's direct gravity */
    PARTICLE *p;
    int i,j,n;
    double dt;
    double add, addd, r2,r3i,r5i,rv, am;
    /*double alpha = 7.0/6.0;*/   

    mdlDiag(pkd->mdl, "Into pkdSunCorrector\n");
            
    p = pkd->pStore;
    n = pkdLocal(pkd);
    for (i=0;i<n;++i) {	
	if (pkdIsActive(pkd,&p[i])) { 		
	      dt =  dTime - p[i].dTime0;
	      r2  = 0.0;
              rv  = 0.0;  
             for (j=0;j<3;++j) {       
	       r2 += p[i].r[j]*p[i].r[j];
               rv += p[i].r[j]*p[i].v[j];     
	       }              
	     r2 = 1.0/r2;
               r3i = r2*sqrt(r2)*dSunMass;
	       r5i = 3.0*r3i*r2*rv;

	    for (j=0;j<3;++j) {   	   
	      p[i].a[j] = -p[i].r[j]*r3i + p[i].app[j];
	      p[i].ad[j] = -(p[i].v[j]*r3i - p[i].r[j]*r5i) + p[i].adpp[j];
	      am = p[i].a0[j]-p[i].a[j];
	      add = 16.0*am+(13.0*p[i].ad0[j]+3.0*p[i].ad[j])*dt;
	      addd = 6.0*am+(5.0*p[i].ad0[j]+p[i].ad[j])*dt;
		p[i].r[j] = p[i].rp[j] - add/120.0*dt*dt;
		p[i].v[j] = p[i].vp[j] - addd/12.0*dt;     
	                        
		}	    
	    }
	}
    mdlDiag(pkd->mdl, "Out of pkdSunCorrector\n");
    }

void
pkdPredictorInactive(PKD pkd,double dTime) 
{
    PARTICLE *p;
    int i,j,n;
    double dt; /* time since last evaluation of force */
    
    mdlDiag(pkd->mdl, "Into pkdPredictorInactive\n");

    if (pkd->param.bVDetails) printf("Into pkdPredictorInactive\n");      

    p = pkd->pStore;
    n = pkdLocal(pkd);
    for (i=0;i<n;++i) {		
	if (pkdIsActive(pkd,&p[i])) continue;	
	dt =  dTime - p[i].dTime0; 
	/*if (pkd->param.bVDetails)
	  {
	  printf("Particle=%d iRung=%d dt=%g \n",p[i].iOrder,p[i].iRung,dt);
	  }*/
	for (j=0;j<3;++j) {                
	    p[i].r[j] = p[i].r0[j] + dt*(p[i].v0[j] + 0.5*dt*(p[i].a0[j]+dt*p[i].ad0[j]/3.0));
	    p[i].v[j] = p[i].v0[j] + dt*(p[i].a0[j]+0.5*dt*p[i].ad0[j]);   
	    p[i].rp[j] = p[i].r[j];
	    p[i].vp[j] = p[i].v[j];  	       
	}	    
    }
    mdlDiag(pkd->mdl, "Out of pkdPredictorInactive\n");
}

void pkdAarsethStep(PKD pkd,double dEta) {
    double dT;
    int i;
    
    for (i=0;i<pkdLocal(pkd);i++) {
	if (pkdIsActive(pkd,&(pkd->pStore[i]))) {
	    mdlassert(pkd->mdl, pkd->pStore[i].dtGrav > 0);
	    dT = dEta*sqrt(pkd->pStore[i].dtGrav);
	    if (dT < pkd->pStore[i].dt)
	    pkd->pStore[i].dt = dT;
	   
	    }
	}
    }

void pkdFirstDt(PKD pkd) {
    int i,j;
    PARTICLE *p ;
    double a0d2,a1d2;
  
    p = pkd->pStore;
    for(i=0;i<pkdLocal(pkd);++i) {
          a0d2 = 0.0;
	  a1d2 = 0.0;
	  if(pkdIsActive(pkd,&p[i]))
	      for (j=0;j<3;++j) { 	                  
		  a0d2 += p[i].a0[j]*p[i].a0[j];
		  a1d2 += p[i].ad0[j]*p[i].ad0[j];
	      }
	  p[i].dtGrav = (a0d2/a1d2);
#ifdef PLANETS
	  if(pkd->param.bCollision){
	      p[i].iColflag = 0; /* initial reset of collision flag */	
	  }
#endif
    }
}
#endif /* Hermite */

/*
 * Stripped down versions of routines from master.c
 */
void pkdKickKDKOpen(PKD pkd,double dTime,double dDelta)
    {
    double H,a;
    double dvFacOne, dvFacTwo;
	
    if (pkd->param.bCannonical) {
	dvFacOne = 1.0;		/* no hubble drag, man! */
	dvFacTwo = csmComoveKickFac(pkd->param.csm,dTime,dDelta);
	}
    else {
	/*
	** Careful! For non-cannonical we want H and a at the 
	** HALF-STEP! This is a bit messy but has to be special
	** cased in some way.
	*/
	dTime += dDelta/2.0;
	a = csmTime2Exp(pkd->param.csm,dTime);
	H = csmTime2Hub(pkd->param.csm,dTime);
	dvFacOne = (1.0 - H*dDelta)/(1.0 + H*dDelta);
	dvFacTwo = dDelta/pow(a,3.0)/(1.0 + H*dDelta);
	}
    pkdKick(pkd, dvFacOne, dvFacTwo, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0);
    }

void pkdKickKDKClose(PKD pkd,double dTime,double dDelta)
    {
    double H,a;
    double dvFacOne, dvFacTwo;
	
    if (pkd->param.bCannonical) {
	dvFacOne = 1.0; /* no hubble drag, man! */
	dvFacTwo = csmComoveKickFac(pkd->param.csm,dTime,dDelta);
	}
    else {
	/*
	** Careful! For non-cannonical we want H and a at the 
	** HALF-STEP! This is a bit messy but has to be special
	** cased in some way.
	*/
	dTime += dDelta/2.0;
	a = csmTime2Exp(pkd->param.csm,dTime);
	H = csmTime2Hub(pkd->param.csm,dTime);
	dvFacOne = (1.0 - H*dDelta)/(1.0 + H*dDelta);
	dvFacTwo = dDelta/pow(a,3.0)/(1.0 + H*dDelta);
	}
    pkdKick(pkd, dvFacOne, dvFacTwo, 0.0, 0.0, 0.0, 0.0, 0, 0.0, 0.0);
    }

void pkdKick(PKD pkd, double dvFacOne, double dvFacTwo, double dvPredFacOne,
	     double dvPredFacTwo, double duDelta, double duPredDelta, int iGasModel,
	     double z, double duDotLimit)
    {
    PARTICLE *p;
    int i,j,n;

    pkdClearTimer(pkd,1);
    pkdStartTimer(pkd,1);

    p = pkd->pStore;
    n = pkdLocal(pkd);
    for (i=0;i<n;++i,++p) {
	if (pkdIsActive(pkd,p)) {
	    for (j=0;j<3;++j) {
		p->v[j] = p->v[j]*dvFacOne + p->a[j]*dvFacTwo;
		}
	    }
	}

    pkdStopTimer(pkd,1);
    mdlDiag(pkd->mdl, "Done pkdkick\n");
    }


void pkdInitStep(PKD pkd, struct parameters *p, CSM csm) {
    pkd->param = *p;
    /*
    ** Need to be careful to correctly copy the cosmo
    ** parameters. This is very ugly!
    */
    csmInitialize(&pkd->param.csm);
    *pkd->param.csm = *csm;
    }


void pkdSetRung(PKD pkd, int iRung) {
    int i;
    
    for(i=0;i<pkdLocal(pkd);++i) {
	pkd->pStore[i].iRung = iRung;
	}
    }

void pkdActiveRung(PKD pkd, int iRung, int bGreater) {
    pkd->uMinRungActive = iRung;
    pkd->uMaxRungActive = bGreater ? 255 : iRung;
    }

int pkdCurrRung(PKD pkd, int iRung) {
    int i;
    int iCurrent;
    
    iCurrent = 0;
    for(i=0;i<pkdLocal(pkd);++i) {
	if(pkd->pStore[i].iRung == iRung) {
	    iCurrent = 1;
	    break;
	    }
	}
    return iCurrent;
    }

void pkdGravStep(PKD pkd,double dEta,double dRhoFac) {
    double dT;
    int i;
    
    for (i=0;i<pkdLocal(pkd);i++) {
	if (pkdIsActive(pkd,&(pkd->pStore[i]))) {
	    mdlassert(pkd->mdl, pkd->pStore[i].dtGrav > 0);
	    dT = dEta/sqrt(pkd->pStore[i].dtGrav*dRhoFac);
	    if (dT < pkd->pStore[i].dt) {
		pkd->pStore[i].dt = dT;	
	
		  mdlassert(pkd->mdl,dT>0);
	        }
	    }
	}
    }

void pkdAccelStep(PKD pkd,double dEta,double dVelFac,double dAccFac,int bDoGravity,int bEpsAcc,int bSqrtPhi,double dhMinOverSoft) {
    int i;
    double vel;
    double acc;
    int j;
    double dT;

    for (i=0;i<pkdLocal(pkd);++i) {
	if (pkdIsActive(pkd,&(pkd->pStore[i]))) {
	    vel = 0;
	    acc = 0;
	    for (j=0;j<3;j++) {
		vel += pkd->pStore[i].v[j]*pkd->pStore[i].v[j];
		acc += pkd->pStore[i].a[j]*pkd->pStore[i].a[j];
		}
	    mdlassert(pkd->mdl,vel >= 0);
	    vel = sqrt(vel)*dVelFac;
	    mdlassert(pkd->mdl,acc >= 0);
	    acc = sqrt(acc)*dAccFac;
	    mdlassert(pkd->mdl,acc > 0);
	    dT = FLOAT_MAXVAL;
	    if (bEpsAcc && acc>0) {
		dT = dEta*sqrt(pkd->pStore[i].fSoft/acc);
		}
	    if (bSqrtPhi) {
		/*
		** NOTE: The factor of 3.5 keeps this criterion in sync
		** with DensityStep. The nominal value of dEta for both
		** cases is then 0.02-0.03.
		*/
		double dtemp =
		    dEta*3.5*sqrt(dAccFac*fabs(pkd->pStore[i].fPot))/acc;
		if (dtemp < dT)
		    dT = dtemp;
		}
	    if (dT < pkd->pStore[i].dt) {
		pkd->pStore[i].dt = dT;
		mdlassert(pkd->mdl,dT>0);
	        }
	    }
	}
    }


void pkdDensityStep(PKD pkd,double dEta,double dRhoFac) {
    int i;
    double dT;
    
    for (i=0;i<pkdLocal(pkd);++i) {
	if (pkdIsActive(pkd,&(pkd->pStore[i]))) {
	    dT = dEta/sqrt(pkd->pStore[i].fDensity*dRhoFac);
	    if (dT < pkd->pStore[i].dt) {
		pkd->pStore[i].dt = dT;
		mdlassert(pkd->mdl,dT>0);
	        }
	    }
	}
    }

int pkdDtToRung(PKD pkd,int iRung,double dDelta,int iMaxRung,int bAll,int *nRungCount) {
    int i;
    int iTempRung;
    int iSteps;
    
    for (i=0;i<iMaxRung;++i) nRungCount[i] = 0;
    for(i=0;i<pkdLocal(pkd);++i) {
	if(pkd->pStore[i].iRung >= iRung) {
	    mdlassert(pkd->mdl,pkdIsActive(pkd,&(pkd->pStore[i])));
	    if(bAll) {          /* Assign all rungs at iRung and above */
		mdlassert(pkd->mdl,pkd->pStore[i].fSoft > 0);
		mdlassert(pkd->mdl,pkd->pStore[i].dt > 0);
		iSteps = dDelta/pkd->pStore[i].dt;
		/* insure that integer boundary goes
		   to the lower rung. */
		if(fmod(dDelta,pkd->pStore[i].dt) == 0.0)
		    iSteps--;
		iTempRung = iRung;
		if(iSteps < 0)
		    iSteps = 0;
		while(iSteps) {
		    ++iTempRung;
		    iSteps >>= 1;
		    }
		if(iTempRung >= iMaxRung) {
		    iTempRung = iMaxRung-1;
		    printf("Maximum Rung %d overshot for particle %d!\n",iMaxRung-1,pkd->pStore[i].iOrder);
		    }
		pkd->pStore[i].iRung = iTempRung;
		}
	    else {
		if(dDelta <= pkd->pStore[i].dt) {
		    pkd->pStore[i].iRung = iRung;
		    }
		else {
		    pkd->pStore[i].iRung = iRung+1;
		    }
                }
	    }
	/*
	** Now produce a count of particles in rungs.
	*/
	nRungCount[pkd->pStore[i].iRung] += 1;
	}
    iTempRung = iMaxRung-1;
    while (nRungCount[iTempRung] == 0 && iTempRung > 0) --iTempRung;
    return iTempRung;
    }

void pkdInitDt(PKD pkd,double dDelta) {
    int i;
    
    for(i=0;i<pkdLocal(pkd);++i) {
	if(pkdIsActive(pkd,&(pkd->pStore[i]))) {
	    pkd->pStore[i].dt = dDelta;
	    mdlassert(pkd->mdl,dDelta>0);
	    }
	}
    }
    

void pkdDeleteParticle(PKD pkd, PARTICLE *p) {
  
#ifdef PLANETS
  /* the following operation is an emergent treatment for 
     to-be-deleted VeryActive particles. */

    int j; 
       p->fMass = 0.0;
	 if(pkd->param.bGravStep){
	   p->dtGrav = 0.00001;
	 }
#ifdef HERMITE
        if(pkd->param.bAarsethStep){
	   p->dtGrav = 10000;
       }
#endif
#ifdef SYMBA
	p->drmin = 1000.0;
	p->n_VA=0;
#endif
       for (j=0;j<3;j++) {
	 p->r[j] = 100.0*p->iOrder;
	 p->v[j] = 0.0;     
#ifdef HERMITE
	 p->r0[j] = 100.0*p->iOrder;
	 p->v0[j] = 0.0;
	 p->a0[j] =  0.000001;
	 p->ad0[j] = 0.000001;
	 p->a[j] =  0.000001;
	 p->ad[j] = 0.000001;
#endif
       }

#endif

    p->iOrder = -2 - p->iOrder;
    TYPEClearACTIVE(p);
    TYPESet(p, TYPE_DELETED);

    }

void pkdNewParticle(PKD pkd, PARTICLE *p) {
    mdlassert(pkd->mdl,pkd->nLocal < pkd->nStore);
    pkd->pStore[pkd->nLocal] = *p;
    pkd->pStore[pkd->nLocal].iOrder = -1;
    pkd->nLocal++;
    }

void pkdColNParts(PKD pkd, int *pnNew, int *nDeltaGas, int *nDeltaDark,
		  int *nDeltaStar) {
    int pi, pj;
    int nNew;
    int ndGas;
    int ndDark;
    int ndStar;
    int newnLocal;
    PARTICLE *p;
    
    nNew = 0;
    ndGas = 0;
    ndDark = 0;
    ndStar = 0;
    newnLocal = pkdLocal(pkd);
    for(pi = 0, pj = 0; pi < pkdLocal(pkd); pi++) {
	if(pj < pi)
	    pkd->pStore[pj] = pkd->pStore[pi];
	p = &pkd->pStore[pi];
	if(p->iOrder == -1) {
	    ++pj;
	    ++nNew;
	    ++ndDark;
	    if(pkdIsActive(pkd,p))
		++pkd->nActive;
	    continue;
	    }
	else if(p->iOrder < -1){
	    --newnLocal;
	    p->iOrder = -2 - p->iOrder;
	    if(pkdIsGas(pkd, p))
		--ndGas;
	    else if(pkdIsDark(pkd, p))
		--ndDark;
	    else if(pkdIsStar(pkd, p))
		--ndStar;
	    else
		mdlassert(pkd->mdl,0);
	    if(pkdIsActive(pkd,p))
		--pkd->nActive;
	    }
	else {
	    ++pj;
	    }
	}

    *pnNew = nNew;
    *nDeltaGas = ndGas;
    *nDeltaDark = ndDark;
    *nDeltaStar = ndStar;
    pkd->nLocal = newnLocal;
    }

void
pkdNewOrder(PKD pkd,int nStart)
    {
    int pi;
    
    for(pi=0;pi<pkdLocal(pkd);pi++) {
	if(pkd->pStore[pi].iOrder == -1) {
	    pkd->pStore[pi].iOrder = nStart++;
	    }
	}
    }

void
pkdSetNParts(PKD pkd,int nGas,int nDark,int nStar,int nMaxOrderGas,
	     int nMaxOrderDark, double dSunMass) {
    pkd->nGas = nGas;
    pkd->nDark = nDark;
    pkd->nStar = nStar;
    pkd->nMaxOrderGas = nMaxOrderGas;
    pkd->nMaxOrderDark = nMaxOrderDark;
#ifdef PLANETS
    pkd->dSunMass;
#endif
    }


void
pkdSetRungVeryActive(PKD pkd, int iRung)
    {
    /* Remember, the first very active particle is at iRungVeryActive + 1 */
    pkd->uRungVeryActive = iRung;
    }

void pkdSetParticleTypes(PKD pkd)
    {
    PARTICLE *p;
    int i, iSetMask;
    
    for(i=0;i<pkdLocal(pkd);++i) {
	p = &pkd->pStore[i];
	iSetMask = 0;
	if (pkdIsGas (pkd,p)) iSetMask |= TYPE_GAS;
	if (pkdIsDark(pkd,p)) iSetMask |= TYPE_DARK;
	if (pkdIsStar(pkd,p)) iSetMask |= TYPE_STAR;

	TYPESet(p,iSetMask);
	}
    }

int pkdIsGas(PKD pkd,PARTICLE *p) {
    if (p->iOrder <= pkd->nMaxOrderGas) return 1;
    else return 0;
    }

int pkdIsDark(PKD pkd,PARTICLE *p) {
    if (p->iOrder > pkd->nMaxOrderGas && p->iOrder <= pkd->nMaxOrderDark)
	return 1;
    else return 0;
    }

int pkdIsStar(PKD pkd,PARTICLE *p) {
    if (p->iOrder > pkd->nMaxOrderDark) return 1;
    else return 0;
    }

#ifdef RELAXATION
void pkdInitRelaxation(PKD pkd)
    {
    int i,j;
    
    for(i=0;i<pkdLocal(pkd);++i) {
	pkd->pStore[i].fRelax = 0.0;
	}
    }

#endif /* RELAXATION */

#ifdef PLANETS
void
pkdReadSS(PKD pkd,char *pszFileName,int nStart,int nLocal)
{
	SSIO ssio;
	SSDATA data;
	PARTICLE *p;
	int i,j, iSetMask;

	pkd->nLocal = nLocal;
	pkd->nActive = nLocal;
	/*
	 ** General initialization (modeled after pkdReadTipsy()).
	 */
	for (i=0;i<nLocal;++i) {
		p = &pkd->pStore[i];
		TYPEClear(p);
		p->iRung = 0;
		p->fWeight = 1.0;
		p->fDensity = 0.0;		
		}
	/*
	 ** Seek past the header and up to nStart.
	 */
	if (ssioOpen(pszFileName,&ssio,SSIO_READ))
		mdlassert(pkd->mdl,0); /* unable to open ss file */
	if (ssioSetPos(&ssio,SSHEAD_SIZE + nStart*SSDATA_SIZE))
		mdlassert(pkd->mdl,0); /* unable to seek in ss file */
	/*
	 ** Read Stuff!
	 */
	for (i=0;i<nLocal;++i) {
		p = &pkd->pStore[i];
		p->iOrder = nStart + i;
		iSetMask = TYPE_DARK;
		if (ssioData(&ssio,&data))
			mdlassert(pkd->mdl,0); /* error during read in ss file */
		p->iOrgIdx = data.org_idx;
		p->fMass = data.mass;
		p->fSoft = data.radius;
		for (j=0;j<3;++j) p->r[j] = data.pos[j];
		for (j=0;j<3;++j) p->v[j] = data.vel[j];
		for (j=0;j<3;++j) p->w[j] = data.spin[j];
		p->iColor = data.color;
#ifdef NEED_VPRED
		for (j=0;j<3;++j) p->vPred[j] = p->v[j];
#endif
		TYPESet(p,iSetMask);
		}
	if (ssioClose(&ssio))
		mdlassert(pkd->mdl,0); /* unable to close ss file */
	}

void
pkdWriteSS(PKD pkd,char *pszFileName,int nStart)
{
	SSIO ssio;
	SSDATA data;
	PARTICLE *p;
	int i,j;

	/*
	 ** Seek past the header and up to nStart.
	 */
	if (ssioOpen(pszFileName,&ssio,SSIO_UPDATE))
		mdlassert(pkd->mdl,0); /* unable to open ss file */
	if (ssioSetPos(&ssio,SSHEAD_SIZE + nStart*SSDATA_SIZE))
		mdlassert(pkd->mdl,0); /* unable to seek in ss file */
	/* 
	 ** Write Stuff!
	 */
	for (i=0;i<pkdLocal(pkd);++i) {
		p = &pkd->pStore[i];
		if (!pkdIsDark(pkd,p))
			mdlassert(pkd->mdl,0); /* only dark particles allowed in ss file */
		data.org_idx = p->iOrgIdx;
		data.mass = p->fMass;
		data.radius = p->fSoft;
		for (j=0;j<3;++j) data.pos[j]  = p->r[j];
		for (j=0;j<3;++j) data.vel[j]  = p->v[j];
		for (j=0;j<3;++j) data.spin[j] = p->w[j];
		data.color = p->iColor;
		if (ssioData(&ssio,&data))
			mdlassert(pkd->mdl,0); /* unable to write in ss file */
		}
	if (ssioClose(&ssio))
		mdlassert(pkd->mdl,0); /* unable to close ss file */
	}

 void pkdSunIndirect(PKD pkd,double aSun[],double adSun[],int iFlag)
 {
     PARTICLE *p;
     FLOAT r2,r1i,r3i,r5i,rv;
     int i,j,n;
   
     for (j=0;j<3;++j)
	       	{                 
                 aSun[j] = 0;
		 adSun[j] = 0;
	       	}
     p = pkd->pStore;
     n = pkdLocal(pkd); 
     for (i=0;i<n;++i) {
       
        if (iFlag == 2){ 
	    if (pkdIsActive(pkd,&p[i])) continue;  /* inactive */
       }else if(iFlag == 1){
	    if (!pkdIsActive(pkd,&p[i])) continue; /* active */
       }

	 r2 = 0;
	 rv = 0;
	 for (j=0;j<3;++j)
	       	{ 
                 r2 += p[i].r[j]*p[i].r[j];
                 rv += p[i].v[j]*p[i].r[j];
	       	}
	 r1i = (r2 == 0 ? 0 : 1/sqrt(r2));
	 r3i = p[i].fMass*r1i*r1i*r1i;
         r5i = 3.0*rv*r3i*r1i*r1i; 
	 for (j=0;j<3;++j)
	       	{ 
                 aSun[j] += p[i].r[j]*r3i;
                 adSun[j] += p[i].v[j]*r3i-p[i].r[j]*r5i;
	       	}	
 }
 } 

void pkdGravSun(PKD pkd,double aSun[],double adSun[],double dSunMass)
{
	PARTICLE *p;
	double r2,v2,r1i,r3i;
        double aai,aa3i,idt2;
        double rv, r5i; 
	int i,j,n;
	double hx,hy,hz,h2,E,e,sum;

	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
	    if (pkdIsActive(pkd,&(p[i]))) {
			r2 = 0;
			v2 = 0;
                        rv = 0;
			for (j=0;j<3;++j)
			{ 
                        r2 += p[i].r[j]*p[i].r[j];
                        v2 += p[i].v[j]*p[i].v[j];
                        rv += p[i].v[j]*p[i].r[j];
			}
		
			r1i = (r2 == 0 ? 0 : 1/sqrt(r2)); /*gravity at origin = zero */
			p[i].fPot -= dSunMass*r1i;
			r3i = dSunMass*r1i*r1i*r1i;
                        r5i = 3.0*rv*r3i*r1i*r1i;
			/* time step is determined by semimajor axis, not the heliocentric distance*/ 
			if(pkd->param.bGravStep && p[i].fMass > 0){
			  /* E and h are normalized by the reduced mass */
			  sum = dSunMass + p[i].fMass;
			  hx =  p[i].r[1]*p[i].v[2] - p[i].r[2]*p[i].v[1];  
			  hy =  p[i].r[2]*p[i].v[0] - p[i].r[0]*p[i].v[2];  
			  hz =  p[i].r[0]*p[i].v[1] - p[i].r[1]*p[i].v[0];  
			  h2 = hx*hx + hy*hy + hz*hz;
			  E = 0.5*v2 - sum*r1i;			  
			  e = sqrt(1.0+2.0*E*h2/sum/sum);
			  aai = -2.0*E/sum/(1.0-e); 
			  aai = aai*aai*aai;
			  idt2 = sum*aai;
			  /*if (p[i].dtSun > p[i].dtGrav) p[i].dtGrav = p[i].dtSun;*/
			  if (idt2 > p[i].dtGrav) p[i].dtGrav = idt2;
			}		   
			for (j=0;j<3;++j) {	
#ifdef HERMITE
				  if(pkd->param.bHermite){
				    p[i].app[j] = p[i].a[j] - aSun[j]; /* perturbation force*/
				    p[i].adpp[j] = p[i].ad[j] - adSun[j];					
				    p[i].ad[j] -= (adSun[j] + p[i].v[j]*r3i-p[i].r[j]*r5i);
				  }
#endif
					p[i].a[j] -= (aSun[j] + p[i].r[j]*r3i);
					}				
			}
		}
	}


#ifdef SYMBA
void
pkdStepVeryActiveSymba(PKD pkd, double dStep, double dTime, double dDelta,
		       int iRung, int iKickRung, int iRungVeryActive,
		       int iAdjust, double diCrit2,
		       int *pnMaxRung, double dSunMass, int multiflag)
{
  int nRungCount[256];
  double dDriftFac;
  double ddDelta, ddStep;

/* get pointa of interacting opponents from thier iOrder's
 multiflag > 0 if close encounters with more than 2 particles exist*/
  if(iRung ==0) {
      pkd->nVeryActive = pkdSortVA(pkd, iRung);
      multiflag = pkdGetPointa(pkd);
      /*printf("Time %e, nVA %d,  multiflag %d \n",dTime, pkd->nVeryActive, multiflag);*/
  }
  if(iRung > 0){
    /* at iRung = 0 we just call pkdStepVeryActiveSymba three times*/  
    
    pkdActiveRung(pkd,iRung,1);
    if(iAdjust && (iRung < pkd->param.iMaxRung-1)) {       
      /* determine KickRung from current position */
	*pnMaxRung = pkdDrminToRungVA(pkd,iRung,pkd->param.iMaxRung,multiflag);	    
	if (pkd->param.bVDetails) {
	    printf("%*cAdjust, iRung: %d\n",2*iRung+2,' ',iRung);
	}
    }
    if(iAdjust == 0){
      /*if(iRung >= 2) pkdSortVA(pkd, iRung); maybe added later */
      pkdGravVA(pkd,iRung);	
      if (pkd->param.bVDetails) {
	    printf("%*cpkdGravVA, iRung: %d\n",2*iRung+2,' ',iRung);
	}
      
    }  
    /* following kick is for particles at iRung and iRung + 1 */     
    pkdKickVA(pkd,0.5*dDelta);  
    
    if (pkd->param.bVDetails) {
      printf("%*cVeryActive pkdKickOpen  at iRung: %d, 0.5*dDelta: %g\n",
	     2*iRung+2,' ',iRung,0.5*dDelta);
    } 
    
    /* if particles exist at the current iRung */
    pkdActiveRung(pkd,iRung,0); 
    pkdKeplerDrift(pkd,dDelta,dSunMass,1); /* 1 means VA */
     if (pkd->param.bVDetails) {
      printf("%*cVeryActive pkdKeplerDrift  at iRung: %d, 0.5*dDelta: %g\n",
	     2*iRung+2,' ',iRung,0.5*dDelta);
    }  
    /* if drmin druing drift is less than R_(iRung+1), 
       this interacting pair is sent to iRung + 1*/
     if(iRung < pkd->param.iMaxRung-1){
	 *pnMaxRung = pkdCheckDrminVA(pkd,iRung,multiflag,*pnMaxRung);
     }
  }	
  
  if (*pnMaxRung > iRung) {
	/*
	** Recurse.
	*/
    ddDelta = dDelta/3.0;
    ddStep = 1.0/pow(3.0, iRung);/* this is currently unnecessary */

    pkdStepVeryActiveSymba(pkd,dStep,dTime,ddDelta,iRung+1,iRung+1,iRungVeryActive,0,
			   diCrit2,pnMaxRung,dSunMass,multiflag);
    dStep += ddStep;
    dTime += ddDelta;
    pkdActiveRung(pkd,iRung,0);
    pkdStepVeryActiveSymba(pkd,dStep,dTime,ddDelta,iRung+1,iRung+1,iRungVeryActive,1,
			   diCrit2,pnMaxRung,dSunMass,multiflag);
    dStep += ddStep;
    dTime += ddDelta;
    pkdActiveRung(pkd,iRung,0);
    pkdStepVeryActiveSymba(pkd,dStep,dTime,ddDelta,iRung+1,iKickRung,iRungVeryActive,1,
			   diCrit2,pnMaxRung,dSunMass,multiflag);
    
    /* move time back to 1 step */
    dTime -= dDelta;			     	
  }
  if(iRung > 0){ 
    dTime += 0.5*dDelta;		
    
    if (pkd->param.bVDetails) {
      printf("%*cVeryActive pkdKickClose at iRung: %d, 0.5*dDelta: %g\n",
	     2*iRung+2,' ',iRung,0.5*dDelta);
    }
    
    /* Kick due to F_iRung for particles at the current or higher iRungs */		
    pkdActiveRung(pkd,iRung,1); 
    pkdGravVA(pkd,iRung);
    pkdKickVA(pkd,0.5*dDelta);
    
    if(pkd->param.bCollision && pkd->iCollisionflag) pkdDoCollisionVeryActive(pkd,dTime);  
    
  } 
}

int pkdSortVA(PKD pkd, int iRung){
    int i,j,n;
    PARTICLE *p = pkd->pStore;
    PARTICLE t;
    n = pkd->nLocal;
    int nSort = 0;
    
    /*
    ** Now start the partitioning of the particles.
    */
    if(iRung == 0){
	i = 0;
    }else{
	i = n-(pkd->nVeryActive);
    }
    j = n - 1;

    while (i <= j) {
	if ( p[i].iRung <= iRung ) ++i;
	else break;
    }
    while (i <= j) {
	if ( p[j].iRung > iRung ) --j;
	else break;
    }

    if (i < j) {
	t = p[i];
	p[i] = p[j];
	p[j] = t;
	while (1) {
	    while ((p[++i].iRung <= iRung));
		   while (p[--j].iRung > iRung);
		   if (i < j) {
		       t = p[i];
		       p[i] = p[j];
		       p[j] = t;
		   }
		   else break;
		   }
	}
    nSort = n - i;
    return(nSort);     
}


void pkdGravVA(PKD pkd, int iRung){
    int i, j, k, n, iStart;
    double x, y, z;
    double fac, d2,a2;
    PARTICLE *p;
    double r_k = 3.0/pow(2.08,iRung-1);
    double r_kk = r_k/2.08;
    double r_kkk = r_kk/2.08;
    double fourh2;

    assert(iRung >= 1);
    iStart = pkd->nLocal - pkd->nVeryActive;
    
    p = pkd->pStore;
    n = pkdLocal(pkd);
    
    /* reset */
    for(i=iStart;i<n;++i) {	
	if(pkdIsActive(pkd,&p[i])){
	    p[i].drmin = 1000.0;
	    for(k=0;k<3;k++){
		p[i].a_VA[k] = 0.0;
	    }
	}
    }
    
    for(i=iStart;i<n;++i) {	
	
	if(pkdIsActive(pkd,&p[i])){
	    for(k=0;k<p[i].n_VA;k++){
		j = p[i].i_VA[k];
		if(i<j){
		    x = p[i].r[0] - p[j].r[0];
		    y = p[i].r[1] - p[j].r[1];
		    z = p[i].r[2] - p[j].r[2];
		    d2 = x*x + y*y +z*z; 
		    a2 = sqrt(d2);
		    
		    if(pkd->param.bCollision){	
		      fourh2 = p[i].fSoft + p[j].fSoft;
		      if(a2 < fourh2){
			  pkd->iCollisionflag = 1;	 
			  p[i].iColflag = 1;
			  p[i].iOrderCol = p[j].iOrder;
			  p[i].dtCol = 1.0*p[i].iOrgIdx;
			  printf("dr = %e, r1+r2 = %e, pi = %i, pj = %i piRung %d iRung %d\n",
				 a2,fourh2,p[i].iOrgIdx,p[j].iOrgIdx, p[i].iRung, iRung);
			  printf("i %d j %d inVA %d, jnVA %d \n",i, j, p[i].n_VA, p[j].n_VA);
		      }
		    }

		    d2 *= a2;
		    a2 /= p[i].hill_VA[k]; /* distance normalized by hill */ 		
		    
		    if (a2 < r_kkk){		 
			fac = 0.0;		  
		    }else if (a2 < r_kk && a2 >r_kkk){
			fac = symfac*(1.0-a2/r_kk);
			fac *= fac*(2.0*fac -3.0);
			fac += 1.0;
		    }else if (a2 < r_k && a2 >r_kk){
			fac = symfac*(1.0-a2/r_k);
			fac *= -fac*(2.0*fac -3.0);
		    } else if (a2 > r_k){
			fac = 0.0;		  
		    }		    

		    p[i].drmin = (a2 < p[i].drmin)?a2:p[i].drmin; 
		    p[j].drmin = (a2 < p[j].drmin)?a2:p[j].drmin; 

		    d2 = fac/d2;
		    /*printf("iOrder %d %d, d2 %e iRung %d\n",
		      p[i].iOrder,p[j].iOrder,d2,iRung);*/
		    p[i].a_VA[0] -= x*d2*p[j].fMass;
		    p[i].a_VA[1] -= y*d2*p[j].fMass;
		    p[i].a_VA[2] -= z*d2*p[j].fMass;
		    p[j].a_VA[0] += x*d2*p[i].fMass;
		    p[j].a_VA[1] += y*d2*p[i].fMass;
		    p[j].a_VA[2] += z*d2*p[i].fMass;
		    
		}
	    }	   
	}
    }
}

int pkdCheckDrminVA(PKD pkd, int iRung, int multiflag, int nMaxRung){
  int i, j, k, n, iStart;
  double x, y, z;
  double d2;
  int iTempRung;
  PARTICLE *p;
  double r_k = 3.0/pow(2.08,iRung);
  
  assert(iRung >= 1);
  p = pkd->pStore;
  iStart = pkd->nLocal - pkd->nVeryActive;
  n = pkdLocal(pkd);

  for(i=iStart;i<n;++i) {	
    if(pkdIsActive(pkd,&p[i])){
      for(k=0;k<p[i].n_VA;k++){
	j = p[i].i_VA[k];
	if(i<j){
	  x = p[i].r[0] - p[j].r[0];
	  y = p[i].r[1] - p[j].r[1];
	  z = p[i].r[2] - p[j].r[2];
	  d2 = x*x + y*y +z*z; 
	  d2 = sqrt(d2);
	  d2 /= p[i].hill_VA[k]; /* distance normalized by hill */ 		
	  
	  /* d2 should be replaced by min.dist. during drift */

	  if(d2 < r_k){
	    p[i].iRung = iRung +1;
	    p[j].iRung = iRung +1;
	    nMaxRung = ((iRung+1) >= nMaxRung)?(iRung+1):nMaxRung;
	    /* retrive */
	    for(k=0;k<3;k++){
	      p[i].r[k] = p[i].rb[k];
	      p[i].v[k] = p[i].vb[k];
	      p[j].r[k] = p[j].rb[k];
	      p[j].v[k] = p[j].vb[k];
	    }
	  } 
	}
      }
    }
  }
  /* If a close encounter with more than 2 particles exists, we synchronize
     iRungs of interacting particles to the highest iRung in them */
  /* Following treatment is not sufficient if p[i].n_VA >= 3,
     in that case we should synchronize from larger p[i].n_VA to smaller */
  if(multiflag){
      for(i=iStart;i<n;++i) {
	  if(pkdIsActive(pkd,&p[i])){
	      if(p[i].n_VA >= 2){		    
		  for(k=0;k<p[i].n_VA;k++){		   
		      iTempRung = p[p[i].i_VA[k]].iRung;
		      if(iTempRung > p[i].iRung){
			  if(iTempRung != (iRung +1)){
			      printf("p_iOrder %d pi_n_VA %d pj_nVA %d iTempRung %d, iRung+1 %d \n",
				     p[i].iOrder,p[i].n_VA,p[p[i].i_VA[k]].n_VA,iTempRung,iRung +1);
			      printf("needs care for four body encounter!! \n");
			  }

			  assert(iTempRung == (iRung +1));
			  p[i].iRung = iRung + 1;
			  for(k=0;k<3;k++){
			      p[i].r[k] = p[i].rb[k];
			      p[i].v[k] = p[i].vb[k];	
			  }
		      } 
		  }
	      }
	  }
      }
  }
  return(nMaxRung);
}

void pkdKickVA(PKD pkd, double dt){  
  int i, k, iStart;
  PARTICLE *p;
  iStart = pkd->nLocal - pkd->nVeryActive;

  p = pkd->pStore;
  for(i=iStart;i<pkdLocal(pkd);++i) {	
    if(pkdIsActive(pkd,&p[i])){
      for(k=0;k<3;k++){
	p[i].v[k] += p[i].a_VA[k]*dt;
      }
    }
  }
}

int pkdDrminToRung(PKD pkd, int iRung, int iMaxRung, int *nRungCount) {
  int i, j, iTempRung;;
  PARTICLE *p ;
  
  /* R_k = 3.0/(2.08)^(k-1) with (k= 1,2, ...)*/ 
  /* note iMaxRung = pkd.param->iMaxRung */
  for (i=0;i<iMaxRung;++i) nRungCount[i] = 0; 
  p = pkd->pStore;
  for(i=0;i<pkdLocal(pkd);++i) {
    if(pkdIsActive(pkd,&p[i])){	  	  
      iTempRung = iRung;
      if(p[i].drmin > 3.0){
	iTempRung = 0;
      }else{
	iTempRung = floor(log(3.0/p[i].drmin)/log(2.08)) + 1;
      }
      if(iTempRung >= iMaxRung) {
	iTempRung = iMaxRung-1;
      }
      
      /* if min. dist. during drift is less than 3 Hill radius, 
	 set iRung = 1 */ 
      if(iTempRung == 0 && p[i].drmin2 < 3.0){	    
	iTempRung = 1;
      }
      
      /* retrive position and velocity of active particle to 
	 those before drift*/
      if(iTempRung > 0){
	for(j=0;j<3;j++){
	  p[i].r[j] = p[i].rb[j];
	  p[i].v[j] = p[i].vb[j];
	}
      }	
      /*
      ** Now produce a count of particles in rungs.
      */
      nRungCount[iTempRung] += 1;
      p[i].iRung = iTempRung;
      
      /*printf("iorder %d, drmin1 %e, drmin2 %e, \n",
	p[i].iOrder, p[i].drmin, p[i].drmin2);*/
    }
    
  }
  iTempRung = iMaxRung-1;
  while (nRungCount[iTempRung] == 0 && iTempRung > 0) --iTempRung;
  return iTempRung;
}

int pkdGetPointa(PKD pkd){
  int i, j, k, n, iStart, iOrderTemp, n_VA;
  int multiflag = 0;
  int iTempRung;
  int bRung = 0;
  PARTICLE *p;
  char c;
  int iMaxRung = pkd->param.iMaxRung; 
  int nRungCount[iMaxRung-1];/* 1 to iMaxRung-1 */

  p = pkd->pStore;
  iStart = pkd->nLocal - pkd->nVeryActive;
  n = pkdLocal(pkd);

/*  for(i=0;i<n;++i) {
      printf("iOrder %d, iRung %d, n_VA %d \n",p[i].iOrder,p[i].iRung,p[i].n_VA);
      }*/
  if(bRung){
  for (i=1;i<iMaxRung;++i) nRungCount[i] = 0; 
  }

  for(i=iStart;i<n;++i) {
      assert(pkdIsActive(pkd,&p[i]));
      if(bRung) nRungCount[p[i].iRung] += 1;
      n_VA = p[i].n_VA;
      assert(n_VA >= 1);
      
      if(n_VA >= 2) multiflag ++;
	        
      for(k=0;k<n_VA;k++){
	  iOrderTemp = p[i].iOrder_VA[k];
	  
	  for(j=iStart;j<n;++j) {
	      if(i==j)continue;
	      if(p[j].iOrder == iOrderTemp){
		  p[i].i_VA[k] = j;
		  break;
	      }
	  }
      }    	
  }
  
  if(bRung){
  for (i=1;i<iMaxRung;++i){
      if (nRungCount[i] == 0) continue;
      else c = ' ';
	    printf(" %c rung:%d %d\n",c,i,nRungCount[i]);
	    }
	printf("\n");
	}
  
  /* If a close encounter with more than 2 particles exists, we synchronize
     iRungs of interacting particles to the highest iRung in them */
  if(multiflag){
      for(i=iStart;i<n;++i) {
	  if(pkdIsActive(pkd,&p[i])){
	      if(p[i].n_VA >= 2){		    
		  for(k=0;k<p[i].n_VA;k++){
		      j = p[i].i_VA[k];
		      iTempRung = p[j].iRung;
		      p[i].iRung = (iTempRung >= p[i].iRung)?iTempRung:p[i].iRung;
		      p[j].iRung = p[i].iRung;
		  }	
	      }
	  }
      }
  }
  return(multiflag);
}

int pkdDrminToRungVA(PKD pkd, int iRung, int iMaxRung, int multiflag) {
  int i, j, k, iTempRung, iStart;
  int nMaxRung = 0; 
  PARTICLE *p;

  /* note iMaxRung = pkd.param->iMaxRung -1 */

  int bRung = 0;
  char c;
  int nRungCount[iMaxRung-1];/* 1 to iMaxRung-1*/
  if(bRung){
  for (i=1;i<iMaxRung;++i) nRungCount[i] = 0; 
  }
  
  iStart = pkd->nLocal - pkd->nVeryActive;
  
  p = pkd->pStore;
  for(i=iStart;i<pkdLocal(pkd);++i) {   
      if(pkdIsActive(pkd,&p[i])){	  	       	    
	/*assert(p[i].n_VA >= 1); n_VA might be 0 after collision */	     
	  iTempRung = floor(log(3.0/p[i].drmin)/log(2.08)) + 1;     	  
	  /*if(iTempRung >= iMaxRung){
	      printf("iRung %d for particle %d larger than Max iRung %d\n",
		     iTempRung, p[i].iOrder,iMaxRung-1);
	      iTempRung = iMaxRung;
	      }*/
	  iTempRung = (iTempRung >= iMaxRung)?(iMaxRung-1):iTempRung;
	  /* iTempRung = (iTempRung >= 0)?iTempRung:0; 
	     p[i].iKickRung = iTempRung; */
      
	  iTempRung = (iTempRung >= iRung)?iTempRung:iRung; 
	  p[i].iRung = iTempRung;
	  nMaxRung = (iTempRung >= nMaxRung)?iTempRung:nMaxRung; 
      }
      if(bRung) nRungCount[p[i].iRung] += 1;
  }
  
  if(bRung){
      printf("nVA %d iRung %d\n",pkd->nVeryActive, iRung);
      for (i=1;i<iMaxRung;++i){
	  if (nRungCount[i] == 0) continue;
	  else c = ' ';
	  printf(" %c rung:%d %d\n",c,i,nRungCount[i]);
      }
      printf("\n");
  }
  /* If a close encounter with more than 2 particles exists, we synchronize
     iRungs of interacting particles to the highest iRung in them */
  if(multiflag){
      for(i=iStart;i<pkdLocal(pkd);++i) {
	  if(pkdIsActive(pkd,&p[i])){
	      if(p[i].n_VA >= 2){	    
		  for(k=0;k<p[i].n_VA;k++){	
		      j = p[i].i_VA[k];
		      iTempRung = p[j].iRung;
		      p[i].iRung = (iTempRung >= p[i].iRung)?iTempRung:p[i].iRung;
		      p[j].iRung = p[i].iRung;
		  }
	      }
	  }
      }
  }
  return(nMaxRung);
}


void pkdMomSun(PKD pkd,double momSun[]){
  PARTICLE *p;
  int i,j,n;
  
  for (j=0;j<3;++j){                 
    momSun[j] = 0.0;
  }
  p = pkd->pStore;
  n = pkdLocal(pkd); 
  for (i=0;i<n;++i) {
    for (j=0;j<3;++j){ 
      momSun[j] -= p[i].fMass*p[i].v[j];     
    }   
  }   
} 

void pkdDriftSun(PKD pkd,double vSun[],double dt){
  PARTICLE *p;
  int i,j,n;
  
  p = pkd->pStore;
  n = pkdLocal(pkd); 
  for (i=0;i<n;++i) {
    for (j=0;j<3;++j){ 
      p[i].r[j] += vSun[j]*dt;     
    } 
  }
} 

void pkdKeplerDrift(PKD pkd,double dt,double mu, int tag_VA) {
  /* 
   * Use the f and g functions to advance an unperturbed orbit.
   */
  PARTICLE *p;	
  int  i,j,n, iStart;
  int iflg = 0;
  double dttmp;
  mdlDiag(pkd->mdl, "Into pkdKeplerDrift \n");
  p = pkd->pStore;
  n = pkdLocal(pkd);
  
  if(tag_VA){
    iStart = pkd->nLocal - pkd->nVeryActive;
  }else{
    iStart = 0; 
  }
  
  for (i=iStart;i<n;++i) {	
    if (pkdIsActive(pkd,&p[i])) {
      
      /* copy r and v before drift */ 
      for (j=0;j<3;++j) {	
        p[i].rb[j] = p[i].r[j];
        p[i].vb[j] = p[i].v[j];
      }
      
      /*printf("before: iorder = %d, (x,y,z)= (%e,%e,%e), (vx,vy,vz)= (%e,%e,%e),ms =%e \n",
	p[i].iOrder, p[i].r[0],p[i].r[1],p[i].r[2],
	p[i].v[0],p[i].v[1],p[i].v[2], mu);*/
      
      iflg = drift_dan(mu,p[i].r,p[i].v,dt); /* see kepler.c */
      
      /* printf("exit drift_dan, iflg = %d, (x,y,z)= (%e,%e,%e),(vx,vy,vz)= (%e,%e,%e), ms =%e \n", 
	 iflg, p[i].r[0],p[i].r[1],p[i].r[2],p[i].v[0],p[i].v[1],p[i].v[2],mu);*/
      
      if(iflg != 0){
	dttmp = 0.1*dt;
	for (j=0;j<10;++j) {	
	  iflg = drift_dan(mu,p[i].r,p[i].v,dttmp);
	  assert(iflg == 0);
	}
      }    
    }
  }
}	

#endif /* SYMBA */
#endif /* PLANETS */


