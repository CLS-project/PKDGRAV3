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

#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h> /* for unlink() */
#endif
#ifdef _MSC_VER
#define mkdir _mkdir
#define unlink _unlink
#endif
#include <stddef.h>
#include <string.h>
#include <ctype.h>
#ifdef HAVE_INTTYPES_H
#include <inttypes.h>
#else
#define PRIu64 "llu"
#endif
#include <limits.h>
#include <stdarg.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <assert.h>
#include <time.h>
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#ifdef HAVE_SYS_STAT_H
#include <sys/stat.h>
#endif
#include <math.h>
#if defined(HAVE_WORDEXP) && defined(HAVE_WORDFREE)
#include <wordexp.h>
#elif defined(HAVE_GLOB) && defined(HAVE_GLOBFREE)
#include <glob.h>
#endif
#include <sys/stat.h>

#ifdef HAVE_SYS_PARAM_H
#include <sys/param.h> /* for MAXHOSTNAMELEN, if available */
#endif

#include "master.h"
#include "illinois.h"
#include "tipsydefs.h"
#include "outtype.h"
#include "smoothfcn.h"
#include "fio.h"

#define LOCKFILE ".lockfile"	/* for safety lock */
#define STOPFILE "STOP"			/* for user interrupt */
#define CHECKFILE "CHECKPOINT"		/* for user interrupt */

void msrprintf(MSR msr, const char *Format, ... ) {
    va_list ap;
    if (msr->param.bVDetails) {
	va_start(ap,Format);
	vprintf(Format,ap);
	va_end(ap);
	}
    }

#ifdef _MSC_VER
double msrTime() {
    FILETIME ft;
    uint64_t clock;
    GetSystemTimeAsFileTime(&ft);
    clock = ft.dwHighDateTime;
    clock <<= 32;
    clock |= ft.dwLowDateTime;
    /* clock is in 100 nano-second units */
    return clock / 10000000.0;
    }
#else
double msrTime() {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return (tv.tv_sec+(tv.tv_usec*1e-6));
    }
#endif

void _msrLeader(void) {
    puts("pkdgrav"PACKAGE_VERSION" Joachim Stadel & Doug Potter Sept 2015");
    puts("USAGE: pkdgrav3 [SETTINGS | FLAGS] [SIM_FILE]");
    puts("SIM_FILE: Configuration file of a particular simulation, which");
    puts("          includes desired settings and relevant input and");
    puts("          output files. Settings specified in this file override");
    puts("          the default settings.");
    puts("SETTINGS");
    puts("or FLAGS: Command line settings or flags for a simulation which");
    puts("          will override any defaults and any settings or flags");
    puts("          specified in the SIM_FILE.");
    }


void _msrTrailer(void) {
    puts("(see man page for more information)");
    }


void _msrExit(MSR msr,int status) {
    MDL mdl=msr->mdl;

    msrFinish(msr);
    exit(status);
    }

char *_BuildName(MSR msr,char *achFile,int iStep,char *defaultPath) {
    char achOutPath[256], *p;
    int n;

    if ( defaultPath[0] ) {
	strcpy( achOutPath, defaultPath );
	p = strstr( achOutPath, "&N" );
	if ( p ) {
	    n = p - achOutPath;
	    strcpy( p, msrOutName(msr) );
	    strcat( p+2, defaultPath + n + 2 );
	    }
	else {
	    n = strlen(achOutPath);
	    if ( !n || achOutPath[n-1]!='/' )
		achOutPath[n++] = '/';
	    strcpy(achOutPath+n,msrOutName(msr));
	    }
	}
    else {
	strcpy(achOutPath,msrOutName(msr));
	}

    p = strstr( achOutPath, "&S" );
    if ( p ) {
	n = p - achOutPath;
	strncpy( achFile, achOutPath, n );
	sprintf( achFile+n, "%05d", iStep );
	strcat( achFile+n, p+2 );
	}
    else {
	char achDigitMask[20];
	sprintf(achDigitMask,"%%s.%%0%ii",msr->param.nDigits);
	sprintf(achFile,achDigitMask,msrOutName(msr),iStep);
	}
    for(p=achFile+1; *p; ++p) {
	if ( *p == '/') {
	    *p = 0;
	    mkdir(achFile,0755);
	    *p = '/';
	    }
	}
    return achFile;
    }

void _msrMakePath(const char *dir,const char *base,char *path) {
    /*
    ** Prepends "dir" to "base" and returns the result in "path". It is the
    ** caller's responsibility to ensure enough memory has been allocated
    ** for "path".
    */

    if (!path) return;
    path[0] = 0;
    if (dir&&dir[0]) {
	strcat(path,dir);
	strcat(path,"/");
	}
    if (!base) return;
    strcat(path,base);
    }

static uint64_t getMemoryModel(MSR msr) {
    uint64_t mMemoryModel = 0;
    /*
    ** Figure out what memory models are in effect.  Configuration flags
    ** can be used to request a specific model, but certain operations
    ** will force these flags to be on.
    */
    if (msr->param.bFindGroups) mMemoryModel |= PKD_MODEL_GROUPS|PKD_MODEL_VELOCITY;
    if (msrDoGravity(msr)) {
	mMemoryModel |= PKD_MODEL_VELOCITY|PKD_MODEL_NODE_MOMENT;
	if (!msr->param.bNewKDK) mMemoryModel |= PKD_MODEL_ACCELERATION;
	}
    if (msr->param.bDoDensity)       mMemoryModel |= PKD_MODEL_DENSITY;
    if (msr->param.bMemUnordered)    mMemoryModel |= PKD_MODEL_UNORDERED;
    if (msr->param.bMemParticleID)   mMemoryModel |= PKD_MODEL_PARTICLE_ID;
    if (msr->param.bTraceRelaxation) mMemoryModel |= PKD_MODEL_RELAXATION;
    if (msr->param.bMemAcceleration || msr->param.bDoAccOutput) mMemoryModel |= PKD_MODEL_ACCELERATION;
    if (msr->param.bMemVelocity)     mMemoryModel |= PKD_MODEL_VELOCITY;
    if (msr->param.bMemPotential || msr->param.bDoPotOutput)    mMemoryModel |= PKD_MODEL_POTENTIAL;
    if (msr->param.bFindHopGroups)   mMemoryModel |= PKD_MODEL_GROUPS | PKD_MODEL_DENSITY | PKD_MODEL_BALL;
    if (msr->param.bMemGroups)       mMemoryModel |= PKD_MODEL_GROUPS;
    if (msr->param.bMemMass)         mMemoryModel |= PKD_MODEL_MASS;
    if (msr->param.bMemSoft)         mMemoryModel |= PKD_MODEL_SOFTENING;
    if (msr->param.bMemRelaxation)   mMemoryModel |= PKD_MODEL_RELAXATION;
    if (msr->param.bMemVelSmooth)    mMemoryModel |= PKD_MODEL_VELSMOOTH;

    if (msr->param.bMemNodeAcceleration) mMemoryModel |= PKD_MODEL_NODE_ACCEL;
    if (msr->param.bMemNodeVelocity)     mMemoryModel |= PKD_MODEL_NODE_VEL;
    if (msr->param.bMemNodeMoment)       mMemoryModel |= PKD_MODEL_NODE_MOMENT;
    if (msr->param.bMemNodeSphBounds)    mMemoryModel |= PKD_MODEL_NODE_SPHBNDS;
    if (msr->param.bDoGas)               mMemoryModel |= (PKD_MODEL_SPH | PKD_MODEL_NODE_SPHBNDS | PKD_MODEL_ACCELERATION);

    if (msr->param.bMemNodeBnd)          mMemoryModel |= PKD_MODEL_NODE_BND;
    if (msr->param.bMemNodeVBnd)         mMemoryModel |= PKD_MODEL_NODE_VBND;

    return mMemoryModel;
    }

void msrInitializePStore(MSR msr, uint64_t *nSpecies) {
    struct inInitializePStore ps;
    double dStorageAmount = (1.0+msr->param.dExtraStore);
    int i;
    for( i=0; i<FIO_SPECIES_LAST; ++i) ps.nSpecies[i] = nSpecies[i];
    /* If we plan to replicate the particles later then reserve space */
    if (msr->param.nInflateReps) {
	i = msr->param.nInflateReps + 1;
	dStorageAmount *= i*i*i;
	}
    ps.nStore = ceil( dStorageAmount * ps.nSpecies[FIO_SPECIES_ALL] / mdlThreads(msr->mdl));
    ps.nTreeBitsLo = msr->param.nTreeBitsLo;
    ps.nTreeBitsHi = msr->param.nTreeBitsHi;
    ps.iCacheSize  = msr->param.iCacheSize;
    ps.iWorkQueueSize  = msr->param.iWorkQueueSize;
    ps.iCUDAQueueSize  = msr->param.iCUDAQueueSize;
    ps.fPeriod[0] = msr->param.dxPeriod;
    ps.fPeriod[1] = msr->param.dyPeriod;
    ps.fPeriod[2] = msr->param.dzPeriod;
    ps.mMemoryModel = getMemoryModel(msr) | PKD_MODEL_VELOCITY;
    ps.bLightCone  = msr->param.bLightCone;
    ps.bLightConeParticles  = msr->param.bLightConeParticles;    

#define SHOW(m) ((ps.mMemoryModel&PKD_MODEL_##m)?" " #m:"")
       printf("Memory Models:%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s\n", 
#ifdef INTEGER_POSITION
	   " INTEGER_POSITION",
#else
	   " DOUBLE_POSITION",
#endif
	   SHOW(UNORDERED),SHOW(VELOCITY),SHOW(ACCELERATION),SHOW(POTENTIAL),
	   SHOW(GROUPS),SHOW(RELAXATION),SHOW(MASS),SHOW(DENSITY),
	   SHOW(BALL),SHOW(SOFTENING),SHOW(VELSMOOTH),SHOW(SPH),
	   SHOW(STAR),SHOW(PARTICLE_ID));
#undef SHOW
    ps.nMinEphemeral = 0;
    ps.nMinTotalStore = 0;

    /* Various features require more or less ephemeral storage */
    ps.nEphemeralBytes = 0;
    if (msr->param.iFofInterval   && ps.nEphemeralBytes < 4) ps.nEphemeralBytes = 4;
    if (msr->param.bFindHopGroups && ps.nEphemeralBytes < 8) ps.nEphemeralBytes = 8;
    if (msr->param.iPkInterval    && ps.nEphemeralBytes < 4) ps.nEphemeralBytes = 4;
    if (msr->param.bGravStep      && ps.nEphemeralBytes < 8) ps.nEphemeralBytes = 8;
    if (msr->param.bDoGas         && ps.nEphemeralBytes < 8) ps.nEphemeralBytes = 8;
    if (msr->param.bDoDensity     && ps.nEphemeralBytes < 12) ps.nEphemeralBytes = 12;
#ifdef MDL_FFTW
    if (msr->param.nGridPk>0) {
	struct inGetFFTMaxSizes inFFTSizes;
	struct outGetFFTMaxSizes outFFTSizes;
	int n;
	inFFTSizes.nx = inFFTSizes.ny = inFFTSizes.nz = msr->param.nGridPk;
	pstGetFFTMaxSizes(msr->pst,&inFFTSizes,sizeof(inFFTSizes),&outFFTSizes,&n);
	/* The new MeasurePk requires two FFTs to eliminate aliasing */
	ps.nMinEphemeral = (msr->param.bPkInterlace?2:1)*outFFTSizes.nMaxLocal*sizeof(FFTW3(real));
	}
    /* 
     * Add some ephemeral memory (if needed) for the linGrid.
     * 3 grids are stored : forceX, forceY, forceZ
     */
    if (strlen(msr->param.csm->val.classData.achLinSpecies)){
	struct inGetFFTMaxSizes inFFTSizes;
	struct outGetFFTMaxSizes outFFTSizes;
	int n;

	inFFTSizes.nx = inFFTSizes.ny = inFFTSizes.nz = msr->param.nGridLin;
	pstGetFFTMaxSizes(msr->pst, &inFFTSizes,sizeof(inFFTSizes), &outFFTSizes, &n);

	if (ps.nMinEphemeral < 3*outFFTSizes.nMaxLocal*sizeof(FFTW3(real)))
	    ps.nMinEphemeral = 3*outFFTSizes.nMaxLocal*sizeof(FFTW3(real));
	}

    if (msr->param.nGrid>0) {
	struct inGetFFTMaxSizes inFFTSizes;
	struct outGetFFTMaxSizes outFFTSizes;
	int n;
	inFFTSizes.nx = inFFTSizes.ny = inFFTSizes.nz = msr->param.nGrid;
	pstGetFFTMaxSizes(msr->pst,&inFFTSizes,sizeof(inFFTSizes),&outFFTSizes,&n);
	ps.nMinTotalStore = 10*outFFTSizes.nMaxLocal*sizeof(FFTW3(real));
	}
#endif
       pstInitializePStore(msr->pst,&ps,sizeof(ps),NULL,NULL);
       PKD pkd = msr->pst->plcl->pkd;
       printf("Allocated %lu MB for particle store on each processor.\n",
	   pkdParticleMemory(pkd)/(1024*1024));
       printf("Particles: %lu bytes (persistent) + %d bytes (ephemeral), Nodes: %lu bytes\n",
	   pkdParticleSize(pkd),ps.nEphemeralBytes,pkdNodeSize(pkd));
    }

static char *formatKey(char *buf,char *fmt,int i) {
    sprintf(buf,fmt,i);
    return buf;
    }

static int readParametersClasses(MSR msr,FILE *fp) {
    int i;
    char key[50];

    if (fscanf(fp,"nClasses=%d\n",&msr->nCheckpointClasses) != 1) return 1;
    for(i=0; i<msr->nCheckpointClasses; ++i) {
	if (fscanf(fp,formatKey(key,"fMass%d=%%g\n",i),&msr->aCheckpointClasses[i].fMass) != 1) return 1;
	if (fscanf(fp,formatKey(key,"fSoft%d=%%g\n",i),&msr->aCheckpointClasses[i].fSoft) != 1) return 1;
	if (fscanf(fp,formatKey(key,"eSpecies%d=%%d\n",i),&msr->aCheckpointClasses[i].eSpecies) != 1) return 1;
	}

    return 0;
    }

static int parseString(char *s, int size, const char *match) {
    char *p;
    if (match) {
	if ( (p = strchr(s,'=')) == NULL) return 0;
	*p++ = 0;
	if (strcmp(s,match)!=0) return 0;
	}
    else p = s;
    if (*p++ != '"') return 0;
    size_t n = strlen(p);
    while(n && isspace(p[n-1])) --n;
    if (n==0) return 0;
    if (p[n-1] != '"') return 0;
    p[n-1] = 0;
    memmove(s,p,n);
    return 1;
    }

static int scanString(char *s, int size, FILE *fp, const char *match) {
    if (fgets(s,size,fp)==NULL) return 0;
    return parseString(s,size,match);
    }

/*
** Read parameters from the checkpoint file. This is extremely STRICT;
** even an extra space will cause problems. This is by design.
*/
int readParameters(MSR msr,const char *fileName) {
    if (fileName == NULL) return 0;
    FILE *fp = fopen(fileName,"r");
    PRM_NODE *pn = NULL;
    char achBuffer[500];
    char achFormat[50];
    char *pScan, *p;
    int bError = 0;
    int n;

    if (fp==NULL) return 0;

    if (scanString(achBuffer,sizeof(achBuffer),fp,"VERSION")!=1 || strcmp(achBuffer,PACKAGE_VERSION)!=0) bError = 1;
    else if (fscanf(fp,"iStep=%d\n",&msr->iCheckpointStep)!=1) bError=1;
    else if (fscanf(fp,"dTime=%lg\n",&msr->dCheckpointTime)!=1) bError=1;
    else if (fscanf(fp,"dEcosmo=%lg\n",&msr->dEcosmo)!=1) bError=1;
    else if (fscanf(fp,"dTimeOld=%lg\n",&msr->dTimeOld)!=1) bError=1;
    else if (fscanf(fp,"dUOld=%lg\n",&msr->dUOld)!=1) bError=1;
    else if (fscanf(fp,"nSpecies0=%"PRIu64"\n",&msr->N)!=1) bError=1;
    else if (fscanf(fp,"nSpecies1=%"PRIu64"\n",&msr->nDark)!=1) bError=1;
    else if (fscanf(fp,"nSpecies2=%"PRIu64"\n",&msr->nGas)!=1) bError=1;
    else if (fscanf(fp,"nSpecies3=%"PRIu64"\n",&msr->nStar)!=1) bError=1;
    else if (readParametersClasses(msr,fp)) bError=1;
    else if (fscanf(fp,"nCheckpointFiles=%d\n",&msr->nCheckpointThreads)!=1) bError=1;
    else if (scanString(msr->achCheckpointName,sizeof(msr->achCheckpointName),fp,"achCheckpointName")==0) bError = 1;
    else for( pn=msr->prm->pnHead; pn!=NULL; pn=pn->pnNext ) {
	if (fgets(achBuffer,sizeof(achBuffer),fp)==NULL) bError=1;
	else {
	    pScan = achBuffer;
	    if (*pScan=='#') ++pScan;
	    else  pn->bFile = 1;
	    sprintf(achFormat,"%s=",pn->pszName);
	    if (strncmp(pScan,achFormat,strlen(achFormat))==0) {
	    	pScan += strlen(achFormat);
		if (pn->pCount) {
		    if (*pScan++ != '[') break;
		    n = strlen(pScan) - 2;
		    if (n<0 || pScan[n]!=']') break;
		    pScan[n] = 0;
		    if (n==0) {
		    	*pn->pCount = 0;
			continue;
			}
		    pScan = strtok(pScan,",");
		    assert(pScan);
		    }
		n = 0;
		while(pScan) {
		    switch (pn->iType) {
		    case 0:
		    case 1:
			assert(pn->iSize == sizeof(int));
			if (sscanf(pScan,"%d\n",(int *)pn->pValue+n)!=1) bError=1;
			break;
		    case 2:
			assert(pn->iSize == sizeof(double));
			if (sscanf(pScan,"%lg\n",(double *)pn->pValue+n)!=1) bError=1;
			break;
		    case 3:
			if (parseString(pScan,256,NULL)!=1) bError = 1;
			else strcpy(pn->pValue,pScan);
			break;
		    case 4:
			sprintf(achFormat,"%s=%%llu\n",pn->pszName);
			if (sscanf(pScan,"%"PRIu64,(uint64_t *)pn->pValue+n)!=1) bError=1;
			break;
			}
		    ++n;
		    if (pn->pCount) pScan = strtok(NULL,",");
		    else pScan = NULL;
		    }
		}
		else bError = 1;
		if (pn->pCount) *pn->pCount = n;
	    }
	if (bError) break;
	}
    fclose(fp);
    return (pn==NULL && bError==0);
    }

double msrRestore(MSR msr) {
    uint64_t nSpecies[FIO_SPECIES_LAST];
    struct inRestore restore;
    int i;
    double sec,dsec;
    BND bnd;

    if (mdlThreads(msr->mdl) != msr->nCheckpointThreads) {
	fprintf(stderr,"ERROR: You must restart a checkpoint with the same number of threads\n");
	fprintf(stderr,"       nThreads=%d, nCheckpointThreads=%d\n",mdlThreads(msr->mdl),msr->nCheckpointThreads);
	fprintf(stderr,"       RESTART WITH %d THREADS\n",msr->nCheckpointThreads);
	_msrExit(msr,1);
	}

    if (msr->param.bVStart)
	printf("Restoring from checkpoint\n");
    sec = msrTime();

    msr->nMaxOrder = msr->N;

    for( i=0; i<FIO_SPECIES_LAST; ++i) nSpecies[i] = 0;
    nSpecies[FIO_SPECIES_ALL]  = msr->N;
    nSpecies[FIO_SPECIES_SPH]  = msr->nGas;
    nSpecies[FIO_SPECIES_DARK] = msr->nDark;
    nSpecies[FIO_SPECIES_STAR] = msr->nStar;
    msrInitializePStore(msr,nSpecies);

    restore.nProcessors = msr->param.bParaRead==0?1:(msr->param.nParaRead<=1 ? msr->nThreads:msr->param.nParaRead);
    strcpy(restore.achInFile,msr->achCheckpointName);
    pstRestore(msr->pst,&restore,sizeof(restore),NULL,NULL);
    pstSetClasses(msr->pst,msr->aCheckpointClasses,msr->nCheckpointClasses*sizeof(PARTCLASS),NULL,NULL);
    pstCalcBound(msr->pst,NULL,0,&bnd,NULL);
    msrCountRungs(msr,NULL);

    dsec = msrTime() - sec;
    PKD pkd = msr->pst->plcl->pkd;
    double dExp = csmTime2Exp(msr->param.csm,msr->dCheckpointTime);
    msrprintf(msr,"Checkpoint Restart Complete @ a=%g, Wallclock: %f secs\n\n",dExp,dsec);

    /* We can indicate that the DD was already done at rung 0 */
    msr->iLastRungRT = 0;
    msr->iLastRungDD = 0;

    return msr->dCheckpointTime;
    }

static void writeParameters(MSR msr,const char *baseName,int iStep,double dTime) {
    PRM_NODE *pn;
    char *p, achOutName[PST_FILENAME_SIZE];
    char szNumber[30];
    double v;
    uint64_t nSpecies[FIO_SPECIES_LAST];
    int i;
    int nBytes;

    pstGetClasses(msr->pst,NULL,0,msr->aCheckpointClasses,&nBytes);
    msr->nCheckpointClasses = nBytes / sizeof(PARTCLASS);
    assert(msr->nCheckpointClasses*sizeof(PARTCLASS)==nBytes);

    strcpy( achOutName, baseName );
    p = strstr( achOutName, "&I" );
    if ( p ) {
	int n = p - achOutName;
	strcpy( p, "chk" );
	strcat( p, baseName + n + 2 );
	}
    else {
	strcat(achOutName,".chk");
	}

    FILE *fp = fopen(achOutName,"w");
    if (fp==NULL) {perror(achOutName); abort(); }

    for(i=0; i<FIO_SPECIES_LAST; ++i) nSpecies[i] = i;
    nSpecies[FIO_SPECIES_ALL]  = msr->N;
    nSpecies[FIO_SPECIES_SPH]  = msr->nGas;
    nSpecies[FIO_SPECIES_DARK] = msr->nDark;
    nSpecies[FIO_SPECIES_STAR] = msr->nStar;

    fprintf(fp,"VERSION=\"%s\"\n",PACKAGE_VERSION);
    fprintf(fp,"iStep=%d\n",iStep);
    fprintf(fp,"dTime=%.17g\n",dTime);
    fprintf(fp,"dEcosmo=%.17g\n",msr->dEcosmo);
    fprintf(fp,"dTimeOld=%.17g\n",msr->dTimeOld);
    fprintf(fp,"dUOld=%.17g\n",msr->dUOld);
    for(i=0; i<FIO_SPECIES_LAST; ++i)
	fprintf(fp,"nSpecies%d=%"PRIu64"\n",i,nSpecies[i]);
    fprintf(fp,"nClasses=%d\n",msr->nCheckpointClasses);
    for(i=0; i<msr->nCheckpointClasses; ++i) {
	fprintf(fp,"fMass%d=%.17g\n",i,msr->aCheckpointClasses[i].fMass);
	fprintf(fp,"fSoft%d=%.17g\n",i,msr->aCheckpointClasses[i].fSoft);
	fprintf(fp,"eSpecies%d=%d\n",i,msr->aCheckpointClasses[i].eSpecies);
	}
    fprintf(fp,"nCheckpointFiles=%d\n",mdlThreads(msr->mdl));
    fprintf(fp,"achCheckpointName=\"%s\"\n",baseName);

    for( pn=msr->prm->pnHead; pn!=NULL; pn=pn->pnNext ) {
    	void *pValue = pn->pValue;
    	int nCount = 1;
	fprintf(fp,"%s%s=",pn->bArg|pn->bFile?"":"#",pn->pszName);
	if (pn->pCount) {
	    fprintf(fp,"[");
	    nCount = *pn->pCount;
	    }
	for(i=0; i<nCount; ++i) {
	    if (i) fprintf(fp, ",");
	    switch (pn->iType) {
	    case 0:
	    case 1:
		assert(pn->iSize == sizeof(int));
		fprintf(fp,"%d",((int *)pn->pValue)[i]);
		break;
	    case 2:
		assert(pn->iSize == sizeof(double));
		/* This is purely cosmetic so 0.15 doesn't turn into 0.14999999... */
		sprintf(szNumber,"%.16g",((double *)pn->pValue)[i]);
		sscanf(szNumber,"%le",&v);
		if ( v!= *(double *)pn->pValue )
		    sprintf(szNumber,"%.17g",((double *)pn->pValue)[i]);
		fprintf(fp,"%s",szNumber);
		break;
	    case 3:
		fprintf(fp,"\"%s\"",(char *)pn->pValue);
		break;
	    case 4:
		assert(pn->iSize == sizeof(uint64_t));
		fprintf(fp,"%"PRIu64,((uint64_t *)pn->pValue)[i]);
		break;
		}
	    }
	fprintf(fp,"%s\n",pn->pCount?"]":"");
	}
    fclose(fp);
    }

void msrCheckpoint(MSR msr,int iStep,double dTime) {
    struct inWrite in;
    double sec,dsec;
    if ( msr->param.achCheckpointPath[0] )
	_BuildName(msr,in.achOutFile,iStep, msr->param.achCheckpointPath);
    else
	_BuildName(msr,in.achOutFile,iStep, msr->param.achOutPath);
    in.nProcessors = msr->param.bParaWrite==0?1:(msr->param.nParaWrite<=1 ? msr->nThreads:msr->param.nParaWrite);
    if (msr->param.csm->val.bComove) {
	double dExp = csmTime2Exp(msr->param.csm,dTime);
	msrprintf(msr,"Writing checkpoint for Step: %d Time:%g Redshift:%g\n",
	    iStep,dTime,(1.0/dExp - 1.0));
	}
    else
	msrprintf(msr,"Writing checkpoint for Step: %d Time:%g\n",iStep,dTime);
    sec = msrTime();

    writeParameters(msr,in.achOutFile,iStep,dTime);

    pstCheckpoint(msr->pst,&in,sizeof(in),NULL,NULL);

    /* This is not necessary, but it means the bounds will be identical upon restore */
    BND bnd;
    pstCalcBound(msr->pst,NULL,0,&bnd,NULL);

    dsec = msrTime() - sec;
    msrprintf(msr,"Checkpoint has been successfully written, Wallclock: %f secs.\n", dsec);
    }

/*
** This routine validates the given parameters and makes any adjustments.
*/
static int validateParameters(MDL mdl,PRM prm,struct parameters *param) {
    
    if (prmSpecified(prm, "dMetalDiffsionCoeff") || prmSpecified(prm,"dThermalDiffusionCoeff")) {
	if (!prmSpecified(prm, "iDiffusion")) param->iDiffusion=1;
	}

    {
	int nCoolingSet=0;
	if (param->bGasIsothermal) nCoolingSet++; 
	if (param->bGasCooling) nCoolingSet++; 
	if (!prmSpecified(prm, "bGasAdiabatic") && nCoolingSet) param->bGasAdiabatic=0;
	else if (param->bGasAdiabatic) nCoolingSet++;

	if (nCoolingSet != 1) {
	    fprintf(stderr,"One of bGasAdiabatic (%d), bGasIsothermal (%d) and bGasCooling (%d) may be set\n", param->bGasAdiabatic, param->bGasIsothermal, param->bGasCooling);
	    assert(0);
	    }
	}



    /* Star parameter checks */

    if (param->bStarForm) {
	param->bAddDelete = 1;
	if (!prmSpecified(prm, "bFeedback")) param->bFeedback=1;
	}

    /* END Gas and Star Parameter Checks */

    if (param->nDigits < 1 || param->nDigits > 9) {
	(void) fprintf(stderr,"Unreasonable number of filename digits.\n");
	return 0;
	}

    /*
    ** Make sure that we have some setting for nReplicas if bPeriodic is set.
    */
    if (param->bPeriodic && !prmSpecified(prm,"nReplicas")) {
	param->nReplicas = 1;
	}
    /*
    ** Warn that we have a setting for nReplicas if bPeriodic NOT set.
    */
    if (!param->bPeriodic && param->nReplicas != 0) {
	printf("WARNING: nReplicas set to non-zero value for non-periodic!\n");
	}

    /*
    ** CUDA likes a larger group size
    */
    if (mdlCudaActive(mdl) && param->iCUDAQueueSize>0 && !prmSpecified(prm,"nGroup") && param->nGroup<256)
	param->nGroup = 256;


#ifndef USE_HDF5
    if (param->bHDF5) {
	printf("WARNING: HDF5 output was requested by is not supported: using Tipsy format\n");
	param->bHDF5 = 0;
	}
#endif

#ifdef MDL_FFTW
    if ( param->nGridPk ) {
	if (prmSpecified(prm,"nBinsPk")) {
	    if (param->nBinsPk > param->nGridPk/2) {
		param->nBinsPk = param->nGridPk/2;
		}
	    }
	else param->nBinsPk = param->nGridPk/2;
	if (param->nBinsPk > PST_MAX_K_BINS)
	    param->nBinsPk = PST_MAX_K_BINS;
	}
    if (param->iPkOrder<1 || param->iPkOrder>4) {
    	puts("ERROR: iPkOrder must be 1 (NGP), 2 (CIC), 3 (TSC) or 4 (PCS)");
    	return 0;
        }
    if ( param->nGrid ) {
	if (param->achInFile[0]) {
	    puts("ERROR: do not specify an input file when generating IC");
	    return 0;
	    }
	if ( param->iSeed == 0 ) {
	    //puts("ERROR: Random seed for IC not specified");
	    param->iSeed = time(NULL);
	    }
	if ( !prmSpecified(prm,"dBoxSize") || param->dBoxSize <= 0 ) {
	    puts("ERROR: Box size for IC not specified");
	    return 0;
	    }
	}
    /* Set the number of bins for the power spectrum measurement of linear species */
    if (param->nGridLin > 0){
        param->nBinsLinPk = param->nGridLin/2;
        if (param->nBinsLinPk > PST_MAX_K_BINS)
            param->nBinsLinPk = PST_MAX_K_BINS;
    }
#endif
    if (param->dTheta <= 0) {
	if (param->dTheta == 0 && param->bVWarnings)
	    fprintf(stderr,"WARNING: Zero opening angle may cause numerical problems\n");
	else if (param->dTheta < 0) {
	    fprintf(stderr,"ERROR: Opening angle must be non-negative\n");
	    return 0;
	    }
	}

    if (!prmSpecified(prm,"dFracNoDomainDecomp")) param->dFracDualTree = param->dFracNoDomainDecomp;
    if ( param->dFracDualTree > param->dFracNoDomainDecomp
	|| param->dFracNoDomainDecomp > param->dFracNoDomainRootFind
	|| param->dFracNoDomainRootFind > param->dFracNoDomainDimChoice
	|| param->dFracNoDomainDecomp<0.0 || param->dFracNoDomainDimChoice > 1.0 ) {
	puts("ERROR: check that 0 <= dFracNoDomainDecomp <= dFracNoDomainRootFind <= dFracNoDomainDimChoice <= 1");
	return 0;
	}

    if (param->bMemUnordered) {
	if (!prmSpecified(prm,"bNewKDK")) {
	    param->bNewKDK = 1;
	    }
	else if (!param->bNewKDK) {
	    puts("ERROR: bNewKDK must be zero when bMemUnordered is enabled!");
	    return 0;
	    }
	}

    /* Make sure that the old behaviour is obeyed. */
    if ( param->nSteps == 0 ) {
	if ( !prmSpecified(prm,"bDoAccOutput") ) param->bDoAccOutput = 1;
	if ( !prmSpecified(prm,"bDoPotOutput") ) param->bDoPotOutput = 1;
	}

    /*
     * Softening
     */
    if (param->bPhysicalSoft ) {
	if (param->bPhysicalSoft && !param->csm->val.bComove) {
	    printf("WARNING: bPhysicalSoft reset to 0 for non-comoving (bComove == 0)\n");
	    param->bPhysicalSoft = 0;
	    }
	}
    /*
    ** Determine the period of the box that we are using.
    ** Set the new d[xyz]Period parameters which are now used instead
    ** of a single dPeriod, but we still want to have compatibility
    ** with the old method of setting dPeriod.
    */
    if (prmSpecified(prm,"dPeriod") &&
	    !prmSpecified(prm,"dxPeriod")) {
	param->dxPeriod = param->dPeriod;
	}
    if (prmSpecified(prm,"dPeriod") &&
	    !prmSpecified(prm,"dyPeriod")) {
	param->dyPeriod = param->dPeriod;
	}
    if (prmSpecified(prm,"dPeriod") &&
	    !prmSpecified(prm,"dzPeriod")) {
	param->dzPeriod = param->dPeriod;
	}
    /*
    ** Periodic boundary conditions can be disabled along any of the
    ** x,y,z axes by specifying a period of zero for the given axis.
    ** Internally, the period is set to infinity (Cf. pkdBucketWalk()
    ** and pkdDrift(); also the INTERSECT() macro in smooth.h).
    */
    if (param->dPeriod  == 0) param->dPeriod  = FLOAT_MAXVAL;
    if (param->dxPeriod == 0) param->dxPeriod = FLOAT_MAXVAL;
    if (param->dyPeriod == 0) param->dyPeriod = FLOAT_MAXVAL;
    if (param->dzPeriod == 0) param->dzPeriod = FLOAT_MAXVAL;
    /*
    ** At the moment, integer positions are only really safe in periodic boxes!Wr
    */
#ifdef INTEGER_POSITION
    if (!param->bPeriodic||param->dxPeriod!=1.0||param->dyPeriod!=1.0||param->dzPeriod!=1.0) {
	fprintf(stderr,"WARNING: Integer coordinates are enabled but the the box is not periodic\n"
	               "       and/or the box size is not 1. Set bPeriodic=1 and dPeriod=1.\n");
	}
#endif

    if (!prmSpecified(prm,"dTheta20")) param->dTheta20 = param->dTheta;
    if (!prmSpecified(prm,"dTheta2")) param->dTheta2 = param->dTheta20;

    /*
    ** Check if fast gas boundaries are needed.
    */
    if (param->bDoGas) {
	param->bMemNodeSphBounds = 1;
    }
    /*
    ** Check timestepping and gravity combinations.
    */
    assert(param->iMaxRung <= IRUNGMAX);
    if (param->bDoGravity) {
	/* Potential is optional, but the default for gravity */
	if (!prmSpecified(prm,"bMemPotential")) param->bMemPotential = 1;
	if (param->iMaxRung < 1) {
	    param->iMaxRung = 0;
	    if (param->bVWarnings) fprintf(stderr,"WARNING: iMaxRung set to 0, SINGLE STEPPING run!\n");
	    /*
	    ** For single stepping we don't need fancy timestepping variables.
	    */
	    param->bMemNodeAcceleration = 0;
	    param->bMemNodeVelocity = 0;
	}
	else {
	    if (param->bEpsAccStep) {
		param->bAccelStep = 1;
	    }
	    if ((param->bAccelStep || param->bDensityStep) && param->bGravStep) {
		/*
		** We cannot combine these 2 types of timestepping criteria, we need to choose one
		** or the other basic timestep criterion, in this case we choose only bGravStep.
		*/
		param->bAccelStep = 0;
		param->bEpsAccStep = 0;
		param->bDensityStep = 0;
		if (param->bVWarnings) fprintf(stderr,"WARNING: bGravStep set in combination with older criteria, now using ONLY bGravStep!\n");
	    }
	    else if (!param->bAccelStep && !param->bGravStep && !param->bDensityStep) {
		param->bGravStep = 1;
		if (param->bVWarnings) fprintf(stderr,"WARNING: none of bAccelStep, bDensityStep, or bGravStep set, now using bGravStep!\n");
	    }
	    /*
	    ** Set the needed memory model based on the chosen timestepping method.
	    */
	    if (param->bGravStep) {
		param->bMemNodeAcceleration = 1;
		if (param->iTimeStepCrit == 1) {
		    param->bMemNodeVelocity = 1;
		}
	    } 
	    else {
		param->bMemNodeAcceleration = 0;
		param->bMemNodeVelocity = 0;
	    }
	}
    }



    return 1;
    }


int msrInitialize(MSR *pmsr,MDL mdl,int argc,char **argv) {
    MSR msr;
    int i,j,ret;
    struct inSetAdd inAdd;
    char ach[256];
    int bDoRestore;

    msr = (MSR)malloc(sizeof(struct msrContext));
    assert(msr != NULL);
    msr->mdl = mdl;
    msr->pst = NULL;
    msr->lcl.pkd = NULL;
    *pmsr = msr;
    csmInitialize(&msr->param.csm);
    /*
    ** Now setup for the input parameters.
    **
    ** NOTE: nThreads & bDiag are parsed here, but the actual values are
    ** read from the command line via mdlInitialize(). This means the
    ** values of nThreads & bDiag read by prmAddParam() are ignored!
    */
    prmInitialize(&msr->prm,_msrLeader,_msrTrailer);
    msr->param.nThreads = 1;
    prmAddParam(msr->prm,"nThreads",1,&msr->param.nThreads,sizeof(int),"sz",
		"<nThreads>");
    msr->param.bDiag = 0;
    prmAddParam(msr->prm,"bDiag",0,&msr->param.bDiag,sizeof(int),"d",
		"enable/disable per thread diagnostic output");
    msr->param.bDedicatedMPI = 0;
    prmAddParam(msr->prm,"bDedicatedMPI",0,&msr->param.bDedicatedMPI,sizeof(int),"dedicated",
		"enable/disable dedicated MPI thread");
    msr->param.bSharedMPI = 0;
    prmAddParam(msr->prm,"bSharedMPI",0,&msr->param.bSharedMPI,sizeof(int),"sharedmpi",
		"enable/disable extra dedicated MPI thread");
    msr->param.bNoGrav = 0;
    prmAddParam(msr->prm,"bNoGrav",0,&msr->param.bNoGrav,sizeof(int),"nograv",
		"enable/disable Gravity calulation for testing = -nograv");
    msr->param.bOverwrite = 0;
    prmAddParam(msr->prm,"bOverwrite",0,&msr->param.bOverwrite,sizeof(int),
		"overwrite","enable/disable overwrite safety lock = -overwrite");
    msr->param.bVWarnings = 1;
    prmAddParam(msr->prm,"bVWarnings",0,&msr->param.bVWarnings,sizeof(int),
		"vwarnings","enable/disable warnings = +vwarnings");
    msr->param.bVStart = 1;
    prmAddParam(msr->prm,"bVStart",0,&msr->param.bVStart,sizeof(int),
		"vstart","enable/disable verbose start = +vstart");
    msr->param.bVStep = 1;
    prmAddParam(msr->prm,"bVStep",0,&msr->param.bVStep,sizeof(int),
		"vstep","enable/disable verbose step = +vstep");
    msr->param.bVRungStat = 1;
    prmAddParam(msr->prm,"bVRungStat",0,&msr->param.bVRungStat,sizeof(int),
		"vrungstat","enable/disable rung statistics = +vrungstat");
    msr->param.bVDetails = 0;
    prmAddParam(msr->prm,"bVDetails",0,&msr->param.bVDetails,sizeof(int),
		"vdetails","enable/disable verbose details = +vdetails");
    msr->param.nDigits = 5;
    prmAddParam(msr->prm,"nDigits",1,&msr->param.nDigits,sizeof(int),"nd",
		"<number of digits to use in output filenames> = 5");
    msr->param.bPeriodic = 0;
    prmAddParam(msr->prm,"bPeriodic",0,&msr->param.bPeriodic,sizeof(int),"p",
		"periodic/non-periodic = -p");
    msr->param.bRestart = 0;
    prmAddParam(msr->prm,"bRestart",0,&msr->param.bRestart,sizeof(int),"restart",
		"restart from checkpoint");
    msr->param.bParaRead = 1;
    prmAddParam(msr->prm,"bParaRead",0,&msr->param.bParaRead,sizeof(int),"par",
		"enable/disable parallel reading of files = +par");
    msr->param.bParaWrite = 0;
    prmAddParam(msr->prm,"bParaWrite",0,&msr->param.bParaWrite,sizeof(int),"paw",
		"enable/disable parallel writing of files = +paw");
    msr->param.nParaRead = 0;
    prmAddParam(msr->prm,"nParaRead",1,&msr->param.nParaRead,sizeof(int),"npar",
		"number of threads to read with during parallel read = 0 (unlimited)");
    msr->param.nParaWrite = 0;
    prmAddParam(msr->prm,"nParaWrite",1,&msr->param.nParaWrite,sizeof(int),"npaw",
		"number of threads to write with during parallel write = 0 (unlimited)");
    msr->param.bDoDensity = 1;
    prmAddParam(msr->prm,"bDoDensity",0,&msr->param.bDoDensity,sizeof(int),
		"den","enable/disable density outputs = +den");
#ifdef USE_PNG
    msr->param.nPNGResolution = 0;
    prmAddParam(msr->prm,"nPNGResolution",1,&msr->param.nPNGResolution,sizeof(int),
		"png","PNG output resolution (zero disables) = 0");
#endif
    msr->param.nBucket = 16;
    prmAddParam(msr->prm,"nBucket",1,&msr->param.nBucket,sizeof(int),"b",
		"<max number of particles in a bucket> = 16");
    msr->param.nGroup = 64;
    prmAddParam(msr->prm,"nGroup",1,&msr->param.nGroup,sizeof(int),"grp",
		"<max number of particles in a group> = 64");
    msr->param.n2min = 50;
    prmAddParam(msr->prm,"n2min",1,&msr->param.n2min,sizeof(int),"nn",
		"<minimum number of p-p interactions for using c-c interactions> = 50");
    msr->param.iStartStep = 0;
    prmAddParam(msr->prm,"iStartStep",1,&msr->param.iStartStep,
		sizeof(int),"nstart","<initial step numbering> = 0");
    msr->param.nSteps = 0;
    prmAddParam(msr->prm,"nSteps",1,&msr->param.nSteps,sizeof(int),"n",
		"<number of timesteps> = 0");
    msr->param.iOutInterval = 0;
    prmAddParam(msr->prm,"iOutInterval",1,&msr->param.iOutInterval,sizeof(int),
		"oi","<number of timesteps between snapshots> = 0");
    msr->param.iFofInterval = 0;
    prmAddParam(msr->prm,"iFofInterval",1,&msr->param.iFofInterval,sizeof(int),
		"fof","<number of timesteps between fof group finding> = 0");
    strcpy(msr->param.achOutTypes,"rvmsp");
    prmAddParam(msr->prm,"achOutTypes",3,msr->param.achOutTypes,256,"ot",
		"<output types for snapshot> = \"rvmsp\"");
    msr->param.iCheckInterval = 0;
    prmAddParam(msr->prm,"iCheckInterval",1,&msr->param.iCheckInterval,sizeof(int),
		"ci","<number of timesteps between checkpoints> = 0");
    strcpy(msr->param.achCheckTypes,"RVMSP");
    prmAddParam(msr->prm,"achCheckTypes",3,msr->param.achCheckTypes,256,"ct",
		"<output types for checkpoints> = \"RVMSP\"");
    msr->param.iLogInterval = 1;
    prmAddParam(msr->prm,"iLogInterval",1,&msr->param.iLogInterval,sizeof(int),
		"ol","<number of timesteps between logfile outputs> = 1");
    msr->param.iPkInterval = 1;
    prmAddParam(msr->prm,"iPkInterval",1,&msr->param.iPkInterval,sizeof(int),
		"opk","<number of timesteps between pk outputs> = 1");
    msr->param.iDeltakInterval = 0;
    prmAddParam(msr->prm,"iDeltakInterval",1,&msr->param.iDeltakInterval,sizeof(int),
		"odk","<number of timesteps between DeltaK outputs> = 0 (off)");
    msr->param.bEwald = 1;
    prmAddParam(msr->prm,"bEwald",0,&msr->param.bEwald,sizeof(int),"ewald",
		"enable/disable Ewald correction = +ewald");
    msr->param.iEwOrder = 4;
    prmAddParam(msr->prm,"iEwOrder",1,&msr->param.iEwOrder,sizeof(int),"ewo",
		"<Ewald multipole expansion order: 1, 2, 3 or 4> = 4");
    msr->param.nReplicas = 0;
    prmAddParam(msr->prm,"nReplicas",1,&msr->param.nReplicas,sizeof(int),"nrep",
		"<nReplicas> = 0 for -p, or 1 for +p");
    msr->param.dSoft = 0.0;
    prmAddParam(msr->prm,"dSoft",2,&msr->param.dSoft,sizeof(double),"e",
		"<gravitational softening length> = 0.0");
    msr->param.dSoftMax = 0.0;
    prmAddParam(msr->prm,"dSoftMax",2,&msr->param.dSoftMax,sizeof(double),"eMax",
		"<maximum comoving gravitational softening length (abs or multiplier)> = 0.0");
    msr->param.bPhysicalSoft = 0;
    prmAddParam(msr->prm,"bPhysicalSoft",0,&msr->param.bPhysicalSoft,sizeof(int),"PhysSoft",
		"<Physical gravitational softening length> -PhysSoft");
    msr->param.bSoftMaxMul = 1;
    prmAddParam(msr->prm,"bSoftMaxMul",0,&msr->param.bSoftMaxMul,sizeof(int),"SMM",
		"<Use maximum comoving gravitational softening length as a multiplier> +SMM");
    msr->param.nSoftNbr = 32;
    prmAddParam(msr->prm,"nSoftNbr",1,&msr->param.nSoftNbr,sizeof(int),"VarSoft",
		"<Neighbours for Variable gravitational softening length> 32");
    msr->param.bSoftByType = 1;
    prmAddParam(msr->prm,"bSoftByType",0,&msr->param.bSoftByType,sizeof(int),"SBT",
		"<Variable gravitational softening length by Type> +SBT");
    msr->param.bDoSoftOutput = 0;
    prmAddParam(msr->prm,"bDoSoftOutput",0,&msr->param.bDoSoftOutput,sizeof(int),
		"softout","enable/disable soft outputs = -softout");
    msr->param.bDoAccOutput = 0;
    prmAddParam(msr->prm,"bDoAccOutput",0,&msr->param.bDoAccOutput,sizeof(int),
		"accout","enable/disable acceleration outputs = -accout");
    msr->param.bDoPotOutput = 0;
    prmAddParam(msr->prm,"bDoPotOutput",0,&msr->param.bDoPotOutput,sizeof(int),
		"potout","enable/disable potential outputs = -potout");
    msr->param.bDoRungOutput = 0;
    prmAddParam(msr->prm,"bDoRungOutput",0,&msr->param.bDoRungOutput,sizeof(int),
		"rungout","enable/disable rung outputs = -rungout");
    msr->param.bDoRungDestOutput = 0;
    prmAddParam(msr->prm,"bDoRungDestOutput",0,&msr->param.bDoRungDestOutput,sizeof(int),
		"rungdestout","enable/disable rung destination outputs = -rungdestout");
    msr->param.dDelta = 0.0;
    prmAddParam(msr->prm,"dDelta",2,&msr->param.dDelta,sizeof(double),"dt",
		"<time step>");
    msr->param.dEta = 0.2;
    prmAddParam(msr->prm,"dEta",2,&msr->param.dEta,sizeof(double),"eta",
		"<time step criterion> = 0.1");
    msr->param.bGravStep = 0;
    prmAddParam(msr->prm,"bGravStep",0,&msr->param.bGravStep,sizeof(int),
		"gs","<Gravity timestepping according to iTimeStep Criterion>");
    msr->param.bEpsAccStep = 0;
    prmAddParam(msr->prm,"bEpsAccStep",0,&msr->param.bEpsAccStep,sizeof(int),
		"ea", "<Sqrt(Epsilon on a) timestepping>");
    msr->param.bDensityStep = 0;
    prmAddParam(msr->prm,"bDensityStep",0,&msr->param.bDensityStep,sizeof(int),
		"isrho", "<Sqrt(1/Rho) timestepping>");
    msr->param.iTimeStepCrit = 0;
    prmAddParam(msr->prm,"iTimeStepCrit",1,&msr->param.iTimeStepCrit,sizeof(int),
		"tsc", "<Criteria for dynamical time-stepping>");
    msr->param.nPartRhoLoc = 32;
    prmAddParam(msr->prm,"nPartRhoLoc",1,&msr->param.nPartRhoLoc,sizeof(int),
		"nprholoc", "<Number of particles for local density in dynamical time-stepping>");
    msr->param.dPreFacRhoLoc = 4.0*M_PI/3.0;
    prmAddParam(msr->prm,"dPreFacRhoLoc",2,&msr->param.dPreFacRhoLoc,sizeof(double),
		"dprefacrholoc", "<Pre-factor for local density in dynamical time-stepping>");
    msr->param.dFacExcludePart = 100;
    prmAddParam(msr->prm,"dFacExcludePart",2,&msr->param.dFacExcludePart,sizeof(double),
		"dfacexclp", "<Pre-factor for exluding far away particles on ILP list>");
    msr->param.dEccFacMax = 3000;
    prmAddParam(msr->prm,"dEccFacMax",2,&msr->param.dEccFacMax,sizeof(double),
		"deccfacmax", "<Maximum correction factor for eccentricity correction>");
    msr->param.nPartColl = 0;
    prmAddParam(msr->prm,"nPartColl",1,&msr->param.nPartColl,sizeof(int),
		"npcoll", "<Number of particles in collisional regime>");
    msr->param.nTruncateRung = 0;
    prmAddParam(msr->prm,"nTruncateRung",1,&msr->param.nTruncateRung,sizeof(int),"nTR",
		"<number of MaxRung particles to delete MaxRung> = 0");
    msr->param.iMaxRung = IRUNGMAX;
    sprintf(ach,"<maximum timestep rung> = %d",msr->param.iMaxRung);
    prmAddParam(msr->prm,"iMaxRung",1,&msr->param.iMaxRung,sizeof(int),
		"mrung",ach);
    msr->param.nRungVeryActive = 31;
    prmAddParam(msr->prm,"nRungVeryActive",1,&msr->param.nRungVeryActive,
		sizeof(int), "nvactrung", "<timestep rung to use very active timestepping>");
    msr->param.nPartVeryActive = 0;
    prmAddParam(msr->prm,"nPartVeryActive",1,&msr->param.nPartVeryActive,
		sizeof(int), "nvactpart", "<number of particles to use very active timestepping>");
    msr->param.bNewKDK = 0;
    prmAddParam(msr->prm,"bNewKDK",0,&msr->param.bNewKDK,
		sizeof(int), "NewKDK", "<Use new implementation of KDK time stepping=no>");
    msr->param.dEwCut = 2.6;
    prmAddParam(msr->prm,"dEwCut",2,&msr->param.dEwCut,sizeof(double),"ew",
		"<dEwCut> = 2.6");
    msr->param.dEwhCut = 2.8;
    prmAddParam(msr->prm,"dEwhCut",2,&msr->param.dEwhCut,sizeof(double),"ewh",
		"<dEwhCut> = 2.8");
    msr->param.dTheta = 0.7;
    msr->param.dTheta2 = msr->param.dTheta;
    msr->param.dTheta20 = msr->param.dTheta;
    prmAddParam(msr->prm,"dTheta",2,&msr->param.dTheta,sizeof(double),"theta",
		"<Barnes opening criterion> = 0.8");
    prmAddParam(msr->prm,"dTheta20",2,&msr->param.dTheta20,sizeof(double),
		"theta20","<Barnes opening criterion for 2 < z <= 20> = 0.8");
    prmAddParam(msr->prm,"dTheta2",2,&msr->param.dTheta2,sizeof(double),
		"theta2","<Barnes opening criterion for z <= 2> = 0.8");
    msr->param.dPeriod = 1.0;
    prmAddParam(msr->prm,"dPeriod",2,&msr->param.dPeriod,sizeof(double),"L",
		"<periodic box length> = 1.0");
    msr->param.dxPeriod = 1.0;
    prmAddParam(msr->prm,"dxPeriod",2,&msr->param.dxPeriod,sizeof(double),"Lx",
		"<periodic box length in x-dimension> = 1.0");
    msr->param.dyPeriod = 1.0;
    prmAddParam(msr->prm,"dyPeriod",2,&msr->param.dyPeriod,sizeof(double),"Ly",
		"<periodic box length in y-dimension> = 1.0");
    msr->param.dzPeriod = 1.0;
    prmAddParam(msr->prm,"dzPeriod",2,&msr->param.dzPeriod,sizeof(double),"Lz",
		"<periodic box length in z-dimension> = 1.0");
    msr->param.achInFile[0] = 0;
    prmAddParam(msr->prm,"achInFile",3,msr->param.achInFile,256,"I",
		"<input file name> (file in TIPSY binary format)");
    strcpy(msr->param.achOutName,"pkdgrav");
    prmAddParam(msr->prm,"achOutName",3,msr->param.achOutName,256,"o",
		"<output name for snapshots and logfile> = \"pkdgrav\"");
    strcpy(msr->param.achOutPath,"");
    prmAddParam(msr->prm,"achOutPath",3,msr->param.achOutPath,256,"op",
		"<output path for snapshots and logfile> = \"\"");
    strcpy(msr->param.achIoPath,"");
    prmAddParam(msr->prm,"achIoPath",3,msr->param.achIoPath,256,"iop",
		"<output path for snapshots and logfile> = \"\"");
    strcpy(msr->param.achCheckpointPath,"");
    prmAddParam(msr->prm,"achCheckpointPath",3,msr->param.achCheckpointPath,256,"cpp",
		"<output path for checkpoints> = \"\"");
    msr->param.csm->val.bComove = 0;
    prmAddParam(msr->prm,"bComove",0,&msr->param.csm->val.bComove,sizeof(int),
		"cm", "enable/disable comoving coordinates = -cm");
    msr->param.csm->val.dHubble0 = 0.0;
    prmAddParam(msr->prm,"dHubble0",2,&msr->param.csm->val.dHubble0,
		sizeof(double),"Hub", "<dHubble0> = 0.0");
    msr->param.csm->val.dOmega0 = 1.0;
    prmAddParam(msr->prm,"dOmega0",2,&msr->param.csm->val.dOmega0,
		sizeof(double),"Om", "<dOmega0> = 1.0");
    msr->param.csm->val.dLambda = 0.0;
    prmAddParam(msr->prm,"dLambda",2,&msr->param.csm->val.dLambda,
		sizeof(double),"Lambda", "<dLambda> = 0.0");
    msr->param.csm->val.dOmegaDE = 0.0;
    prmAddParam(msr->prm,"dOmegaDE",2,&msr->param.csm->val.dOmegaDE,
		sizeof(double),"OmDE", "Omega for Dark Energy using w0 and wa parameters: <dOmegaDE> = 0.0");
    msr->param.csm->val.w0 = -1.0;
    prmAddParam(msr->prm,"w0",2,&msr->param.csm->val.w0,
		sizeof(double),"w0", "w0 parameter for Dark Energy <w0> = -1.0 (pure Lambda)");
    msr->param.csm->val.wa = 0.0;
    prmAddParam(msr->prm,"wa",2,&msr->param.csm->val.wa,
		sizeof(double),"wa", "wa parameter for Dark Energy <wa> = 0.0 (pure Lambda)");
    msr->param.csm->val.dOmegaRad = 0.0;
    prmAddParam(msr->prm,"dOmegaRad",2,&msr->param.csm->val.dOmegaRad,
		sizeof(double),"Omrad", "<dOmegaRad> = 0.0");
    msr->param.csm->val.dOmegab = 0.0;
    prmAddParam(msr->prm,"dOmegab",2,&msr->param.csm->val.dOmegab,
		sizeof(double),"Omb", "<dOmegab> = 0.0");
    msr->param.csm->val.dSigma8 = 0.0;
    prmAddParam(msr->prm,"dSigma8",2,&msr->param.csm->val.dSigma8,
		sizeof(double),"S8", "<dSimga8> = 0.0");
    msr->param.csm->val.dNormalization = 0.0;
    prmAddParam(msr->prm,"dNormalization",2,&msr->param.csm->val.dNormalization,
		sizeof(double),"As", "<dNormalization> = 0.0");
    msr->param.csm->val.dSpectral = 0.0;
    prmAddParam(msr->prm,"dSpectral",2,&msr->param.csm->val.dSpectral,
		sizeof(double),"ns", "<dSpectral> = 0.0");
    msr->param.csm->val.dRunning = 0.0;
    prmAddParam(msr->prm,"dRunning",2,&msr->param.csm->val.dRunning,
        sizeof(double), "alphas", "Primordial tilt running: <dRunning> = 0.0");
    msr->param.csm->val.dPivot = 0.05;
    prmAddParam(msr->prm,"dPivot",2,&msr->param.csm->val.dPivot,
        sizeof(double), "kpivot", "Primordial pivot scale in 1/Mpc (not h/Mpc): <dPivot> = 0.05");
    strcpy(msr->param.achDataSubPath,"");
    prmAddParam(msr->prm,"achDataSubPath",3,msr->param.achDataSubPath,256,
		NULL,NULL);
    msr->param.dExtraStore = 0.1;
    prmAddParam(msr->prm,"dExtraStore",2,&msr->param.dExtraStore,
		sizeof(double),NULL,NULL);
    msr->param.bDualTree = 0;
    prmAddParam(msr->prm,"bDualTree",0,&msr->param.bDualTree,sizeof(int),"2tree",
		"enable/disable second tree for active rungs = -2tree");
    msr->param.nTreeBitsLo = 14;
    prmAddParam(msr->prm,"nTreeBitsLo",1,&msr->param.nTreeBitsLo,
	sizeof(int),"treelo",
	"<number of low bits for tree> = 14");
    msr->param.nTreeBitsHi = 18;
    prmAddParam(msr->prm,"nTreeBitsHi",1,&msr->param.nTreeBitsHi,
	sizeof(int),"treehi",
	"<number of high bits for tree> = 18");
#ifdef MDL_CACHE_SIZE
    msr->param.iCacheSize = MDL_CACHE_SIZE;
#else
    msr->param.iCacheSize = 0;
#endif
    prmAddParam(msr->prm,"iCacheSize",1,&msr->param.iCacheSize,sizeof(int),"cs",
		"<size of the MDL cache (0=default)> = 0");
    msr->param.iWorkQueueSize = 0;
    prmAddParam(msr->prm,"iWorkQueueSize",1,&msr->param.iWorkQueueSize,sizeof(int),"wqs",
		"<size of the MDL work queue> = 0");
    msr->param.iCUDAQueueSize = 8;
    prmAddParam(msr->prm,"iCUDAQueueSize",1,&msr->param.iCUDAQueueSize,sizeof(int),"cqs",
		"<size of the CUDA work queue> = 8");
    msr->param.nSmooth = 64;
    prmAddParam(msr->prm,"nSmooth",1,&msr->param.nSmooth,sizeof(int),"s",
		"<number of particles to smooth over> = 64");
    msr->param.bStandard = 1;
    prmAddParam(msr->prm,"bStandard",0,&msr->param.bStandard,sizeof(int),"std",
		"output in standard TIPSY binary format = -std");
    msr->param.iCompress = 0;
    prmAddParam(msr->prm,"iCompress",1,&msr->param.iCompress,sizeof(int),NULL,
		"compression format, 0=none, 1=gzip, 2=bzip2");
    msr->param.bHDF5 = 0;
    prmAddParam(msr->prm,"bHDF5",0,&msr->param.bHDF5,sizeof(int),"hdf5",
		"output in HDF5 format = -hdf5");
    msr->param.bDoublePos = 0;
    prmAddParam(msr->prm,"bDoublePos",0,&msr->param.bDoublePos,sizeof(int),"dp",
		"input/output double precision positions (standard format only) = -dp");
    msr->param.bDoubleVel = 0;
    prmAddParam(msr->prm,"bDoubleVel",0,&msr->param.bDoubleVel,sizeof(int),"dv",
		"input/output double precision velocities (standard format only) = -dv");
    msr->param.bLightCone = 0;
    prmAddParam(msr->prm,"bLightCone",0,&msr->param.bLightCone,sizeof(int),"lc",
		"output light cone data = -lc");
    msr->param.nSideHealpix = 8192;
    prmAddParam(msr->prm,"nSideHealpix",1,&msr->param.nSideHealpix,
		sizeof(int),"healpix",
		"<Number per side of the healpix map> = 8192");
    msr->param.bLightConeParticles = 0;
    prmAddParam(msr->prm,"bLightConeParticles",0,&msr->param.bLightConeParticles,sizeof(int),"lcp",
		"output light cone particles = -lcp");
    msr->param.bInFileLC = 0;
    prmAddParam(msr->prm,"bInFileLC",0,&msr->param.bInFileLC,sizeof(int),"lcin",
		"input light cone data = -lcin");
    msr->param.dRedshiftLCP = 0;
    prmAddParam(msr->prm,"dRedshiftLCP",2,&msr->param.dRedshiftLCP,sizeof(double),"zlcp",
		"starting redshift to output light cone particles = 0");
    msr->param.bCenterOfMassExpand = 1;
    prmAddParam(msr->prm,"bCenterOfMassExpand",0,&msr->param.bCenterOfMassExpand,sizeof(int),"CoM",
		"use multipole expansions about the center of mass = +CoM");
    msr->param.dRedTo = 0.0;
    prmAddParam(msr->prm,"dRedTo",2,&msr->param.dRedTo,sizeof(double),"zto",
		"specifies final redshift for the simulation");
    msr->param.dRedFrom = 0.0;
    prmAddParam(msr->prm,"dRedFrom",2,&msr->param.dRedFrom,sizeof(double),"z",
		"specifies initial redshift for the simulation");
    msr->param.dGrowDeltaM = 0.0;
    prmAddParam(msr->prm,"dGrowDeltaM",2,&msr->param.dGrowDeltaM,
		sizeof(double),"gmdm","<Total growth in mass/particle> = 0.0");
    msr->param.dGrowStartT = 0.0;
    prmAddParam(msr->prm,"dGrowStartT",2,&msr->param.dGrowStartT,
		sizeof(double),"gmst","<Start time for growing mass> = 0.0");
    msr->param.dGrowEndT = 1.0;
    prmAddParam(msr->prm,"dGrowEndT",2,&msr->param.dGrowEndT,
		sizeof(double),"gmet","<End time for growing mass> = 1.0");
    msr->param.dFracDualTree = 0.05;
    prmAddParam(msr->prm,"dFracDualTree",2,&msr->param.dFracDualTree,
		sizeof(double),"fndt",
		"<Fraction of Active Particles for to use a dual tree> = dFracNoDomainDecomp ");
    msr->param.dFracNoDomainDecomp = 0.1;
    prmAddParam(msr->prm,"dFracNoDomainDecomp",2,&msr->param.dFracNoDomainDecomp,
		sizeof(double),"fndd",
		"<Fraction of Active Particles for no DD> = 0.1");
    msr->param.dFracNoDomainRootFind = 0.1;
    prmAddParam(msr->prm,"dFracNoDomainRootFind",2,&msr->param.dFracNoDomainRootFind,
		sizeof(double),"fndrf",
		"<Fraction of Active Particles for no DD root finding> = 0.1");
    msr->param.dFracNoDomainDimChoice = 0.1;
    prmAddParam(msr->prm,"dFracNoDomainDimChoice",2,&msr->param.dFracNoDomainDimChoice,
		sizeof(double),"fnddc",
		"<Fraction of Active Particles for no DD dimension choice> = 0.1");
    msr->param.bDoGravity = 1;
    prmAddParam(msr->prm,"bDoGravity",0,&msr->param.bDoGravity,sizeof(int),"g",
		"enable/disable interparticle gravity = +g");
    msr->param.bAarsethStep = 0;
    prmAddParam(msr->prm,"bAarsethStep",0,&msr->param.bAarsethStep,sizeof(int),
		"aas","<Aarseth timestepping>");
    msr->param.iWallRunTime = 0;
    prmAddParam(msr->prm,"iWallRunTime",1,&msr->param.iWallRunTime,
		sizeof(int),"wall",
		"<Maximum Wallclock time (in minutes) to run> = 0 = infinite");
    msr->param.iSignalSeconds = 0;
    prmAddParam(msr->prm,"iSignalSeconds",1,&msr->param.iSignalSeconds,
		sizeof(int),"signal",
		"<Time (in seconds) that USR1 is sent before termination> = 0 = immediate");
    msr->param.bFindGroups = 0;
    prmAddParam(msr->prm,"bFindGroups",0,&msr->param.bFindGroups,sizeof(int),
		"groupfinder","<enable/disable group finder> = -groupfinder");
    msr->param.bFindHopGroups = 0;
    prmAddParam(msr->prm,"bFindHopGroups",0,&msr->param.bFindHopGroups,sizeof(int),
		"hop","<enable/disable phase-space group finder> = -hop");
    msr->param.dHopTau = -4.0;
    prmAddParam(msr->prm,"dHopTau",2,&msr->param.dHopTau,sizeof(double),"hoptau",
		"<linking length for Gasshopper (negative for multiples of softening)> = -4.0");
    msr->param.nMinMembers = 10;
    prmAddParam(msr->prm,"nMinMembers",1,&msr->param.nMinMembers,sizeof(int),
		"nMinMembers","<minimum number of group members> = 10");
    msr->param.dTau = 0.164;
    prmAddParam(msr->prm,"dTau",2,&msr->param.dTau,sizeof(double),"dTau",
		"<linking length for FOF in units of mean particle separation> = 0.164");
    msr->param.dEnvironment0 = -1.0;
    prmAddParam(msr->prm,"dEnvironment0",2,&msr->param.dEnvironment0,sizeof(double),"dEnv0",
		"<first radius for density environment about a group> = -1.0 (disabled)");
    msr->param.dEnvironment1 = -1.0;
    prmAddParam(msr->prm,"dEnvironment1",2,&msr->param.dEnvironment1,sizeof(double),"dEnv1",
		"<second radius for density environment about a group> = -1.0 (disabled)");
    msr->param.nBins = 0;
    prmAddParam(msr->prm,"nBins",1,&msr->param.nBins,sizeof(int),"nBins",
		"<number of bin in profiles, no profiles if 0 or negative> = 0");
    msr->param.iCenterType = 2;
    prmAddParam(msr->prm,"iCenterType",1,&msr->param.iCenterType,sizeof(int),"iCenterType",
		"<sets center type for group finder: 0 com; 1 potmin; 2 denmax> = 2");
    msr->param.binFactor = 0.2;
    prmAddParam(msr->prm,"binFactor",2,&msr->param.binFactor,sizeof(double),"binFactor",
		"<ratio of largest spherical bin to fof determined group radius> = 0.2");
    msr->param.fMinRadius = 1.0e-5;
    prmAddParam(msr->prm,"fMinRadius",2,&msr->param.fMinRadius,sizeof(double),"fMinRadius",
                "<radius of first, smallest spherical bin in the group profiles> = 1.0e-5");
    msr->param.bLogBins = 1;
    prmAddParam(msr->prm,"bLogBins",0,&msr->param.bLogBins,
		sizeof(int),"bLogBins","use logaritmic bins instead of linear = +bLogBins");
    msr->param.bTraceRelaxation = 0;
    prmAddParam(msr->prm,"bTraceRelaxation",0,&msr->param.bTraceRelaxation,sizeof(int),
		"rtrace","<enable/disable relaxation tracing> = -rtrace");

#ifdef MDL_FFTW
    msr->param.nBinsPk = 0;
    prmAddParam(msr->prm,"nBinsPk",1,&msr->param.nBinsPk,
		sizeof(int),"npk","<Number of log bins for P(k)> = nGridPk/2");
    msr->param.nGridPk = 0;
    prmAddParam(msr->prm,"nGridPk",1,&msr->param.nGridPk,
		sizeof(int),"pk","<Grid size for measure P(k) 0=disabled> = 0");
    msr->param.bPkInterlace = 1;
    prmAddParam(msr->prm,"bPkInterlace",0,&msr->param.bPkInterlace,
		sizeof(int),"pkinterlace","<Use interlacing to measure P(k)> = +pkinterlace");
    msr->param.iPkOrder = 4;
    prmAddParam(msr->prm,"iPkOrder",1,&msr->param.iPkOrder,
		sizeof(int),"pko","<Mass assignment order for measuring P(k) = 3");
    msr->param.bFixedAmpIC = 0;
    prmAddParam(msr->prm,"bFixedAmpIC",0,&msr->param.bFixedAmpIC,
		sizeof(int),"fixedamp","<Use fixed amplitude of 1 for ICs> = -fixedamp");
    msr->param.dFixedAmpPhasePI = 0.0;
    prmAddParam(msr->prm,"dFixedAmpPhasePI",2,&msr->param.dFixedAmpPhasePI,
		sizeof(double),"fixedphase","<Phase shift for fixed amplitude in units of PI> = 0.0");
    msr->param.nGridLin = 0;
    prmAddParam(msr->prm, "nGridLin", 1, &msr->param.nGridLin,
        sizeof(int), "lingrid", "<Grid size for linear species 0=disabled> =0");
    msr->param.bDoLinPkOutput = 0;
    prmAddParam(msr->prm, "bDoLinPkOutput", 0, &msr->param.bDoLinPkOutput,
        sizeof(int), "linPk", "<enable/disable power spectrum output for linear species> = 0");
#endif

    msr->param.iInflateStep = 0;
    prmAddParam(msr->prm,"iInflateStep",1,&msr->param.iInflateStep,
		sizeof(int),"iis","<Step when to inflate the number of particles> = 0");
    msr->param.nInflateReps = 0;
    prmAddParam(msr->prm,"nInflateReps",1,&msr->param.nInflateReps,
		sizeof(int),"nir","<Number of replicas when inflating> = 0");

    /* IC Generation */
    msr->param.csm->val.classData.bClass = 0;
    prmAddParam(msr->prm,"bClass",0,&msr->param.csm->val.classData.bClass,
		sizeof(int),"class","<Enable/disable the use of CLASS> = -class");
    msr->param.csm->val.classData.achFilename[0] = 0;
    prmAddParam(msr->prm, "achClassFilename", 3, msr->param.csm->val.classData.achFilename,
		256, "class_filename", "<Name of hdf5 file containing the CLASS data> -class_filename");
    msr->param.csm->val.classData.achLinSpecies[0] = 0;
    prmAddParam(msr->prm, "achLinSpecies", 3, msr->param.csm->val.classData.achLinSpecies,
                128, "lin_species",
                "<plus-separated string of linear species, e.g. \"ncdm[0]+g+metric\"> -lin_species");
    msr->param.h = 0.0;
    prmAddParam(msr->prm,"h",2,&msr->param.h,
		sizeof(double),"h","<hubble parameter h> = 0");
    msr->param.dBoxSize = 1.0;
    prmAddParam(msr->prm,"dBoxSize",2,&msr->param.dBoxSize,
		sizeof(double),"mpc","<Simulation Box size in Mpc> = 1.0");
    msr->param.nGrid = 0;
    prmAddParam(msr->prm,"nGrid",1,&msr->param.nGrid,
		sizeof(int),"grid","<Grid size for IC 0=disabled> = 0");
    msr->param.achTfFile[0] = 0;
    prmAddParam(msr->prm,"achTfFile",3,msr->param.achTfFile,256,"tf",
		"<transfer file name> (file in CMBFAST format)");
    msr->param.iSeed = 0;
    prmAddParam(msr->prm,"iSeed",1,&msr->param.iSeed,
		sizeof(int),"seed","<Random seed for IC> = 0");
    msr->param.b2LPT = 1;
    prmAddParam(msr->prm,"b2LPT",0,&msr->param.b2LPT,
		sizeof(int),"2lpt","<Enable/disable 2LPT> = 1");
#ifdef USE_PYTHON
    strcpy(msr->param.achScriptFile,"");
    prmAddParam(msr->prm,"achScript",3,msr->param.achScriptFile,256,"script",
		"<Python script for analysis> = \"\"");
#endif
    msr->param.bWriteIC = 0;
    prmAddParam(msr->prm,"bWriteIC",0,&msr->param.bWriteIC,
		sizeof(int),"wic","<Write IC after generating> = 0");

    /* Memory models */
    msr->param.bMemUnordered = 0;
    prmAddParam(msr->prm,"bMemUnordered",0,&msr->param.bMemUnordered,
		sizeof(int),"unordered","<Particles have no specific order> = -unordered");
    msr->param.bMemParticleID = 0;
    prmAddParam(msr->prm,"bMemParticleID",0,&msr->param.bMemParticleID,
		sizeof(int),"pid","<Particles have a unique identifier> = -pid");
    msr->param.bMemAcceleration = 0;
    prmAddParam(msr->prm,"bMemAcceleration",0,&msr->param.bMemAcceleration,
		sizeof(int),"Ma","<Particles have acceleration> = -Ma");
    msr->param.bMemVelocity = 0;
    prmAddParam(msr->prm,"bMemVelocity",0,&msr->param.bMemVelocity,
		sizeof(int),"Mv","<Particles have velocity> = -Mv");
    msr->param.bMemPotential = 0;
    prmAddParam(msr->prm,"bMemPotential",0,&msr->param.bMemPotential,
		sizeof(int),"Mp","<Particles have potential> = -Mp");
    msr->param.bMemGroups = 0;
    prmAddParam(msr->prm,"bMemGroups",0,&msr->param.bMemGroups,
		sizeof(int),"Mg","<Particles support group finding> = -Mg");
    msr->param.bMemMass = 0;
    prmAddParam(msr->prm,"bMemMass",0,&msr->param.bMemMass,
		sizeof(int),"Mm","<Particles have individual masses> = -Mm");
    msr->param.bMemSoft = 0;
    prmAddParam(msr->prm,"bMemSoft",0,&msr->param.bMemSoft,
		sizeof(int),"Ms","<Particles have individual softening> = -Ms");
    msr->param.bMemRelaxation = 0;
    prmAddParam(msr->prm,"bMemRelaxation",0,&msr->param.bMemRelaxation,
		sizeof(int),"Mr","<Particles have relaxation> = -Mr");
    msr->param.bMemVelSmooth = 0;
    prmAddParam(msr->prm,"bMemVelSmooth",0,&msr->param.bMemVelSmooth,
		sizeof(int),"Mvs","<Particles support velocity smoothing> = -Mvs");
    msr->param.bMemNodeMoment = 0;
    prmAddParam(msr->prm,"bMemNodeMoment",0,&msr->param.bMemNodeMoment,
		sizeof(int),"MNm","<Tree nodes support multipole moments> = 0");
    msr->param.bMemNodeAcceleration = 0;
    prmAddParam(msr->prm,"bMemNodeAcceleration",0,&msr->param.bMemNodeAcceleration,
		sizeof(int),"MNa","<Tree nodes support acceleration (for bGravStep)> = 0");
    msr->param.bMemNodeVelocity = 0;
    prmAddParam(msr->prm,"bMemNodeVelocity",0,&msr->param.bMemNodeVelocity,
		sizeof(int),"MNv","<Tree nodes support velocity (for iTimeStepCrit = 1)> = 0");
    msr->param.bMemNodeSphBounds = 0;
    prmAddParam(msr->prm,"bMemNodeSphBounds",0,&msr->param.bMemNodeSphBounds,
		sizeof(int),"MNsph","<Tree nodes support fast-gas bounds> = 0");

    msr->param.bMemNodeBnd = 1;
    /*prmAddParam(msr->prm,"bMemNodeBnd",1,&msr->param.bMemNodeBnd,
      sizeof(int),"MNbnd","<Tree nodes support 3D bounds> = 1");*/

    msr->param.bMemNodeVBnd = 0;
    prmAddParam(msr->prm,"bMemNodeVBnd",0,&msr->param.bMemNodeVBnd,
		sizeof(int),"MNvbnd","<Tree nodes support velocity bounds> = 0");

/* Gas Parameters */
    msr->param.bDoGas = 0;
    prmAddParam(msr->prm,"bDoGas",0,&msr->param.bDoGas,sizeof(int),"gas",
		"calculate gas/don't calculate gas = +gas");
    msr->param.bGasAdiabatic = 1;
    prmAddParam(msr->prm,"bGasAdiabatic",0,&msr->param.bGasAdiabatic,
		sizeof(int),"GasAdiabatic",
		"<Gas is Adiabatic> = +GasAdiabatic");
    msr->param.bGasIsothermal = 0;
    prmAddParam(msr->prm,"bGasIsothermal",0,&msr->param.bGasIsothermal,
		sizeof(int),"GasIsothermal",
		"<Gas is Isothermal> = +GasIsothermal");
    msr->param.bGasCooling = 0;
    prmAddParam(msr->prm,"bGasCooling",0,&msr->param.bGasCooling,
		sizeof(int),"GasCooling",
		"<Gas is Cooling> = +GasCooling");
    msr->param.bInitTFromCooling = 0;
    prmAddParam(msr->prm,"bInitTFromCooling",0,&msr->param.bInitTFromCooling,
		sizeof(int),"bInitTFromCooling",
		"set T (also E, Y, etc..) using Cooling initialization value = +bInitTFromCooling");
    msr->param.iRungCoolTableUpdate = 0;
    prmAddParam(msr->prm,"iRungCoolTableUpdate",1,&msr->param.iRungCoolTableUpdate,
		sizeof(int),"iRungCoolTableUpdate",
		"<Rung on which to update cool tables, def. 0>");
    msr->param.bInitTFromCooling = 0;
    msr->param.dEtaCourant = 0.4;
    prmAddParam(msr->prm,"dEtaCourant",2,&msr->param.dEtaCourant,sizeof(double),"etaC",
				"<Courant criterion> = 0.4");
    msr->param.dEtaUDot = 0.25;
    prmAddParam(msr->prm,"dEtauDot",2,&msr->param.dEtaUDot,sizeof(double),"etau",
		"<uDot timestep criterion> = 0.25");
    msr->param.dConstAlpha = 1.0; 
    prmAddParam(msr->prm,"dConstAlpha",2,&msr->param.dConstAlpha,
		sizeof(double),"alpha",
		"<Alpha constant in viscosity> = 1.0 or 0.5 (bBulkViscosity)");
    msr->param.dConstBeta = 2.0; 
    prmAddParam(msr->prm,"dConstBeta",2,&msr->param.dConstBeta,
		sizeof(double),"beta",
		"<Beta constant in viscosity> = 2.0 or 0.5 (bBulkViscosity)");
    msr->param.dConstGamma = 5.0/3.0;
    prmAddParam(msr->prm,"dConstGamma",2,&msr->param.dConstGamma,
		sizeof(double),"gamma",
		"<Ratio of specific heats> = 5/3");
    msr->param.dMeanMolWeight = 1.0;
    prmAddParam(msr->prm,"dMeanMolWeight",2,&msr->param.dMeanMolWeight,
		sizeof(double),"mmw",
		"<Mean molecular weight in amu> = 1.0");
    msr->param.dGasConst = 1.0;
    prmAddParam(msr->prm,"dGasConst",2,&msr->param.dGasConst,
		sizeof(double),"gcnst",
		"<Gas Constant>");
    msr->param.dKBoltzUnit = 1.0;
    prmAddParam(msr->prm,"dKBoltzUnit",2,&msr->param.dKBoltzUnit,
		sizeof(double),"kb",
		"<Boltzmann Constant in System Units>");
    msr->param.dhMinOverSoft = 0.0;
    prmAddParam(msr->prm,"dhMinOverSoft",2,&msr->param.dhMinOverSoft,
		sizeof(double),"hmin",
		"<Minimum h as a fraction of Softening> = 0.0");
    msr->param.dMetalDiffusionCoeff = 0;
    prmAddParam(msr->prm,"dMetalDiffusionCoeff",2,&msr->param.dMetalDiffusionCoeff,
		sizeof(double),"metaldiff",
		"<Coefficient in Metal Diffusion> = 0.0");
    msr->param.dThermalDiffusionCoeff = 0;
    prmAddParam(msr->prm,"dThermalDiffusionCoeff",2,&msr->param.dThermalDiffusionCoeff,
		sizeof(double),"thermaldiff",
		"<Coefficient in Thermal Diffusion> = 0.0");
    msr->param.dMsolUnit = 1.0;
    prmAddParam(msr->prm,"dMsolUnit",2,&msr->param.dMsolUnit,
		sizeof(double),"msu",
		"<Solar mass/system mass unit>");
    msr->param.dKpcUnit = 1000.0;
    prmAddParam(msr->prm,"dKpcUnit",2,&msr->param.dKpcUnit,
		sizeof(double),"kpcu",
		"<Kiloparsec/system length unit>");
    msr->param.ddHonHLimit = 0.1;
    prmAddParam(msr->prm,"ddHonHLimit",2,&msr->param.ddHonHLimit,
		sizeof(double),"dhonh",
		"<|dH|/H Limiter> = 0.1");
    msr->param.iViscosityLimiter = 1;
    prmAddParam(msr->prm,"iViscosityLimiter",1,&msr->param.iViscosityLimiter,sizeof(int),
		"vlim","<iViscosity Limiter> = 1");
    msr->param.iDiffusion = 0;
    prmAddParam(msr->prm,"iDiffusion",1,&msr->param.iDiffusion,sizeof(int),
		"idiff","<iDiffusion> = 0");
    msr->param.bAddDelete = 0;
    prmAddParam(msr->prm,"bAddDelete",0,&msr->param.bAddDelete,sizeof(int),
		"adddel","<Add Delete Particles> = 0");
    /* Star Form parameters */
    msr->param.bStarForm = 0;
    prmAddParam(msr->prm,"bStarForm",0,&msr->param.bStarForm,sizeof(int),
		"stfm","<Star Forming> = 0");
    msr->param.bFeedback = 0;
    prmAddParam(msr->prm,"bFeedback",0,&msr->param.bFeedback,sizeof(int),
		"fdbk","<Stars provide feedback> = 0");
    msr->param.SFdComovingDenMin = 0; /* RAMSES DEFAULT */
    prmAddParam(msr->prm,"SFdComovingDenMin", 2, &msr->param.SFdComovingDenMin,
		sizeof(double), "stODmin",
		"<Minimum overdensity for forming stars> = 2");
    msr->param.SFdPhysDenMin =  0.1/.76*1.66e-24; /* 0.1 nH/cc; RAMSES DEFAULT */
    prmAddParam(msr->prm,"SFdPhysDenMin", 2, &msr->param.SFdPhysDenMin,
		sizeof(double), "stPDmin",
		"<Minimum physical density for forming stars (gm/cc)> =  7e-26");
/* Duplicatee below.
    msr->param.SFdInitStarMass = 0;
    prmAddParam(msr->prm,"SFdInitStarMass", 2, &msr->param.SFdInitStarMass,
		sizeof(double), "stm0",
		"<Initial star mass> = 0");
*/
    msr->param.SFdESNPerStarMass = 5.0276521e+15; /* RAMSES DEFAULT */
    prmAddParam(msr->prm,"SFdESNPerStarMass", 2, &msr->param.SFdESNPerStarMass,
		sizeof(double), "ESNPerStarMass",
		"<ESN per star mass, erg per g of stars> = 1.25e16");
    msr->param.SFdTMax = 1e20; /* RAMSES DEFAULT */
    prmAddParam(msr->prm,"SFdTMax", 2, &msr->param.SFdTMax,
		sizeof(double), "SFdTMax",
		"<Maximum temperature for forming stars, K> = 3e4");
    msr->param.SFdEfficiency = 0.01; /* RAMSES DEFAULT */
    prmAddParam(msr->prm,"SFdEfficiency", 2, &msr->param.SFdEfficiency,
		sizeof(double), "SFdEfficiency",
		"<SF Efficiency> = 0.1");
    msr->param.SFdtCoolingShutoff = 30e6; /* RAMSES DEFAULT */
    prmAddParam(msr->prm,"SFdtCoolingShutoff", 2, &msr->param.SFdtCoolingShutoff,
		sizeof(double), "SFdtCoolingShutoff",
		"<SF Cooling Shutoff duration> = 30e6");

    msr->param.SFdtFeedbackDelay = 10e6; /* RAMSES DEFAULT */
    prmAddParam(msr->prm,"SFdtFeedbackDelay", 2, &msr->param.SFdtFeedbackDelay,
		sizeof(double), "SFdtFBD",
		"<SF FB delay> = 10e6");
    msr->param.SFdMassLossPerStarMass = 0.1; /* RAMSES DEFAULT */
    prmAddParam(msr->prm,"SFdMassLossPerStarMass", 2, &msr->param.SFdMassLossPerStarMass,
		sizeof(double), "SFMLPSM",
		"<SFMSPSM > = 0.1");
    msr->param.SFdZMassPerStarMass = 0.01; /* RAMSES DEFAULT */
    prmAddParam(msr->prm,"SFdZMassPerStarMass", 2, &msr->param.SFdZMassPerStarMass,
		sizeof(double), "SFZMPSM",
		"<SF ZMPSM> = 0.01");
    msr->param.SFdInitStarMass = 0; 
    prmAddParam(msr->prm,"SFdInitStarMass", 2, &msr->param.SFdInitStarMass,
		sizeof(double), "SFISM",
		"<SF ISM> = ?");
    msr->param.SFdMinGasMass = 0;
    prmAddParam(msr->prm,"SFdMinGasMass", 2, &msr->param.SFdMinGasMass,
		sizeof(double), "SFMGM",
		"<SF MGM> = ?");
    msr->param.SFdvFB = 100;
    prmAddParam(msr->prm,"SFdvFB", 2, &msr->param.SFdvFB,
		sizeof(double), "SFVFB",
		"<SF dvFB sound speed in FB region expected, km/s> = 100");

    msr->param.SFbdivv = 0;
    prmAddParam(msr->prm,"SFbdivv", 0, &msr->param.SFbdivv,
		sizeof(int), "SFbdivv",
		"<SF Use div v for star formation> = 1");
    /* END Gas/Star Parameters */
    msr->param.nOutputParticles = 0;
    prmAddArray(msr->prm,"lstOrbits",4,&msr->param.iOutputParticles,sizeof(uint64_t),&msr->param.nOutputParticles);

    /* New params added by IA for the hydrodynamics */
    msr->param.dCFLacc = 0.01;
    prmAddParam(msr->prm,"dCFLacc",2,&msr->param.dCFLacc,sizeof(double),"CFLacc",
				"<CFL for the acceleration criteria> = 0.01");
    msr->param.bMeshlessHydro = 0;
    prmAddParam(msr->prm,"bMeshlessHydro", 0, &msr->param.bMeshlessHydro,
		sizeof(int), "meshlessHydro",
		"Use the new implementation of the hydrodynamics");

    msr->param.bConservativeReSmooth = 1;
    prmAddParam(msr->prm,"bConservativeReSmooth", 0, &msr->param.bConservativeReSmooth,
		sizeof(int), "resmooth",
		"Use re-smoothing for the fluxes and gradients computation");

    msr->param.bGlobalDt = 0;
    prmAddParam(msr->prm,"bGlobalDt", 0, &msr->param.bGlobalDt,
		sizeof(int), "globaldt",
		"Force all particles to the same rung");

    msr->param.dFixedDelta = 0.0;
    prmAddParam(msr->prm,"dFixedDelta",2,&msr->param.dFixedDelta,sizeof(double),"FixedDelta",
				"<Fixed time step of the simulation (requires bGlobalDt=1)> = 0.0");

    msr->param.bIterativeSmoothingLength = 1;
    prmAddParam(msr->prm,"bIterativeSmoothingLength", 0, &msr->param.bIterativeSmoothingLength,
		sizeof(int), "iterh",
		"Use an iterative scheme to obtain h");

    msr->param.dNeighborsStd0 = 1;
    prmAddParam(msr->prm,"dNeighborsStd", 2, &msr->param.dNeighborsStd0,
		sizeof(double), "neighstd",
		"Maximum deviation from desired number of neighbors");

    msr->param.dMinDt = 0.;
    prmAddParam(msr->prm,"dMinDt", 2, &msr->param.dMinDt,
		sizeof(double), "minDt",
		"Minimum allowed timestep for the particles (in code units)");

#ifdef COOLING
    prmAddParam(msr->prm,"strCoolingTables",3,msr->param.strCoolingTables,256,"coolingtables",
		"Path to cooling tables");

    /// Hydrogen reionization
    msr->param.fH_reion_z = 11.5;
    prmAddParam(msr->prm,"fH_reion_z", 2, &msr->param.fH_reion_z,
		sizeof(double), "H_reion_z",
		"Redshift of Hydrogen reionization");

    msr->param.fH_reion_eV_p_H = 2.0;
    prmAddParam(msr->prm,"fH_reion_eV_p_H", 2, &msr->param.fH_reion_eV_p_H,
		sizeof(double), "H_reion_eV_p_H",
		"Energy (in eV) injected per proton during H reionization");

    /// Helium reionization
    msr->param.fHe_reion_eV_p_H = 2.0;
    prmAddParam(msr->prm,"fHe_reion_eV_p_H", 2, &msr->param.fHe_reion_eV_p_H,
		sizeof(double), "He_reion_eV_p_H",
		"Energy (in eV) injected per proton during He reionization");

    msr->param.fHe_reion_z_centre = 3.5;
    prmAddParam(msr->prm,"fHe_reion_z_centre", 2, &msr->param.fHe_reion_z_centre,
		sizeof(double), "He_reion_z_centre",
		"Mean redshift of Helium reionization");

    msr->param.fHe_reion_z_sigma = 0.5;
    prmAddParam(msr->prm,"fHe_reion_z_sigma", 2, &msr->param.fHe_reion_z_sigma,
		sizeof(double), "He_reion_z_sigma",
		"Redshift interval for Helium reionization");

    /// Relatives abundances of Ca and S
    msr->param.fCa_over_Si_in_Solar = 1.;
    msr->param.fS_over_Si_in_Solar = 1.;

    /// Temperature of the CMB at z=0
    msr->param.fT_CMB_0 = 2.725;
    prmAddParam(msr->prm,"fT_CMB_0", 2, &msr->param.fT_CMB_0,
		sizeof(float), "fT_CMB_0",
		"Temperature of the CMB at z=0");
#endif
#ifdef STAR_FORMATION
    msr->param.SFdMinOverDensity = 0.;
    prmAddParam(msr->prm,"SFdMinOverDensity", 2, &msr->param.SFdMinOverDensity,
		sizeof(double), "SFdMinOverDensity",
		"Minimium overdensity for allowing star formation");

    msr->param.SFdGasFraction = 0.;
    prmAddParam(msr->prm,"SFdGasFraction", 2, &msr->param.SFdGasFraction,
		sizeof(double), "SFdGasFraction",
		"Gas fraction (assumed constant) for the star formation");
#endif
    /* END of new params */

    msr->param.bAccelStep = 0;

#ifndef USE_DIAPOLE
    msr->param.bCenterOfMassExpand = 1;
#endif

    /*
    ** Set the box center to (0,0,0) for now!
    */
    for (j=0;j<6;++j) msr->fCenter[j] = 0.0;
    /*
    ** Process command line arguments.
    */
    ret = prmArgProc(msr->prm,argc,argv);
    if (!ret) {
	_msrExit(msr,1);
	}

    /* This was a checkpoint file! */
    bDoRestore = readParameters(msr,msr->prm->pszFilename);
    if (!bDoRestore) {
	/*
	** Now read parameter file if one was specified.
	** NOTE: command line argument take precedence.
	*/
	if (!prmParseParam(msr->prm,msr)) {
	    _msrExit(msr,1);
	    }
	if (!validateParameters(mdl,msr->prm,&msr->param)) _msrExit(msr,1);
	}

    /*
    ** Convert kboltz/mhydrogen to system units, assuming that
    ** G == 1.
    */
    if(prmSpecified(msr->prm, "dMsolUnit") &&
       prmSpecified(msr->prm, "dKpcUnit")) {
	msr->param.dGasConst = msr->param.dKpcUnit*KPCCM*KBOLTZ
	    /MHYDR/GCGS/msr->param.dMsolUnit/MSOLG;
	/* code energy per unit mass --> erg per g */
	msr->param.dErgPerGmUnit = GCGS*msr->param.dMsolUnit*MSOLG/(msr->param.dKpcUnit*KPCCM);
      /* code energy --> erg IA: TODO be sure about this*/
	msr->param.dErgUnit = GCGS*pow(msr->param.dMsolUnit*MSOLG,2.0)/(msr->param.dKpcUnit*KPCCM);
	/* code density --> g per cc */
	msr->param.dGmPerCcUnit = (msr->param.dMsolUnit*MSOLG)/pow(msr->param.dKpcUnit*KPCCM,3.0);
	/* code time --> seconds */
	msr->param.dSecUnit = sqrt(1/(msr->param.dGmPerCcUnit*GCGS));
	/* code speed --> km/s */
	msr->param.dKmPerSecUnit = sqrt(GCGS*msr->param.dMsolUnit*MSOLG/(msr->param.dKpcUnit*KPCCM))/1e5;
	/* code comoving density --> g per cc = msr->param.dGmPerCcUnit (1+z)^3 */
	msr->param.dComovingGmPerCcUnit = msr->param.dGmPerCcUnit;
	}
    else {
	msr->param.dSecUnit = 1;
	msr->param.dKmPerSecUnit = 1;
	msr->param.dComovingGmPerCcUnit = 1;
	msr->param.dGmPerCcUnit = 1;
	msr->param.dErgPerGmUnit = 1;
      msr->param.dErgUnit = 1;
	}
    msr->param.dTuFac = msr->param.dGasConst/(msr->param.dConstGamma - 1)/
		msr->param.dMeanMolWeight;

    // We convert the minimum dt to rung
    msr->param.dMinDt = ceil( log2(msr->param.dDelta/msr->param.dMinDt) );
#ifdef COOLING
    // We convert the parameters of the entropy floor into code units
    msr->param.dJeansFloorIndex -= 1.; 
    msr->param.dJeansFloorDen *=  MHYDR / (msr->param.dMsolUnit * MSOLG ) / 0.75 ; // Now in rho
    msr->param.dJeansFlooru *= msr->param.dGasConst/(msr->param.dConstGamma - 1.)/1.22; // Now in internal energy per unit mass (assuming neutral gas here with primordial abundances)

    msr->param.dCoolingFloorDen *= MHYDR / (msr->param.dMsolUnit * MSOLG ) / 0.75;
    msr->param.dCoolingFlooru *= msr->param.dGasConst/(msr->param.dConstGamma - 1.)/1.22;

#endif

    /* Gas parameter checks */
#ifdef CLASSICAL_FOPEN
    fprintf(stderr,"WARNING: CLASSICAL_FOPEN\n");
#endif

    /* Determine current opening angle  */
    msr->dThetaMin = msr->param.dTheta;
    if ( !prmSpecified(msr->prm,"nReplicas") && msr->param.nReplicas>=1 ) {
	if ( msr->dThetaMin < 0.52 ) msr->param.nReplicas = 2;
	else msr->param.nReplicas = 1;
	}

    /*
    ** Initialize comove variables.
    */
    msr->nMaxOuts = 100;
    msr->pdOutTime = malloc(msr->nMaxOuts*sizeof(double));
    assert(msr->pdOutTime != NULL);
    msr->nOuts = msr->iOut = 0;


    pstInitialize(&msr->pst,msr->mdl,&msr->lcl);
    pstAddServices(msr->pst,msr->mdl);

    msr->nThreads = mdlThreads(mdl);

    /* Make sure that parallel read and write are sane */
    if (msr->param.nParaRead>msr->nThreads) msr->param.nParaRead = msr->nThreads;
    if (msr->param.nParaWrite>msr->nThreads) msr->param.nParaWrite = msr->nThreads;

    /*
    ** Create the processor subset tree.
    */
    inAdd.idLower = 0;
    inAdd.idUpper = msr->nThreads;
    if (msr->nThreads > 1)
	msrprintf(msr,"Adding %d through %d to the PST\n",
		  inAdd.idLower+1,inAdd.idUpper-1);
    pstSetAdd(msr->pst,&inAdd,sizeof(inAdd),NULL,NULL);

    msr->iCurrMaxRung = 0;
    /*
    ** Mark the Domain Decompositon as not done
    */
    msr->iRungDD = 0;
    msr->iRungDT = 0;
    msr->iLastRungRT = -1;
    msr->iLastRungDD = -1;
    msr->nRung = malloc((MAX_RUNG+1)*sizeof(uint64_t));
    assert(msr->nRung != NULL);
    for (i=0;i<=MAX_RUNG;++i) msr->nRung[i] = 0;

    msr->iRungVeryActive = msr->param.iMaxRung; /* No very active particles */
    msr->bSavePending = 0;                      /* There is no pending save */

    if (msr->param.csm->val.classData.bClass){
        csmClassRead(msr->param.csm, msr->param.dBoxSize, msr->param.h);
        csmClassGslInitialize(msr->param.csm);
    }
    if (strlen(msr->param.csm->val.classData.achLinSpecies) && msr->param.nGridLin == 0){
        fprintf(stderr, "ERROR: you must specify nGridLin when running with linear species\n");
        abort();
    }
    if (msr->param.csm->val.classData.bClass && msr->param.b2LPT){
        fprintf(stderr, "ERROR: 2LPT not yet implemented for Class ICs\n");
        abort();
    }

    return bDoRestore;
    }

void msrLogParams(MSR msr,FILE *fp) {
#if defined(MAXHOSTNAMELEN) && defined(HAVE_GETHOSTNAME)
    char hostname[MAXHOSTNAMELEN];
#endif
    double z;
    int i;

#ifdef __DATE__
#ifdef __TIME__
    fprintf(fp,"# Compiled: %s %s\n",__DATE__,__TIME__);
#endif
#endif
#if defined(__AVX__)
    fprintf(fp,"# with AVX support\n");
#elif defined(__SSE3__)
    fprintf(fp,"# with SSE3 support\n");
#elif defined(__SSE2__)
    fprintf(fp,"# with SSE2 support\n");
#elif defined(__SSE__)
    fprintf(fp,"# with SSE support\n");
#endif

#ifdef USE_BT
#endif
    fprintf(fp,"# Preprocessor macros:");
#ifdef DEBUG
    fprintf(fp," DEBUG");
#endif
#ifdef _REENTRANT
    fprintf(fp," _REENTRANT");
#endif
#if defined(MAXHOSTNAMELEN) && defined(HAVE_GETHOSTNAME)
    fprintf(fp,"\n# Master host: ");
    if (gethostname(hostname,MAXHOSTNAMELEN))
	fprintf(fp,"unknown");
    else
	fprintf(fp,"%s",hostname);
#endif
    fprintf(fp,"\n# N: %"PRIu64,msr->N);
    fprintf(fp," ngas: %"PRIu64,msr->nGas);
    fprintf(fp," nstar: %"PRIu64,msr->nStar);
    fprintf(fp," nThreads: %d",msr->nThreads);
    fprintf(fp," bDiag: %d",msr->param.bDiag);
    fprintf(fp," Verbosity flags: (%d,%d,%d,%d,%d)",msr->param.bVWarnings,
	    msr->param.bVStart,msr->param.bVStep,msr->param.bVRungStat,
	    msr->param.bVDetails);
    fprintf(fp,"\n# bPeriodic: %d",msr->param.bPeriodic);
    fprintf(fp," bComove: %d",msr->param.csm->val.bComove);
    fprintf(fp,"\n# bRestart: %d",msr->param.bRestart);
    fprintf(fp," bParaRead: %d",msr->param.bParaRead);
    fprintf(fp," nParaRead: %d",msr->param.nParaRead);
    fprintf(fp," bParaWrite: %d",msr->param.bParaWrite);
    fprintf(fp," nParaWrite: %d",msr->param.nParaWrite);
    fprintf(fp," bStandard: %d",msr->param.bStandard);
    fprintf(fp," iCompress: %d",msr->param.iCompress);
    fprintf(fp," bHDF5: %d",msr->param.bHDF5);
    fprintf(fp," nBucket: %d",msr->param.nBucket);
    fprintf(fp," nGroup: %d",msr->param.nGroup);
    fprintf(fp," n2min: %d",msr->param.n2min);
    fprintf(fp,"\n# iOutInterval: %d",msr->param.iOutInterval);
    fprintf(fp," iCheckInterval: %d",msr->param.iCheckInterval);
    fprintf(fp," iLogInterval: %d",msr->param.iLogInterval);
    fprintf(fp," iEwOrder: %d",msr->param.iEwOrder);
    fprintf(fp," nReplicas: %d",msr->param.nReplicas);
    fprintf(fp,"\n# dEwCut: %f",msr->param.dEwCut);
    fprintf(fp," dEwhCut: %f",msr->param.dEwhCut);
    fprintf(fp,"\n# iStartStep: %d",msr->param.iStartStep);
    fprintf(fp," nSteps: %d",msr->param.nSteps);
    fprintf(fp," nSmooth: %d",msr->param.nSmooth);
    fprintf(fp," dExtraStore: %f",msr->param.dExtraStore);
    fprintf(fp," nTreeBitsLo: %d",msr->param.nTreeBitsLo);
    fprintf(fp," nTreeBitsHi: %d",msr->param.nTreeBitsHi);
    fprintf(fp," iCacheSize: %d",msr->param.iCacheSize);
    fprintf(fp," iWorkQueueSize: %d",msr->param.iWorkQueueSize);
    fprintf(fp," iCUDAQueueSize: %d",msr->param.iCUDAQueueSize);
    if (prmSpecified(msr->prm,"dSoft"))
	fprintf(fp," dSoft: %g",msr->param.dSoft);
    else
	fprintf(fp," dSoft: input");
    fprintf(fp,"\n# bPhysicalSoft: %d",msr->param.bPhysicalSoft);
    fprintf(fp," nSoftNbr: %d",msr->param.nSoftNbr);
    fprintf(fp," bSoftByType: %d",msr->param.bSoftByType);
    fprintf(fp," bSoftMaxMul: %d",msr->param.bSoftMaxMul);
    fprintf(fp," dSoftMax: %g",msr->param.dSoftMax);
    fprintf(fp," bDoSoftOutput: %d",msr->param.bDoSoftOutput);
    fprintf(fp," bDoAccOutput: %d",msr->param.bDoAccOutput);
    fprintf(fp," bDoPotOutput: %d",msr->param.bDoPotOutput);
    fprintf(fp,"\n# dDelta: %g",msr->param.dDelta);
    fprintf(fp," dEta: %g",msr->param.dEta);
    fprintf(fp," iMaxRung: %d",msr->param.iMaxRung);
    fprintf(fp," nRungVeryActive: %d",msr->param.nRungVeryActive);
    fprintf(fp," bDoRungOutput: %d",msr->param.bDoRungOutput);
    fprintf(fp," bDoRungDestOutput: %d",msr->param.bDoRungDestOutput);
    fprintf(fp,"\n# bGravStep: %d",msr->param.bGravStep);
    fprintf(fp," bEpsAccStep: %d",msr->param.bEpsAccStep);
    fprintf(fp," bDensityStep: %d",msr->param.bDensityStep);
    fprintf(fp," nTruncateRung: %d",msr->param.nTruncateRung);
    fprintf(fp,"\n# iTimeStepCrit: %d",msr->param.iTimeStepCrit);
    fprintf(fp," nPartRhoLoc: %d", msr->param.nPartRhoLoc);
    fprintf(fp," dPreFacRhoLoc: %g", msr->param.dPreFacRhoLoc);
    fprintf(fp," dFacExcludePart: %g", msr->param.dFacExcludePart);
    fprintf(fp," dEccFacMax: %g", msr->param.dEccFacMax);
    fprintf(fp," nPartColl: %d", msr->param.nPartColl);
    fprintf(fp,"\n# bDoGravity: %d",msr->param.bDoGravity);
    fprintf(fp," bAarsethStep: %d",msr->param.bAarsethStep);
    fprintf(fp,"\n# dFracDualTree: %g",msr->param.dFracDualTree);
    fprintf(fp,"dFracNoDomainDecomp: %g",msr->param.dFracNoDomainDecomp);
    fprintf(fp," dFracNoDomainRootFind: %g",msr->param.dFracNoDomainRootFind);
    fprintf(fp," dFracNoDomainDimChoice: %g",msr->param.dFracNoDomainDimChoice);
    fprintf(fp,"\n# nTruncateRung: %d",msr->param.nTruncateRung);
    fprintf(fp," dGrowDeltaM: %g",msr->param.dGrowDeltaM);
    fprintf(fp," dGrowStartT: %g",msr->param.dGrowStartT);
    fprintf(fp," dGrowEndT: %g",msr->param.dGrowEndT);

    fprintf(fp,"\n# SPH: bDoGas: %d",msr->param.bDoGas);	
    fprintf(fp," bGasAdiabatic: %d",msr->param.bGasAdiabatic);	
    fprintf(fp," bGasIsothermal: %d",msr->param.bGasIsothermal);	
    fprintf(fp," bGasCooling: %d",msr->param.bGasCooling);	
    fprintf(fp," bInitTFromCooling: %d",msr->param.bInitTFromCooling);	
    fprintf(fp," iRungCoolTableUpdate: %d",msr->param.iRungCoolTableUpdate);
    fprintf(fp," iViscosityLimiter: %d",msr->param.iViscosityLimiter);	
    fprintf(fp," iDiffusion: %d",msr->param.iDiffusion);	
    fprintf(fp,"\n# dConstAlpha: %g",msr->param.dConstAlpha);
    fprintf(fp," dConstBeta: %g",msr->param.dConstBeta);
    fprintf(fp," dConstGamma: %g",msr->param.dConstGamma);
    fprintf(fp," dMeanMolWeight: %g",msr->param.dMeanMolWeight);
    fprintf(fp," dGasConst: %g",msr->param.dGasConst);
    fprintf(fp,"\n# dEtaCourant: %g",msr->param.dEtaCourant);
    fprintf(fp," dEtaUDot: %g",msr->param.dEtaUDot);
    fprintf(fp," dTuFac: %g",msr->param.dTuFac);
    fprintf(fp," dhMinOverSoft: %g",msr->param.dhMinOverSoft);
    fprintf(fp," dMetalDiffusionCoeff: %g",msr->param.dMetalDiffusionCoeff);
    fprintf(fp," dThermalDiffusionCoeff: %g",msr->param.dThermalDiffusionCoeff);
    fprintf(fp,"\n# UNITS: dKBoltzUnit: %g",msr->param.dKBoltzUnit);
    fprintf(fp," dMsolUnit: %g",msr->param.dMsolUnit);
    fprintf(fp," dKpcUnit: %g",msr->param.dKpcUnit);
    if(prmSpecified(msr->prm, "dMsolUnit") &&
       prmSpecified(msr->prm, "dKpcUnit")) {
	fprintf(fp," dErgPerGmUnit: %g", msr->param.dErgPerGmUnit );
	fprintf(fp," dGmPerCcUnit (z=0): %g", msr->param.dGmPerCcUnit );
	fprintf(fp," dSecUnit: %g", msr->param.dSecUnit );
	fprintf(fp," dKmPerSecUnit (z=0): %g", msr->param.dKmPerSecUnit );
	}
    fprintf(fp,"\n# STARFORM: bStarForm %d",msr->param.bStarForm);
    fprintf(fp," bFeedback %d",msr->param.bFeedback);
    fprintf(fp," SFdEfficiency %g",msr->param.SFdEfficiency);
    fprintf(fp," SFdTMax %g",msr->param.SFdTMax);
    fprintf(fp," SFdPhysDenMin %g",msr->param.SFdPhysDenMin);
    fprintf(fp," SFdComovingDenMin %g",msr->param.SFdComovingDenMin);
    fprintf(fp," SFdESNPerStarMass %g",msr->param.SFdESNPerStarMass);
    fprintf(fp,"\n# SFdtCoolingShutoff %g",msr->param.SFdtCoolingShutoff);
    fprintf(fp," SFdtFeedbackDelay %g",msr->param.SFdtFeedbackDelay);
    fprintf(fp," SFdMassLossPerStarMass %g",msr->param.SFdMassLossPerStarMass);
    fprintf(fp," SFdZMassPerStarMass %g",msr->param.SFdZMassPerStarMass);
    fprintf(fp," SFdInitStarMass %g",msr->param.SFdInitStarMass);
    fprintf(fp," SFdMinGasMass %g",msr->param.SFdMinGasMass);
    fprintf(fp," SFbdivv %d",msr->param.SFbdivv);
    /* -- */
    fprintf(fp,"\n# Group Find: bFindHopGroups: %d",msr->param.bFindHopGroups);
    fprintf(fp," dHopTau: %g",msr->param.dHopTau);
    /* -- */
    fprintf(fp,"\n# Group Find: bFindGroups: %d",msr->param.bFindGroups);
    fprintf(fp," dTau: %g",msr->param.dTau);
    fprintf(fp," nMinMembers: %d",msr->param.nMinMembers);
    fprintf(fp," nBins: %d",msr->param.nBins);
    fprintf(fp," iCenterType: %d",msr->param.iCenterType);
    fprintf(fp," binFactor: %g",msr->param.binFactor);
    fprintf(fp," fMinRadius: %g",msr->param.fMinRadius);
    fprintf(fp," bLogBins: %d",msr->param.bLogBins);
    fprintf(fp,"\n# Relaxation estimate: bTraceRelaxation: %d",msr->param.bTraceRelaxation);
    fprintf(fp," dTheta: %f",msr->param.dTheta);
    fprintf(fp,"\n# dPeriod: %g",msr->param.dPeriod);
    fprintf(fp," dxPeriod: %g",
	    msr->param.dxPeriod >= FLOAT_MAXVAL ? 0 : msr->param.dxPeriod);
    fprintf(fp," dyPeriod: %g",
	    msr->param.dyPeriod >= FLOAT_MAXVAL ? 0 : msr->param.dyPeriod);
    fprintf(fp," dzPeriod: %g",
	    msr->param.dzPeriod >= FLOAT_MAXVAL ? 0 : msr->param.dzPeriod);
    fprintf(fp,"\n# dHubble0: %g",msr->param.csm->val.dHubble0);
    fprintf(fp," dOmega0: %g",msr->param.csm->val.dOmega0);
    fprintf(fp," dLambda: %g",msr->param.csm->val.dLambda);
    fprintf(fp," dOmegaDE: %g",msr->param.csm->val.dOmegaDE);
    fprintf(fp," w0: %g",msr->param.csm->val.w0);
    fprintf(fp," wa: %g",msr->param.csm->val.wa);
    fprintf(fp," dOmegaRad: %g",msr->param.csm->val.dOmegaRad);
    fprintf(fp," dOmegab: %g",msr->param.csm->val.dOmegab);
    fprintf(fp,"\n# achInFile: %s",msr->param.achInFile);
    fprintf(fp,"\n# achOutName: %s",msr->param.achOutName);
    fprintf(fp,"\n# achOutPath: %s",msr->param.achOutPath);
    fprintf(fp,"\n# achIoPath: %s",msr->param.achIoPath);
    fprintf(fp,"\n# achDataSubPath: %s",msr->param.achDataSubPath);
    if (msr->param.csm->val.bComove) {
	fprintf(fp,"\n# RedOut:");
	if (msr->nOuts == 0) fprintf(fp," none");
	for (i=0;i<msr->nOuts;i++) {
	    if (i%5 == 0) fprintf(fp,"\n#   ");
	    z = 1.0/csmTime2Exp(msr->param.csm, msr->pdOutTime[i]) - 1.0;
	    fprintf(fp," %f",z);
	    }
	fprintf(fp,"\n");
	}
    else {
	fprintf(fp,"\n# TimeOut:");
	if (msr->nOuts == 0) fprintf(fp," none");
	for (i=0;i<msr->nOuts;i++) {
	    if (i%5 == 0) fprintf(fp,"\n#   ");
	    fprintf(fp," %f",msr->pdOutTime[i]);
	    }
	fprintf(fp,"\n");
	}
    }

int
msrGetLock(MSR msr) {
    /*
    ** Attempts to lock run directory to prevent overwriting. If an old lock
    ** is detected with the same achOutName, an abort is signaled. Otherwise
    ** a new lock is created. The bOverwrite parameter flag can be used to
    ** suppress lock checking.
    */

    FILE *fp = NULL;
    char achTmp[256],achFile[256];

    _msrMakePath(msr->param.achDataSubPath,LOCKFILE,achFile);
    if (!msr->param.bOverwrite && (fp = fopen(achFile,"r"))) {
	if (fscanf(fp,"%s",achTmp) != 1) achTmp[0] = '\0';
	(void) fclose(fp);
	if (!strcmp(msr->param.achOutName,achTmp)) {
	    (void) printf("ABORT: %s detected.\nPlease ensure data is safe to "
			  "overwrite. Delete lockfile and try again.\n",achFile);
	    return 0;
	    }
	}
    if (!(fp = fopen(achFile,"w"))) {
	if (msr->param.bOverwrite && msr->param.bVWarnings) {
	    (void) printf("WARNING: Unable to create %s...ignored.\n",achFile);
	    return 1;
	    }
	else {
	    (void) printf("Unable to create %s\n",achFile);
	    return 0;
	    }
	}
    (void) fprintf(fp,"%s",msr->param.achOutName);
    (void) fclose(fp);
    return 1;
    }

int msrCheckForStop(MSR msr,const char *achStopFile) {
    /*
    ** Checks for existence of STOPFILE in run directory. If found, the file
    ** is removed and the return status is set to 1, otherwise 0.
    */

    char achFile[256];
    FILE *fp = NULL;
    _msrMakePath(msr->param.achDataSubPath,achStopFile,achFile);
    if ((fp = fopen(achFile,"r"))) {
	(void) printf("User interrupt detected.\n");
	(void) fclose(fp);
	(void) unlink(achFile);
	return 1;
	}
    return 0;
    }

void msrFinish(MSR msr) {
   int id;
   
   for (id=1;id<msr->nThreads;++id) {
	int rID;
	rID = mdlReqService(msr->mdl,id,SRV_STOP,NULL,0);
	mdlGetReply(msr->mdl,rID,NULL,NULL);
	}
    pstFinish(msr->pst);
    csmFinish(msr->param.csm);
    /*
    ** finish with parameter stuff, deallocate and exit.
    */
    prmFinish(msr->prm);
    free(msr->nRung);
    free(msr->pdOutTime);
    free(msr);
    }

static int CmpPC(const void *v1,const void *v2) {
    PARTCLASS *pClass1 = (PARTCLASS*)v1;
    PARTCLASS *pClass2 = (PARTCLASS*)v2;
    if ( pClass1->fMass < pClass2->fMass ) return -1;
    else if ( pClass1->fMass > pClass2->fMass ) return 1;
    else if ( pClass1->fSoft < pClass2->fSoft ) return -1;
    else if ( pClass1->fSoft > pClass2->fSoft ) return 1;
    else return 0;
    }

void msrSetClasses(MSR msr) {
    PARTCLASS *pClass;
    int n, nClass;
    pClass = malloc(PKD_MAX_CLASSES*sizeof(PARTCLASS));
    assert(pClass!=NULL);
    pstGetClasses(msr->pst,NULL,0,pClass,&nClass);
    n = nClass / sizeof(PARTCLASS);
    assert(n*sizeof(PARTCLASS)==nClass);
    qsort(pClass,n,sizeof(PARTCLASS),CmpPC);
    pstSetClasses(msr->pst,pClass,nClass,NULL,NULL);
    free(pClass);
    }

static void _SwapClasses(MSR msr, int id) {
    LCL *plcl = msr->pst->plcl;
    PST pst0 = msr->pst;
    PARTCLASS *pClass;
    int n;
    int rID;

    pClass = malloc(PKD_MAX_CLASSES*sizeof(PARTCLASS));
    assert(pClass!=NULL);

    n = pkdGetClasses( plcl->pkd, PKD_MAX_CLASSES, pClass );
    rID = mdlReqService(pst0->mdl,id,PST_SWAPCLASSES,pClass,n*sizeof(PARTCLASS));
    mdlGetReply(pst0->mdl,rID,pClass,&n);
    n = n / sizeof(PARTCLASS);
    pkdSetClasses( plcl->pkd, n, pClass, 0 );
    free(pClass);
    }

void msrOneNodeRead(MSR msr, struct inReadFile *in, FIO fio) {
    int id;
    int *nParts;		/* number of particles for each processor */
    uint64_t nStart;
    PST pst0;
    LCL *plcl;
    char achInFile[PST_FILENAME_SIZE];
    int nid;
    int inswap;
    int rID;

    nParts = malloc(msr->nThreads*sizeof(*nParts));
    for (id=0;id<msr->nThreads;++id) {
	nParts[id] = -1;
	}

    pstOneNodeReadInit(msr->pst, in, sizeof(*in), nParts, &nid);
    assert((size_t)nid == msr->nThreads*sizeof(*nParts));
    for (id=0;id<msr->nThreads;++id) {
	assert(nParts[id] > 0);
	}

    pst0 = msr->pst;
    while (pst0->nLeaves > 1)
	pst0 = pst0->pstLower;
    plcl = pst0->plcl;

    nStart = nParts[0];
    for (id=1;id<msr->nThreads;++id) {
	/*
	 * Read particles into the local storage.
	 */
	assert(plcl->pkd->nStore >= nParts[id]);
	pkdReadFIO(plcl->pkd, fio, nStart, nParts[id], in->dvFac,in->dTuFac);
	nStart += nParts[id];
	/*
	 * Now shove them over to the remote processor.
	 */
	_SwapClasses(msr,id);
	inswap = 0;
	rID = mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
	pkdSwapAll(plcl->pkd, id);
	mdlGetReply(pst0->mdl,rID,NULL,NULL);
	}
    assert(nStart == msr->N);
    /*
     * Now read our own particles.
     */
    pkdReadFIO(plcl->pkd, fio, 0, nParts[0], in->dvFac, in->dTuFac);

    free(nParts);
    }

double getTime(MSR msr, double dExpansion, double *dvFac) {
    double dTime,aTo,tTo,z;
    if (msr->param.csm->val.bComove) {
	if (msr->param.csm->val.dHubble0 == 0.0) {
	    printf("No hubble constant specified\n");
	    _msrExit(msr,1);
	    }
	dTime = csmExp2Time(msr->param.csm,dExpansion);
	z = 1.0/dExpansion - 1.0;
	if (msr->param.bVStart)
	    printf("Input file, Time:%g Redshift:%g Expansion factor:%g iStartStep:%d\n",
		   dTime,z,dExpansion,msr->param.iStartStep);
	if (prmSpecified(msr->prm,"dRedTo")) {
	    if (msr->param.dRedTo <= -1.0) {
		printf("Badly specified final redshift (zTo <= -1.0), check -zto parameter.\n");
		_msrExit(msr,1);
		}
	    if (!prmArgSpecified(msr->prm,"nSteps") &&
		    prmArgSpecified(msr->prm,"dDelta")) {
		aTo = 1.0/(msr->param.dRedTo + 1.0);
		tTo = csmExp2Time(msr->param.csm,aTo);
		if (msr->param.bVStart)
		    printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
			   tTo,1.0/aTo-1.0,aTo);
		if (tTo < dTime) {
		    printf("Badly specified final redshift, check -zto parameter.\n");
		    _msrExit(msr,1);
		    }
		msr->param.nSteps = (int)ceil((tTo-dTime)/msr->param.dDelta);
		msr->param.dDelta =
		    (tTo-dTime)/(msr->param.nSteps -
				 msr->param.iStartStep);
		}
	    else if (!prmArgSpecified(msr->prm,"dDelta") &&
		     prmArgSpecified(msr->prm,"nSteps")) {
		aTo = 1.0/(msr->param.dRedTo + 1.0);
		tTo = csmExp2Time(msr->param.csm,aTo);
		if (msr->param.bVStart)
		    printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
			   tTo,1.0/aTo-1.0,aTo);
		if (tTo < dTime) {
		    printf("Badly specified final redshift, check -zto parameter.\n");
		    _msrExit(msr,1);
		    }
		if (msr->param.nSteps != 0)
		    msr->param.dDelta =
			(tTo-dTime)/(msr->param.nSteps -
				     msr->param.iStartStep);

		else
		    msr->param.dDelta = 0.0;
		}
	    else if (!prmSpecified(msr->prm,"nSteps") &&
		     prmFileSpecified(msr->prm,"dDelta")) {
		aTo = 1.0/(msr->param.dRedTo + 1.0);
		tTo = csmExp2Time(msr->param.csm,aTo);
		if (msr->param.bVStart)
		    printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
			   tTo,1.0/aTo-1.0,aTo);
		if (tTo < dTime) {
		    printf("Badly specified final redshift, check -zto parameter.\n");
		    _msrExit(msr,1);
		    }
		msr->param.nSteps = (int)ceil((tTo-dTime)/msr->param.dDelta);
		msr->param.dDelta =
		    (tTo-dTime)/(msr->param.nSteps -
				 msr->param.iStartStep);
		}
	    else if (!prmSpecified(msr->prm,"dDelta") &&
		     prmFileSpecified(msr->prm,"nSteps")) {
		aTo = 1.0/(msr->param.dRedTo + 1.0);
		tTo = csmExp2Time(msr->param.csm,aTo);
		if (msr->param.bVStart)
		    printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
			   tTo,1.0/aTo-1.0,aTo);
		if (tTo < dTime) {
		    printf("Badly specified final redshift, check -zto parameter.\n");
		    _msrExit(msr,1);
		    }
		if (msr->param.nSteps != 0)
		    msr->param.dDelta =	(tTo-dTime)/(msr->param.nSteps
						     - msr->param.iStartStep);
		else
		    msr->param.dDelta = 0.0;
		}
	    }
	else {
	    tTo = dTime + msr->param.nSteps*msr->param.dDelta;
	    aTo = csmTime2Exp(msr->param.csm,tTo);
	    if (msr->param.bVStart)
		printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
		       tTo,1.0/aTo-1.0,aTo);
	    }
	if (msr->param.csm->val.bComove) {
	    *dvFac = dExpansion*dExpansion;
	    }
	else {
	    *dvFac = 1.0;
	    }
	}
    else {
	dTime = dExpansion;
	if (msr->param.bVStart) printf("Input file, Time:%g iStartStep:%d\n",dTime,msr->param.iStartStep);
	tTo = dTime + (msr->param.nSteps - msr->param.iStartStep)*msr->param.dDelta;
	if (msr->param.bVStart) {
	    printf("Simulation to Time:%g\n",tTo);
	    }
	*dvFac = 1.0;
	}

    return dTime;
    }

/*
** This function makes some potentially problematic assumptions!!!
** Main problem is that it calls pkd level routines, bypassing the
** pst level. It uses plcl pointer which is not desirable.
*/
void msrAllNodeWrite(MSR msr, const char *pszFileName, double dTime, double dvFac, int bDouble) {
    int nProcessors;
    PST pst0;
    LCL *plcl;
    struct inWrite in;
    int j;

    pst0 = msr->pst;
    while (pst0->nLeaves > 1)
	pst0 = pst0->pstLower;
    plcl = pst0->plcl;


    /*
    ** Add Data Subpath for local and non-local names.
    */
    _msrMakePath(msr->param.achDataSubPath,pszFileName,in.achOutFile);

    in.bStandard = msr->param.bStandard;
    /*
    ** If bParaWrite is 0, then we write serially; if it is 1, then we write
    ** in parallel using all available threads, otherwise we write in parallel
    ** using the specified number of threads.  The latter option will reduce
    ** the total amount of simultaneous I/O for file systems that cannot
    ** handle it.
    */
    nProcessors = msr->param.bParaWrite==0?1:(msr->param.nParaWrite<=1 ? msr->nThreads:msr->param.nParaWrite);
    in.iIndex = 0;

    if (msr->param.csm->val.bComove) {
	in.dTime = csmTime2Exp(msr->param.csm,dTime);
	in.dvFac = 1.0/(in.dTime*in.dTime);
	}
    else {
	in.dTime = dTime;
	in.dvFac = 1.0;
	}

    /* We need to enforce periodic boundaries (when applicable) */
    if (msr->param.bPeriodic &&
	    msr->param.dxPeriod < FLOAT_MAXVAL &&
	    msr->param.dyPeriod < FLOAT_MAXVAL &&
	    msr->param.dzPeriod < FLOAT_MAXVAL) {
	for (j=0;j<3;++j) {
	    in.bnd.fCenter[j] = msr->fCenter[j];
	    }
	in.bnd.fMax[0] = 0.5*msr->param.dxPeriod;
	in.bnd.fMax[1] = 0.5*msr->param.dyPeriod;
	in.bnd.fMax[2] = 0.5*msr->param.dzPeriod;
	}
    else {
	for (j=0;j<3;++j) {
	    in.bnd.fCenter[j] = 0.0;
	    in.bnd.fMax[j] = FLOAT_MAXVAL;
	    }
	}

    in.dEcosmo    = msr->dEcosmo;
    in.dTimeOld   = msr->dTimeOld;
    in.dUOld      = msr->dUOld;
    in.dBoxSize   = msr->param.dBoxSize;
    in.Omega0     = msr->param.csm->val.dOmega0;
    in.OmegaLambda= msr->param.csm->val.dLambda;
    in.HubbleParam= msr->param.h;

    in.nDark = msr->nDark;
    in.nSph  = msr->nGas;
    in.nStar = msr->nStar;

    in.bHDF5 = msr->param.bHDF5;
    in.mFlags = FIO_FLAG_POTENTIAL | FIO_FLAG_DENSITY
	| (bDouble?FIO_FLAG_CHECKPOINT:0)
	| (msr->param.bDoublePos?FIO_FLAG_DOUBLE_POS:0)
	| (msr->param.bDoubleVel?FIO_FLAG_DOUBLE_VEL:0)
	| (msr->param.bMemMass?0:FIO_FLAG_COMPRESS_MASS)
	| (msr->param.bMemSoft?0:FIO_FLAG_COMPRESS_SOFT);

    if (!msr->param.bHDF5 && strstr(in.achOutFile,"&I")==0) {
	FIO fio;
	fio = fioTipsyCreate(in.achOutFile,
			     in.mFlags&FIO_FLAG_CHECKPOINT,
			     in.bStandard,in.dTime,
			     in.nSph, in.nDark, in.nStar);
	fioClose(fio);
	}
    in.iLower = 0;
    in.iUpper = msr->nThreads;
    in.iIndex = 0;
    in.nProcessors = nProcessors;
    pstWrite(msr->pst,&in,sizeof(in),NULL,NULL);
    }


uint64_t msrCalcWriteStart(MSR msr) {
    struct outSetTotal out;
    struct inSetWriteStart in;

    pstSetTotal(msr->pst,NULL,0,&out,NULL);
    assert(out.nTotal <= msr->N);
    in.nWriteStart = 0;
    pstSetWriteStart(msr->pst,&in,sizeof(in),NULL,NULL);
    return out.nTotal;
    }

void msrWrite(MSR msr,const char *pszFileName,double dTime,int bCheckpoint) {
    char achOutFile[PST_FILENAME_SIZE];
    LCL *plcl = msr->pst->plcl;
    int nProcessors;
    double dvFac, dExp;
    double sec,dsec;
    uint64_t N;

#ifdef NAIVE_DOMAIN_DECOMP
    msrReorder(msr);
#else
    if ( msr->iLastRungRT >= 0 ) msrReorder(msr);
#endif

    /*
    ** Calculate where to start writing.
    ** This sets plcl->nWriteStart.
    */
    N = msrCalcWriteStart(msr);
    /*
    ** Add Data Subpath for local and non-local names.
    */
    _msrMakePath(msr->param.achDataSubPath,pszFileName,achOutFile);

    /*
    ** If bParaWrite is 0, then we write serially; if it is 1, then we write
    ** in parallel using all available threads, otherwise we write in parallel
    ** using the specified number of threads.  The latter option will reduce
    ** the total amount of simultaneous I/O for file systems that cannot
    ** handle it.
    */
    nProcessors = msr->param.bParaWrite==0?1:(msr->param.nParaWrite<=1 ? msr->nThreads:msr->param.nParaWrite);

    if (msr->param.csm->val.bComove) {
	dExp = csmTime2Exp(msr->param.csm,dTime);
	dvFac = 1.0/(dExp*dExp);
	}
    else {
	dExp = dTime;
	dvFac = 1.0;
	}
    if ( !msr->param.bParaWrite ) {
	msrprintf(msr,"Writing %s in %s format serially ...\n",
		  achOutFile, (msr->param.bHDF5?"HDF5":"Tipsy"));
	}
    else {
	if ( msr->param.nParaWrite > 1 )
	    msrprintf(msr,"Writing %s in %s format in parallel (but limited to %d processors) ...\n",
		      achOutFile, (msr->param.bHDF5?"HDF5":"Tipsy"), nProcessors);
	else
	    msrprintf(msr,"Writing %s in %s format in parallel ...\n",
		      achOutFile, (msr->param.bHDF5?"HDF5":"Tipsy"));
	}

    if (msr->param.csm->val.bComove)
	msrprintf(msr,"Time:%g Redshift:%g\n",dTime,(1.0/dExp - 1.0));
    else
	msrprintf(msr,"Time:%g\n",dTime);

    sec = msrTime();
    msrAllNodeWrite(msr, achOutFile, dTime, dvFac, bCheckpoint);
    dsec = msrTime() - sec;

    msrprintf(msr,"Output file has been successfully written, Wallclock: %f secs.\n", dsec);
    }


void msrSetSoft(MSR msr,double dSoft) {
    struct inSetSoft in;

    in.dSoft = dSoft;
    pstSetSoft(msr->pst,&in,sizeof(in),NULL,NULL);
    }

// IA: If the initial condition do not provide information on the smoothing
//    length of the gas particles, we set it assuming an equipartition of the whole
//    volume. This will be the first estimate before the iterative computation of fBall
void msrSetSmooth(MSR msr) {
    struct inSetSoft in;  // For simplicty we reuse this struct

    double V = msr->param.dxPeriod * msr->param.dyPeriod * msr->param.dzPeriod;
    double dSmooth = pow( 3.*V/(4.*M_PI*msr->N), 1./3. );

    printf("Setting initial smoothing length to %e \n", dSmooth);

    in.dSoft = dSmooth;
    pstSetSmooth(msr->pst,&in,sizeof(in),NULL,NULL);
    }

void msrDomainDecompOld(MSR msr,int iRung,int bSplitVA) {
    struct inDomainDecomp in;
    uint64_t nActive;
    const uint64_t nDT = d2u64(msr->N*msr->param.dFracDualTree);
    const uint64_t nDD = d2u64(msr->N*msr->param.dFracNoDomainDecomp);
    const uint64_t nRT = d2u64(msr->N*msr->param.dFracNoDomainRootFind);
    const uint64_t nSD = d2u64(msr->N*msr->param.dFracNoDomainDimChoice);
    double sec,dsec;
    int iRungDT, iRungDD=0,iRungRT,iRungSD;
    int i,j;
    int bRestoreActive = 0;

    in.bDoRootFind = 1;
    in.bDoSplitDimFind = 1;
    if (iRung >= 0) {
	/*
	** All of this could be calculated once for the case that the number
	** of particles don't change. Or calculated every time the number of
	** particles does change.
	*/
	nActive = 0;
	iRungDT = 0;
	iRungDD = 0;
	iRungRT = 0;
	iRungSD = 0;
	for (i=msr->iCurrMaxRung;i>=0;--i) {
	    nActive += msr->nRung[i];
	    if (nActive > nDT && !iRungDT) iRungDT = i;
	    if (nActive > nDD && !iRungDD) iRungDD = i;
	    if (nActive > nRT && !iRungRT) iRungRT = i;
	    if (nActive > nSD && !iRungSD) iRungSD = i;
	    }
	assert(iRungDD >= iRungRT);
	assert(iRungRT >= iRungSD);
	msr->iRungDD = iRungDD;
	msr->iRungDT = iRungDT;
#ifdef NAIVE_DOMAIN_DECOMP
	if (msr->iLastRungRT < 0) {
	    /*
	    ** We need to do a full domain decompotition with iRungRT particles being active.
	    ** However, since I am not sure what the exact state of the domains can be at this point
	    ** I had better do a split dim find as well.
	    */
	    msr->iLastRungRT = 0;
	    msrActiveRung(msr,msr->iLastRungRT,1);
	    bRestoreActive = 1;
	    in.bDoRootFind = 1;
	    in.bDoSplitDimFind = 1;
	    if (msr->param.bVRungStat) {
		printf("Doing Domain Decomposition (nActive = %"PRIu64"/%"PRIu64", iRung:%d iRungRT:%d)\n",
		    msr->nActive,msr->N,iRung,iRungRT);
		}
	    }
	else if (iRung <= iRungRT) {
	    /*
	    ** We need to do a full domain decomposition with *ALL* particles being active.
	    */
	    in.bDoRootFind = 1;
	    if (iRung <= iRungSD) {
		if (msr->param.bVRungStat) {
		    printf("Doing Domain Decomposition (nActive = %"PRIu64"/%"PRIu64", iRung:%d iRungRT:%d)\n",
			msr->nActive,msr->N,iRung,iRungRT);
		    }
		in.bDoSplitDimFind = 1;
		}
	    else { 
		if (msr->param.bVRungStat) {
		    printf("Skipping Domain Dim Choice (nActive = %"PRIu64"/%"PRIu64", iRung:%d iRungSD:%d)\n",
			msr->nActive,msr->N,iRung,iRungSD);
		    }
		in.bDoSplitDimFind = 0;
		}
	    msrActiveRung(msr,0,1); /* Here we activate all particles. */
	    bRestoreActive = 1;
	    }	    
	else if (iRung <= iRungDD) {
	    if (msr->param.bVRungStat) {
		printf("Skipping Root Finder (nActive = %"PRIu64"/%"PRIu64", iRung:%d iRungRT:%d iRungDD:%d)\n",
		    msr->nActive,msr->N,iRung,iRungRT,iRungDD);
		}
	    in.bDoRootFind = 0;
	    in.bDoSplitDimFind = 0;
	    bRestoreActive = 0;	    
	    }
	else {
	    if (msr->param.bVRungStat) {
		printf("Skipping Domain Decomposition (nActive = %"PRIu64"/%"PRIu64", iRung:%d iRungDD:%d)\n",
		    msr->nActive,msr->N,iRung,iRungDD);
		}
	    return; /* do absolutely nothing! */
	    }
#else
	if (msr->iLastRungRT < 0) {
	    /*
	    ** We need to do a full domain decompotition with iRungRT particles being active.
	    ** However, since I am not sure what the exact state of the domains can be at this point
	    ** I had better do a split dim find as well.
	    */
	    msr->iLastRungRT = iRungRT;
	    msrActiveRung(msr,iRungRT,1);
	    bRestoreActive = 1;
	    in.bDoRootFind = 1;
	    in.bDoSplitDimFind = 1;
	    }
	else if (iRung == msr->iLastRungDD) {
	    if (msr->param.bVRungStat) {
		printf("Skipping Domain Decomposition (nActive = %"PRIu64"/%"PRIu64", iRung:%d iRungDD:%d iLastRungRT:%d)\n",
		    msr->nActive,msr->N,iRung,iRungDD,msr->iLastRungRT);
		}
	    return;  /* do absolutely nothing! */
	    }
	else if (iRung >= iRungDD && !bSplitVA) {
	    if (msr->iLastRungRT < iRungRT) {
		msr->iLastRungRT = iRungRT;
		msrActiveRung(msr,iRungRT,1);
		bRestoreActive = 1;
		in.bDoRootFind = 1;
		in.bDoSplitDimFind = 0;
		}
	    else {
		if (msr->param.bVRungStat) {
		    printf("Skipping Domain Decomposition (nActive = %"PRIu64"/%"PRIu64", iRung:%d iRungDD:%d iLastRungRT:%d)\n",
			msr->nActive,msr->N,iRung,iRungDD,msr->iLastRungRT);
		    }
		return;  /* do absolutely nothing! */
		}
	    }
	else if (iRung > iRungRT) {
	    if (msr->iLastRungRT < iRungRT) {
		msr->iLastRungRT = iRungRT;
		msrActiveRung(msr,iRungRT,1);
		bRestoreActive = 1;
		in.bDoRootFind = 1;
		in.bDoSplitDimFind = 0;
		}
	    else {
		if (msr->param.bVRungStat) {
		    printf("Skipping Root Finder (nActive = %"PRIu64"/%"PRIu64", iRung:%d iRungRT:%d iRungDD:%d iLastRungRT:%d)\n",
			msr->nActive,msr->N,iRung,iRungRT,iRungDD,msr->iLastRungRT);
		    }
		in.bDoRootFind = 0;
		in.bDoSplitDimFind = 0;
		}
	    }
	else if (iRung > iRungSD) {
	    if (msr->iLastRungRT == iRung) {
		if (msr->param.bVRungStat) {
		    printf("Skipping Root Finder (nActive = %"PRIu64"/%"PRIu64", iRung:%d iRungRT:%d iRungDD:%d iLastRungRT:%d)\n",
			msr->nActive,msr->N,iRung,iRungRT,iRungDD,msr->iLastRungRT);
		    }
		in.bDoRootFind = 0;
		in.bDoSplitDimFind = 0;
		}
	    else {
		if (msr->param.bVRungStat) {
		    printf("Skipping Domain Dim Choice (nActive = %"PRIu64"/%"PRIu64", iRung:%d iRungSD:%d iLastRungRT:%d)\n",
			msr->nActive,msr->N,iRung,iRungSD,msr->iLastRungRT);
		    }
		msr->iLastRungRT = iRung;
		in.bDoRootFind = 1;
		in.bDoSplitDimFind = 0;
		}
	    }
	else {
	    if (msr->iLastRungRT == iRung) {
		in.bDoRootFind = 0;
		in.bDoSplitDimFind = 0;
		}
	    else {
		msr->iLastRungRT = iRung;
		in.bDoRootFind = 1;
		in.bDoSplitDimFind = 1;
		}
	    }
#endif
	}
    msr->iLastRungDD = msr->iLastRungRT;
    in.nActive = msr->nActive;
    in.nTotal = msr->N;

    in.nBndWrap[0] = 0;
    in.nBndWrap[1] = 0;
    in.nBndWrap[2] = 0;
    /*
    ** If we are dealing with a nice periodic volume in all
    ** three dimensions then we can set the initial bounds
    ** instead of calculating them.
    */
    if (msr->param.bPeriodic &&
	msr->param.dxPeriod < FLOAT_MAXVAL &&
	msr->param.dyPeriod < FLOAT_MAXVAL &&
	msr->param.dzPeriod < FLOAT_MAXVAL) {
	for (j=0;j<3;++j) {
	    in.bnd.fCenter[j] = msr->fCenter[j];
	    }
	in.bnd.fMax[0] = 0.5*msr->param.dxPeriod;
	in.bnd.fMax[1] = 0.5*msr->param.dyPeriod;
	in.bnd.fMax[2] = 0.5*msr->param.dzPeriod;

	pstEnforcePeriodic(msr->pst,&in.bnd,sizeof(BND),NULL,NULL);
	}
    else {
	pstCombineBound(msr->pst,NULL,0,&in.bnd,NULL);
	}
    /*
    ** If we are doing SPH we need to make absolutely certain to clear
    ** all neighbor lists here since the pointers in the particles will
    ** only be valid on the node where it was malloc'ed!
    */
#ifdef FAST_GAS
    if (msr->param.bDoGas) {
	pstFastGasCleanup(msr->pst,NULL,0,NULL,NULL);
	}
#endif
    in.bSplitVA = bSplitVA;
    msrprintf(msr,"Domain Decomposition: nActive (Rung %d) %"PRIu64" SplitVA:%d\n",
	msr->iLastRungRT,msr->nActive,bSplitVA);
    msrprintf(msr,"Domain Decomposition... \n");
    sec = msrTime();

    pstDomainDecomp(msr->pst,&in,sizeof(in),NULL,NULL);
    dsec = msrTime() - sec;
    printf("Domain Decomposition complete, Wallclock: %f secs\n\n",dsec);
    if (bRestoreActive) {
	/* Restore Active data */
	msrActiveRung(msr,iRung,1);
	}
    }

void msrInflate(MSR msr,int iStep) {
    if (msr->param.nInflateReps>0 && iStep==msr->param.iInflateStep) {
	struct inInflate inflate;
	int i = msr->param.nInflateReps + 1;
	printf("Inflating number of particles by a factor of %d\n",i*i*i);
	inflate.nInflateReps = msr->param.nInflateReps;
	pstInflate(msr->pst,&inflate,sizeof(inflate),NULL,NULL);
	msr->N *= i*i*i;
	msr->nGas *= i*i*i;
	msr->nDark *= i*i*i;
	msr->nStar *= i*i*i;
	msr->nMaxOrder = msr->N;
	}
    }

void msrDomainDecomp(MSR msr,int iRung,int bOthers,int bSplitVA) {
    msrDomainDecompOld(msr,iRung,bSplitVA);
    }

/*
** This the meat of the tree build, but will be called by differently named
** functions in order to implement special features without recoding...
*/
static void BuildTree(MSR msr,int bNeedEwald,uint32_t uRoot,uint32_t utRoot) {
    struct inBuildTree in;
    struct inCalcRoot calc;
    struct outCalcRoot root;
    struct inDistribTopTree *pDistribTop;
    int i;
    PST pst0;
    LCL *plcl;
    PKD pkd;
    KDN *pkdn;
    int iDum;
    int nTopTree;
    double sec,dsec;

    pst0 = msr->pst;
    while (pst0->nLeaves > 1)
	pst0 = pst0->pstLower;
    plcl = pst0->plcl;
    pkd = plcl->pkd;

    nTopTree = pkdNodeSize(pkd) * (2*msr->nThreads-1);
    pDistribTop = malloc( sizeof(struct inDistribTopTree) + nTopTree );
    assert(pDistribTop != NULL);
    pDistribTop->uRoot = uRoot;
    pkdn = (KDN *)(pDistribTop + 1);

    in.nBucket = msr->param.nBucket;
    in.nGroup = msr->param.nGroup;
    in.uRoot = uRoot;
    in.utRoot = utRoot;
    sec = msrTime();
    pstBuildTree(msr->pst,&in,sizeof(in),pkdn,&nTopTree);
    pDistribTop->nTop = nTopTree / pkdNodeSize(pkd);
    assert(pDistribTop->nTop == (2*msr->nThreads-1));
    pstDistribTopTree(msr->pst,pDistribTop,sizeof(struct inDistribTopTree) + nTopTree,NULL,NULL);
    dsec = msrTime() - sec;
    printf("Tree built, Wallclock: %f secs\n\n",dsec);

    if (bNeedEwald) {
	/*
	** For simplicity we will skip calculating the Root for all particles
	** with exclude very active since there are missing particles which
	** could add to the mass and because it probably is not important to
	** update the root so frequently.
	*/
	double kdn_r[3];
	pkdNodeGetPos(pkd,pkdn,kdn_r);
	calc.com[0] = kdn_r[0];
	calc.com[1] = kdn_r[1];
	calc.com[2] = kdn_r[2];
	calc.uRoot = uRoot;
	pstCalcRoot(msr->pst,&calc,sizeof(calc),&root,&iDum);
	msr->momTreeRoot[uRoot] = root.momc;
	msr->momTreeCom[uRoot][0] = kdn_r[0];
	msr->momTreeCom[uRoot][1] = kdn_r[1];
	msr->momTreeCom[uRoot][2] = kdn_r[2];
	}

    free(pDistribTop);
    }

void msrBuildTree(MSR msr,double dTime,int bNeedEwald) {
    msrprintf(msr,"Building local trees...\n\n");

    struct inDumpTrees dump;
    dump.bOnlyVA = 0;
    dump.uRungDD = IRUNGMAX;
    pstDumpTrees(msr->pst,&dump,sizeof(dump),NULL,NULL);
    BuildTree(msr,bNeedEwald,ROOT,0);

    if (bNeedEwald) {
	struct ioDistribRoot droot;
	droot.momc = msr->momTreeRoot[ROOT];
	droot.r[0] = msr->momTreeCom[ROOT][0];
	droot.r[1] = msr->momTreeCom[ROOT][1];
	droot.r[2] = msr->momTreeCom[ROOT][2];
	pstDistribRoot(msr->pst,&droot,sizeof(struct ioDistribRoot),NULL,NULL);
	}
    }

/*
** Separates the particles into two trees, and builds the "fixed" tree.
*/
void msrBuildTreeFixed(MSR msr,double dTime,int bNeedEwald,uint8_t uRungDD) {
    msrprintf(msr,"Building fixed local trees...\n\n");
    BuildTree(msr,bNeedEwald,FIXROOT,0);
    }

void msrBuildTreeActive(MSR msr,double dTime,int bNeedEwald,uint8_t uRungDD) {
   /*
    ** The trees reset/removed. This does the following:
    **   1. Closes any open cell cache (it will be subsequently invalid)
    **   2. Resets the number of used nodes to zero (or more if we keep the main tree)
    **   3. Sets up the ROOT and FIXROOT node (either of which may have zero particles).
    */

    msrprintf(msr,"Building active local trees...\n\n");

    struct inDumpTrees dump;
    dump.bOnlyVA = 1;
    dump.uRungDD = uRungDD;
    pstDumpTrees(msr->pst,&dump,sizeof(dump),NULL,NULL);

    /* New build the very active tree */
    BuildTree(msr,bNeedEwald,ROOT,FIXROOT);

    /* For ewald we have to shift and combine the individual tree moments */
    if (bNeedEwald) {
	struct ioDistribRoot droot;
	MOMC momc;
	double *com1 = msr->momTreeCom[FIXROOT];
	double    m1 = msr->momTreeRoot[FIXROOT].m;
	double *com2 = msr->momTreeCom[ROOT];
	double    m2 = msr->momTreeRoot[ROOT].m;
	double ifMass = 1.0 / (m1 + m2);
	double x, y, z;
	int j;

	/* New Center of Mass, then shift and scale the moments */
	for (j=0;j<3;++j) droot.r[j] = ifMass*(m1*com1[j] + m2*com2[j]);

	droot.momc = msr->momTreeRoot[FIXROOT];
	x = com1[0] - droot.r[0];
	y = com1[1] - droot.r[1];
	z = com1[2] - droot.r[2];
	momShiftMomc(&droot.momc,x,y,z);

	momc = msr->momTreeRoot[ROOT];
	x = com2[0] - droot.r[0];
	y = com2[1] - droot.r[1];
	z = com2[2] - droot.r[2];
	momShiftMomc(&momc,x,y,z);

	momAddMomc(&droot.momc, &momc);

	pstDistribRoot(msr->pst,&droot,sizeof(struct ioDistribRoot),NULL,NULL);
	}
    }

void msrBuildTreeByRung(MSR msr,double dTime,int bNeedEwald,int iRung) {
    assert(0);
//    TREESPEC spec[2];
//    int nTrees = 1;
//    spec[0].uRoot = ROOT;
//    spec[0].uRungFirst = iRung;
//    spec[0].uRungLast = iRung;
//    if (msrCurrMaxRung(msr) > iRung) {
//	spec[1].uRoot = ROOT+1;
//	spec[1].uRungFirst = iRung+1;
//	spec[1].uRungLast = MAX_RUNG;
//	++nTrees;
//	}
//    BuildTree(msr,bNeedEwald,nTrees,spec);
    }

void msrBuildTreeExcludeVeryActive(MSR msr,double dTime) {
    assert(0);
//    TREESPEC spec[2];
//    spec[0].uRoot = ROOT;
//    spec[0].uRungFirst = 0;
//    spec[0].uRungLast = msr->iRungVeryActive;
//    spec[1].uRoot = ROOT+1;
//    spec[1].uRungFirst = spec[0].uRungLast+1;
//    spec[1].uRungLast = MAX_RUNG;
//    BuildTree(msr,0,2,spec);
    }

void msrBuildTreeMarked(MSR msr,double dTime) {
    pstTreeInitMarked(msr->pst,NULL,0,NULL,NULL);
    BuildTree(msr,0,ROOT,0);
    }

void msrReorder(MSR msr) { //IA TODO FIXME figure out why this crashes and creates pathological configurations (probably sph field not being copied...)
//  if (!msr->param.bMemUnordered) {
//    struct inDomainOrder in;
//    double sec,dsec;

//    msrprintf(msr,"Ordering...\n");
//    sec = msrTime();
//    in.iMinOrder = 0;
//    in.iMaxOrder = msrMaxOrder(msr)-1;
//    pstDomainOrder(msr->pst,&in,sizeof(in),NULL,NULL);
//    in.iMinOrder = 0;
//    in.iMaxOrder = msrMaxOrder(msr)-1;
//    pstLocalOrder(msr->pst,&in,sizeof(in),NULL,NULL);
//    dsec = msrTime() - sec;
//    msrprintf(msr,"Order established, Wallclock: %f secs\n\n",dsec);

//  /*
//  ** Mark domain decomp as not done.
//  */
//    msr->iLastRungRT = -1;
//    msr->iLastRungDD = -1;
//    }
  }

void msrOutASCII(MSR msr,const char *pszFile,int iType,int nDims) {

    char achOutFile[PST_FILENAME_SIZE];
    LCL *plcl;
    PST pst0;
    int id,iDim;
    int inswap;
    PKDOUT pkdout;
    const char *arrayOrVector;
    struct outSetTotal total;
    int rID;

    switch(nDims) {
    case 1: arrayOrVector = "vector"; break; /* JW -- seems like a bug */
    case 3: arrayOrVector = "array";  break;
    default:arrayOrVector = NULL;
	    assert(nDims==1 || nDims==3);
	}

    pst0 = msr->pst;
    while (pst0->nLeaves > 1)
	pst0 = pst0->pstLower;
    plcl = pst0->plcl;

    pstSetTotal(msr->pst,NULL,0,&total,NULL);

    if (pszFile) {
	/*
	** Add Data Subpath for local and non-local names.
	*/
	_msrMakePath(msr->param.achDataSubPath,pszFile,achOutFile);

	switch(msr->param.iCompress) {
#ifdef HAVE_LIBBZ2
	case PKDOUT_TYPE_BZIP2:
	    strcat(achOutFile,".bz2");
	    break;
#endif
#ifdef HAVE_LIBZ
	case PKDOUT_TYPE_ZLIB:
	    strcat(achOutFile,".gz");
	    break;
#endif
	default:
	    break;
	    }

	msrprintf(msr, "Writing %s to %s\n", arrayOrVector, achOutFile );
	}
    else {
	printf("No %s Output File specified\n", arrayOrVector);
	_msrExit(msr,1);
	return;
	}

    if (msr->param.bParaWrite && msr->param.iCompress) {
	struct inCompressASCII in;
	struct outCompressASCII out;
	struct inWriteASCII inWrite;
	int nOut;
	FILE *fp;

	fp = fopen(achOutFile,"wb");
	if ( fp==NULL) {
	    printf("Could not create %s Output File:%s\n",arrayOrVector, achOutFile);
	    _msrExit(msr,1);
	    }
	fclose(fp);

	inWrite.nFileOffset = 0;
	for( iDim=0; iDim<nDims; iDim++ ) {
	    in.nTotal = total.nTotal;
	    in.iFile = msr->param.iCompress;
	    in.iType = iType;
	    in.iDim = iDim;
	    pstCompressASCII(msr->pst,&in,sizeof(in),&out,&nOut);
	    strcpy(inWrite.achOutFile,achOutFile);
	    pstWriteASCII(msr->pst,&inWrite,sizeof(inWrite),NULL,NULL);
	    inWrite.nFileOffset += out.nBytes;
	    }
	}
    else {
	pkdout = pkdOpenOutASCII(plcl->pkd,achOutFile,"wb",msr->param.iCompress,iType);
	if (!pkdout) {
	    printf("Could not open %s Output File:%s\n",arrayOrVector,achOutFile);
	    _msrExit(msr,1);
	    }

	pkdOutHdr(plcl->pkd,pkdout,total.nTotal);

	/*
	 * First write our own particles.
	 */
	for (iDim=0;iDim<nDims;++iDim) {
	    pkdOutASCII(plcl->pkd,pkdout,iType,iDim);
	    for (id=1;id<msr->nThreads;++id) {
		/*
		 * Swap particles with the remote processor.
		 */
		inswap = 0;
		rID = mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
		pkdSwapAll(plcl->pkd, id);
		mdlGetReply(pst0->mdl,rID,NULL,NULL);
		/*
		 * Write the swapped particles.
		 */
		pkdOutASCII(plcl->pkd,pkdout,iType,iDim);
		/*
		 * Swap them back again.
		 */
		inswap = 0;
		rID = mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
		pkdSwapAll(plcl->pkd,id);
		mdlGetReply(pst0->mdl,rID,NULL,NULL);
		}
	    }
	pkdCloseOutASCII(plcl->pkd,pkdout);
	}
    }

void msrOutArray(MSR msr,const char *pszFile,int iType) {
    msrOutASCII(msr,pszFile,iType,1);
    }

void msrOutVector(MSR msr,const char *pszFile,int iType) {
    msrOutASCII(msr,pszFile,iType,3);
    }

void msrSmoothSetSMF(MSR msr, SMF *smf, double dTime) {
#ifdef SYMBA
    smf->dSunMass = msr->dSunMass;
#endif
    smf->dTime = dTime;
    if (msrComove(msr)) {
	smf->bComove = 1;
	smf->H = csmTime2Hub(msr->param.csm,dTime);
	smf->a = csmTime2Exp(msr->param.csm,dTime);
	}
    else {
	smf->bComove = 0;
	smf->H = 0.0;
	smf->a = 1.0;
	}
    smf->alpha = msr->param.dConstAlpha;
    smf->beta = msr->param.dConstBeta;
    smf->gamma = msr->param.dConstGamma;
    smf->dDelta = msr->param.dDelta;
    smf->dEtaCourant = msr->param.dEtaCourant;
    smf->iViscosityLimiter = msr->param.iViscosityLimiter;
    smf->iDiffusion = msr->param.iDiffusion;
    smf->dMetalDiffusionCoeff = msr->param.dMetalDiffusionCoeff;
    smf->dThermalDiffusionCoeff = msr->param.dThermalDiffusionCoeff;
    /* For SF & FB in code units */
#define SECONDSPERYEAR   31557600.
    if (msr->param.bGasIsothermal) smf->SFdESNPerStarMass = 0;
    else smf->SFdESNPerStarMass = msr->param.SFdESNPerStarMass/msr->param.dErgPerGmUnit;
    smf->SFdtCoolingShutoff = msr->param.SFdtCoolingShutoff*SECONDSPERYEAR/msr->param.dSecUnit;
    /* avoid alignment in steps + feedback time -- increase FB time by small tweak */
    smf->SFdtFeedbackDelay = msr->param.SFdtFeedbackDelay*1.0000013254678*SECONDSPERYEAR/msr->param.dSecUnit;
    smf->SFdMassLossPerStarMass = msr->param.SFdMassLossPerStarMass;
    smf->SFdZMassPerStarMass = msr->param.SFdZMassPerStarMass;
    smf->SFdFBFac = 0.5/((1+0.6*smf->alpha)/(smf->a*smf->dEtaCourant))
	/(msr->param.SFdvFB/msr->param.dKmPerSecUnit);
    smf->FirstHydroLoop = msrFirstHydroLoop(msr);
    smf->dCFLacc = msr->param.dCFLacc;
    smf->dNeighborsStd = msrNeighborsStd(msr);
    }

void msrSmooth(MSR msr,double dTime,int iSmoothType,int bSymmetric,int nSmooth) {
    struct inSmooth in;

    in.nSmooth = nSmooth;
    in.bPeriodic = msr->param.bPeriodic;
    in.bSymmetric = bSymmetric;
    in.iSmoothType = iSmoothType;
    msrSmoothSetSMF(msr, &(in.smf), dTime);
    if (msr->param.bVStep) {
	double sec,dsec;
	printf("Smoothing...\n");
	sec = msrTime();
	pstSmooth(msr->pst,&in,sizeof(in),NULL,NULL);
	dsec = msrTime() - sec;
	printf("Smooth Calculated, Wallclock: %f secs\n\n",dsec);
	}
    else {
	pstSmooth(msr->pst,&in,sizeof(in),NULL,NULL);
	}
    }


#ifdef FAST_GAS
void msrFastGasPhase1(MSR msr,double dTime,int iSmoothType) {
    struct inSmooth in;

    in.nSmooth = msr->param.nSmooth;
    in.bPeriodic = msr->param.bPeriodic;
    in.bSymmetric = 0;
    in.iSmoothType = iSmoothType;
    msrSmoothSetSMF(msr, &(in.smf), dTime);
    if (msr->param.bVStep) {
	double sec,dsec;
	printf("FastGas Phase 1 Smoothing...\n");
	sec = msrTime();
	pstFastGasPhase1(msr->pst,&in,sizeof(in),NULL,NULL);
	dsec = msrTime() - sec;
	printf("FastGas Phase 1 Smooth Calculated, Wallclock: %f secs\n\n",dsec);
	}
    else {
	pstFastGasPhase1(msr->pst,&in,sizeof(in),NULL,NULL);
	}
    }


void msrFastGasPhase2(MSR msr,double dTime,int iSmoothType) {
    struct inSmooth in;

    in.nSmooth = msr->param.nSmooth;
    in.bPeriodic = msr->param.bPeriodic;
    in.bSymmetric = 0;
    in.iSmoothType = iSmoothType;
    msrSmoothSetSMF(msr, &(in.smf), dTime);
    if (msr->param.bVStep) {
	double sec,dsec;
	printf("FastGas Phase 2 Smoothing...\n");
	sec = msrTime();
	pstFastGasPhase2(msr->pst,&in,sizeof(in),NULL,NULL);
	dsec = msrTime() - sec;
	printf("FastGas Phase 2 Smooth Calculated, Wallclock: %f secs\n\n",dsec);
	}
    else {
	pstFastGasPhase2(msr->pst,&in,sizeof(in),NULL,NULL);
	}
    }
#endif

int msrReSmooth(MSR msr,double dTime,int iSmoothType,int bSymmetric, int bFirstStep) {
    struct inSmooth in;
    struct outSmooth out;
    int nOut;

    in.nSmooth = msr->param.nSmooth;
    in.bPeriodic = msr->param.bPeriodic;
    in.bSymmetric = bSymmetric;
    in.iSmoothType = iSmoothType;
    msrSmoothSetSMF(msr, &(in.smf), dTime);
    if (bFirstStep) in.smf.dDelta = 0.0; // Avoid adding fluxes and doing the spatial extrapolation

    if (msr->param.bVStep) {
	double sec,dsec;
	//printf("ReSmoothing...\n");
	sec = msrTime();
	pstReSmooth(msr->pst,&in,sizeof(in),&out,&nOut);
	dsec = msrTime() - sec;
	//printf("ReSmooth Calculated on %d particles, Wallclock: %f secs\n\n", out.nSmoothed, dsec);
	}
    else {
	pstReSmooth(msr->pst,&in,sizeof(in),&out,&nOut);
	}
    return out.nSmoothed;
    }

void msrUpdateSoft(MSR msr,double dTime) {
    if (!(msr->param.bPhysicalSoft)) return;
    if (msr->param.bPhysicalSoft) {
	struct inPhysicalSoft in;

	in.dFac = 1./csmTime2Exp(msr->param.csm,dTime);
	in.bSoftMaxMul = msr->param.bSoftMaxMul;
	in.dSoftMax = msr->param.dSoftMax;

	if (msr->param.bSoftMaxMul && in.dFac > in.dSoftMax) in.dFac = in.dSoftMax;

	pstPhysicalSoft(msr->pst,&in,sizeof(in),NULL,NULL);
	}
    }

#define PRINTGRID(w,FRM,VAR) {						\
    printf("      % *d % *d % *d % *d % *d % *d % *d % *d % *d % *d\n",\
	   w,0,w,1,w,2,w,3,w,4,w,5,w,6,w,7,w,8,w,9);		       \
    for (i=0;i<msr->nThreads/10;++i) {\
	printf("%4d: "FRM" "FRM" "FRM" "FRM" "FRM" "FRM" "FRM" "FRM" "FRM" "FRM"\n",i*10,\
	       out[i*10+0].VAR,out[i*10+1].VAR,out[i*10+2].VAR,out[i*10+3].VAR,out[i*10+4].VAR,\
	       out[i*10+5].VAR,out[i*10+6].VAR,out[i*10+7].VAR,out[i*10+8].VAR,out[i*10+9].VAR);\
	}\
    switch (msr->nThreads%10) {\
    case 0: break;\
    case 1: printf("%4d: "FRM"\n",i*10,\
		   out[i*10+0].VAR); break;\
    case 2: printf("%4d: "FRM" "FRM"\n",i*10,\
		   out[i*10+0].VAR,out[i*10+1].VAR); break;\
    case 3: printf("%4d: "FRM" "FRM" "FRM"\n",i*10,\
		   out[i*10+0].VAR,out[i*10+1].VAR,out[i*10+2].VAR); break;\
    case 4: printf("%4d: "FRM" "FRM" "FRM" "FRM"\n",i*10,\
		   out[i*10+0].VAR,out[i*10+1].VAR,out[i*10+2].VAR,out[i*10+3].VAR); break;\
    case 5: printf("%4d: "FRM" "FRM" "FRM" "FRM" "FRM"\n",i*10,\
		   out[i*10+0].VAR,out[i*10+1].VAR,out[i*10+2].VAR,out[i*10+3].VAR,out[i*10+4].VAR); break;\
    case 6: printf("%4d: "FRM" "FRM" "FRM" "FRM" "FRM" "FRM"\n",i*10,\
		   out[i*10+0].VAR,out[i*10+1].VAR,out[i*10+2].VAR,out[i*10+3].VAR,out[i*10+4].VAR,\
		   out[i*10+5].VAR); break;\
    case 7: printf("%4d: "FRM" "FRM" "FRM" "FRM" "FRM" "FRM" "FRM"\n",i*10,\
		   out[i*10+0].VAR,out[i*10+1].VAR,out[i*10+2].VAR,out[i*10+3].VAR,out[i*10+4].VAR,\
		   out[i*10+5].VAR,out[i*10+6].VAR); break;\
    case 8: printf("%4d: "FRM" "FRM" "FRM" "FRM" "FRM" "FRM" "FRM" "FRM"\n",i*10,\
		   out[i*10+0].VAR,out[i*10+1].VAR,out[i*10+2].VAR,out[i*10+3].VAR,out[i*10+4].VAR,\
		   out[i*10+5].VAR,out[i*10+6].VAR,out[i*10+7].VAR); break;\
    case 9: printf("%4d: "FRM" "FRM" "FRM" "FRM" "FRM" "FRM" "FRM" "FRM" "FRM"\n",i*10,\
		   out[i*10+0].VAR,out[i*10+1].VAR,out[i*10+2].VAR,out[i*10+3].VAR,out[i*10+4].VAR,\
		   out[i*10+5].VAR,out[i*10+6].VAR,out[i*10+7].VAR,out[i*10+8].VAR); break;\
    }\
}

void msrHostname(MSR msr) {
    struct outHostname *out;
    int i,iDum;
    out = malloc(msr->nThreads*sizeof(struct outHostname));
    assert(out != NULL);
    pstHostname(msr->pst,0,0,out,&iDum);
    printf("Host Names:\n");
    PRINTGRID(12,"%12.12s",szHostname);
    printf("MPI Rank:\n");
    PRINTGRID(8,"% 8d",iMpiID);
    free(out);
    }

void msrMemStatus(MSR msr) {
    struct outMemStatus *out;
    int i,iDum;
    if (msr->param.bVDetails) {
	out = malloc(msr->nThreads*sizeof(struct outMemStatus));
	assert(out != NULL);
	pstMemStatus(msr->pst,0,0,out,&iDum);
#ifdef __linux__
	printf("Resident (MB):\n");
	PRINTGRID(8,"%8"PRIu64,rss);
	printf("Free Memory (MB):\n");
	PRINTGRID(8,"%8"PRIu64,freeMemory);
#endif
	printf("Tree size (MB):\n");
	PRINTGRID(8,"%8"PRIu64,nBytesTree/1024/1024);
	printf("Checklist size (KB):\n");
	PRINTGRID(8,"%8"PRIu64,nBytesCl/1024);
	printf("Particle List size (KB):\n");
	PRINTGRID(8,"%8"PRIu64,nBytesIlp/1024);
	printf("Cell List size (KB):\n");
	PRINTGRID(8,"%8"PRIu64,nBytesIlc/1024);
	free(out);
	}
    }

void msrPrintStat(STAT *ps,char *pszPrefix,int p) {
    double dSum = ps->dSum;
    double dMax = ps->dMax;
    const char *minmax = "max";

    if (dSum<0) {
	dSum = -dSum;
	dMax = -dMax;
	minmax = "min";
	}
    if (ps->n > 1) {
	printf("%s %s=%8.*f @%5d avg=%8.*f of %5d std-dev=%8.*f\n",pszPrefix,minmax,
	    p,dMax,ps->idMax,p,dSum/ps->n,ps->n,p,sqrt((ps->dSum2 - ps->dSum*ps->dSum/ps->n)/(ps->n-1)));
	}
    else if (ps->n == 1) {
	printf("%s %s=%8.*f @%5d\n",pszPrefix,minmax,p,dMax,ps->idMax);
	}
    else {
	printf("%s no data\n",pszPrefix);
	}	
    }


void msrLightCone(MSR msr,double dTime,uint8_t uRungLo,uint8_t uRungHi) {
    struct inLightCone in;
    double sec,dsec,dt;
    double dTimeLCP;
    int i;

    if (!msr->param.bLightCone) return;
    in.dLookbackFac = csmComoveKickFac(msr->param.csm,dTime,(csmExp2Time(msr->param.csm,1.0) - dTime));
    dTimeLCP = csmExp2Time(msr->param.csm,1.0/(1.0+msr->param.dRedshiftLCP));
    in.dLookbackFacLCP = csmComoveKickFac(msr->param.csm,dTimeLCP,(csmExp2Time(msr->param.csm,1.0) - dTimeLCP));
    in.uRungLo = uRungLo;
    in.uRungHi = uRungHi;
    for (i=0,dt=msr->param.dDelta;i<=msr->param.iMaxRung;++i,dt*=0.5) {
	in.dtLCDrift[i] = 0.0;
	in.dtLCKick[i] = 0.0;
	if (i>=uRungLo) {
	    if (msr->param.csm->val.bComove) {
		in.dtLCDrift[i] = csmComoveDriftFac(msr->param.csm,dTime,dt);
		in.dtLCKick[i] = csmComoveKickFac(msr->param.csm,dTime,dt);
		}
	    else {
	        in.dtLCDrift[i] = dt;
	        in.dtLCKick[i] = dt;
		}
	    }
	}
    sec = msrTime();
    pstLightCone(msr->pst,&in,sizeof(in),NULL,NULL);
    dsec = msrTime() - sec;
    if (msr->param.bVStep) {
	/*
	** Output some info...
	*/
	printf("Light-Cone Calculated, Wallclock: %f secs\n",dsec);
	}
    }


uint8_t msrGravity(MSR msr,uint8_t uRungLo, uint8_t uRungHi,int iRoot1,int iRoot2,
    double dTime, double dStep,int bKickClose,int bKickOpen,int bEwald,int nGroup,int *piSec,uint64_t *pnActive) {
    struct inGravity in;
    struct outGravityPerProc *out;
    struct outGravityPerProc *outend;
    struct outGravityReduct *outr;
    uint64_t nRungSum[IRUNGMAX+1];
    int i,id,iDum;
    double sec,dsec,dTotFlop,dt,a;
    double dTimeLCP;
    uint8_t uRungMax=0;
    char c;

    if (msr->param.bVStep) printf("Calculating Gravity, Step:%f (rung %d)\n",dStep,uRungLo);
    in.dTime = dTime;
    in.nReps = msr->param.nReplicas;
    in.bPeriodic = msr->param.bPeriodic;
    in.bEwald = bEwald;
    in.nGroup = nGroup;
    in.dEwCut = msr->param.dEwCut;
    in.dEwhCut = msr->param.dEwhCut;
    in.uRungLo = uRungLo;
    in.uRungHi = uRungHi;
    in.dThetaMin = msr->dThetaMin;
    in.iRoot1 = iRoot1;
    in.iRoot2 = iRoot2;
    if (msr->param.bLightCone) {
	in.dLookbackFac = csmComoveKickFac(msr->param.csm,dTime,(csmExp2Time(msr->param.csm,1.0) - dTime));
	dTimeLCP = csmExp2Time(msr->param.csm,1.0/(1.0+msr->param.dRedshiftLCP));
	in.dLookbackFacLCP = csmComoveKickFac(msr->param.csm,dTimeLCP,(csmExp2Time(msr->param.csm,1.0) - dTimeLCP));
	}
    else {
	in.dLookbackFac = 0.0;
	in.dLookbackFacLCP = 0.0;
	}
    /*
    ** Now calculate the timestepping factors for kick close and open if the
    ** gravity should kick the particles. If the code uses bKickClose and 
    ** bKickOpen it no longer needs to store accelerations per particle.
    */
    in.bKickClose = bKickClose;
    in.bKickOpen = bKickOpen;
    if (msr->param.csm->val.bComove) {
	a = csmTime2Exp(msr->param.csm,dTime);
	in.dAccFac = 1.0/(a*a*a);
#if 0
	// erf2: in.dThetaMin = 0.4 + 0.3*erf(a*10.0);
	in.dThetaMin = 0.35 + 0.35*erf(a*8.0); // erf3
	if ( !prmSpecified(msr->prm,"nReplicas") && msr->param.nReplicas>=1 ) {
	    in.nReps = in.dThetaMin < 0.52 ? 2 : 1;
	    }
	printf("z=%.1f theta=%.3f nReps=%d\n", 1.0/a - 1.0, in.dThetaMin,in.nReps);
#endif
	}
    else {
	in.dAccFac = 1.0;
	}
    for (i=0,dt=0.5*msr->param.dDelta;i<=msr->param.iMaxRung;++i,dt*=0.5) {
	in.dtClose[i] = 0.0;
	in.dtOpen[i] = 0.0;
	if (i>=uRungLo) {
	    if (msr->param.csm->val.bComove) {
		if (bKickClose) {
		    in.dtClose[i] = csmComoveKickFac(msr->param.csm,dTime-dt,dt);
		    }
		if (bKickOpen) {
		    in.dtOpen[i] = csmComoveKickFac(msr->param.csm,dTime,dt);
		    }
		}
	    else {
		if (bKickClose) in.dtClose[i] = dt;
		if (bKickOpen) in.dtOpen[i] = dt;
		}
	    }
	}
    /*
    ** Note that in this loop we initialize dt with the full step, not a half step!
    */
    for (i=0,dt=msr->param.dDelta;i<=msr->param.iMaxRung;++i,dt*=0.5) {
	in.dtLCDrift[i] = 0.0;
	in.dtLCKick[i] = 0.0;
	if (i>=uRungLo) {
	    if (msr->param.csm->val.bComove) {
		in.dtLCDrift[i] = csmComoveDriftFac(msr->param.csm,dTime,dt);
		in.dtLCKick[i] = csmComoveKickFac(msr->param.csm,dTime,dt);
		}
	    else {
	        in.dtLCDrift[i] = dt;
	        in.dtLCKick[i] = dt;
		}
	    }
	}
    in.bLinearSpecies = (strlen(msr->param.csm->val.classData.achLinSpecies) > 0);
    out = malloc(msr->nThreads*sizeof(struct outGravityPerProc) + sizeof(struct outGravityReduct));
    assert(out != NULL);
    outend = out + msr->nThreads;
    outr = (struct outGravityReduct *)outend;

    sec = msrTime();
    pstGravity(msr->pst,&in,sizeof(in),out,&iDum);
    dsec = msrTime() - sec;

    *piSec = d2i(dsec);
    *pnActive = outr->nActive;

    if (bKickOpen) {
	for (i=IRUNGMAX;i>=uRungLo;--i) {
	    if (outr->nRung[i]) break;
	    }
	assert(i >= uRungLo);
	uRungMax = i;
	msr->iCurrMaxRung = uRungMax;   /* this assignment shouldn't be needed */
	/*
	** Update only the active rung counts in the master rung counts.
	** We need to go all the way to IRUNGMAX to clear any prior counts at rungs 
	** deeper than the current uRungMax!
	*/
	for (i=uRungLo;i<=IRUNGMAX;++i) msr->nRung[i] = outr->nRung[i];

	const uint64_t nDT = d2u64(msr->N*msr->param.dFracDualTree);
	const uint64_t nDD = d2u64(msr->N*msr->param.dFracNoDomainDecomp);
	uint64_t nActive = 0;
	msr->iRungDD = 0;
	msr->iRungDT = 0;
	for (i=msr->iCurrMaxRung;i>=0;--i) {
	    nActive += msr->nRung[i];
	    if (nActive > nDT && !msr->iRungDT) msr->iRungDT = i;
	    if (nActive > nDD && !msr->iRungDD) msr->iRungDD = i;
	    }
	}
    if (msr->param.bVStep) {
	/*
	** Output some info...
	*/
	dTotFlop = outr->sFlop.dSum;
	if (dsec > 0.0) {
	    double dGFlops = dTotFlop/dsec;
	    printf("Gravity Calculated, Wallclock: %f secs, Gflops:%.1f, Total Gflop:%.3g\n",
		   dsec,dGFlops,dTotFlop);
	    printf("  Gflops: CPU:%.1f,%.1f GPU:%.1f,%.1f",
		outr->dFlopSingleCPU/dsec,outr->dFlopDoubleCPU/dsec,
		outr->dFlopSingleGPU/dsec,outr->dFlopDoubleGPU/dsec);
	    }
	else {
	    printf("Gravity Calculated, Wallclock: %f secs, Gflops:unknown, Total Gflop:%.3g\n",
		   dsec,dTotFlop);
	    }
	printf("  Gflop: CPU:%.3g,%.3g GPU:%.3g,%.3g\n",
	    outr->dFlopSingleCPU,outr->dFlopDoubleCPU,
	    outr->dFlopSingleGPU,outr->dFlopDoubleGPU);
	msrPrintStat(&outr->sLocal,         "  particle  load:",0);
	msrPrintStat(&outr->sActive,        "  actives   load:",0);
	msrPrintStat(&outr->sFlop,          "  Gflop     load:",1);
	msrPrintStat(&outr->sPart,          "  P-P per active:",2);
	msrPrintStat(&outr->sCell,          "  P-C per active:",2);
#ifdef INSTRUMENT
	msrPrintStat(&outr->sComputing,     "     % computing:",3);
	msrPrintStat(&outr->sWaiting,       "     %   waiting:",3);
	msrPrintStat(&outr->sSynchronizing, "     %   syncing:",3);
#endif
#ifdef __linux__
	msrPrintStat(&outr->sFreeMemory, "free memory (GB):", 3);
	msrPrintStat(&outr->sRSS,           "   resident size:",3);
#endif
	printf("  (cache access statistics are given per active particle)\n");
	msrPrintStat(&outr->sPartNumAccess, "  P-cache access:",1);
	msrPrintStat(&outr->sCellNumAccess, "  C-cache access:",1);
	msrPrintStat(&outr->sPartMissRatio, "  P-cache miss %:",2);
	msrPrintStat(&outr->sCellMissRatio, "  C-cache miss %:",2);
	/*
	** Now comes the really verbose output for each processor.
	*/
	if (msr->param.bVDetails) {
	    printf("Walk Timings:\n");
	    PRINTGRID(8,"% 8.2f",dWalkTime);
	    }
	}
    if (msr->param.bVRungStat && bKickOpen) {
	printf("Rung distribution:\n");
	printf("\n");
	nRungSum[uRungMax] = msr->nRung[uRungMax];
	for (i=uRungMax-1;i>=0;--i) {
	    nRungSum[i] = nRungSum[i+1] + msr->nRung[i];
	    }
	for (i=0;i<=uRungMax;++i) {
	    if (msr->nRung[i]) break;
	    }
	if (nRungSum[0]>0) for (;i<=uRungMax;++i) {
	    c = ' ';
	    printf(" %c rung:%d %14"PRIu64"    %14"PRIu64"  %3.0f %%\n",
		c,i,msr->nRung[i],nRungSum[i],
		ceil(100.0 * nRungSum[i] / nRungSum[0]));
	    }
	printf("\n");
	}
    free(out);
    return(uRungMax);
    }


void msrCalcEandL(MSR msr,int bFirst,double dTime,double *E,double *T,double *U,double *Eth,double *L,double *F,double *W) {
    struct outCalcEandL out;
    double a;
    int k;

    pstCalcEandL(msr->pst,NULL,0,&out,NULL);
    *T = out.T;
    *U = out.U;
    *Eth = out.Eth;
    for (k=0;k<3;k++) L[k] = out.L[k];
    for (k=0;k<3;k++) F[k] = out.F[k];
    *W = out.W;
    /*
    ** Do the comoving coordinates stuff.
    ** Currently L is not adjusted for this. Should it be?
    */
    a = csmTime2Exp(msr->param.csm,dTime);
    if (!msr->param.csm->val.bComove) *T *= pow(a,4.0);
    /*
     * Estimate integral (\dot a*U*dt) over the interval.
     * Note that this is equal to integral (W*da) and the latter
     * is more accurate when a is changing rapidly.
     */
    if (msr->param.csm->val.bComove && !bFirst) {
	msr->dEcosmo += 0.5*(a - csmTime2Exp(msr->param.csm, msr->dTimeOld))
			*((*U) + msr->dUOld);
	}
    else {
	msr->dEcosmo = 0.0;
	}
    msr->dTimeOld = dTime;
    msr->dUOld = *U;
    *U *= a;
    *E = (*T) + (*U) - msr->dEcosmo + a*a*(*Eth);
    }


void msrDrift(MSR msr,double dTime,double dDelta,int iRoot) {
    struct inDrift in;

    if (msr->param.csm->val.bComove) {
	in.dDelta = csmComoveDriftFac(msr->param.csm,dTime,dDelta);
	in.dDeltaVPred = csmComoveKickFac(msr->param.csm,dTime,dDelta);
	}
    else {
	in.dDelta = dDelta;
	in.dDeltaVPred = dDelta;
	}
    in.dTime = dTime;
    in.dDeltaUPred = dDelta;
    in.iRoot = iRoot;
    pstDrift(msr->pst,&in,sizeof(in),NULL,NULL);
    }
/* IA
 * We add the gravitational work following a second-order scheme (e.g., equations H1 and H2 of Hopkins 2015)
 * To that end, we need the gradient of the gravitational potential at the begining and end of the timestep.
 * This gradient is just the acceleration due to gravity.
 *
 */
void  msrApplyGravWork(MSR msr,double dTime,double dDelta,uint8_t uRungLo,uint8_t uRungHi) {
    struct inKick in;


    in.dTime = dTime;
    if (msr->param.csm->val.bComove) {
	in.dDelta = csmComoveKickFac(msr->param.csm,dTime,dDelta);
	in.dDeltaVPred = 0;
    }
    else {
	in.dDelta = dDelta;
	in.dDeltaVPred = 0;
    }
    in.dDeltaU = dDelta;
    in.dDeltaUPred = 0;
    in.uRungLo = uRungLo;
    in.uRungHi = uRungHi;

    pstApplyGravWork(msr->pst,&in,sizeof(in),NULL,NULL);
    }



void msrUpdateConsVars(MSR msr,double dTime,double dDelta,int iRoot) {
    struct inDrift in;

    assert(iRoot!=0); //IA: Dual tree not implemented

    printf("(msrUpdateConsVars) Begin \n");
    if (msr->param.csm->val.bComove) {
	in.dDelta = csmComoveDriftFac(msr->param.csm,dTime,dDelta);
	in.dDeltaVPred = csmComoveKickFac(msr->param.csm,dTime,dDelta);
	}
    else {
	in.dDelta = dDelta;
	in.dDeltaVPred = dDelta;
	}
    in.dTime = dTime;
    in.dDeltaUPred = dDelta;
    in.iRoot = iRoot;
    pstUpdateConsVars(msr->pst,&in,sizeof(in),NULL,NULL);
    printf("(msrUpdateConsVars) End \n");
    }



void msrMeshlessGradients(MSR msr,double dTime,double dDelta,int iRoot){
    printf("Computing gradients... \n");
    if (msr->param.bConservativeReSmooth){
       msrReSmooth(msr,dTime,SMX_SECONDHYDROLOOP,0,0);
    }else{
       msrSmooth(msr,dTime,SMX_SECONDHYDROLOOP,0, msr->param.nSmooth);
    }
}


void msrMeshlessFluxes(MSR msr,double dTime,double dDelta,int iRoot){
#ifdef MAKE_GLASS
    return;
#endif
    printf("Computing fluxes... \n");
    if (msr->param.bConservativeReSmooth){
       if (dDelta==0.0){
          msrReSmooth(msr,dTime,SMX_THIRDHYDROLOOP,1,1); 
       }else{
          msrReSmooth(msr,dTime,SMX_THIRDHYDROLOOP,1,0); 
       }
    }else{
       msrSmooth(msr,dTime,SMX_THIRDHYDROLOOP,0,msr->param.nSmooth);
    }
}


void msrUpdatePrimVars(MSR msr,double dTime,double dDelta,int iRoot){
    struct inDrift in; //IA: TODO new struct for this, as I am using more space than needed
    int nSmoothed = 1, it=0, maxit = 100;  
    in.iRoot = iRoot;
    in.dTime = dTime;
    in.dDelta = dDelta;

    // IA: We also update the particles densities here (i.e., first hydro loop)
    printf("Computing density... \n");
    if (msr->param.bIterativeSmoothingLength){
       msrSelAll(msr); // We set all particles as "not converged"
       msrResetNeighborsStd(msr);

       //pstPredictSmoothing(msr->pst,&in,sizeof(in),NULL,NULL);

       while (nSmoothed>0 && it <= maxit){
          msrSetFirstHydroLoop(msr, 1); // 1-> we care if the particle is marked ; 0-> we dont
          nSmoothed = msrReSmooth(msr,dTime,SMX_FIRSTHYDROLOOP,0,0);
          msrSetFirstHydroLoop(msr, 0);
          it++;
          //if (it>4)
          //   msrIncreaseNeighborsStd(msr);
       }
       if (nSmoothed >0) { /* IA: If we reach the maximum iteration number without full convergence, we compute the smoothing length
                            * as the mean of the surrounding particles which have already converged.
                            * If there are none, the fBall is unmodified 
                            */
                              
         // msrSetFirstHydroLoop(msr, 1); // 1-> we care if the particle is marked ; 0-> we dont
         // nSmoothed = msrReSmooth(msr,dTime,SMX_MEANSMOOTHING,0);
         // msrSetFirstHydroLoop(msr, 0);
         printf("Smoothing length did not converge for %d particles\n", nSmoothed);
       } 
       printf("Computing h took %d iterations \n", it);
    }else{
       msrSetFirstHydroLoop(msr, 1); // 1-> we update the particle's h ; 0-> we dont
       msrSmooth(msr,dTime,SMX_FIRSTHYDROLOOP,0,msr->param.nSmooth);
       msrSetFirstHydroLoop(msr, 0);
    }

    printf("Computing primitive variables... \n");
    pstComputePrimVars(msr->pst,&in,sizeof(in),NULL,NULL); 
}



void msrScaleVel(MSR msr,double dvFac) {
    struct inScaleVel in;

    in.dvFac = dvFac;
    pstScaleVel(msr->pst,&in,sizeof(in),NULL,NULL);
    }

static double ddplus(double a,double omegam,double omegav) {
    double eta;
    if ( a == 0.0 ) return 0.0;
    eta = sqrt(omegam/a + omegav*a*a + 1.0 - omegam - omegav);
    return 2.5/(eta*eta*eta);
    }

#define MAXLEV 20
static double Romberg(double a1, double b1, double eps/* = 1e-8*/,
               double omegam,double omegav) {
    double tllnew;
    double tll;
    double tlk[MAXLEV+1];
    int n = 1;
    int nsamples = 1;

    tlk[0] = tllnew = (b1-a1)*ddplus(0.5*(b1+a1),omegam,omegav);
    if(a1 == b1) return tllnew;

    eps*=0.5;

    do {
        /*
         * midpoint rule.
         */
        double deltax;
        double tlktmp;
        int i;

        nsamples *= 3;
        deltax = (b1-a1)/nsamples;
        tlktmp = tlk[0];
        tlk[0] = tlk[0]/3.0;
        
        for(i=0;i<nsamples/3;i++) {
            tlk[0] += deltax*ddplus(a1 + (3*i + 0.5)*deltax,omegam,omegav);
            tlk[0] += deltax*ddplus(a1 + (3*i + 2.5)*deltax,omegam,omegav);
        }
    
        /*
         * Romberg extrapolation.
         */

        for(i=0;i<n;i++) {
            double tlknew = (pow(9.0, i+1.)*tlk[i] - tlktmp)
                /(pow(9.0, i+1.) - 1.0);
            
            tlktmp = tlk[i+1];
            tlk[i+1] = tlknew;
        }
        tll = tllnew;
        tllnew = tlk[n];
        n++;

    } while((fabs((tllnew-tll)/(tllnew+tll)) > eps) && (n < MAXLEV));
    assert((fabs((tllnew-tll)/(tllnew+tll)) < eps));
    return tllnew;
}

static double dplus(double a,double omegam,double omegav) {
    double eta;
    eta = sqrt(omegam/a + omegav*a*a + 1.0 - omegam - omegav);
    return eta/a * Romberg(0,a,1e-8,omegam,omegav);
}

static double fomega(double a,double omegam,double omegav) {
    double eta, omegak;
    if ( omegam == 1.0 && omegav == 0.0 ) return 1.0;
    omegak=1.0-omegam-omegav;
    eta=sqrt(omegam/a+omegav*a*a+omegak);
    return (2.5/dplus(a,omegam,omegav)-1.5*omegam/a-omegak)/(eta*eta);
}

static double dladt( double a, double omegam, double omegav ) {
    double eta;
    eta=sqrt(omegam/a+omegav*a*a+1.0-omegam-omegav);
    return a * eta;
    }

/*
** Assumes that particles have been drifted using the Zel'dovich approximation.
** This function will drift them to a new time.
*/
double msrAdjustTime(MSR msr, double aOld, double aNew) {
    struct inDrift in;
    double dOmegaM = msr->param.csm->val.dOmega0;
    double dOmegaV = msr->param.csm->val.dLambda;
    double dvFac;
    double dOld, dNew;

    dOld = fomega(aOld,dOmegaM,dOmegaV) * dladt(aOld,dOmegaM,dOmegaV);

    in.dDeltaVPred = 0.0;
    in.dDeltaUPred = 0.0;
    in.iRoot = ROOT;
    msrprintf(msr,"Drifing particles from Time:%g Redshift:%g to Time:%g Redshift:%g ...\n",
	      aOld, 1.0/aOld-1.0, aNew, aNew>0 ? 1.0/aNew-1.0 : 999999.0);
    msrprintf(msr,"WARNING: This only works if the input file is a Zel'dovich perturbed grid\n");
    in.dDelta = -1.0 / (sqrt(8.0/3.0*M_PI)*dOld );
    pstDrift(msr->pst,&in,sizeof(in),NULL,NULL);

    if (aNew > 0.0) {
	dNew = fomega(aNew,dOmegaM,dOmegaV) * dladt(aNew,dOmegaM,dOmegaV);
	dvFac = dOld / dNew * pow(dplus(aNew,dOmegaM,dOmegaV)/dplus(aOld,dOmegaM,dOmegaV),2);
	msrScaleVel(msr,dvFac);
	in.dDelta = 1.0 / (sqrt(8.0/3.0*M_PI)*dNew );
	pstDrift(msr->pst,&in,sizeof(in),NULL,NULL);
	}
    else aNew = 1e-5;

    return getTime(msr,aNew,&dvFac);
    }

/*
 * For gas, updates predicted velocities to beginning of timestep.
 */
void msrKickKDKOpen(MSR msr,double dTime,double dDelta,uint8_t uRungLo,uint8_t uRungHi) {
    struct inKick in;
    struct outKick out;

    in.dTime = dTime;
    if (msr->param.csm->val.bComove) {
	in.dDelta = csmComoveKickFac(msr->param.csm,dTime,dDelta);
	in.dDeltaVPred = 0;
    }
    else {
	in.dDelta = dDelta;
	in.dDeltaVPred = 0;
    }
    in.dDeltaU = dDelta;
    in.dDeltaUPred = 0;
    in.uRungLo = uRungLo;
    in.uRungHi = uRungHi;
    pstKick(msr->pst,&in,sizeof(in),&out,NULL);
    msrprintf(msr,"KickOpen: Avg Wallclock %f, Max Wallclock %f\n",
	      out.SumTime/out.nSum,out.MaxTime);
    }

/*
 * For gas, updates predicted velocities to end of timestep.
 */
void msrKickKDKClose(MSR msr,double dTime,double dDelta,uint8_t uRungLo,uint8_t uRungHi) {
    struct inKick in;
    struct outKick out;

    in.dTime = dTime;
    if (msr->param.csm->val.bComove) {
	in.dDelta = csmComoveKickFac(msr->param.csm,dTime,dDelta);
	in.dDeltaVPred = in.dDelta;
    }
    else {
	in.dDelta = dDelta;
	in.dDeltaVPred = in.dDelta;
    }
    in.dDeltaU = dDelta;
    in.dDeltaUPred = in.dDeltaU;
    in.uRungLo = uRungLo;
    in.uRungHi = uRungHi;
    pstKick(msr->pst,&in,sizeof(in),&out,NULL);
    msrprintf(msr,"KickClose: Avg Wallclock %f, Max Wallclock %f\n",
	      out.SumTime/out.nSum,out.MaxTime);
    }

int msrOutTime(MSR msr,double dTime) {
    if (msr->iOut < msr->nOuts) {
	if (dTime >= msr->pdOutTime[msr->iOut]) {
	    ++msr->iOut;
	    return(1);
	    }
	else return(0);
	}
    else return(0);
    }

int cmpTime(const void *v1,const void *v2) {
    double *d1 = (double *)v1;
    double *d2 = (double *)v2;

    if (*d1 < *d2) return(-1);
    else if (*d1 == *d2) return(0);
    else return(1);
    }

void msrReadOuts(MSR msr,double dTime) {
    char achFile[PST_FILENAME_SIZE];
    char ach[PST_FILENAME_SIZE];
    LCL *plcl = &msr->lcl;
    FILE *fp;
    int i,ret;
    double z,a,n;
    char achIn[80];

    /*
    ** Add Data Subpath for local and non-local names.
    */
    _msrMakePath(msr->param.achDataSubPath,msr->param.achOutName,achFile);
    strcat(achFile,".red");

    fp = fopen(achFile,"r");
    if (!fp) {
	msr->nOuts = 0;
	return;
	}
    i = 0;
    while (1) {
	if (!fgets(achIn,80,fp)) goto NoMoreOuts;
	switch (achIn[0]) {
	case 'z':
	    ret = sscanf(&achIn[1],"%lf",&z);
	    if (ret != 1) goto NoMoreOuts;
	    a = 1.0/(z+1.0);
	    msr->pdOutTime[i] = csmExp2Time(msr->param.csm,a);
	    break;
	case 'a':
	    ret = sscanf(&achIn[1],"%lf",&a);
	    if (ret != 1) goto NoMoreOuts;
	    msr->pdOutTime[i] = csmExp2Time(msr->param.csm,a);
	    break;
	case 't':
	    ret = sscanf(&achIn[1],"%lf",&msr->pdOutTime[i]);
	    if (ret != 1) goto NoMoreOuts;
	    break;
	case 'n':
	    ret = sscanf(&achIn[1],"%lf",&n);
	    if (ret != 1) goto NoMoreOuts;
	    msr->pdOutTime[i] = dTime + (n-0.5)*msrDelta(msr);
	    break;
	default:
	    ret = sscanf(achIn,"%lf",&z);
	    if (ret != 1) goto NoMoreOuts;
	    a = 1.0/(z+1.0);
	    msr->pdOutTime[i] = csmExp2Time(msr->param.csm,a);
	    }
	++i;
	if (i > msr->nMaxOuts) {
	    msr->nMaxOuts *= 2;
	    msr->pdOutTime = realloc(msr->pdOutTime,
				     msr->nMaxOuts*sizeof(double));
	    assert(msr->pdOutTime != NULL);
	    }
	}
NoMoreOuts:
    msr->nOuts = i;
    /*
    ** Now sort the array of output times into ascending order.
    */
    qsort(msr->pdOutTime,msr->nOuts,sizeof(double),cmpTime);
    fclose(fp);
    }


int msrSteps(MSR msr) {
    return(msr->param.nSteps);
    }

char *msrOutName(MSR msr) {
    return(msr->param.achOutName);
    }

char *msrBuildName(MSR msr,char *achFile,int iStep) {
    return _BuildName(msr,achFile,iStep, msr->param.achOutPath);
    }

char *msrBuildIoName(MSR msr,char *achFile,int iStep) {
    if ( msr->param.achIoPath[0] )
	return _BuildName(msr,achFile,iStep, msr->param.achIoPath);
    else
	return msrBuildName(msr,achFile,iStep);
    }



double msrDelta(MSR msr) {
    return(msr->param.dDelta);
    }


int msrLogInterval(MSR msr) {
    return(msr->param.iLogInterval);
    }


int msrOutInterval(MSR msr) {
    return(msr->param.iOutInterval);
    }


const char *msrOutTypes(MSR msr) {
    return(msr->param.achOutTypes);
    }


int msrCheckInterval(MSR msr) {
    return(msr->param.iCheckInterval);
    }


const char *msrCheckTypes(MSR msr) {
    return(msr->param.achCheckTypes);
    }


int msrComove(MSR msr) {
    return(msr->param.csm->val.bComove);
    }


double msrSoft(MSR msr) {
    return(msr->param.dSoft);
    }

/*
** Theta switch. Default is to use dTheta, then switch:
**   at z=20 to dTheta20
**   at z=2 to dTheta2
** We also adjust the number of replicas if the accuracy warrants it.
*/
void msrSwitchTheta(MSR msr,double dTime) {
    double a = csmTime2Exp(msr->param.csm,dTime);
    double dNewTheta;
    if (a < (1.0/21.0)) dNewTheta = msr->param.dTheta;
    else if (a < (1.0/3.0)) dNewTheta = msr->param.dTheta20;
    else dNewTheta = msr->param.dTheta2;
    if (msr->dThetaMin != dNewTheta) {
	msr->dThetaMin = dNewTheta;
	if (msr->param.bVDetails)
	    printf("Theta changed to %g at z=%g\n",dNewTheta,1.0/a - 1.0);
	}
    if ( !prmSpecified(msr->prm,"nReplicas") && msr->param.nReplicas>=1 ) {
	if ( msr->dThetaMin < 0.52 ) msr->param.nReplicas = 2;
	else msr->param.nReplicas = 1;
	}
    }

void msrInitStep(MSR msr) {
    struct inSetRung insr;
    static/*FIXME: this is a HACK: message is too large for stack */ struct inInitStep in;

    /*
    ** Here we can pass down all parameters of the simulation
    ** before any timestepping takes place. This should happen
    ** just after the file has been read and the PKD structure
    ** initialized for each processor.
    */
    in.param = msr->param;
    in.cosmo = msr->param.csm->val;
    pstInitStep(msr->pst, &in, sizeof(in), NULL, NULL);

    /*
    ** Initialize particles to lowest rung. (what for? JW: Seconded and removed)
    */
//    insr.uRung = 0; /* msr->param.iMaxRung - 1; */
//    insr.uRungLo = 0;
//    insr.uRungHi = MAX_RUNG;
//    pstSetRung(msr->pst, &insr, sizeof(insr), NULL, NULL);
//    msr->iCurrMaxRung = insr.uRung;
    }


void msrSetRung(MSR msr, uint8_t uRungLo, uint8_t uRungHi, int uRung) {
    struct inSetRung in;

    in.uRung = uRung;
    in.uRungLo = uRungLo;
    in.uRungHi = uRungHi;
    pstSetRung(msr->pst, &in, sizeof(in), NULL, NULL);
    msr->iCurrMaxRung = in.uRung;
    }


void msrZeroNewRung(MSR msr, uint8_t uRungLo, uint8_t uRungHi, int uRung) {
    struct inZeroNewRung in;

    in.uRung = uRung;
    in.uRungLo = uRungLo;
    in.uRungHi = uRungHi;
    pstZeroNewRung(msr->pst, &in, sizeof(in), NULL, NULL);
    }


int msrMaxRung(MSR msr) {
    return msr->param.iMaxRung;
    }


int msrCurrMaxRung(MSR msr) {
    return msr->iCurrMaxRung;
    }


double msrEta(MSR msr) {
    return msr->param.dEta;
    }

/*
 * bGreater = 1 => activate all particles at this rung and greater.
 */
void msrActiveRung(MSR msr, int iRung, int bGreater) {
    struct inActiveRung in;

    in.iRung = iRung;
    in.bGreater = bGreater;
    pstActiveRung(msr->pst, &in, sizeof(in), NULL, NULL);

    if ( iRung==0 && bGreater )
	msr->nActive = msr->N;
    else {
	int i;

	assert( msr->nRung != NULL );

	msr->nActive = 0;
	for ( i=iRung; i<= (bGreater?msr->param.iMaxRung:iRung); i++ )
	    msr->nActive += msr->nRung[i];
	}
    }

void msrActiveOrder(MSR msr) {
    pstActiveOrder(msr->pst,NULL,0,&(msr->nActive),NULL);
    }

void msrSetRungVeryActive(MSR msr, int iRung) {
    struct inSetRung in;

    msr->iRungVeryActive = iRung;

    in.uRung = iRung;
    in.uRungLo = 0;
    in.uRungHi = MAX_RUNG;
    pstSetRungVeryActive(msr->pst,&in,sizeof(in),NULL,NULL);
    }

int msrCountRungs(MSR msr, uint64_t *nRungs) {
    struct outCountRungs out;
    int i, iMaxRung=0;
    pstCountRungs(msr->pst, NULL, 0, &out, NULL);
    for(i=0; i<=MAX_RUNG; ++i) {
	msr->nRung[i] = out.nRungs[i];
	if (msr->nRung[i]) iMaxRung = i;
	if (nRungs) nRungs[i] = msr->nRung[i];
	}
    msr->iCurrMaxRung = iMaxRung;

    const uint64_t nDT = d2u64(msr->N*msr->param.dFracDualTree);
    const uint64_t nDD = d2u64(msr->N*msr->param.dFracNoDomainDecomp);
    uint64_t nActive = 0;
    msr->iRungDD = 0;
    msr->iRungDT = 0;
    for (i=msr->iCurrMaxRung;i>=0;--i) {
	nActive += msr->nRung[i];
	if (nActive > nDT && !msr->iRungDT) msr->iRungDT = i;
	if (nActive > nDD && !msr->iRungDD) msr->iRungDD = i;
	}

    return iMaxRung;
    }

void msrAccelStep(MSR msr,uint8_t uRungLo,uint8_t uRungHi,double dTime) {
    struct inAccelStep in;
    double a;

    in.dEta = msrEta(msr);
    a = csmTime2Exp(msr->param.csm,dTime);
    if (msr->param.csm->val.bComove) {
	in.dVelFac = 1.0/(a*a);
	}
    else {
	in.dVelFac = 1.0;
	}
    in.dAccFac = 1.0/(a*a*a);
    in.bDoGravity = msrDoGravity(msr);
    in.bEpsAcc = msr->param.bEpsAccStep;
    in.uRungLo = uRungLo;
    in.uRungHi = uRungHi;
    pstAccelStep(msr->pst,&in,sizeof(in),NULL,NULL);
    }

/* Requires full forces and full udot (i.e. sph and cooling both done) */
void msrSphStep(MSR msr,uint8_t uRungLo,uint8_t uRungHi,double dTime) {
    struct inSphStep in;
    double a;

    a = csmTime2Exp(msr->param.csm,dTime);
    in.dAccFac = 1.0/(a*a*a);
    in.uRungLo = uRungLo;
    in.uRungHi = uRungHi;
    pstSphStep(msr->pst,&in,sizeof(in),NULL,NULL);
    }

/* IA: Computes the dt criteria being known the fluxes and the signal velocities */
void msrHydroStep(MSR msr,uint8_t uRungLo,uint8_t uRungHi,double dTime) {
    int bSymmetric = 0; 

    printf("Computing hydro time step... \n");
    msrReSmooth(msr,dTime,SMX_HYDROSTEP,1, 0);

    if (msr->param.bGlobalDt){
       if (msr->param.dFixedDelta != 0.0){
          msrSetGlobalDt(msr, msr->param.dFixedDelta);
       }else{
          uint8_t minDt;
          minDt = msrGetMinDt(msr);
          msrSetGlobalDt(msr, minDt);
       }
    }
}

uint8_t msrGetMinDt(MSR msr){
    struct outGetMinDt out;
    int nOut;

    pstGetMinDt(msr->pst, NULL, 0, &out, &nOut);
    return out.uMinDt;
    }

void msrSetGlobalDt(MSR msr, uint8_t minDt){
   struct outGetMinDt in;
   in.uMinDt = minDt;

   pstSetGlobalDt(msr->pst, &in, sizeof(in), NULL, NULL); 

   }

void msrDensityStep(MSR msr,uint8_t uRungLo,uint8_t uRungHi,double dTime) {
    struct inDensityStep in;
    double expand;
    int bSymmetric;

    msrprintf(msr,"Calculating Rung Densities...\n");
    bSymmetric = 0;
    msrSmooth(msr,dTime,SMX_DENSITY,bSymmetric,msr->param.nSmooth);
    in.dEta = msrEta(msr);
    expand = csmTime2Exp(msr->param.csm,dTime);
    in.dRhoFac = 1.0/(expand*expand*expand);
    in.uRungLo = uRungLo;
    in.uRungHi = uRungHi;
    pstDensityStep(msr->pst,&in,sizeof(in),NULL,NULL);
    }

/*
** Updates all particles in the tree at iRoot.
** Particles may not go higer (lower rung number) than uRung.
*/
void msrUpdateRungByTree(MSR msr, uint8_t uRung, int iRoot) {
    struct inUpdateRungByTree in;
    struct outUpdateRung out;
    int iTempRung;

    in.iRoot = iRoot;
    in.uMinRung = uRung;
    in.uMaxRung = msrMaxRung(msr);
    pstUpdateRungByTree(msr->pst, &in, sizeof(in), &out, NULL);
    for (iTempRung=uRung;iTempRung < msrMaxRung(msr);++iTempRung) msr->nRung[iTempRung] = out.nRungCount[iTempRung];
    if (msr->param.bVRungStat) {
	printf("Rung distribution:\n");
	printf("\n");
	for (iTempRung=0;iTempRung <= msr->iCurrMaxRung;++iTempRung) {
	    if (msr->nRung[iTempRung] == 0) continue;
	    printf("   rung:%d %"PRIu64"\n",iTempRung,msr->nRung[iTempRung]);
	    }
	printf("\n");
	}


    }

/*
** Returns the Very Active rung based on the number of very active particles desired,
** or the fixed rung that was specified in the parameters.
*/
int msrUpdateRung(MSR msr, uint8_t uRung) {
    struct inUpdateRung in;
    struct outUpdateRung out;
    int iTempRung,iOutMaxRung,iRungVeryActive;
    uint64_t sum;
    char c;

    /* If we are called, it is a mistake -- this happens in analysis mode */
    if (msr->param.bMemUnordered) return 0;

    in.uRungLo = uRung;
    in.uRungHi = msrMaxRung(msr);
    in.uMinRung = uRung;
    in.uMaxRung = msrMaxRung(msr);

    pstUpdateRung(msr->pst, &in, sizeof(in), &out, NULL);

    iTempRung =msrMaxRung(msr)-1;
    while (out.nRungCount[iTempRung] == 0 && iTempRung > 0) --iTempRung;
    iOutMaxRung = iTempRung;

    while (out.nRungCount[iOutMaxRung] <= msr->param.nTruncateRung && iOutMaxRung > uRung) {
	msrprintf(msr,"n_CurrMaxRung = %"PRIu64"  (iCurrMaxRung = %d):  Promoting particles to iCurrMaxrung = %d\n",
		  out.nRungCount[iOutMaxRung],iOutMaxRung,iOutMaxRung-1);

	in.uMaxRung = iOutMaxRung; /* Note this is the forbidden rung so no -1 here */
	pstUpdateRung(msr->pst, &in, sizeof(in), &out, NULL);

	iTempRung =msrMaxRung(msr)-1;
	while (out.nRungCount[iTempRung] == 0 && iTempRung > 0) --iTempRung;
	iOutMaxRung = iTempRung;
	}
    /*
    ** Now copy the rung distribution to the msr structure!
    */
    for (iTempRung=0;iTempRung < msrMaxRung(msr);++iTempRung) msr->nRung[iTempRung] = out.nRungCount[iTempRung];
    /*
    ** Now we want to make a suggestion for the current very active rung based on the number of
    ** particles in the deepest rungs.
    */
    if (msr->param.nPartVeryActive > 0) {
	iRungVeryActive = msrMaxRung(msr);
	sum = 0;
	while (sum < msr->param.nPartVeryActive && iRungVeryActive > 0) {
	    sum += out.nRungCount[--iRungVeryActive];
	    }
	}
    else {
	iRungVeryActive = msr->param.nRungVeryActive;
	}

    msr->iCurrMaxRung = iOutMaxRung;

    if (msr->param.bVRungStat) {
	printf("Rung distribution:\n");
	printf("\n");
	for (iTempRung=0;iTempRung <= msr->iCurrMaxRung;++iTempRung) {
	    if (out.nRungCount[iTempRung] == 0) continue;
	    if (iTempRung > iRungVeryActive) c = 'v';
	    else c = ' ';
	    printf(" %c rung:%d %"PRIu64"\n",c,iTempRung,out.nRungCount[iTempRung]);
	    }
	printf("\n");
	}

    /*
     * Set VeryActive particles
     */
    msrSetRungVeryActive(msr, iRungVeryActive);
    return(iRungVeryActive);
    }


/*
** Open the healpix output file, and also the particles files if requested.
*/
void msrLightConeOpen(MSR msr, int iStep) {
    if (msr->param.bLightCone) {
	struct inLightConeOpen lc;
	if (msr->param.bLightConeParticles )
	    msrBuildName(msr,lc.achOutFile,iStep);
	else lc.achOutFile[0] = 0;
	lc.nSideHealpix = msr->param.nSideHealpix;
	pstLightConeOpen(msr->pst,&lc,sizeof(lc),NULL,NULL);
	}
    }


/*
** Close the files for this step.
*/
void msrLightConeClose(MSR msr,int iStep) {
    if (msr->param.bLightCone) {
	struct inLightConeClose lc;
	msrBuildName(msr,lc.achOutFile,iStep);
	pstLightConeClose(msr->pst,&lc,sizeof(lc),NULL,NULL);
	}
    }

/*
** Correct velocities from a^2 x_dot to a x_dot (physical peculiar velocities) using the 
** position dependent scale factor within the light cone. This could be expensive.
*/
void msrLightConeVel(MSR msr) {
    double sec,dsec;

    sec = msrTime();
    pstLightConeVel(msr->pst,NULL,0,NULL,NULL);
    dsec = msrTime() - sec;
    printf("Converted lightcone velocities to physical, Wallclock: %f secs.\n", dsec);
    }



void msrCheckForOutput(MSR msr,int iStep,double dTime,int *pbDoCheckpoint,int *pbDoOutput) {
    int iStop, iCheck;
    long lSec;

    /*
    ** Check for user interrupt.
    */
    iStop = msrCheckForStop(msr,STOPFILE);
    iCheck = msrCheckForStop(msr,CHECKFILE);

    /*
    ** Check to see if the runtime has been exceeded.
    */
    if (!iStop && msr->param.iWallRunTime > 0) {
	if (msr->param.iWallRunTime*60 - (time(0)-msr->lStart) < ((int) (lSec*1.5)) ) {
	    printf("RunTime limit exceeded.  Writing checkpoint and exiting.\n");
	    printf("    iWallRunTime(sec): %d   Time running: %ld   Last step: %ld\n",
		msr->param.iWallRunTime*60,time(0)-msr->lStart,lSec);
	    iStop = 1;
	    }
	}

    /* Check to see if there should be an output */
    if (!iStop && timeGlobalSignalTime>0) { /* USR1 received */
	if ( (time(0)+(lSec*1.5)) > timeGlobalSignalTime+msr->param.iSignalSeconds) {
	    printf("RunTime limit exceeded.  Writing checkpoint and exiting.\n");
	    printf("    iSignalSeconds: %d   Time running: %ld   Last step: %ld\n",
		msr->param.iSignalSeconds,time(0)-msr->lStart,lSec);
	    iStop = 1;
	    }
	}

    /*
    ** Output if 1) we've hit an output time
    **           2) We are stopping
    **           3) we're at an output interval
    */
    if (iCheck || (msrCheckInterval(msr)>0 &&
		    (bGlobalOutput
		    || iStop
			|| (iStep%msrCheckInterval(msr) == 0) )) ) {
	bGlobalOutput = 0;
	*pbDoCheckpoint = 1;
	}

    if (msrOutTime(msr,dTime) 
	|| (msrOutInterval(msr) > 0 &&
		(bGlobalOutput
		|| iStop
		|| iStep == msrSteps(msr)
		|| (iStep%msrOutInterval(msr) == 0))) ) {
	bGlobalOutput = 0;
	*pbDoOutput = 1;
	}
    }


int msrNewTopStepKDK(MSR msr,
    int bDualTree,      /* Should be zero at rung 0! */
    uint8_t uRung,	/* Rung level */
    double *pdStep,	/* Current step */
    double *pdTime,	/* Current time */
    uint8_t *puRungMax,
    int *piSec,int *pbDoCheckpoint,int *pbDoOutput) {
    uint64_t nActive;
    double dDelta,dTimeFixed;
    uint32_t uRoot2=0;
    int iRungDT = msr->iRungDT;
    char achFile[256];
    int bKickOpen;
    /*
    ** The iStep variable serves only to give a number to the lightcone and group output files.
    ** We define this to be the output number of the final radius of the lightcone surface.
    */
    int iStep = (int)(*pdStep) + 1;

    if (uRung == iRungDT+1) {
	if ( msr->param.bDualTree && uRung < *puRungMax) {
	    /* HACK: FIXME: Don't use the dual tree before z=2; the overlap region is too large */
	    /* better would be to construct the tree matching remote processor shape as well as local */
	    double a = csmTime2Exp(msr->param.csm,*pdTime);
	    if (a < (1.0/3.0)) bDualTree = 0;
	    else {
		bDualTree = 1;
		struct inDumpTrees dump;
		dump.bOnlyVA = 0;
		dump.uRungDD = iRungDT;
		pstDumpTrees(msr->pst,&dump,sizeof(dump),NULL,NULL);
		msrprintf(msr,"Half Drift, uRung: %d\n",iRungDT);
		dDelta = msr->param.dDelta/(1 << iRungDT); // Main tree step
		msrDrift(msr,*pdTime,0.5 * dDelta,FIXROOT);
		dTimeFixed = *pdTime + 0.5 * dDelta;
		msrBuildTreeFixed(msr,*pdTime,msr->param.bEwald,iRungDT);
		}
	    }
	else bDualTree = 0;
	}
    if (uRung < *puRungMax) {
	bDualTree = msrNewTopStepKDK(msr,bDualTree,uRung+1,pdStep,pdTime,puRungMax,piSec,pbDoCheckpoint,pbDoOutput);
	}

    /* Drift the "ROOT" (active) tree or all particle */
    dDelta = msr->param.dDelta/(1 << *puRungMax);
    if (bDualTree) {
	msrprintf(msr,"Drift very actives, uRung: %d\n",*puRungMax);
	msrDrift(msr,*pdTime,dDelta,ROOT);
	}
    else {
	msrprintf(msr,"Drift, uRung: %d\n",*puRungMax);
	msrDrift(msr,*pdTime,dDelta,-1);
	}
    *pdTime += dDelta;
    *pdStep += 1.0/(1 << *puRungMax);

    msrActiveRung(msr,uRung,1);
    msrUpdateSoft(msr,*pdTime);
    if (bDualTree && uRung > iRungDT) {
	uRoot2 = FIXROOT;
	msrBuildTreeActive(msr,*pdTime,msr->param.bEwald,iRungDT);
	}
    else {
	if (uRung==0) msrInflate(msr,round(*pdStep));
	msrDomainDecomp(msr,uRung,0,0);
	uRoot2 = 0;
	msrBuildTree(msr,*pdTime,msr->param.bEwald);
	}

    if (!uRung && msr->param.iPkInterval && iStep%msr->param.iPkInterval == 0) {
	msrOutputPk(msr,iStep,*pdTime);
	}

    if (!uRung) {
	msrCheckForOutput(msr,iStep,*pdTime,pbDoCheckpoint,pbDoOutput);	
	if (*pbDoCheckpoint || *pbDoOutput) bKickOpen = 0;
	else bKickOpen = 1;
	}
    else bKickOpen = 1;
    
    /*
    ** We need to write all light cone files (healpix and LCP) at this point before the last
    ** gravity is called since it will advance the particles in the light cone as part of the
    ** opening kick! We also need to open
    */
    if (!uRung) {
	msrLightConeClose(msr,iStep);
	if (bKickOpen) msrLightConeOpen(msr,iStep+1);
	}

    /* Compute the grids of linear species at main timesteps, before gravity is called */
    if (!uRung && strlen(msr->param.csm->val.classData.achLinSpecies) && msr->param.nGridLin){
	msrGridCreateFFT(msr,msr->param.nGridLin);
        msrSetLinGrid(msr, *pdTime, msr->param.nGridLin,1,bKickOpen);
        if (msr->param.bDoLinPkOutput)
            msrOutputLinPk(msr, *pdStep, *pdTime);
	msrLinearKick(msr,*pdTime,1,bKickOpen);
	msrGridDeleteFFT(msr);
	}

    if (!uRung && msr->param.bFindGroups) msrNewFof(msr,*pdTime);

    // We need to make sure we descend all the way to the bucket with the
    // active tree, or we can get HUGE group cells, and hence too much P-P/P-C
    int nGroup = (bDualTree && uRung > iRungDT) ? 1 : msr->param.nGroup;
    *puRungMax = msrGravity(msr,uRung,msrMaxRung(msr),ROOT,uRoot2,*pdTime,
	*pdStep,1,bKickOpen,msr->param.bEwald,nGroup,piSec,&nActive);

    if (!uRung && msr->param.bFindGroups) {
	msrGroupStats(msr);
	msrBuildName(msr,achFile,iStep);
	strncat(achFile,".fofstats",256);
	msrHopWrite(msr,achFile);
	}

    if (uRung && uRung < *puRungMax) bDualTree = msrNewTopStepKDK(msr,bDualTree,uRung+1,pdStep,pdTime,puRungMax,piSec,pbDoCheckpoint,pbDoOutput);
    if (bDualTree && uRung==iRungDT+1) {
	msrprintf(msr,"Half Drift, uRung: %d\n",iRungDT);
	dDelta = msr->param.dDelta/(1 << iRungDT);
	msrDrift(msr,dTimeFixed,0.5 * dDelta,FIXROOT);
	}
    return bDualTree;
    }

void msrTopStepKDK(MSR msr,
		   double dStep,	/* Current step */
		   double dTime,	/* Current time */
		   double dDelta,	/* Time step */
		   int iRung,		/* Rung level */
		   int iKickRung,	/* Gravity on all rungs from iRung
					   to iKickRung */
		   int iRungVeryActive,  /* current setting for iRungVeryActive */
		   /*
		   ** Note that iRungVeryActive is one less than the first rung with VA particles!
		   */
		   int iAdjust,		/* Do an adjust? */
		   double *pdActiveSum,
		   int *piSec) {
    uint64_t nActive;
    int bSplitVA;

    if (iAdjust && (iRung < msrMaxRung(msr)-1)) {
	msrprintf(msr,"%*cAdjust, iRung: %d\n",2*iRung+2,' ',iRung);
	/* JW: Note -- can't trash uRungNew here! Force calcs set values for it! */
	msrActiveRung(msr, iRung, 1);
	if (msr->param.bAccelStep) {
	    msrAccelStep(msr,iRung,MAX_RUNG,dTime);
	    }
	if (msrDoGas(msr)) {
          if (!msrMeshlessHydro(msr)) {
	       msrSphStep(msr,iRung,MAX_RUNG,dTime);
	       }else{
             msrHydroStep(msr,iRung, MAX_RUNG, dTime);
          }
      }    
	if (msr->param.bDensityStep) {
	    bSplitVA = 0;
	    msrDomainDecomp(msr,iRung,0,bSplitVA);
	    msrActiveRung(msr,iRung,1);
	    msrBuildTree(msr,dTime,0);
	    msrDensityStep(msr,iRung,MAX_RUNG,dTime);
	    }
	iRungVeryActive = msrUpdateRung(msr,iRung);
        }

    msrprintf(msr,"%*cmsrKickOpen  at iRung: %d 0.5*dDelta: %g\n",
	      2*iRung+2,' ',iRung,0.5*dDelta);
    msrKickKDKOpen(msr,dTime,0.5*dDelta,iRung,iRung);
    //msrUpdateConsVars(msr, dTime, dDelta, ROOT);
      if (dTime == 1 && msrDoGas(msr) && msrMeshlessHydro(msr) && msrDoGravity(msr)){
    //      msrApplyGravWork(msr,dTime,0.0,iRung,iRung); 
      }
    if ((msrCurrMaxRung(msr) > iRung) && (iRungVeryActive > iRung)) {
	/*
	** Recurse.
	*/
	msrTopStepKDK(msr,dStep,dTime,0.5*dDelta,iRung+1,iRung+1,iRungVeryActive,0,
		      pdActiveSum,piSec);
	dTime += 0.5*dDelta;
	dStep += 1.0/(2 << iRung);

	msrActiveRung(msr,iRung,0); /* is this call even needed? */

	msrTopStepKDK(msr,dStep,dTime,0.5*dDelta,iRung+1,iKickRung,iRungVeryActive,1,
		      pdActiveSum,piSec);
	}
    else if (msrCurrMaxRung(msr) == iRung) {
      if (msrDoGas(msr) && msrMeshlessHydro(msr)){
         msrActiveRung(msr,iKickRung,1); //IA: The repeated call after msrDrift would not be needed
         if (msr->param.bVStep) printf("Step:%f (iKickRung %d) (iRung %d) \n",dStep,iKickRung, iRung);
         msrUpdateConsVars(msr, dTime, dDelta, ROOT);
         msrMeshlessFluxes(msr, dTime, dDelta, ROOT);
      }


	/* This Drifts everybody */
	msrprintf(msr,"%*cDrift, iRung: %d\n",2*iRung+2,' ',iRung);
	msrDrift(msr,dTime,dDelta,ROOT);
	dTime += dDelta;
	dStep += 1.0/(1 << iRung);
      if (msr->param.csm->val.bComove){
         const float a = csmTime2Exp(msr->param.csm,dTime);
#ifdef COOLING
         const float z = 1./a - 1.;

         msrCoolingUpdate(msr, z);
#endif

      }

      printf("dTime %e \n", dTime);
	msrActiveRung(msr,iKickRung,1);
	bSplitVA = 0;
	msrDomainDecomp(msr,iKickRung,0,bSplitVA);

	/* JW: Good place to zero uNewRung  IA: not really... This is done in the init call of smSmooth*/ 
	msrZeroNewRung(msr,iKickRung,MAX_RUNG,iKickRung); /* brute force */


      if (msrDoGas(msr) && msrMeshlessHydro(msr)){
         //msrUpdateConsVars(msr, dTime, dDelta, ROOT);
      }

	if (msrDoGravity(msr) || msrDoGas(msr)) {
	    msrActiveRung(msr,iKickRung,1);
	    if (msrDoGravity(msr)) msrUpdateSoft(msr,dTime);
	    msrprintf(msr,"%*cForces, iRung: %d to %d\n",2*iRung+2,' ',iKickRung,iRung);
	    msrBuildTree(msr,dTime,msr->param.bEwald);
	    }
	if (msrDoGravity(msr)) {
	    msrGravity(msr,iKickRung,MAX_RUNG,ROOT,0,dTime,dStep,0,0,msr->param.bEwald,msr->param.nGroup,piSec,&nActive);
	    *pdActiveSum += (double)nActive/msr->N;
	    }
	
	if (msrDoGas(msr) && !msrMeshlessHydro(msr)) {
	    msrSph(msr,dTime,dStep);  /* dTime = Time at end of kick */
	    msrCooling(msr,dTime,dStep,0,
		       (iKickRung<=msr->param.iRungCoolTableUpdate ? 1:0),0);
	}


      if (msrDoGas(msr) && msrMeshlessHydro(msr)){
         msrUpdatePrimVars(msr, dTime, dDelta, ROOT);
         msrMeshlessGradients(msr, dTime, dDelta, ROOT);
      }

      if (msrDoGas(msr) && msrMeshlessHydro(msr) && msrDoGravity(msr)){
          //msrApplyGravWork(msr,dTime,dDelta,iKickRung,iRung); // IA: dDelta is the step for iRung, but iKickRung *may* be lower
      }
	/*
	 * move time back to 1/2 step so that KickClose can integrate
	 * from 1/2 through the timestep to the end.
	 */
	dTime -= 0.5*dDelta;
	}
    else { /* Very active particles */
	double dDeltaTmp;
	int i;

	assert(!msrDoGas(msr)); /* Gas doesn't handle VA yet */
	/*
	 * We have more rungs to go, but we've hit the very active limit.
	 */

	/*
	 * Drift the non-VeryActive particles forward 1/2 timestep
	 */
	msrprintf(msr,"%*cInActiveDrift at iRung: %d, 0.5*dDelta: %g\n",
		  2*iRung+2,' ',iRung,0.5*dDelta);
	msrDrift(msr,dTime,0.5*dDelta,ROOT);
	/*
	 * Build a tree out of them for use by the VeryActives
	 */
	if (msrDoGravity(msr)) {
	    msrUpdateSoft(msr,dTime + 0.5*dDelta);
	    /*
	    ** Domain decomposition for parallel exclude very active is going to be
	    ** placed here shortly.
	    */
	    bSplitVA = 1;
	    msrDomainDecomp(msr,iRung,0,bSplitVA);

	    msrprintf(msr,"%*cBuilding exclude very active tree: iRung: %d\n",
		      2*iRung+2,' ',iRung);
	    /*
	     * Activate VeryActives, this is needed for the BuildTreeExcludeVeryActive below?
	     */
	    msrActiveRung(msr,msr->iRungVeryActive+1,1);
	    msrBuildTreeExcludeVeryActive(msr,dTime + 0.5*dDelta);
	    }
	/*
	 * Perform timestepping on individual processors.
	 */
	msrStepVeryActiveKDK(msr, dStep, dTime, dDelta, iRung);
	dTime += dDelta;
	dStep += 1.0/(1 << iRung);
	/*
	 * Move Inactives to the end of the step.
	 */
	msrprintf(msr,"%*cInActiveDrift at iRung: %d, 0.5*dDelta: %g\n",
		  2*iRung+2,' ',iRung,0.5*dDelta);
	/*
	** The inactives are half time step behind the actives.
	** Move them a half time step ahead to synchronize everything again.
	*/
	msrDrift(msr,dTime-0.5*dDelta,0.5*dDelta,ROOT);

	/*
	 * Regular Tree gravity
	 */
	msrActiveRung(msr,iKickRung,1);
	bSplitVA = 0;
	msrDomainDecomp(msr,iKickRung,0,bSplitVA);

	if (msrDoGravity(msr)) {
	    msrActiveRung(msr,iKickRung,1);
	    msrUpdateSoft(msr,dTime);
	    msrprintf(msr,"%*cGravity, iRung: %d to %d\n",
		      2*iRung+2,' ',iKickRung,msrCurrMaxRung(msr));
	    msrBuildTree(msr,dTime,msr->param.bEwald);
	    msrGravity(msr,iKickRung,MAX_RUNG,ROOT,0,dTime,dStep,0,0,msr->param.bEwald,msr->param.nGroup,piSec,&nActive);
	    *pdActiveSum += (double)nActive/msr->N;
	    }

	dDeltaTmp = dDelta;
	for (i = msrCurrMaxRung(msr); i > iRung; i--)
	    dDeltaTmp *= 0.5;

	for (i = msrCurrMaxRung(msr); i > iRung; i--) { /* close off all
							 the VeryActive Kicks
						      */
	    msrprintf(msr,"%*cVeryActive msrKickClose at iRung: %d, 0.5*dDelta: %g\n",
		      2*iRung+2,' ',i, 0.5*dDeltaTmp);
	    msrKickKDKClose(msr,dTime-0.5*dDeltaTmp,0.5*dDeltaTmp,i,i);
	    dDeltaTmp *= 2.0;
	    }
	/*
	 * move time back to 1/2 step so that KickClose can integrate
	 * from 1/2 through the timestep to the end.
	 */
	dTime -= 0.5*dDelta;
	 } /* END of very active particles */

    msrprintf(msr,"%*cKickClose, iRung: %d, 0.5*dDelta: %g\n",
	      2*iRung+2,' ',iRung, 0.5*dDelta);
    msrKickKDKClose(msr,dTime,0.5*dDelta,iRung,iRung); /* uses dTime-0.5*dDelta */
//      if (msrDoGas(msr) && msrMeshlessHydro(msr) && msrDoGravity(msr)){
//          msrApplyGravWork(msr,dTime,0.5*dDelta,iRung,iRung);  
//      }


    dTime += 0.5*dDelta; /* Important to have correct time at step end for SF! */

//      if (msrDoGas(msr) && msrMeshlessHydro(msr)){
//	   msrActiveRung(msr,iRung,0);
//         msrUpdatePrimVars(msr, dTime, dDelta, ROOT);
//         msrMeshlessGradients(msr, dTime, dDelta, ROOT);
//      }



    /* JW: Creating/Deleting/Merging is best done outside (before or after) KDK cycle 
       -- Tree should still be valid from last force eval (only Drifts + Deletes invalidate it) */
//    if (iKickRung == iRung) { /* Do all co-incident kicked p's in one go */
//	msrStarForm(msr, dTime, iKickRung); /* re-eval timesteps as needed for next step 
//					    -- apply to all neighbours */
//	}

#ifdef STAR_FORMATION
      msrStarForm(msr, dTime, iKickRung);
#endif

    }



//void msrStarForm(MSR msr, double dTime, int iRung)
//    {
//    struct inStarForm in;
//    struct outStarForm out;
//    double a,d1,d2;
//    double sec,sec1,dsec;
//
//    if(!msr->param.bStarForm) return;
//    sec = msrTime();
//
//    a = csmTime2Exp(msr->param.csm,dTime);
//    in.dDelta = 0;
//    /* Convert input parameters to code units */
//    in.dRateCoeff = msr->param.SFdEfficiency*sqrt(32/(3*M_PI)/pow(a,3)); /* G=1 */
//    in.dTMax = msr->param.SFdTMax;
//    d1 = msr->param.SFdComovingDenMin;
//    d2 = msr->param.SFdPhysDenMin/msr->param.dGmPerCcUnit*pow(a,3);
//    in.dDenMin = (d1>d2 ? d1 : d2);
//    
//    in.dTime = dTime;
//    in.dInitStarMass = msr->param.SFdInitStarMass;
//    in.dESNPerStarMass = msr->param.SFdESNPerStarMass/msr->param.dErgPerGmUnit;
//#define SECONDSPERYEAR   31557600.
//    in.dtCoolingShutoff = msr->param.SFdtCoolingShutoff*SECONDSPERYEAR/msr->param.dSecUnit;
//    in.dtFeedbackDelay = msr->param.SFdtFeedbackDelay*1.0000013254678*SECONDSPERYEAR/msr->param.dSecUnit;
//    in.dMassLossPerStarMass = msr->param.SFdMassLossPerStarMass;
//    in.dZMassPerStarMass = msr->param.SFdZMassPerStarMass;
//    in.dInitStarMass = msr->param.SFdInitStarMass;
//    in.dMinGasMass = msr->param.SFdMinGasMass;
//    in.bdivv = msr->param.SFbdivv;
//    
//    if (msr->param.bVDetails) printf("Star Form ... ");
//    
//    msrActiveRung(msr,iRung,1); /* important to limit to active gas only */
//    pstStarForm(msr->pst, &in, sizeof(in), &out, NULL);
//    if (msr->param.bVDetails)
//	printf("%d Stars formed with mass %g, %d gas deleted\n",
//	       out.nFormed, out.dMassFormed, out.nDeleted);
//    
//    if (out.nDeleted) {
//	//msrSelSrcGas(msr); /* Not really sure what the setting here needs to be */
//	//msrSelDstDeleted(msr); /* Select only deleted particles */
//	msrActiveRung(msr,0,1); /* costs nothing -- may be redundant */
///*	msrBuildTree(msr,dTime,msr->param.bEwald);*/
//	msrSmooth(msr, dTime, SMX_DIST_DELETED_GAS, 1,msr->param.nSmooth); /* use full smooth to account for deleted */
//	//msrSelSrcAll(msr);
//	//msrSelDstAll(msr);
//	}
//
//    /* Strictly speaking adding/deleting particles invalidates the tree 
//       In practice we can find deleted particles (and ignore them) and we aren't looking for 
//       any new star particles so we should be able to function with an old tree for FB 
//       Could do FB before the AddDel call if careful to exclude deleted particles */
//    if (out.nDeleted || out.nFormed) msrAddDelParticles(msr);
//    
//    sec1 = msrTime();
//    dsec = sec1 - sec;
//    printf("Star Formation Calculated, Wallclock: %f secs\n\n",dsec);
//
//    if (msr->param.bFeedback) {
//	msrActiveRung(msr,iRung,1); /* costs nothing -- important to limit to active stars only */
// 	//msrSelSrcGas(msr); /* Not really sure what the setting here needs to be */
//	//msrSelDstStar(msr,1,dTime); /* Select only stars that have FB to do */ 
///*	msrBuildTree(msr,dTime,msr->param.bEwald);*/
//	msrSmooth(msr, dTime, SMX_DIST_SN_ENERGY, 1, msr->param.nSmooth); /* full smooth for stars */
//	//msrSelSrcAll(msr);
//	//msrSelDstAll(msr);
//
//	dsec = msrTime() - sec1;
//	printf("Feedback Calculated, Wallclock: %f secs\n\n",dsec);
//	}
//    }

void msrStepVeryActiveKDK(MSR msr, double dStep, double dTime, double dDelta,
		     int iRung) {
    struct inStepVeryActive in;
    struct outStepVeryActive out;

    in.dStep = dStep;
    in.dTime = dTime;
    in.dDelta = dDelta;
    in.iRung = iRung;
    in.nMaxRung = msrCurrMaxRung(msr);
    /*TODO: Are the next two lines okay?  nMaxRung and iRung needed? */
    in.uRungLo = iRung;
    in.uRungHi = msrCurrMaxRung(msr);
    /*
     * Start Particle Cache on all nodes (could be done as part of
     * tree build)
     */
    pstROParticleCache(msr->pst, NULL, 0, NULL, NULL);

    pstStepVeryActiveKDK(msr->pst, &in, sizeof(in), &out, NULL);
    /*
     * Finish Particle Cache on all nodes
     */
    pstParticleCacheFinish(msr->pst, NULL, 0, NULL, NULL);
    msr->iCurrMaxRung = out.nMaxRung;
    }

uint64_t msrMaxOrder(MSR msr) {
    return msr->nMaxOrder;
    }

void msrGetNParts(MSR msr) { /* JW: Not pretty -- may be better way via fio */
    struct outGetNParts outget;

    pstGetNParts(msr->pst,NULL,0,&outget,NULL);
    assert(outget.nGas == msr->nGas);
    assert(outget.nDark == msr->nDark);
    assert(outget.nStar == msr->nStar);
    msr->nMaxOrder = outget.nMaxOrder;
#if 0
    if (outget.iMaxOrderGas > msr->nMaxOrder) {
	msr->nMaxOrder = outget.iMaxOrderGas;
	fprintf(stderr,"WARNING: Largest iOrder of gas > Largest iOrder of star\n");
	}
    if (outget.iMaxOrderDark > msr->nMaxOrder) {
	msr->nMaxOrder = outget.iMaxOrderDark;
	fprintf(stderr,"WARNING: Largest iOrder of dark > Largest iOrder of star\n");
	}
#endif
    }

void
msrAddDelParticles(MSR msr) {
    struct outColNParts *pColNParts;
    uint64_t *pNewOrder;
    struct inSetNParts in;
    int iOut;
    int i;

    msrprintf(msr,"Changing Particle number\n");
    pColNParts = malloc(msr->nThreads*sizeof(*pColNParts));
    assert(pColNParts!=NULL);
    pstColNParts(msr->pst, NULL, 0, pColNParts, &iOut);
    /*
     * Assign starting numbers for new particles in each processor.
     */
    pNewOrder = malloc(msr->nThreads*sizeof(*pNewOrder));
    assert(pNewOrder!=NULL);
    for (i=0;i<msr->nThreads;i++) {
	/*
	 * Detect any changes in particle number, and force a tree
	 * build.
	 */
	if (pColNParts[i].nNew != 0 || pColNParts[i].nDeltaGas != 0 ||
		pColNParts[i].nDeltaDark != 0 || pColNParts[i].nDeltaStar != 0) {
	    /*printf("Particle assignments have changed!\n");
	      printf("need to rebuild tree, code in msrAddDelParticles()\n");
	      printf("needs to be updated. Bailing out for now...\n");
	      exit(-1); */
	    pNewOrder[i] = msr->nMaxOrder+1; /* JW: +1 was missing for some reason */
	    msr->nMaxOrder += pColNParts[i].nNew;
	    msr->nGas += pColNParts[i].nDeltaGas;
	    msr->nDark += pColNParts[i].nDeltaDark;
	    msr->nStar += pColNParts[i].nDeltaStar;
	    }
	}
    msr->N = msr->nGas + msr->nDark + msr->nStar;

    /*msr->nMaxOrderDark = msr->nMaxOrder;*/

    pstNewOrder(msr->pst,pNewOrder,(int)sizeof(*pNewOrder)*msr->nThreads,NULL,NULL);

    msrprintf(msr,"New numbers of particles: %"PRIu64" gas %"PRIu64" dark %"PRIu64" star\n",
	      msr->nGas, msr->nDark, msr->nStar);

    in.nGas = msr->nGas;
    in.nDark = msr->nDark;
    in.nStar = msr->nStar;
    pstSetNParts(msr->pst,&in,sizeof(in),NULL,NULL);

    free(pNewOrder);
    free(pColNParts);
    }

int msrDoDensity(MSR msr) {
    return(msr->param.bDoDensity);
    }

#ifdef USE_PNG
int msrPNGResolution(MSR msr) {
    return(msr->param.nPNGResolution);
    }
#endif

int msrDoGravity(MSR msr) {
    return(msr->param.bDoGravity);
    }

/* Gas routines */

int msrDoGas(MSR msr) {
    return(msr->param.bDoGas);
    }

int msrMeshlessHydro(MSR msr) {
    return(msr->param.bMeshlessHydro);
    }

double msrNeighborsStd(MSR msr) {
    return(msr->param.dNeighborsStd);
    }
void msrIncreaseNeighborsStd(MSR msr) {
    msr->param.dNeighborsStd *= 1.4;
    if (msr->param.dNeighborsStd > 0.25*msr->param.nSmooth) msr->param.dNeighborsStd = 0.25*msr->param.nSmooth;
    }
void msrResetNeighborsStd(MSR msr) {
    msr->param.dNeighborsStd  = msr->param.dNeighborsStd0;
    }

int msrFirstHydroLoop(MSR msr) {
    return(msr->param.bFirstHydroLoop);
    }
void msrSetFirstHydroLoop(MSR msr, int value) {
    msr->param.bFirstHydroLoop = value;
    }

void msrInitSph(MSR msr,double dTime)
    {
    if (!msrMeshlessHydro(msr)){
    /* Init gas, internal energy -- correct estimate from dTuFac */
    msrActiveRung(msr,0,1);
    /* Very important NOT to do this if starting from a checkpoint */
    if (msr->param.bGasCooling && !msr->param.bRestart) {
	struct inCorrectEnergy in;
	double a;

 	//msrSelSrcGas(msr); /* Not really sure what the setting here needs to be */
	//msrSelDstGas(msr);  
	msrSmooth(msr,dTime,SMX_DENDVDX,0,msr->param.nSmooth);  
	//msrSelSrcAll(msr);
	//msrSelDstAll(msr);

	in.dTuFac = msr->param.dTuFac;
	a = csmTime2Exp(msr->param.csm,dTime);
	in.z = 1/a - 1;
	in.dTime = dTime;
	if (msr->param.bInitTFromCooling) {
	    fprintf(stderr,"INFO: Resetting thermal energies to special value in cooling routines\n");
	    in.iDirection = CORRECTENERGY_SPECIAL;
	    }
	else {
	    fprintf(stderr,"INFO: Correcting (dTuFac) thermal energies using cooling routines\n");
	    in.iDirection = CORRECTENERGY_IN;
	    }
	   
	pstCorrectEnergy(msr->pst, &in, sizeof(in), NULL, NULL);
	}
    
    /* Init forces, ... */
    msrActiveRung(msr,0,1);
    msrSph(msr,dTime,0);
    msrSphStep(msr,0,MAX_RUNG,dTime); /* Requires SPH */
    msrCooling(msr,dTime,0,0,1,1); /* Interate cooling for consistent dt */
    }else{ //IA: we set the initial rungs of the particles
        msrActiveRung(msr,0,1);
      if (msrDoGravity(msr)){ //IA: We need this for the acceleration time step criteria in the kepler ring!
//          msrApplyGravWork(msr, -1, 0.0, 0, MAX_RUNG);  
      }
        msrSetSmooth(msr);
        msrUpdatePrimVars(msr, dTime, 0.0, ROOT);
        msrMeshlessGradients(msr, dTime, 0.0, ROOT);
        msrMeshlessFluxes(msr, dTime, 0.0, ROOT);
	msrZeroNewRung(msr,0,MAX_RUNG,0); 
        msrHydroStep(msr,0,MAX_RUNG,dTime); // We do this twice because we need to have uNewRung for the time limiter
        msrHydroStep(msr,0,MAX_RUNG,dTime);  // of Durier & Dalla Vecchia
        msrUpdateConsVars(msr, dTime, 0.0, ROOT); // Reset the fluxes
    }
}

void msrSph(MSR msr,double dTime, double dStep) {
    double sec,dsec;

    if (msr->param.bVStep) printf("Calculating Sph, Step:%f\n",dStep);
    sec = msrTime();

/* JW: Is the tree aware of this -- does it need to be? 
       Will smooth behave correctly? */

    //msrSelSrcGas(msr); /* Not really sure what the setting here needs to be */
    //msrSelDstGas(msr);  
    msrSmooth(msr,dTime,SMX_DENDVDX,0,msr->param.nSmooth);  
    msrSmooth(msr,dTime,SMX_SPHFORCES,1,msr->param.nSmooth); /* Should be a resmooth */
    //msrSelSrcAll(msr);
    //msrSelDstAll(msr);

    dsec = msrTime() - sec;
    if (msr->param.bVStep) {
	printf("SPH Calculated, Wallclock: %f secs\n",  dsec);
	}
    }

void msrCoolSetup(MSR msr, double dTime) {
    }

void msrCooling(MSR msr,double dTime,double dStep,int bUpdateState, int bUpdateTable, int bIterateDt) {
    }
#ifdef COOLING
void msrCoolingUpdate(MSR msr,float redshift) {
   printf("Updating cooling.. \n");
   struct inCoolUpdate in;
   in.redshift = redshift;
   pstCoolingUpdate(msr->pst, &in, sizeof(in), NULL, NULL);
    }
void msrCoolingInit(MSR msr) {
    pstCoolingInit(msr->pst,NULL,0,NULL,NULL);
    }
#endif

/* END Gas routines */

void msrHopWrite(MSR msr, const char *fname) {
    FILE *fp;
    LCL *plcl;
    PST pst0;
    int id;
    double sec,dsec;

    pst0 = msr->pst;
    while (pst0->nLeaves > 1)
	pst0 = pst0->pstLower;
    plcl = pst0->plcl;


    if (msr->param.bVStep)
	printf("Writing group statistics to %s\n", fname );
    sec = msrTime();

    /* This is the new parallel binary format */
    struct inOutput out;
    out.eOutputType = OUT_TINY_GROUP;
    out.iPartner = -1;
    out.nPartner = -1;
    out.iProcessor = 0;
    out.nProcessor = msr->param.bParaWrite==0?1:(msr->param.nParaWrite<=1 ? msr->nThreads:msr->param.nParaWrite);
    strcpy(out.achOutFile,fname);
    pstOutput(msr->pst,&out,sizeof(out),NULL,NULL);
    dsec = msrTime() - sec;
    if (msr->param.bVStep)
	printf("Written statistics, Wallclock: %f secs\n",dsec);

    }

void msrHop(MSR msr, double dTime) {
    struct inSmooth in;
    struct inHopLink h;
    struct outHopJoin j;
    struct inHopFinishUp inFinish;
    struct inHopTreeBuild inTreeBuild;
    struct inHopGravity inGravity;
    struct inHopUnbind inUnbind;
    struct outHopUnbind outUnbind;
    struct inGroupStats inGroupStats;
    int i;
    uint64_t nGroups;
    double sec,dsec,ssec;

    ssec = msrTime();

    h.nSmooth    = in.nSmooth = 80;
    h.bPeriodic  = in.bPeriodic = msr->param.bPeriodic;
    h.bSymmetric = in.bSymmetric = 0;
    h.dHopTau    = msr->param.dHopTau<0 ? msr->param.dHopTau : msr->param.dHopTau;
    h.smf.a      = in.smf.a = dTime;
    h.smf.dTau2  = in.smf.dTau2 = 0.0;
    h.smf.nMinMembers = in.smf.nMinMembers = msr->param.nMinMembers;
    msrSmoothSetSMF(msr, &(in.smf), dTime);
    msrSmoothSetSMF(msr, &(h.smf), dTime);

    if (msr->param.bVStep) {
	if (h.dHopTau<0.0)
	    printf("Running Grasshopper with adaptive linking length (%g times softening)\n", -h.dHopTau );
	else
	    printf("Running Grasshopper with fixed linking length %g\n", h.dHopTau );
	}

    in.iSmoothType = SMX_DENSITY_M3;
    sec = msrTime();
    pstSmooth(msr->pst,&in,sizeof(in),NULL,NULL);
    dsec = msrTime() - sec;
    if (msr->param.bVStep)
	printf("Density calculation complete in %f secs, finding chains...\n",dsec);

    h.iSmoothType = SMX_GRADIENT_M3;
    sec = msrTime();
    nGroups = 0;
    pstHopLink(msr->pst,&h,sizeof(h),&nGroups,NULL);
    dsec = msrTime() - sec;
    if (msr->param.bVStep)
	printf("Chain search complete in %f secs, building minimal tree...\n",dsec);

    /* Build a new tree with only marked particles */
    sec = msrTime();
    msrBuildTreeMarked(msr,dTime);
    dsec = msrTime() - sec;
    if (msr->param.bVStep)
	printf("Tree build complete in %f secs, merging %"PRIu64" chains...\n",dsec,nGroups);

    h.iSmoothType = SMX_HOP_LINK;
    sec = msrTime();
    i = 0;
    do {
	++i;
	assert(i<100);
	pstHopJoin(msr->pst,&h,sizeof(h),&j,NULL);
	if (msr->param.bVStep)
	    printf("... %d iteration%s, %"PRIu64" chains remain\n",i,i==1?"":"s",j.nGroups);
	} while( !j.bDone );
    nGroups = j.nGroups;
    dsec = msrTime() - sec;
    if (msr->param.bVStep)
	printf("Chain merge complete in %f secs, %"PRIu64" groups\n",dsec,nGroups);
    inFinish.nMinGroupSize = msr->param.nMinMembers;
    inFinish.bPeriodic = msr->param.bPeriodic;
    inFinish.fPeriod[0] = msr->param.dxPeriod;
    inFinish.fPeriod[1] = msr->param.dyPeriod;
    inFinish.fPeriod[2] = msr->param.dzPeriod;
    pstHopFinishUp(msr->pst,&inFinish,sizeof(inFinish),&nGroups,NULL);
    if (msr->param.bVStep)
	printf("Removed groups with fewer than %d particles, %"PRIu64" remain\n",
	    inFinish.nMinGroupSize, nGroups);
#if 0
    if (msr->param.bVStep)
	printf("Unbinding\n");

    inUnbind.dTime = dTime;
    inUnbind.bPeriodic = msr->param.bPeriodic;
    inUnbind.fPeriod[0] = msr->param.dxPeriod;
    inUnbind.fPeriod[1] = msr->param.dyPeriod;
    inUnbind.fPeriod[2] = msr->param.dzPeriod;
    inUnbind.nMinGroupSize = msr->param.nMinMembers;
    inUnbind.iIteration = 0;
    inGravity.dTime = dTime;
    inGravity.bPeriodic = msr->param.bPeriodic;
    inGravity.nGroup = msr->param.nGroup;
    inGravity.dEwCut = msr->param.dEwCut;
    inGravity.dEwhCut = msr->param.dEwhCut;
    inGravity.uRungLo = 0;
    inGravity.uRungHi = MAX_RUNG;
    inGravity.dThetaMin = msr->dThetaMin;

    inUnbind.iIteration=0;
    do {
	sec = msrTime();
	inTreeBuild.nBucket = msr->param.nBucket;
	inTreeBuild.nGroup = msr->param.nGroup;
	pstHopTreeBuild(msr->pst,&inTreeBuild,sizeof(inTreeBuild),NULL,NULL);
	dsec = msrTime() - sec;
	if (msr->param.bVStep)
	    printf("... group trees built, Wallclock: %f secs\n",dsec);

	sec = msrTime();
	pstHopGravity(msr->pst,&inGravity,sizeof(inGravity),NULL,NULL);
	dsec = msrTime() - sec;
	if (msr->param.bVStep)
	    printf("... gravity complete, Wallclock: %f secs\n",dsec);
	
	sec = msrTime();
	pstHopUnbind(msr->pst,&inUnbind,sizeof(inUnbind),&outUnbind,NULL);
	nGroups = outUnbind.nGroups;
	dsec = msrTime() - sec;
	if (msr->param.bVStep)
	    printf("Unbinding completed in %f secs, %"PRIu64" particles evaporated, %"PRIu64" groups remain\n",
		dsec,outUnbind.nEvaporated, nGroups);
	} while(++inUnbind.iIteration < 100 && outUnbind.nEvaporated);
#endif
    /*
    ** This should be done as a separate msr function.
    */
    inGroupStats.bPeriodic = msr->param.bPeriodic;
    inGroupStats.dPeriod[0] = msr->param.dxPeriod;
    inGroupStats.dPeriod[1] = msr->param.dyPeriod;
    inGroupStats.dPeriod[2] = msr->param.dzPeriod;
    inGroupStats.rEnvironment[0] = msr->param.dEnvironment0;
    inGroupStats.rEnvironment[1] = msr->param.dEnvironment1;
    if ( prmSpecified(msr->prm,"dBoxSize") && msr->param.dBoxSize > 0.0 ) {
	inGroupStats.rEnvironment[0] /= msr->param.dBoxSize;
	inGroupStats.rEnvironment[1] /= msr->param.dBoxSize;
	}
    pstGroupStats(msr->pst,&inGroupStats,sizeof(inGroupStats),NULL,NULL);

    dsec = msrTime() - ssec;
    if (msr->param.bVStep)
	printf("Grasshopper complete, Wallclock: %f secs\n\n",dsec);
    }


void msrNewFof(MSR msr, double dTime) {
    struct inNewFof in;
    struct outFofPhases out;
    struct inFofFinishUp inFinish;
    int i;
    uint64_t nGroups;
    double sec,dsec,ssec;

    ssec = msrTime();

    in.dTau2 = msr->param.dTau*msr->param.dTau;
    in.nMinMembers = msr->param.nMinMembers;
    if (msr->param.bVStep) {
	printf("Running FoF with fixed linking length %g\n", msr->param.dTau );
	}

    sec = msrTime();
    pstNewFof(msr->pst,&in,sizeof(in),NULL,NULL);
    dsec = msrTime() - sec;
    if (msr->param.bVStep)
	printf("Initial FoF calculation complete in %f secs\n",dsec);

    sec = msrTime();
    i = 0;
    do {
	++i;
	assert(i<100);
	pstFofPhases(msr->pst,NULL,0,&out,NULL);
	if (msr->param.bVStep)
	    printf("... %d iteration%s\n",i,i==1?"":"s");
	} while( out.bMadeProgress );
    dsec = msrTime() - sec;
    if (msr->param.bVStep)
	printf("Global merge complete in %f secs\n",dsec);

    inFinish.nMinGroupSize = msr->param.nMinMembers;
    pstFofFinishUp(msr->pst,&inFinish,sizeof(inFinish),&nGroups,NULL);
    if (msr->param.bVStep)
	printf("Removed groups with fewer than %d particles, %"PRIu64" remain\n",
	    inFinish.nMinGroupSize, nGroups);
//    pstGroupRelocate(msr->pst,NULL,0,NULL,NULL);
    dsec = msrTime() - ssec;
    if (msr->param.bVStep)
	printf("FoF complete, Wallclock: %f secs\n",dsec);
    }


void msrGroupStats(MSR msr) {
    struct inGroupStats inGroupStats;
    double sec,dsec,ssec;

    if (msr->param.bVStep)
	printf("Generating Group statistics\n");
    sec = msrTime();
    inGroupStats.bPeriodic = msr->param.bPeriodic;
    inGroupStats.dPeriod[0] = msr->param.dxPeriod;
    inGroupStats.dPeriod[1] = msr->param.dyPeriod;
    inGroupStats.dPeriod[2] = msr->param.dzPeriod;
    inGroupStats.rEnvironment[0] = msr->param.dEnvironment0;
    inGroupStats.rEnvironment[1] = msr->param.dEnvironment1;
    if ( prmSpecified(msr->prm,"dBoxSize") && msr->param.dBoxSize > 0.0 ) {
	inGroupStats.rEnvironment[0] /= msr->param.dBoxSize;
	inGroupStats.rEnvironment[1] /= msr->param.dBoxSize;
	}
    pstGroupStats(msr->pst,&inGroupStats,sizeof(inGroupStats),NULL,NULL);
    dsec = msrTime() - sec;
    if (msr->param.bVStep)
	printf("Group statistics complete, Wallclock: %f secs\n\n",dsec);
    }

void msrDeleteGroups(MSR msr) {

    LCL *plcl;
    PST pst0;

    pst0 = msr->pst;
    while (pst0->nLeaves > 1)
	pst0 = pst0->pstLower;
    plcl = pst0->plcl;

    if (plcl->pkd->groupBin)free(plcl->pkd->groupBin);
    plcl->pkd->nBins = 0;
    plcl->pkd->nGroups = 0;
    }

void msrInitRelaxation(MSR msr) {
    pstInitRelaxation(msr->pst,NULL,0,NULL,NULL);
    }

void msrRelaxation(MSR msr,double dTime,double deltaT,int iSmoothType,int bSymmetric) {
    struct inSmooth in;
    in.nSmooth = msr->param.nSmooth;
    in.bPeriodic = msr->param.bPeriodic;
    in.bSymmetric = bSymmetric;
    in.iSmoothType = iSmoothType;
#if 0
    in.dfBall2OverSoft2 = (msr->param.bLowerSoundSpeed ? 0 :
			   4.0*msr->param.dhMinOverSoft*msr->param.dhMinOverSoft);
#endif
    if (msrComove(msr)) {
	in.smf.H = csmTime2Hub(msr->param.csm,dTime);
	in.smf.a = csmTime2Exp(msr->param.csm,dTime);
	}
    else {
	in.smf.H = 0.0;
	in.smf.a = 1.0;
	}
    in.smf.dDeltaT = deltaT;
    if (msr->param.bVStep) {
	double sec,dsec;
	printf("Smoothing for relaxation...dDeltaT = %f \n",deltaT);
	sec = msrTime();
	pstSmooth(msr->pst,&in,sizeof(in),NULL,NULL);
	dsec = msrTime() - sec;
	printf("Relaxation Calculated, Wallclock: %f secs\n\n",dsec);
	}
    else {
	pstSmooth(msr->pst,&in,sizeof(in),NULL,NULL);
	}
    }

#ifdef MDL_FFTW
double msrGenerateIC(MSR msr) {
    static/*FIXME: this is a HACK: message is too large for stack */ struct inGenerateIC in;
    struct outGenerateIC out;
    struct inGetFFTMaxSizes inFFTSizes;
    struct outGetFFTMaxSizes outFFTSizes;
    uint64_t nSpecies[FIO_SPECIES_LAST];
    int nOut;
    double sec,dsec;
    double dvFac;
    double mean, rms;
    uint64_t nTotal;
    int j;

    in.dBoxSize = msr->param.dBoxSize;
    in.iSeed = msr->param.iSeed;
    in.bFixed = msr->param.bFixedAmpIC;
    in.fPhase = msr->param.dFixedAmpPhasePI * M_PI;
    in.nGrid = msr->param.nGrid;
    in.b2LPT = msr->param.b2LPT;
    in.bClass = msr->param.csm->val.classData.bClass;
    in.cosmo = msr->param.csm->val;
    in.nInflateFactor = msr->param.nInflateReps + 1;
    in.nInflateFactor *= in.nInflateFactor * in.nInflateFactor;

    nTotal  = in.nGrid; /* Careful: 32 bit integer cubed => 64 bit integer */
    nTotal *= in.nGrid;
    nTotal *= in.nGrid;
    in.dBoxMass = msr->param.csm->val.dOmega0 / nTotal;

    for( j=0; j<FIO_SPECIES_LAST; j++) nSpecies[j] = 0;
    nSpecies[FIO_SPECIES_ALL] = nSpecies[FIO_SPECIES_DARK] = nTotal;
    msrInitializePStore(msr,nSpecies);

    if (prmSpecified(msr->prm,"dRedFrom")) {
	assert(msr->param.dRedFrom >= 0.0 );
	in.dExpansion = 1.0 / (1.0 + msr->param.dRedFrom);
	}
    else in.dExpansion = 0.0;

    msr->N     = nSpecies[FIO_SPECIES_ALL];
    msr->nGas  = nSpecies[FIO_SPECIES_SPH];
    msr->nDark = nSpecies[FIO_SPECIES_DARK];
    msr->nStar = nSpecies[FIO_SPECIES_STAR];
    msr->nMaxOrder = msr->N;

    if (msr->param.bVStart)
	printf("Generating IC...\nN:%"PRIu64" nDark:%"PRIu64
	       " nGas:%"PRIu64" nStar:%"PRIu64"\n",
	       msr->N, msr->nDark,msr->nGas,msr->nStar);

    /* Read the transfer function */
    in.nTf = 0;
    if (prmSpecified(msr->prm,"achTfFile")) {
	FILE *fp = fopen(msr->param.achTfFile,"r");
	char buffer[256];

	if (msr->param.bVStart)
	    printf("Reading transfer function from %s\n", msr->param.achTfFile);
	if (fp == NULL) {
	    perror(msr->param.achTfFile);
	    _msrExit(msr,1);
	    }
	while(fgets(buffer,sizeof(buffer),fp)) {
	    assert(in.nTf < MAX_TF);
	    if (sscanf(buffer," %lg %lg\n",&in.k[in.nTf],&in.tf[in.nTf])==2) {
		in.k[in.nTf] = log(in.k[in.nTf]);
		in.tf[in.nTf] = log(in.tf[in.nTf]);
		++in.nTf;
		}
	    }
	fclose(fp);
	if (msr->param.bVStart)
	    printf("Transfer function : %d lines kmin %g kmax %g\n",
		in.nTf, exp(in.k[0]), exp(in.k[in.nTf-1]));

	}

    sec = msrTime();

    /* Figure out the minimum number of particles */
    inFFTSizes.nx = inFFTSizes.ny = inFFTSizes.nz = in.nGrid;
    pstGetFFTMaxSizes(msr->pst,&inFFTSizes,sizeof(inFFTSizes),&outFFTSizes,&nOut);
    printf("Grid size %d x %d x %d, per node %d x %d x %d and %d x %d x %d\n",
	inFFTSizes.nx, inFFTSizes.ny, inFFTSizes.nz,
	inFFTSizes.nx, inFFTSizes.ny, outFFTSizes.nMaxZ,
	inFFTSizes.nx, outFFTSizes.nMaxY, inFFTSizes.nz);

    msrprintf(msr,"IC Generation @ a=%g with seed %d\n",in.dExpansion,msr->param.iSeed);
    in.nPerNode = outFFTSizes.nMaxLocal;
    pstGenerateIC(msr->pst,&in,sizeof(in),&out,&nOut);
    mean = 2*out.noiseMean / msr->N;
    rms = sqrt(2*out.noiseCSQ / msr->N);

    msrprintf(msr,"Transferring particles between/within nodes\n");
    pstMoveIC(msr->pst,&in,sizeof(in),NULL,NULL);

    msrSetClasses(msr);
    dsec = msrTime() - sec;
    PKD pkd = msr->pst->plcl->pkd;
    msrprintf(msr,"IC Generation Complete @ a=%g, Wallclock: %f secs\n\n",out.dExpansion,dsec);
    msrprintf(msr,"Mean of noise same is %g, RMS %g.\n",mean,rms);

    return getTime(msr,out.dExpansion,&dvFac);
    }
#endif

double msrRead(MSR msr, const char *achInFile) {
    double dTime,dExpansion;
    FIO fio;
    int j;
    double sec,dsec;
    struct inReadFile *read;
    uint64_t nSpecies[FIO_SPECIES_LAST];
    size_t nBytes;
    inReadFileFilename achFilename;
    uint64_t mMemoryModel = 0;
    LCL *plcl;
    PST pst0;

    pst0 = msr->pst;
    while (pst0->nLeaves > 1)
	pst0 = pst0->pstLower;
    plcl = pst0->plcl;

    mMemoryModel = getMemoryModel(msr);

    sec = msrTime();

    nBytes = PST_MAX_FILES*(sizeof(fioSpeciesList)+PST_FILENAME_SIZE);
    read = malloc(sizeof(struct inReadFile) + nBytes);
    assert(read != NULL);

    /* Add Data Subpath for local and non-local names. */
    _msrMakePath(msr->param.achDataSubPath,achInFile,achFilename);
//    strcpy(read.achFilename,achFilename);
    fio = fioOpen(achFilename,msr->param.csm->val.dOmega0,msr->param.csm->val.dOmegab);
    if (fio==NULL) {
	fprintf(stderr,"ERROR: unable to open input file\n");
	perror(achFilename);
	_msrExit(msr,1);
	}
    nBytes = fioDump(fio,nBytes,read+1);

    if (!fioGetAttr(fio,"dTime",FIO_TYPE_DOUBLE,&dExpansion)) dExpansion = 0.0;
    if (!fioGetAttr(fio,"dEcosmo",FIO_TYPE_DOUBLE,&msr->dEcosmo)) msr->dEcosmo = 0.0;
    if (!fioGetAttr(fio,"dTimeOld",FIO_TYPE_DOUBLE,&msr->dTimeOld)) msr->dTimeOld = 0.0;
    if (!fioGetAttr(fio,"dUOld",FIO_TYPE_DOUBLE,&msr->dUOld)) msr->dUOld = 0.0;

    if(!prmSpecified(msr->prm, "dOmega0")) fioGetAttr(fio,"dOmega0",FIO_TYPE_DOUBLE,&msr->param.csm->val.dOmega0);
    if(!prmSpecified(msr->prm, "dLambda")) fioGetAttr(fio,"dLambda",FIO_TYPE_DOUBLE,&msr->param.csm->val.dLambda);
    if(!prmSpecified(msr->prm, "dBoxSize")) fioGetAttr(fio,"dBoxSize",FIO_TYPE_DOUBLE,&msr->param.dBoxSize);
    if(!prmSpecified(msr->prm, "h")) fioGetAttr(fio,"h",FIO_TYPE_DOUBLE,&msr->param.h);

    msr->N     = fioGetN(fio,FIO_SPECIES_ALL);
    msr->nGas  = fioGetN(fio,FIO_SPECIES_SPH);
    msr->nDark = fioGetN(fio,FIO_SPECIES_DARK);
    msr->nStar = fioGetN(fio,FIO_SPECIES_STAR);
    msr->nMaxOrder = msr->N;

    read->nProcessors = msr->param.bParaRead==0?1:(msr->param.nParaRead<=1 ? msr->nThreads:msr->param.nParaRead);

    if (!fioGetAttr(fio,"nFiles",FIO_TYPE_UINT32,&j)) j = 1;
    printf("Reading %"PRIu64" particles from %d file%s using %d processor%s\n",
	msr->N, j, (j==1?"":"s"), read->nProcessors, (read->nProcessors==1?"":"s") );

    dTime = getTime(msr,dExpansion,&read->dvFac);
    if (msr->param.bInFileLC) read->dvFac = 1.0;
    read->dTuFac = msr->param.dTuFac;
    
    if (msr->nGas && !prmSpecified(msr->prm,"bDoGas")) msr->param.bDoGas = 1;
    if (msrDoGas(msr) || msr->nGas) mMemoryModel |= (PKD_MODEL_SPH|PKD_MODEL_ACCELERATION|PKD_MODEL_VELOCITY|PKD_MODEL_NODE_SPHBNDS);		
    if (msr->param.bStarForm || msr->nStar) mMemoryModel |= (PKD_MODEL_SPH|PKD_MODEL_ACCELERATION|PKD_MODEL_VELOCITY|PKD_MODEL_MASS|PKD_MODEL_SOFTENING|PKD_MODEL_STAR);
    
    read->nNodeStart = 0;
    read->nNodeEnd = msr->N - 1;


    for( j=0; j<FIO_SPECIES_LAST; j++) nSpecies[j] = fioGetN(fio,j);
    msrInitializePStore(msr, nSpecies);

    read->dOmega0 = msr->param.csm->val.dOmega0;
    read->dOmegab = msr->param.csm->val.dOmegab;

    /*
    ** If bParaRead is 0, then we read serially; if it is 1, then we read
    ** in parallel using all available threads, otherwise we read in parallel
    ** using the specified number of threads.  The latter option will reduce
    ** the total amount of simultaneous I/O for file systems that cannot
    ** handle it.
    */

    if (msr->param.bParaRead) {
	fioClose(fio);
	pstReadFile(msr->pst,read,sizeof(struct inReadFile)+nBytes,NULL,NULL);
	}
    else {
	msrOneNodeRead(msr,read,fio);
	fioClose(fio);
	}

    dsec = msrTime() - sec;
    msrSetClasses(msr);
    printf("Input file has been successfully read, Wallclock: %f secs.\n", dsec);

    free(read);

    /*
    ** If this is a non-periodic box, then we must precalculate the bounds.
    ** We throw away the result, but PKD will keep track for later.
    */
    if (!msr->param.bPeriodic ||
	    msr->param.dxPeriod >= FLOAT_MAXVAL ||
	    msr->param.dyPeriod >= FLOAT_MAXVAL ||
	    msr->param.dzPeriod >= FLOAT_MAXVAL) {
	BND bnd;
	msrCalcBound(msr,&bnd);
	}

    /*
    ** Now read in the output points, passing the initial time.
    ** We do this only if nSteps is not equal to zero.
    */
    if (msrSteps(msr) > 0) msrReadOuts(msr,dTime);

    return dTime;
    }


void msrCalcBound(MSR msr,BND *pbnd) {
    /*
    ** This sets the local pkd->bnd.
    */
    pstCalcBound(msr->pst,NULL,0,pbnd,NULL);
    }

void msrCalcVBound(MSR msr,BND *pbnd) {
    /*
    ** This sets the local pkd->bnd.
    */
    pstCalcVBound(msr->pst,NULL,0,pbnd,NULL);
    }

#ifdef MDL_FFTW
void msrOutputPk(MSR msr,int iStep,double dTime) {
#ifdef _MSC_VER
    char achFile[MAX_PATH];
#else
    char achFile[PATH_MAX];
#endif
    float *fK, *fPk;
    uint64_t *nPk;
    double a, vfact, kfact;
    FILE *fp;
    int i;

    if (msr->param.nGridPk == 0) return;

    msrGridCreateFFT(msr,msr->param.nGridPk);


    fK = malloc(sizeof(float)*(msr->param.nBinsPk));
    assert(fK != NULL);
    fPk = malloc(sizeof(float)*(msr->param.nBinsPk));
    assert(fPk != NULL);
    nPk = malloc(sizeof(uint64_t)*(msr->param.nBinsPk));
    assert(nPk != NULL);

    msrMeasurePk(msr,msr->param.iPkOrder,msr->param.bPkInterlace,msr->param.nGridPk,msr->param.nBinsPk,nPk,fK,fPk);

    msrBuildName(msr,achFile,iStep);
    strncat(achFile,".pk",256);

    if (!msr->param.csm->val.bComove) a = 1.0;
    else a = csmTime2Exp(msr->param.csm,dTime);

    /* If the Box Size (in mpc/h) was specified, then we can scale the output power spectrum measurement */
    if ( prmSpecified(msr->prm,"dBoxSize") && msr->param.dBoxSize > 0.0 ) kfact = msr->param.dBoxSize;
    else kfact = 1.0;
    vfact = kfact * kfact * kfact;
    kfact = 1.0 / kfact;

    fp = fopen(achFile,"w");
    if ( fp==NULL) {
	printf("Could not create P(k) File:%s\n",achFile);
	_msrExit(msr,1);
	}
    for(i=0; i<msr->param.nBinsPk; ++i) {
	if (fPk[i] > 0.0) fprintf(fp,"%g %g %" PRIu64 "\n",
 	    kfact * fK[i] * 2.0 * M_PI,vfact * fPk[i], nPk[i]);
	}
    fclose(fp);
    free(fK);
    free(fPk);
    free(nPk);
    /* Output the k-grid if requested */
    if (msr->param.iDeltakInterval && (iStep % msr->param.iDeltakInterval == 0)) {
	struct inOutput out;
	double dsec, sec = msrTime();
	out.eOutputType = OUT_KGRID;
	out.iPartner = -1;
	out.nPartner = -1;
	out.iProcessor = 0;
	out.nProcessor = msr->param.bParaWrite==0?1:(msr->param.nParaWrite<=1 ? msr->nThreads:msr->param.nParaWrite);
	if (out.nProcessor > mdlProcs(msr->mdl)) out.nProcessor = mdlProcs(msr->mdl); 
	msrBuildName(msr,out.achOutFile,iStep);
	strncat(out.achOutFile,".deltak",256);
        printf("Writing Delta(k) to %s ...\n",out.achOutFile);
	pstOutput(msr->pst,&out,sizeof(out),NULL,NULL);
	dsec = msrTime() - sec;
	msrprintf(msr,"Delta(k) has been successfully written, Wallclock: %f secs.\n\n", dsec);
	}
    msrGridDeleteFFT(msr);
    }


void msrOutputLinPk(MSR msr,int iStep,double dTime) {
#ifdef _MSC_VER
    char achFile[MAX_PATH];
#else
    char achFile[PATH_MAX];
#endif
    float *fK, *fPk;
    uint64_t *nPk;
    double a, vfact, kfact;
    FILE *fp;
    int i;

    if (msr->param.nGridLin == 0) return;
    if (!msr->param.csm->val.bComove) return;
    if (!prmSpecified(msr->prm, "dBoxSize")) return; 
    fK = malloc(sizeof(float)*(msr->param.nBinsLinPk));
    assert(fK != NULL);
    fPk = malloc(sizeof(float)*(msr->param.nBinsLinPk));
    assert(fPk != NULL);
    nPk = malloc(sizeof(uint64_t)*(msr->param.nBinsLinPk));
    assert(nPk != NULL);

    a = csmTime2Exp(msr->param.csm, dTime);

    msrMeasureLinPk(msr,msr->param.nGridLin,a,msr->param.dBoxSize,nPk,fK,fPk);

    msrBuildName(msr,achFile,iStep);
    strncat(achFile,".lin_pk",256);

    if (!msr->param.csm->val.bComove) a = 1.0;
    else a = csmTime2Exp(msr->param.csm,dTime);

    if ( msr->param.dBoxSize > 0.0 ) kfact = msr->param.dBoxSize;
    else kfact = 1.0;
    vfact = kfact * kfact * kfact;
    kfact = 1.0 / kfact;

    fp = fopen(achFile,"w");
    if ( fp==NULL) {
	printf("Could not create P_lin(k) File:%s\n",achFile);
	_msrExit(msr,1);
	}
    for(i=0; i<msr->param.nBinsLinPk; ++i) {
	if (fPk[i] > 0.0) fprintf(fp,"%g %g %" PRIu64 "\n",
 	    kfact * fK[i] * 2.0 * M_PI,vfact * fPk[i], nPk[i]);
	}
    fclose(fp);
    free(fK);
    free(fPk);
    free(nPk);
    }

#endif

/*
**  This routine will output all requested files and fields
*/

void msrOutput(MSR msr, int iStep, double dTime, int bCheckpoint) {
#ifdef _MSC_VER
    char achFile[MAX_PATH];
#else
    char achFile[PATH_MAX];
#endif
    int bSymmetric;
    int nFOFsDone;
    int i,iSec=0;
    uint64_t nActive;

    printf( "Writing output for step %d\n", iStep );
    msrBuildIoName(msr,achFile,iStep);

    // IA: I want to do this even at step 0
    /*if ( iStep )*/ msrWrite(msr,achFile,dTime,bCheckpoint );

    if (msrDoGas(msr) && !msr->param.nSteps) {  /* Diagnostic Gas */ 
	msrReorder(msr);
	msrBuildName(msr,achFile,iStep);
	strncat(achFile,".uDot",256);
	msrOutArray(msr,achFile,OUT_UDOT_ARRAY);
	msrBuildName(msr,achFile,iStep);
	strncat(achFile,".u",256);
	msrOutArray(msr,achFile,OUT_U_ARRAY);
	msrBuildName(msr,achFile,iStep);
	strncat(achFile,".c",256);
	msrOutArray(msr,achFile,OUT_C_ARRAY);
	msrBuildName(msr,achFile,iStep);
	strncat(achFile,".hsph",256);
	msrOutArray(msr,achFile,OUT_HSPH_ARRAY);
	}

    if (msrDoDensity(msr)) {
#ifdef FAST_GAS
	msrActiveRung(msr,3,1); /* Activate some particles */
	msrDomainDecomp(msr,0,0,0);
	msrBuildTree(msr,dTime,0);

	//msrSelSrcGas(msr);  /* FOR TESTING!! of gas active particles */
	msrFastGasPhase1(msr,dTime,SMX_DENSITY);
	msrFastGasPhase2(msr,dTime,SMX_PRINTNN);
	//msrSelSrcAll(msr);  /* FOR TESTING!! of gas active particles */
#else
	msrActiveRung(msr,0,1); /* Activate all particles */
	msrDomainDecomp(msr,-1,0,0);
	msrBuildTree(msr,dTime,0);
	bSymmetric = 0;  /* should be set in param file! */
	msrSmooth(msr,dTime,SMX_DENSITY,bSymmetric,msr->param.nSmooth);
#endif
	}
    if ( msr->param.bFindGroups ) {
	msrReorder(msr);
	//sprintf(achFile,"%s.fof",msrOutName(msr));
	//msrOutArray(msr,achFile,OUT_GROUP_ARRAY);
	msrBuildName(msr,achFile,iStep);
	strncat(achFile,".fofstats",256);
	msrHopWrite(msr,achFile);
	}
    if ( msr->param.bFindHopGroups ) {
	msrActiveRung(msr,0,1); /* Activate all particles */
	msrDomainDecomp(msr,-1,0,0);
	msrBuildTree(msr,dTime,0);
	msrHop(msr,dTime);
	msrReorder(msr);

	msrBuildName(msr,achFile,iStep);
	//strncat(achFile,".hopgrp",256);
	//msrOutArray(msr,achFile,OUT_GROUP_ARRAY);

	msrBuildName(msr,achFile,iStep);
	strncat(achFile,".hopstats",256);
	msrHopWrite(msr,achFile);
	}

    if (msr->param.bDoAccOutput) {
	msrReorder(msr);
	msrBuildName(msr,achFile,iStep);
	strncat(achFile,".acc",256);
	msrOutVector(msr,achFile,OUT_ACCEL_VECTOR);
	}
    if (msr->param.bDoPotOutput) {
	msrReorder(msr);
	msrBuildName(msr,achFile,iStep);
	strncat(achFile,".pot",256);
	msrOutArray(msr,achFile,OUT_POT_ARRAY);
	}

    if ( msr->param.bTraceRelaxation) {
	msrReorder(msr);
	msrBuildName(msr,achFile,iStep);
	strncat(achFile,".relax",256);
	msrOutArray(msr,achFile,OUT_RELAX_ARRAY);
	}
    if ( msrDoDensity(msr) ) {
	msrReorder(msr);
	msrBuildName(msr,achFile,iStep);
	strncat(achFile,".den",256);
	msrOutArray(msr,achFile,OUT_DENSITY_ARRAY);
	}
    if (msr->param.bDoRungOutput) {
	msrReorder(msr);
	msrBuildName(msr,achFile,iStep);
	strncat(achFile,".rung",256);
	msrOutArray(msr,achFile,OUT_RUNG_ARRAY);
	}
    if (msr->param.bDoRungDestOutput) {
	msrReorder(msr);
	msrBuildName(msr,achFile,iStep);
	strncat(achFile,".rd",256);
	msrOutArray(msr,achFile,OUT_RUNGDEST_ARRAY);
	}
    if (msr->param.bDoSoftOutput) {
	msrReorder(msr);
	msrBuildName(msr,achFile,iStep);
	strncat(achFile,".soft",256);
	msrOutArray(msr,achFile,OUT_SOFT_ARRAY);
	}
    /*
    ** Don't allow duplicate outputs.
    */
    while (msrOutTime(msr,dTime));
    }

void msrSelAll(MSR msr) {
    pstSelAll(msr->pst, NULL, 0, NULL, NULL );
    }
void msrSelGas(MSR msr) {
    pstSelGas(msr->pst, NULL, 0, NULL, NULL );
    }
void msrSelStar(MSR msr) {
    pstSelStar(msr->pst, NULL, 0, NULL, NULL );
    }
void msrSelDeleted(MSR msr) {
    pstSelDeleted(msr->pst, NULL, 0, NULL, NULL );
    }
uint64_t msrSelById(MSR msr,uint64_t idStart,uint64_t idEnd,int setIfTrue,int clearIfFalse) {
    struct inSelById in;
    struct outSelById out;
    int nOut;

    in.idStart = idStart;
    in.idEnd = idEnd;
    in.setIfTrue = setIfTrue;
    in.clearIfFalse = clearIfFalse;
    pstSelById(msr->pst, &in, sizeof(in), &out, &nOut);
    return out.nSelected;
    }
uint64_t msrSelMass(MSR msr,double dMinMass,double dMaxMass,int setIfTrue,int clearIfFalse) {
    struct inSelMass in;
    struct outSelMass out;
    int nOut;

    in.dMinMass = dMinMass;
    in.dMaxMass = dMaxMass;
    in.setIfTrue = setIfTrue;
    in.clearIfFalse = clearIfFalse;
    pstSelMass(msr->pst, &in, sizeof(in), &out, &nOut);
    return out.nSelected;
    }
uint64_t msrSelPhaseDensity(MSR msr,double dMinPhaseDensity,double dMaxPhaseDensity,int setIfTrue,int clearIfFalse) {
    struct inSelPhaseDensity in;
    struct outSelPhaseDensity out;
    int nOut;

    in.dMinDensity = dMinPhaseDensity;
    in.dMaxDensity = dMaxPhaseDensity;
    in.setIfTrue = setIfTrue;
    in.clearIfFalse = clearIfFalse;
    pstSelPhaseDensity(msr->pst, &in, sizeof(in), &out, &nOut);
    return out.nSelected;
    }
uint64_t msrSelBox(MSR msr,double *dCenter, double *dSize,int setIfTrue,int clearIfFalse) {
    struct inSelBox in;
    struct outSelBox out;
    int nOut;

    in.dCenter[0] = dCenter[0];
    in.dCenter[1] = dCenter[1];
    in.dCenter[2] = dCenter[2];
    in.dSize[0] = dSize[0];
    in.dSize[1] = dSize[1];
    in.dSize[2] = dSize[2];
    in.setIfTrue = setIfTrue;
    in.clearIfFalse = clearIfFalse;
    pstSelBox(msr->pst, &in, sizeof(in), &out, &nOut);
    return out.nSelected;
    }
uint64_t msrSelSphere(MSR msr,double *r, double dRadius,int setIfTrue,int clearIfFalse) {
    struct inSelSphere in;
    struct outSelSphere out;
    int nOut;

    in.r[0] = r[0];
    in.r[1] = r[1];
    in.r[2] = r[2];
    in.dRadius = dRadius;
    in.setIfTrue = setIfTrue;
    in.clearIfFalse = clearIfFalse;
    pstSelSphere(msr->pst, &in, sizeof(in), &out, &nOut);
    return out.nSelected;
    }
uint64_t msrSelCylinder(MSR msr,double *dP1, double *dP2, double dRadius,
			   int setIfTrue, int clearIfFalse ) {
    struct inSelCylinder in;
    struct outSelCylinder out;
    int nOut,j;

    for(j=0;j<3;j++) {
	in.dP1[j] = dP1[j];
	in.dP2[j] = dP2[j];
	}
    in.dRadius = dRadius;
    in.setIfTrue = setIfTrue;
    in.clearIfFalse = clearIfFalse;
    pstSelCylinder(msr->pst, &in, sizeof(in), &out, &nOut);
    return out.nSelected;
    }

double msrTotalMass(MSR msr) {
    struct outTotalMass out;
    int nOut;

    pstTotalMass(msr->pst, NULL, 0, &out, &nOut);
    return out.dMass;
    }

void msrDeleteProfile(MSR msr) {
    LCL *plcl;
    PST pst0;

    pst0 = msr->pst;
    while (pst0->nLeaves > 1)
	pst0 = pst0->pstLower;
    plcl = pst0->plcl;

    if (plcl->pkd->profileBins) mdlFree(msr->mdl,plcl->pkd->profileBins);
    plcl->pkd->profileBins = NULL;
    }

void msrCalcDistance(MSR msr,const double *dCenter, double dRadius ) {
    struct inCalcDistance in;
    int j;

    for(j=0;j<3;j++) in.dCenter[j] = dCenter[j];
    in.dRadius = dRadius;
    pstCalcDistance(msr->pst, &in, sizeof(in), NULL, NULL);
    }

void msrCalcCOM(MSR msr,const double *dCenter, double dRadius,
		double *com, double *vcm, double *L, double *M) {
    struct inCalcCOM in;
    struct outCalcCOM out;
    int nOut;
    int j;
    double T[3];

    for(j=0;j<3;j++) in.dCenter[j] = dCenter[j];
    in.dRadius = dRadius;
    pstCalcCOM(msr->pst, &in, sizeof(in), &out, &nOut);
    assert( nOut == sizeof(out) );

    *M = out.M;
    if ( out.M > 0.0 ) {
	for( j=0; j<3; j++ ) {
	    com[j] = out.com[j] / out.M;
	    vcm[j] = out.vcm[j] / out.M;
	    }
	cross_product(T, com, vcm);
	vec_add_const_mult(L,out.L,-out.M,T);
	for( j=0; j<3; j++ ) L[j] /= out.M;
	}
    }

uint64_t msrCountDistance(MSR msr,double dRadius2Inner, double dRadius2Outer) {
    struct inCountDistance in;
    struct outCountDistance out;
    int nOut;
    in.dRadius2Inner = dRadius2Inner;
    in.dRadius2Outer = dRadius2Outer;
    pstCountDistance(msr->pst, &in, sizeof(in), &out, &nOut);
    assert( nOut == sizeof(out) );
    return out.nCount;
    }

typedef struct {
    double dFrac;       /* Fraction of particles in each bin */
    uint64_t nTotal;    /* Total number of particles in the range */
    uint64_t nInner;    /* Number inside minimum radius */
    uint64_t nTarget;   /* Target number of particles */
    uint64_t nSelected;
    MSR msr;
    } SPHERECTX;

static double countSphere(double r,void *vctx) {
    SPHERECTX *ctx = vctx;
    ctx->nSelected = msrCountDistance(ctx->msr,0.0,r*r);
    return 1.0*ctx->nSelected - 1.0*ctx->nTarget;
    }

static void profileRootFind( double *dBins, int lo, int hi, int nAccuracy, SPHERECTX *ctx ) {
    int nIter;
    int iBin = (lo+hi) / 2;
    if ( lo == iBin ) return;

    ctx->nTarget = d2u64((ctx->nTotal-ctx->nInner) * ctx->dFrac * iBin + ctx->nInner);
    dBins[iBin] = illinois( countSphere, ctx, dBins[lo], dBins[hi], 0.0, 1.0*nAccuracy, &nIter );
    profileRootFind(dBins,lo,iBin,nAccuracy,ctx);
    profileRootFind(dBins,iBin,hi,nAccuracy,ctx);
    }


typedef struct {
    double rMiddle;
    total_t nTarget;   /* Target number of particles */
    MSR msr;
    } SHELLCTX;

static double countShell(double rInner,void *vctx) {
    SHELLCTX *ctx = vctx;
    double rOuter;
    local_t nSelected;

    if ( rInner == ctx->rMiddle ) nSelected = 0;
    else {
	rOuter = pow(10,2.0*log10(ctx->rMiddle)-log10(rInner));
	nSelected = msrCountDistance(ctx->msr,rInner*rInner,rOuter*rOuter);
	}
    return 1.0*nSelected - 1.0*ctx->nTarget;
    }


/*
** Calculate a profile.
** Bins are of equal size (same number of particles) between dMinRadius and dLogRadius.
** From dLogRadius to dMaxRadius, the binning is done logarithmicly.
** Setting dLogRadius to dMinRadius results in purely logarithmic binning, while
** setting dLogRadius to dMaxRadius results in purely equal sized binning.
*/
void msrProfile( MSR msr, const PROFILEBIN **ppBins, int *pnBins,
		 double *r, double dMinRadius, double dLogRadius, double dMaxRadius,
		 int nPerBin, int nBins, int nAccuracy ) {
    SPHERECTX ctxSphere;
    SHELLCTX ctxShell;
    PROFILEBIN *pBins;
    double sec, dsec;
    double com[3], vcm[3], L[3], M;
    struct inProfile *in;
    size_t inSize;
    int i,j;
    int nBinsInner;
    total_t N,n;
    LCL *plcl;
    PST pst0;

    assert(dMinRadius<=dLogRadius);
    assert(dLogRadius<=dMaxRadius);
    assert(dLogRadius==dMinRadius || nPerBin>0);
    assert(dLogRadius==dMaxRadius || nBins>0);

    if ( dLogRadius == dMaxRadius ) nBins = 0;

    pst0 = msr->pst;
    while (pst0->nLeaves > 1)
        pst0 = pst0->pstLower;
    plcl = pst0->plcl;

    msrCalcDistance(msr,r,dMaxRadius);
    msrCalcCOM(msr,r,dMaxRadius,com,vcm,L,&M);

    if ( dLogRadius > dMinRadius ) {
	/*
	** The inner radius is calculated such that the logarithmic mid-point
	** falls on dMinRadius.  This is done so that the profile is plotted
	** all the way to the inner radius.  The correct radius must be between
	** dMinRadius and the logrithmic difference between dMinRadius and
	** dMaxRadius below dMinRadius.
	*/
	ctxShell.rMiddle = dMinRadius;
	ctxShell.nTarget = nPerBin;
	ctxShell.msr = msr;
	dMinRadius = illinois( countShell, &ctxShell,
			       pow(10,2.0*log10(dMinRadius)-log10(dMaxRadius)), dMinRadius,
			       0.0, 0.0, NULL );
	N = msrCountDistance(msr,dMinRadius*dMinRadius,dLogRadius*dLogRadius);
	nBinsInner = (N+nPerBin/2) / nPerBin;
	}
    else {
	double dOuter;

	nBinsInner = 0;

	/*
	** Calculate the logarithmic mid-point and verify that there are enough particles
	** in the first bin.  If not, invoke the root finder.
	*/
	ctxShell.rMiddle = dMinRadius;
	ctxShell.nTarget = nPerBin;
	ctxShell.msr = msr;
	dMinRadius = pow(10,(2.0*(nBins+1)*log10(dMinRadius)-log10(dMaxRadius))/(2*nBins));
	dOuter = pow(10,2.0*log10(ctxShell.rMiddle)-log10(dMinRadius));
	N = msrCountDistance(msr, dMinRadius*dMinRadius,dOuter*dOuter);
	if ( N < nPerBin-nAccuracy ) {
	    dMinRadius = illinois( countShell, &ctxShell,
				   pow(10,2.0*log10(dMinRadius)-log10(dMaxRadius)), dMinRadius,
				   0.0, 0.0, NULL );
	    }
	dLogRadius = dMinRadius;
	}

    inSize = sizeof(struct inProfile)-sizeof(in->dRadii[0])*(sizeof(in->dRadii)/sizeof(in->dRadii[0])-nBins-nBinsInner-1);
    in = malloc(inSize);
    assert(in!=NULL);

    in->dRadii[0] = dMinRadius;

    /*
    ** Inner, fixed size bins
    */
    if ( nBinsInner ) {
	sec = msrTime();
	msrprintf(msr, "Root finding for %d bins\n", nBinsInner );
	ctxSphere.nTotal = msrCountDistance(msr,0.0,dLogRadius*dLogRadius);
	ctxSphere.nInner = msrCountDistance(msr,0.0,dMinRadius*dMinRadius);
	ctxSphere.msr = msr;
	ctxSphere.dFrac = 1.0 / nBinsInner;
	in->dRadii[nBinsInner] = dLogRadius;
	profileRootFind( in->dRadii, 0, nBinsInner, nAccuracy, &ctxSphere );
	dsec = msrTime() - sec;
	msrprintf(msr,"Root finding complete, Wallclock: %f secs\n\n",dsec);
	}

    /*
    ** Now logarithmic binning for the outer region.  We still obey nPerBin
    ** as the minimum number of particles to include in each bin.
    */
    if ( nBins ) {
	double dLogMin;
	double dLogMax = log10(dMaxRadius);
	double dRadius;

	ctxSphere.nTotal = msrSelSphere(msr,r,dMaxRadius,1,1);
	ctxSphere.msr = msr;

	N = msrCountDistance(msr,0.0,dLogRadius*dLogRadius);
	for( i=1; i<nBins; i++ ) {
	    int nBinsRem = nBins - i + 1;

	    dLogMin = log10(in->dRadii[nBinsInner+i-1]);
	    dRadius = pow(10,(dLogMax-dLogMin)/nBinsRem + dLogMin);
	    n = msrCountDistance(msr,0.0,dRadius*dRadius);
	    if ( n-N < nPerBin-nAccuracy ) {
		ctxSphere.nTarget = N + nPerBin;
		dRadius = illinois( countSphere, &ctxSphere, 0.0, dMaxRadius,
				    0.0, 1.0*nAccuracy, NULL );
		n = ctxSphere.nSelected;
		}
	    in->dRadii[nBinsInner+i] = dRadius;
	    N = n;
	    }
	}

    nBins = nBins+nBinsInner;

    in->dRadii[nBins] = dMaxRadius;

    sec = msrTime();
    msrprintf( msr, "Profiling\n" );
    for(i=0; i<3; i++) {
	in->dCenter[i] = r[i];
	in->com[i] = com[i];
	in->vcm[i] = vcm[i];
	in->L[i] = L[i];
	}
    in->nBins = nBins+1;
    in->uRungLo = 0;
    in->uRungHi = msrMaxRung(msr)-1;
    pstProfile(msr->pst, in, inSize, NULL, NULL);
    free(in);

    /*
    ** Finalize bin values
    */
    pBins = plcl->pkd->profileBins;
    for( i=0; i<nBins+1; i++ ) {
	if ( pBins[i].dMassInBin > 0.0 ) {
	    pBins[i].vel_radial /= pBins[i].dMassInBin;
	    pBins[i].vel_radial_sigma /= pBins[i].dMassInBin;
	    pBins[i].vel_tang_sigma = sqrt(pBins[i].vel_tang_sigma / pBins[i].dMassInBin);
	    if (pBins[i].vel_radial_sigma > pBins[i].vel_radial*pBins[i].vel_radial)
		pBins[i].vel_radial_sigma = sqrt(pBins[i].vel_radial_sigma-pBins[i].vel_radial*pBins[i].vel_radial);
	    else
		pBins[i].vel_radial_sigma = 0.0;
	    for(j=0; j<3;j++) {
		pBins[i].L[j] /= pBins[i].dMassInBin;
		}
	    }
	}

    dsec = msrTime() - sec;
    msrprintf(msr,"Profiling complete, Wallclock: %f secs\n\n",dsec);

    if ( ppBins ) *ppBins = plcl->pkd->profileBins;
    if ( pnBins ) *pnBins = nBins+1;
    }

void msrInitGrid(MSR msr,int x,int y,int z) {
    struct inInitGrid in;
    in.n1 = x;
    in.n2 = y;
    in.n3 = z;
    in.a1 = x;
    in.s = 0;
    in.n = z;

    pstInitGrid(msr->pst, &in, sizeof(in), NULL, NULL);
    }

void msrGridProject(MSR msr,double x,double y,double z) {
    struct inGridProject in;
    in.r[0] = x;
    in.r[1] = y;
    in.r[2] = z;
    pstGridProject(msr->pst, &in, sizeof(in), NULL, NULL);
    }

#ifdef MDL_FFTW
void msrGridCreateFFT(MSR msr, int nGrid) {
    struct inGridCreateFFT in;
    in.nGrid = nGrid;
    pstGridCreateFFT(msr->pst, &in, sizeof(in), NULL, NULL);
    }

void msrGridDeleteFFT(MSR msr) {
    pstGridDeleteFFT(msr->pst, NULL, 0, NULL, NULL);
    }

/* Important: call msrGridCreateFFT() before, and msrGridDeleteFFT() after */
void msrMeasurePk(MSR msr,int iAssignment,int bInterlace,int nGrid,int nBins,uint64_t *nPk,float *fK,float *fPk) {
    struct inMeasurePk in;
    struct outMeasurePk *out;
    int nOut;
    int i;
    double fftNormalize;
    double sec,dsec;

    if (nGrid/2 < nBins) nBins = nGrid/2;
    assert(nBins <= PST_MAX_K_BINS);

    sec = msrTime();
    printf("Measuring P(k) with grid size %d (%d bins)...\n",nGrid,nBins);

    in.iAssignment = iAssignment;
    in.bInterlace = bInterlace;
    in.nGrid = nGrid;
    in.nBins = nBins;
    in.dTotalMass = msrTotalMass(msr);

    out = malloc(sizeof(struct outMeasurePk));
    assert(out != NULL);
    pstMeasurePk(msr->pst, &in, sizeof(in), out, &nOut);
    for( i=0; i<nBins; i++ ) {
	if ( out->nPower[i] == 0 ) fK[i] = fPk[i] = 0;
	else {
	    if (nPk) nPk[i] = out->nPower[i];
	    fK[i] = exp(out->fK[i]/out->nPower[i]);
	    fPk[i] = out->fPower[i]/out->nPower[i];
	    }
	}
    /* At this point, dPk[] needs to be corrected by the box size */
    free(out);

    dsec = msrTime() - sec;
    printf("P(k) Calculated, Wallclock: %f secs\n\n",dsec);
    }

void msrMeasureLinPk(MSR msr,int nGrid, double dA, double dBoxSize,
                    uint64_t *nPk,float *fK,float *fPk) {
    struct inMeasureLinPk in;
    struct outMeasureLinPk *out;
    int nOut;
    int i;
    double sec,dsec;

    sec = msrTime();

    in.nGrid = nGrid;
    in.nBins = msr->param.nBinsLinPk;
    in.dBoxSize = dBoxSize;
    in.dA = dA;
    in.iSeed = msr->param.iSeed;
    in.bFixed = msr->param.bFixedAmpIC;
    in.fPhase = msr->param.dFixedAmpPhasePI * M_PI;

    out = malloc(sizeof(struct outMeasureLinPk));
    assert(out != NULL);
    printf("Measuring P_lin(k) with grid size %d (%d bins)...\n",in.nGrid,in.nBins);
    pstMeasureLinPk(msr->pst, &in, sizeof(in), out, &nOut);
    for( i=0; i<in.nBins; i++ ) {
	if ( out->nPower[i] == 0 ) fK[i] = fPk[i] = 0;
	else {
	    if (nPk) nPk[i] = out->nPower[i];
	    fK[i] = out->fK[i]/out->nPower[i];
	    fPk[i] = out->fPower[i]/out->nPower[i];
	    }
	}
    /* At this point, dPk[] needs to be corrected by the box size */

    dsec = msrTime() - sec;
    printf("P_lin(k) Calculated, Wallclock: %f secs\n\n",dsec);
    }

void msrSetLinGrid(MSR msr,double dTime, int nGrid, int bKickClose, int bKickOpen){
    printf("Setting force grids of linear species with nGridLin = %d \n", nGrid);
    double sec, dsec;
    sec = msrTime();

    struct inSetLinGrid in;
    in.nGrid = nGrid;

    int do_DeltaRho_lin_avg = 1;
    in.a0 = in.a1 = in.a = csmTime2Exp(msr->param.csm, dTime);
    if (do_DeltaRho_lin_avg) {
	if (bKickClose) in.a0 = csmTime2Exp(msr->param.csm, dTime - 0.5*msrDelta(msr));
	if (bKickOpen)  in.a1 = csmTime2Exp(msr->param.csm, dTime + 0.5*msrDelta(msr));
	}

    in.dBSize = msr->param.dBoxSize;
    /* Parameters for the grid realization */
    in.iSeed = msr->param.iSeed;
    in.bFixed = msr->param.bFixedAmpIC;
    in.fPhase = msr->param.dFixedAmpPhasePI*M_PI;
    pstSetLinGrid(msr->pst, &in, sizeof(in), NULL, NULL);

    dsec = msrTime() - sec;
    printf("Force from linear species calculated, Wallclock: %f, secs\n\n", dsec);
    }

/* First call msrSetLinGrid() to setup the grid */
void msrLinearKick(MSR msr, double dTime, int bKickClose, int bKickOpen) {
    struct inLinearKick in;
    double dt = 0.5*msr->param.dDelta;
    double sec, dsec;

    printf("Applying Linear Kick...\n");
    sec = msrTime();
    in.dtOpen = in.dtClose = 0.0;
    if (msr->param.csm->val.bComove) {
	if (bKickClose) in.dtClose = csmComoveKickFac(msr->param.csm,dTime-dt,dt);
	if (bKickOpen) in.dtOpen = csmComoveKickFac(msr->param.csm,dTime,dt);
	}
    else {
	if (bKickClose) in.dtClose = dt;
	if (bKickOpen) in.dtOpen = dt;
	}
    pstLinearKick(msr->pst, &in, sizeof(in), NULL, NULL);
    dsec = msrTime() - sec;
    printf("Linear Kick Applied, Wallclock: %f secs\n\n",dsec);
    }
#endif

int msrGetParticles(MSR msr, int nIn, uint64_t *ID, struct outGetParticles *out) {
    int nOut;
    pstGetParticles(msr->pst, ID, sizeof(uint64_t)*nIn, out, &nOut);
    return nOut / sizeof(struct outGetParticles);
    }

void msrOutputOrbits(MSR msr,int iStep,double dTime) {
#ifdef _MSC_VER
    char achFile[MAX_PATH];
#else
    char achFile[PATH_MAX];
#endif
    FILE *fp;
    int i;

    if (msr->param.nOutputParticles) {
	struct outGetParticles particles[GET_PARTICLES_MAX];
	double dExp, dvFac;

	if (msr->param.csm->val.bComove) {
	    dExp = csmTime2Exp(msr->param.csm,dTime);
	    dvFac = 1.0/(dExp*dExp);
	    }
	else {
	    dExp = dTime;
	    dvFac = 1.0;
	    }

	msrGetParticles(msr,msr->param.nOutputParticles,msr->param.iOutputParticles,particles);
	msrBuildName(msr,achFile,iStep);
	strncat(achFile,".orb",256);

	fp = fopen(achFile,"w");
	if ( fp==NULL) {
	    printf("Could not create orbit File:%s\n",achFile);
	    _msrExit(msr,1);
	    }
	fprintf(fp,"%d %f\n",msr->param.nOutputParticles,dExp);
	for(i=0; i<msr->param.nOutputParticles; ++i) {
	    fprintf(fp,"%"PRIu64" %.8e %.16e %.16e %.16e %.8e %.8e %.8e %.8e\n",
		particles[i].id, particles[i].mass,
		particles[i].r[0], particles[i].r[1], particles[i].r[2],
		particles[i].v[0]*dvFac, particles[i].v[1]*dvFac, particles[i].v[2]*dvFac,
		particles[i].phi);
	    }
	fclose(fp);
	}
    }
