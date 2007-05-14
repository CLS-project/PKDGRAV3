#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#define _LARGEFILE_SOURCE 
#define _FILE_OFFSET_BITS 64 
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include "mdl.h"
#ifdef USE_BT
#include "bt.h"
#endif
#include "master.h"
#include "outtype.h"
#include "smoothfcn.h"
#ifdef USE_MDL_IO
#include "io.h"
#endif

#ifdef USE_MDL_IO
static void main_io(MDL mdl) {
    IO io;

    ioInitialize(&io,mdl);
    ioAddServices(io,mdl);

    if ( mdlSelf(mdl) == 0 ) {
	mdlSetComm(mdl,1);
    }

    mdlprintf( mdl, "I/O thread %d started\n", mdlSelf(mdl) );
    mdlHandler(mdl);
    mdlprintf( mdl, "I/O thread %d terminated\n", mdlSelf(mdl) );

    if ( mdlSelf(mdl) == 0 ) {
	int id;
	mdlSetComm(mdl,0);
	for( id=1; id<mdlIO(mdl); ++id ) {
	    mdlReqService(mdl,id,SRV_STOP,NULL,0);
	    mdlGetReply(mdl,id,NULL,NULL);
	}
    }
}
#endif

void main_ch(MDL mdl) {
    PST pst;
    LCL lcl;

    lcl.pszDataPath = (char *)getenv("PTOOLS_DATA_PATH");
    lcl.pkd = NULL;
    pstInitialize(&pst,mdl,&lcl);

    pstAddServices(pst,mdl);

    mdlHandler(mdl);

    pstFinish(pst);
    }

/*DEBUG Should opaque nature of "msr" be enforced in main()? -- DCR*/ 

int main(int argc,char **argv) {
    MDL mdl;
    MSR msr;
    FILE *fpLog = NULL;
    char achFile[256];			/*DEBUG use MAXPATHLEN here (& elsewhere)? -- DCR*/
    double dTime;
    double E=0,T=0,U=0,Eth=0,L[3]={0,0,0};
    double dWMax=0,dIMax=0,dEMax=0,dMass=0,dMultiEff=0;
    long lSec=0,lStart;
    int i,iStep,iSec=0,iStop=0,nActive;
    int bGasOnly,bSymmetric;
    char achBaseMask[256];
    int nFOFsDone;
#ifdef TINY_PTHREAD_STACK
    static int first = 1;
    static char **save_argv;

    /*
     * Hackery to get around SGI's tiny pthread stack.
     * Main will be called twice.  The second time, argc and argv
     * will be garbage, so we have to save them from the first.
     * Another way to do this would involve changing the interface
     * to mdlInitialize(), so that this hackery could be hidden
     * down there.
     */
    if(first) {
	save_argv = malloc((argc+1)*sizeof(*save_argv));
	for(i = 0; i < argc; i++)
	    save_argv[i] = strdup(argv[i]);
	save_argv[argc] = NULL;
	}
    else {
	argv = save_argv;
	}
    first = 0;
#endif /* TINY_PTHREAD_STACK */

    printf("%s\n", PACKAGE_STRING );

#ifdef USE_BT
    bt_initialize();
#endif
#ifndef CCC
    /* no stdout buffering */
    setbuf(stdout,(char *) NULL);
#endif

    lStart=time(0);
#ifdef USE_MDL_IO
    mdlInitialize(&mdl,argv,main_ch,main_io);
#else
    mdlInitialize(&mdl,argv,main_ch);
#endif
    for(argc = 0; argv[argc]; argc++); /* some MDLs can trash argv */
    msrInitialize(&msr,mdl,argc,argv);

    (void) strncpy(achBaseMask,msr->param.achDigitMask,256);

    /*
    ** Establish safety lock.
    */
    if (!msrGetLock(msr)) {
	msrFinish(msr);
	mdlFinish(mdl);
	return 1;
	}
    /*
    ** Read in the binary file, this may set the number of timesteps or
    ** the size of the timestep when the zto parameter is used.
    */
      if(msr->param.bHeliocentric){
	dTime = msrReadSS(msr); /* must use "Solar System" (SS) I/O format... */
    }else{
      dTime = msrReadTipsy(msr);  /*read initial conditions */
    } 
    msrInitStep(msr);
    if (prmSpecified(msr->prm,"dSoft")) msrSetSoft(msr,msrSoft(msr));
    /*
    ** If the simulation is periodic make sure to wrap all particles into
    ** the "unit" cell. Doing a drift of 0.0 will always take care of this.
    */
    msrDrift(msr,dTime,0.0);

    if (msrSteps(msr) > 0) {
	if (msrComove(msr)) {
	    msrSwitchTheta(msr,dTime);
	    }
	/*
	** Now we have all the parameters for the simulation we can make a 
	** log file entry.
	*/
	if (msrLogInterval(msr)) {
	    sprintf(achFile,"%s.log",msrOutName(msr));
	    fpLog = fopen(achFile,"w");
	    assert(fpLog != NULL);
	    setbuf(fpLog,(char *) NULL); /* no buffering */
	    /*
	    ** Include a comment at the start of the log file showing the
	    ** command line options.
	    */
	    fprintf(fpLog,"# ");
	    for (i=0;i<argc;++i) fprintf(fpLog,"%s ",argv[i]);
	    fprintf(fpLog,"\n");
	    msrLogParams(msr,fpLog);
	    }
	/*
	** Build tree, activating all particles first (just in case).
	*/
	msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE);
	msrDomainDecomp(msr,0,1,0);
	msrUpdateSoft(msr,dTime);
	msrBuildTree(msr,dMass,dTime);
	if (msrDoGravity(msr)) {
	  msrGravity(msr,dTime,msr->param.iStartStep,&iSec,&dWMax,&dIMax,&dEMax,&nActive);
	  if (msr->param.bGravStep && msr->param.iTimeStepCrit == -1) {
	    msrGravity(msr,dTime,msr->param.iStartStep,&iSec,&dWMax,&dIMax,&dEMax,&nActive);
	  }
	    msrCalcEandL(msr,MSR_INIT_E,dTime,&E,&T,&U,&Eth,L);
	    dMultiEff = 1.0;
	    if (msrLogInterval(msr)) {
		(void) fprintf(fpLog,"%e %e %.16e %e %e %e %.16e %.16e %.16e "
			       "%i %e %e %e %e\n",dTime,
			       1.0/csmTime2Exp(msr->param.csm,dTime)-1.0,
			       E,T,U,Eth,L[0],L[1],L[2],iSec,dWMax,dIMax,dEMax,
			       dMultiEff);
		}
	    }
#ifdef RELAXATION
	msrInitRelaxation(msr);
#endif
	for (iStep=msr->param.iStartStep+1;iStep<=msrSteps(msr);++iStep) {
	    if (msrComove(msr)) {
		msrSwitchTheta(msr,dTime);
		}
	    dMultiEff = 0.0;
	    lSec = time(0);

	    msrTopStepKDK(msr,iStep-1,dTime,
			  msrDelta(msr),0,0,msrMaxRung(msr),1,
			  &dMultiEff,&dWMax,&dIMax,
			  &dEMax,&iSec);

	    dTime += msrDelta(msr);
	    /*
	    ** Output a log file line if requested.
	    ** Note: no extra gravity calculation required.
	    */
	    if (msrLogInterval(msr) && iStep%msrLogInterval(msr) == 0) {
		msrCalcEandL(msr,MSR_STEP_E,dTime,&E,&T,&U,&Eth,L);
		lSec = time(0) - lSec;
		(void) fprintf(fpLog,"%e %e %.16e %e %e %e %.16e %.16e "
			       "%.16e %li %e %e %e %e\n",dTime,
			       1.0/csmTime2Exp(msr->param.csm,dTime)-1.0,
			       E,T,U,Eth,L[0],L[1],L[2],lSec,dWMax,dIMax,
			       dEMax,dMultiEff);
		}
#ifdef RELAXATION 
	    if ( msr->param.bTraceRelaxation) {
		msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
		msrDomainDecomp(msr,0,1,0);
		msrBuildTree(msr,dMass,dTime);
		msrRelaxation(msr,dTime,msrDelta(msr),SMX_RELAXATION,0);
		}
#endif	
	    /*
	    ** Check for user interrupt.
	    */
	    iStop = msrCheckForStop(msr);
	    /*
	    ** Output if 1) we've hit an output time
	    **           2) We are stopping
	    **           3) we're at an output interval
	    */
	    if (msrOutTime(msr,dTime) || iStep == msrSteps(msr) || iStop ||
		(msrOutInterval(msr) > 0 &&	iStep%msrOutInterval(msr) == 0)) {

		msrReorder(msr);
		sprintf(achFile,msr->param.achDigitMask,msrOutName(msr),iStep);
		msrWriteTipsy(msr,achFile,dTime);
		if (msrDoDensity(msr)) {
		    msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
		    msrDomainDecomp(msr,0,1,0);
		    msrBuildTree(msr,dMass,dTime);
		    bGasOnly = 0;
		    bSymmetric = 0; /* FOR TESTING!!*/
		    msrSmooth(msr,dTime,SMX_DENSITY,bGasOnly,bSymmetric);
		    }
		nFOFsDone = 0;
		while( msr->param.nFindGroups > nFOFsDone) {
		    /*
		    ** Build tree, activating all particles first (just in case).
		    */	
		    msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE);
		    msrDomainDecomp(msr,0,1,0);
		    msrBuildTree(msr,dMass,dTime);
		    msrFof(msr,nFOFsDone,SMX_FOF,0);
		    msrGroupMerge(msr);
		    if(msr->param.nBins > 0) msrGroupProfiles(msr,nFOFsDone,SMX_FOF,0);
		    msrReorder(msr);
		    sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
		    for(i=0;i<=nFOFsDone;++i)strncat(achFile,".fof",256);
		    msrOutArray(msr,achFile,OUT_GROUP_ARRAY);
		    sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
		    for(i=0;i<=nFOFsDone;++i)strncat(achFile,".stats",256);
		    msrOutGroups(msr,achFile,OUT_GROUP_STATS,dTime);
		    if( msr->nBins > 0){
			sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
			for(i=0;i<=nFOFsDone;++i)strncat(achFile,".pros",256);
			msrOutGroups(msr,achFile,OUT_GROUP_PROFILES,dTime);
			}
		    sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
		    for(i=0;i<=nFOFsDone;++i)strncat(achFile,".grps",256);
		    if(	msr->param.bStandard) msrOutGroups(msr,achFile,OUT_GROUP_TIPSY_STD,dTime);
		    else msrOutGroups(msr,achFile,OUT_GROUP_TIPSY_NAT,dTime);			  
		    nFOFsDone++;
		    }
		if (msr->param.nFindGroups > 0) msrDeleteGroups(msr);
#ifdef RELAXATION	
		if ( msr->param.bTraceRelaxation) {
		    msrReorder(msr);
		    sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
		    strncat(achFile,".relax",256);
		    msrOutArray(msr,achFile,OUT_RELAX_ARRAY);
		    }	
#endif 
		if (msrDoDensity(msr)) {
		    msrReorder(msr);
		    sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
		    strncat(achFile,".den",256);
		    msrOutArray(msr,achFile,OUT_DENSITY_ARRAY);
		    }
		if (msr->param.bDoRungOutput) {
		    msrReorder(msr);
		    sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
		    strncat(achFile,".rung",256);
		    msrOutArray(msr,achFile,OUT_RUNG_ARRAY);
		    }
		if (msr->param.bDoSoftOutput) {
		    msrReorder(msr);
		    sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
		    strncat(achFile,".soft",256);
		    msrOutArray(msr,achFile,OUT_SOFT_ARRAY);
		    }
		if (msr->param.bDodtOutput) {
		    msrReorder(msr);
		    sprintf(achFile,achBaseMask,msrOutName(msr),iStep);
		    strncat(achFile,".dt",256);
		    msrOutArray(msr,achFile,OUT_DT_ARRAY);
		    }
		/*
		** Don't allow duplicate outputs.
		*/
		while (msrOutTime(msr,dTime));
		}
	    if (!iStop && msr->param.iWallRunTime > 0) {
		if (msr->param.iWallRunTime*60 - (time(0)-lStart) < ((int) (lSec*1.5)) ) {
		    printf("RunTime limit exceeded.  Writing checkpoint and exiting.\n");
		    printf("    iWallRunTime(sec): %d   Time running: %ld   Last step: %ld\n",
			   msr->param.iWallRunTime*60,time(0)-lStart,lSec);
		    iStop = 1;
		    }
		}
	    if (iStop || iStep == msrSteps(msr)) {
		/*
		** Want to write an output file here really...
		*/
		}
	    if (iStop) break;
	    }
	if (msrLogInterval(msr)) (void) fclose(fpLog);
	}
	
    else {
	struct inInitDt in;
	msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);

	in.dDelta = 1e37;		/* large number */
	pstInitDt(msr->pst,&in,sizeof(in),NULL,NULL);
    
	fprintf(stderr,"Initialized Accel and dt\n");

	msrDomainDecomp(msr,0,1,0);
	msrUpdateSoft(msr,dTime);
	msrBuildTree(msr,dMass,dTime);

	if (msrDoGravity(msr)) {
	  msrGravity(msr,dTime,msr->param.iStartStep,&iSec,&dWMax,&dIMax,&dEMax,&nActive);
	  if (msr->param.bGravStep && msr->param.iTimeStepCrit == -1) {
	    msrGravity(msr,dTime,msr->param.iStartStep,&iSec,&dWMax,&dIMax,&dEMax,&nActive);
	  }
	    msrReorder(msr);
	    sprintf(achFile,"%s.accg",msrOutName(msr));
	    msrOutVector(msr,achFile,OUT_ACCEL_VECTOR);
	    sprintf(achFile,"%s.pot",msrOutName(msr));
	    msrReorder(msr);
	    msrOutArray(msr,achFile,OUT_POT_ARRAY);
	    }
	if (msrDoDensity(msr) || msr->param.bDensityStep) {
	    bGasOnly = 0;
	    bSymmetric = 1;
	    msrSmooth(msr,dTime,SMX_DENSITY,bGasOnly,bSymmetric);
	    msrReorder(msr);
	    sprintf(achFile,"%s.den",msrOutName(msr));
	    msrOutArray(msr,achFile,OUT_DENSITY_ARRAY);
	    } 
	nFOFsDone = 0;
	while ( msr ->param.nFindGroups > nFOFsDone) {
	    /*
	    ** Build tree, activating all particles first (just in case).
	    */
	    msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE);
	    msrDomainDecomp(msr,0,1,0);
	    msrBuildTree(msr,dMass,dTime);
	    msrFof(msr,nFOFsDone,SMX_FOF,0);
	    msrGroupMerge(msr);
	    if(msr->param.nBins > 0) msrGroupProfiles(msr,nFOFsDone,SMX_FOF,0);
	    msrReorder(msr);
	    sprintf(achFile,"%s.%i.fof",msrOutName(msr),nFOFsDone);
	    msrOutArray(msr,achFile,OUT_GROUP_ARRAY);
	    sprintf(achFile,"%s.%i.stats",msrOutName(msr),nFOFsDone);
	    msrOutGroups(msr,achFile,OUT_GROUP_STATS,dTime);			
	    sprintf(achFile,"%s.%i.grps",msrOutName(msr),nFOFsDone);
	    if(	msr->param.bStandard) msrOutGroups(msr,achFile,OUT_GROUP_TIPSY_STD,dTime);
	    else msrOutGroups(msr,achFile,OUT_GROUP_TIPSY_NAT,dTime);
	    if( msr->nBins > 0){
		sprintf(achFile,"%s.%i.pros",msrOutName(msr),nFOFsDone);
		msrOutGroups(msr,achFile,OUT_GROUP_PROFILES,dTime);
		}			
	    nFOFsDone++;	
	    }	
	if (msr->param.nFindGroups > 0) msrDeleteGroups(msr);
	if (msrDoDensity(msr)) {
	    msrReorder(msr);
	    sprintf(achFile,"%s.den",msrOutName(msr));
	    msrOutArray(msr,achFile,OUT_DENSITY_ARRAY);
	    }
	if (msrDoGravity(msr)) {
	    if (msr->param.bGravStep) {
		fprintf(stderr,"Adding GravStep dt\n");
		msrGravStep(msr,dTime);
		}
	    if (msr->param.bAccelStep) {
		fprintf(stderr,"Adding AccelStep dt\n");
		msrAccelStep(msr,dTime);
		}
	    }
	if (msr->param.bDensityStep) {
	    fprintf(stderr,"Adding DensStep dt\n");
	    msrDensityStep(msr,dTime);
	    }
	msrReorder(msr);
	sprintf(achFile,"%s.dt",msrOutName(msr));
	msrOutArray(msr,achFile,OUT_DT_ARRAY);
	if(msr->param.bDensityStep || msrDoGravity(msr)) {
	    msrDtToRung(msr,0,msrDelta(msr),1);
	    }
	}
	
    msrFinish(msr);
    mdlFinish(mdl);
    return 0;
    }
