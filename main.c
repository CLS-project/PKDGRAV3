#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#define _LARGEFILE_SOURCE
#define _FILE_OFFSET_BITS 64
#ifdef ENABLE_FE
#include <fenv.h>
#endif
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
#ifdef USE_PYTHON
#include "pkdpython.h"
#endif
#include <signal.h>

static time_t timeGlobalSignalTime = 0;
static void USR1_handler(int signo) {
    signal(SIGUSR1,USR1_handler);
    timeGlobalSignalTime = time(0);
    }

static int bGlobalOutput = 0;
static void USR2_handler(int signo) {
    signal(SIGUSR2,USR2_handler);
    bGlobalOutput = 1;
    }

void * main_ch(MDL mdl) {
    PST pst;
    LCL lcl;

    /* a USR1 signal indicates that the queue wants us to exit */
    timeGlobalSignalTime = 0;
    signal(SIGUSR1,NULL);
    signal(SIGUSR1,USR1_handler);

    /* a USR2 signal indicates that we should write an output when convenient */
    bGlobalOutput = 0;
    signal(SIGUSR2,NULL);
    signal(SIGUSR2,USR2_handler);

    lcl.pkd = NULL;
    pstInitialize(&pst,mdl,&lcl);

    pstAddServices(pst,mdl);

    mdlHandler(mdl);

    pstFinish(pst);
    return NULL;
    }

/*
** This is invoked instead of main_ch for the "master" process.
*/
void * master_ch(MDL mdl) {
    int argc = mdl->base.argc;
    char **argv = mdl->base.argv;
    MSR msr;
    FILE *fpLog = NULL;
    char achFile[256];			/*DEBUG use MAXPATHLEN here (& elsewhere)? -- DCR*/
    double dTime;
    double E=0,T=0,U=0,Eth=0,L[3]={0,0,0},F[3]={0,0,0},W=0;
    double dMultiEff=0;
    long lSec=0,lStart;
    int i,iStep,iStartStep,iSec=0,iStop=0;
    uint64_t nActive;
    int bKickOpen = 0;
    int64_t nRungs[MAX_RUNG+1];
    uint8_t uRungMax;
    double diStep;
    double ddTime;
    int bRestore;

    lStart=time(0);

    printf("%s\n", PACKAGE_STRING );

    bRestore = msrInitialize(&msr,mdl,argc,argv);

    /*
    ** Establish safety lock.
    */
    if (!msrGetLock(msr)) {
	msrFinish(msr);
	return NULL;
	}

    /* a USR1 signal indicates that the queue wants us to exit */
    timeGlobalSignalTime = 0;
    signal(SIGUSR1,NULL);
    signal(SIGUSR1,USR1_handler);

    /* a USR2 signal indicates that we should write an output when convenient */
    bGlobalOutput = 0;
    signal(SIGUSR2,NULL);
    signal(SIGUSR2,USR2_handler);

    /*
    ** Output the host names to make troubleshooting easier
    */
    msrHostname(msr);

    /* Restore from a previous checkpoint */
    if (bRestore) {
	dTime = msrRestore(msr);
	iStartStep = msr->iCheckpointStep;
	msrInitStep(msr);
	if (prmSpecified(msr->prm,"dSoft")) msrSetSoft(msr,msrSoft(msr));
	iStep = msrSteps(msr); /* 0=analysis, >1=simulate, <0=python */
	uRungMax = msr->iCurrMaxRung;
	}

    /* Generate the initial particle distribution */
    else if (prmSpecified(msr->prm,"nGrid")) {
	iStartStep = msr->param.iStartStep; /* Should be zero */
#ifdef MDL_FFTW
	dTime = msrGenerateIC(msr); /* May change nSteps/dDelta */
	if ( msr->param.bWriteIC ) {
	    msrBuildIoName(msr,achFile,0);
	    msrWrite( msr,achFile,dTime,msr->param.bWriteIC-1);
	    }
#else
	printf("To generate initial conditions, compile with FFTW\n");
	msrFinish(msr);
	return NULL;
#endif
	msrInitStep(msr);
	if (prmSpecified(msr->prm,"dSoft")) msrSetSoft(msr,msrSoft(msr));
	iStep = msrSteps(msr); /* 0=analysis, >1=simulate, <0=python */
	}

    /* Read in a binary file */
    else if ( msr->param.achInFile[0] ) {
	iStartStep = msr->param.iStartStep; /* Should be zero */
	dTime = msrRead(msr,msr->param.achInFile); /* May change nSteps/dDelta */
	msrInitStep(msr);
	if (msr->param.bAddDelete) msrGetNParts(msr);
	if (prmSpecified(msr->prm,"dRedFrom")) {
	    double aOld, aNew;
	    aOld = csmTime2Exp(msr->param.csm,dTime);
	    aNew = 1.0 / (1.0 + msr->param.dRedFrom);
	    dTime = msrAdjustTime(msr,aOld,aNew);
	    /* Seriously, we shouldn't need to send parameters *again*.
	       When we remove sending parameters, we should remove this. */
	    msrInitStep(msr);
	    }
	if (prmSpecified(msr->prm,"dSoft")) msrSetSoft(msr,msrSoft(msr));
	iStep = msrSteps(msr); /* 0=analysis, >1=simulate, <0=python */
	}
#ifdef USE_PYTHON
    else if ( msr->param.achScriptFile[0] ) {
	iStep = -1;
	}
#endif
    else {
	printf("No input file specified\n");
	msrFinish(msr);
	return NULL;
	}

    /* Adjust theta for gravity calculations. */
    if (msrComove(msr)) msrSwitchTheta(msr,dTime);

    /* Analysis mode */
    if (iStep==0) {
	if (msrDoGravity(msr) ||msrDoGas(msr)) {
	    msrActiveRung(msr,0,1); /* Activate all particles */
	    msrDomainDecomp(msr,0,0,0);
	    msrUpdateSoft(msr,dTime);
	    msrBuildTree(msr,dTime,msr->param.bEwald);
#ifdef MDL_FFTW
	    if (msr->param.nGridPk>0) {
		msrOutputPk(msr,iStartStep,dTime);
		}
#endif
	    if (msrDoGravity(msr)) {
		msrGravity(msr,0,MAX_RUNG,ROOT,0,dTime,iStartStep,0,0,
		    msr->param.bEwald,msr->param.nGroup,&iSec,&nActive);
		msrMemStatus(msr);
		if (msr->param.bGravStep) {
		    msrBuildTree(msr,dTime,msr->param.bEwald);
		    msrGravity(msr,0,MAX_RUNG,ROOT,0,dTime,iStartStep,0,0,
			msr->param.bEwald,msr->param.nGroup,&iSec,&nActive);
		    msrMemStatus(msr);
		    }
		}
		
	    if (msrDoGas(msr)) {
		/* Initialize SPH, Cooling and SF/FB and gas time step */
		msrCoolSetup(msr,dTime);
		/* Fix dTuFac conversion of T in InitSPH */
		msrInitSph(msr,dTime);
		}
		   
	    msrUpdateRung(msr,0); /* set rungs for output */
	    }
	msrOutput(msr,0,dTime,0);  /* JW: Will trash gas density */
	}

#ifdef USE_PYTHON
    else if ( iStep < 0 ) {
	PPY ppy;
	ppyInitialize(&ppy,msr,dTime);
	msr->prm->script_argv[0] = msr->param.achScriptFile;
	ppyRunScript(ppy,msr->prm->script_argc,msr->prm->script_argv);
	ppyFinish(ppy);
	}
#endif

    /* Simulation mode */
    else {
	/*
	** Now we have all the parameters for the simulation we can make a
	** log file entry.
	*/
	if (msrLogInterval(msr)) {
	    sprintf(achFile,"%s.log",msrOutName(msr));
	    fpLog = fopen(achFile,"a");
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

	if (bRestore) goto test_continue;

	if (msr->param.bLightCone && msrComove(msr)) {
	    printf("LightCone output will begin at z=%.10g\n",
		1.0/csmComoveLookbackTime2Exp(msr->param.csm,1.0 / dLightSpeedSim(msr->param.dBoxSize)) - 1.0 );;
	    }

	/*
	** Build tree, activating all particles first (just in case).
	*/
	msrActiveRung(msr,0,1); /* Activate all particles */
	msrDomainDecomp(msr,0,0,0);
	msrUpdateSoft(msr,dTime);
	msrBuildTree(msr,dTime,msr->param.bEwald);
#ifdef MDL_FFTW
	if (msr->param.nGridPk>0) {
	    msrOutputPk(msr,iStartStep,dTime);
	    }
#endif
	if (msrDoGravity(msr) && !msr->param.bHSDKD) {
	    if (msr->param.bNewKDK) bKickOpen = 1;
	    else bKickOpen = 0;
	    uRungMax = msrGravity(msr,0,MAX_RUNG,ROOT,0,dTime,iStartStep,0,bKickOpen,msr->param.bEwald,msr->param.nGroup,&iSec,&nActive);
	    msrMemStatus(msr);
	    if (msr->param.bGravStep) {
		assert(msr->param.bNewKDK == 0);    /* for now! */
		msrBuildTree(msr,dTime,msr->param.bEwald);
		msrGravity(msr,0,MAX_RUNG,ROOT,0,dTime,iStartStep,0,0,msr->param.bEwald,msr->param.nGroup,&iSec,&nActive);
		msrMemStatus(msr);
		if (msr->param.bHSDKD) {
		    msrAccelStep(msr,0,MAX_RUNG,dTime);
		    msrUpdateRung(msr,0);
		    }
		}
	    }
	if (msrDoGas(msr)) {
	    /* Initialize SPH, Cooling and SF/FB and gas time step */
	    msrCoolSetup(msr,dTime);
	    /* Fix dTuFac conversion of T in InitSPH */
	    msrInitSph(msr,dTime);
	    }

	msrCalcEandL(msr,MSR_INIT_E,dTime,&E,&T,&U,&Eth,L,F,&W);
	dMultiEff = 1.0;
	if (msrLogInterval(msr)) {
		(void) fprintf(fpLog,"%e %e %.16e %e %e %e %.16e %.16e %.16e "
			       "%.16e %.16e %.16e %.16e %i %e\n",dTime,
			       1.0/csmTime2Exp(msr->param.csm,dTime)-1.0,
			       E,T,U,Eth,L[0],L[1],L[2],F[0],F[1],F[2],W,iSec,dMultiEff);
	    }
	if ( msr->param.bTraceRelaxation) {
	    msrInitRelaxation(msr);
	    }
    test_continue:

	for (iStep=iStartStep+1;iStep<=msrSteps(msr)&&!iStop;++iStep) {
	    if (msrComove(msr)) msrSwitchTheta(msr,dTime);
	    dMultiEff = 0.0;
	    lSec = time(0);
	    msrLightConeOpen(msr,iStep);
	    if (msr->param.bHSDKD) {
		/* Perform select */
		msrActiveRung(msr,0,1); /* Activate all particles */
		msrBuildTree(msr,dTime,msr->param.bEwald);
		msrGravity(msr,0,MAX_RUNG,ROOT,0,dTime,iStep-1,0,0,msr->param.bEwald,msr->param.nGroup,&iSec,&nActive);
		msrMemStatus(msr);
		msrAccelStep(msr,0,MAX_RUNG,dTime);
		msrUpdateRung(msr,0);

		msrTopStepHSDKD(msr,iStep-1,dTime,
		    msrDelta(msr),0,0,msrMaxRung(msr),1,
		    &dMultiEff,&iSec);
		}
	    else if (msr->param.bNewKDK) {
		diStep = (double)(iStep-1);
		ddTime = dTime;
		msrNewTopStepKDK(msr,0,0,&diStep,&ddTime,&uRungMax,&iSec);
		}
	    else {
		msrTopStepKDK(msr,iStep-1,dTime,
		    msrDelta(msr),0,0,msrMaxRung(msr),1,
		    &dMultiEff,&iSec);
		}
	    dTime += msrDelta(msr);
	    lSec = time(0) - lSec;
	    msrMemStatus(msr);

	    msrLightConeClose(msr,iStep);

#ifdef MDL_FFTW
	    if (msr->param.iPkInterval && iStep%msr->param.iPkInterval == 0) {
		msrOutputPk(msr,iStep,dTime);
		}
#endif
	    if (msr->param.iFofInterval<0 || msr->param.iFofInterval>0 && iStep%msr->param.iFofInterval==0) {
		msrNewFof(msr,dTime);
		if (msr->param.iFofInterval<0) iStep = iStep<<(-msr->param.iFofInterval);
		msrBuildName(msr,achFile,iStep);
		strncat(achFile,".fofstats",256);
		msrHopWrite(msr,achFile);
		}
	    /*
	    ** Output a log file line if requested.
	    ** Note: no extra gravity calculation required.
	    */
	    if (msrLogInterval(msr) && iStep%msrLogInterval(msr) == 0) {
		msrCalcEandL(msr,MSR_STEP_E,dTime,&E,&T,&U,&Eth,L,F,&W);
		(void) fprintf(fpLog,"%e %e %.16e %e %e %e %.16e %.16e "
			       "%.16e %.16e %.16e %.16e %.16e %li %e\n",dTime,
			       1.0/csmTime2Exp(msr->param.csm,dTime)-1.0,
			       E,T,U,Eth,L[0],L[1],L[2],F[0],F[1],F[2],W,lSec,dMultiEff);
		}
	    if ( msr->param.bTraceRelaxation) {
		msrActiveRung(msr,0,1); /* Activate all particles */
		msrDomainDecomp(msr,0,0,0);
		msrBuildTree(msr,dTime,0);
		msrRelaxation(msr,dTime,msrDelta(msr),SMX_RELAXATION,0);
		}
	    /*
	    ** Check for user interrupt.
	    */
	    iStop = msrCheckForStop(msr);

	    /*
	    ** Check to see if the runtime has been exceeded.
	    */
	    if (!iStop && msr->param.iWallRunTime > 0) {
		if (msr->param.iWallRunTime*60 - (time(0)-lStart) < ((int) (lSec*1.5)) ) {
		    printf("RunTime limit exceeded.  Writing checkpoint and exiting.\n");
		    printf("    iWallRunTime(sec): %d   Time running: %ld   Last step: %ld\n",
			   msr->param.iWallRunTime*60,time(0)-lStart,lSec);
		    iStop = 1;
		    }
		}

	    /* Check to see if there should be an output */
	    if (!iStop && timeGlobalSignalTime>0) { /* USR1 received */
		if ( (time(0)+(lSec*1.5)) > timeGlobalSignalTime+msr->param.iSignalSeconds) {
		    printf("RunTime limit exceeded.  Writing checkpoint and exiting.\n");
		    printf("    iSignalSeconds: %d   Time running: %ld   Last step: %ld\n",
			msr->param.iSignalSeconds,time(0)-lStart,lSec);
		    iStop = 1;
		    }
		}

	    /*
	    ** Output if 1) we've hit an output time
	    **           2) We are stopping
	    **           3) we're at an output interval
	    */
	    if (msrCheckInterval(msr)>0 &&
		(iStop || bGlobalOutput
		    || (iStep%msrCheckInterval(msr) == 0) ) ) {
		bGlobalOutput = 0;
		msrCheckpoint(msr,iStep,dTime);
		}

	    if (bGlobalOutput || msrOutTime(msr,dTime) || iStep == msrSteps(msr) || iStop ||
		    (msrOutInterval(msr) > 0 && iStep%msrOutInterval(msr) == 0)) {
		bGlobalOutput = 0;
		msrOutput(msr,iStep,dTime, 0);
		}
	    }
	if (msrLogInterval(msr)) (void) fclose(fpLog);
	}

    msrFinish(msr);
    return NULL;
    }

#ifdef FC_DUMMY_MAIN
int FC_DUMMY_MAIN() { return 1; }
#endif

int FC_MAIN(int argc,char **argv) {
#ifdef USE_BT
    bt_initialize();
#endif
#ifdef ENABLE_FE
    feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
#endif
#ifndef CCC
    /* no stdout buffering */
    setbuf(stdout,(char *) NULL);
#endif

    mdlLaunch(argc,argv,master_ch,main_ch);

    return 0;
    }
