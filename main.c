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

void * main_ch(MDL mdl) {
    PST pst;
    LCL lcl;

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
    int i,iStep,iSec=0,iStop=0;
    uint64_t nActive;

    lStart=time(0);

    printf("%s\n", PACKAGE_STRING );

    msrInitialize(&msr,mdl,argc,argv);

    /*
    ** Establish safety lock.
    */
    if (!msrGetLock(msr)) {
	msrFinish(msr);
	return NULL;
	}

    /*
    ** Output the host names to make troubleshooting easier
    */
    msrHostname(msr);
    /*
    ** Read in the binary file, this may set the number of timesteps or
    ** the size of the timestep when the zto parameter is used.
    */
#ifdef USE_GRAFIC
    if (prmSpecified(msr->prm,"nGrid")) {
	dTime = msrGenerateIC(msr);
	msrInitStep(msr);
	if (prmSpecified(msr->prm,"dSoft")) msrSetSoft(msr,msrSoft(msr));
	}
    else {
#endif
	if ( msr->param.achInFile[0] ) {
	    dTime = msrRead(msr,msr->param.achInFile);
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
	    }
	else {
#ifdef USE_PYTHON
	    if ( !msr->param.achScriptFile[0] ) {
#endif
		printf("No input file specified\n");
		return NULL;
		}
#ifdef USE_PYTHON
	    }
#endif
#ifdef USE_GRAFIC
	}
#endif

    if ( msr->param.bWriteIC ) {
	msrBuildIoName(msr,achFile,0);
	msrWrite( msr,achFile,dTime,msr->param.bWriteIC-1);
	}

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

#ifdef USE_PYTHON
    /* If a script file was specified, enter analysis mode */
    if ( msr->param.achScriptFile[0] ) iStep = 0;
    else
#endif
	iStep = msrSteps(msr);
    if (iStep > 0) {
	if (msrComove(msr)) {
	    msrSwitchTheta(msr,dTime);
	    }
	/*
	** Build tree, activating all particles first (just in case).
	*/
	msrActiveRung(msr,0,1); /* Activate all particles */
	msrDomainDecomp(msr,0,0,0);
	msrUpdateSoft(msr,dTime);
	msrBuildTree(msr,dTime,msr->param.bEwald);
	if (msrDoGravity(msr) && !msr->param.bHSDKD) {
	    msrGravity(msr,0,MAX_RUNG,ROOT,0,dTime,msr->param.iStartStep,msr->param.bEwald,msr->param.nGroup,&iSec,&nActive);
	    msrMemStatus(msr);
	    if (msr->param.bGravStep) {
		msrBuildTree(msr,dTime,msr->param.bEwald);
		msrGravity(msr,0,MAX_RUNG,ROOT,0,dTime,msr->param.iStartStep,msr->param.bEwald,msr->param.nGroup,&iSec,&nActive);
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
	for (iStep=msr->param.iStartStep+1;iStep<=msrSteps(msr)&&!iStop;++iStep) {
	    if (msrComove(msr)) {
		msrSwitchTheta(msr,dTime);
		}
	    dMultiEff = 0.0;
	    lSec = time(0);
	    if (msr->param.bHSDKD) {

		/* Perform select */
		msrActiveRung(msr,0,1); /* Activate all particles */
		msrBuildTree(msr,dTime,msr->param.bEwald);
		msrGravity(msr,0,MAX_RUNG,ROOT,0,dTime,iStep-1,msr->param.bEwald,msr->param.nGroup,&iSec,&nActive);
		msrMemStatus(msr);
		msrAccelStep(msr,0,MAX_RUNG,dTime);
		msrUpdateRung(msr,0);

		msrTopStepHSDKD(msr,iStep-1,dTime,
		    msrDelta(msr),0,0,msrMaxRung(msr),1,
		    &dMultiEff,&iSec);
		}
	    else {
		msrTopStepKDK(msr,iStep-1,dTime,
		    msrDelta(msr),0,0,msrMaxRung(msr),1,
		    &dMultiEff,&iSec);
		}

	    dTime += msrDelta(msr);
	    lSec = time(0) - lSec;

	    msrMemStatus(msr);

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

	    /*
	    ** Output if 1) we've hit an output time
	    **           2) We are stopping
	    **           3) we're at an output interval
	    */
	    if (msrOutTime(msr,dTime) || iStep == msrSteps(msr) || iStop ||
		    (msrOutInterval(msr) > 0 && iStep%msrOutInterval(msr) == 0) ||
		    (msrCheckInterval(msr) > 0 && iStep%msrCheckInterval(msr) == 0)) {
		msrOutput(msr,iStep,dTime, msrCheckInterval(msr)>0
			  && (iStep%msrCheckInterval(msr) == 0
			      || iStep == msrSteps(msr) || iStop));
		}
	    }
	}

    /* No steps were requested */
    else {
#ifdef USE_PYTHON
	if ( msr->param.achScriptFile[0] ) {
	    PPY ppy;
	    ppyInitialize(&ppy,msr,dTime);
	    msr->prm->script_argv[0] = msr->param.achScriptFile;
	    ppyRunScript(ppy,msr->prm->script_argc,msr->prm->script_argv);
	    ppyFinish(ppy);
	    }
	else {

#endif
	    if (msrDoGravity(msr) ||msrDoGas(msr)) {
		msrActiveRung(msr,0,1); /* Activate all particles */
		msrDomainDecomp(msr,0,0,0);
		msrUpdateSoft(msr,dTime);
		msrBuildTree(msr,dTime,msr->param.bEwald);
		
		if (msrDoGravity(msr)) {
		    msrGravity(msr,0,MAX_RUNG,ROOT,0,dTime,msr->param.iStartStep,msr->param.bEwald,msr->param.nGroup,&iSec,&nActive);
		    msrMemStatus(msr);
		    if (msr->param.bGravStep) {
			msrBuildTree(msr,dTime,msr->param.bEwald);
			msrGravity(msr,0,MAX_RUNG,ROOT,0,dTime,msr->param.iStartStep,msr->param.bEwald,msr->param.nGroup,&iSec,&nActive);
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
 
		msrCalcEandL(msr,MSR_INIT_E,dTime,&E,&T,&U,&Eth,L,F,&W);
		dMultiEff = 1.0;
		if (msrLogInterval(msr)) {
		    (void) fprintf(fpLog,"%e %e %.16e %e %e %e %.16e %.16e "
			       "%.16e %.16e %.16e %.16e %.16e %i %e\n",dTime,
			       1.0/csmTime2Exp(msr->param.csm,dTime)-1.0,
			       E,T,U,Eth,L[0],L[1],L[2],F[0],F[1],F[2],W,iSec,dMultiEff);
		    }
		}

	    msrOutput(msr,0,dTime,0);  /* JW: Will trash gas density */

#ifdef USE_PYTHON
	    }
#endif

	}

    if (msrLogInterval(msr)) (void) fclose(fpLog);

#ifdef PP_SIMD_BENCHMARK
    PPDumpStats();
#endif

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
