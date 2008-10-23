#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
const char *pkdgrav2_module_id = PACKAGE_STRING
    " ($Id$)";

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
#ifdef USE_MDL_IO
#include "io.h"
#endif

#ifdef USE_PYTHON
#include "python.h"
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
	for ( id=1; id<mdlIO(mdl); ++id ) {
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
    double dMultiEff=0;
    long lSec=0,lStart;
    int i,iStep,iSec=0,iStop=0;
    uint64_t nActive;
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
    if (first) {
	save_argv = malloc((argc+1)*sizeof(*save_argv));
	for (i = 0; i < argc; i++)
	    save_argv[i] = strdup(argv[i]);
	save_argv[argc] = NULL;
	}
    else {
	argv = save_argv;
	}
    first = 0;
#endif /* TINY_PTHREAD_STACK */
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

    lStart=time(0);
#ifdef USE_MDL_IO
    mdlInitialize(&mdl,argv,main_ch,main_io);
#else
    mdlInitialize(&mdl,argv,main_ch,0);
#endif
    for (argc = 0; argv[argc]; argc++); /* some MDLs can trash argv */

    printf("%s\n", PACKAGE_STRING );

    msrInitialize(&msr,mdl,argc,argv);

    /*
    ** Establish safety lock.
    */
    if (!msrGetLock(msr)) {
	msrFinish(msr);
	mdlFinish(mdl);
	return 1;
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
	    if (prmSpecified(msr->prm,"dSoft")) msrSetSoft(msr,msrSoft(msr));
	    }
	else {
#ifdef USE_PYTHON
	    if ( !msr->param.achScriptFile[0] ) {
#endif
		printf("No input file specified\n");
		return 1;
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
	msrDomainDecomp(msr,0,1,0);
	msrUpdateSoft(msr,dTime);
	msrBuildTree(msr,dTime,msr->param.bEwald);
	if (msrDoGravity(msr)) {
	    msrGravity(msr,0,MAX_RUNG,dTime,msr->param.iStartStep,msr->param.bEwald,&iSec,&nActive);
	    msrMemStatus(msr);
	    if (msr->param.bGravStep) {
		msrBuildTree(msr,dTime,msr->param.bEwald);
		msrGravity(msr,0,MAX_RUNG,dTime,msr->param.iStartStep,msr->param.bEwald,&iSec,&nActive);
		msrMemStatus(msr);
		}
	    msrCalcEandL(msr,MSR_INIT_E,dTime,&E,&T,&U,&Eth,L);
	    dMultiEff = 1.0;
	    if (msrLogInterval(msr)) {
		(void) fprintf(fpLog,"%e %e %.16e %e %e %e %.16e %.16e %.16e "
			       "%i %e\n",dTime,
			       1.0/csmTime2Exp(msr->param.csm,dTime)-1.0,
			       E,T,U,Eth,L[0],L[1],L[2],iSec,dMultiEff);
		}
	    }
#ifdef PLANETS
	if (msr->param.bHeliocentric) {
	    msrGravSun(msr);
	    }
#ifdef SYMBA
	msrDriftSun(msr,dTime,0.5*msrDelta(msr));
#endif
#endif
	if ( msr->param.bTraceRelaxation) {
	    msrInitRelaxation(msr);
	    }
#ifdef HERMITE
	if (msr->param.bHermite) {
	    msrActiveRung(msr,0,1); /* Activate all particles */
	    msrCopy0(msr, dTime);
	    if (msr->param.bAarsethStep) {
		msrFirstDt(msr);
		}
	    }
#endif
	for (iStep=msr->param.iStartStep+1;iStep<=msrSteps(msr)&&!iStop;++iStep) {
	    if (msrComove(msr)) {
		msrSwitchTheta(msr,dTime);
		}
	    dMultiEff = 0.0;
	    lSec = time(0);
#ifdef HERMITE
	    if (msr->param.bHermite) {
		msrTopStepHermite(msr,iStep-1,dTime,
				  msrDelta(msr),0,0,msrMaxRung(msr),1,
				  &dMultiEff,&iSec);
		}
	    else
#endif
#ifdef SYMBA
		if (msr->param.bSymba) {
		    msrTopStepSymba(msr,iStep-1,dTime,
				    msrDelta(msr),0,0,msrMaxRung(msr),1,
				    &dMultiEff,&iSec);
		    }
		else
#endif

		    {
		    msrTopStepKDK(msr,iStep-1,dTime,
				  msrDelta(msr),0,0,msrMaxRung(msr),1,
				  &dMultiEff,&iSec);
		    }

	    dTime += msrDelta(msr);

	    msrMemStatus(msr);

	    /*
	    ** Output a log file line if requested.
	    ** Note: no extra gravity calculation required.
	    */
	    if (msrLogInterval(msr) && iStep%msrLogInterval(msr) == 0) {
		msrCalcEandL(msr,MSR_STEP_E,dTime,&E,&T,&U,&Eth,L);
		lSec = time(0) - lSec;
		(void) fprintf(fpLog,"%e %e %.16e %e %e %e %.16e %.16e "
			       "%.16e %li %e\n",dTime,
			       1.0/csmTime2Exp(msr->param.csm,dTime)-1.0,
			       E,T,U,Eth,L[0],L[1],L[2],lSec,dMultiEff);
		}
	    if ( msr->param.bTraceRelaxation) {
		msrActiveRung(msr,0,1); /* Activate all particles */
		msrDomainDecomp(msr,0,1,0);
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
	    ppyRunScript(ppy,msr->param.achScriptFile);
	    ppyFinish(ppy);
	    }
	else {

#endif
	    if (msrDoGravity(msr)) {
		msrActiveRung(msr,0,1); /* Activate all particles */
		msrDomainDecomp(msr,0,1,0);
		msrUpdateSoft(msr,dTime);
		msrBuildTree(msr,dTime,msr->param.bEwald);

		msrGravity(msr,0,MAX_RUNG,dTime,msr->param.iStartStep,msr->param.bEwald,&iSec,&nActive);
		if (msr->param.bGravStep && msr->param.iTimeStepCrit == -1) {
		    msrGravity(msr,0,MAX_RUNG,dTime,msr->param.iStartStep,msr->param.bEwald,&iSec,&nActive);
		    }

		msrCalcEandL(msr,MSR_INIT_E,dTime,&E,&T,&U,&Eth,L);
		dMultiEff = 1.0;
		if (msrLogInterval(msr)) {
			(void) fprintf(fpLog,"%e %e %.16e %e %e %e %.16e %.16e %.16e "
				       "%i %e\n",dTime,
				       1.0/csmTime2Exp(msr->param.csm,dTime)-1.0,
				       E,T,U,Eth,L[0],L[1],L[2],iSec,dMultiEff);
		    }
		}

	    msrOutput(msr,0,dTime,0);
#ifdef USE_PYTHON
	    }
#endif

	}

    if (msrLogInterval(msr)) (void) fclose(fpLog);

#ifdef PP_SIMD_BENCHMARK
    PPDumpStats();
#endif

    msrFinish(msr);
    mdlFinish(mdl);
    return 0;
    }
