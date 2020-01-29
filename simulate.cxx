/*  This file is part of PKDGRAV3 (http://www.pkdgrav.org/).
 *  Copyright (c) 2001-2020 Douglas Potter & Joachim Stadel
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
#include "pkd_config.h"
#include <string>
#include <cmath>
#include "master.h"


/******************************************************************************\
*   Simulation Mode: normal method of populating the simulation data
\******************************************************************************/
double msrLoadOrGenerateIC(MSR msr) {
    double dTime = -HUGE_VAL;
    if (prmSpecified(msr->prm,"nGrid")) {
	dTime = msrGenerateIC(msr); /* May change nSteps/dDelta */
	if ( msr->param.bWriteIC ) {
	    char achFile[256];			/*DEBUG use MAXPATHLEN here (& elsewhere)? -- DCR*/
	    msrBuildIoName(msr,achFile,0);
	    msrWrite( msr,achFile,dTime,msr->param.bWriteIC-1);
	    }
	}

    /* Read in a binary file */
    else if ( msr->param.achInFile[0] ) {
	dTime = msrRead(msr,msr->param.achInFile); /* May change nSteps/dDelta */
	if (msr->param.bAddDelete) msrGetNParts(msr);
	if (prmSpecified(msr->prm,"dRedFrom")) {
	    double aOld, aNew;
	    aOld = csmTime2Exp(msr->csm,dTime);
	    aNew = 1.0 / (1.0 + msr->param.dRedFrom);
	    dTime = msrAdjustTime(msr,aOld,aNew);
	    /* Seriously, we shouldn't need to send parameters *again*.
	    When we remove sending parameters, we should remove this. */
	    msrSetParameters(msr);
	    }
	}
    else {
	printf("No input file specified\n");
	}
    return dTime;
    }

/******************************************************************************\
*   Simulation Mode: normal operation mode for pkdgrav
\******************************************************************************/
void msrSimulate(MSR msr,double dTime) {
    return msrSimulate(msr,dTime,msr->param.iStartStep);
    }
void msrSimulate(MSR msr,double dTime,int iStartStep) {
    FILE *fpLog = NULL;

    msrSetParameters(msr);
    msrInitCosmology(msr);
    if (prmSpecified(msr->prm,"dSoft")) msrSetSoft(msr,msrSoft(msr));
    if (msrComove(msr)) msrSwitchTheta(msr,dTime); // Adjust theta for gravity calculations.
    /*
    ** Now we have all the parameters for the simulation we can make a
    ** log file entry.
    */
    if (msrLogInterval(msr)) {
	char achFile[256];			/*DEBUG use MAXPATHLEN here (& elsewhere)? -- DCR*/
    	sprintf(achFile,"%s.log",msrOutName(msr));
	fpLog = fopen(achFile,"a");
	assert(fpLog != NULL);
	setbuf(fpLog,(char *) NULL); /* no buffering */
	// fprintf(fpLog,"# ");
	// for (auto i=0;i<argc;++i) fprintf(fpLog,"%s ",argv[i]);
	// fprintf(fpLog,"\n");
	// msrLogParams(msr,fpLog);
	}

    if (msr->param.bLightCone && msrComove(msr)) {
	printf("One, Two, Three replica depth is z=%.10g, %.10g, %.10g\n",
	    1.0/csmComoveLookbackTime2Exp(msr->csm,1.0 / dLightSpeedSim(1*msr->param.dBoxSize)) - 1.0,
	    1.0/csmComoveLookbackTime2Exp(msr->csm,1.0 / dLightSpeedSim(2*msr->param.dBoxSize)) - 1.0,
	    1.0/csmComoveLookbackTime2Exp(msr->csm,1.0 / dLightSpeedSim(3*msr->param.dBoxSize)) - 1.0 );
	}

    /*
    ** Build tree, activating all particles first (just in case).
    */
    msrInflate(msr,iStartStep);
    msrActiveRung(msr,0,1); /* Activate all particles */
    msrDomainDecomp(msr,0,0);
    msrUpdateSoft(msr,dTime);
    msrBuildTree(msr,dTime,msr->param.bEwald);
    msrOutputOrbits(msr,iStartStep,dTime);
    if (msr->param.nGridPk>0) msrOutputPk(msr,iStartStep,dTime);


    int bKickOpen, bKickClose;
    uint8_t uRungMax;
    uint64_t nActive;
    int iSec;

    if (msrDoGravity(msr)) {
	msrSwitchDelta(msr,dTime,iStartStep);
	msrSetParameters(msr);
	if (msr->param.bNewKDK) {
	    msrLightConeOpen(msr,iStartStep + 1);
	    bKickOpen = 1;
	    }
	else bKickOpen = 0;

        /* Compute the grids of the linear species before doing gravity */
        if (strlen(msr->param.achLinearSpecies) && msr->param.nGridLin > 0){
	    msrGridCreateFFT(msr,msr->param.nGridLin);
            msrSetLinGrid(msr,dTime, msr->param.nGridLin,bKickClose,bKickOpen);
            if (msr->param.bDoLinPkOutput)
                msrOutputLinPk(msr, iStartStep, dTime);
	    msrLinearKick(msr,dTime,bKickClose,bKickOpen);
	    msrGridDeleteFFT(msr);
        }
	uRungMax = msrGravity(msr,0,MAX_RUNG,ROOT,0,dTime,iStartStep,0,bKickOpen,
	        msr->param.bEwald,msr->param.bGravStep,msr->param.nPartRhoLoc,msr->param.iTimeStepCrit,msr->param.nGroup,&iSec,&nActive);
	msrMemStatus(msr);
	if (msr->param.bGravStep) {
	    assert(msr->param.bNewKDK == 0);    /* for now! */
	    msrBuildTree(msr,dTime,msr->param.bEwald);
	    msrGravity(msr,0,MAX_RUNG,ROOT,0,dTime,iStartStep,0,0,
		    msr->param.bEwald,msr->param.bGravStep,msr->param.nPartRhoLoc,msr->param.iTimeStepCrit,msr->param.nGroup,&iSec,&nActive);
	    msrMemStatus(msr);
	    }
	}
    if (msrDoGas(msr)) {
	/* Initialize SPH, Cooling and SF/FB and gas time step */
	msrCoolSetup(msr,dTime);
	/* Fix dTuFac conversion of T in InitSPH */
	msrInitSph(msr,dTime);
	}

    double E=0,T=0,U=0,Eth=0,L[3]={0,0,0},F[3]={0,0,0},W=0,dMultiEff = 1.0;
    msrCalcEandL(msr,MSR_INIT_E,dTime,&E,&T,&U,&Eth,L,F,&W);
    if (msrLogInterval(msr)) {
	(void) fprintf(fpLog,"%e %e %.16e %e %e %e %.16e %.16e %.16e "
			"%.16e %.16e %.16e %.16e %i %e\n",dTime,
			1.0/csmTime2Exp(msr->csm,dTime)-1.0,
			E,T,U,Eth,L[0],L[1],L[2],F[0],F[1],F[2],W,iSec,dMultiEff);
	}
    if ( msr->param.bTraceRelaxation) {
	msrInitRelaxation(msr);
	}

    bKickOpen = 0;
    int iStop=0, bDoCheckpoint=0, bDoOutput=0;
    for (auto iStep=iStartStep+1;iStep<=msrSteps(msr)&&!iStop;++iStep) {
	msrSwitchDelta(msr,dTime,iStep-1);
	msrSetParameters(msr);
	if (msrComove(msr)) msrSwitchTheta(msr,dTime);
	dMultiEff = 0.0;
	msr->lPrior = time(0);
	if (msr->param.bNewKDK) {
	    double diStep = (double)(iStep-1);
	    double ddTime = dTime;
	    if (bKickOpen) {
		msrBuildTree(msr,dTime,0);
                msrLightConeOpen(msr,iStep);  /* open the lightcone */
		uRungMax = msrGravity(msr,0,MAX_RUNG,ROOT,0,ddTime,diStep,0,1,
		        msr->param.bEwald,msr->param.bGravStep,msr->param.nPartRhoLoc,msr->param.iTimeStepCrit,msr->param.nGroup,&iSec,&nActive);
                /* Set the grids of the linear species */
                if (strlen(msr->param.achLinearSpecies) && msr->param.nGridLin > 0){
		    msrGridCreateFFT(msr,msr->param.nGridLin);
		    msrSetLinGrid(msr, dTime, msr->param.nGridLin,bKickClose,bKickOpen);
                    if (msr->param.bDoLinPkOutput)
                        msrOutputLinPk(msr, iStartStep, dTime);
		    msrLinearKick(msr,dTime,bKickClose,bKickOpen);
		    msrGridDeleteFFT(msr);
                    }
		bKickOpen = 0; /* clear the opening kicking flag */
		}
	    msrNewTopStepKDK(msr,0,0,&diStep,&ddTime,&uRungMax,&iSec,&bDoCheckpoint,&bDoOutput,&bKickOpen);
	    }
	else {
	    msrTopStepKDK(msr,iStep-1,dTime,msrDelta(msr),0,0,1,&dMultiEff,&iSec);
	    }
	dTime += msrDelta(msr);
	auto lSec = time(0) - msr->lPrior;
	msrMemStatus(msr);

	msrOutputOrbits(msr,iStep,dTime);

	/*
	** Output a log file line if requested.
	** Note: no extra gravity calculation required.
	*/
	if (msrLogInterval(msr) && iStep%msrLogInterval(msr) == 0) {
	    msrCalcEandL(msr,MSR_STEP_E,dTime,&E,&T,&U,&Eth,L,F,&W);
	    (void) fprintf(fpLog,"%e %e %.16e %e %e %e %.16e %.16e "
			    "%.16e %.16e %.16e %.16e %.16e %li %e\n",dTime,
			    1.0/csmTime2Exp(msr->csm,dTime)-1.0,
			    E,T,U,Eth,L[0],L[1],L[2],F[0],F[1],F[2],W,lSec,dMultiEff);
	    }
	if ( msr->param.bTraceRelaxation) {
	    msrActiveRung(msr,0,1); /* Activate all particles */
	    msrDomainDecomp(msr,0,0);
	    msrBuildTree(msr,dTime,0);
	    msrRelaxation(msr,dTime,msrDelta(msr),SMX_RELAXATION,0);
	    }
	if (!msr->param.bNewKDK) {
	    msrCheckForOutput(msr,iStep,dTime,&bDoCheckpoint,&bDoOutput);
	    }
	iStop = (bDoCheckpoint&2) || (bDoOutput&2);
	if (bDoCheckpoint) {
	    msrCheckpoint(msr,iStep,dTime);
	    bDoCheckpoint = 0;
	    }
	if (bDoOutput) {
	    msrOutput(msr,iStep,dTime,0);
	    bDoOutput = 0;
	    if (msr->param.bNewKDK) {
		msrDomainDecomp(msr,0,0);
		msrBuildTree(msr,dTime,msr->param.bEwald);
		}
	    }
	}
    if (msrLogInterval(msr)) (void) fclose(fpLog);

    }

/******************************************************************************\
*   Parameter validation
\******************************************************************************/

#define MAX_CSM_SPECIES 20
static int parseSpeciesNames(const char *aSpecies[], char *achSpecies) {
    if (achSpecies==NULL || achSpecies[0]==0 ) return 0;
    char *p, *stringp = achSpecies;
    int nSpecies = 0;
    while ((p = strsep(&stringp, "+")) != NULL) {
	assert(nSpecies<MAX_CSM_SPECIES);
        if (p[0]) aSpecies[nSpecies++] = p;
	}
    return nSpecies;
    }

/*
** This routine validates the given parameters and makes any adjustments.
*/
int msrValidateParameters(MSR msr) {
    auto mdl = msr->mdl;
    auto csm = msr->csm;
    auto prm = msr->prm;
    auto param = &msr->param;
    
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
	if (param->bPhysicalSoft && !csm->val.bComove) {
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
    if (param->bMemIntegerPosition && (!param->bPeriodic||param->dxPeriod!=1.0||param->dyPeriod!=1.0||param->dzPeriod!=1.0)) {
	fprintf(stderr,"WARNING: Integer coordinates are enabled but the the box is not periodic\n"
	               "       and/or the box size is not 1. Set bPeriodic=1 and dPeriod=1.\n");
	}

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

    /* Make sure that parallel read and write are sane */
    int nThreads = mdlThreads(mdl);
    if (param->nParaRead  > nThreads) param->nParaRead  = nThreads;
    if (param->nParaWrite > nThreads) param->nParaWrite = nThreads;


    /**********************************************************************\
    * The following "parameters" are derived from real parameters.
    \**********************************************************************/

    msr->dTuFac = param->dGasConst/(param->dConstGamma - 1)/param->dMeanMolWeight;
#define KBOLTZ	1.38e-16     /* bolzman constant in cgs */
#define MHYDR 1.67e-24       /* mass of hydrogen atom in grams */
#define MSOLG 1.99e33        /* solar mass in grams */
#define GCGS 6.67e-8         /* G in cgs */
#define KPCCM 3.085678e21    /* kiloparsec in centimeters */
#define SIGMAT 6.6524e-25    /* Thompson cross-section (cm^2) */
#define LIGHTSPEED 2.9979e10 /* Speed of Light cm/s */
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
	}

    /* Determine current opening angle  */
    msr->dThetaMin = msr->param.dTheta;
    if ( !prmSpecified(msr->prm,"nReplicas") && msr->param.nReplicas>=1 ) {
	if ( msr->dThetaMin < 0.52 ) msr->param.nReplicas = 2;
	else msr->param.nReplicas = 1;
	}

    if (msr->csm->val.classData.bClass){
	const char *aLinear[MAX_CSM_SPECIES];
	const char *aPower[MAX_CSM_SPECIES];
	char *achLinearSpecies = strdup(msr->param.achLinearSpecies);
	char *achPowerSpecies = strdup(msr->param.achPowerSpecies);
	int nLinear = parseSpeciesNames(aLinear,achLinearSpecies);
	int nPower = parseSpeciesNames(aPower,achPowerSpecies);
        if (!prmSpecified(msr->prm,"dOmega0")) msr->csm->val.dOmega0 = 0.0;
        csmClassRead(msr->csm, msr->param.achClassFilename, msr->param.dBoxSize, msr->param.h, nLinear, aLinear, nPower, aPower);
        free(achLinearSpecies);
        free(achPowerSpecies);
        csmClassGslInitialize(msr->csm);
	}
    if (strlen(msr->param.achLinearSpecies) && msr->param.nGridLin == 0){
        fprintf(stderr, "ERROR: you must specify nGridLin when running with linear species\n");
        abort();
	}
    return 1;
    }
