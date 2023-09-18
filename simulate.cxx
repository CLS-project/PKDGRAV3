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
#include "units.h"
#include "fmt/format.h"
#include "SPH/SPHOptions.h"

/******************************************************************************\
*   Simulation Mode: normal method of populating the simulation data
\******************************************************************************/
double MSR::LoadOrGenerateIC() {
    double dTime = -HUGE_VAL;
    if (parameters.has_nGrid()) {
        dTime = GenerateIC(
                    parameters.get_nGrid(),
                    parameters.get_iSeed(),
                    parameters.get_dRedFrom(),
                    parameters.get_dBoxSize(),
                    csm); /* May change nSteps/dDelta */
        if ( param.bWriteIC ) {
            Write(BuildIoName(0).c_str(),dTime,param.bWriteIC-1);
        }
    }

    /* Read in a binary file */
    else if ( parameters.has_achInFile() ) {
        dTime = Read(parameters.get_achInFile()); /* May change nSteps/dDelta */
    }
    else {
        printf("No input file specified\n");
    }
    if (param.bAddDelete) GetNParts();
    return dTime;
}

/******************************************************************************\
*   Calculate dDelta and nSteps from: iStartStep, nSteps, dRedTo, dDelta
\******************************************************************************/

bool MSR::getDeltaSteps(double dTime,int iStartStep,double &dDelta,int &nSteps) {
    dDelta = param.dDelta;
    nSteps = parameters.get_nSteps();
    if (!csm->val.bComove) {
        if (prmSpecified(prm,"dRedTo")) {
            printf("WARNING: dRedTo is meaningless for non-cosmological simulations, ignoring.\n");
        }
    }
    else if (!prmSpecified(prm,"dRedTo")) {}
    else if (prmSpecified(prm,"dDelta") && parameters.has_nSteps()) {
        printf("Specify at most two of: dDelta, nSteps, dRedTo -- all three were specified\n");
        return false;
    }
    else if (prmSpecified(prm,"dDelta")) {
        auto aTo = 1.0/(param.dRedTo + 1.0);
        auto tTo = csmExp2Time(csm,aTo);
        if (tTo < dTime) {
            printf("Badly specified final redshift, check -zto parameter.\n");
            return false;
        }
        nSteps = (int)ceil((tTo-dTime)/param.dDelta);
        dDelta = (tTo-dTime)/(nSteps - iStartStep);
    }
    else if (parameters.has_nSteps()) {
        auto aTo = 1.0/(param.dRedTo + 1.0);
        auto tTo = csmExp2Time(csm,aTo);
        if (tTo < dTime) {
            printf("Badly specified final redshift, check -zto parameter.\n");
            return false;
        }
        if (parameters.get_nSteps() == 0) dDelta = 0.0;
        else dDelta = (tTo-dTime) / (nSteps - iStartStep);
    }
    return true;
}

/******************************************************************************\
*   Simulation Mode: normal operation mode for pkdgrav
\******************************************************************************/
void MSR::Simulate(double dTime) {
    double dDelta;
    int nSteps;
    const auto iStartStep = parameters.get_iStartStep();
    if (prmSpecified(prm,"achOutTimes")) {
        nSteps = ReadOuts(dTime);
    }
    else {
        getDeltaSteps(dTime,iStartStep, /* OUTPUT -> */ dDelta,nSteps);
    }
    return Simulate(dTime,dDelta,iStartStep,nSteps);
}
void MSR::Simulate(double dTime,double dDelta,int iStartStep,int nSteps, bool bRestart) {
    FILE *fpLog = NULL;
    const auto bEwald = parameters.get_bEwald();
    const auto bGravStep = parameters.get_bGravStep();
    const auto nPartRhoLoc = parameters.get_nPartRhoLoc();
    const auto iTimeStepCrit = parameters.get_iTimeStepCrit();

    InitCosmology(csm);
    auto dTheta = set_dynamic(iStartStep,dTime);

    if (parameters.has_dSoft()) SetSoft(Soft());

    /*
    ** Now we have all the parameters for the simulation we can make a
    ** log file entry.
    */
    if (LogInterval()) {
        std::string filename = std::string(OutName()) + ".log";
        fpLog = fopen(filename.c_str(),"a");
        assert(fpLog != NULL);
        setbuf(fpLog,(char *) NULL); /* no buffering */
        // fprintf(fpLog,"# ");
        // for (auto i=0;i<argc;++i) fprintf(fpLog,"%s ",argv[i]);
        // fprintf(fpLog,"\n");
        // msrLogParams(msr,fpLog);
    }

    TimerHeader();

    if (param.bLightCone && Comove()) {
        auto dBoxSize = parameters.get_dBoxSize();
        printf("One, Two, Three replica depth is z=%.10g, %.10g, %.10g\n",
               1.0/csmComoveLookbackTime2Exp(csm,1.0 / dLightSpeedSim(1*dBoxSize)) - 1.0,
               1.0/csmComoveLookbackTime2Exp(csm,1.0 / dLightSpeedSim(2*dBoxSize)) - 1.0,
               1.0/csmComoveLookbackTime2Exp(csm,1.0 / dLightSpeedSim(3*dBoxSize)) - 1.0 );
    }

    if (DoGas() && MeshlessHydro()) {
        ChemCompInit();
    }

#ifdef COOLING
    float redshift;
    if (csm->val.bComove)
        redshift = 1./csmTime2Exp(csm,dTime) - 1.;
    else
        redshift = 0.0;

    CoolingInit(redshift);
    CoolingUpdate(redshift, 1);
#endif

#ifdef GRACKLE
    if ((csm->val.bComove)) {
        GrackleInit(1, csmTime2Exp(csm,dTime));
    }
    else {
        GrackleInit(0, 1.0);
    }
#endif

#if defined(STAR_FORMATION) || defined(FEEDBACK)
    StarFormInit(dTime);
#ifdef STAR_FORMATION
    starFormed = 0.;
    massFormed = 0.;
#endif
#endif

#ifdef STELLAR_EVOLUTION
    StellarEvolutionInit(dTime);
#endif

    OutputFineStatistics(0.0, -1);
    /*
    ** Build tree, activating all particles first (just in case).
    */
    ActiveRung(0,1); /* Activate all particles */
    DomainDecomp();
    UpdateSoft(dTime);
    BuildTree(bEwald);
    runAnalysis(iStartStep,dTime); // Run any registered Python analysis tasks
    OutputOrbits(iStartStep,dTime);
    if (param.nGridPk>0) OutputPk(iStartStep,dTime);

    int bKickOpen, bKickClose=0;
    uint8_t uRungMax;
    int iSec = time(0);

    dDelta = SwitchDelta(dTime,dDelta,iStartStep,parameters.get_nSteps());
    if (param.bNewKDK) {
        LightConeOpen(iStartStep + 1);
        bKickOpen = 1;
    }
    else bKickOpen = 0;
    if (DoGravity()) {
        /* Compute the grids of the linear species before doing gravity */
        if (csm->val.classData.bClass && strlen(param.achLinSpecies) && param.nGridLin > 0) {
            GridCreateFFT(param.nGridLin);
            SetLinGrid(dTime,dDelta,param.nGridLin,bKickClose,bKickOpen);
            if (param.bDoLinPkOutput)
                OutputLinPk( iStartStep, dTime);
            LinearKick(dTime,dDelta,bKickClose,bKickOpen);
            GridDeleteFFT();
        }
    }
    if (DoGas() && NewSPH()) {
        // Calculate Density
        SelAll(-1,1);
        SPHOptions SPHoptions = initializeSPHOptions(parameters,param,csm,dTime);
        SPHoptions.doGravity = 0;
        SPHoptions.doDensity = 1;
        SPHoptions.doSPHForces = 0;
        uRungMax = Gravity(0,MAX_RUNG,ROOT,0,dTime,dDelta,iStartStep,dTheta,0,bKickOpen,
                           bEwald,bGravStep,nPartRhoLoc,iTimeStepCrit,SPHoptions);
        MemStatus();
        if (SPHoptions.doInterfaceCorrection) {
            SPHoptions.doDensity = 0;
            SPHoptions.doDensityCorrection = 1;
            SPHoptions.dofBallFactor = 0;
            TreeUpdateFlagBounds(bEwald,ROOT,0,SPHoptions);
            uRungMax = Gravity(0,MAX_RUNG,ROOT,0,dTime,dDelta,iStartStep,dTheta,0,bKickOpen,
                               bEwald,bGravStep,nPartRhoLoc,iTimeStepCrit,SPHoptions);
            UpdateGasValues(0,dTime,dDelta,iStartStep,0,bKickOpen,SPHoptions);
            SPHoptions.doDensityCorrection = 0;
        }
        // Calculate Forces
        SelAll(-1,1);
        SPHoptions.doGravity = param.bDoGravity;
        SPHoptions.doDensity = 0;
        SPHoptions.doSPHForces = 1;
        SPHoptions.dofBallFactor = 0;
        TreeUpdateFlagBounds(bEwald,ROOT,0,SPHoptions);
        uRungMax = Gravity(0,MAX_RUNG,ROOT,0,dTime,dDelta,iStartStep,dTheta,0,bKickOpen,
                           bEwald,bGravStep,nPartRhoLoc,iTimeStepCrit,SPHoptions);
        MemStatus();
    }
    else if (DoGravity()) {
        SPHOptions SPHoptions = initializeSPHOptions(parameters,param,csm,dTime);
        SPHoptions.doGravity = param.bDoGravity;
        uRungMax = Gravity(0,MAX_RUNG,ROOT,0,dTime,dDelta,iStartStep,dTheta,0,bKickOpen,
                           bEwald,bGravStep,nPartRhoLoc,iTimeStepCrit,SPHoptions);
        MemStatus();
    }
    if (DoGravity() && bGravStep) {
        assert(param.bNewKDK == 0);    /* for now! */
        BuildTree(bEwald);
        SPHOptions SPHoptions = initializeSPHOptions(parameters,param,csm,dTime);
        SPHoptions.doGravity = param.bDoGravity;
        Gravity(0,MAX_RUNG,ROOT,0,dTime,dDelta,iStartStep,dTheta,0,0,
                bEwald,bGravStep,nPartRhoLoc,iTimeStepCrit,SPHoptions);
        MemStatus();
    }
    if (DoGas() && MeshlessHydro()) {
        InitSph(dTime, dDelta,bRestart);
    }
#ifdef BLACKHOLES
    uRungMax = GetMinDt();
#ifndef DEBUG_BH_ONLY
    BlackholeInit(uRungMax);
#endif
    if (parameters.get_bFindGroups() && param.bBHPlaceSeed) {
        NewFof(parameters.get_dTau(),parameters.get_nMinMembers());
        GroupStats();
    }
#endif

    double E=0,T=0,U=0,Eth=0,L[3]= {0,0,0},F[3]= {0,0,0},W=0;
    CalcEandL(MSR_INIT_E,dTime,&E,&T,&U,&Eth,L,F,&W);
    iSec = time(0) - iSec;
    if (LogInterval()) {
        (void) fprintf(fpLog,"%e %e %.16e %e %e %e %.16e %.16e %.16e "
                       "%.16e %.16e %.16e %.16e %i\n",dTime,
                       1.0/csmTime2Exp(csm,dTime)-1.0,
                       E,T,U,Eth,L[0],L[1],L[2],F[0],F[1],F[2],W,iSec);
    }

    if (param.bWriteIC && !parameters.has_nGrid() && !NewSPH()) {
#ifndef BLACKHOLES
        if (parameters.get_bFindGroups()) {
            NewFof(parameters.get_dTau(),parameters.get_nMinMembers());
            GroupStats();
        }
#endif
        Output(iStartStep,dTime,dDelta,0);
    }

    // Make sure that the tree is usable before the start of the simulation
    if (parameters.get_bFindGroups() || param.bWriteIC) {
        DomainDecomp();
        BuildTree(bEwald);
    }

    bKickOpen = 0;
    int iStop=0, bDoCheckpoint=0, bDoOutput=0;
    for (auto iStep=iStartStep+1; iStep<=nSteps&&!iStop; ++iStep) {
        auto dTheta = set_dynamic(iStep-1,dTime);
        dDelta = SwitchDelta(dTime,dDelta,iStep-1,parameters.get_nSteps());
        lPrior = time(0);
        TimerRestart();
        if (param.bNewKDK) {
            double diStep = (double)(iStep-1);
            double ddTime = dTime;
            if (bKickOpen) {
                BuildTree(0);
                LightConeOpen(iStep);  /* open the lightcone */
                if (DoGas() && NewSPH()) {
                    // Calculate Density
                    SelAll(-1,1);
                    SPHOptions SPHoptions = initializeSPHOptions(parameters,param,csm,dTime);
                    SPHoptions.doGravity = 0;
                    SPHoptions.doDensity = 1;
                    SPHoptions.doSPHForces = 0;
                    uRungMax = Gravity(0,MAX_RUNG,ROOT,0,ddTime,dDelta,diStep,dTheta,0,1,
                                       bEwald,bGravStep,nPartRhoLoc,iTimeStepCrit,SPHoptions);
                    if (SPHoptions.doInterfaceCorrection) {
                        SPHoptions.doDensity = 0;
                        SPHoptions.doDensityCorrection = 1;
                        SPHoptions.dofBallFactor = 0;
                        TreeUpdateFlagBounds(bEwald,ROOT,0,SPHoptions);
                        uRungMax = Gravity(0,MAX_RUNG,ROOT,0,ddTime,dDelta,diStep,dTheta,0,1,
                                           bEwald,bGravStep,nPartRhoLoc,iTimeStepCrit,SPHoptions);
                        UpdateGasValues(0,ddTime,dDelta,diStep,0,1,SPHoptions);
                        SPHoptions.doDensityCorrection = 0;
                    }
                    // Calculate Forces
                    SelAll(-1,1);
                    SPHoptions.doGravity = param.bDoGravity;
                    SPHoptions.doDensity = 0;
                    SPHoptions.doSPHForces = 1;
                    SPHoptions.dofBallFactor = 0;
                    TreeUpdateFlagBounds(bEwald,ROOT,0,SPHoptions);
                    uRungMax = Gravity(0,MAX_RUNG,ROOT,0,ddTime,dDelta,diStep,dTheta,0,1,
                                       bEwald,bGravStep,nPartRhoLoc,iTimeStepCrit,SPHoptions);
                }
                else {
                    SPHOptions SPHoptions = initializeSPHOptions(parameters,param,csm,dTime);
                    SPHoptions.doGravity = param.bDoGravity;
                    uRungMax = Gravity(0,MAX_RUNG,ROOT,0,ddTime,dDelta,diStep,dTheta,0,1,
                                       bEwald,bGravStep,nPartRhoLoc,iTimeStepCrit,SPHoptions);
                }
                /* Set the grids of the linear species */
                if (csm->val.classData.bClass && strlen(param.achLinSpecies) && param.nGridLin > 0) {
                    GridCreateFFT(param.nGridLin);
                    SetLinGrid(dTime,dDelta,param.nGridLin,bKickClose,bKickOpen);
                    if (param.bDoLinPkOutput)
                        OutputLinPk( iStartStep, dTime);
                    LinearKick(dTime,dDelta,bKickClose,bKickOpen);
                    GridDeleteFFT();
                }
                bKickOpen = 0; /* clear the opening kicking flag */
            }
            NewTopStepKDK(ddTime,dDelta,dTheta,nSteps,0,0,&diStep,&uRungMax,&bDoCheckpoint,&bDoOutput,&bKickOpen);
        }
        else {
            TopStepKDK(iStep-1,dTime,dDelta,dTheta,0,0,1);
            runAnalysis(iStep,dTime); // Run any registered Python analysis tasks
        }
        dTime += dDelta;
        auto lSec = time(0) - lPrior;
        MemStatus();

        OutputOrbits(iStep,dTime);

        /*
        ** Output a log file line if requested.
        ** Note: no extra gravity calculation required.
        */
        if (LogInterval() && iStep%LogInterval() == 0) {
            CalcEandL(MSR_STEP_E,dTime,&E,&T,&U,&Eth,L,F,&W);
            (void) fprintf(fpLog,"%e %e %.16e %e %e %e %.16e %.16e "
                           "%.16e %.16e %.16e %.16e %.16e %li\n",dTime,
                           1.0/csmTime2Exp(csm,dTime)-1.0,
                           E,T,U,Eth,L[0],L[1],L[2],F[0],F[1],F[2],W,lSec);
        }
        if (!param.bNewKDK) {
            CheckForOutput(iStep,nSteps,dTime,&bDoCheckpoint,&bDoOutput);
        }
        iStop = (bDoCheckpoint&2) || (bDoOutput&2);
        if (bDoCheckpoint) {
            Checkpoint(iStep,nSteps,dTime,dDelta);
            bDoCheckpoint = 0;
        }
        if (bDoOutput) {
            if (DoGas() && NewSPH() && param.bCentrifugal) {
                ResetCOM();
            }
            Output(iStep,dTime,param.dDelta,0);
            bDoOutput = 0;
            DomainDecomp();
            BuildTree(bEwald);
        }
        TimerDump(iStep);
    }
    if (LogInterval()) (void) fclose(fpLog);

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

// Check to see if a path specification is valid. Paths are constructed using
// the Python/Format string formats and may contain the following fields.
//   {name} - replaced with the name of the simulation
//   {step} - replaced with the current step number
//   {type} - replaced with the file type, including the leading dot
// The default is:
//   {name}.{step:05d}{type}
static void validate_path(const char *name,const char *path) {
    if (path && path[0]) {
        using namespace fmt::literals;
        auto r1 = fmt::format(path,"name"_a="name","step"_a=1,"type"_a=".dat");
        if (r1 == fmt::format(path,"name"_a="name","step"_a=1,"type"_a=".XXX"))
            throw fmt::format_error(std::string(name) + ": " + path);
        if (r1 == fmt::format(path,"name"_a="name","step"_a=9,"type"_a=".dat"))
            throw fmt::format_error(std::string(name) + ": " + path);
        if (r1 == fmt::format(path,"name"_a="XXXX","step"_a=1,"type"_a=".dat"))
            throw fmt::format_error(std::string(name) + ": " + path);
    }
}

/*
** This routine validates the given parameters and makes any adjustments.
*/
int MSR::ValidateParameters() {
    try {
        validate_path("achOutPath",       param.achOutPath);
        validate_path("achCheckpointPath",param.achCheckpointPath);
        validate_path("achIoPath",        param.achIoPath);
    }
    catch (fmt::format_error &e) {
        fprintf(stderr,"ERROR: %s\n",e.what());
        fprintf(stderr,"       When specified must contain {name}, {step} and {type} and no other fields\n");
        fprintf(stderr,"       Default: {name}.{step:05d}{type}\n");
        fprintf(stderr,"       Example: /path/to/output/{step:05d}/{name}.{step:05d}{type}\n");
        return 0;
    }

    if (prmSpecified(prm, "achOutTimes") && parameters.has_nSteps() ) {
        fprintf(stderr, "ERROR: achOutTimes and nSteps can not given at the same time.\n");
        return 0;
    }
    if (prmSpecified(prm, "achOutTimes") && prmSpecified(prm, "dRedTo") ) {
        fprintf(stderr, "ERROR: achOutTimes and dRedTo can not be given at the same time.\n");
        fprintf(stderr, "       Add your final redshift (dRedTo) to the end of achOutTimes file.\n");
        return 0;
    }

    if (!parameters.has_bDoGas()) parameters.set_bDoGas(param.bMeshlessHydro||param.bNewSPH);
    if (parameters.get_bDoGas() && !(param.bMeshlessHydro||param.bNewSPH) ) {
        fprintf(stderr,"ERROR: Please provide an hydrodynamic solver to be used: bMeshlessHydro or bNewSPH.\n");
        return 0;
    }
    if ((param.bMeshlessHydro||param.bNewSPH) && !parameters.get_bDoGas()) {
        fprintf(stderr,"ERROR: An hydro scheme is selected but bDoGas=0! Did you forget to add bDoGas=1?\n");
        return 0;
    }
    if (param.bMeshlessHydro && param.bNewSPH) {
        fprintf(stderr,"ERROR: Only one hydrodynamic scheme can be used.\n");
        return 0;
    }

#ifdef BLACKHOLES
    if  (!ValidateBlackholeParam()) return 0;
#endif

#ifdef STAR_FORMATION
    if  (!ValidateStarFormationParam()) return 0;
#endif

    if (param.bGasInterfaceCorrection && param.bGasOnTheFlyPrediction) {
        fprintf(stderr,"Warning: On-the-fly prediction is not compatible with interface correction, disabled\n");
        param.bGasOnTheFlyPrediction = 0;
    }

#ifndef NN_FLAG_IN_PARTICLE
    if (param.bNewSPH && param.bGasInterfaceCorrection && param.dFastGasFraction > 0.0f) {
        fprintf(stderr,"ERROR: Interface correction and FastGas is active, but the NN flag is not compiled in. Set NN_FLAG_IN_PARTICLE to ON in CMakeLists.txt and recompile.\n");
        return 0;
    }
#endif

    /*
    ** CUDA likes a larger group size
    */
    if ( (mdl->isCudaActive() || mdl->isMetalActive()) && !parameters.has_nGroup() && parameters.get_nGroup()<256)
        parameters.set_nGroup(256);

#ifndef USE_HDF5
    if (param.bHDF5) {
        printf("WARNING: HDF5 output was requested by is not supported: using Tipsy format\n");
        param.bHDF5 = 0;
    }
#endif

#ifdef MDL_FFTW
    if ( param.nGridPk ) {
        if (prmSpecified(prm,"nBinsPk")) {
            if (param.nBinsPk > param.nGridPk/2) {
                param.nBinsPk = param.nGridPk/2;
            }
        }
        else param.nBinsPk = param.nGridPk/2;
        if (param.nBinsPk > PST_MAX_K_BINS)
            param.nBinsPk = PST_MAX_K_BINS;
    }
    if (param.iPkOrder<1 || param.iPkOrder>4) {
        puts("ERROR: iPkOrder must be 1 (NGP), 2 (CIC), 3 (TSC) or 4 (PCS)");
        return 0;
    }
    if ( parameters.get_nGrid() ) {
        if (parameters.has_achInFile()) {
            puts("ERROR: do not specify an input file when generating IC");
            return 0;
        }
        if ( parameters.get_iSeed() == 0 ) {
            //puts("ERROR: Random seed for IC not specified");
            parameters.set(parameters.str_iSeed,time(NULL));
        }
        if ( !parameters.has_dBoxSize() || parameters.get_dBoxSize() <= 0 ) {
            puts("ERROR: Box size for IC not specified");
            return 0;
        }
        if ( !csm->val.classData.bClass ) {
            if ( ( !prmSpecified(prm,"dSigma8") || csm->val.dSigma8 <= 0 ) &&
                    ( !prmSpecified(prm,"dNormalization") || csm->val.dNormalization <= 0 ) ) {
                puts("ERROR: Either dSigma8 or dNormalization should be specified for generating IC");
                return 0;
            }
            if ( !prmSpecified(prm,"dSpectral") || csm->val.dSpectral <= 0 ) {
                puts("ERROR: dSpectral for IC not specified");
                return 0;
            }
        }
        if ( param.bICgas ) {
            if ( !prmSpecified(prm,"dOmegab") || csm->val.dOmegab <= 0 ) {
                puts("ERROR: Can not generate IC with gas if dOmegab is not specified");
                return 0;
            }
            if ( !parameters.get_bDoGas() ) {
                puts("ERROR: Can not generate gas if bDoGas=0");
                return 0;
            }
        }
    }
    if ( csm->val.bComove && !csm->val.classData.bClass ) {
        if ( !prmSpecified(prm,"h") ) {
            fprintf(stderr, "WARNING: Running with bComove without specifying a Hubble parameter, h\n");
        }
    }
    /* Set the number of bins for the power spectrum measurement of linear species */
    if (param.nGridLin > 0) {
        param.nBinsLinPk = param.nGridLin/2;
        if (param.nBinsLinPk > PST_MAX_K_BINS)
            param.nBinsLinPk = PST_MAX_K_BINS;
    }
#endif
    if (!prmSpecified(prm,"dFracNoDomainRootFind") && param.dFracNoDomainRootFind > param.dFracNoDomainDimChoice) param.dFracNoDomainRootFind = param.dFracNoDomainDimChoice;
    if (!prmSpecified(prm,"dFracNoDomainDecomp") && param.dFracNoDomainDecomp > param.dFracNoDomainRootFind) param.dFracNoDomainDecomp = param.dFracNoDomainRootFind;
    if (!prmSpecified(prm,"dFracDualTree") && param.dFracDualTree > param.dFracNoDomainDecomp) param.dFracDualTree = param.dFracNoDomainDecomp;
    if (!prmSpecified(prm,"dFracNoDomainDecomp") && param.dFracNoDomainDecomp < param.dFracDualTree) param.dFracNoDomainDecomp = param.dFracDualTree;
    if (!prmSpecified(prm,"dFracNoDomainRootFind") && param.dFracNoDomainRootFind < param.dFracNoDomainDecomp) param.dFracNoDomainRootFind = param.dFracNoDomainDecomp;
    if (!prmSpecified(prm,"dFracNoDomainDimChoice") && param.dFracNoDomainDimChoice < param.dFracNoDomainRootFind) param.dFracNoDomainDimChoice = param.dFracNoDomainRootFind;
    if ( param.dFracDualTree > param.dFracNoDomainDecomp
            || param.dFracNoDomainDecomp > param.dFracNoDomainRootFind
            || param.dFracNoDomainRootFind > param.dFracNoDomainDimChoice
            || param.dFracNoDomainDecomp<0.0 || param.dFracNoDomainDimChoice > 1.0 ) {
        puts("ERROR: check that 0 <= dFracNoDomainDecomp <= dFracNoDomainRootFind <= dFracNoDomainDimChoice <= 1");
        return 0;
    }

    /* Make sure that the old behaviour is obeyed. */
    if ( parameters.has_nSteps() && parameters.get_nSteps() == 0 ) {
        if ( !prmSpecified(prm,"bDoAccOutput") ) param.bDoAccOutput = 1;
        if ( !prmSpecified(prm,"bDoPotOutput") ) param.bDoPotOutput = 1;
    }

    /*
     * Softening
     */
    auto bPhysicalSoft = parameters.get_bPhysicalSoft();
    const auto dMaxPhysicalSoft = parameters.get_dMaxPhysicalSoft();
    if (bPhysicalSoft && !csm->val.bComove) {
        printf("WARNING: bPhysicalSoft reset to 0 for non-comoving (bComove == 0)\n");
        parameters.set_bPhysicalSoft(bPhysicalSoft = false);
    }
    if (bPhysicalSoft && dMaxPhysicalSoft>0) {
        fprintf(stderr, "ERROR: Setting both bPhysicalSoft and dMaxPhysicalSoft "
                "is not allowed.\n Did you mean to limit the physical softening"
                "with bPhysicalSoft and dSoftMax? or just limit the comoving "
                "softening with dMaxPhysicalSoft?\n");
        return 0;
    }
    if ( dMaxPhysicalSoft>0 && parameters.get_dSoft()==0.0 && !parameters.get_bSoftMaxMul()) {
        fprintf(stderr, "ERROR: Trying to limit individual softenings setting a "
                "maximum physical softening rather than a factor...\nThis is "
                "not supported.\n Did you mean to use dSoft for a global softening? "
                "or bSoftMaxMul for setting the limit as a factor?\n");
        return 0;
    }
    if ( bPhysicalSoft && parameters.get_dSoftMax()==0.0) {
        fprintf(stderr, "ERROR: If setting bPhysicalSoft, dSoftMax should be "
                "provided to avoid divergences in the early universe.\n");
        return 0;
    }
    /*
    ** Determine the period of the box that we are using.
    ** Set the new d[xyz]Period parameters which are now used instead
    ** of a single dPeriod, but we still want to have compatibility
    ** with the old method of setting dPeriod.
    */
    if (prmSpecified(prm,"dPeriod") &&
            !prmSpecified(prm,"dxPeriod")) {
        param.dxPeriod = param.dPeriod;
    }
    if (prmSpecified(prm,"dPeriod") &&
            !prmSpecified(prm,"dyPeriod")) {
        param.dyPeriod = param.dPeriod;
    }
    if (prmSpecified(prm,"dPeriod") &&
            !prmSpecified(prm,"dzPeriod")) {
        param.dzPeriod = param.dPeriod;
    }
    /*
    ** Periodic boundary conditions can be disabled along any of the
    ** x,y,z axes by specifying a period of zero for the given axis.
    ** Internally, the period is set to infinity (Cf. pkdBucketWalk()
    ** and pkdDrift(); also the INTERSECT() macro in smooth.h).
    */
    if (param.dPeriod  == 0) param.dPeriod  = FLOAT_MAXVAL;
    if (param.dxPeriod == 0) param.dxPeriod = FLOAT_MAXVAL;
    if (param.dyPeriod == 0) param.dyPeriod = FLOAT_MAXVAL;
    if (param.dzPeriod == 0) param.dzPeriod = FLOAT_MAXVAL;
    /*
    ** At the moment, integer positions are only really safe in periodic boxes!Wr
    */
    if (parameters.get_bMemIntegerPosition() && (!parameters.get_bPeriodic()||param.dxPeriod!=1.0||param.dyPeriod!=1.0||param.dzPeriod!=1.0)) {
        fprintf(stderr,"WARNING: Integer coordinates are enabled but the the box is not periodic\n"
                "       and/or the box size is not 1. Set bPeriodic=1 and dPeriod=1.\n");
    }

    /*
    ** Check if fast gas boundaries are needed.
    */
    if (parameters.get_bDoGas() && !NewSPH()) {
        parameters.set_bMemNodeSphBounds(1);
    }
    /*
    ** Check timestepping and gravity combinations.
    */
    const auto iMaxRung = parameters.get_iMaxRung();
    assert(iMaxRung <= IRUNGMAX);
    if (parameters.get_bEpsAccStep()) param.bAccelStep = 1;
    if (param.bDoGravity) {
        /* Potential is optional, but the default for gravity */
        if (!parameters.has_bMemPotential()) parameters.set_bMemPotential(1);
        if (iMaxRung < 1) {
            parameters.set_iMaxRung(0);
            if (parameters.get_bVWarnings()) fprintf(stderr,"WARNING: iMaxRung set to 0, SINGLE STEPPING run!\n");
            /*
            ** For single stepping we don't need fancy timestepping variables.
            */
            parameters.set_bMemNodeAcceleration(0);
            parameters.set_bMemNodeVelocity(0);
        }
        else {
            if ((param.bAccelStep || parameters.get_bDensityStep()) && parameters.get_bGravStep()) {
                /*
                ** We cannot combine these 2 types of timestepping criteria, we need to choose one
                ** or the other basic timestep criterion, in this case we choose only bGravStep.
                */
                param.bAccelStep = 0;
                parameters.set_bEpsAccStep(false);
                parameters.set_bDensityStep(false);
                if (parameters.get_bVWarnings()) fprintf(stderr,"WARNING: bGravStep set in combination with older criteria, now using ONLY bGravStep!\n");
            }
            else if (!param.bAccelStep && !parameters.get_bGravStep() && !parameters.get_bDensityStep()) {
                parameters.set_bGravStep(true);
                if (parameters.get_bVWarnings()) fprintf(stderr,"WARNING: none of bAccelStep, bDensityStep, or bGravStep set, now using bGravStep!\n");
            }
            /*
            ** Set the needed memory model based on the chosen timestepping method.
            */
            if (parameters.get_bGravStep()) {
                parameters.set_bMemNodeAcceleration(1);
                if (parameters.get_iTimeStepCrit()) {
                    parameters.set_bMemNodeVelocity(1);
                }
            }
            else {
                parameters.set_bMemNodeAcceleration(0);
                parameters.set_bMemNodeVelocity(0);
            }
        }
    }

    if (parameters.get_bFindGroups() && !parameters.has_dTau()) {
        fprintf(stderr, "ERROR: you must specify dTau when FOF is to be run\n");
        return 0;
    }

    if (csm->val.classData.bClass) {
        const char *aLinear[MAX_CSM_SPECIES];
        const char *aPower[MAX_CSM_SPECIES];
        char *achLinSpecies = strdup(param.achLinSpecies);
        char *achPkSpecies = strdup(param.achPkSpecies);
        int nLinear = parseSpeciesNames(aLinear,achLinSpecies);
        int nPower = parseSpeciesNames(aPower,achPkSpecies);
        if (!prmSpecified(prm,"dOmega0")) csm->val.dOmega0 = 0.0;
        auto achClassFilename = parameters.get_achClassFilename();
        csmClassRead(csm, achClassFilename.data(), parameters.get_dBoxSize(), parameters.get_h(), nLinear, aLinear, nPower, aPower);
        free(achLinSpecies);
        free(achPkSpecies);
        csmClassGslInitialize(csm);
    }
    if (strlen(param.achLinSpecies) && param.nGridLin == 0) {
        fprintf(stderr, "ERROR: you must specify nGridLin when running with linear species\n");
        abort();
    }
    return 1;
}
