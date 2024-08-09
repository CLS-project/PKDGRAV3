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
#include <fstream>
#include "master.h"
#include "units.h"
#include "fmt/format.h"
#include <fmt/ranges.h>
#include "SPH/SPHOptions.h"
using namespace fmt::literals; // Gives us ""_a and ""_format literals

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
        if ( parameters.get_bWriteIC() ) {
            Write(BuildIoName(0).c_str(),dTime,0);
        }
    }

    /* Read in a binary file */
    else if ( parameters.has_achInFile() ) {
        dTime = Read(parameters.get_achInFile()); /* May change nSteps/dDelta */
    }
    else {
        print_error("No input file specified\n");
    }
    if (parameters.get_bAddDelete()) CountSpecies();
    return dTime;
}

/******************************************************************************\
*   Calculate dDelta and nSteps from: iStartStep, nSteps, dRedTo, dDelta
\******************************************************************************/

bool MSR::getDeltaSteps(double dTime,int iStartStep) {
    dDelta_list.clear();    // Step with the dDelta
    iStep_list.clear();     // to this step number

    auto nSteps = parameters.get_nSteps();
    auto dDelta = parameters.get_dDelta();
    auto dRedTo = parameters.get_dRedTo();

    // nSteps gives the number of steps for each interval. Convert this to the
    // ending step number for each interval.
    if (nSteps.size()) {
        decltype(nSteps)::value_type sum = 0;
        for (auto nStep : nSteps) iStep_list.push_back(sum += nStep);
    }

    // If dDelta is used, then we just need to copy the values from the parameters
    if (dDelta.size()) {
        if (dDelta.size() != iStep_list.size()) {
            print_error("ERROR: dDelta must have the same number of elements as nSteps\n");
            print_error("dDelta: {}\n", dDelta);
            print_error("nSteps: {}\n", nSteps);
            return false;
        }
        dDelta_list = dDelta;
    }

    // Non-cosmolicical simulations. nSteps and dDelta should have been set above.
    if (!csm->val.bComove) {
        if (dRedTo.size()) {
            print_warning("WARNING: dRedTo is meaningless for non-cosmological simulations, ignoring.\n");
        }
        if (nSteps.size() == 0 || dDelta.size() == 0) {
            print_error("ERROR: nSteps and dDelta must be specified for non-cosmological simulations.\n");
            return false;
        }
    }
    // Cosmological simulations
    else {
        if (dRedTo.size() > 0) {
            if (dDelta.size() > 0) {
                print_error("ERROR: dRedTo and dDelta can not be given at the same time for cosmological simulations.\n");
                return false;
            }
            if (dRedTo.size() != iStep_list.size()) {
                print_error("ERROR: dRedTo must have the same number of elements as nSteps\n");
                print_error("dRedTo: {}\n", dRedTo);
                print_error("nSteps: {}\n", nSteps);
                return false;
            }
            // Catchup. dTime corresponds to the current step (iStartStep)
            while (iStep_list.size() && iStep_list.front() <= iStartStep) {
                dRedTo.erase(dRedTo.begin());
                iStep_list.erase(iStep_list.begin());
            }

            auto tFrom = dTime;
            auto sFrom = iStartStep;
            for (auto i = 0; i<dRedTo.size(); ++i) {
                auto z = dRedTo[i];
                auto sTo = iStep_list[i];
                auto tTo = csmExp2Time(csm,1.0/(z + 1.0));
                dDelta_list.push_back((tTo-tFrom)/(sTo - sFrom));
                tFrom = tTo;
                sFrom = sTo;
            }
        }
    }
    assert(dDelta_list.size() == iStep_list.size());
    while (iStep_list.size()>0 && iStartStep >= iStep_list.front()) {
        dDelta_list.erase(dDelta_list.begin());
        iStep_list.erase(iStep_list.begin());
    }
    parameters.set_dynamic("delta",dDelta_list.front());
    print_notice (" Step From     Step To  Delta-T\n");
    for (auto i=0, j=iStartStep; i<dDelta_list.size(); j=iStep_list[i],++i) {
        print_notice("{:10d}  {:10d}  {:.8g}\n",j,iStep_list[i],dDelta_list[i]);

    }
    return true;
}

/******************************************************************************\
*   Simulation Mode: normal operation mode for pkdgrav
\******************************************************************************/
void MSR::Simulate(double dTime) {
    const auto iStartStep = parameters.get_iStartStep();
    if (getDeltaSteps(dTime,iStartStep)) Simulate(dTime,iStartStep);
}
void MSR::Simulate(double dTime,int iStartStep,bool bRestart) {
    std::ofstream log;
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
        log.open(filename,std::ios::app);
        if (!log.is_open()) {
            print_error("Failed to open log file: {}\n",filename);
            abort();
        }
        // fprintf(fpLog,"# ");
        // for (auto i=0;i<argc;++i) fprintf(fpLog,"%s ",argv[i]);
        // fprintf(fpLog,"\n");
        // msrLogParams(msr,fpLog);
    }

    TimerHeader();

    if (parameters.get_bLightCone() && Comove()) {
        auto dBoxSize = parameters.get_dBoxSize();
        print("One, Two, Three replica depth is z={:.10g}, {:.10g}, {:.10g}\n",
              1.0/csmComoveLookbackTime2Exp(csm,1.0 / dLightSpeedSim(1*dBoxSize)) - 1.0,
              1.0/csmComoveLookbackTime2Exp(csm,1.0 / dLightSpeedSim(2*dBoxSize)) - 1.0,
              1.0/csmComoveLookbackTime2Exp(csm,1.0 / dLightSpeedSim(3*dBoxSize)) - 1.0 );
    }

    if (MeshlessHydro()) {
        ChemCompInit();
    }

#ifdef COOLING
    float redshift;
    if (csm->val.bComove)
        redshift = 1./csmTime2Exp(csm,dTime) - 1.;
    else
        redshift = 0.0;

    CoolingInit(redshift);
    CoolingUpdate(redshift);
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
    if (parameters.get_nGridPk()>0) OutputPk(iStartStep,dTime);

    int bKickOpen, bKickClose=0;
    uint8_t uRungMax;
    int iSec = time(0);

    double dDelta = SwitchDelta(dTime,iStartStep);
    if (parameters.get_bNewKDK()) {
        LightConeOpen(iStartStep + 1);
        bKickOpen = 1;
    }
    else bKickOpen = 0;
    if (DoGravity()) {
        /* Compute the grids of the linear species before doing gravity */
        auto nGridLin = parameters.get_nGridLin();
        if (csm->val.classData.bClass && parameters.get_achLinSpecies().length() && nGridLin > 0) {
            GridCreateFFT(nGridLin);
            SetLinGrid(dTime,dDelta,nGridLin,bKickClose,bKickOpen);
            if (parameters.get_bDoLinPkOutput())
                OutputLinPk( iStartStep, dTime);
            LinearKick(dTime,dDelta,bKickClose,bKickOpen);
            GridDeleteFFT();
        }
    }

    const bool bDoStartOutput = parameters.get_bWriteIC() && !parameters.has_nGrid() && !NewSPH();
    const bool bDoStartFof =  parameters.get_bFindGroups() && ( bDoStartOutput || parameters.get_bBHPlaceSeed() );
    if (bDoStartFof) {
        NewFof(parameters.get_dTau(),parameters.get_nMinMembers());
    }

    if (NewSPH()) {
        // Calculate Density
        SelAll(-1,1);
        SPHOptions SPHoptions = initializeSPHOptions(parameters,csm,dTime,dDelta);
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
        SPHoptions.doGravity = parameters.get_bDoGravity();
        SPHoptions.doDensity = 0;
        SPHoptions.doSPHForces = 1;
        SPHoptions.dofBallFactor = 0;
        TreeUpdateFlagBounds(bEwald,ROOT,0,SPHoptions);
        uRungMax = Gravity(0,MAX_RUNG,ROOT,0,dTime,dDelta,iStartStep,dTheta,0,bKickOpen,
                           bEwald,bGravStep,nPartRhoLoc,iTimeStepCrit,SPHoptions);
        MemStatus();
    }
    else if (DoGravity()) {
        SPHOptions SPHoptions = initializeSPHOptions(parameters,csm,dTime,dDelta);
        SPHoptions.doGravity = parameters.get_bDoGravity();
        uRungMax = Gravity(0,MAX_RUNG,ROOT,0,dTime,dDelta,iStartStep,dTheta,0,bKickOpen,
                           bEwald,bGravStep,nPartRhoLoc,iTimeStepCrit,SPHoptions);
        MemStatus();
    }
    if (DoGravity() && bGravStep) {
        assert(parameters.get_bNewKDK() == false);    /* for now! */
        BuildTree(bEwald);
        SPHOptions SPHoptions = initializeSPHOptions(parameters,csm,dTime,dDelta);
        SPHoptions.doGravity = parameters.get_bDoGravity();
        Gravity(0,MAX_RUNG,ROOT,0,dTime,dDelta,iStartStep,dTheta,0,0,
                bEwald,bGravStep,nPartRhoLoc,iTimeStepCrit,SPHoptions);
        MemStatus();
    }
    if (MeshlessHydro()) {
        InitSph(dTime, dDelta,bRestart);
    }
#ifdef BLACKHOLES
    uRungMax = GetMinDt();
#ifndef DEBUG_BH_ONLY
    BlackholeInit(uRungMax);
#endif
#endif

    if (bDoStartFof) {
        GroupStats();
    }

    double E=0,T=0,U=0,Eth=0,L[3]= {0,0,0},F[3]= {0,0,0},W=0;
    CalcEandL(MSR_INIT_E,dTime,&E,&T,&U,&Eth,L,F,&W);
    iSec = time(0) - iSec;
    if (LogInterval()) {
        fmt::print(log,"{time:.10e} {z:.10e} {E:.10e} {T:.10e} {U:.10e} {Eth:.10e} {L0:.10e} {L1:.10e} {L2:.10e} {F0:.10e} {F1:.10e} {F2:.10e} {W:.10e} {elapsed:d}\n",
                   "time"_a=dTime,
                   "z"_a=1.0/csmTime2Exp(csm,dTime)-1.0,
                   "E"_a=E,"T"_a=T,"U"_a=U,"Eth"_a=Eth,
                   "L0"_a=L[0],"L1"_a=L[1],"L2"_a=L[2],
                   "F0"_a=F[0],"F1"_a=F[1],"F2"_a=F[2],
                   "W"_a=W,"elapsed"_a=iSec);
        log.rdbuf()->pubsync();
    }

    if (bDoStartOutput) {
        Output(iStartStep,dTime,dDelta,0);
    }

    // Make sure that the tree is usable before the start of the simulation
    if (bDoStartFof || bDoStartOutput) {
        DomainDecomp();
        BuildTree(bEwald);
    }

    bKickOpen = 0;
    int iStop=0, bDoCheckpoint=0, bDoOutput=0;
    auto nSteps = iStep_list.back();
    for (auto iStep=iStartStep+1; iStep<=nSteps&&!iStop; ++iStep) {
        dTheta = set_dynamic(iStep,dTime);
        dDelta = SwitchDelta(dTime,iStep-1);
        lPrior = time(0);
        TimerRestart();
        if (parameters.get_bNewKDK()) {
            double diStep = (double)(iStep-1);
            double ddTime = dTime;
            if (bKickOpen) {
                BuildTree(0);
                LightConeOpen(iStep);  /* open the lightcone */
                if (NewSPH()) {
                    // Calculate Density
                    SelAll(-1,1);
                    SPHOptions SPHoptions = initializeSPHOptions(parameters,csm,dTime,dDelta);
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
                    SPHoptions.doGravity = parameters.get_bDoGravity();
                    SPHoptions.doDensity = 0;
                    SPHoptions.doSPHForces = 1;
                    SPHoptions.dofBallFactor = 0;
                    TreeUpdateFlagBounds(bEwald,ROOT,0,SPHoptions);
                    uRungMax = Gravity(0,MAX_RUNG,ROOT,0,ddTime,dDelta,diStep,dTheta,0,1,
                                       bEwald,bGravStep,nPartRhoLoc,iTimeStepCrit,SPHoptions);
                }
                else {
                    SPHOptions SPHoptions = initializeSPHOptions(parameters,csm,dTime,dDelta);
                    SPHoptions.doGravity = parameters.get_bDoGravity();
                    uRungMax = Gravity(0,MAX_RUNG,ROOT,0,ddTime,dDelta,diStep,dTheta,0,1,
                                       bEwald,bGravStep,nPartRhoLoc,iTimeStepCrit,SPHoptions);
                }
                /* Set the grids of the linear species */
                auto nGridLin = parameters.get_nGridLin();
                if (csm->val.classData.bClass && parameters.get_achLinSpecies().length() && nGridLin > 0) {
                    GridCreateFFT(nGridLin);
                    SetLinGrid(dTime,dDelta,nGridLin,bKickClose,bKickOpen);
                    if (parameters.get_bDoLinPkOutput())
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
            fmt::print(log,"{time:.10e} {z:.10e} {E:.10e} {T:.10e} {U:.10e} {Eth:.10e} {L0:.10e} {L1:.10e} {L2:.10e} {F0:.10e} {F1:.10e} {F2:.10e} {W:.10e} {elapsed:d}\n",
                       "time"_a=dTime,
                       "z"_a=1.0/csmTime2Exp(csm,dTime)-1.0,
                       "E"_a=E,"T"_a=T,"U"_a=U,"Eth"_a=Eth,
                       "L0"_a=L[0],"L1"_a=L[1],"L2"_a=L[2],
                       "F0"_a=F[0],"F1"_a=F[1],"F2"_a=F[2],
                       "W"_a=W,"elapsed"_a=lSec);
            log.rdbuf()->pubsync();
        }
        if (!parameters.get_bNewKDK()) {
            CheckForOutput(iStep,nSteps,dTime,&bDoCheckpoint,&bDoOutput);
        }
        iStop = (bDoCheckpoint&2) || (bDoOutput&2);
        if (bDoCheckpoint) {
            Checkpoint(iStep,nSteps,dTime,dDelta);
            bDoCheckpoint = 0;
        }
        if (bDoOutput) {
            if (NewSPH() && parameters.get_bCentrifugal()) {
                ResetCOM();
            }
            Output(iStep,dTime,dDelta,0);
            bDoOutput = 0;
            DomainDecomp();
            BuildTree(bEwald);
        }
        TimerDump(iStep);
    }
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
static void validate_path(const char *name,const std::string_view &path) {
    if (path.length()) {
        using namespace fmt::literals;
        auto r1 = fmt::format(path,"name"_a="name","step"_a=1,"type"_a=".dat");
        if (r1 == fmt::format(path,"name"_a="name","step"_a=1,"type"_a=".XXX"))
            throw fmt::format_error(std::string(name) + ": " + std::string(path));
        if (r1 == fmt::format(path,"name"_a="name","step"_a=9,"type"_a=".dat"))
            throw fmt::format_error(std::string(name) + ": " + std::string(path));
        if (r1 == fmt::format(path,"name"_a="XXXX","step"_a=1,"type"_a=".dat"))
            throw fmt::format_error(std::string(name) + ": " + std::string(path));
    }
}

/*
** This routine validates the given parameters and makes any adjustments.
*/
int MSR::ValidateParameters() {
    int success = 1;
    try {
        validate_path("achOutPath",       parameters.get_achOutPath());
        validate_path("achCheckpointPath",parameters.get_achCheckpointPath());
        validate_path("achIoPath",        parameters.get_achIoPath());
    }
    catch (fmt::format_error &e) {
        print_error("ERROR: {}\n{}",e.what(),
                    "       When specified must contain {name}, {step} and {type} and no other fields\n"
                    "       Default: {name}.{step:05d}{type}\n"
                    "       Example: /path/to/output/{step:05d}/{name}.{step:05d}{type}\n");
        return 0;
    }

    //**************************************************************************
    // Verify the hydro model
    //**************************************************************************
    auto model = parameters.get_hydro_model();
    switch (model) {
    case HYDRO_MODEL::NONE:
    case HYDRO_MODEL::SPH:
    case HYDRO_MODEL::MFM:
    case HYDRO_MODEL::MFV:
        break;
    default:
        print_error("ERROR: Unknown hydro_model\n");
        return 0;
    }

    if (model == HYDRO_MODEL::MFM || model == HYDRO_MODEL::MFV) {
        if (parameters.get_bNewKDK()) {
            parameters.set_bNewKDK(false);
            print_warning("WARNING: Meshless hydrodynamics does not support bNewKDK. Setting bNewKDK to false\n");
        }
        if (parameters.get_bDualTree()) {
            parameters.set_bDualTree(false);
            print_warning("WARNING: Meshless hydrodynamics does not support bDualTree. "
                          "Setting bDualTree to false\n");
        }
        if (parameters.get_bMemIntegerPosition()) {
            parameters.set_bMemIntegerPosition(false);
            print_warning("WARNING: Meshless hydrodynamics does not support bMemIntegerPosition. "
                          "Setting bMemIntegerPosition to false\n");
        }
        if (!parameters.get_bMemUnordered()) {
            parameters.set_bMemUnordered(true);
            print_warning("WARNING: Meshless hydrodynamics requires bMemUnordered. "
                          "Setting bMemUnordered to true\n");
        }
        if (model == HYDRO_MODEL::MFV && !parameters.get_bMemMass()) {
            parameters.set_bMemMass(true);
            print_warning("WARNING: Meshless Finite Volume scheme requires bMemMass. "
                          "Setting bMemMass to true\n");
        }
#ifdef COOLING
        if (parameters.get_bComove()) {
            auto fH_reion_z = parameters.get_fH_reion_z();
            auto dRedTo = parameters.get_dRedTo();
            if (!std::any_of(dRedTo.begin(), dRedTo.end(), [fH_reion_z](double z) {
            return std::fabs(z - fH_reion_z) < 1e-10;
            })) {
                print_error("ERROR: Meshless hydrodynamics with cooling requires an element of dRedTo "
                            "to be set to the redshift of Hydrogen reionization, as set in parameter "
                            "fH_reion_z. Include the value of fH_reion_z in dRedTo, and choose a "
                            "value for the corresponding element of nSteps\n");
                return 0;
            }
        }
#endif
    }

#ifdef BLACKHOLES
    if (!ValidateBlackholeParam()) return 0;
#endif

#ifdef STAR_FORMATION
    if (!ValidateStarFormationParam()) return 0;
#endif

#ifdef STELLAR_EVOLUTION
    if (parameters.get_bChemEnrich() && !parameters.get_bMemMass()) {
        parameters.set_bMemMass(true);
        print_warning("WARNING: Chemical enrichment requires bMemMass. "
                      "Setting bMemMass to true\n");
    }
#endif

    if (parameters.get_bGasInterfaceCorrection() && parameters.get_bGasOnTheFlyPrediction()) {
        print_warning("Warning: On-the-fly prediction is not compatible with interface correction, disabled\n");
        parameters.set_bGasOnTheFlyPrediction(false);
    }

    if ((parameters.get_dVelocityDamper() > 0.0) && ((parameters.get_dVelocityDamperEnd() > 0.0) || (parameters.get_dVelocityDamperEndTime() > 0.0))) {
        if (parameters.get_dVelocityDamperEnd() <= 0.0) {
            print_error("ERROR: dVelocityDamper and dVelocityDamperEndTime specified, but not dVelocityDamperEnd.\n");
            return 0;
        }
        if (parameters.get_dVelocityDamperEndTime() <= 0.0) {
            print_error("ERROR: dVelocityDamper and dVelocityDamperEnd specified, but not dVelocityDamperEndTime.\n");
            return 0;
        }
        if (!parameters.get_bGasConsistentPrediction()) {
            print_error("ERROR: dVelocityDamper, dVelocityDamperEnd and dVelocityDamperEndTime specified, but not bGasConsistentPrediction.\n");
            return 0;
        }
    }

#ifndef NN_FLAG_IN_PARTICLE
    if (NewSPH() && parameters.get_bGasInterfaceCorrection() && parameters.get_dFastGasFraction() > 0.0f) {
        print_error("ERROR: Interface correction and FastGas is active, but the NN flag is not compiled in. Set NN_FLAG_IN_PARTICLE to ON in CMakeLists.txt and recompile.\n");
        return 0;
    }
#endif

    /*
    ** CUDA likes a larger group size
    */
    if ( (mdl->isCudaActive() || mdl->isMetalActive()) && !parameters.has_nGroup() && parameters.get_nGroup()<256)
        parameters.set_nGroup(256);

#ifndef USE_HDF5
    if (parameters.get_bHDF5()) {
        print_warning("WARNING: HDF5 output was requested but it is not supported: using Tipsy format\n");
        parameters.set_bHDF5(false);
    }
#endif

#ifdef MDL_FFTW
    auto nGridPk = parameters.get_nGridPk();
    if ( nGridPk ) {
        parameters.set_nBinsPk(std::min(parameters.has_nBinsPk() ? parameters.get_nBinsPk() : PST_MAX_K_BINS, nGridPk/2));
    }
    auto iPkOrder = parameters.get_iPkOrder();
    if (iPkOrder<ASSIGNMENT_ORDER::NGP || iPkOrder>ASSIGNMENT_ORDER::PCS) {
        print_error("ERROR: iPkOrder must be 0 (NGP), 1 (CIC), 2 (TSC) or 3 (PCS)\n");
        return 0;
    }
    if ( parameters.get_nGrid() ) {
        if (parameters.has_achInFile()) {
            print_error("ERROR: do not specify an input file when generating IC\n");
            return 0;
        }
        if ( parameters.get_iSeed() == 0 ) {
            //print_error("ERROR: Random seed for IC not specified"\n);
            parameters.set(parameters.str_iSeed,time(NULL));
        }
        if ( !parameters.has_dBoxSize() || parameters.get_dBoxSize() <= 0 ) {
            print_error("ERROR: Box size for IC not specified\n");
            return 0;
        }
        if ( !parameters.get_bClass() ) {
            if ( ( !parameters.has_dSigma8() || parameters.get_dSigma8() <= 0 ) &&
                    ( !parameters.has_dNormalization() || parameters.get_dNormalization() <= 0 ) ) {
                print_error("ERROR: Either dSigma8 or dNormalization should be specified for generating IC\n");
                return 0;
            }
            if ( !parameters.has_dSpectral() || parameters.get_dSpectral() <= 0 ) {
                print_error("ERROR: dSpectral for IC not specified\n");
                return 0;
            }
        }
        if ( parameters.get_bICgas() ) {
            if ( !parameters.has_dOmegab() || parameters.get_dOmegab() <= 0 ) {
                print_error("ERROR: Can not generate IC with gas if dOmegab is not specified\n");
                return 0;
            }
            if ( !DoGas() ) {
                print_error("ERROR: Can not generate gas if a hydrodynamic solver is not selected\n");
                return 0;
            }
        }
    }
    if ( parameters.get_bComove() && !parameters.get_bClass() ) {
        if ( !parameters.has_h() ) {
            print_warning("WARNING: Running with bComove without specifying a Hubble parameter, h\n");
        }
    }
#endif
    auto dFracNoDomainRootFind = parameters.get_dFracNoDomainRootFind();
    auto dFracNoDomainDecomp = parameters.get_dFracNoDomainDecomp();
    auto dFracNoDomainDimChoice = parameters.get_dFracNoDomainDimChoice();
    auto dFracDualTree = parameters.get_dFracDualTree();
    if (!parameters.has_dFracNoDomainRootFind() && dFracNoDomainRootFind  > dFracNoDomainDimChoice)parameters.set_dFracNoDomainRootFind(dFracNoDomainRootFind=dFracNoDomainDimChoice);
    if (!parameters.has_dFracNoDomainDecomp()   && dFracNoDomainDecomp    > dFracNoDomainRootFind) parameters.set_dFracNoDomainDecomp(dFracNoDomainDecomp=dFracNoDomainRootFind);
    if (!parameters.has_dFracDualTree()         && dFracDualTree          > dFracNoDomainDecomp)   parameters.set_dFracDualTree(dFracDualTree=dFracNoDomainDecomp);
    if (!parameters.has_dFracNoDomainDecomp()   && dFracNoDomainDecomp    < dFracDualTree)         parameters.set_dFracNoDomainDecomp(dFracNoDomainDecomp=dFracDualTree);
    if (!parameters.has_dFracNoDomainRootFind() && dFracNoDomainRootFind  < dFracNoDomainDecomp)   parameters.set_dFracNoDomainRootFind(dFracNoDomainRootFind=dFracNoDomainDecomp);
    if (!parameters.has_dFracNoDomainDimChoice()&& dFracNoDomainDimChoice < dFracNoDomainRootFind) parameters.set_dFracNoDomainDimChoice(dFracNoDomainDimChoice=dFracNoDomainRootFind);
    if ( dFracDualTree > dFracNoDomainDecomp
            || dFracNoDomainDecomp > dFracNoDomainRootFind
            || dFracNoDomainRootFind > dFracNoDomainDimChoice
            || dFracNoDomainDecomp<0.0 || dFracNoDomainDimChoice > 1.0 ) {
        print_error("ERROR: check that 0 <= dFracNoDomainDecomp <= dFracNoDomainRootFind <= dFracNoDomainDimChoice <= 1\n");
        return 0;
    }

    /* Make sure that the old behaviour is obeyed. */
    if ( NoSteps() ) {
        if ( !parameters.has_bDoAccOutput() ) parameters.set_bDoAccOutput(true);
        if ( !parameters.has_bDoPotOutput() ) parameters.set_bDoPotOutput(true);
    }

    /*
     * Softening
     */
    auto bPhysicalSoft = parameters.get_bPhysicalSoft();
    const auto dMaxPhysicalSoft = parameters.get_dMaxPhysicalSoft();
    if (bPhysicalSoft && !parameters.get_bComove()) {
        print_warning("WARNING: bPhysicalSoft reset to 0 for non-comoving (bComove == 0)\n");
        parameters.set_bPhysicalSoft(bPhysicalSoft = false);
    }
    if (bPhysicalSoft && dMaxPhysicalSoft>0) {
        print_error("ERROR: Setting both bPhysicalSoft and dMaxPhysicalSoft "
                    "is not allowed.\n Did you mean to limit the physical softening"
                    "with bPhysicalSoft and dSoftMax? or just limit the comoving "
                    "softening with dMaxPhysicalSoft?\n");
        return 0;
    }
    if ( dMaxPhysicalSoft>0 && parameters.get_dSoft()==0.0 && !parameters.get_bSoftMaxMul()) {
        print_error("ERROR: Trying to limit individual softenings setting a "
                    "maximum physical softening rather than a factor...\nThis is "
                    "not supported.\n Did you mean to use dSoft for a global softening? "
                    "or bSoftMaxMul for setting the limit as a factor?\n");
        return 0;
    }
    if ( bPhysicalSoft && parameters.get_dSoftMax()==0.0) {
        print_error("ERROR: If setting bPhysicalSoft, dSoftMax should be "
                    "provided to avoid divergences in the early universe.\n");
        return 0;
    }
    /*
    ** Periodic boundary conditions can be disabled along any of the
    ** x,y,z axes by specifying a period of zero for the given axis.
    ** Internally, the period is set to infinity (Cf. pkdBucketWalk()
    ** and pkdDrift(); also the INTERSECT() macro in smooth.h).
    */
    auto period = parameters.get_dPeriod();
    if (period[0] == 0) period[0] = FLOAT_MAXVAL;
    if (period[1] == 0) period[1] = FLOAT_MAXVAL;
    if (period[2] == 0) period[2] = FLOAT_MAXVAL;
    parameters.set_dPeriod(period);
    /*
    ** At the moment, integer positions are only really safe in periodic boxes!Wr
    */
    if (parameters.get_bMemIntegerPosition() && (!parameters.get_bPeriodic()||blitz::any(period!=1.0))) {
        print_warning("WARNING: Integer coordinates are enabled but the the box is not periodic\n"
                      "       and/or the box size is not 1. Set bPeriodic=1 and dPeriod=1.\n");
    }

    /*
    ** Check timestepping and gravity combinations.
    */
    const auto iMaxRung = parameters.get_iMaxRung();
    assert(iMaxRung <= IRUNGMAX);
    if (parameters.get_bEpsAccStep()) parameters.set_bAccelStep(true);
    if (parameters.get_bDoGravity()) {
        /* Potential is optional, but the default for gravity */
        if (!parameters.has_bMemPotential()) parameters.set_bMemPotential(1);
        if (iMaxRung < 1) {
            parameters.set_iMaxRung(0);
            if (parameters.get_bVWarnings()) print_warning("WARNING: iMaxRung set to 0, SINGLE STEPPING run!\n");
            /*
            ** For single stepping we don't need fancy timestepping variables.
            */
            parameters.set_bMemNodeAcceleration(0);
            parameters.set_bMemNodeVelocity(0);
        }
        else {
            if ((parameters.get_bAccelStep() || parameters.get_bDensityStep()) && parameters.get_bGravStep()) {
                /*
                ** We cannot combine these 2 types of timestepping criteria, we need to choose one
                ** or the other basic timestep criterion, in this case we choose only bGravStep.
                */
                parameters.set_bAccelStep(false);
                parameters.set_bEpsAccStep(false);
                parameters.set_bDensityStep(false);
                if (parameters.get_bVWarnings()) print_warning("WARNING: bGravStep set in combination with older criteria, now using ONLY bGravStep!\n");
            }
            else if (!parameters.get_bAccelStep() && !parameters.get_bGravStep() && !parameters.get_bDensityStep()) {
                parameters.set_bGravStep(true);
                if (parameters.get_bVWarnings()) print_warning("WARNING: none of bAccelStep, bDensityStep, or bGravStep set, now using bGravStep!\n");
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
        print_error("ERROR: you must specify dTau when FOF is to be run\n");
        return 0;
    }

    csmInitialize(&csm);
    csm->val.bComove = parameters.get_bComove();
    csm->val.dHubble0 = parameters.get_dHubble0();
    csm->val.h = parameters.get_h();
    csm->val.dOmega0 = parameters.get_dOmega0();
    csm->val.dLambda = parameters.get_dLambda();
    csm->val.dOmegaDE = parameters.get_dOmegaDE();
    csm->val.w0 = parameters.get_w0();
    csm->val.wa = parameters.get_wa();
    csm->val.dOmegaRad = parameters.get_dOmegaRad();
    csm->val.dOmegab = parameters.get_dOmegab();
    csm->val.dSigma8 = parameters.get_dSigma8();
    csm->val.dNormalization = parameters.get_dNormalization();
    csm->val.dSpectral = parameters.get_dSpectral();
    csm->val.dRunning = parameters.get_dRunning();
    csm->val.dPivot = parameters.get_dPivot();
    csm->val.classData.bClass = parameters.get_bClass();

    if (csm->val.classData.bClass) {
        const char *aLinear[MAX_CSM_SPECIES];
        const char *aPower[MAX_CSM_SPECIES];
        char *achLinSpecies = strdup(parameters.get_achLinSpecies().data());
        char *achPkSpecies = strdup(parameters.get_achPkSpecies().data());
        int nLinear = parseSpeciesNames(aLinear,achLinSpecies);
        int nPower = parseSpeciesNames(aPower,achPkSpecies);
        if (!parameters.has_dOmega0()) csm->val.dOmega0 = 0.0;
        auto achClassFilename = parameters.get_achClassFilename();
        csmClassRead(csm, achClassFilename.data(), parameters.get_dBoxSize(), parameters.get_h(), nLinear, aLinear, nPower, aPower);
        free(achLinSpecies);
        free(achPkSpecies);
        csmClassGslInitialize(csm);
    }
    if (parameters.get_achLinSpecies().length() && parameters.get_nGridLin() == 0) {
        print_error("ERROR: you must specify nGridLin when running with linear species\n");
        abort();
    }
    return success;
}
