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

#include "pkd_config.h"
#include "starformation/init.h"

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
#include "master.h"
#include "io/outtype.h"
#include "smooth/smoothfcn.h"

#include "core/setadd.h"
#include "core/swapall.h"
#include "core/hostname.h"
#include "core/initcosmology.h"
#include "core/calcroot.h"
#include "core/select.h"
#include "core/fftgrid.h"
#include "core/countspecies.h"
#include "core/removedeleted.h"
#include "core/ignoresigbus.h"
#include "core/decorator_reorder.h"
#ifdef INSTRUMENT_MPI
    #include "core/cachestats.h"
#endif
#include "io/service.h"
#include "io/restore.h"

#include "analysis/analysis.h"
#include "analysis/orbits.h"
#include "analysis/measurepk.h"
#include "analysis/assignmass.h"
#include "analysis/bispectrum.h"
#ifdef HAVE_ROCKSTAR
    #include "analysis/rshalocount.h"
    #include "analysis/rshaloloadids.h"
#endif
#include "analysis/rsloadids.h"
#include "analysis/rssaveids.h"
#include "analysis/rsreorder.h"
#include "analysis/rsextract.h"
#include "group/fof.h"
#include "group/hop.h"

#include "analysis/lightcone.h"

#include "io/output.h"
#include "io/classes.h"

#include "domains/buildtree.h"
#include "domains/calcbound.h"
#include "domains/combinebound.h"
#include "domains/distribtoptree.h"
#include "domains/distribroot.h"
#include "domains/dumptrees.h"
#include "domains/enforceperiodic.h"
#include "domains/freestore.h"
#include "domains/olddd.h"
#include "domains/reorder.h"
#include "domains/getordsplits.h"

#include "gravity/setsoft.h"
#include "gravity/activerung.h"
#include "gravity/countrungs.h"
#include "gravity/pmforces.h"
#include "gravity/updaterung.h"
#include "gravity/zeronewrung.h"
#include "gravity/timestep.h"

#include "SPH/gasvalues.h"
#include "SPH/resetcom.h"

#include "meshless/hydro/wakeup.h"
#include "meshless/hydro/endstep.h"

#include "meshless/cooling/services.h"
#include "meshless/starformation/init.h"
#include "meshless/starformation/starformation.h"

time_t timeGlobalSignalTime = 0;
int bGlobalOutput = 0;

#ifndef _MSC_VER
static inline void USR1_handler(int signo) {
    signal(SIGUSR1, USR1_handler);
    timeGlobalSignalTime = time(0);
}

static inline void USR2_handler(int signo) {
    signal(SIGUSR2, USR2_handler);
    bGlobalOutput = 1;
}
#endif

/*
** This function is called at the very start by every thread.
** It returns the "worker context"; in this case the PST.
*/
void* worker_init(MDL vmdl) {
    auto mdl = reinterpret_cast<mdl::mdlClass*>(vmdl);
    LCL* plcl = new LCL;
    plcl->pkd = NULL;
    auto pst = new pstNode(mdl, plcl);
    pstAddServices(pst, mdl);
    mdl->AddService<ServiceStarFormInit>(pst);
    mdl->AddService<ServiceStarForm>(pst);
    mdl->AddService<ServiceCoolingInit>(pst);
    mdl->AddService<ServiceCoolingUpdate>(pst);
    mdl->AddService<ServiceCoolingUpdateZ>(pst);
    mdl->AddService<ServiceCoolingHydReionization>(pst);
    mdl->AddService<ServiceWakeParticles>(pst);
    mdl->AddService<ServiceEndTimestepIntegration>(pst);
    mdl->AddService<ServiceSetAdd>(pst);
    mdl->AddService<ServiceSwapAll>(pst);
    mdl->AddService<ServiceHostname>(pst);
    mdl->AddService<ServiceInitCosmology>(pst);
    mdl->AddService<ServiceInitLightcone>(pst);
    mdl->AddService<ServiceCalcRoot>(pst);
    mdl->AddService<ServiceCountSelected>(pst);
    mdl->AddService<ServiceSelBlackholes>(pst);
    mdl->AddService<ServiceSelSpecies>(pst);
    mdl->AddService<ServiceSelActives>(pst);
    mdl->AddService<ServiceSelGroup>(pst);
    mdl->AddService<ServiceSelMass>(pst);
    mdl->AddService<ServiceSelPhaseDensity>(pst);
    mdl->AddService<ServiceSelById>(pst);
    mdl->AddService<ServiceSelBox>(pst);
    mdl->AddService<ServiceSelSphere>(pst);
    mdl->AddService<ServiceSelCylinder>(pst);
    mdl->AddService<ServiceFftSizes>(pst);
    mdl->AddService<ServiceIgnoreSIGBUS>(pst);
    mdl->AddService<ServiceCountSpecies>(pst);
    mdl->AddService<ServiceTotalMass>(pst);
    mdl->AddService<ServiceResetCOM>(pst);
    mdl->AddService<ServiceRemoveDeleted>(pst);
    mdl->AddService<ServiceFileSizes>(pst);
    mdl->AddService<ServiceRestore>(pst);
    mdl->AddService<ServiceRestoreDecorator>(pst);
    mdl->AddService<ServiceBuildTree>(pst);
    mdl->AddService<ServiceTreeUpdateFlagBounds>(pst);
    mdl->AddService<ServiceCalcBound>(pst);
    mdl->AddService<ServiceCombineBound>(pst);
    mdl->AddService<ServiceDistribTopTree>(pst);
    mdl->AddService<ServiceDistribRoot>(pst);
    mdl->AddService<ServiceDumpTrees>(pst);
    mdl->AddService<ServiceTreeInitMarked>(pst);
    mdl->AddService<ServiceEnforcePeriodic>(pst);
    mdl->AddService<ServiceFreeStore>(pst);
    mdl->AddService<ServiceSetSoft>(pst);
    mdl->AddService<ServiceActiveRung>(pst);
    mdl->AddService<ServiceCountRungs>(pst);
    mdl->AddService<ServiceUpdateRung>(pst);
    mdl->AddService<ServiceGridCreateFFT>(pst);
    mdl->AddService<ServiceGridDeleteFFT>(pst);
    mdl->AddService<ServiceLinearKick>(pst);
    mdl->AddService<ServiceMeasureLinearPk>(pst);
    mdl->AddService<ServiceAssignMass>(pst);
    mdl->AddService<ServiceWindowCorrection>(pst);
    mdl->AddService<ServiceDensityContrast>(pst);
    mdl->AddService<ServiceInterlace>(pst);
    mdl->AddService<ServiceGridBinK>(pst);
    mdl->AddService<ServiceSetLinearGrid>(pst);
    mdl->AddService<ServiceAddLinearSignal>(pst);
    mdl->AddService<ServiceBispectrumSelect>(pst);
    mdl->AddService<ServiceBispectrumCalculate>(pst);
    mdl->AddService<ServiceNewFof>(pst);
    mdl->AddService<ServiceFofPhases>(pst);
    mdl->AddService<ServiceFofFinishUp>(pst);
    mdl->AddService<ServiceHopLink>(pst);
    mdl->AddService<ServiceHopJoin>(pst);
    mdl->AddService<ServiceHopFinishUp>(pst);
    mdl->AddService<ServiceHopTreeBuild>(pst);
    mdl->AddService<ServiceHopGravity>(pst);
    mdl->AddService<ServiceHopUnbind>(pst);
    mdl->AddService<ServiceProfile>(pst);
    mdl->AddService<ServiceCalcDistance>(pst);
    mdl->AddService<ServiceCalcCOM>(pst);
    mdl->AddService<ServiceZeroNewRung>(pst);
    mdl->AddService<ServiceGetOrdSplits>(pst);
    mdl->AddService<ServiceSaveGroupStats>(pst);
    mdl->AddService<ServiceCountDistance>(pst);
    mdl->AddService<ServiceGetOrbits>(pst);
    mdl->AddService<ServiceSwapClasses>(pst);
    mdl->AddService<ServiceGetClasses>(pst);
    mdl->AddService<ServiceSetClasses>(pst);
    mdl->AddService<ServiceMaxRung>(pst);
    mdl->AddService<ServiceSetGlobalRung>(pst);
#ifdef HAVE_ROCKSTAR
    mdl->AddService<ServiceRsHaloCount>(pst);
    mdl->AddService<ServiceRsHaloLoadIds>(pst);
#endif
    mdl->AddService<ServiceRsLoadIds>(pst);
    mdl->AddService<ServiceRsSaveIds>(pst);
    mdl->AddService<ServiceRsReorderIds>(pst);
    mdl->AddService<ServiceRsExtract>(pst);
    mdl->AddService<ServiceReorderDecorator>(pst);
    mdl->AddService<OldDD::ServiceDomainDecomp>(pst);
#ifdef NEW_REORDER
    mdl->AddService<NewDD::ServiceReorder>(pst);
#else
    mdl->AddService<OldDD::ServiceDomainOrder>(pst);
    mdl->AddService<OldDD::ServiceLocalOrder>(pst);
#endif
    mdl->AddService<OldDD::ServiceColRejects>(pst);
    mdl->AddService<OldDD::ServiceSwapRejects>(pst);
    mdl->AddService<OldDD::ServiceColOrdRejects>(pst);
    mdl->AddService<OldDD::ServiceWeight>(pst);
    mdl->AddService<OldDD::ServiceWeightWrap>(pst);
    mdl->AddService<OldDD::ServiceOrdWeight>(pst);
    mdl->AddService<OldDD::ServiceActiveOrder>(pst);
    mdl->AddService<ServiceUpdateGasValues>(pst);
    mdl->AddService<ServiceLightConeOpen>(pst);
    mdl->AddService<ServiceLightConeClose>(pst);
    mdl->AddService<ServiceLightConeVel>(pst);
#ifdef INSTRUMENT_MPI
    mdl->AddService<ServiceCacheStats>(pst);
#endif
    return pst;
}

/*
** This function is called at the very end for every thread.
** It needs to destroy the worker context (PST).
*/
void worker_done(MDL mdl, void* ctx) {
    PST pst = reinterpret_cast<PST>(ctx);
    LCL* plcl = pst->lcl();

    while (pst) {
        auto pstKill = pst;
        if (pst->AmCore() && pst->lcl()->pkd) delete pst->lcl()->pkd;
        pst = pst->lower();
        delete pstKill;
    }
    delete plcl;
}

/*
** This is invoked for the "master" process after the worker has been setup.
*/
int master(MDL vmdl, void* vpst) {
    auto mdl = reinterpret_cast<mdl::mdlClass*>(vmdl);
    auto pst = reinterpret_cast<PST>(vpst);
    int argc = mdl->GetArgc();
    char** argv = mdl->GetArgv();

    /* a USR1 signal indicates that the queue wants us to exit */
#ifndef _MSC_VER
    timeGlobalSignalTime = 0;
    signal(SIGUSR1, NULL);
    signal(SIGUSR1, USR1_handler);

    /* a USR2 signal indicates that we should write an output when convenient */
    bGlobalOutput = 0;
    signal(SIGUSR2, NULL);
    signal(SIGUSR2, USR2_handler);
#endif

    MSR msr(mdl, pst);
    auto rc = msr.Python(argc, argv);
    if (rc < 0) {
        printf("%s using Python %d.%d.%d\n", PACKAGE_STRING, PY_MAJOR_VERSION, PY_MINOR_VERSION, PY_MICRO_VERSION);
        if (!msr.ValidateParameters()) return 2;

        /* Establish safety lock. */
        if (!msr.GetLock()) {
            return 1;
        }

        if (msr.parameters.get_bIgnoreSIGBUS()) msr.IgnoreSIGBUS();
        msr.Hostname(); // List all host names

        auto dTime = msr.LoadOrGenerateIC();
        if (dTime != -HUGE_VAL) msr.Simulate(dTime);
        rc = 0;
    }
    return rc;
}

int main(int argc, char** argv) {
#ifdef ENABLE_FE
    feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
#endif
#ifndef CCC
    /* no stdout buffering */
    setvbuf(stdout, NULL, _IONBF, 0);
#endif

    return mdl::Launch(master, worker_init, worker_done, argc, argv);
}
