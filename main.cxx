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

#include "io/service.h"
#include "io/restore.h"

#ifdef HAVE_ROCKSTAR
    #include "analysis/rshalocount.h"
    #include "analysis/rshaloloadids.h"
#endif
#include "analysis/rsloadids.h"
#include "analysis/rssaveids.h"
#include "analysis/rsreorder.h"
#include "analysis/rsextract.h"

#include "initlightcone.h"

#include "domains/calcbound.h"
#include "domains/combinebound.h"
#include "domains/distribtoptree.h"
#include "domains/distribroot.h"
#include "domains/dumptrees.h"
#include "domains/enforceperiodic.h"
#include "domains/freestore.h"
#include "domains/olddd.h"
#include "domains/getordsplits.h"

#include "gravity/setsoft.h"
#include "gravity/activerung.h"
#include "gravity/countrungs.h"
#include "gravity/zeronewrung.h"

time_t timeGlobalSignalTime = 0;
int bGlobalOutput = 0;

#ifndef _MSC_VER
static inline void USR1_handler(int signo) {
    signal(SIGUSR1,USR1_handler);
    timeGlobalSignalTime = time(0);
}

static inline void USR2_handler(int signo) {
    signal(SIGUSR2,USR2_handler);
    bGlobalOutput = 1;
}
#endif

/*
** This function is called at the very start by every thread.
** It returns the "worker context"; in this case the PST.
*/
void *worker_init(MDL vmdl) {
    auto mdl = reinterpret_cast<mdl::mdlClass *>(vmdl);
    PST pst;
    LCL *plcl = new LCL;
    plcl->pkd = NULL;
    pstInitialize(&pst,mdl,plcl);
    pstAddServices(pst,vmdl);
    mdl->AddService(std::make_unique<ServiceSetAdd>(pst));
    mdl->AddService(std::make_unique<ServiceSwapAll>(pst));
    mdl->AddService(std::make_unique<ServiceHostname>(pst));
    mdl->AddService(std::make_unique<ServiceInitCosmology>(pst));
    mdl->AddService(std::make_unique<ServiceInitLightcone>(pst));
    mdl->AddService(std::make_unique<ServiceCalcRoot>(pst));
    mdl->AddService(std::make_unique<ServiceCountSelected>(pst));
    mdl->AddService(std::make_unique<ServiceSelBlackholes>(pst));
    mdl->AddService(std::make_unique<ServiceSelSpecies>(pst));
    mdl->AddService(std::make_unique<ServiceSelActives>(pst));
    mdl->AddService(std::make_unique<ServiceSelGroup>(pst));
    mdl->AddService(std::make_unique<ServiceSelMass>(pst));
    mdl->AddService(std::make_unique<ServiceSelPhaseDensity>(pst));
    mdl->AddService(std::make_unique<ServiceSelById>(pst));
    mdl->AddService(std::make_unique<ServiceSelBox>(pst));
    mdl->AddService(std::make_unique<ServiceSelSphere>(pst));
    mdl->AddService(std::make_unique<ServiceSelCylinder>(pst));
    mdl->AddService(std::make_unique<ServiceFileSizes>(pst));
    mdl->AddService(std::make_unique<ServiceRestore>(pst));
    mdl->AddService(std::make_unique<ServiceCalcBound>(pst));
    mdl->AddService(std::make_unique<ServiceCombineBound>(pst));
    mdl->AddService(std::make_unique<ServiceDistribTopTree>(pst));
    mdl->AddService(std::make_unique<ServiceDistribRoot>(pst));
    mdl->AddService(std::make_unique<ServiceDumpTrees>(pst));
    mdl->AddService(std::make_unique<ServiceEnforcePeriodic>(pst));
    mdl->AddService(std::make_unique<ServiceFreeStore>(pst));
    mdl->AddService(std::make_unique<ServiceSetSoft>(pst));
    mdl->AddService(std::make_unique<ServiceActiveRung>(pst));
    mdl->AddService(std::make_unique<ServiceCountRungs>(pst));
    mdl->AddService(std::make_unique<ServiceZeroNewRung>(pst));
    mdl->AddService(std::make_unique<ServiceGetOrdSplits>(pst));
#ifdef HAVE_ROCKSTAR
    mdl->AddService(std::make_unique<ServiceRsHaloCount>(pst));
    mdl->AddService(std::make_unique<ServiceRsHaloLoadIds>(pst));
#endif
    mdl->AddService(std::make_unique<ServiceRsLoadIds>(pst));
    mdl->AddService(std::make_unique<ServiceRsSaveIds>(pst));
    mdl->AddService(std::make_unique<ServiceRsReorderIds>(pst));
    mdl->AddService(std::make_unique<ServiceRsExtract>(pst));
    mdl->AddService(std::make_unique<OldDD::ServiceDomainDecomp>(pst));
    mdl->AddService(std::make_unique<OldDD::ServiceDomainOrder>(pst));
    mdl->AddService(std::make_unique<OldDD::ServiceLocalOrder>(pst));
    mdl->AddService(std::make_unique<OldDD::ServiceColRejects>(pst));
    mdl->AddService(std::make_unique<OldDD::ServiceSwapRejects>(pst));
    mdl->AddService(std::make_unique<OldDD::ServiceColOrdRejects>(pst));
    mdl->AddService(std::make_unique<OldDD::ServiceWeight>(pst));
    mdl->AddService(std::make_unique<OldDD::ServiceWeightWrap>(pst));
    mdl->AddService(std::make_unique<OldDD::ServiceOrdWeight>(pst));
    return pst;
}

/*
** This function is called at the very end for every thread.
** It needs to destroy the worker context (PST).
*/
void worker_done(MDL mdl, void *ctx) {
    PST pst = reinterpret_cast<PST>(ctx);
    LCL *plcl = pst->plcl;
    pstFinish(pst);
    delete plcl;
}

/*
** This is invoked for the "master" process after the worker has been setup.
*/
int master(MDL mdl,void *vpst) {
    auto pst = reinterpret_cast<PST>(vpst);
    int argc = mdlGetArgc(mdl);
    char **argv = mdlGetArgv(mdl);

    /* a USR1 signal indicates that the queue wants us to exit */
#ifndef _MSC_VER
    timeGlobalSignalTime = 0;
    signal(SIGUSR1,NULL);
    signal(SIGUSR1,USR1_handler);

    /* a USR2 signal indicates that we should write an output when convenient */
    bGlobalOutput = 0;
    signal(SIGUSR2,NULL);
    signal(SIGUSR2,USR2_handler);
#endif

    MSR msr(mdl,pst);
    auto rc = msr.Python(argc,argv);
    if (rc < 0) {
        printf("%s using Python %d.%d.%d\n", PACKAGE_STRING, PY_MAJOR_VERSION, PY_MINOR_VERSION, PY_MICRO_VERSION );
        if (!msr.ValidateParameters())
            return 2;

        /* Establish safety lock. */
        if (!msr.GetLock()) {
            return 1;
        }

        msr.Hostname(); // List all host names

        auto dTime = msr.LoadOrGenerateIC();
        if (dTime != -HUGE_VAL) msr.Simulate(dTime);
        rc = 0;
    }
    return rc;
}

int main(int argc,char **argv) {
#ifdef ENABLE_FE
    feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
#endif
#ifndef CCC
    /* no stdout buffering */
    setbuf(stdout,(char *) NULL);
#endif

    return mdlLaunch(argc,argv,master,worker_init,worker_done);
}
