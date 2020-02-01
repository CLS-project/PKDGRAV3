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
void *worker_init(MDL mdl) {
    PST pst;
    LCL *plcl = new LCL;
    plcl->pkd = NULL;
    pstInitialize(&pst,mdl,plcl);
    pstAddServices(pst,mdl);
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
void master(MDL mdl,void *vpst) {
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
    if (!msr.Python(argc,argv)) {
	printf("%s using Python %d.%d.%d\n", PACKAGE_STRING, PY_MAJOR_VERSION, PY_MINOR_VERSION, PY_MICRO_VERSION );
	msr.ValidateParameters();

	/* Establish safety lock. */
	if (!msr.GetLock()) {
	    return;
	    }

	msr.Hostname(); // List all host names

	auto dTime = msr.LoadOrGenerateIC();
	if (dTime != -HUGE_VAL) msr.Simulate(dTime);
	}
    }

int main(int argc,char **argv) {
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

    mdlLaunch(argc,argv,master,worker_init,worker_done);

    return 0;
    }
