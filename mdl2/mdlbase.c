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

#include "mdlbase.h"
//#include <WinSock2.h> /* gethostname */
#ifdef __linux__
#include <sys/resource.h>
#else
#include <time.h>
#endif
#include <assert.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include <stdarg.h>

#define MDL_DEFAULT_SERVICES	120

static void _srvNull(void *p1, void *vin, int nIn, void *vout, int *pnOut) {
    return;
    }

void mdlBaseInitialize(mdlBASE *base,int argc,char **argv) {
    int i;

#ifdef _MSC_VER
    WSADATA wsaData;
    WSAStartup(0x202, &wsaData);
#endif
    if (gethostname(base->nodeName, sizeof(base->nodeName)))
        base->nodeName[0] = 0;
    else
        base->nodeName[sizeof(base->nodeName) - 1] = 0;

    base->bDiag = 0;
    base->fpDiag = NULL;

    base->argc = argc;
    base->argv = argv;

    /* Some sensible defaults */
    base->nThreads = 1;
    base->idSelf = 0;
    base->nProcs = 1;
    base->iProc = 0;
    base->nCores = 1;
    base->iCore = 0;
    base->iProcToThread = NULL;
    base->nTicks = 0;

    /*
    ** Set default "maximums" for structures. These are NOT hard
    ** maximums, as the structures will be realloc'd when these
    ** values are exceeded.
    */
    base->nMaxServices = MDL_DEFAULT_SERVICES;
    base->nMaxInBytes  = 0;
    base->nMaxOutBytes = 0;
    /*
    ** Now allocate the initial service slots.
    */
    base->psrv = malloc(base->nMaxServices*sizeof(SERVICE));
    assert(base->psrv != NULL);
    /*
    ** Provide a 'null' service for sid = 0, so that stopping the
    ** service handler is well defined!
    */
    base->psrv[0].p1 = NULL;
    base->psrv[0].nInBytes = 0;
    base->psrv[0].nOutBytes = 0;
    base->psrv[0].fcnService = _srvNull;
    /*
    ** Initialize the remaining new service slots.
    */
    for (i = 1; i<base->nMaxServices; ++i) {
        base->psrv[i].p1 = NULL;
        base->psrv[i].nInBytes = 0;
        base->psrv[i].nOutBytes = 0;
        base->psrv[i].fcnService = NULL;
        }
    }

/* O(1): Given a process id, return the first global thread id */
int mdlBaseProcToThread(mdlBASE *base, int iProc) {
    assert(iProc >= 0 && iProc <= base->nProcs);
    return base->iProcToThread[iProc];
    }

/* O(l2(nProc)): Given a global thread id, return the process to which it belongs */
int mdlBaseThreadToProc(mdlBASE *base, int iThread) {
    int l=0, u=base->nProcs;
    assert(iThread >= 0 && iThread <= base->nThreads);
    assert(base->nThreads == base->iProcToThread[base->nProcs]);
    while (l <= u) {
	int m = (u + l) / 2;
	if (iThread < base->iProcToThread[m]) u = m - 1;
	else l = m+1;
	}
    return l-1;
    }

void mdlBaseFinish(mdlBASE *base) {
    free(base->psrv);

#ifdef _MSC_VER
    WSACleanup();
#endif
    }

void mdlBaseAddService(mdlBASE *base, int sid, void *p1,
    void(*fcnService)(void *, void *, int, void *, int *),
    int nInBytes, int nOutBytes) {
    int i, nMaxServices;

    assert(sid >= 0); /* We can replace SRV_STOP to do something at the end */
    if (sid >= base->nMaxServices) {
        /*
        ** reallocate service buffer, adding space for 8 new services
        ** including the one just defined.
        */
        nMaxServices = sid + 9;
        base->psrv = realloc(base->psrv, nMaxServices*sizeof(SERVICE));
        assert(base->psrv != NULL);
        /*
        ** Initialize the new service slots.
        */
        for (i = base->nMaxServices; i<nMaxServices; ++i) {
            base->psrv[i].p1 = NULL;
            base->psrv[i].nInBytes = 0;
            base->psrv[i].nOutBytes = 0;
            base->psrv[i].fcnService = NULL;
        }
        base->nMaxServices = nMaxServices;
    }
    if (nInBytes  > base->nMaxInBytes)  base->nMaxInBytes  = nInBytes;
    if (nOutBytes > base->nMaxOutBytes) base->nMaxOutBytes = nOutBytes;
    base->psrv[sid].p1 = p1;
    base->psrv[sid].nInBytes = nInBytes;
    base->psrv[sid].nOutBytes = nOutBytes;
    base->psrv[sid].fcnService = fcnService;
}

const char *mdlName(void *mdl) {
    mdlBASE *base = mdl;
    return base->nodeName;
    }

double mdlCpuTimer(void * mdl) {
#ifdef __linux__
    struct rusage ru;

    getrusage(RUSAGE_SELF, &ru);
    return((double)ru.ru_utime.tv_sec + 1e-6*(double)ru.ru_utime.tv_usec);
#elif defined(_MSC_VER)
    FILETIME createTime;
    FILETIME exitTime;
    FILETIME kernelTime;
    FILETIME userTime;
    if (GetProcessTimes(GetCurrentProcess(),
        &createTime, &exitTime, &kernelTime, &userTime) != -1)
        {
        SYSTEMTIME userSystemTime;
        if (FileTimeToSystemTime(&userTime, &userSystemTime) != -1)
            return (double)userSystemTime.wHour * 3600.0 +
            (double)userSystemTime.wMinute * 60.0 +
            (double)userSystemTime.wSecond +
            (double)userSystemTime.wMilliseconds / 1000.0;
        }
    return 0.0;
#else
    return(((double)clock()) / CLOCKS_PER_SEC);
#endif
    }

void mdlDiag(void *mdl, char *psz) {
    mdlBASE *base = mdl;
    if (base->bDiag) {
        fputs(psz, base->fpDiag);
        fflush(base->fpDiag);
        }
    }


void mdlprintf(void *mdl, const char *format, ...) {
    mdlBASE *base = mdl;
    if (base->bDiag) {
        va_list args;
        va_start(args, format);
        vfprintf(base->fpDiag, format, args);
        fflush(base->fpDiag);
        va_end(args);
        }
    }

/*
* MDL debug and Timer functions
*/
#ifdef MDLTIMER
void mdlZeroTimer(void * mdl, mdlTimer *t) {
#ifdef _MSC_VER
    FILETIME ft;
    uint64_t clock;
    GetSystemTimeAsFileTime(&ft);
    clock = ft.dwHighDateTime;
    clock <<= 32;
    clock |= ft.dwLowDateTime;
    /* clock is in 100 nano-second units */
    t->wallclock = 1.0 * (clock / 10000000UL);
#else
    struct timezone tz;
    struct timeval tv;
    struct rusage ru;
    tz.tz_minuteswest = 0;
    tz.tz_dsttime = 0;
    gettimeofday(&tv, &tz);
    t->wallclock = tv.tv_sec + 1e-6*(double)tv.tv_usec;
    getrusage(RUSAGE_SELF, &ru);
    t->cpu = (double)ru.ru_utime.tv_sec + 1e-6*(double)ru.ru_utime.tv_usec;
    t->system = (double)ru.ru_stime.tv_sec + 1e-6*(double)ru.ru_stime.tv_usec;
#endif
    }

void mdlGetTimer(void *mdl, mdlTimer *t0, mdlTimer *t) {
#ifdef _MSC_VER
    FILETIME ft;
    uint64_t clock;
    GetSystemTimeAsFileTime(&ft);
    clock = ft.dwHighDateTime;
    clock <<= 32;
    clock |= ft.dwLowDateTime;
    /* clock is in 100 nano-second units */
    t->wallclock = clock / 10000000UL - t0->wallclock;
#else
    struct timezone tz;
    struct timeval tv;
    struct rusage ru;

    getrusage(RUSAGE_SELF, &ru);
    t->cpu = (double)ru.ru_utime.tv_sec + 1e-6*(double)ru.ru_utime.tv_usec - t0->cpu;
    t->system = (double)ru.ru_stime.tv_sec + 1e-6*(double)ru.ru_stime.tv_usec - t0->system;
    tz.tz_minuteswest = 0;
    tz.tz_dsttime = 0;
    gettimeofday(&tv, &tz);
    t->wallclock = tv.tv_sec + 1e-6*(double)tv.tv_usec - t0->wallclock;
#endif
    }

void mdlPrintTimer(void *mdl, char *message, mdlTimer *t0) {
    mdlBASE *base = mdl;
    mdlTimer lt;

    if (base->bDiag) {
        mdlGetTimer(mdl, t0, &lt);
        mdlprintf(mdl, "%s %f %f %f\n", message, lt.wallclock, lt.cpu, lt.system);
        }
    }


#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
void mdlTimeReset(void *mdl) {
    mdlBASE *base = mdl;
    base->dWaiting = base->dComputing = base->dSynchronizing = 0.0;
    base->nTicks = getticks();
    }

void mdlTimeAddComputing(void *mdl) {
    mdlBASE *base = mdl;
    ticks nTicks = getticks();
    base->dComputing += elapsed(nTicks, base->nTicks);
    base->nTicks = nTicks;
    }

void mdlTimeAddSynchronizing(void *mdl) {
    mdlBASE *base = mdl;
    ticks nTicks = getticks();
    base->dSynchronizing += elapsed(nTicks, base->nTicks);
    base->nTicks = nTicks;
    }

void mdlTimeAddWaiting(void *mdl) {
    mdlBASE *base = mdl;
    ticks nTicks = getticks();
    base->dWaiting += elapsed(nTicks, base->nTicks);
    base->nTicks = nTicks;
    }

static double TimeFraction(void * mdl) {
    mdlBASE *base = mdl;
    double dTotal = base->dComputing + base->dWaiting + base->dSynchronizing;
    if (dTotal <= 0.0) return 0.0;
    return 100.0 / dTotal;
    }

double mdlTimeComputing(void * mdl) {
    mdlBASE *base = mdl;
    return base->dComputing * TimeFraction(mdl);
    }

double mdlTimeSynchronizing(void * mdl) {
    mdlBASE *base = mdl;
    return base->dSynchronizing * TimeFraction(mdl);
    }

double mdlTimeWaiting(void * mdl) {
    mdlBASE *base = mdl;
    return base->dWaiting * TimeFraction(mdl);
    }
#endif





#endif

