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

static int _srvNull(void *p1, void *vin, int nIn, void *vout, int nOut) {
    return 0;
    }



int mdlBASE::SERVICE::operator()(int nIn, char *pszIn, char *pszOut) {
    assert(nIn <= nInBytes);
    assert(fcnService != NULL);
    int nOut = (*fcnService)(p1, pszIn, nIn, pszOut, nOutBytes);
    assert(nOut <= nOutBytes);
    return nOut;
    }


mdlBASE::mdlBASE(int argc,char **argv) {
    int i;

#ifdef _MSC_VER
    WSADATA wsaData;
    WSAStartup(0x202, &wsaData);
#endif
    if (gethostname(nodeName, sizeof(nodeName)))
        nodeName[0] = 0;
    else
        nodeName[sizeof(nodeName) - 1] = 0;

    bDiag = 0;
    fpDiag = NULL;

    this->argc = argc;
    this->argv = argv;

    /* Some sensible defaults */
    nThreads = 1;
    idSelf = 0;
    nProcs = 1;
    iProc = 0;
    nCores = 1;
    iCore = 0;
    nTicks = 0;

    /*
    ** Set default "maximums" for structures. These are NOT hard
    ** maximums, as the structures will be realloc'd when these
    ** values are exceeded.
    */
    nMaxInBytes  = 0;
    nMaxOutBytes = 0;
    /*
    ** Now allocate the initial service slots.
    */
    services.resize(MDL_DEFAULT_SERVICES);
    /*
    ** Provide a 'null' service for sid = 0, so that stopping the
    ** service handler is well defined!
    */
    services[0] = SERVICE(_srvNull);
    }

mdlBASE::~mdlBASE() {
#ifdef _MSC_VER
    WSACleanup();
#endif
    }

/* O(1): Given a process id, return the first global thread id */
int32_t mdlBASE::ProcToThread(int32_t iProc) const {
    assert(iProc >= 0 && iProc <= nProcs);
    return iProcToThread[iProc];
    }

/* O(l2(nProc)): Given a global thread id, return the process to which it belongs */
int32_t mdlBASE::ThreadToProc(int32_t iThread) const {
    int l=0, u=nProcs;
    assert(iThread >= 0 && iThread <= nThreads);
    assert(nThreads == iProcToThread[nProcs]);
    while (l <= u) {
	int m = (u + l) / 2;
	if (iThread < iProcToThread[m]) u = m - 1;
	else l = m+1;
	}
    return l-1;
    }

void mdlBASE::AddService(int sid, void *p1,
    fcnService_t *fcnService,
    int nInBytes, int nOutBytes) {
    if (nInBytes  > nMaxInBytes)  nMaxInBytes  = nInBytes;
    if (nOutBytes > nMaxOutBytes) nMaxOutBytes = nOutBytes;
    if (sid >= services.size()) services.resize(sid+9);
    services[sid] = SERVICE(fcnService,p1,nInBytes,nOutBytes);
}

void mdlBASE::yield() {
#ifdef _MSC_VER
    SwitchToThread();
#else
    sched_yield();
#endif
    }

const char *mdlName(void *mdl) {
    mdlBASE *base = reinterpret_cast<mdlBASE *>(mdl);
    return base->nodeName;
    }

#if 0
int mdlThreads(void *mdl) {  return reinterpret_cast<mdlBASE *>(mdl)->nThreads; }
int mdlSelf(void *mdl)    {  return reinterpret_cast<mdlBASE *>(mdl)->idSelf; }
int mdlCore(void *mdl)    {  return reinterpret_cast<mdlBASE *>(mdl)->iCore; }
int mdlCores(void *mdl)   {  return reinterpret_cast<mdlBASE *>(mdl)->nCores; }
int mdlProc(void *mdl)    {  return reinterpret_cast<mdlBASE *>(mdl)->iProc; }
int mdlProcs(void *mdl)   {  return reinterpret_cast<mdlBASE *>(mdl)->nProcs; }

int mdlGetArgc(void *mdl) {  return reinterpret_cast<mdlBASE *>(mdl)->argc; }
char **mdlGetArgv(void *mdl) {  return reinterpret_cast<mdlBASE *>(mdl)->argv; }
#endif

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
#if 0
    mdlBASE *base = reinterpret_cast<mdlBASE *>(mdl);
    if (base->bDiag) {
        fputs(psz, base->fpDiag);
        fflush(base->fpDiag);
        }
#endif
    }

void mdlBASE::mdl_vprintf(const char *format, va_list ap) {
    if (bDiag) {
        vfprintf(fpDiag, format, ap);
        fflush(fpDiag);
        }
    }
void mdlBASE::mdl_printf(const char *format, ...) {
    va_list args;
    va_start(args, format);
    mdl_vprintf(format,args);
    va_end(args);
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
#if 0
    mdlBASE *base = reinterpret_cast<mdlBASE *>(mdl);
    mdlTimer lt;

    if (base->bDiag) {
        mdlGetTimer(mdl, t0, &lt);
        base->mdl_printf("%s %f %f %f\n", message, lt.wallclock, lt.cpu, lt.system);
        }
#endif
    }


void mdlBASE::TimeReset() {
    dWaiting = dComputing = dSynchronizing = 0.0;
#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
    nTicks = getticks();
#endif
    }

void mdlBASE::TimeAddComputing() {
#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
    ticks nTicks = getticks();
    dComputing += elapsed(nTicks, this->nTicks);
    this->nTicks = nTicks;
#endif
    }

void mdlBASE::TimeAddSynchronizing() {
#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
    ticks nTicks = getticks();
    dSynchronizing += elapsed(nTicks, this->nTicks);
    this->nTicks = nTicks;
#endif
    }

void mdlBASE::TimeAddWaiting() {
#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
    ticks nTicks = getticks();
    dWaiting += elapsed(nTicks, this->nTicks);
    this->nTicks = nTicks;
#endif
    }

double mdlBASE::TimeFraction() const {
    double dTotal = dComputing + dWaiting + dSynchronizing;
    if (dTotal <= 0.0) return 0.0;
    return 100.0 / dTotal;
    }

double mdlBASE::TimeComputing() const {
    return dComputing * TimeFraction();
    }

double mdlBASE::TimeSynchronizing() const {
    return dSynchronizing * TimeFraction();
    }

double mdlBASE::TimeWaiting() const {
    return dWaiting * TimeFraction();
    }
#endif

