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
using namespace mdl;

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
#include <string.h>
#include <numeric>

/*****************************************************************************\
* LegacyService
*
* This is the legacy MDL service support where a function pointer and parameter
* are used. RunService() just calls fcnService() with p1 as the first parameter.
\*****************************************************************************/
class LegacyService : public BasicService {
protected:
    fcnService_t *fcnService;
    void *p1;
public:
    explicit LegacyService(fcnService_t *fcnService,void *p1=nullptr, int service_id=0,
                           int nInBytes=0, int nOutBytes=0, const char *service_name="")
        : BasicService(service_id, nInBytes, nOutBytes, service_name),
          fcnService(fcnService), p1(p1) {}
    virtual ~LegacyService() = default;
protected:
    virtual int operator()(int nIn, void *pIn, void *pOut) override;
};

int LegacyService::operator()(int nIn, void *pIn, void *pOut) {
    assert(nIn <= getMaxBytesIn());
    assert(fcnService != NULL);
    auto nOut = (*fcnService)(p1, pIn, nIn, pOut, getMaxBytesOut());
    assert(nOut <= getMaxBytesOut());
    return nOut;
}

/*****************************************************************************\
* mdlBASE
\*****************************************************************************/

#define MDL_DEFAULT_SERVICES    120

static int _srvNull(void *p1, void *vin, int nIn, void *vout, int nOut) {
    return 0;
}

mdlBASE::mdlBASE(int argc,char **argv) {
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
    layout.nThreads = 1;
    layout.idSelf = 0;
    layout.nProcs = 1;
    layout.iProc = 0;
    layout.nCores = 1;
    layout.iCore = 0;
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
    services[0] = std::make_unique<LegacyService>(_srvNull);
}

mdlBASE::~mdlBASE() {
#ifdef _MSC_VER
    WSACleanup();
#endif
}

/* O(1): Given a process id, return the first global thread id */
int32_t mdlBASE::ProcToThread(int32_t iProc) const {
    assert(iProc >= 0 && iProc <= layout.nProcs);
    return iProcToThread[iProc];
}

/* O(l2(nProc)): Given a global thread id, return the process to which it belongs */
int32_t mdlBASE::ThreadToProc(int32_t iThread) const {
    int l=0, u=layout.nProcs;
    assert(iThread >= 0 && iThread <= layout.nThreads);
    assert(layout.nThreads == iProcToThread[layout.nProcs]);
    while (l <= u) {
        int m = (u + l) / 2;
        if (iThread < iProcToThread[m]) u = m - 1;
        else l = m+1;
    }
    return l-1;
}

void mdlBASE::AddService(int sid, void *p1,
                         fcnService_t *fcnService,
                         int nInBytes, int nOutBytes, const char *name) {
    if (nInBytes  > nMaxInBytes)  nMaxInBytes  = nInBytes;
    if (nOutBytes > nMaxOutBytes) nMaxOutBytes = nOutBytes;
    if (sid >= services.size()) services.resize(sid+9);
    assert(services[sid]==nullptr);
    services[sid] = std::make_unique<LegacyService>(fcnService,p1,sid,nInBytes,nOutBytes,name);
}

void mdlBASE::AddService(std::unique_ptr<BasicService> &&service) {
    auto sid = service->getServiceID();
    if (service->getMaxBytesIn()  > nMaxInBytes)  nMaxInBytes  = service->getMaxBytesIn();
    if (service->getMaxBytesOut() > nMaxOutBytes) nMaxOutBytes = service->getMaxBytesOut();
    if (sid >= services.size()) services.resize(sid+9);
    assert(services[sid]==nullptr);
    services[sid] = std::move(service);
}

int mdlBASE::RunService(int sid,int nIn, void *pIn, void *pOut) {
    assert(sid < services.size());
    assert(services[sid] != nullptr);
    return (*services[sid])(nIn,pIn,pOut);
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

double mdlCpuTimer(void *mdl) {
#ifdef __linux__
    struct rusage ru;

    getrusage(RUSAGE_SELF, &ru);
    return ((double)ru.ru_utime.tv_sec + 1e-6*(double)ru.ru_utime.tv_usec);
#elif defined(_MSC_VER)
    FILETIME createTime;
    FILETIME exitTime;
    FILETIME kernelTime;
    FILETIME userTime;
    if (GetProcessTimes(GetCurrentProcess(),
                        &createTime, &exitTime, &kernelTime, &userTime) != -1) {
        SYSTEMTIME userSystemTime;
        if (FileTimeToSystemTime(&userTime, &userSystemTime) != -1)
            return (double)userSystemTime.wHour * 3600.0 +
                   (double)userSystemTime.wMinute * 60.0 +
                   (double)userSystemTime.wSecond +
                   (double)userSystemTime.wMilliseconds / 1000.0;
    }
    return 0.0;
#else
    return (((double)clock()) / CLOCKS_PER_SEC);
#endif
}

void mdlDiag(void *mdl, const char *psz) {
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
void mdlZeroTimer(void *mdl, mdlTimer *t) {
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

void mdlPrintTimer(void *mdl, const char *message, mdlTimer *t0) {
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
    std::fill(dTimer,dTimer+TIME_COUNT,0.0);
    //dWaiting = dComputing = dSynchronizing = 0.0;
#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
    nTicks = getticks();
#endif
}

void mdlBASE::TimeAddComputing() {
#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
    ticks nTicks = getticks();
    dTimer[TIME_COMPUTING] += elapsed(nTicks, this->nTicks);
    this->nTicks = nTicks;
#endif
}

void mdlBASE::TimeAddSynchronizing() {
#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
    ticks nTicks = getticks();
    dTimer[TIME_SYNCHRONIZING] += elapsed(nTicks, this->nTicks);
    this->nTicks = nTicks;
#endif
}

void mdlBASE::TimeAddWaiting() {
#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
    ticks nTicks = getticks();
    dTimer[TIME_WAITING] += elapsed(nTicks, this->nTicks);
    this->nTicks = nTicks;
#endif
}

double mdlBASE::TimeFraction() const {
    double dTotal = std::accumulate(dTimer,dTimer+TIME_COUNT,0.0);
    if (dTotal <= 0.0) return 0.0;
    return 100.0 / dTotal;
}

double mdlBASE::TimeComputing() const {
    return dTimer[TIME_COMPUTING] * TimeFraction();
}

double mdlBASE::TimeSynchronizing() const {
    return dTimer[TIME_SYNCHRONIZING] * TimeFraction();
}

double mdlBASE::TimeWaiting() const {
    return dTimer[TIME_WAITING] * TimeFraction();
}
#endif

