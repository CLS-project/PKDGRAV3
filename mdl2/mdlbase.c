#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "mdlbase.h"
//#include <WinSock2.h> /* gethostname */
#ifdef __linux__
#include <sys/resource.h>
#else
#include <time.h>
#endif
#include <assert.h>
#if !defined(HAVE_CONFIG_H) || defined(HAVE_MALLOC_H)
#include <malloc.h>
#endif
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#include <stdarg.h>

#define MDL_DEFAULT_SERVICES	120

static void _srvNull(void *p1, void *vin, int nIn, void *vout, int *pnOut) {
    return;
    }

void mdlBaseInitialize(mdlBASE *base) {
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

    assert(sid > 0);
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

    getrusage(0, &ru);
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
    getrusage(0, &ru);
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

    getrusage(0, &ru);
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
#endif

