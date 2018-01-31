#ifndef MDLBASE_H
#define MDLBASE_H
#ifdef HAVE_CONFIG_H
#include "config.h"
#else
#include "mdl_config.h"
#endif
#ifdef INSTRUMENT
#include "cycle.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "opa_queue.h"

typedef struct {
    OPA_Queue_element_hdr_t hdr;
    uint32_t iServiceID;
    uint32_t iCoreFrom;
    } MDLserviceElement;

/*
** A MDL Key must be large enough to hold the largest unique particle key.
** It must be (several times) larger than the total number of particles.
** An "unsigned long" is normally 32 bits on a 32 bit machine and
** 64 bits on a 64 bit machine.  We use uint64_t to be sure.
*/
#ifndef MDLKEY
#define MDLKEY uint64_t
#endif
typedef MDLKEY mdlkey_t;
static const mdlkey_t MDL_INVALID_KEY = (mdlkey_t)(-1);

/*
** The purpose of this routine is to safely cast a size_t to an int.
** If this were done inline, the type of "v" would not be checked.
**
** We cast for the following reasons:
** - We have masked a mdlkey_t (to an local id or processor), so we
**   know it will fit in an integer.
** - We are taking a part of a memory block to send or receive
**   through MPI.  MPI takes an "int" parameter.
** - We have calculated a memory size with pointer subtraction and
**   know that it must be smaller than 4GB.
**
** The compiler will cast automatically, but warnings can be generated unless
** the cast is done explicitly.  Putting it inline is NOT type safe as in:
**   char *p;
**   int i;
**   i = (int)p;
** "works", while:
**   i = size_t_to_int(p);
** would fail.
*/
static inline int size_t_to_int(size_t v) {
    return (int)v;
    }

static inline int mdlkey_t_to_int(mdlkey_t v) {
    return (int)v;
    }

/*
* Compile time mdl debugging options
*
* mdl asserts: define MDLASSERT
* Probably should always be on unless you want no mdlDiag output at all
*
* NB: defining NDEBUG turns off all asserts so MDLASSERT will not assert
* however it will output uding mdlDiag and the code continues.
*/
#define MDLASSERT

typedef struct serviceRec {
    int nInBytes;
    int nOutBytes;
    void *p1;
    void(*fcnService)(void *, void *, int, void *, int *);
} SERVICE;

#define MAX_NODE_NAME_LENGTH      256
typedef struct {
    int32_t nThreads; /* Global number of threads (total) */
    int32_t idSelf;   /* Global index of this thread */
    int32_t nProcs;   /* Number of global processes (e.g., MPI ranks) */
    int32_t iProc;    /* Index of this process (MPI rank) */
    int16_t nCores;   /* Number of threads in this process */
    int16_t iCore;    /* Local core id */
    int bDiag;        /* When true, debug output is enabled */
    int argc;

    FILE *fpDiag;
    char **argv;

    /* Services information */
    int nMaxServices;
    int nMaxInBytes;
    int nMaxOutBytes;
    SERVICE *psrv;

    /* Maps a give process (Proc) to the first global thread ID */
    int *iProcToThread; /* [0,nProcs] (note inclusive extra element) */

    char nodeName[MAX_NODE_NAME_LENGTH];

#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
    ticks nTicks;
    double dWaiting;
    double dComputing;
    double dSynchronizing;
#endif


    } mdlBASE;

void mdlBaseInitialize(mdlBASE *base,int argc,char **argv);
void mdlBaseFinish(mdlBASE *base);
void mdlBaseAddService(mdlBASE *base, int sid, void *p1,
    void(*fcnService)(void *, void *, int, void *, int *),
    int nInBytes, int nOutBytes);

#define mdlThreads(mdl) ((mdl)->base.nThreads)
#define mdlSelf(mdl) ((mdl)->base.idSelf)
#define mdlCore(mdl) ((mdl)->base.iCore)
#define mdlCores(mdl) ((mdl)->base.nCores)
#define mdlProc(mdl) ((mdl)->base.iProc)
#define mdlProcs(mdl) ((mdl)->base.nProcs)
const char *mdlName(void *mdl);

int mdlBaseProcToThread(mdlBASE *base, int iProc);
int mdlBaseThreadToProc(mdlBASE *base, int iThread);
#define mdlProcToThread(mdl,iProc) mdlBaseProcToThread(&(mdl)->base,iProc)
#define mdlThreadToProc(mdl,iThread) mdlBaseThreadToProc(&(mdl)->base,iThread)

void mdlDiag(void *mdl, char *psz);
void mdlprintf(void *mdl, const char *format, ...);
#ifdef MDLASSERT
#ifndef __STRING
#define __STRING( arg )   (("arg"))
#endif
#define mdlassert(mdl,expr) \
    { \
    if (!(expr)) { \
    mdlprintf( mdl, "%s:%d Assertion `%s' failed.\n", __FILE__, __LINE__, __STRING(expr) ); \
    assert( expr ); \
        } \
    }
#else
#define mdlassert(mdl,expr)  assert(expr)
#endif

double mdlCpuTimer(void * mdl);


#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
void mdlTimeReset(void *mdl);
void mdlTimeAddComputing(void *mdl);
void mdlTimeAddSynchronizing(void *mdl);
void mdlTimeAddWaiting(void *mdl);
double mdlTimeComputing(void *mdl);
double mdlTimeSynchronizing(void *mdl);
double mdlTimeWaiting(void *mdl);
#else
#define mdlTimeReset(mdl)
#define mdlTimeAddComputing(mdl)
#define mdlTimeAddSynchronizing(mdl)
#define mdlTimeAddWaiting(mdl)
#define mdlTimeComputing(mdl) 0.0
#define mdlTimeSynchronizing(mdl) 0.0
#define mdlTimeWaiting(mdl) 0.0
#endif



/*
* Timer functions active: define MDLTIMER
* Makes mdl timer functions active
*/
#ifndef _CRAYMPP
#define MDLTIMER
#endif

typedef struct {
    double wallclock;
    double cpu;
    double system;
    } mdlTimer;

#ifdef MDLTIMER
void mdlZeroTimer(void * mdl, mdlTimer *);
void mdlGetTimer(void * mdl, mdlTimer *, mdlTimer *);
void mdlPrintTimer(void *mdl, char *message, mdlTimer *);
#else
#define mdlZeroTimer
#define mdlGetTimer
#define mdlPrintTimer
#endif

#endif
