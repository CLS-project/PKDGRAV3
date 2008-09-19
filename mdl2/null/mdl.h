#ifndef MDL_HINCLUDED
#define MDL_HINCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdio.h>
#include <assert.h>
#ifdef INSTRUMENT
#include "cycle.h"
#endif

#if defined(__osf__) || defined(__sgi)
#define vsnprintf(a,b,c,d) vsprintf((a),(c),(d))
#endif

#define MAX_PROCESSOR_NAME      256

#define SRV_STOP		0

typedef struct cacheSpace {
    int iType;
    void *pData;
    int iDataSize;
    int nData;
    void *ctx;
    void (*init)(void *,void *);
    void (*combine)(void *,void *,void *);
    void * (*getElt)(void *pData,int i,int iDataSize);
    /*
     ** Statistics stuff.
     */
    int nAccess;
    int nAccHigh;
    long nMiss;
    long nColl;
    long nMin;
    } CACHE;


typedef struct serviceRec {
    int nInBytes;
    int nOutBytes;
    void *p1;
    void (*fcnService)(void *,void *,int,void *,int *);
    } SERVICE;


typedef struct mdlContext {
    int nThreads;
    int idSelf;
    int bDiag;
    int *atid;
    FILE *fpDiag;
    /*
     ** Services stuff!
     */
    int nMaxServices;
    int nMaxInBytes;
    int nMaxOutBytes;
    SERVICE *psrv;
    char *pszIn;
    char *pszOut;
    /*
     ** Caching stuff!
     */
    unsigned long uRand;
    int iMaxDataSize;
    int nMaxCacheIds;
    CACHE *cache;
    char nodeName[MAX_PROCESSOR_NAME];
#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
    ticks nTicks;
    double dComputing;
#endif
    } * MDL;


/*
 * MDL debug and Timer macros and prototypes
 */
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
/*
 * Debug functions active: define MDLDEBUG
 * Adds debugging mdldebug prints and mdldebugassert asserts
 */
#define MDLDEBUG
/*
 * Timer functions active: define MDLTIMER
 * Makes mdl timer functions active
 */
#ifndef _CRAYMPP
#define MDLTIMER
#endif


void mdlprintf( MDL mdl, const char *format, ... );

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

#ifdef MDLDEBUG
#define mdldebugassert(mdl,expr)   mdlassert(mdl,expr)
void mdldebug( MDL mdl, const char *format, ... );
#else
#define mdldebug
#define mdldebugassert
#endif

typedef struct {
    double wallclock;
    double cpu;
    double system;
    } mdlTimer;

#ifdef MDLTIMER
void mdlZeroTimer(MDL mdl,mdlTimer *);
void mdlGetTimer(MDL mdl,mdlTimer *,mdlTimer *);
void mdlPrintTimer(MDL mdl,char *message,mdlTimer *);
#else
#define mdlZeroTimer
#define mdlGetTimer
#define mdlPrintTimer
#endif

/*
 ** General Functions
 */
double mdlCpuTimer(MDL);
int mdlInitialize(MDL *pmdl,char **argv,void (*fcnChild)(MDL),void (*fcnIOChild)(MDL));
void mdlFinish(MDL);
int mdlThreads(MDL);
int mdlSelf(MDL);
int mdlOldSelf(MDL);
const char *mdlName(MDL);
int mdlSwap(MDL,int,size_t,void *,size_t,size_t *,size_t *);
void mdlDiag(MDL,char *);
void mdlAddService(MDL,int,void *,void (*)(void *,void *,int,void *,int *),
		   int,int);
void mdlReqService(MDL,int,int,void *,int);
void mdlGetReply(MDL,int,void *,int *);
void mdlHandler(MDL);
/*
 ** Caching functions.
 */
void *mdlMalloc(MDL,int);
void mdlFree(MDL,void *);
void mdlROcache(MDL mdl,int cid,
                void * (*getElt)(void *pData,int i,int iDataSize),
                void *pData,int iDataSize,int nData);
void mdlCOcache(MDL mdl,int cid,
                void * (*getElt)(void *pData,int i,int iDataSize),
                void *pData,int iDataSize,int nData,
                void *ctx,void (*init)(void *,void *),void (*combine)(void *,void *,void *));
void mdlFinishCache(MDL,int);
void mdlCacheCheck(MDL);
void mdlCacheBarrier(MDL,int);
void *mdlAquire(MDL,int,int,int);
void mdlPrefetch(MDL,int,int,int);
void mdlRelease(MDL,int,void *);
/*
 ** Cache statistics functions.
 */
double mdlNumAccess(MDL,int);
double mdlMissRatio(MDL,int);
double mdlCollRatio(MDL,int);
double mdlMinRatio(MDL,int);

#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
void mdlTimeReset(MDL);
double mdlTimeComputing(MDL);
double mdlTimeSynchronizing(MDL);
double mdlTimeWaiting(MDL);
#endif

#endif

