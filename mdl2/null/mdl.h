#ifndef MDL_HINCLUDED
#define MDL_HINCLUDED
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#ifdef MDL_FFTW
#include <srfftw.h>
#endif
#ifdef INSTRUMENT
#include "cycle.h"
#endif

#if defined(__osf__) || defined(__sgi)
#define vsnprintf(a,b,c,d) vsprintf((a),(c),(d))
#endif

#ifdef __cplusplus
extern "C" {
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
void *mdlMalloc(MDL,size_t);
void mdlFree(MDL,void *);
void mdlSetCacheSize(MDL,int);
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

/*
** Collective operations
*/
typedef int MDL_Op;

#define MDL_BAND ((MDL_Op)0x40000001)
#define MDL_BOR ((MDL_Op)0x40000002)
#define MDL_BXOR ((MDL_Op)0x40000003)
#define MDL_LAND ((MDL_Op)0x40000004)
#define MDL_LOR ((MDL_Op)0x40000005)
#define MDL_LXOR ((MDL_Op)0x40000006)
#define MDL_MAX ((MDL_Op)0x40000007)
#define MDL_MAXLOC ((MDL_Op)0x40000008)
#define MDL_MIN ((MDL_Op)0x40000009)
#define MDL_MINLOC ((MDL_Op)0x4000000a)
#define MDL_PROD ((MDL_Op)0x4000000b)
#define MDL_REPLACE ((MDL_Op)0x4000000c)
#define MDL_SUM ((MDL_Op)0x4000000d)

typedef int MDL_Datatype;
#define MDL_FLOAT ((MDL_Datatype)0x50000001)
#define MDL_DOUBLE ((MDL_Datatype)0x50000002)
#define MDL_BYTE ((MDL_Datatype)0x50000003)
#define MDL_INT ((MDL_Datatype)0x50000004)
#define MDL_LONG_LONG ((MDL_Datatype)0x50000005)

int mdlReduce ( MDL mdl, void *sendbuf, void *recvbuf, int count,
		MDL_Datatype datatype, MDL_Op op, int root );
int mdlAllreduce( MDL mdl, void *sendbuf, void *recvbuf, int count,
		  MDL_Datatype datatype, MDL_Op op );
int mdlAlltoall( MDL mdl, void *sendbuf, int scount, MDL_Datatype stype,
    void *recvbuf, int rcount, MDL_Datatype rtype);
int mdlAlltoallv( MDL mdl, void *sendbuf, int *sendcnts, int *sdispls, MDL_Datatype sendtype,
    void *recvbuf, int *recvcnts, int *rdispls, MDL_Datatype recvtype);
int mdlAllGather( MDL mdl, void *sendbuf, int scount, MDL_Datatype stype,
    void *recvbuf, int rcount, MDL_Datatype recvtype);
int mdlReduceScatter( MDL mdl, void* sendbuf, void* recvbuf, int *recvcounts,
    MDL_Datatype datatype, MDL_Op op);
int mdlTypeContiguous(MDL mdl,int count, MDL_Datatype old_type, MDL_Datatype *newtype);
int mdlTypeCommit(MDL mdl, MDL_Datatype *datatype );
int mdlTypeFree (MDL mdl, MDL_Datatype *datatype );

/*
** Grid Operations
*/

typedef struct mdlGridContext {
    uint32_t n1, n2, n3; /* Real dimensions */
    uint32_t a1;         /* Actual size of dimension 1 */
    uint32_t nlocal;     /* Number of local elements */
    } * MDLGRID;


/*
** Allocate a MDLGRID context.  This has no actual data, but only describes
** the grid geometry.  The global geometry is set.
*/
void mdlGridInitialize(MDL mdl,MDLGRID *pgrid,int n1,int n2,int n3,int a1);
/*
** Free all memory associated with a MDLGRID context.
*/
void mdlGridFinish(MDL mdl, MDLGRID grid);
/*
** Sets the local geometry (i.e., what is on this processor) of this grid.
*/
void mdlGridSetLocal(MDL mdl,MDLGRID grid,int s, int n, int nlocal);
/*
** Share the local geometry between processors.
*/
void mdlGridShare(MDL mdl,MDLGRID grid);
/*
** Allocate the local elements.  The size of a single element is
** given and the local GRID information is consulted to determine
** how many to allocate.
*/
void *mdlGridMalloc(MDL mdl,MDLGRID grid,int nEntrySize);
void mdlGridFree( MDL mdl, MDLGRID grid, void *p );
/*
** This gives the processor on which the given slab can be found.
*/
static inline int mdlGridId(MDLGRID grid, uint32_t x, uint32_t y, uint32_t z) {
    return 0; /* Always processor 0 */
    }
/*
** This returns the index into the array on the appropriate processor.
*/
static inline int mdlGridIdx(MDLGRID grid, uint32_t x, uint32_t y, uint32_t z) {
    assert(x<grid->a1&&y<grid->n2&&z<grid->n3);
    return x + grid->a1*(y + grid->n2*z); /* Local index */
    }

/*
** FFT Operations
*/
#ifdef MDL_FFTW
typedef struct mdlFFTContext {
    MDLGRID rgrid;
    MDLGRID kgrid;
    rfftwnd_plan fplan;
    rfftwnd_plan iplan;
    } * MDLFFT;

size_t mdlFFTInitialize(MDL mdl,MDLFFT *fft,
			int nx,int ny,int nz,int bMeasure);
void mdlFFTFinish( MDL mdl, MDLFFT fft );
fftw_real *mdlFFTMAlloc( MDL mdl, MDLFFT fft );
void mdlFFTFree( MDL mdl, MDLFFT fft, void *p );
void mdlFFT( MDL mdl, MDLFFT fft, fftw_real *data, int bInverse );

/* Grid accessors: r-space */
#define mdlFFTrId(fft,x,y,z) mdlGridId((fft)->rgrid,x,y,z)
#define mdlFFTrIdx(fft,x,y,z) mdlGridIdx((fft)->rgrid,x,y,z)

/* Grid accessors: k-space (note permuted indices) */
#define mdlFFTkId(fft,x,y,z) mdlGridId((fft)->kgrid,x,z,y)
#define mdlFFTkIdx(fft,x,y,z) mdlGridIdx((fft)->kgrid,x,z,y)
#endif

void mdlSetWorkQueueSize(MDL,int);
void mdlAddWork(MDL mdl, int (*doWork)(void *ctx), void *ctx);

#ifdef __cplusplus
    }
#endif

#endif

