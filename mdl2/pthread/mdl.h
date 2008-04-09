#ifndef MDL_HINCLUDED
#define MDL_HINCLUDED
#include <stdio.h>
#include <pthread.h>
#include <assert.h>
#include <unistd.h>

#ifdef MDL_FFTW
#include <srfftw_threads.h>
#endif

#ifdef __osf__
#define vsnprintf(a,b,c,d) vsprintf((a),(c),(d))
#endif

#define MDL_CACHE_SIZE		2000000
#define MDL_CACHELINE_BITS	3
#define MDL_CACHELINE_ELTS	(1<<MDL_CACHELINE_BITS)
#define MDL_CACHE_MASK		(MDL_CACHELINE_ELTS-1)
#define MDL_INDEX_MASK		(~MDL_CACHE_MASK)

#define MDL_MBX_RING_SZ	8

#define MAX_PROCESSOR_NAME      256

#define SRV_STOP		0

#ifndef MDLKEY
#define MDLKEY unsigned long
#endif
typedef MDLKEY mdlkey_t;
static const mdlkey_t MDL_INVALID_KEY = (mdlkey_t)(-1);

typedef struct cacheTag {
    mdlkey_t iKey;
    int nLock;
    int nLast;
    int iLink;
    } CTAG;

typedef struct cacheHeader {
    int cid;
    int mid;
    int id;
    int iLine;
    int iSeq;
    } CAHEAD;

typedef struct mbxStruct {
    pthread_mutex_t mux;
    pthread_cond_t sigReq;
    pthread_cond_t sigRpl;
    pthread_cond_t sigRel;
    int bReq;
    int bRpl;
    int bRel;
    int sid;
    int nBytes;
    char *pszIn;
    char *pszOut;
    } MBX;

typedef struct swxStruct {
    pthread_mutex_t mux;
    pthread_cond_t sigRel;
    pthread_cond_t sigRec;
    pthread_cond_t sigSnd;
    int bRel;
    int bRec;
    size_t nInBytes;
    size_t nOutBufBytes;
    char *pszBuf;
    } SWX;

typedef struct cacheSpace {
    int iType;
    char *pData;
    int iDataSize;
    int nData;
    size_t pDataMax;
    int iLineSize;
    int nLines;
    int nTrans;
    mdlkey_t iTransMask;
    mdlkey_t iIdMask;
    int iKeyShift;
    int iInvKeyShift;
    int *pTrans;
    CTAG *pTag;
    char *pLine;
    int nCheckOut;
    void (*init)(void *);
    void (*combine)(void *,void *);
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
    pthread_t *pt;
    struct mdlContext **pmdl;
    /*
     ** Services stuff!
     */
    int nMaxServices;
    int nMaxInBytes;
    int nMaxOutBytes;
    SERVICE *psrv;
    MBX mbxOwn;
    /*
     ** Swapping Box.
     */
    SWX swxOwn;
    /*
     ** Caching stuff!
     */
    MBX mbxCache[MDL_MBX_RING_SZ];
    pthread_mutex_t muxRing;
    int iRingHd;
    int iRingTl;
    int iRecSeq[128];
    int iSndSeq[128];
    unsigned long uRand;
    int iCaBufSize;
    int nMaxCacheIds;
    int iMaxDataSize;
    CACHE *cache;
    char nodeName[MAX_PROCESSOR_NAME];
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
 * however it will output using mdlDiag and the code continues.
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
#define MDLTIMER


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
int mdlInitialize(MDL *,char **,void (*)(MDL),void (*)(MDL));
void mdlFinish(MDL);
int mdlThreads(MDL);
int mdlSelf(MDL);
int mdlOldSelf(MDL);
int mdlSwap(MDL,int,size_t,void *,size_t,size_t *,size_t *);
void mdlDiag(MDL,char *);
void mdlAddService(MDL,int,void *,void (*)(void *,void *,int,void *,int *),
		   int,int);
void mdlReqService(MDL,int,int,void *,int);
void mdlGetReply(MDL,int,void *,int *);
void mdlHandler(MDL);
/*
** Collective operations
*/
/*typedef MPI_Op MDL_Op;
typedef MPI_Datatype MDL_Datatype;
int mdlReduce ( MDL mdl, void *sendbuf, void *recvbuf, int count,
MDL_Datatype datatype, MDL_Op op, int root );*/
/*
** FFT Operations
*/
#ifdef MDL_FFTW
typedef struct mdlFFTContext {
    MDL mdl;
    rfftwnd_plan fplan;
    rfftwnd_plan iplan;

    int rx, ry, rz; /* Real dimensions */
    int sz, nz;     /* Start z and number of z */
    int sy, ny;     /* Transposed start and number */
    int nlocal;     /* Number of local elements */
    } * MDLFFT;

size_t mdlFFTInitialize(MDL mdl,MDLFFT *fft,
			int nx,int ny,int nz,int bMeasure);
void mdlFFTFinish( MDLFFT fft );
void mdlFFT( MDLFFT fft, fftw_real *data, int bInverse );
#endif
/*
 ** Caching functions.
 */
void *mdlMalloc(MDL,size_t);
void mdlFree(MDL,void *);
void mdlROcache(MDL,int,void *,int,int);
void mdlCOcache(MDL,int,void *,int,int,
		void (*)(void *),void (*)(void *,void *));
void mdlFinishCache(MDL,int);
void *mdlAquire(MDL,int,int,int);
void mdlPrefetch(MDL,int,int,int);
void mdlRelease(MDL,int,void *);
void mdlCacheBarrier(MDL,int);
void mdlCacheCheck(MDL);
/*
 ** Cache statistics functions.
 */
double mdlNumAccess(MDL,int);
double mdlMissRatio(MDL,int);
double mdlCollRatio(MDL,int);
double mdlMinRatio(MDL,int);

#endif
