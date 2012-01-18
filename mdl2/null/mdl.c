/*
 ** NULL mdl.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <stdarg.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#ifdef _MSC_VER
#include <WinSock2.h> /* gethostname */
#else
#include <sys/resource.h>
#endif
#include "mdl.h"
const char *null_mdl_module_id = "NULL ($Id$)";

#define MDL_NOCACHE			0
#define MDL_ROCACHE			1
#define MDL_COCACHE			2


#define MDL_DEFAULT_BYTES		4096
#define MDL_DEFAULT_SERVICES	50
#define MDL_DEFAULT_CACHEIDS	5

void _srvNull(void *p1,void *vin,int nIn,void *vout,int *pnOut) {
    return;
    }


double mdlCpuTimer(MDL mdl) {
#ifdef __linux__
    struct rusage ru;

    getrusage(0,&ru);
    return((double)ru.ru_utime.tv_sec + 1e-6*(double)ru.ru_utime.tv_usec);
#else
    return 0.0;
#endif
    }

/*
 * MDL debug and Timer functions
 */
#define MDLPRINTF_STRING_MAXLEN 256
void mdlprintf( MDL mdl, const char *format, ... ) {
    static char ach[MDLPRINTF_STRING_MAXLEN];
    va_list args;

    if (mdl->bDiag) {
	va_start( args, format);
	vsnprintf( ach, MDLPRINTF_STRING_MAXLEN, format, args);
	mdlDiag( mdl, ach);
	va_end( args);
	}
    }

#ifdef MDLDEBUG
void mdldebug( MDL mdl, const char *format, ... ) {
    static char ach[MDLPRINTF_STRING_MAXLEN];
    va_list args;

    if (mdl->bDiag) {
	va_start( args, format);
	vsnprintf( ach, MDLPRINTF_STRING_MAXLEN, format, args);
	mdlDiag( mdl, ach);
	va_end( args);
	}
    }
#endif

#ifdef MDLTIMER
void mdlZeroTimer(MDL mdl, mdlTimer *t) {
#ifdef _MSC_VER
    FILETIME ft;
    uint64_t clock;
    GetSystemTimeAsFileTime(&ft);
    clock = ft.dwHighDateTime;
    clock <<= 32;
    clock |= ft.dwLowDateTime;
    /* clock is in 100 nano-second units */
    t->wallclock = clock / 10000000.0;
#else
    struct timezone tz;
    struct timeval tv;
    struct rusage ru;
    tz.tz_minuteswest = 0;
    tz.tz_dsttime = 0;
    gettimeofday(&tv,&tz);
    t->wallclock = tv.tv_sec + 1e-6*(double) tv.tv_usec;
    getrusage(0,&ru);
    t->cpu = (double)ru.ru_utime.tv_sec + 1e-6*(double)ru.ru_utime.tv_usec;
    t->system = (double)ru.ru_stime.tv_sec + 1e-6*(double)ru.ru_stime.tv_usec;
#endif
    }

void mdlGetTimer(MDL mdl, mdlTimer *t0, mdlTimer *t) {
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

    getrusage(0,&ru);
    t->cpu = (double)ru.ru_utime.tv_sec + 1e-6*(double)ru.ru_utime.tv_usec - t0->cpu;
    t->system = (double)ru.ru_stime.tv_sec + 1e-6*(double)ru.ru_stime.tv_usec - t0->system;
    tz.tz_minuteswest = 0;
    tz.tz_dsttime = 0;
    gettimeofday(&tv,&tz);
    t->wallclock = tv.tv_sec + 1e-6*(double) tv.tv_usec - t0->wallclock;
#endif
}

void mdlPrintTimer(MDL mdl,char *message, mdlTimer *t0) {
    mdlTimer lt;

    if (mdl->bDiag) {
	mdlGetTimer(mdl,t0,&lt);
	mdlprintf(mdl,"%s %f %f %f\n",message,lt.wallclock,lt.cpu,lt.system);
	}
    }
#endif

int mdlInitialize(MDL *pmdl,char **argv,void (*fcnChild)(MDL),void (*fcnIOChild)(MDL)) {
    MDL mdl;
    int i,nThreads,bDiag,bThreads;
    char *p,ach[256],achDiag[256];

    assert(fcnIOChild==NULL);

    *pmdl = NULL;
    mdl = malloc(sizeof(struct mdlContext));
    assert(mdl != NULL);
    /*
     ** Set default "maximums" for structures. These are NOT hard
     ** maximums, as the structures will be realloc'd when these
     ** values are exceeded.
     */
    mdl->nMaxServices = MDL_DEFAULT_SERVICES;
    mdl->nMaxInBytes = MDL_DEFAULT_BYTES;
    mdl->nMaxOutBytes = MDL_DEFAULT_BYTES;
    mdl->nMaxCacheIds = MDL_DEFAULT_CACHEIDS;
    /*
     ** Now allocate the initial service slots.
     */
    mdl->psrv = malloc(mdl->nMaxServices*sizeof(SERVICE));
    assert(mdl->psrv != NULL);
    /*
     ** Initialize the new service slots.
     */
    for (i=0;i<mdl->nMaxServices;++i) {
	mdl->psrv[i].p1 = NULL;
	mdl->psrv[i].nInBytes = 0;
	mdl->psrv[i].nOutBytes = 0;
	mdl->psrv[i].fcnService = NULL;
	}
    /*
     ** Provide a 'null' service for sid = 0, so that stopping the
     ** service handler is well defined!
     */
    mdl->psrv[0].p1 = NULL;
    mdl->psrv[0].nInBytes = 0;
    mdl->psrv[0].nOutBytes = 0;
    mdl->psrv[0].fcnService = _srvNull;
    /*
     ** Allocate service buffers.
     */
    mdl->pszIn = malloc(mdl->nMaxInBytes);
    assert(mdl->pszIn != NULL);
    mdl->pszOut = malloc(mdl->nMaxOutBytes);
    assert(mdl->pszOut != NULL);
    /*
     ** Allocate initial cache spaces.
     */
    mdl->cache = malloc(mdl->nMaxCacheIds*sizeof(CACHE));
    assert(mdl->cache != NULL);
    /*
     ** Initialize caching spaces.
     */
    for (i=0;i<mdl->nMaxCacheIds;++i) {
	mdl->cache[i].iType = MDL_NOCACHE;
	}
    /*
     ** Do some low level argument parsing for number of threads, and
     ** diagnostic flag!
     */
    bDiag = 0;
    bThreads = 0;
    i = 1;
    while (argv[i]) {
	if (!strcmp(argv[i],"-sz") && !bThreads) {
	    ++i;
	    if (argv[i]) {
		nThreads = atoi(argv[i]);
		bThreads = 1;
		}
	    }
	if (!strcmp(argv[i],"+d") && !bDiag) {
	    p = getenv("MDL_DIAGNOSTIC");
	    if (!p) p = getenv("HOME");
	    if (!p) sprintf(ach,"/tmp");
	    else sprintf(ach,"%s",p);
	    bDiag = 1;
	    }
	++i;
	}
    nThreads = 1;
    mdl->bDiag = bDiag;
    mdl->nThreads = nThreads;
    *pmdl = mdl;

    gethostname(mdl->nodeName,sizeof(mdl->nodeName));
    mdl->nodeName[sizeof(mdl->nodeName)-1] = 0;

#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
    mdl->dComputing = 0.0;
    mdl->nTicks = getticks();
#endif

    /*
     ** A unik!
     */
    mdl->idSelf = 0;
    if (mdl->bDiag) {
	char *tmp = strrchr(argv[0],'/');
	if (!tmp) tmp = argv[0];
	else ++tmp;
	sprintf(achDiag,"%s/%s.%d",ach,tmp,mdl->idSelf);
	mdl->fpDiag = fopen(achDiag,"w");
	assert(mdl->fpDiag != NULL);
	}
    return(nThreads);
    }


void mdlFinish(MDL mdl) {
    /*
     ** Close Diagnostic file.
     */
    if (mdl->bDiag) {
	fclose(mdl->fpDiag);
	}
    /*
     ** Deregister from PVM and deallocate storage.
     */
    free(mdl->psrv);
    free(mdl->pszIn);
    free(mdl->pszOut);
    free(mdl->cache);
    free(mdl);
    }


/*
 ** This function returns the number of threads in the set of
 ** threads.
 */
int mdlThreads(MDL mdl) {
    return(mdl->nThreads);
    }


/*
 ** This function returns this threads 'id' number within the specified
 ** MDL Context. Parent thread always has 'id' of 0, where as children
 ** have 'id's ranging from 1..(nThreads - 1).
 */
int mdlSelf(MDL mdl) {
    return(mdl->idSelf);
    }

/*
** This function returns the absolute ID, when using IO threads. This is 
** only really useful for debugging. For now it is equivalent to mdlSelf
** since IO threads are not supported yet.
*/
int mdlOldSelf(MDL mdl) {
    return(mdl->idSelf);
    }

size_t typeSize(MDL_Datatype type) {
    switch(type) {
    case MDL_FLOAT: return sizeof(float);
    case MDL_DOUBLE:return sizeof(double);
    default:
	assert(0);
	return 0;
	}
    }

/*
** Collective operations
*/
int mdlReduce ( MDL mdl, void *sendbuf, void *recvbuf, int count,
		MDL_Datatype datatype, MDL_Op op, int root ) {
    memcpy(recvbuf,sendbuf,count*typeSize(datatype));
    return 0;
    }

int mdlAllreduce ( MDL mdl, void *sendbuf, void *recvbuf, int count,
		   MDL_Datatype datatype, MDL_Op op ) {
    memcpy(recvbuf,sendbuf,count*typeSize(datatype));
    return 0;
    }

int mdlAlltoall( MDL mdl, void *sendbuf, int scount, MDL_Datatype stype,
		 void *recvbuf, int rcount, MDL_Datatype rtype) {
    assert(rcount==scount && stype==rtype);
    memcpy(recvbuf,sendbuf,scount*typeSize(stype));
    return 0;
    }




/*
 ** This is a tricky function. It initiates a bilateral transfer between
 ** two threads. Both threads MUST be expecting this transfer. The transfer
 ** occurs between idSelf <---> 'id' or 'id' <---> idSelf as seen from the
 ** opposing thread. It is designed as a high performance non-local memory
 ** swapping primitive and implementation will vary in non-trivial ways
 ** between differing architectures and parallel paradigms (eg. message
 ** passing and shared address space). A buffer is specified by 'pszBuf'
 ** which is 'nBufBytes' in size. Of this buffer the LAST 'nOutBytes' are
 ** transfered to the opponent, in turn, the opponent thread transfers his
 ** nBufBytes to this thread's buffer starting at 'pszBuf'.
 ** If the transfer completes with no problems the function returns 1.
 ** If the function returns 0 then one of the players has not received all
 ** of the others memory, however he will have successfully transfered all
 ** of his memory.
 */
int mdlSwap(MDL mdl,int id,size_t nBufBytes,void *vBuf,size_t nOutBytes,
	    size_t *pnSndBytes,size_t *pnRcvBytes) {
    assert(0);
    return 0;
    }

const char *mdlName(MDL mdl) {
    return mdl->nodeName;
    }

void mdlDiag(MDL mdl,char *psz) {
    if (mdl->bDiag) {
	fputs(psz,mdl->fpDiag);
	fflush(mdl->fpDiag);
	}
    }


void mdlAddService(MDL mdl,int sid,void *p1,
		   void (*fcnService)(void *,void *,int,void *,int *),
		   int nInBytes,int nOutBytes) {
    int i,nMaxServices;

    assert(sid > 0);
    if (sid >= mdl->nMaxServices) {
	/*
	 ** reallocate service buffer, adding space for 8 new services
	 ** including the one just defined.
	 */
	nMaxServices = sid + 9;
	mdl->psrv = realloc(mdl->psrv,nMaxServices*sizeof(SERVICE));
	assert(mdl->psrv != NULL);
	/*
	 ** Initialize the new service slots.
	 */
	for (i=mdl->nMaxServices;i<nMaxServices;++i) {
	    mdl->psrv[i].p1 = NULL;
	    mdl->psrv[i].nInBytes = 0;
	    mdl->psrv[i].nOutBytes = 0;
	    mdl->psrv[i].fcnService = NULL;
	    }
	mdl->nMaxServices = nMaxServices;
	}
    /*
     ** Make sure the service buffers are big enough!
     */
    if (nInBytes > mdl->nMaxInBytes) {
	mdl->pszIn = realloc(mdl->pszIn,nInBytes);
	assert(mdl->pszIn != NULL);
	mdl->nMaxInBytes = nInBytes;
	}
    if (nOutBytes > mdl->nMaxOutBytes) {
	mdl->pszOut = realloc(mdl->pszOut,nOutBytes);
	assert(mdl->pszOut != NULL);
	mdl->nMaxOutBytes = nOutBytes;
	}
    mdl->psrv[sid].p1 = p1;
    mdl->psrv[sid].nInBytes = nInBytes;
    mdl->psrv[sid].nOutBytes = nOutBytes;
    mdl->psrv[sid].fcnService = fcnService;
    }


void mdlReqService(MDL mdl,int id,int sid,void *vin,int nInBytes) {
    assert(0);
    }


void mdlGetReply(MDL mdl,int id,void *vout,int *pnOutBytes) {
    assert(0);
    }


void mdlHandler(MDL mdl) {
    assert(0);
    }


#define BILLION				1000000000

/*
 ** Special MDL memory allocation functions for allocating memory
 ** which must be visible to other processors thru the MDL cache
 ** functions.
 ** mdlMalloc() is defined to return a pointer to AT LEAST iSize bytes
 ** of memory. This pointer will be passed to either mdlROcache or
 ** mdlCOcache as the pData parameter.
 ** For PVM and most machines these functions are trivial, but on the
 ** T3D and perhaps some future machines these functions are required.
 */
void *mdlMalloc(MDL mdl,size_t iSize) {
    return(malloc(iSize));
    }


void mdlFree(MDL mdl,void *p) {
    free(p);
    }

/*
** This is the default element fetch routine.  It impliments the old behaviour
** of a single large array.  New data structures need to be more clever.
*/
static void *getArrayElement(void *vData,int i,int iDataSize) {
    char *pData = vData;
    return pData + i*iDataSize;
    }

void mdlSetCacheSize(MDL mdl,int cacheSize) {
    }

/*
 ** Common initialization for all types of caches.
 */
CACHE *CacheInitialize(MDL mdl,int cid,
    void * (*getElt)(void *pData,int i,int iDataSize),
    void *pData,int iDataSize,int nData) {
    CACHE *c;
    int i,nMaxCacheIds;

    /*
     ** Allocate more cache spaces if required!
     */
    assert(cid >= 0);
    if (cid >= mdl->nMaxCacheIds) {
	/*
	 ** reallocate cache spaces, adding space for 2 new cache spaces
	 ** including the one just defined.
	 */
	nMaxCacheIds = cid + 3;
	mdl->cache = realloc(mdl->cache,nMaxCacheIds*sizeof(CACHE));
	assert(mdl->cache != NULL);
	/*
	 ** Initialize the new cache slots.
	 */
	for (i=mdl->nMaxCacheIds;i<nMaxCacheIds;++i) {
	    mdl->cache[i].iType = MDL_NOCACHE;
	    }
	mdl->nMaxCacheIds = nMaxCacheIds;
	}
    c = &mdl->cache[cid];
    assert(c->iType == MDL_NOCACHE);
    c->getElt = getElt==NULL ? getArrayElement : getElt;
    c->pData = pData;
    c->iDataSize = iDataSize;
    c->nData = nData;
    c->nAccess = 0;
    c->nAccHigh = 0;
    c->nMiss = 0;
    c->nColl = 0;
    c->nMin = 0;
    return(c);
    }


/*
 ** Initialize a Read-Only caching space.
 */

void mdlROcache(MDL mdl,int cid,
    void * (*getElt)(void *pData,int i,int iDataSize),
    void *pData,int iDataSize,int nData) {
    CACHE *c;

    c = CacheInitialize(mdl,cid,getElt,pData,iDataSize,nData);
    c->iType = MDL_ROCACHE;
    c->init = NULL;
    c->combine = NULL;
    }


/*
 ** Initialize a Combiner caching space.
 */
void mdlCOcache(MDL mdl,int cid,
    void * (*getElt)(void *pData,int i,int iDataSize),
    void *pData,int iDataSize,int nData,
    void *ctx,void (*init)(void *,void *),void (*combine)(void *,void *,void *)) {
    CACHE *c;

    c = CacheInitialize(mdl,cid,getElt,pData,iDataSize,nData);
    c->iType = MDL_COCACHE;
    c->init = init;
    c->combine = combine;
    c->ctx = ctx;
    }


void mdlFinishCache(MDL mdl,int cid) {
    CACHE *c = &mdl->cache[cid];

#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
	{
	ticks nTicks = getticks();
	mdl->dComputing += elapsed( nTicks, mdl->nTicks );
	mdl->nTicks = nTicks;
	}
#endif

    /*
     ** Free up storage and finish.
     */
    c->iType = MDL_NOCACHE;
    }

void mdlCacheBarrier(MDL mdl,int cid) {
    }

void mdlCacheCheck(MDL mdl) {
    }

void mdlPrefetch(MDL mdl,int cid,int iIndex,int id) {
    /* The data is already local so there is nothing to do */
    }

void *mdlAquire(MDL mdl,int cid,int iIndex,int id) {
    CACHE *c = &mdl->cache[cid];

    ++c->nAccess;
    assert(id == mdl->idSelf);
    return (*c->getElt)(c->pData,iIndex,c->iDataSize);
    }


void mdlRelease(MDL mdl,int cid,void *p) {
    }


double mdlNumAccess(MDL mdl,int cid) {
    CACHE *c = &mdl->cache[cid];

    return(c->nAccHigh*1e9 + c->nAccess);
    }


double mdlMissRatio(MDL mdl,int cid) {
    CACHE *c = &mdl->cache[cid];
    double dAccess = c->nAccHigh*1e9 + c->nAccess;

    if (dAccess > 0.0) return(c->nMiss/dAccess);
    else return(0.0);
    }


double mdlCollRatio(MDL mdl,int cid) {
    CACHE *c = &mdl->cache[cid];
    double dAccess = c->nAccHigh*1e9 + c->nAccess;

    if (dAccess > 0.0) return(c->nColl/dAccess);
    else return(0.0);
    }


double mdlMinRatio(MDL mdl,int cid) {
    CACHE *c = &mdl->cache[cid];
    double dAccess = c->nAccHigh*1e9 + c->nAccess;

    if (dAccess > 0.0) return(c->nMin/dAccess);
    else return(0.0);
    }

#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
void mdlTimeReset(MDL mdl) {
    mdl->dComputing = 0.0;
    mdl->nTicks = getticks();
    }

static double TimeFraction(MDL mdl) {
    double dTotal = mdl->dComputing;
    if ( dTotal <= 0.0 ) return 0.0;
    return 100.0 / dTotal;
    }

double mdlTimeComputing(MDL mdl) {
    return mdl->dComputing * TimeFraction(mdl);
    }

double mdlTimeSynchronizing(MDL mdl) {
    return 0.0;
    }

double mdlTimeWaiting(MDL mdl) {
    return 0.0;
    }
#endif

/*
** GRID Geometry information.  The basic process is as follows:
** - Initialize: Create a MDLGRID giving the global geometry information (total grid size)
** - SetLocal:   Set the local grid geometry (which slabs are on this processor)
** - GridShare:  Share this information between processors
** - Malloc:     Allocate one or more grid instances
** - Free:       Free the memory for all grid instances
** - Finish:     Free the GRID geometry information.
*/
void mdlGridInitialize(MDL mdl,MDLGRID *pgrid,int n1,int n2,int n3,int a1) {
    MDLGRID grid;
    assert(n1>0&&n2>0&&n3>0);
    assert(a1<=n1);
    *pgrid = grid = malloc(sizeof(struct mdlGridContext)); assert(grid!=NULL);
    grid->n1 = n1;
    grid->n2 = n2;
    grid->n3 = n3;
    grid->a1 = a1;
    }

void mdlGridFinish(MDL mdl, MDLGRID grid) {
    free(grid);
    }

void mdlGridSetLocal(MDL mdl,MDLGRID grid,int s, int n, int nlocal) {
    assert(s==0);
    assert(s+n==grid->n3);
    grid->nlocal = nlocal;
    }

/*
** Share the local GRID information with other processors by,
**   - finding the starting slab and number of slabs on each processor
**   - building a mapping from slab to processor id.
*/
void mdlGridShare(MDL mdl,MDLGRID grid) {
    }

/*
** Allocate the local elements.  The size of a single element is
** given and the local GRID information is consulted to determine
** how many to allocate.
*/
void *mdlGridMalloc(MDL mdl,MDLGRID grid,int nEntrySize) {
    return mdlMalloc(mdl,nEntrySize*grid->nlocal);
    }

void mdlGridFree( MDL mdl, MDLGRID grid, void *p ) {
    mdlFree(mdl,p);
    }

#ifdef MDL_FFTW
size_t mdlFFTInitialize(MDL mdl,MDLFFT *pfft,
			int n1,int n2,int n3,int bMeasure) {
    MDLFFT fft;
    *pfft = NULL;
    fft = malloc(sizeof(struct mdlFFTContext));
    assert(fft != NULL);

    fft->fplan = rfftw3d_create_plan( n3, n2, n1,
				      FFTW_REAL_TO_COMPLEX,
				      FFTW_IN_PLACE | (bMeasure ? FFTW_MEASURE : FFTW_ESTIMATE) );

    fft->iplan = rfftw3d_create_plan(/* dim.'s of REAL data --> */ n3, n2, n1,
				     FFTW_COMPLEX_TO_REAL,
				     FFTW_IN_PLACE | (bMeasure ? FFTW_MEASURE : FFTW_ESTIMATE));

    mdlGridInitialize(mdl,&fft->rgrid,n1,n2,n3,2*(n1/2+1));
    mdlGridInitialize(mdl,&fft->kgrid,n1/2+1,n3,n2,n1/2+1);



    *pfft = fft;
    return n1 * n2 * n3;
    }

void mdlFFTFinish( MDL mdl, MDLFFT fft ) {
    rfftwnd_destroy_plan(fft->fplan);
    rfftwnd_destroy_plan(fft->iplan);
    free(fft);
    }

void mdlFFT( MDL mdl, MDLFFT fft, fftw_real *data, int bInverse ) {
    rfftwnd_plan plan = bInverse ? fft->iplan : fft->fplan;
    rfftwnd_one(plan,data,0); /* ,FFTW_TRANSPOSED_ORDER); */
    }
#endif

void mdlSetWorkQueueSize(MDL mdl,int wqSize,int cudaSize) {
    }

/* Just do the work immediately */
void mdlAddWork(MDL mdl, void *ctx, mdlWorkFunction initWork, mdlWorkFunction checkWork, mdlWorkFunction doWork, mdlWorkFunction doneWork) {
    while( doWork(ctx) != 0 ) {}
    doneWork(ctx);
    }
