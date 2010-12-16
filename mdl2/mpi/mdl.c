/*
 ** MPI version of MDL.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <stdarg.h>
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#ifdef __linux__
#include <sys/resource.h>
#endif
#include "mpi.h"
#include "mdl.h"
#ifdef USE_BSC
#include "mpitrace_user_events.h"
#endif

const char *mpi_mdl_module_id = "MPI ($Id$)";

#define MDL_NOCACHE			0
#define MDL_ROCACHE			1
#define MDL_COCACHE			2

#define MDL_DEFAULT_BYTES		80000
#define MDL_DEFAULT_SERVICES	120
#define MDL_DEFAULT_CACHEIDS	5

#define MDL_TRANS_SIZE		5000000
#define MDL_TAG_INIT 		1
#define MDL_TAG_SWAPINIT 	2
#define MDL_TAG_SWAP		3
#define MDL_TAG_REQ	   	4
#define MDL_TAG_RPL		5
#define MDL_TAG_SEND            6

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
static inline int size_t_to_int( size_t v ) {
    return (int)v;
    }

static inline int mdlkey_t_to_int( mdlkey_t v ) {
    return (int)v;
    }


/*
 ** This structure should be "maximally" aligned, with 4 ints it
 ** should align up to at least QUAD word, which should be enough.
 */
typedef struct srvHeader {
    int idFrom;
    int sid;
    int nInBytes;
    int nOutBytes;
    } SRVHEAD;

void _srvNull(void *p1,void *vin,int nIn,void *vout,int *pnOut) {
    return;
    }


double mdlCpuTimer(MDL mdl) {
#ifdef __linux__
    struct rusage ru;

    getrusage(0,&ru);
    return((double)ru.ru_utime.tv_sec + 1e-6*(double)ru.ru_utime.tv_usec);
#else
    return( ((double) clock())/CLOCKS_PER_SEC);
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
    t->wallclock = clock / 10000000UL;
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

#if 0
static inline unsigned int swap( unsigned int v, int N ) {
    unsigned int r = 0;
    do {
	r <<= 1;
	r |= ( v&1 );
	v >>= 1;
	}
    while ( N >>= 1 );
    return r;
    }

static int order( int N, int T ) {
    int K, V, R=N;
    int rank = 0;

    assert( T < N );

    for ( K=0; R; K++ ) {

	V = swap(K,N);
	if ( V < N ) {
	    if ( V == T ) return rank;
	    rank++;
	    R--;
	    }

	}
    assert(0);
    }
#endif

/*
** This will remove nProcs processors from communicator zero
** and reconfigure the remaining processors.  commList[] is:
**   [0] - Our private communicator
**   [1]
*/
int  mdlSplitComm(MDL mdl, int nProcs) {
    int stride, offset, oldRank, isIO, peerID, iComm;
    MPI_Comm newComm;

    if ( nProcs == 0 ) return 0;

    assert( nProcs < mdl->nThreads );
    assert( mdl->commCount < MDL_MAX_COMM );

    oldRank = mdl->idSelf;

    isIO = 0;
    stride = mdl->nThreads / nProcs;
    offset = mdl->nThreads - nProcs * stride;
    mdl->idSelf = (oldRank - offset) / stride;  /* I/O processor rank */
    if ( mdl->idSelf < 0 )
	mdl->idSelf = oldRank;
    else if ( (oldRank - offset) % stride != stride-1 )
	mdl->idSelf = oldRank - mdl->idSelf;
    else
	isIO = 1;
    mdl->nThreads -= nProcs;
    peerID = isIO ? 0 : stride+offset-1;

    MPI_Comm_split(mdl->commList[0], isIO, mdl->idSelf, &newComm );
    MPI_Intercomm_create( newComm, 0, mdl->commList[0],
			  peerID, 10, &mdl->commList[mdl->commCount++]);
    if ( mdl->commList[0] != MPI_COMM_WORLD )
	MPI_Comm_free(&mdl->commList[0]);
    mdl->commMDL = mdl->commList[0] = newComm;
    MPI_Comm_rank(mdl->commMDL, &mdl->idSelf);

    iComm = isIO ? mdl->commCount-1 : 0;
    mdl->commMDL = mdl->commList[iComm];
    return iComm;
    }

int mdlInitialize(MDL *pmdl,char **argv,
		  void (*fcnChild)(MDL),
		  void (*fcnIO)(MDL)) {
    MDL mdl;
    int iComm;
    int i,j,bDiag,bThreads;
    char *p,ach[256],achDiag[256];
    int argc;

    *pmdl = NULL;
    mdl = malloc(sizeof(struct mdlContext));
    assert(mdl != NULL);
    /*
     ** Set default "maximums" for structures. These are NOT hard
     ** maximums, as the structures will be realloc'd when these
     ** values are exceeded.
     */
    mdl->nMaxServices = MDL_DEFAULT_SERVICES;
    mdl->nMaxSrvBytes = MDL_DEFAULT_BYTES;
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
    mdl->pszIn = malloc(mdl->nMaxSrvBytes+sizeof(SRVHEAD));
    assert(mdl->pszIn != NULL);
    mdl->pszOut = malloc(mdl->nMaxSrvBytes+sizeof(SRVHEAD));
    assert(mdl->pszOut != NULL);
    mdl->pszBuf = malloc(mdl->nMaxSrvBytes+sizeof(SRVHEAD));
    assert(mdl->pszBuf != NULL);
    /*
     ** Allocate swapping transfer buffer. This buffer remains fixed.
     */
    mdl->pszTrans = malloc(MDL_TRANS_SIZE);
    assert(mdl->pszTrans != NULL);
    /*
     ** Allocate initial cache spaces.
     */
    mdl->cache = malloc(mdl->nMaxCacheIds*sizeof(CACHE));
    assert(mdl->cache != NULL);
    mdl->freeCacheLines = NULL;
    /*
     ** Initialize caching spaces.
     */
    mdl->cacheSize = MDL_CACHE_SIZE;
    for (i=0;i<mdl->nMaxCacheIds;++i) {
	mdl->cache[i].iType = MDL_NOCACHE;
	}

    for (argc = 0; argv[argc]; argc++);

    MPI_Init(&argc, &argv);

#ifdef USE_BSC_trace
    MPItrace_shutdown();
#endif

    /* We start with a single communicator */
    mdl->commCount = 1;
    mdl->commList[0] = MPI_COMM_WORLD;

    /*
     ** Do some low level argument parsing for number of threads, and
     ** diagnostic flag!
     */
    bDiag = 0;
    bThreads = 0;
    mdl->nIO = 0;
    i = 1;
    while (argv[i]) {
	if (!strcmp(argv[i],"-io") && fcnIO!=0) {
	    if (!argv[i+1]) {
		fprintf(stderr, "-io requires an argument\n");
		exit(0);
		}
	    mdl->nIO = atoi(argv[i+1]);
	    j = i;
	    do {
		argv[j] = argv[j+2];
		}
	    while ( argv[j++] );
	    continue;
	    }
	if (!strcmp(argv[i],"-sz") && !bThreads) {
	    ++i;
	    if (argv[i]) bThreads = 1;
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
    if (bThreads) {
	fprintf(stderr,"Warning: -sz parameter ignored, using as many\n");
	fprintf(stderr,"         processors as specified in environment.\n");
	fflush(stderr);
	}

    MPI_Comm_size(mdl->commList[0], &mdl->nThreads);
    MPI_Comm_rank(mdl->commList[0], &mdl->idSelf);

    mdlassert( mdl, mdl->nIO >= 0 && mdl->nIO <= MDL_MAX_IO_PROCS );

    /*
    ** It is important that every processor know the IO processor count.
    ** MPI only passes parameters to rank 0, so we sent it on.
    */
    MPI_Bcast( &(mdl->nIO), 1, MPI_INT, 0, mdl->commList[0]);

    iComm = mdlSplitComm(mdl,mdl->nIO);
    mdl->commMDL = mdl->commList[iComm];

    /*
     ** Allocate caching buffers, with initial data size of 0.
     ** We need one reply buffer for each thread, to deadlock situations.
     */
    mdl->iMaxDataSize = 0;
    mdl->iCaBufSize = sizeof(CAHEAD);
    mdl->pszRcv = malloc(mdl->iCaBufSize);
    assert(mdl->pszRcv != NULL);
    mdl->ppszRpl = malloc(mdl->nThreads*sizeof(char *));
    assert(mdl->ppszRpl != NULL);
    mdl->pmidRpl = malloc(mdl->nThreads*sizeof(int));
    assert(mdl->pmidRpl != NULL);
    for (i=0;i<mdl->nThreads;++i)
	mdl->pmidRpl[i] = -1;
    mdl->pReqRpl = malloc(mdl->nThreads*sizeof(MPI_Request));
    assert(mdl->pReqRpl != NULL);
    for (i=0;i<mdl->nThreads;++i) {
	mdl->ppszRpl[i] = malloc(mdl->iCaBufSize);
	assert(mdl->ppszRpl[i] != NULL);
	}
    mdl->pszFlsh = malloc(mdl->iCaBufSize);
    assert(mdl->pszFlsh != NULL);
    mdl->bDiag = bDiag;
    *pmdl = mdl;
    if (mdl->bDiag) {
	char *tmp = strrchr(argv[0],'/');
	if (!tmp) tmp = argv[0];
	else ++tmp;
	sprintf(achDiag,"%s/%s.%d",ach,tmp,mdlOldSelf(mdl));
	mdl->fpDiag = fopen(achDiag,"w");
	assert(mdl->fpDiag != NULL);
	}

    MPI_Get_processor_name( mdl->nodeName, &i );
    mdl->nodeName[i] = 0;

    /*if ( oldRank >= mdl->nThreads ) {*/
    if ( iComm ) {
	/*
	 ** I/O thread.
	 */
	mdlassert( mdl, mdl->nIO > 0 );
	(*fcnIO)(mdl);
	mdlFinish(mdl);
	exit(0);
	}
    else if (mdl->nThreads > 1 && mdlOldSelf(mdl)) {
	/*
	 ** Child thread.
	 */
	(*fcnChild)(mdl);
	mdlFinish(mdl);
	exit(0);
	}
    return(mdl->nThreads);
    }


void mdlFinish(MDL mdl) {
    int i;
    FREECACHELINES *fcl, *fcn;

    MPI_Barrier(mdl->commMDL);
    for ( i=0; i<mdl->commCount; i++ )
	if ( mdl->commList[i] != MPI_COMM_WORLD )
	    MPI_Comm_free(&mdl->commList[i]);

    MPI_Finalize();
    /*
     ** Close Diagnostic file.
     */
    if (mdl->bDiag) {
	fclose(mdl->fpDiag);
	}
    /*
     ** Deallocate storage.
     */
    for(fcl=mdl->freeCacheLines;fcl!=NULL;fcl=fcn) {
	fcn = fcl->next;
	free(fcl);
	}
    free(mdl->psrv);
    free(mdl->pszIn);
    free(mdl->pszOut);
    free(mdl->pszBuf);
    free(mdl->pszTrans);
    free(mdl->cache);
    free(mdl->pszRcv);
    free(mdl->pszFlsh);
    for (i=0;i<mdl->nThreads;++i) free(mdl->ppszRpl[i]);
    free(mdl->ppszRpl);
    free(mdl->pmidRpl);
    free(mdl->pReqRpl);
    free(mdl);
    }


void mdlSetComm(MDL mdl, int iComm) {
    assert( iComm >=0 && iComm < mdl->commCount );
    mdl->commMDL = mdl->commList[iComm];
    }


/*
 ** This function returns the number of threads in the set of
 ** threads.
 */
int mdlThreads(MDL mdl) {
    return(mdl->nThreads);
    }


/*
 ** This function returns the number of threads dedicated to IO.
 */
int mdlIO(MDL mdl) {
    return(mdl->nIO);
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
** This function returns the absolute ID (the MPI rank), useful for
** debugging any problems that might occur.
*/
int mdlOldSelf(MDL mdl) {
    int oldSelf;
    MPI_Comm_rank(MPI_COMM_WORLD, &oldSelf);
    return oldSelf;
    }


const char *mdlName(MDL mdl) {
    return mdl->nodeName;
    }


/*
** This needs to be improved by abstracting away more of the MPI functionality
*/
int mdlReduce ( MDL mdl, void *sendbuf, void *recvbuf, int count,
		MDL_Datatype datatype, MDL_Op op, int root ) {
    return MPI_Reduce( sendbuf, recvbuf, count, datatype, op, root, mdl->commMDL );
    }

int mdlAllreduce ( MDL mdl, void *sendbuf, void *recvbuf, int count,
		MDL_Datatype datatype, MDL_Op op ) {
    return MPI_Allreduce( sendbuf, recvbuf, count, datatype, op, mdl->commMDL );
    }

int mdlAlltoall( MDL mdl, void *sendbuf, int scount, MDL_Datatype stype,
		 void *recvbuf, int rcount, MDL_Datatype rtype) {
    return MPI_Alltoall(sendbuf,scount,stype,
			recvbuf,rcount,rtype,mdl->commMDL);
    }
/*
** This function will transfer a block of data using a pack function.
** The corresponding node must call mdlRecv.
 */

#define SEND_BUFFER_SIZE (1*1024*1024)

void mdlSend(MDL mdl,int id,mdlPack pack, void *ctx) {
    size_t nBuff;
    char *vOut;

    vOut = malloc(SEND_BUFFER_SIZE);
    mdlassert(mdl,vOut!=NULL);

    do {
	nBuff = (*pack)(ctx,&id,SEND_BUFFER_SIZE,vOut);
	if ( nBuff != 0 ) {
	    MPI_Ssend(vOut,nBuff,MPI_BYTE,id,MDL_TAG_SEND,mdl->commMDL);
	    }
	}
    while ( nBuff != 0 );

    free(vOut);
    }

void mdlRecv(MDL mdl,int id,mdlPack unpack, void *ctx) {
    void *vIn;
    size_t nUnpack;
    int nBytes;
    MPI_Status status;
    int inid;

    if ( id < 0 ) id = MPI_ANY_SOURCE;

    vIn = malloc(SEND_BUFFER_SIZE);
    mdlassert(mdl,vIn!=NULL);

    do {
	MPI_Recv(vIn,SEND_BUFFER_SIZE,MPI_BYTE,id,MDL_TAG_SEND,
		 mdl->commMDL,&status);
	MPI_Get_count(&status, MPI_BYTE, &nBytes);
	inid = status.MPI_SOURCE;
	nUnpack = (*unpack)(ctx,&inid,nBytes,vIn);
	}
    while (nUnpack>0);

    free(vIn);
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
    size_t nInBytes,nOutBufBytes;
    size_t i;
    int nInMax,nOutMax,nBytes;
    int iTag, pid;
    char *pszBuf = vBuf;
    char *pszIn,*pszOut;
    struct swapInit {
	size_t nOutBytes;
	size_t nBufBytes;
	} swi,swo;
    MPI_Status status;
    MPI_Request request;

    *pnRcvBytes = 0;
    *pnSndBytes = 0;
    /*
     **	Send number of rejects to target thread amount of free space
     */
    swi.nOutBytes = nOutBytes;
    swi.nBufBytes = nBufBytes;
    MPI_Isend(&swi,sizeof(swi),MPI_BYTE,id,MDL_TAG_SWAPINIT,
	      mdl->commMDL, &request);
    /*
     ** Receive the number of target thread rejects and target free space
     */
    iTag = MDL_TAG_SWAPINIT;
    pid = id;
    MPI_Recv(&swo,sizeof(swo),MPI_BYTE,pid,iTag,mdl->commMDL, &status);
    MPI_Get_count(&status, MPI_BYTE, &nBytes);
    assert(nBytes == sizeof(swo));
    MPI_Wait(&request, &status);
    nInBytes = swo.nOutBytes;
    nOutBufBytes = swo.nBufBytes;
    /*
     ** Start bilateral transfers. Note: One processor is GUARANTEED to
     ** complete all its transfers.
     */
    assert(nBufBytes >= nOutBytes);
    pszOut = &pszBuf[nBufBytes-nOutBytes];
    pszIn = pszBuf;
    while (nOutBytes && nInBytes) {
	/*
	 ** nOutMax is the maximum number of bytes allowed to be sent
	 ** nInMax is the number of bytes which will be received.
	 */
	nOutMax = size_t_to_int((nOutBytes < MDL_TRANS_SIZE)?nOutBytes:MDL_TRANS_SIZE);
	nOutMax = size_t_to_int((nOutMax < nOutBufBytes)?nOutMax:nOutBufBytes);
	nInMax = size_t_to_int((nInBytes < MDL_TRANS_SIZE)?nInBytes:MDL_TRANS_SIZE);
	nInMax = size_t_to_int((nInMax < nBufBytes)?nInMax:nBufBytes);
	/*
	 ** Copy to a temp buffer to be safe.
	 */
	for (i=0;i<nOutMax;++i) mdl->pszTrans[i] = pszOut[i];
	MPI_Isend(mdl->pszTrans,nOutMax,MPI_BYTE,id,MDL_TAG_SWAP,
		  mdl->commMDL, &request);
	iTag = MDL_TAG_SWAP;
	pid = id;
	MPI_Recv(pszIn,nInMax,MPI_BYTE,pid,iTag,mdl->commMDL,
		 &status);
	MPI_Get_count(&status, MPI_BYTE, &nBytes);
	assert(nBytes == nInMax);
	MPI_Wait(&request, &status);
	/*
	 ** Adjust pointers and counts for next itteration.
	 */
	pszOut = &pszOut[nOutMax];
	nOutBytes -= nOutMax;
	nOutBufBytes -= nOutMax;
	*pnSndBytes += nOutMax;
	pszIn = &pszIn[nInMax];
	nInBytes -= nInMax;
	nBufBytes -= nInMax;
	*pnRcvBytes += nInMax;
	}
    /*
     ** At this stage we perform only unilateral transfers, and here we
     ** could exceed the opponent's storage capacity.
     ** Note: use of Ssend is mandatory here, also because of this we
     ** don't need to use the intermediate buffer mdl->pszTrans.
     */
    while (nOutBytes && nOutBufBytes) {
	nOutMax = size_t_to_int((nOutBytes < MDL_TRANS_SIZE)?nOutBytes:MDL_TRANS_SIZE);
	nOutMax = size_t_to_int((nOutMax < nOutBufBytes)?nOutMax:nOutBufBytes);
	MPI_Ssend(pszOut,nOutMax,MPI_BYTE,id,MDL_TAG_SWAP,
		  mdl->commMDL);
	pszOut = &pszOut[nOutMax];
	nOutBytes -= nOutMax;
	nOutBufBytes -= nOutMax;
	*pnSndBytes += nOutMax;
	}
    while (nInBytes && nBufBytes) {
	nInMax = size_t_to_int((nInBytes < MDL_TRANS_SIZE)?nInBytes:MDL_TRANS_SIZE);
	nInMax = size_t_to_int((nInMax < nBufBytes)?nInMax:nBufBytes);
	iTag = MDL_TAG_SWAP;
	MPI_Recv(pszIn,nInMax,MPI_BYTE,id,iTag,mdl->commMDL,
		 &status);
	MPI_Get_count(&status, MPI_BYTE, &nBytes);
	assert(nBytes == nInMax);
	pszIn = &pszIn[nInMax];
	nInBytes -= nInMax;
	nBufBytes -= nInMax;
	*pnRcvBytes += nInMax;
	}
    if (nOutBytes) return(0);
    else if (nInBytes) return(0);
    else return(1);
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
    int i,nMaxServices,nMaxBytes;

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
    nMaxBytes = (nInBytes > nOutBytes)?nInBytes:nOutBytes;
    if (nMaxBytes > mdl->nMaxSrvBytes) {
	mdl->pszIn = realloc(mdl->pszIn,nMaxBytes+sizeof(SRVHEAD));
	assert(mdl->pszIn != NULL);
	mdl->pszOut = realloc(mdl->pszOut,nMaxBytes+sizeof(SRVHEAD));
	assert(mdl->pszOut != NULL);
	mdl->pszBuf = realloc(mdl->pszBuf,nMaxBytes+sizeof(SRVHEAD));
	assert(mdl->pszBuf != NULL);
	mdl->nMaxSrvBytes = nMaxBytes;
	}
    mdl->psrv[sid].p1 = p1;
    mdl->psrv[sid].nInBytes = nInBytes;
    mdl->psrv[sid].nOutBytes = nOutBytes;
    mdl->psrv[sid].fcnService = fcnService;
    }


void mdlReqService(MDL mdl,int id,int sid,void *vin,int nInBytes) {
    char *pszIn = vin;
    /*
     ** If this looks like dangerous magic, it's because it is!
     */
    SRVHEAD *ph = (SRVHEAD *)mdl->pszBuf;
    char *pszOut = &mdl->pszBuf[sizeof(SRVHEAD)];
    int i;

    ph->idFrom = mdl->idSelf;
    ph->sid = sid;
    if (!pszIn) ph->nInBytes = 0;
    else ph->nInBytes = nInBytes;
    if (nInBytes > 0 && pszIn != NULL) {
	for (i=0;i<nInBytes;++i) pszOut[i] = pszIn[i];
	}
    MPI_Send(mdl->pszBuf,nInBytes+(int)sizeof(SRVHEAD),MPI_BYTE,id,MDL_TAG_REQ,
	     mdl->commMDL);
    }


void mdlGetReply(MDL mdl,int id,void *vout,int *pnOutBytes) {
    char *pszOut = vout;
    SRVHEAD *ph = (SRVHEAD *)mdl->pszBuf;
    char *pszIn = &mdl->pszBuf[sizeof(SRVHEAD)];
    int i,iTag,nBytes;
    MPI_Status status;

    iTag = MDL_TAG_RPL;
    MPI_Recv(mdl->pszBuf,mdl->nMaxSrvBytes+(int)sizeof(SRVHEAD),MPI_BYTE,
	     id,iTag,mdl->commMDL, &status);
    MPI_Get_count(&status, MPI_BYTE, &nBytes);
    assert(nBytes == ph->nOutBytes + sizeof(SRVHEAD));
    if (ph->nOutBytes > 0 && pszOut != NULL) {
	for (i=0;i<ph->nOutBytes;++i) pszOut[i] = pszIn[i];
	}
    if (pnOutBytes) *pnOutBytes = ph->nOutBytes;
    }

void mdlStop(MDL mdl) {
    int id;

    /* Stop the worker processes */
    for ( id=1; id<mdl->nThreads; ++id ) {
	/*if (msr->param.bVDetails)*/
	printf("Stopping worker thread %d\n",id);
	mdlReqService(mdl,id,SRV_STOP,NULL,0);
	mdlGetReply(mdl,id,NULL,NULL);
	}
    if ( mdl->nIO ) {
	printf("Stopping I/O threads\n");
	mdlSetComm(mdl,1);
	mdlReqService(mdl,0,SRV_STOP,NULL,0);
	mdlGetReply(mdl,0,NULL,NULL);
	mdlSetComm(mdl,0);
	}

    printf( "MDL terminated\n" );
    }

static int mdlHandleOne(MDL mdl) {
    SRVHEAD *phi = (SRVHEAD *)mdl->pszIn;
    SRVHEAD *pho = (SRVHEAD *)mdl->pszOut;
    char *pszIn = &mdl->pszIn[sizeof(SRVHEAD)];
    char *pszOut = &mdl->pszOut[sizeof(SRVHEAD)];
    int sid,iTag,id,nOutBytes,nBytes;
    MPI_Status status;
    MPI_Comm   comm;

    /* Save this communicator... reply ALWAYS goes here */
    comm = mdl->commMDL;
    iTag = MDL_TAG_REQ;
    id = MPI_ANY_SOURCE;
    MPI_Recv(mdl->pszIn,mdl->nMaxSrvBytes+sizeof(SRVHEAD),
	     MPI_BYTE, id,iTag,comm,&status);
    /*
    ** Quite a few sanity checks follow.
    */
    id = status.MPI_SOURCE;
    MPI_Get_count(&status, MPI_BYTE, &nBytes);
    assert(nBytes == phi->nInBytes + sizeof(SRVHEAD));
    assert(id == phi->idFrom);
    sid = phi->sid;
    assert(sid < mdl->nMaxServices);
    if (phi->nInBytes > mdl->psrv[sid].nInBytes) {
	printf( "ERROR: pid=%d, sid=%d, nInBytes=%d, sid.nInBytes=%d\n",
		mdlSelf(mdl), sid, phi->nInBytes, mdl->psrv[sid].nInBytes );
	}
    assert(phi->nInBytes <= mdl->psrv[sid].nInBytes);
    nOutBytes = 0;
    assert(mdl->psrv[sid].fcnService != NULL);
    (*mdl->psrv[sid].fcnService)(mdl->psrv[sid].p1,pszIn,phi->nInBytes,
				 pszOut,&nOutBytes);
    assert(nOutBytes <= mdl->psrv[sid].nOutBytes);
    pho->idFrom = mdl->idSelf;
    pho->sid = sid;
    pho->nInBytes = phi->nInBytes;
    pho->nOutBytes = nOutBytes;
    MPI_Send(mdl->pszOut,nOutBytes+sizeof(SRVHEAD),
	     MPI_BYTE, id,MDL_TAG_RPL, comm);

    return sid;
    }



void mdlHandler(MDL mdl) {
    int sid;

    /*
    ** We must choose the correct communicator.  All nodes initially
    ** communicate only to their peers, except for the split masters;
    ** those wait for instructions from the main process group.
    */
    if ( mdlSelf(mdl) == 0 && mdlOldSelf(mdl) != 0 ) {
	assert( mdl->commCount == 2 );
	mdl->commMDL = mdl->commList[mdl->commCount-1];
	}
    else {
	mdl->commMDL = mdl->commList[0];
	}

    do {
	sid = mdlHandleOne(mdl);
	}
    while (sid != SRV_STOP);

    }

#define MDL_TAG_CACHECOM	10
#define MDL_MID_CACHEIN		1
#define MDL_MID_CACHEREQ	2
#define MDL_MID_CACHERPL	3
#define MDL_MID_CACHEOUT	4
#define MDL_MID_CACHEFLSH	5
#define MDL_MID_CACHEDONE	6

#define BILLION				1000000000

int mdlCacheReceive(MDL mdl,char *pLine) {
    CACHE *c;
    CAHEAD *ph = (CAHEAD *)mdl->pszRcv;
    char *pszRcv = &mdl->pszRcv[sizeof(CAHEAD)];
    CAHEAD *phRpl;
    char *pszRpl;
    char *t;
    int id, iTag;
    int s,n,i;
    MPI_Status status;
    int ret;
    int iLineSize;
#if 0
    char achDiag[256];
#endif

    ret = MPI_Wait(&mdl->ReqRcv, &status);
    assert(ret == MPI_SUCCESS);

#if 0
    /* Doesn't work...  MPI_SOURCE is -2 */
    if ( ph->id != status.MPI_SOURCE ) {
	printf( "MDL ERROR: id=%d, nThreads=%d, MPI_SOURCE=%d\n",
		ph->id, mdl->nThreads, status.MPI_SOURCE );
	printf("%d: cache %d, message %d, from %d, rec top\n",
	       mdl->idSelf, ph->cid, ph->mid, ph->id);
	assert( ph->id == status.MPI_SOURCE );
	}
#endif
#if 0
    sprintf(achDiag, "%d: cache %d, message %d, from %d, rec top\n",
	    mdl->idSelf, ph->cid, ph->mid, ph->id);
    mdlDiag(mdl, achDiag);
#endif

    c = &mdl->cache[ph->cid];
    assert(c->iType != MDL_NOCACHE);

    switch (ph->mid) {
    case MDL_MID_CACHEIN:
	++c->nCheckIn;
	ret = 0;
	break;
    case MDL_MID_CACHEOUT:
	++c->nCheckOut;
	ret = 0;
	break;
    case MDL_MID_CACHEREQ:
	/*
	 ** This is the tricky part! Here is where the real deadlock
	 ** difficulties surface. Making sure to have one buffer per
	 ** thread solves those problems here.
	 */
	pszRpl = &mdl->ppszRpl[ph->id][sizeof(CAHEAD)];
	phRpl = (CAHEAD *)mdl->ppszRpl[ph->id];
	phRpl->cid = ph->cid;
	phRpl->mid = MDL_MID_CACHERPL;
	phRpl->id = mdl->idSelf;

	s = ph->iLine*MDL_CACHELINE_ELTS;
	n = s + MDL_CACHELINE_ELTS;
	if ( n > c->nData ) n = c->nData;
	iLineSize = (n-s) * c->iDataSize;
	for(i=s; i<n; i++ ) {
	    t = (*c->getElt)(c->pData,i,c->iDataSize);
	    memcpy(pszRpl+(i-s)*c->iDataSize,t,c->iDataSize);
	    }
	if (mdl->pmidRpl[ph->id] != -1) {
	    MPI_Wait(&mdl->pReqRpl[ph->id], &status);
	    }
	mdl->pmidRpl[ph->id] = 0;
	MPI_Isend(phRpl,(int)sizeof(CAHEAD)+iLineSize,MPI_BYTE,
		  ph->id, MDL_TAG_CACHECOM, mdl->commMDL,
		  &mdl->pReqRpl[ph->id]);
	ret = 0;
	break;
    case MDL_MID_CACHEFLSH:
	assert(c->iType == MDL_COCACHE);
	s = ph->iLine*MDL_CACHELINE_ELTS;
	n = s + MDL_CACHELINE_ELTS;
	if (n > c->nData) n = c->nData;
	for(i=s;i<n;i++) {
		(*c->combine)(c->ctx,(*c->getElt)(c->pData,i,c->iDataSize),
			      &pszRcv[(i-s)*c->iDataSize]);
	    }
	ret = 0;
	break;
    case MDL_MID_CACHERPL:
	/*
	 ** For now assume no prefetching!
	 ** This means that this WILL be the reply to this Aquire
	 ** request.
	 */
	assert(pLine != NULL);
	iLineSize = c->iLineSize;
	for (i=0;i<iLineSize;++i) pLine[i] = pszRcv[i];
	if (c->iType == MDL_COCACHE && c->init) {
	    /*
	     ** Call the initializer function for all elements in
	     ** the cache line.
	     */
	    for (i=0;i<c->iLineSize;i+=c->iDataSize) {
		    (*c->init)(c->ctx,&pLine[i]);
		}
	    }
	ret = 1;
	break;
    case MDL_MID_CACHEDONE:
	/*
	 * No more caches, shouldn't get here.
	 */
	assert(0);
	break;
    default:
	assert(0);
	}

#if 0
    sprintf(achDiag, "%d: cache %d, message %d rec bottom\n", mdl->idSelf,
	    ph->cid, ph->mid);
    mdlDiag(mdl, achDiag);
#endif
    /*
     * Fire up next receive
     */
    id = MPI_ANY_SOURCE;
    iTag = MDL_TAG_CACHECOM;
    MPI_Irecv(mdl->pszRcv,mdl->iCaBufSize, MPI_BYTE, id,
	      iTag, mdl->commMDL, &mdl->ReqRcv);

    return ret;
    }


void AdjustDataSize(MDL mdl) {
    int i,iMaxDataSize;

    /*
     ** Change buffer size?
     */
    iMaxDataSize = 0;
    for (i=0;i<mdl->nMaxCacheIds;++i) {
	if (mdl->cache[i].iType == MDL_NOCACHE) continue;
	if (mdl->cache[i].iDataSize > iMaxDataSize) {
	    iMaxDataSize = mdl->cache[i].iDataSize;
	    }
	}
    if (iMaxDataSize != mdl->iMaxDataSize) {
	/*
	 ** Create new buffer with realloc?
	 ** Be very careful when reallocing buffers in other libraries
	 ** (not PVM) to be sure that the buffers are not in use!
	 ** A pending non-blocking receive on a buffer which is realloced
	 ** here will cause problems, make sure to take this into account!
	 ** This is certainly true in using the MPL library.
	 */
	MPI_Status status;
	CAHEAD caOut;

	/* cancel outstanding receive by sending a message to
	   myself */

	caOut.cid = 0;
	caOut.mid = MDL_MID_CACHEDONE;
	caOut.id = mdl->idSelf;
	MPI_Send(&caOut,sizeof(CAHEAD),MPI_BYTE, mdl->idSelf,
		 MDL_TAG_CACHECOM, mdl->commMDL);
	MPI_Wait(&mdl->ReqRcv, &status);

	mdl->iMaxDataSize = iMaxDataSize;
	mdl->iCaBufSize = (int)sizeof(CAHEAD) +
			  iMaxDataSize*(1 << MDL_CACHELINE_BITS);
	mdl->pszRcv = realloc(mdl->pszRcv,mdl->iCaBufSize);
	assert(mdl->pszRcv != NULL);
	for (i=0;i<mdl->nThreads;++i) {
	    mdl->ppszRpl[i] = realloc(mdl->ppszRpl[i],mdl->iCaBufSize);
	    assert(mdl->ppszRpl[i] != NULL);
	    }
	mdl->pszFlsh = realloc(mdl->pszFlsh,mdl->iCaBufSize);
	assert(mdl->pszFlsh != NULL);

	/*
	 * Fire up receive again.
	 */
	MPI_Irecv(mdl->pszRcv,mdl->iCaBufSize, MPI_BYTE,
		  MPI_ANY_SOURCE, MDL_TAG_CACHECOM,
		  mdl->commMDL, &mdl->ReqRcv);
	}
    }

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
    mdl->cacheSize = cacheSize;
    }

/*
 ** Common initialization for all types of caches.
 */
CACHE *CacheInitialize(MDL mdl,int cid,
		       void * (*getElt)(void *pData,int i,int iDataSize),
		       void *pData,int iDataSize,int nData) {
    CACHE *c;
    int i,nMaxCacheIds;
    int first;

    /*
     ** Allocate more cache spaces if required!
     */
    assert(cid >= 0);
    /*
     * first cache?
     */
    first = 1;
    for (i = 0; i < mdl->nMaxCacheIds; ++i) {
	if (mdl->cache[i].iType != MDL_NOCACHE) {
	    first = 0;
	    break;
	    }
	}
    if (first) {
	/*
	 * Fire up first receive
	 */
	MPI_Irecv(mdl->pszRcv,mdl->iCaBufSize, MPI_BYTE, MPI_ANY_SOURCE,
		  MDL_TAG_CACHECOM, mdl->commMDL, &mdl->ReqRcv);

#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
	mdl->dWaiting = mdl->dComputing = mdl->dSynchronizing = 0.0;
	mdl->nTicks = getticks();
#endif
	}

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
    c->iLineSize = MDL_CACHELINE_ELTS*c->iDataSize;

    c->iKeyShift = 0;
    while ((1 << c->iKeyShift) < mdl->nThreads) ++c->iKeyShift;
    c->iIdMask = (1 << c->iKeyShift) - 1;

    if (c->iKeyShift < MDL_CACHELINE_BITS) {
	/*
	 * Key will be (index & MDL_INDEX_MASK) | id.
	 */
	c->iInvKeyShift = MDL_CACHELINE_BITS;
	c->iKeyShift = 0;
	}
    else {
	/*
	 * Key will be (index & MDL_INDEX_MASK) << KeyShift | id.
	 */
	c->iInvKeyShift = c->iKeyShift;
	c->iKeyShift -= MDL_CACHELINE_BITS;
	}

    /*
     ** Determine the number of cache lines to be allocated.
     */
    assert( mdl->cacheSize >= (1<<MDL_CACHELINE_BITS)*10*c->iDataSize);
    c->nLines = (mdl->cacheSize/c->iDataSize) >> MDL_CACHELINE_BITS;
    c->iLine = 1;
    c->nTrans = 1;
    while (c->nTrans < c->nLines) c->nTrans *= 2;
    c->nTrans *= 2;
    c->iTransMask = c->nTrans-1;
    /*
     **	Set up the translation table.
     */
    c->pTrans = malloc(c->nTrans*sizeof(int));
    assert(c->pTrans != NULL);
    for (i=0;i<c->nTrans;++i) c->pTrans[i] = 0;
    /*
     ** Set up the tags. Note pTag[0] is a Sentinel!
     */
    c->pTag = malloc(c->nLines*sizeof(CTAG));
    assert(c->pTag != NULL);
    for (i=0;i<c->nLines;++i) {
	c->pTag[i].iKey = MDL_INVALID_KEY;
	c->pTag[i].nLock = 0;
	c->pTag[i].iLink = 0;
	}
    c->pTag[0].nLock = 1;     /* always locked */
    c->iLastVictim = 0;       /* This makes sure we have the first iVictim = 1 */
    c->nAccess = 0;
    c->nAccHigh = 0;
    c->nMiss = 0;				/* !!!, not NB */
    c->nColl = 0;				/* !!!, not NB */
    c->nMin = 0;				/* !!!, not NB */
    /*
     ** Allocate cache data lines.
     */
    if (mdl->freeCacheLines) {
	c->pLine = (char *)mdl->freeCacheLines;
	mdl->freeCacheLines = mdl->freeCacheLines->next;
	}
    else {
	c->pLine = malloc(mdl->cacheSize);
	}
    /*c->pLine = malloc(c->nLines*c->iLineSize);*/
    assert(c->pLine != NULL);
    c->nCheckOut = 0;

    /*
     ** Set up the request message as much as possible!
     */
    c->caReq.cid = cid;
    c->caReq.mid = MDL_MID_CACHEREQ;
    c->caReq.id = mdl->idSelf;
    return(c);
    }

/*
 ** Initialize a Read-Only caching space.
 */
void mdlROcache(MDL mdl,int cid,
		void * (*getElt)(void *pData,int i,int iDataSize),
		void *pData,int iDataSize,int nData) {
    CACHE *c;
    int id;
    CAHEAD caIn;
    char achDiag[256];

    c = CacheInitialize(mdl,cid,getElt,pData,iDataSize,nData);
    c->iType = MDL_ROCACHE;
    /*
     ** For an ROcache these two functions are not needed.
     */
    c->init = NULL;
    c->combine = NULL;
    sprintf(achDiag, "%d: before CI, cache %d\n", mdl->idSelf, cid);
    mdlDiag(mdl, achDiag);
    /*
     ** THIS IS A SYNCHRONIZE!!!
     */
    caIn.cid = cid;
    caIn.mid = MDL_MID_CACHEIN;
    caIn.id = mdl->idSelf;
    if (mdl->idSelf == 0) {
	c->nCheckIn = 1;
	while (c->nCheckIn < mdl->nThreads) {
	    mdlCacheReceive(mdl, NULL);
	    }
	}
    else {
	/*
	 ** Must use non-blocking sends here, we will never wait
	 ** for these sends to complete, but will know for sure
	 ** that they have completed.
	 */
	MPI_Send(&caIn,sizeof(CAHEAD),MPI_BYTE, 0,
		 MDL_TAG_CACHECOM, mdl->commMDL);
	}
    sprintf(achDiag, "%d: In CI, cache %d\n", mdl->idSelf, cid);
    mdlDiag(mdl, achDiag);
    if (mdl->idSelf == 0) {
	for (id = 1; id < mdl->nThreads; id++) {
	    MPI_Send(&caIn,sizeof(CAHEAD),MPI_BYTE, id,
		     MDL_TAG_CACHECOM, mdl->commMDL);
	    }
	}
    else {
	c->nCheckIn = 0;
	while (c->nCheckIn == 0) {
	    mdlCacheReceive(mdl,NULL);
	    }
	}
    sprintf(achDiag, "%d: After CI, cache %d\n", mdl->idSelf, cid);
    mdlDiag(mdl, achDiag);
    AdjustDataSize(mdl);
    MPI_Barrier(mdl->commMDL);
    }

/*
 ** Initialize a Combiner caching space.
 */
void mdlCOcache(MDL mdl,int cid,
		void * (*getElt)(void *pData,int i,int iDataSize),
		void *pData,int iDataSize,int nData,
		void *ctx,void (*init)(void *,void *),void (*combine)(void *,void *,void *)) {
    CACHE *c;
    int id;
    CAHEAD caIn;

    c = CacheInitialize(mdl,cid,getElt,pData,iDataSize,nData);
    c->iType = MDL_COCACHE;
    assert(init);
    c->init = init;
    assert(combine);
    c->combine = combine;
    c->ctx = ctx;
    /*
     ** THIS IS A SYNCHRONIZE!!!
     */
    caIn.cid = cid;
    caIn.mid = MDL_MID_CACHEIN;
    caIn.id = mdl->idSelf;
    if (mdl->idSelf == 0) {
	c->nCheckIn = 1;
	while (c->nCheckIn < mdl->nThreads) {
	    mdlCacheReceive(mdl, NULL);
	    }
	}
    else {
	/*
	 ** Must use non-blocking sends here, we will never wait
	 ** for these sends to complete, but will know for sure
	 ** that they have completed.
	 */
	MPI_Send(&caIn,sizeof(CAHEAD),MPI_BYTE, 0,
		 MDL_TAG_CACHECOM, mdl->commMDL);
	}
    if (mdl->idSelf == 0) {
	for (id = 1; id < mdl->nThreads; id++) {
	    MPI_Send(&caIn,sizeof(CAHEAD),MPI_BYTE, id,
		     MDL_TAG_CACHECOM, mdl->commMDL);
	    }
	}
    else {
	c->nCheckIn = 0;
	while (c->nCheckIn == 0) {
	    mdlCacheReceive(mdl,NULL);
	    }
	}
    AdjustDataSize(mdl);
    MPI_Barrier(mdl->commMDL);
    }

void mdlFinishCache(MDL mdl,int cid) {
    CACHE *c = &mdl->cache[cid];
    CAHEAD caOut;
    CAHEAD *caFlsh = (CAHEAD *)mdl->pszFlsh;
    char *pszFlsh = &mdl->pszFlsh[sizeof(CAHEAD)];
    mdlkey_t iKey;
    int i,id;
    char *t;
    int j;
    int last;
    MPI_Status status;
    MPI_Request reqFlsh;
    MPI_Request reqBoth[2];
    int index;

#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
	{
	ticks nTicks = getticks();
	mdl->dComputing += elapsed( nTicks, mdl->nTicks );
	mdl->nTicks = nTicks;
	}
#endif
    if (c->iType == MDL_COCACHE) {
	/*
	 * Extra checkout to let everybody finish before
	 * flushes start.
	* I think this makes for bad synchronizes --trq
	caOut.cid = cid;
	caOut.mid = MDL_MID_CACHEOUT;
	caOut.id = mdl->idSelf;
	for(id = 0; id < mdl->nThreads; id++) {
	    if(id == mdl->idSelf)
		continue;
	    MPI_Send(&caOut,sizeof(CAHEAD),MPI_BYTE, id,
		     MDL_TAG_CACHECOM, mdl->commMDL);
	    }
	++c->nCheckOut;
	while(c->nCheckOut < mdl->nThreads) {
	  mdlCacheReceive(mdl, NULL);
	    }
	c->nCheckOut = 0;
	 */
	/*
	 ** Must flush all valid data elements.
	 */
	caFlsh->cid = cid;
	caFlsh->mid = MDL_MID_CACHEFLSH;
	caFlsh->id = mdl->idSelf;
	for (i=1;i<c->nLines;++i) {
	    iKey = c->pTag[i].iKey;
	    if (iKey != MDL_INVALID_KEY) {
		/*
		 ** Flush element since it is valid!
		 */
		id = mdlkey_t_to_int(iKey & c->iIdMask);
		caFlsh->iLine = mdlkey_t_to_int(iKey >> c->iInvKeyShift);
		t = &c->pLine[i*c->iLineSize];
		for (j = 0; j < c->iLineSize; ++j)
		    pszFlsh[j] = t[j];
		/*
		 * Use Synchronous send so as not to
		 * overwhelm the receiver.
		 */
		MPI_Issend(caFlsh, (int)sizeof(CAHEAD)+c->iLineSize,
			   MPI_BYTE, id, MDL_TAG_CACHECOM,
			   mdl->commMDL, &reqFlsh);
		/*
		 * Wait for the Flush to complete, but
		 * also service any incoming cache requests.
		*/
		reqBoth[0] = mdl->ReqRcv;
		reqBoth[1] = reqFlsh;

		while (1) {
		    MPI_Waitany(2, reqBoth, &index, &status);
		    assert(!(index != 0 && reqBoth[0] ==
			     MPI_REQUEST_NULL));
		    mdl->ReqRcv = reqBoth[0];
		    if (index == 1) /* Flush has completed */
			break;
		    else if (index == 0) {
			mdlCacheReceive(mdl, NULL);
			reqBoth[0] = mdl->ReqRcv;
			}
		    else
			assert(0);
		    }
		}
	    }
	}
    /*
     ** THIS IS A SYNCHRONIZE!!!
     */
    caOut.cid = cid;
    caOut.mid = MDL_MID_CACHEOUT;
    caOut.id = mdl->idSelf;
    if (mdl->idSelf == 0) {
	++c->nCheckOut;
	while (c->nCheckOut < mdl->nThreads) {
	    mdlCacheReceive(mdl, NULL);
	    }
	}
    else {
	MPI_Send(&caOut,sizeof(CAHEAD),MPI_BYTE, 0,
		 MDL_TAG_CACHECOM, mdl->commMDL);
	}
    if (mdl->idSelf == 0) {
	for (id = 1; id < mdl->nThreads; id++) {
	    MPI_Send(&caOut,sizeof(CAHEAD),MPI_BYTE, id,
		     MDL_TAG_CACHECOM, mdl->commMDL);
	    }
	}
    else {
	c->nCheckOut = 0;
	while (c->nCheckOut == 0) {
	    mdlCacheReceive(mdl,NULL);
	    }
	}
    /*
     ** Free up storage and finish.
     */
    free(c->pTrans);
    free(c->pTag);
    {
        FREECACHELINES *fcl = (FREECACHELINES *)c->pLine;
	fcl->next = mdl->freeCacheLines;
	mdl->freeCacheLines = fcl;
    }
    /*free(c->pLine);*/
    c->iType = MDL_NOCACHE;

    AdjustDataSize(mdl);
    /*
     * last cache?
     */
    last = 1;
    for (i = 0; i < mdl->nMaxCacheIds; ++i) {
	if (mdl->cache[i].iType != MDL_NOCACHE) {
	    last = 0;
	    break;
	    }
	}
    /*
     * shut down CacheReceive.
     * Note: I'm sending a message to myself.
     */
    if (last) {
	MPI_Status status;

	caOut.cid = cid;
	caOut.mid = MDL_MID_CACHEDONE;
	caOut.id = mdl->idSelf;
	MPI_Send(&caOut,sizeof(CAHEAD),MPI_BYTE, mdl->idSelf,
		 MDL_TAG_CACHECOM, mdl->commMDL);
	MPI_Wait(&mdl->ReqRcv, &status);
	}
    MPI_Barrier(mdl->commMDL);
#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
	{
	ticks nTicks = getticks();
	mdl->dSynchronizing += elapsed( nTicks, mdl->nTicks );
	mdl->nTicks = nTicks;
	}
#endif
    }


void mdlCacheBarrier(MDL mdl,int cid) {
    CACHE *c = &mdl->cache[cid];
    CAHEAD caOut;
    int id;

#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
	{
	ticks nTicks = getticks();
	mdl->dComputing += elapsed( nTicks, mdl->nTicks );
	mdl->nTicks = nTicks;
	}
#endif

    /*
    ** THIS IS A SYNCHRONIZE!!!
    */
    caOut.cid = cid;
    caOut.mid = MDL_MID_CACHEOUT;
    caOut.id = mdl->idSelf;
    if (mdl->idSelf == 0) {
	++c->nCheckOut;
	while (c->nCheckOut < mdl->nThreads) {
	    mdlCacheReceive(mdl, NULL);
	    }
	}
    else {
	MPI_Send(&caOut,sizeof(CAHEAD),MPI_BYTE, 0,
		 MDL_TAG_CACHECOM, mdl->commMDL);
	}
    if (mdl->idSelf == 0) {
	for (id = 1; id < mdl->nThreads; id++) {
	    MPI_Send(&caOut,sizeof(CAHEAD),MPI_BYTE, id,
		     MDL_TAG_CACHECOM, mdl->commMDL);
	    }
	}
    else {
	c->nCheckOut = 0;
	while (c->nCheckOut == 0) {
	    mdlCacheReceive(mdl,NULL);
	    }
	}
    c->nCheckOut = 0;
    MPI_Barrier(mdl->commMDL);
#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
	{
	ticks nTicks = getticks();
	mdl->dSynchronizing += elapsed( nTicks, mdl->nTicks );
	mdl->nTicks = nTicks;
	}
#endif
    }


void mdlCacheCheck(MDL mdl) {
    int flag;
    MPI_Status status;

    while (1) {
	MPI_Test(&mdl->ReqRcv, &flag, &status);
	if (flag == 0)
	    break;

	mdlCacheReceive(mdl,NULL);
	}
    }


void *mdlDoMiss(MDL mdl, int cid, int iIndex, int id, mdlkey_t iKey, int lock) {
    CACHE *c = &mdl->cache[cid];
    char *pLine;
    int iElt,i,iLine;
    mdlkey_t iKeyVic;
    int idVic;
    int iVictim,*pi;
    char ach[80];
    CAHEAD *caFlsh;
    char *pszFlsh;
    MPI_Status status;
    MPI_Request reqFlsh;


#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
	{
	ticks nTicks = getticks();
	mdl->dComputing += elapsed( nTicks, mdl->nTicks );
	mdl->nTicks = nTicks;
	}
#endif

    /*
    ** Cache Miss.
    */
    iLine = iIndex >> MDL_CACHELINE_BITS;
    c->caReq.cid = cid;
    c->caReq.mid = MDL_MID_CACHEREQ;
    c->caReq.id = mdl->idSelf;
    c->caReq.iLine = iLine;
    MPI_Send(&c->caReq,sizeof(CAHEAD),MPI_BYTE,
	     id,MDL_TAG_CACHECOM, mdl->commMDL);
    ++c->nMiss;

    if (++c->iLastVictim == c->nLines) c->iLastVictim = 1;
    iVictim = c->iLastVictim;

    iElt = iIndex & MDL_CACHE_MASK;
    for (i=1;i<c->nLines;++i) {
	if (!c->pTag[iVictim].nLock) {
	    /*
	     ** Found victim.
	     */
	    iKeyVic = c->pTag[iVictim].iKey;
	    /*
	     ** 'pLine' will point to the actual data line in the cache.
	     */
	    pLine = &c->pLine[iVictim*c->iLineSize];
	    caFlsh = NULL;
	    if (iKeyVic != MDL_INVALID_KEY) {
		if (c->iType == MDL_COCACHE) {
		    /*
		    ** Flush element since it is valid!
		    */
		    idVic = mdlkey_t_to_int(iKeyVic&c->iIdMask);
		    caFlsh = (CAHEAD *)mdl->pszFlsh;
		    pszFlsh = &mdl->pszFlsh[sizeof(CAHEAD)];
		    caFlsh->cid = cid;
		    caFlsh->mid = MDL_MID_CACHEFLSH;
		    caFlsh->id = mdl->idSelf;
		    caFlsh->iLine = mdlkey_t_to_int(iKeyVic >> c->iInvKeyShift);
		    for (i = 0; i < c->iLineSize; ++i)
			pszFlsh[i] = pLine[i];
		    MPI_Isend(caFlsh, (int)sizeof(CAHEAD)+c->iLineSize,
			      MPI_BYTE, idVic,
			      MDL_TAG_CACHECOM, mdl->commMDL, &reqFlsh);
		    }
		/*
		 ** If valid iLine then "unlink" it from the cache.
		 */
		pi = &c->pTrans[iKeyVic & c->iTransMask];
		while (*pi != iVictim) pi = &c->pTag[*pi].iLink;
		*pi = c->pTag[iVictim].iLink;
		}
	    c->pTag[iVictim].iKey = iKey;
	    if (lock) c->pTag[iVictim].nLock = 1;
	    /*
	     **	Add the modified victim tag back into the cache.
	     ** Note: the new element is placed at the head of the chain.
	     */
	    pi = &c->pTrans[iKey & c->iTransMask];
	    c->pTag[iVictim].iLink = *pi;
	    *pi = iVictim;
	    goto Await;
	    }
	if (++iVictim == c->nLines) iVictim = 1;
	}
    /*
     ** Cache Failure!
     */
    sprintf(ach,"MDL CACHE FAILURE: cid == %d, no unlocked lines!\n",cid);
    mdlDiag(mdl,ach);
    exit(1);
Await:
    c->iLastVictim = iVictim;
    /*
     ** At this point 'pLine' is the recipient cache line for the
     ** data requested from processor 'id'.
     */
    while (1) {
	if (mdlCacheReceive(mdl,pLine)) {
	    if (caFlsh)
		MPI_Wait(&reqFlsh, &status);
#if defined(INSTRUMENT) && defined(HAVE_TICK_COUNTER)
		{
		ticks nTicks = getticks();
		mdl->dWaiting += elapsed( nTicks, mdl->nTicks );
		mdl->nTicks = nTicks;
		}
#endif
	    return(&pLine[iElt*c->iDataSize]);
	    }
	}
    }


void mdlRelease(MDL mdl,int cid,void *p) {
    CACHE *c = &mdl->cache[cid];
    int iLine;

    iLine = ((char *)p - c->pLine) / c->iLineSize;
    /*
     ** Check if the pointer fell in a cache line, otherwise it
     ** must have been a local pointer.
     */
    if (iLine > 0 && iLine < c->nLines) {
	--c->pTag[iLine].nLock;
	assert(c->pTag[iLine].nLock >= 0);
	}
#ifdef OLD_CACHE
    else {
	iData = ((char *)p - c->pData) / c->iDataSize;
	assert(iData < c->nData);
	}
#endif
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
    mdl->dWaiting = mdl->dComputing = mdl->dSynchronizing = 0.0;
    mdl->nTicks = getticks();
    }

static double TimeFraction(MDL mdl) {
    double dTotal = mdl->dComputing + mdl->dWaiting + mdl->dSynchronizing;
    if ( dTotal <= 0.0 ) return 0.0;
    return 100.0 / dTotal;
    }

double mdlTimeComputing(MDL mdl) {
    return mdl->dComputing * TimeFraction(mdl);
    }

double mdlTimeSynchronizing(MDL mdl) {
    return mdl->dSynchronizing * TimeFraction(mdl);
    }

double mdlTimeWaiting(MDL mdl) {
    return mdl->dWaiting * TimeFraction(mdl);
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
    assert(n1<=a1);
    *pgrid = grid = malloc(sizeof(struct mdlGridContext)); assert(grid!=NULL);
    grid->n1 = n1;
    grid->n2 = n2;
    grid->n3 = n3;
    grid->a1 = a1;

    /* This will be shared later (see mdlGridShare) */
    grid->id = malloc(sizeof(*grid->id)*(grid->n3));    assert(grid->id!=NULL);
    grid->rs = mdlMalloc(mdl,sizeof(*grid->rs)*mdl->nThreads); assert(grid->rs!=NULL);
    grid->rn = mdlMalloc(mdl,sizeof(*grid->rn)*mdl->nThreads); assert(grid->rn!=NULL);

    /* The following need to be set to appropriate values still. */
    grid->s = grid->n = grid->nlocal = 0;
    }

void mdlGridFinish(MDL mdl, MDLGRID grid) {
    if (grid->rs) free(grid->rs);
    if (grid->rn) free(grid->rn);
    if (grid->id) free(grid->id);
    free(grid);
    }

void mdlGridSetLocal(MDL mdl,MDLGRID grid,int s, int n, int nlocal) {
    assert( s>=0 && s<grid->n3);
    assert( n>=0 && s+n<=grid->n3);
    grid->s = s;
    grid->n = n;
    grid->nlocal = nlocal;
    }

/*
** Share the local GRID information with other processors by,
**   - finding the starting slab and number of slabs on each processor
**   - building a mapping from slab to processor id.
*/
void mdlGridShare(MDL mdl,MDLGRID grid) {
    int i, id;

    MPI_Allgather(&grid->s,sizeof(*grid->rs),MPI_BYTE,
		  grid->rs,sizeof(*grid->rs),MPI_BYTE,
		  mdl->commMDL);
    MPI_Allgather(&grid->n,sizeof(*grid->rn),MPI_BYTE,
		  grid->rn,sizeof(*grid->rn),MPI_BYTE,
		  mdl->commMDL);
    /* Calculate on which processor each slab can be found. */
    for(id=0; id<mdl->nThreads; id++ ) {
	for( i=grid->rs[id]; i<grid->rs[id]+grid->rn[id]; i++ ) grid->id[i] = id;
	}
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
    int sz,nz,sy,ny,nlocal;

    *pfft = NULL;
    fft = malloc(sizeof(struct mdlFFTContext));
    assert(fft != NULL);

    /* Contruct the FFTW plans */
    fft->fplan = rfftw3d_mpi_create_plan(mdl->commMDL,
					 n3, n2, n1,
					 FFTW_REAL_TO_COMPLEX,
					 (bMeasure ? FFTW_MEASURE : FFTW_ESTIMATE) );
    fft->iplan = rfftw3d_mpi_create_plan(mdl->commMDL,
					 /* dim.'s of REAL data --> */ n3, n2, n1,
					 FFTW_COMPLEX_TO_REAL,
					 (bMeasure ? FFTW_MEASURE : FFTW_ESTIMATE));
    rfftwnd_mpi_local_sizes( fft->fplan, &nz, &sz, &ny, &sy,&nlocal);

    /*
    ** Dimensions of k-space and r-space grid.  Note transposed order.
    ** Note also that the "actual" dimension 1 side of the r-space array
    ** can be (and usually is) larger than "n1" because of the inplace FFT.
    */
    mdlGridInitialize(mdl,&fft->rgrid,n1,n2,n3,2*(n1/2+1));
    mdlGridInitialize(mdl,&fft->kgrid,n1/2+1,n3,n2,n1/2+1);

    mdlGridSetLocal(mdl,fft->rgrid,sz,nz,nlocal);
    mdlGridSetLocal(mdl,fft->kgrid,sy,ny,nlocal/2);
    mdlGridShare(mdl,fft->rgrid);
    mdlGridShare(mdl,fft->kgrid);

    *pfft = fft;
    return nlocal;
    }

void mdlFFTFinish( MDL mdl, MDLFFT fft ) {
    rfftwnd_mpi_destroy_plan(fft->fplan);
    rfftwnd_mpi_destroy_plan(fft->iplan);
    mdlGridFinish(mdl,fft->kgrid);
    mdlGridFinish(mdl,fft->rgrid);
    free(fft);
    }

fftw_real *mdlFFTMAlloc( MDL mdl, MDLFFT fft ) {
    return mdlGridMalloc(mdl,fft->rgrid,sizeof(fftw_real));
    }

void mdlFFTFree( MDL mdl, MDLFFT fft, void *p ) {
    mdlGridFree(mdl,fft->rgrid,p);
    }

void mdlFFT( MDL mdl, MDLFFT fft, fftw_real *data, int bInverse ) {
    rfftwnd_mpi_plan plan = bInverse ? fft->iplan : fft->fplan;
    rfftwnd_mpi(plan,1,data,0,FFTW_TRANSPOSED_ORDER);
    }
#endif
