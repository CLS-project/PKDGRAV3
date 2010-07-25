#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <rpc/xdr.h>

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#ifdef HAVE_ALLOCA_H
#include <alloca.h>
#endif
#ifdef USE_PNG
#include "png.h"
#endif
#include <math.h>
#include <inttypes.h>

#include "fio.h"
#include "pst.h"
#include "io.h"

const char *io_c_module_id = "$Id$";
const char *io_h_module_id = IO_H_MODULE_ID;

#define CHUNKSIZE (32*1024)
#define MINVALUE (-1e20)

static void makeName( IO io, char *achOutName, const char *inName, int iIndex ) {
    char *p;

    strcpy( achOutName, inName );
    p = strstr( achOutName, "&I" );
    if ( p ) {
	int n = p - achOutName;
	sprintf( p, "%03d", iIndex );
	strcat( p, inName + n + 2 );
	}
    else {
	p = achOutName + strlen(achOutName);
	sprintf(p,".%03d", iIndex);
	}
    }


#ifdef USE_HDF5
/* Create an HDF5 file for output */
hid_t ioCreate( const char *filename ) {
    hid_t fileID;

    /* Create the output file */
    fileID = H5Fcreate( filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
    H5assert(fileID);

    return fileID;
    }

hid_t ioOpen( const char *filename ) {
    hid_t fileID;

    /* Create the output file */
    fileID = H5Fopen( filename, H5F_ACC_RDONLY, H5P_DEFAULT );
    H5assert(fileID);

    return fileID;
    }
#endif

static void writeToFIO(IO io, FIO fio) {
    local_t i;

    for ( i=0; i<io->N; i++ ) {
	fioWriteDark(fio,io->iMinOrder+i,
		     io->r[i].v, io->v[i].v,
		     io->ioClasses[io->vClass[i]].dMass,
		     io->ioClasses[io->vClass[i]].dSoft,
		     io->p[i] );
	}
    }

static void ioSave(IO io, const char *filename, total_t N,
		   double dTime,
		   double dEcosmo, double dTimeOld, double dUOld,
		   int bDouble, int iStandard, int bHDF5 ) {
    FIO fio;
    char *tname;
    local_t i;
    int nSph, nDark, nStar;

    nSph = nStar = 0;
    nDark = N;

/* Make sure we save something */
#ifndef USE_HDF5
    if ( iStandard < 0 ) iStandard = 1;
#else
    if ( !bHDF5 && iStandard < 0 ) iStandard = 1;
#endif
    if ( iStandard >= 0 ) {
	tname = malloc(strlen(filename)+5);
	strcpy(tname,filename);
	if ( iStandard > 0 ) strcat(tname,".std");
	else strcat(tname,".bin");
	fio = fioTipsyCreatePart(tname,0,bDouble,iStandard, dTime, 
				 nSph, nDark, nStar, io->iMinOrder);
	free(tname);
	writeToFIO(io,fio);
	fioClose(fio);
	}

#ifdef USE_HDF5
    if ( bHDF5 ) {
	fio = fioHDF5Create(filename, (bDouble?FIO_FLAG_CHECKPOINT:0) | FIO_FLAG_POTENTIAL);
	fioSetAttr( fio, "dTime",   FIO_TYPE_DOUBLE, &dTime);
	fioSetAttr( fio, "dEcosmo", FIO_TYPE_DOUBLE, &dEcosmo );
	fioSetAttr( fio, "dTimeOld",FIO_TYPE_DOUBLE, &dTimeOld );
	fioSetAttr( fio, "dUOld",   FIO_TYPE_DOUBLE, &dUOld );
	writeToFIO(io,fio);
	fioClose(fio);
	}
#endif
    }

void ioInitialize(IO *pio,MDL mdl) {
    IO io;
    io = (IO)malloc(sizeof(struct ioContext));
    mdlassert(mdl,io != NULL);
    io->mdl = mdl;

    io->nAllocated = 0;
    io->N = 0;
    io->r = NULL;
    io->v = NULL;
    io->d = NULL;
    io->p = NULL;

    io->nClasses = 0;

    *pio = io;
    }

void ioAddServices(IO io,MDL mdl) {
    mdlAddService(mdl,IO_SETUP,io,
		  (void (*)(void *,void *,int,void *,int *)) ioSetup,
		  sizeof(struct inIOSetup),0);
    mdlAddService(mdl,IO_START_SAVE,io,
		  (void (*)(void *,void *,int,void *,int *)) ioStartSave,
		  sizeof(struct inStartSave),0);

    mdlAddService(mdl,IO_ALLOCATE,io,
		  (void (*)(void *,void *,int,void *,int *)) ioAllocate,
		  sizeof(struct inIOAllocate),0);
    mdlAddService(mdl,IO_START_RECV,io,
		  (void (*)(void *,void *,int,void *,int *)) ioStartRecv,
		  sizeof(struct inStartRecv),0);
#ifdef USE_PNG
    mdlAddService(mdl,IO_MAKE_PNG,io,
		  (void (*)(void *,void *,int,void *,int *)) ioMakePNG,
		  sizeof(struct inMakePNG),0);
#endif
    }

/*
**  This service is called from the peer root node (processor 0) and initiates
**  a save process.
*/
void ioStartSave(IO io,void *vin,int nIn,void *vout,int *pnOut) {
    int id;
    float scale;
    struct inStartSave *save = vin;
    struct inStartRecv recv;
    total_t iCount;
#ifdef USE_PNG
    struct inMakePNG png;
#endif

    mdlassert(io->mdl,sizeof(struct inStartSave)==nIn);
    mdlassert(io->mdl,mdlSelf(io->mdl)==0);

    mdlSetComm(io->mdl,0); /* Talk to our peers */
    recv.dTime     = save->dTime;
    recv.dEcosmo   = save->dEcosmo;
    recv.dTimeOld  = save->dTimeOld;
    recv.dUOld     = save->dUOld;
    strcpy(recv.achOutName,save->achOutName);
    recv.bCheckpoint = save->bCheckpoint;
    recv.N         = save->N;
    recv.iStandard = save->iStandard;
    recv.bHDF5     = save->bHDF5;
    iCount = save->N / mdlIO(io->mdl);

    printf( "Starting to save %"PRIu64" particles (~%"PRIu64" per I/O node)\n",
	    save->N, iCount );

    for ( id=1; id<mdlIO(io->mdl); id++ ) {
	recv.iIndex = iCount * id;
	recv.nCount = iCount;
	if ( id+1 == mdlIO(io->mdl) )
	    recv.nCount = save->N - recv.iIndex;
	mdlReqService(io->mdl,id,IO_START_RECV,&recv,sizeof(recv));
	}

    recv.iIndex = 0;
    recv.nCount = iCount;
    ioStartRecv(io,&recv,sizeof(recv),NULL,0);

    for ( id=1; id<mdlIO(io->mdl); id++ ) {
	mdlGetReply(io->mdl,id,NULL,NULL);
	}

#ifdef USE_PNG
    for ( scale=1.0; scale<=8.0+1e-5; scale*=2.0 ) {
	png.iResolution = 10240;
	png.minValue = -1;
	png.maxValue = 5;
	png.scale = scale;
	strcpy(png.achOutName,save->achOutName);
	for ( id=1; id<mdlIO(io->mdl); id++ ) {
	    mdlReqService(io->mdl,id,IO_MAKE_PNG,&png,sizeof(png));
	    }
	ioMakePNG(io,&png,sizeof(png),NULL,0);

	for ( id=1; id<mdlIO(io->mdl); id++ ) {
	    mdlGetReply(io->mdl,id,NULL,NULL);
	    }
	}
#endif
    mdlSetComm(io->mdl,1);
    }

static int ioUnpackIO(void *vctx, int *id, size_t nSize, void *vBuff) {
    IO io = vctx;
    PIO *pio = vBuff;
    uint_fast32_t nIO = nSize / sizeof(PIO);
    uint_fast32_t i, d;

    mdlassert(io->mdl,nIO*sizeof(PIO) == nSize);
    mdlassert(io->mdl,nIO<=io->nExpected);

    for ( i=0; i<nIO; i++ ) {
	total_t iOrder = pio[i].iOrder;
	size_t iLocal = iOrder - io->iMinOrder;
	size_t iClass;

	mdlassert(io->mdl,iOrder>=io->iMinOrder);
	mdlassert(io->mdl,iOrder<io->iMaxOrder);
	mdlassert(io->mdl,iLocal<io->N);
	mdlassert(io->mdl,pio[i].fMass > 0.0);
	mdlassert(io->mdl,pio[i].fSoft > 0.0);

	for ( d=0; d<3; d++ ) {
	    io->r[iLocal].v[d] = pio[i].r[d];
	    io->v[iLocal].v[d] = pio[i].v[d];
	    }

	/*FIXME: linear search - MAX 256, <10 typical */
	for ( iClass=0; iClass<io->nClasses; iClass++ )
	    if ( io->ioClasses[iClass].dMass == pio[i].fMass && io->ioClasses[iClass].dSoft == pio[i].fSoft )
		break;
	if ( iClass != io->nClasses ) {
	    if ( io->ioClasses[iClass].iMaxOrder < iOrder )
		io->ioClasses[iClass].iMinOrder = iOrder;
	    else if ( io->ioClasses[iClass].iMinOrder > iOrder )
		io->ioClasses[iClass].iMinOrder = iOrder;
	    }
	else {
	    assert( iClass<MAX_IO_CLASSES);
	    io->nClasses++;
	    io->ioClasses[iClass].dMass = pio[i].fMass;
	    io->ioClasses[iClass].dSoft = pio[i].fSoft;
	    io->ioClasses[iClass].iMinOrder = iOrder;
	    io->ioClasses[iClass].iMaxOrder = iOrder;
	    }
	io->vClass[iLocal] = iClass;
	io->d[iLocal] = pio[i].fDensity;
	io->p[iLocal] = pio[i].fPot;

	}

    io->nExpected -= nIO;
    return io->nExpected;
    }

static int ioPackIO(void *vctx, int *id, size_t nSize, void *vBuff) {
    IO io = vctx;
    PIO *pio = vBuff;
    uint_fast32_t nIO = nSize / sizeof(PIO);
    uint_fast32_t i, d;
    total_t n;
    local_t nPerThread;
    total_t maxOrder;

    /* Calculate how many particles to send, and to whom they should be sent */
    nPerThread = io->nTotal / mdlThreads(io->mdl);
    *id = io->iOrder / nPerThread;
    maxOrder = (*id+1)==mdlThreads(io->mdl) ? io->iMaxOrder : (*id+1)*nPerThread;
    if ( maxOrder > io->iMaxOrder ) maxOrder = io->iMaxOrder;
    n = maxOrder - io->iOrder;
    if ( n > nIO ) n = nIO;

    /* Okay, send "n" particles to "id" -- n can be zero meaning we are done */
    for ( i=0; i<n; i++ ) {
	int I = io->iOrder - io->iMinOrder;
	assert( I < io->N );

	pio->iOrder = io->iOrder;
	for ( d=0; d<3; d++ ) {
	    pio->r[d] = io->r[I].v[d];
	    pio->v[d] = io->v[I].v[d];
	    }
	pio->fDensity = 0.0;
	pio->fPot = 0.0;
	pio->fMass = io->ioClasses[io->vClass[I]].dMass;
	pio->fSoft = io->ioClasses[io->vClass[I]].dSoft;

	io->iOrder++;
	pio++;
	}
    assert( n>0 || io->iOrder == io->iMaxOrder );

    return n * sizeof(PIO);
    }



void ioAllocate(IO io,void *vin,int nIn,void *vout,int *pnOut) {
    struct inIOAllocate *alloc = vin;

    mdlassert(io->mdl,sizeof(struct inIOAllocate)==nIn);

    if ( alloc->nCount > io->nAllocated ) {
	if ( io->nAllocated ) {
	    free(io->vClass);
	    free(io->p);
	    free(io->d);
	    free(io->v);
	    free(io->r);
	    }
	io->nAllocated = alloc->nCount + 100; /* Room to grow... */
	io->r = malloc(io->nAllocated*sizeof(ioV3));  assert(io->r != NULL );
	io->v = malloc(io->nAllocated*sizeof(ioV3));  assert(io->v != NULL );
	io->d = malloc(io->nAllocated*sizeof(float));  assert(io->d != NULL );
	io->p = malloc(io->nAllocated*sizeof(float));  assert(io->p != NULL );
	io->vClass = malloc(io->nAllocated*sizeof(uint8_t)); assert( io->vClass != NULL );
	}
    }

void ioSetup(IO io,void *vin,int nIn,void *vout,int *pnOut) {
    struct inIOSetup *setup = vin;
    struct inIOAllocate alloc;
    total_t iCount;
    int id;

    mdlassert(io->mdl,sizeof(struct inIOSetup)==nIn);

    iCount = setup->N / mdlIO(io->mdl);
    alloc.nCount = iCount;
    ioAllocate(io,&alloc, sizeof(alloc),0,0);

    mdlSetComm(io->mdl,0); /* Talk to our peers */
    for ( id=1; id<mdlIO(io->mdl); id++ ) {
	alloc.nCount = iCount;
	if ( id+1 == mdlIO(io->mdl) )
	    alloc.nCount = setup->N - id*iCount;
	mdlReqService(io->mdl,id,IO_ALLOCATE,&alloc,sizeof(alloc));
	}
    for ( id=1; id<mdlIO(io->mdl); id++ ) {
	mdlGetReply(io->mdl,id,NULL,NULL);
	}
    mdlSetComm(io->mdl,1);
    }

/*
**  Here we actually wait for the data from the Work nodes
*/
void ioStartRecv(IO io,void *vin,int nIn,void *vout,int *pnOut) {
    struct inStartRecv *recv = vin;
    struct inIOAllocate alloc;
    char achOutName[256];

    mdlassert(io->mdl,sizeof(struct inStartRecv)==nIn);

    alloc.nCount = recv->nCount;
    ioAllocate(io,&alloc, sizeof(alloc),0,0);

    io->N = io->nExpected = recv->nCount;
    io->iMinOrder = recv->iIndex;
    io->iMaxOrder = recv->iIndex + recv->nCount;

    mdlSetComm(io->mdl,1); /* Talk to the work process */
    mdlRecv(io->mdl,-1,ioUnpackIO,io);
    mdlSetComm(io->mdl,0);

    makeName( io, achOutName, recv->achOutName, mdlSelf(io->mdl) );

    ioSave(io, achOutName, recv->N, recv->dTime, recv->dEcosmo,
	   recv->dTimeOld, recv->dUOld, recv->bCheckpoint,
	   recv->iStandard, recv->bHDF5 );
    }


#ifdef USE_PNG
void ioMakePNG(IO io,void *vin,int nIn,void *vout,int *pnOut) {
    struct inMakePNG *make = vin;
    char achOutName[256];
    float *limg, *slice, X, Y;
    uint_fast32_t R, N, i;
    uint_fast32_t Ns, Is, Ss;
    int x, y;
    PNG png;
    FILE *fp;

    R = make->iResolution;
    N = R * R;
    Ns = 1024*1024/sizeof(float);

    limg = malloc( N*sizeof(float) );
    assert( limg != NULL );

    for ( i=0; i<N; i++ ) limg[i] = MINVALUE;

    // Project the density onto a grid and find the maximum (ala Tipsy)
    for ( i=0; i<io->N; i++ ) {
	// Scale, crop and adjust range to [0,1)
	X = io->r[i].v[0] * make->scale + 0.5;
	Y = io->r[i].v[1] * make->scale + 0.5;
	if ( X < 0 || X >= 1 || Y < 0 || Y >= 1 ) continue;

	x = floor(X * R);
	y = floor(Y * R);
	assert( x>=0 && x<R && y>=0 && y<R );
	if ( io->d[i] > limg[x+R*y] )
	    limg[x+R*y] = io->d[i];
	}

    if ( mdlSelf(io->mdl) == 0 ) {
	slice = malloc( Ns*sizeof(float) );
	assert( slice != NULL );
	for ( Is=0; Is<N; Is+=Ss ) {
	    Ss = (N-Is) > Ns ? Ns : (N-Is);
	    mdlReduce(io->mdl,limg+Is,slice,Ss,MDL_FLOAT,MDL_MAX,0);
	    memcpy(limg+Is,slice,Ss*sizeof(float));
	    }
	makeName( io, achOutName, make->achOutName, mdlSelf(io->mdl) );
	sprintf(achOutName+strlen(achOutName),"-%.1f.png", make->scale );
	fp = fopen( achOutName, "wb" );
	if ( fp != NULL ) {
	    png = pngInitialize( make->iResolution, make->minValue, make->maxValue );
	    pngWrite( png, fp, limg );
	    pngFinish(png);
	    fclose(fp);
	    }
	free(slice);
	}
    else {
	for ( Is=0; Is<N; Is+=Ss ) {
	    Ss = (N-Is) > Ns ? Ns : (N-Is);
	    mdlReduce(io->mdl,limg+Is,0,Ss,MDL_FLOAT,MDL_MAX,0);
	    }
	}

    free(limg);
    }
#endif
