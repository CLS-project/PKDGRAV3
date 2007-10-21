#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#define USE_IO_TIPSY
#ifdef USE_IO_TIPSY
#include <rpc/xdr.h>
#endif

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

#include "iohdf5.h"
#include "pst.h"
#include "io.h"

#define CHUNKSIZE (32*1024)
#define EPSILON (-1e20)

static void makeName( IO io, char *achOutName, const char *inName, int iIndex )
{
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



static void ioSave(IO io, const char *filename, total_t N,
		   double dTime,
		   double dEcosmo, double dTimeOld, double dUOld,
		   int bDouble )
{
    hid_t fileID;
    IOHDF5 iohdf5;
    IOHDF5V ioDen;
    IOHDF5V ioPot;
    local_t i;
#ifdef USE_IO_TIPSY
    char *tname = malloc(strlen(filename)+5);
    FILE *tfp;
    XDR xdr;
    int j;
    float fTmp;
#endif

#ifdef USE_IO_TIPSY
    strcpy(tname,filename);
    strcat(tname,".std");
    tfp = fopen(tname,"wb");
    assert( tfp != NULL );
    free(tname);
    xdrstdio_create(&xdr,tfp,XDR_ENCODE);
    if ( mdlSelf(io->mdl) == 0 ) {
	uint32_t nBodies = N;
	uint32_t nDims = 3;
	uint32_t nZero = 0;
	assert(xdr_double(&xdr,&dTime));
	assert(xdr_u_int(&xdr,&nBodies));
	assert(xdr_u_int(&xdr,&nDims));
	assert(xdr_u_int(&xdr,&nZero));
	assert(xdr_u_int(&xdr,&nBodies));
	assert(xdr_u_int(&xdr,&nZero));
	assert(xdr_u_int(&xdr,&nZero));
    }
#endif

    /* Create the output file */
    fileID = ioCreate(filename);
    iohdf5 = ioHDF5Initialize( fileID, CHUNKSIZE, bDouble );
    ioDen  = ioHDFF5NewVector( iohdf5, "density",  IOHDF5_SINGLE );
    ioPot  = ioHDFF5NewVector( iohdf5, "potential",IOHDF5_SINGLE );

    ioHDF5WriteAttribute( iohdf5, "dTime",   H5T_NATIVE_DOUBLE, &dTime );
    ioHDF5WriteAttribute( iohdf5, "dEcosmo", H5T_NATIVE_DOUBLE, &dEcosmo );
    ioHDF5WriteAttribute( iohdf5, "dTimeOld",H5T_NATIVE_DOUBLE, &dTimeOld );
    ioHDF5WriteAttribute( iohdf5, "dUOld",   H5T_NATIVE_DOUBLE, &dUOld );

    for( i=0; i<io->N; i++ ) {
/*	ioHDF5AddDark(iohdf5, io->iMinOrder+i,
		      io->r[i].v, io->v[i].v,
		      io->m[i], io->s[i], io->p[i] );*/
	ioHDF5AddDark(iohdf5, io->iMinOrder+i,
		      io->r[i].v, io->v[i].v,
		      io->ioClasses[io->vClass[i]].dMass,
		      io->ioClasses[io->vClass[i]].dSoft, io->p[i] );
	ioHDF5AddVector( ioDen, io->iMinOrder+i, io->d[i] );
	ioHDF5AddVector( ioPot, io->iMinOrder+1, io->p[i] );

#ifdef USE_IO_TIPSY
	fTmp = io->ioClasses[io->vClass[i]].dMass;
	xdr_float(&xdr,&fTmp);
	for(j=0;j<3;j++) {
	    if (bDouble)
		xdr_double(&xdr,&io->r[i].v[j]);
	    else {
		fTmp = io->r[i].v[j];
		xdr_float(&xdr,&fTmp);
	    }
	}
	for(j=0;j<3;j++) {
	    fTmp = io->v[i].v[j];
	    xdr_float(&xdr,&fTmp);
	}
	fTmp = io->ioClasses[io->vClass[i]].dSoft;
	xdr_float(&xdr,&fTmp);
	fTmp = io->p[i];
	xdr_float(&xdr,&fTmp);
#endif
    }
    ioHDF5Finish(iohdf5);

    H5assert(H5Fflush(fileID,H5F_SCOPE_GLOBAL));
    H5assert(H5Fclose(fileID));

#ifdef USE_IO_TIPSY
    xdr_destroy(&xdr);
    fclose(tfp);
#endif
}


static void ioLoad(IO io, const char *filename,
		   total_t iIndex, total_t N,
		   double *dTime, double *dEcosmo,
		   double *dTimeOld, double *dUOld )
{
    hid_t fileID;
    IOHDF5 iohdf5;
    IOHDF5V ioDen;
    IOHDF5V ioPot;
    PINDEX iOrder = 0;
    local_t iOffset;
    local_t i, iClass;

    /* Open the output file */
    fileID = ioOpen(filename);
    iohdf5 = ioHDF5Initialize( fileID, CHUNKSIZE, 0 );

    ioHDF5ReadAttribute( iohdf5, "dTime",   H5T_NATIVE_DOUBLE, dTime );
    ioHDF5ReadAttribute( iohdf5, "dEcosmo", H5T_NATIVE_DOUBLE, dEcosmo );
    ioHDF5ReadAttribute( iohdf5, "dTimeOld",H5T_NATIVE_DOUBLE, dTimeOld );
    ioHDF5ReadAttribute( iohdf5, "dUOld",   H5T_NATIVE_DOUBLE, dUOld );

    if ( N == 0 ) N = ioHDF5DarkCount( iohdf5 ) - iIndex;
    iOffset = io->N;
    io->N += N;

    assert( io->N <= io->nAllocated );

    if ( iIndex ) ioHDF5SeekDark( iohdf5, iIndex );

    /*printf( "%s: %lu -> %lu\n", filename, iIndex, iIndex+N );*/

    for( i=0; i<N; i++ ) {
	local_t iLocal = iOffset + i;
	FLOAT dMass, dSoft, dPot;

	ioHDF5GetDark(iohdf5, &iOrder,
		   io->r[iLocal].v, io->v[iLocal].v,
		   &dMass, &dSoft, &dPot );

	/*FIXME: linear search - MAX 256, <10 typical */
	for( iClass=0; iClass<io->nClasses; iClass++ )
	    if ( io->ioClasses[iClass].dMass == dMass && io->ioClasses[iClass].dSoft == dSoft )
		break;
	if ( iClass != io->nClasses ) {
	    io->ioClasses[iClass].iMaxOrder = iOrder;
	}
	else {
	    assert( iClass<MAX_IO_CLASSES);
	    io->nClasses++;
	    io->ioClasses[iClass].dMass = dMass;
	    io->ioClasses[iClass].dSoft = dSoft;
	    io->ioClasses[iClass].iMinOrder = iOrder;
	    io->ioClasses[iClass].iMaxOrder = iOrder;
	}
	io->vClass[iLocal] = iClass;
    }

    ioHDF5Finish(iohdf5);

    H5assert(H5Fclose(fileID));
}

void ioInitialize(IO *pio,MDL mdl)
{
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

void ioAddServices(IO io,MDL mdl)
{
    mdlAddService(mdl,IO_SETUP,io,
		  (void (*)(void *,void *,int,void *,int *)) ioSetup,
		  sizeof(struct inIOSetup),0);
    mdlAddService(mdl,IO_START_SAVE,io,
		  (void (*)(void *,void *,int,void *,int *)) ioStartSave,
		  sizeof(struct inStartSave),0);
    mdlAddService(mdl,IO_START_LOAD,io,
		  (void (*)(void *,void *,int,void *,int *)) ioStartLoad,
		  sizeof(struct inStartLoad),0);
    mdlAddService(mdl,IO_PLAN_LOAD,io,
		  (void (*)(void *,void *,int,void *,int *)) ioPlanLoad,
		  sizeof(struct inPlanLoad),sizeof(struct outPlanLoad));

    mdlAddService(mdl,IO_ALLOCATE,io,
		  (void (*)(void *,void *,int,void *,int *)) ioAllocate,
		  sizeof(struct inIOAllocate),0);
    mdlAddService(mdl,IO_START_RECV,io,
		  (void (*)(void *,void *,int,void *,int *)) ioStartRecv,
		  sizeof(struct inStartRecv),0);
    mdlAddService(mdl,IO_START_SEND,io,
		  (void (*)(void *,void *,int,void *,int *)) ioStartSend,
		  sizeof(struct inStartSend),0);
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
void ioStartSave(IO io,void *vin,int nIn,void *vout,int *pnOut)
{
    int id;
    struct inStartSave *save = vin;
    struct inStartRecv recv;
    total_t iCount;
#ifdef USE_PNG
	struct inMakePNG png;
#endif

    mdlassert(io->mdl,sizeof(struct inStartSave)==nIn);
    mdlassert(io->mdl,mdlSelf(io->mdl)==0);

    mdlSetComm(io->mdl,0); /* Talk to our peers */
    recv.dTime   = save->dTime;
    recv.dEcosmo = save->dEcosmo;
    recv.dTimeOld= save->dTimeOld;
    recv.dUOld   = save->dUOld;
    strcpy(recv.achOutName,save->achOutName);
    recv.bCheckpoint = save->bCheckpoint;
    recv.N       = save->N;
    iCount = save->N / mdlIO(io->mdl);

    printf( "Starting to save %"PRIu64" particles (~%"PRIu64" per I/O node)\n",
	    save->N, iCount );

    for( id=1; id<mdlIO(io->mdl); id++ ) {
	recv.iIndex = iCount * id;
	recv.nCount = iCount;
	if ( id+1 == mdlIO(io->mdl) )
	    recv.nCount = save->N - recv.iIndex;
	mdlReqService(io->mdl,id,IO_START_RECV,&recv,sizeof(recv));
    }

    recv.iIndex = 0;
    recv.nCount = iCount;
    ioStartRecv(io,&recv,sizeof(recv),NULL,0);

    for( id=1; id<mdlIO(io->mdl); id++ ) {
	mdlGetReply(io->mdl,id,NULL,NULL);
    }

#ifdef USE_PNG
    png.iResolution = 1024;
    png.minValue = -1;
    png.maxValue = 5;
    png.scale = 1;
    strcpy(png.achOutName,save->achOutName);
    for( id=1; id<mdlIO(io->mdl); id++ ) {
	mdlReqService(io->mdl,id,IO_MAKE_PNG,&png,sizeof(png));
    }
    ioMakePNG(io,&png,sizeof(png),NULL,0);

    for( id=1; id<mdlIO(io->mdl); id++ ) {
	mdlGetReply(io->mdl,id,NULL,NULL);
    }
#endif
    mdlSetComm(io->mdl,1);
}

/*
** Here we figure out how many input file there are.
*/
void ioPlanLoad(IO io,void *vin,int nIn,void *vout,int *pnOut)
{
    int id,i;
    int me;
    H5E_auto_t oldAutoFunc;
    void *     oldAutoData;
    char achInName[256];
    hid_t fileID;
    IOHDF5 iohdf5;

    struct inPlanLoad  *in  = vin;
    struct outPlanLoad *out = vout;

    me = mdlSelf(io->mdl);

    if ( me == 0 ) {
	mdlSetComm(io->mdl,0); /* Talk to our peers */
	for( id=1; id<mdlIO(io->mdl); id++ ) {
	    mdlReqService(io->mdl,id,IO_PLAN_LOAD,vin,nIn);
	}
	mdlSetComm(io->mdl,1);
    }

    for( i=0; i<MDL_MAX_IO_PROCS; i++ )
	out->nCount[i] = 0;

    /* Check only files for which I am responsible */
    H5Eget_auto(&oldAutoFunc, &oldAutoData );
    H5Eset_auto(0,0);
    for( i=me; i<MDL_MAX_IO_PROCS; i += mdlIO(io->mdl) ) {
	makeName( io, achInName, in->achInName, i );
	if ( H5Fis_hdf5(achInName) > 0 ) {
	    fileID = ioOpen(achInName);
	    iohdf5 = ioHDF5Initialize( fileID, CHUNKSIZE, 0 );
	    out->nCount[i] = ioHDF5DarkCount(iohdf5);
	    assert(ioHDF5ReadAttribute(
		       iohdf5, "dTime", H5T_NATIVE_DOUBLE, &out->dExpansion ));
	    ioHDF5Finish(iohdf5);
	    H5assert(H5Fclose(fileID));
	}
	else break;
    }
    H5Eset_auto(oldAutoFunc, oldAutoData );

    if ( me == 0 ) {
	struct outPlanLoad outChild;
	int outN;
	mdlSetComm(io->mdl,0); /* Talk to our peers */
	for( id=1; id<mdlIO(io->mdl); id++ ) {
	    mdlGetReply(io->mdl,id,&outChild,&outN);
	    assert(outN == sizeof(outChild) );

	    for( i=0; i<MDL_MAX_IO_PROCS; i++ )
		out->nCount[i] += outChild.nCount[i];
	}
	mdlSetComm(io->mdl,1);
    }
    if ( pnOut ) *pnOut = sizeof(struct outPlanLoad);
}


/*
**  This service is called from the peer root node (processor 0) and initiates
**  a load process.
*/
void ioStartLoad(IO io,void *vin,int nIn,void *vout,int *pnOut)
{
    int id, i;
    struct inStartLoad *load = vin;
    struct inStartSend send, send0;
    struct inIOAllocate alloc;
    total_t N, nLocal, n;

    mdlassert(io->mdl,sizeof(struct inStartLoad)==nIn);
    mdlassert(io->mdl,mdlSelf(io->mdl)==0);

    N = 0;
    for( i=0; i<load->nFiles; i++ )
	N += load->nCount[i];
    nLocal = N / mdlIO(io->mdl);

    mdlSetComm(io->mdl,0); /* Talk to our peers */
    for( id=1; id<mdlIO(io->mdl); id++ ) {
	if ( id+1 == mdlIO(io->mdl) )
	    alloc.nCount = N - nLocal*(mdlIO(io->mdl)-1);
	else
	    alloc.nCount = nLocal;
	mdlReqService(io->mdl,id,IO_ALLOCATE,&alloc,sizeof(alloc));
    }
    mdlSetComm(io->mdl,1);

    alloc.nCount = nLocal;
    ioAllocate(io,&alloc, sizeof(alloc),0,0);

    mdlSetComm(io->mdl,0); /* Talk to our peers */
    for( id=1; id<mdlIO(io->mdl); id++ ) {
	mdlGetReply(io->mdl,id,NULL,NULL);
    }
    mdlSetComm(io->mdl,1);


    mdlSetComm(io->mdl,0); /* Talk to our peers */

    strcpy( send.achInName, load->achInName );
    send.N = N;
    send.iFirstFile = send.iLastFile = 0;
    send.iFirstOffset = send.iLastOffset = 0;

    for( id=0; id<mdlIO(io->mdl); id++ ) {
	send.iMinOrder = id * nLocal;
	if ( id+1 == mdlIO(io->mdl) )
	    nLocal = N - nLocal*(mdlIO(io->mdl)-1);
	send.iMaxOrder = send.iMinOrder + nLocal;
	send.iFirstFile   = send.iLastFile;
	send.iFirstOffset = send.iLastOffset;
	if ( send.iFirstOffset == load->nCount[send.iFirstFile] ) {
	    send.iFirstOffset = 0;
	    send.iFirstFile++;
	}

	send.iLastFile    = send.iFirstFile;
	send.iLastOffset  = send.iFirstOffset + nLocal;

	if ( send.iLastOffset > load->nCount[send.iLastFile] ) {
	    n = nLocal - (load->nCount[send.iLastFile] - send.iFirstOffset);
	    send.iLastFile = send.iFirstFile + 1;

	    while ( n > load->nCount[send.iLastFile] ) {
		n -= load->nCount[send.iLastFile];
		send.iLastFile++;
	    }
	    send.iLastOffset = n;
	}
	/*printf( "%d %d.%lu -> %d.%lu\n",
		id, send.iFirstFile, send.iFirstOffset,
		send.iLastFile, send.iLastOffset );*/
	if ( id == 0 ) {
	    strcpy( send0.achInName, send.achInName );
	    send0.N = N;
	    send0.iMinOrder = send.iMinOrder;
	    send0.iMaxOrder = send.iMaxOrder;
	    send0.iFirstFile = send.iFirstFile;
	    send0.iFirstOffset = send.iFirstOffset;
	    send0.iLastFile = send.iLastFile;
	    send0.iLastOffset = send.iLastOffset;
	}
	else {
	    mdlReqService(io->mdl,id,IO_START_SEND,&send,sizeof(send));
	}
    }
    mdlSetComm(io->mdl,1);

    ioStartSend(io,&send0,sizeof(send0),NULL,0);

    mdlSetComm(io->mdl,0); /* Talk to our peers */
    for( id=1; id<mdlIO(io->mdl); id++ ) {
	mdlGetReply(io->mdl,id,NULL,NULL);
    }
    mdlSetComm(io->mdl,1);

}

static int ioUnpackIO(void *vctx, int *id, size_t nSize, void *vBuff)
{
    IO io = vctx;
    PIO *pio = vBuff;
    uint_fast32_t nIO = nSize / sizeof(PIO);
    uint_fast32_t i, k, d;

    mdlassert(io->mdl,nIO*sizeof(PIO) == nSize);
    mdlassert(io->mdl,nIO<=io->nExpected);

    for( i=0; i<nIO; i++ ) {
	total_t iOrder = pio[i].iOrder;
	size_t iLocal = iOrder - io->iMinOrder;
	size_t iClass;

	mdlassert(io->mdl,iOrder>=io->iMinOrder);
	mdlassert(io->mdl,iOrder<io->iMaxOrder);
	mdlassert(io->mdl,iLocal<io->N);
	mdlassert(io->mdl,pio[i].fMass > 0.0);
	mdlassert(io->mdl,pio[i].fSoft > 0.0);

	for( d=0; d<3; d++ ) {
	    io->r[iLocal].v[d] = pio[i].r[d];
	    io->v[iLocal].v[d] = pio[i].v[d];
	}

	/*FIXME: linear search - MAX 256, <10 typical */
	for( iClass=0; iClass<io->nClasses; iClass++ )
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

static int ioPackIO(void *vctx, int *id, size_t nSize, void *vBuff)
{
    IO io = vctx;
    PIO *pio = vBuff;
    uint_fast32_t nIO = nSize / sizeof(PIO);
    uint_fast32_t i, d;
    total_t n;
    local_t nPerThread;
    total_t iOffset, maxOrder;

    /* Calculate how many particles to send, and to whom they should be sent */
    nPerThread = io->nTotal / mdlThreads(io->mdl);
    *id = io->iOrder / nPerThread;
    maxOrder = (*id+1)==mdlThreads(io->mdl) ? io->iMaxOrder : (*id+1)*nPerThread;
    if ( maxOrder > io->iMaxOrder ) maxOrder = io->iMaxOrder;
    n = maxOrder - io->iOrder;
    if ( n > nIO ) n = nIO;

    /* Okay, send "n" particles to "id" -- n can be zero meaning we are done */
    for( i=0; i<n; i++ ) {
	int I = io->iOrder - io->iMinOrder;
	assert( I < io->N );

	pio->iOrder = io->iOrder;
	for( d=0; d<3; d++ ) {
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



void ioAllocate(IO io,void *vin,int nIn,void *vout,int *pnOut)
{
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

void ioSetup(IO io,void *vin,int nIn,void *vout,int *pnOut)
{
    struct inIOSetup *setup = vin;
    struct inIOAllocate alloc;
    total_t iCount;
    int id;

    mdlassert(io->mdl,sizeof(struct inIOSetup)==nIn);

    iCount = setup->N / mdlIO(io->mdl);
    alloc.nCount = iCount;
    ioAllocate(io,&alloc, sizeof(alloc),0,0);

    mdlSetComm(io->mdl,0); /* Talk to our peers */
    for( id=1; id<mdlIO(io->mdl); id++ ) {
	alloc.nCount = iCount;
	if ( id+1 == mdlIO(io->mdl) )
	    alloc.nCount = setup->N - id*iCount;
	mdlReqService(io->mdl,id,IO_ALLOCATE,&alloc,sizeof(alloc));
    }
    for( id=1; id<mdlIO(io->mdl); id++ ) {
	mdlGetReply(io->mdl,id,NULL,NULL);
    }
    mdlSetComm(io->mdl,1);
}

/*
**  Here we actually wait for the data from the Work nodes
*/
void ioStartRecv(IO io,void *vin,int nIn,void *vout,int *pnOut)
{
    struct inStartRecv *recv = vin;
    struct inIOAllocate alloc;
    char achOutName[256];
    total_t iOrder;

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
	   recv->dTimeOld, recv->dUOld,
	   recv->bCheckpoint ? IOHDF5_DOUBLE : IOHDF5_SINGLE );
}


void ioStartSend(IO io,void *vin,int nIn,void *vout,int *pnOut)
{
    struct inStartSend *send = vin;
    char achInName[256];
    double dTime, dEcosmo, dTimeOld, dUOld;
    int i;
    local_t N, n, O;
    total_t F, L;

    io->iMinOrder = send->iMinOrder;
    io->iMaxOrder = send->iMaxOrder;
    io->nTotal = send->N;
    io->N = 0;
    for( i=send->iFirstFile; i<=send->iLastFile; i++ ) {
	O = i==send->iFirstFile ? send->iFirstOffset : 0;
	n = i==send->iLastFile ? send->iLastOffset-O : 0;
	makeName( io, achInName, send->achInName, i );
	ioLoad(io, achInName, O, n, &dTime, &dEcosmo, &dTimeOld, &dUOld );
    }
    assert( io->N == io->iMaxOrder - io->iMinOrder );

    n = send->N / mdlThreads(io->mdl);

    io->iOrder = io->iMinOrder;
    mdlSetComm(io->mdl,1); /* Talk to the work process */
    mdlSend(io->mdl,-1,ioPackIO,io);
    mdlSetComm(io->mdl,0);

}

#ifdef USE_PNG
void ioMakePNG(IO io,void *vin,int nIn,void *vout,int *pnOut)
{
    struct inMakePNG *make = vin;
    char achOutName[256];
    float *limg, *img, X, Y;
    uint_fast32_t R, N, i;
    int x, y;
    PNG png;
    FILE *fp;

    R = make->iResolution;
    N = R * R;

    limg = malloc( N*sizeof(float) );
    assert( limg != NULL );

    for( i=0; i<N; i++ ) limg[i] = EPSILON;

    // Project the density onto a grid and find the maximum (ala Tipsy)
    for( i=0; i<io->N; i++ ) {
	X = (io->r[i].v[0]+0.5) * make->scale;
	Y = (io->r[i].v[1]+0.5) * make->scale;
	if ( X < 0 || X >= 1 || Y < 0 || Y >= 1 ) continue;

	x = X * R;
	y = Y * R;
	assert( x>=0 && x<R && y>=0 && y<R );
	if ( io->d[i] > limg[x+R*y] )
	    limg[x+R*y] = io->d[i];
    }

    if ( mdlSelf(io->mdl) == 0 ) {
	img = malloc( N*sizeof(float) );
	assert( img != NULL );
	mdlReduce(io->mdl,limg,img,N,MPI_FLOAT,MPI_MAX,0);

	makeName( io, achOutName, make->achOutName, mdlSelf(io->mdl) );
	sprintf(achOutName+strlen(achOutName),"-%.1f.png", make->scale );
	fp = fopen( achOutName, "wb" );
	if ( fp != NULL ) {
	    png = pngInitialize( make->iResolution, make->minValue, make->maxValue );
	    pngWrite( png, fp, img );
	    pngFinish(png);
	    fclose(fp);
	}
	free(img);
    }
    else {
	mdlReduce(io->mdl,limg,0,N,MPI_FLOAT,MPI_MAX,0);
    }

    free(limg);
}
#endif
