#include <malloc.h>
#include <alloca.h>
#include <hdf5.h>
#include "pst.h"
#include "io.h"

#define H5assert(rc) assert( rc >= 0 )
#define CHUNKSIZE (16*1024)

static void writeAttribute( hid_t groupID, const char *name,
                            hid_t dataType, void *data ) {
    hid_t dataSpace, attrID;
    herr_t rc;
    hsize_t dims = 1;

    dataSpace = H5Screate_simple( 1, &dims, NULL ); H5assert(dataSpace);
    attrID = H5Acreate( groupID, name, dataType, dataSpace, H5P_DEFAULT );
    H5assert(attrID);
    rc = H5Awrite( attrID, dataType, data ); H5assert(rc);
    H5assert(H5Aclose( attrID ));
    H5assert(H5Sclose(dataSpace));
}

static hid_t newSet(hid_t fileID, const char *group, const char *name, 
                    uint64_t count, int nDims, hid_t dataType )
{
    hid_t dataProperties, dataSpace, dataSet;
    hsize_t iDims[2];
    char *nameBuffer;

    /* Construct the name of this object */
    if ( group && group[0] ) {
        nameBuffer = alloca( strlen(group) + strlen(name) + 4 );
        strcpy( nameBuffer, "/" );
        strcat( nameBuffer, group );
        strcat( nameBuffer, "/" );
    }
    else {
        nameBuffer = alloca( strlen(name) + 2 );
        strcpy( nameBuffer, "/" );
    }
    strcat( nameBuffer, name );

    /* Create a dataset property so we can set the chunk size */
    dataProperties = H5Pcreate( H5P_DATASET_CREATE );
    H5assert( dataProperties );
    iDims[0] = CHUNKSIZE;
    iDims[1] = 1;
    H5assert( H5Pset_chunk( dataProperties, nDims>1?2:1, iDims ));

    /* And the dataspace */
    iDims[0] = count;
    iDims[1] = nDims;
    dataSpace = H5Screate_simple( nDims>1?2:1, iDims, iDims );
    H5assert( dataSpace );

    /* Create the data set */
    dataSet = H5Dcreate( fileID, nameBuffer, 
                         dataType, dataSpace, dataProperties );
    H5assert( dataSet );
    H5assert( H5Pclose( dataProperties ) );
    H5assert( H5Sclose( dataSpace ) );

    return dataSet;
}


static void writeSet(
    hid_t fileID,           /* HDF5 into which to add this set */
    hid_t fileType,         /* Data type in the file */
    const char *groupName,  /* Name of the group ("dark", etc.) */
    const char *setName,    /* Name of the dataset ("pos", "vel", etc.) */
    const void *pBuffer,    /* Buffer containing the data */
    hid_t memType,          /* Data type in memory */
    hsize_t nRows,          /* Number of rows to write */
    hsize_t nCols )         /* Number of columns (normally 1 or 3) */
{
    hid_t dataSet;
    herr_t rc;

    dataSet = newSet( fileID, groupName, setName, nRows, nCols, fileType );
    rc = H5Dwrite( dataSet, memType, H5S_ALL, H5S_ALL, H5P_DEFAULT, pBuffer );
    H5Dclose( dataSet );
}


static hid_t newGroup( hid_t fileID, const char *groupName ) {
    hid_t groupID;
    char *nameBuffer;

    nameBuffer = alloca( strlen(groupName)+2 );
    strcpy( nameBuffer, "/" );
    strcat( nameBuffer, groupName );
    groupID = H5Gcreate( fileID, nameBuffer, 0 ); H5assert(groupID);
    return groupID;
}


static void writeType(
    IO io,
    double dTime,
    int bSingle,
    hid_t fileID,
    int typeID,
    const char *typeName)
{
    hid_t groupID;

    /* This is the native type of FLOAT values - normally double */
    hid_t memType = sizeof(FLOAT)==sizeof(float) ? H5T_NATIVE_FLOAT : H5T_NATIVE_DOUBLE;
    hid_t fileType = bSingle ? H5T_NATIVE_FLOAT : memType;

    groupID = newGroup( fileID, typeName ); H5assert( groupID );

    /* Save the current simulation time (expansion factor) */
    writeAttribute( groupID, "time", H5T_NATIVE_DOUBLE, &dTime );

    writeSet( fileID, fileType, typeName, "pos",  io->r, memType, io->N, 3 );
    writeSet( fileID, fileType, typeName, "vel",  io->v, memType, io->N, 3 );
    writeSet( fileID, fileType, typeName, "mass", io->m, memType, io->N, 1 );
}

static void ioSave(IO io, const char *filename, double dTime, int bSingle)
{
    hid_t fileID;

    /* Create the output file */
    fileID = H5Fcreate( filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
    H5assert(fileID);

    writeType( io, dTime, bSingle, fileID, 0, "dark" );
}

void ioInitialize(IO *pio,MDL mdl)
{
    IO io;
    io = (IO)malloc(sizeof(struct ioContext));
    mdlassert(mdl,io != NULL);
    io->mdl = mdl;

    io->N = 0;
    io->r = NULL;
    io->v = NULL;
    io->m = NULL;

    *pio = io;
}

void ioAddServices(IO io,MDL mdl)
{
    mdlAddService(mdl,IO_START_SAVE,io,
		  (void (*)(void *,void *,int,void *,int *)) ioStartSave,
		  sizeof(struct inStartSave),0);
    mdlAddService(mdl,IO_START_RECV,io,
		  (void (*)(void *,void *,int,void *,int *)) ioStartRecv,
		  sizeof(struct inStartRecv),0);
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

    mdlassert(io->mdl,sizeof(struct inStartSave)==nIn);
    mdlassert(io->mdl,mdlSelf(io->mdl)==0);

    mdlSetComm(io->mdl,0); /* Talk to our peers */
    recv.dTime = save->dTime;
    for( id=1; id<mdlIO(io->mdl); id++ ) {
	recv.nCount = save->nCount[id];
	mdlReqService(io->mdl,id,IO_START_RECV,&recv,sizeof(recv));
    }

    recv.nCount = save->nCount[0];
    ioStartRecv(io,&recv,sizeof(recv),NULL,0);

    for( id=1; id<mdlIO(io->mdl); id++ ) {
	mdlGetReply(io->mdl,id,NULL,NULL);
    }
    mdlSetComm(io->mdl,1);

}


static int ioUnpackIO(void *vctx, int nSize, void *vBuff)
{
    IO io = vctx;
    PIO *pio = vBuff;
    int nIO = nSize / sizeof(PIO);
    int i, d;

    mdlassert(io->mdl,nIO<=io->nExpected);

    for( i=0; i<nIO; i++ ) {
	for( d=0; d<3; d++ ) {
	    io->r[io->nReceived].v[d] = pio[i].r[d];
	    io->v[io->nReceived].v[d] = pio[i].v[d];
	}
	io->m[io->nReceived] = pio[i].fMass;
	io->nReceived++;
    }

    io->nExpected -= nIO;
    return io->nExpected;
}

/*
**  Here we actually wait for the data from the Work nodes
*/
void ioStartRecv(IO io,void *vin,int nIn,void *vout,int *pnOut)
{
    struct inStartRecv *recv = vin;
    char testFilename[200];

    mdlassert(io->mdl,sizeof(struct inStartRecv)==nIn);
    io->nExpected = recv->nCount;
    io->nReceived = 0;

    if ( io->nExpected > io->N ) {
	if ( io->N ) {
	    free(io->m);
	    free(io->v);
	    free(io->r);
	}
	io->N = io->nExpected;
	io->r = malloc(io->N*sizeof(ioV3));
	io->v = malloc(io->N*sizeof(ioV3));
	io->m = malloc(io->N*sizeof(FLOAT));
    }

    mdlSetComm(io->mdl,1); /* Talk to the work process */
    mdlRecv(io->mdl,MPI_ANY_SOURCE,ioUnpackIO,io);
    mdlSetComm(io->mdl,0);


    sprintf(testFilename, "testout.%02d.%05d.h5", mdlSelf(io->mdl), 12 );
    ioSave(io, testFilename, recv->dTime, 0 );

}
