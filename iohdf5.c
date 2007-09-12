/******************************************************************************
 *  iohdf5.c - Read and Write particle data in HDF5 format.
 *
 *  - Call H5Fcreate or H5Fopen to create or open an HDF5 file
 *  - Call ioHDF5Initialize
 *
 *  Write:
 *    - Call ioHDF5AddDark (or Gas or Star) for each particle
 *  Read:
 *    - Call ioHDF5GetDark (or Gas or Star) for each particle
 *
 *  - Call ioHDF5Finish.
 *  - Call H5Fclose
 *
 *****************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <assert.h>
#include <stdlib.h>

#include "iohdf5.h"

#define GROUP_PARAMETERS "parameters"
#define GROUP_DARK       "dark"
#define GROUP_GAS        "gas"
#define GROUP_STAR       "star"

#define FIELD_POSITION   "position"
#define FIELD_VELOCITY   "velocity"
#define FIELD_ORDER      "order"
#define FIELD_CLASS      "class"
#define FIELD_CLASSES    "classes"

#define ATTR_IORDER      "iOrder"

/* Create a data group: dark, gas, star, etc. */
static hid_t CreateGroup( hid_t fileID, const char *groupName ) {
    hid_t groupID;

    groupID = H5Gcreate( fileID, groupName, 0 ); H5assert(groupID);

    return groupID;
}

static hid_t openSet(hid_t fileID, const char *name )
{
    hid_t dataSet;
    dataSet = H5Dopen( fileID, name );
    return dataSet;
}

static hsize_t getSetSize(hid_t setID) {
    hid_t spaceID;
    hsize_t dims[2], maxs[2];

    spaceID = H5Dget_space(setID); H5assert(spaceID);
    assert( H5Sis_simple(spaceID) > 0 );
    assert( H5Sget_simple_extent_ndims(spaceID) <= 3 );
    H5Sget_simple_extent_dims(spaceID,dims,maxs);
    H5Sclose(spaceID);
    return dims[0];
}

/* Create a dataset inside a group (positions, velocities, etc. */
static hid_t newSet(hid_t locID, const char *name, uint64_t chunk,
		    uint64_t count, int nDims, hid_t dataType )
{
    hid_t dataProperties, dataSpace, dataSet;
    hsize_t iDims[2], iMax[2];

    /* Create a dataset property so we can set the chunk size */
    dataProperties = H5Pcreate( H5P_DATASET_CREATE );
    H5assert( dataProperties );
    iDims[0] = chunk;
    iDims[1] = 1;
    H5assert( H5Pset_chunk( dataProperties, nDims>1?2:1, iDims ));

    /* Also request the FLETCHER checksum */
    H5assert( H5Pset_filter(dataProperties,H5Z_FILTER_FLETCHER32,0,0,NULL) );

    /* And the dataspace */
    iDims[0] = count;
    iDims[1] = nDims;
    iMax[0] = H5S_UNLIMITED;
    iMax[1] = nDims;

    dataSpace = H5Screate_simple( nDims>1?2:1, iDims, iMax );
    H5assert( dataSpace );

    /* Create the data set */
    dataSet = H5Dcreate( locID, name, 
                         dataType, dataSpace, dataProperties );
    H5assert( dataSet );
    H5assert( H5Pclose( dataProperties ) );
    H5assert( H5Sclose( dataSpace ) );

    return dataSet;
}

/* Read part of a set from disk into memory */
static void readSet(
    hid_t set_id,           /* set from which to read the data */
    void *pBuffer,          /* Buffer for the data */
    hid_t memType,          /* Data type in memory */
    hsize_t iOffset,        /* Offset into the set */
    hsize_t nRows,          /* Number of rows to read */
    hsize_t nCols )         /* Number of columns (normally 1 or 3) */
{
    hid_t memSpace, diskSpace;
    hsize_t dims[2], start[2];
    dims[0] = nRows;
    dims[1] = nCols;
    memSpace = H5Screate_simple(nCols>1?2:1,dims,0);
    diskSpace = H5Dget_space(set_id);
    start[0] = iOffset;
    start[1] = 0;

    H5Sselect_hyperslab(diskSpace,H5S_SELECT_SET,start,0,dims,0);
    H5Dread(set_id,memType,memSpace,diskSpace,H5P_DEFAULT,pBuffer);
    H5Sclose(memSpace);
    H5Sclose(diskSpace);
}


/* Add an attribute to a group */
static void writeAttribute( hid_t groupID, const char *name,
			    hid_t dataType, void *data ) {
    hid_t dataSpace, attrID;
    hsize_t dims = 1;

    dataSpace = H5Screate_simple( 1, &dims, NULL ); H5assert(dataSpace);
    attrID = H5Acreate( groupID,name,dataType,dataSpace,H5P_DEFAULT );
    H5assert(attrID);
    H5assert(H5Awrite( attrID, dataType, data ));
    H5assert(H5Aclose( attrID ));
    H5assert(H5Sclose(dataSpace));
}

/* Read an attribute from a group */
static int readAttribute( hid_t groupID, const char *name,
			  hid_t dataType, void *data ) {
    hid_t attrID;

    attrID = H5Aopen_name( groupID,name );
    if ( attrID == H5I_INVALID_HID ) return 0;
    H5assert(H5Aread( attrID, dataType, data ));
    H5assert(H5Aclose( attrID ));
    return 1;
}



/* Write part of a set from memory */
static void writeSet(
    hid_t set_id,           /* set into which to write the data */
    const void *pBuffer,    /* Buffer containing the data */
    hid_t memType,          /* Data type in memory */
    hsize_t iOffset,        /* Offset into the set */
    hsize_t nRows,          /* Number of rows to write */
    hsize_t nCols )         /* Number of columns (normally 1 or 3) */
{
    hid_t memSpace, diskSpace;
    hsize_t dims[2], start[2], size[2];

    dims[0] = nRows;
    dims[1] = nCols;
    memSpace = H5Screate_simple(nCols>1?2:1,dims,0);
    size[0] = iOffset + nRows;
    size[1] = nCols;
    H5Dextend(set_id,size);
    diskSpace = H5Dget_space(set_id);
    start[0] = iOffset;
    start[1] = 0;
    H5Sselect_hyperslab(diskSpace,H5S_SELECT_SET,start,0,dims,0);
    H5Dwrite(set_id,memType,memSpace,diskSpace,H5P_DEFAULT,pBuffer);
    H5Sclose(memSpace);
    H5Sclose(diskSpace);
}

static void flushVector( IOHDF5V iov )
{
    if ( iov->nBuffered ) {

	if ( iov->set_id == H5I_INVALID_HID ) {
	    iov->set_id = newSet(
		iov->io->darkBase.group_id, iov->name,
		iov->io->iChunkSize, 0, 1, iov->diskFloat );
	}

	if ( iov->diskFloat == H5T_NATIVE_FLOAT )
	    writeSet(iov->set_id, iov->s, H5T_NATIVE_FLOAT, iov->iOffset, iov->nBuffered, 1 );
	else
	    writeSet(iov->set_id, iov->d, H5T_NATIVE_DOUBLE, iov->iOffset, iov->nBuffered, 1 );
	iov->iOffset += iov->nBuffered;
	iov->nBuffered = 0;
    }
}

/* Writes the fields common to all particles */
static void flushBase( IOHDF5 io, IOBASE *Base,
		       void (*flushExtra)( IOHDF5 io, IOBASE *Base ) )
{
    IOHDF5V iov;

    assert( io != NULL );
    if ( Base->nBuffered ) {
	hid_t dataType = sizeof(PINDEX)==4
	    ? H5T_NATIVE_UINT32 : H5T_NATIVE_UINT64;

	if ( Base->group_id == H5I_INVALID_HID ) {
	    Base->group_id = CreateGroup(io->fileID,Base->szGroupName);


	    writeAttribute( Base->group_id, ATTR_IORDER,
			    dataType, &Base->Order.iStart );
	}

	if ( Base->setR_id == H5I_INVALID_HID ) {
	    Base->setR_id = newSet(
		Base->group_id, FIELD_POSITION,
		io->iChunkSize, 0, 3, io->diskFloat_R );
	    Base->setV_id = newSet(
		Base->group_id, FIELD_VELOCITY,
		io->iChunkSize, 0, 3, io->diskFloat_V );
	    if ( Base->Order.iOrder != NULL ) {
		//assert( sizeof(PINDEX) == 4 );
		Base->Order.setOrder_id = newSet(
		    Base->group_id, FIELD_ORDER,
		    io->iChunkSize, 0, 1, dataType);
	    }
	    if ( Base->Class.piClass != NULL ) {
		Base->Class.setClass_id = newSet(
		    Base->group_id, FIELD_CLASS,
		    io->iChunkSize, 0, 1, H5T_NATIVE_UINT8);
	    }
	}

	writeSet( Base->setR_id, Base->R, io->memFloat,
		  Base->iOffset, Base->nBuffered, 3 );
	writeSet( Base->setV_id, Base->V, io->memFloat,
		  Base->iOffset, Base->nBuffered, 3 );

	if ( Base->Order.iOrder != NULL ) {
	    writeSet( Base->Order.setOrder_id, Base->Order.iOrder,
		      dataType, Base->iOffset, Base->nBuffered, 1 );
	}
	if ( Base->Class.piClass != NULL ) {
	    writeSet( Base->Class.setClass_id, Base->Class.piClass,
		      H5T_NATIVE_UINT8, Base->iOffset, Base->nBuffered, 1 );
	}

	if ( flushExtra != NULL ) flushExtra(io,Base);
	Base->iOffset += Base->nBuffered;
	Base->nBuffered = 0;
    }

    for( iov = io->vectorList; iov!=NULL; iov = iov->next ) {
	flushVector(iov);
    }


}

/* Write any accumulated dark particles */
static void flushDark( IOHDF5 io )
{
    flushBase( io, &io->darkBase, 0 );
}

/* Write any accumulated gas particles */
static void flushGas( IOHDF5 io )
{
    flushBase( io, &io->gasBase, 0 );
}

/* Write any accumulated star particles */
static void flushStar( IOHDF5 io )
{
    flushBase( io, &io->starBase, 0 );
}

static hid_t makeClassType(hid_t floatType, int bStart) {
    hid_t tid;
    hid_t dataType = sizeof(PINDEX)==4
	? H5T_NATIVE_UINT32 : H5T_NATIVE_UINT64;

    tid = H5Tcreate (H5T_COMPOUND, sizeof(classEntry));
    H5assert(tid);
    H5Tinsert(tid,"class",HOFFSET(classEntry,iClass), H5T_NATIVE_UINT8);
    H5Tinsert(tid,"mass", HOFFSET(classEntry,fMass), floatType);
    H5Tinsert(tid,"soft", HOFFSET(classEntry,fSoft), floatType);
    if ( bStart )
	H5Tinsert(tid,"start",HOFFSET(classEntry,iOrderStart), dataType);
    return tid;
}

static void writeClassTable(IOHDF5 io, IOBASE *Base ) {
    hid_t tid, set;

    if ( Base->nTotal > 0 ) {
	tid = makeClassType( io->memFloat, Base->Class.piClass==NULL );

	set = newSet(Base->group_id, FIELD_CLASSES,
		     io->iChunkSize, Base->Class.nClasses, 1, tid );

	writeSet( set, Base->Class.Class, tid,
		  0, Base->Class.nClasses, 1 );

	H5Dclose(set);
	H5Tclose(tid);
    }
}

static void readClassTable( IOHDF5 io, IOBASE *Base ) {
    hid_t tid, set;

    set = openSet( Base->group_id, FIELD_CLASSES );
    if ( set != H5I_INVALID_HID ) {

	if ( Base->Class.setClass_id != H5I_INVALID_HID 
	     && Base->Class.piClass == NULL ) {
	    Base->Class.piClass = (uint8_t *)malloc( io->iChunkSize * sizeof(uint8_t) );
	    assert(Base->Class.piClass != NULL );
	}
	tid = makeClassType( io->memFloat, Base->Class.piClass==NULL );
	Base->Class.nClasses = getSetSize(set);
	readSet( set, Base->Class.Class, tid,
		 0, Base->Class.nClasses, 1 );
	H5Dclose(set);
	H5Tclose(tid);
    }

}

/* Read an attribute from a group */
int ioHDF5ReadAttribute( IOHDF5 io, const char *name,
			 hid_t dataType, void *data )
{
    return readAttribute( io->parametersID, name, dataType, data );
}

/* Add an attribute to a group */
void ioHDF5WriteAttribute( IOHDF5 io, const char *name,
			   hid_t dataType, void *data )
{
    writeAttribute( io->parametersID, name, dataType, data );
}

void ioHDF5Flush( IOHDF5 io )
{
    flushDark(io);
    flushGas(io);
    flushStar(io);
    writeClassTable(io,&io->darkBase);
    writeClassTable(io,&io->gasBase);
    writeClassTable(io,&io->starBase);
}


static void baseInitialize( IOHDF5 io, IOBASE *Base,
			    hid_t locID, const char *group)
{
    Base->iOffset = 0;
    Base->iIndex = 0;
    Base->nBuffered = 0;


    Base->Order.iStart = Base->Order.iNext = 0;
    Base->Order.iOrder = NULL;
    Base->R = Base->V = NULL;
    Base->Class.fMass = Base->Class.fSoft = NULL;
    Base->Class.piClass = NULL;
    Base->Class.nClasses=0;

    Base->group_id = H5Gopen( locID, group );

    /* Oh.  We are reading this file. */
    if ( Base->group_id != H5I_INVALID_HID ) {
	Base->setR_id = openSet(Base->group_id,FIELD_POSITION);
	H5assert(Base->setR_id);
	Base->nTotal = getSetSize(Base->setR_id);
	Base->setV_id = openSet(Base->group_id,FIELD_VELOCITY);
	H5assert(Base->setV_id);
	assert( Base->nTotal == getSetSize(Base->setV_id) );
	Base->Order.setOrder_id = openSet(Base->group_id,FIELD_ORDER);
	Base->Class.setClass_id = openSet(Base->group_id,FIELD_CLASS);
	//if ( Base->Class.setClass_id != H5I_INVALID_HID ) {
	    readClassTable( io, Base );
	    //}
    }

    /* No group: we have to create this later.  We delay until later because
       it is possible that we won't have any of this type of particle. */
    else {
	/*Base->group_id = H5I_INVALID_HID;*/
	Base->setR_id = Base->setV_id = H5I_INVALID_HID;
	Base->Order.setOrder_id = H5I_INVALID_HID;
	Base->Class.setClass_id = H5I_INVALID_HID;
	Base->nTotal = 0;
    }

    assert( strlen(group) < sizeof(Base->szGroupName) );
    strcpy( Base->szGroupName, group );

}

IOHDF5 ioHDF5Initialize( hid_t fileID, hid_t iChunkSize, int bDouble )
{
    IOHDF5 io;
    H5E_auto_t save_func;
    void *     save_data;

    io = (IOHDF5)malloc( sizeof(struct ioHDF5) ); assert(io!=NULL);
    io->fileID = fileID;
    io->iChunkSize = iChunkSize;

    io->bRead = io->bWrite = 0;

    io->vectorList = NULL;

    /* This is the native type of FLOAT values - normally double */
    io->memFloat = sizeof(FLOAT)==sizeof(float)
	? H5T_NATIVE_FLOAT : H5T_NATIVE_DOUBLE;
    io->diskFloat_R = io->diskFloat_V = bDouble ? io->memFloat : H5T_NATIVE_FLOAT;

    /* The group might already exist */
    H5Eget_auto(&save_func,&save_data);
    H5Eset_auto(0,0);

    io->parametersID = H5Gopen( fileID, GROUP_PARAMETERS );
    if ( io->parametersID == H5I_INVALID_HID ) {
	io->parametersID = CreateGroup( fileID, GROUP_PARAMETERS );
    }

    baseInitialize( io, &io->darkBase, fileID, GROUP_DARK );
    baseInitialize( io, &io->gasBase,  fileID, GROUP_GAS  );
    baseInitialize( io, &io->starBase, fileID, GROUP_STAR );

    H5Eset_auto(save_func,save_data);

    return io;
}

static void baseFinish( IOBASE *Base )
{
    if ( Base->group_id != H5I_INVALID_HID ) H5Gclose(Base->group_id);
    if ( Base->setR_id != H5I_INVALID_HID ) H5Dclose(Base->setR_id);
    if ( Base->setV_id != H5I_INVALID_HID ) H5Dclose(Base->setV_id);
    if ( Base->Order.setOrder_id != H5I_INVALID_HID )
	H5Dclose(Base->Order.setOrder_id);
    if ( Base->Class.setClass_id != H5I_INVALID_HID )
	H5Dclose(Base->Class.setClass_id);

    if ( Base->Order.iOrder != NULL ) free( Base->Order.iOrder );
    if ( Base->R != NULL ) free( Base->R );
    if ( Base->V != NULL ) free( Base->V );
    if ( Base->Class.fMass != NULL ) free( Base->Class.fMass );
    if ( Base->Class.fSoft != NULL ) free( Base->Class.fSoft );
}

void ioHDF5Finish( IOHDF5 io )
{
    IOHDF5V iov,nextv;
    assert( io != NULL );

    if ( io->bWrite )
	ioHDF5Flush(io);

    for( iov = io->vectorList; iov!=NULL; iov = nextv ) {
	nextv = iov->next;
	if ( iov->s ) free(iov->s);
	if ( iov->d ) free(iov->d);
	if ( iov->set_id != H5I_INVALID_HID )
	    H5Dclose( iov->set_id );
	free(iov);
    }
    io->vectorList = NULL;

    baseFinish( &io->darkBase );
    baseFinish( &io->gasBase );
    baseFinish( &io->starBase );

    free( io );
}


/* Create the class set if it doesn't already exist and flush out the table. */
static void createClass(IOHDF5 io, IOBASE *Base)
{
    PINDEX i, j, n, o;
    IOCLASS *Class = &Base->Class;

    /* We already created the set */
    if ( Class->piClass != NULL ) return;

    Class->piClass = (uint8_t *)malloc( io->iChunkSize * sizeof(uint8_t) );
    assert(Class->piClass!=NULL);

    /* If the group exists, we will have to write */
    if ( Class->setClass_id == H5I_INVALID_HID 
	 && Base->group_id != H5I_INVALID_HID ) {
	Class->setClass_id = newSet(
	    Base->group_id, FIELD_CLASS,
	    io->iChunkSize, 0, 1, H5T_NATIVE_UINT8);
    }

    for( i=n=o=0; i<Class->nClasses; i++ ) {
	PINDEX s, e;

	s = Class->Class[i].iOrderStart;
	e = i==Class->nClasses-1 ? Base->nTotal : Class->Class[i+1].iOrderStart;

	for( j=s; j<e; j++ ) {
	    Class->piClass[n++] = i;
	    if ( n == io->iChunkSize ) {
		writeSet( Base->Class.setClass_id, Base->Class.piClass,
		      H5T_NATIVE_UINT8, o, n, 1 );
		o += n;
		n = 0;
	    }
	}
    }
}

static void addClass( IOHDF5 io, IOBASE *Base,
		      PINDEX iOrder, FLOAT fMass, FLOAT fSoft )
{
    IOCLASS *Class = &Base->Class;
    uint_fast32_t i;

    /* See if we already have this class: Mass/Softening pair */
    for( i=0; i<Class->nClasses; i++ ) {
	if ( Class->Class[i].fMass == fMass && Class->Class[i].fSoft == fSoft )
	    break;
    }

    /* Case 1: This is a new class */
    if ( i == Class->nClasses ) {
	assert( Class->nClasses < 256 ); /*TODO: handle this case */
	Class->Class[i].iClass = i;
	Class->Class[i].iOrderStart = iOrder;
	Class->Class[i].fMass = fMass;
	Class->Class[i].fSoft = fSoft;
	Class->Class[i].nCount= 0;
	Class->nClasses++;
	if ( Class->piClass != NULL )
	    Class->piClass[Base->nBuffered] = i;
    }

    /* Case 2: This was the last class, and we might be compressing */
    else if ( i == Class->nClasses - 1 && Class->piClass==NULL ) {
    }

    /* Case 3: A match, but a prior class */
    else {
	createClass(io,Base);
	Class->piClass[Base->nBuffered] = i;
    }
    Class->Class[i].nCount++;
}

static void addOrder( IOHDF5 io, IOBASE *Base,
		      PINDEX iOrder )
{
    IOORDER *Order = &Base->Order;
    PINDEX i;
    uint_fast32_t n;

    /* If we still think that the particles are in order */
    if ( Order->iOrder==NULL ) {
	if ( Base->nTotal == 1 )
	    Order->iStart = Order->iNext = iOrder;
	/* So far, so good.  Bump up the "next expected" value */
	if ( iOrder == Order->iNext )
	    Order->iNext++;
	/* Darn...  we have out of order particles, so write the array as
	   it should be up to this point. */
	else {
	    int iOffset = 0;
	    hid_t dataType = sizeof(PINDEX)==4
		? H5T_NATIVE_UINT32 : H5T_NATIVE_UINT64;

	    Order->iOrder = (PINDEX*)malloc( io->iChunkSize * sizeof(PINDEX) );
	    assert(Order->iOrder!=NULL);

	    if ( Order->setOrder_id == H5I_INVALID_HID 
		&& Base->group_id != H5I_INVALID_HID ) {
		assert( sizeof(PINDEX) == 4 );
		Order->setOrder_id = newSet(
		    Base->group_id, FIELD_ORDER,
		    io->iChunkSize, 0, 1, dataType);
	    }

	    n = 0;
	    for( i=Order->iStart; i<Order->iNext; i++ ) {
		Order->iOrder[n++] = i;
		if ( n == io->iChunkSize ) {
		    writeSet( Order->setOrder_id, Order->iOrder,
			      dataType,iOffset,n,1);
		    iOffset += n;
		    n = 0;
		}
	    }
	    assert( n == Base->nBuffered );
	    Order->iOrder[Base->nBuffered] = iOrder;
	}
    }

    /* Particles are out of order for sure: just buffer iOrder for writing */
    else {
	Order->iOrder[Base->nBuffered] = iOrder;
    }
}

/* If the structures have not been allocated, do so now */
static void allocateBase( IOHDF5 io, IOBASE *Base )
{
    if ( Base->R == NULL ) {
	Base->R    = (ioV3*)malloc( io->iChunkSize * sizeof(ioV3) );
	assert( Base->R != NULL );
	Base->V    = (ioV3*)malloc( io->iChunkSize * sizeof(ioV3) );
	assert( Base->V != NULL );
    }
}


static int getBase( IOHDF5 io, IOBASE *Base, PINDEX *iOrder,
		    FLOAT *r, FLOAT *v,
		    FLOAT *fMass, FLOAT *fSoft )
{
    uint_fast32_t iClass;

    assert(Base->setR_id!=H5I_INVALID_HID);
    assert(Base->setV_id!=H5I_INVALID_HID);
    allocateBase(io,Base);

    assert( io->bWrite == 0 );
    io->bRead = 1;

    /* If we have to read more from the file */
    if ( Base->nBuffered == Base->iIndex ) {
	hid_t dataType = sizeof(PINDEX)==4
	    ? H5T_NATIVE_UINT32 : H5T_NATIVE_UINT64;
	hsize_t N;
	Base->iOffset += Base->iIndex;
	Base->iIndex = Base->nBuffered = 0;
	if ( Base->iOffset >= Base->nTotal ) return 0;

 	N = Base->nTotal - Base->iOffset;
	if ( N > io->iChunkSize )
	    Base->nBuffered = io->iChunkSize;
	else
	    Base->nBuffered = N;

	readSet( Base->setR_id, Base->R, io->memFloat,
		 Base->iOffset, Base->nBuffered, 3 );
	readSet( Base->setV_id, Base->V, io->memFloat,
		 Base->iOffset, Base->nBuffered, 3 );

	if ( Base->Order.setOrder_id != H5I_INVALID_HID ) {
	    if ( Base->Order.iOrder == NULL ) {
		Base->Order.iOrder = (PINDEX*)malloc( io->iChunkSize * sizeof(PINDEX) );
		assert( Base->Order.iOrder != NULL );
	    }
	    readSet( Base->Order.setOrder_id, Base->Order.iOrder,
		     dataType,
		     Base->iOffset, Base->nBuffered, 1 );
	}
	if ( Base->Class.setClass_id != H5I_INVALID_HID ) {
	    if ( Base->Class.piClass == NULL ) {
		Base->Class.piClass=(uint8_t*)malloc( io->iChunkSize * sizeof(uint8_t) );
		assert( Base->Class.piClass != NULL );
	    }
	    readSet( Base->Class.setClass_id, Base->Class.piClass,
		     H5T_NATIVE_UINT8, Base->iOffset, Base->nBuffered, 1 );
	}
    }
    *iOrder = Base->iOffset + Base->iIndex; /*FIXME: */
    r[0] = Base->R[Base->iIndex].v[0];
    r[1] = Base->R[Base->iIndex].v[1];
    r[2] = Base->R[Base->iIndex].v[2];
    v[0] = Base->V[Base->iIndex].v[0];
    v[1] = Base->V[Base->iIndex].v[1];
    v[2] = Base->V[Base->iIndex].v[2];

    if ( Base->Class.piClass == NULL ) {
	assert(Base->Class.nClasses>=1);
	for( iClass=0; iClass<Base->Class.nClasses; iClass++ )
	    if ( Base->Class.Class[iClass].iOrderStart > *iOrder )
		break;
	assert( iClass>0 );
	--iClass;
    }
    else {
	iClass = Base->Class.piClass[Base->iIndex];
	assert( iClass < Base->Class.nClasses );
    }
    *fMass = Base->Class.Class[iClass].fMass;
    *fSoft = Base->Class.Class[iClass].fSoft;

    Base->iIndex++;
    return 1;
}

static void addBase( IOHDF5 io, IOBASE *Base, PINDEX iOrder,
		     void (*flush)(IOHDF5 io),
		     const FLOAT *r, const FLOAT *v,
		     FLOAT fMass, FLOAT fSoft )
{
    int i;

    assert( io != NULL );
    assert( Base != NULL );

    assert( io->bRead == 0 );
    io->bWrite = 1;
    Base->nTotal++;

    allocateBase(io,Base);

    /* If the particles are not in order, then we have to store iOrder */
    addOrder( io, Base, iOrder );

    for( i=0; i<3; i++ ) {
	Base->R[Base->nBuffered].v[i] = r[i];
	Base->V[Base->nBuffered].v[i] = v[i];
    }
    addClass( io, Base, iOrder, fMass, fSoft );

    if ( ++Base->nBuffered == io->iChunkSize ) {
	(*flush)(io);
    }
}

static void seekBase( IOHDF5 io, IOBASE *Base, PINDEX Offset ) {
    /*TODO: (optional) seek to chunk boundary */
    Base->nBuffered = Base->iIndex = 0;
    Base->iOffset = Offset;
}


PINDEX ioHDF5DarkCount( IOHDF5 io )
{
    return io->darkBase.nTotal;
}

PINDEX ioHDF5GasCount( IOHDF5 io )
{
    return io->starBase.nTotal;
}

PINDEX ioHDF5StarCount( IOHDF5 io )
{
    return io->starBase.nTotal;
}

IOHDF5V ioHDFF5NewVector( IOHDF5 io, const char *name, int bDouble )
{
    IOHDF5V iov;

    /* Make the new vector and link it onto the vector list */
    iov = malloc( sizeof(struct ioHDF5v) ); assert( iov != NULL );
    iov->next = io->vectorList;
    io->vectorList = iov;
    iov->io = io;
    iov->diskFloat = bDouble ? io->memFloat : H5T_NATIVE_FLOAT;
    iov->set_id = H5I_INVALID_HID;

    assert( strlen(name) < sizeof(iov->name) );
    strcpy( iov->name, name );

    // Filled with cheese - FIXME

    if ( bDouble ) {
	iov->d = (double*)malloc( io->iChunkSize * sizeof(double) );
	iov->s = NULL;
    }
    else {
	iov->s = (float*)malloc( io->iChunkSize * sizeof(float) );
	iov->d = NULL;
    }
    iov->nBuffered = 0;
    iov->iOffset = 0;

    return iov;
}

void ioHDF5AddVector( IOHDF5V iov, PINDEX iOrder, FLOAT v )
{
    if ( iov->diskFloat == H5T_NATIVE_FLOAT )
	iov->s[iov->nBuffered] = v;
    else
	iov->d[iov->nBuffered] = v;
    if ( ++iov->nBuffered == iov->io->iChunkSize )
	flushVector(iov);
}

void ioHDF5AddDark( IOHDF5 io, PINDEX iOrder,
		    const FLOAT *r, const FLOAT *v,
		    FLOAT fMass, FLOAT fSoft, FLOAT fPot )
{
    addBase( io, &io->darkBase, iOrder, flushDark,
	     r, v, fMass, fSoft );
    /*TODO: Save Potential as well */
}

int  ioHDF5GetDark( IOHDF5 io, PINDEX *iOrder,
		    FLOAT *r, FLOAT *v,
		    FLOAT *fMass, FLOAT *fSoft, FLOAT *fPot )
{
    return getBase( io, &io->darkBase, iOrder, r, v, fMass, fSoft );
    *fPot = 0.0;
    /*TODO: Load Potential as well */
}

void ioHDF5SeekDark( IOHDF5 io, PINDEX Offset ) {
    seekBase(io,&io->darkBase,Offset);
}

void ioHDF5AddGas(IOHDF5 io, PINDEX iOrder,
		  const FLOAT *r, const FLOAT *v,
		  FLOAT fMass, FLOAT fSoft, FLOAT fPot,
		  FLOAT fTemp, FLOAT fMetals)
{
    addBase( io, &io->gasBase, iOrder, flushGas,
	     r, v, fMass, fSoft );
    /*TODO: Save fPot, fTemp, fMetals */
}

int ioHDF5GetGas(IOHDF5 io, PINDEX *iOrder,
		  FLOAT *r, FLOAT *v,
		  FLOAT *fMass, FLOAT *fSoft, FLOAT *fPot,
		  FLOAT *fTemp, FLOAT *fMetals)
{
    return getBase( io, &io->gasBase, iOrder, r, v, fMass, fSoft );
    *fPot = *fTemp = *fMetals = 0.0;
    /*TODO: Load fPot, fTemp, fMetals */
}

void ioHDF5SeekGas( IOHDF5 io, PINDEX Offset ) {
    seekBase(io,&io->gasBase,Offset);
}

void ioHDF5AddStar(IOHDF5 io, PINDEX iOrder,
	       const FLOAT *r, const FLOAT *v,
	       FLOAT fMass, FLOAT fSoft, FLOAT fPot,
	       FLOAT fMetals, FLOAT fTForm)
{
    addBase( io, &io->starBase, iOrder, flushStar,
	     r, v, fMass, fSoft );
    /*TODO: Save fPot, fMetals, fTForm */
}

int ioHDF5GetStar( IOHDF5 io, PINDEX *iOrder,
		   FLOAT *r, FLOAT *v,
		   FLOAT *fMass, FLOAT *fSoft, FLOAT *fPot,
		   FLOAT *fMetals, FLOAT *fTForm)
{
    return getBase( io, &io->starBase, iOrder, r, v, fMass, fSoft );
    *fPot = *fMetals = *fTForm = 0.0;
    /*TODO: Load fPot, fMetals, fTForm */
}

void ioHDF5SeekStar( IOHDF5 io, PINDEX Offset ) {
    seekBase(io,&io->starBase,Offset);
}
