#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <assert.h>
#include <rpc/xdr.h>
#include "iohdf5.h"

int main( int argc, char *argv[] )
{
    XDR xdr;
    FILE *fin;
    hid_t fileID;
    IOHDF5 io;
    IOHDF5V ioPot;

    int i, iOrder;
    float fTemp;
    FLOAT fMass, fSoft, fPot;
    FLOAT r[3];
    FLOAT v[3];

    double dTime;
    int nBodies, nDims, nGas, nDark, nStar, nPad;


    if ( argc != 3 ) {
	fprintf( stderr, "Usage: %s <instd> <outhdf5>\n", argv[0] );
	exit(0);
    }

    fin = fopen( argv[1], "rb" );
    if ( fin == NULL ) {
	fprintf( stderr, "Unable to open %s for reading\n", argv[1] );
	exit(1);
    }
    xdrstdio_create(&xdr,fin,XDR_DECODE);

    /* Header */
    assert( xdr_double(&xdr,&dTime) );
    assert( xdr_int(&xdr,&nBodies) );
    assert( xdr_int(&xdr,&nDims) );
    assert( xdr_int(&xdr,&nGas) );
    assert( xdr_int(&xdr,&nDark) );
    assert( xdr_int(&xdr,&nStar) );
    assert( xdr_int(&xdr,&nPad) );

    fprintf( stderr, "Converting %d particles, %d gas, %d dark, %d star\n",
	     nBodies, nGas, nDark, nStar );

    fileID=H5Fcreate(argv[2], H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if ( fileID < 0 ) {
	fprintf( stderr, "Unable to create %s\n", argv[2] );
	exit(2);
    }

    io = ioHDF5Initialize( fileID, 32768, IOHDF5_SINGLE );
    ioPot  = ioHDFF5NewVector( io, "potential",IOHDF5_SINGLE );
    ioHDF5WriteAttribute( io, "dTime", H5T_NATIVE_DOUBLE, &dTime );

    iOrder = 0;
    for( i=0; i<nDark; i++ ) {
	xdr_float(&xdr,&fTemp); fMass = fTemp;
	xdr_float(&xdr,&fTemp); r[0] = fTemp;
	xdr_float(&xdr,&fTemp); r[1] = fTemp;
	xdr_float(&xdr,&fTemp); r[2] = fTemp;
	xdr_float(&xdr,&fTemp); v[0] = fTemp;
	xdr_float(&xdr,&fTemp); v[1] = fTemp;
	xdr_float(&xdr,&fTemp); v[2] = fTemp;
	xdr_float(&xdr,&fTemp); fSoft = fTemp;
	xdr_float(&xdr,&fTemp); fPot  = fTemp;
	ioHDF5AddDark( io, iOrder, r, v, fMass, fSoft, fPot );
	ioHDF5AddVector( ioPot, iOrder, fPot );
	iOrder++;
    }

    ioHDF5Finish(io);

    H5Fclose( fileID );
    xdr_destroy(&xdr);

    return 0;
}
