#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <assert.h>
#include <getopt.h>
#include <rpc/xdr.h>
#include "iohdf5.h"

#define OPT_HELP   'h'
#define OPT_DOUBLE 'd'

static void usage(const char *argv0) {
    printf( "Usage: %s [-d,--double] <infile> <outfile>\n", argv0 );
    exit(0);
}

int main( int argc, char *argv[] )
{
    int bHelp = 0;
    int bDoublePos = 0;
    const char *inName, *outName;

    XDR xdr;
    FILE *fin;
    hid_t fileID;
    IOHDF5 io;
    IOHDF5V ioPot;

    int i, iOrder;
    float fTemp;
    double dTemp;
    FLOAT fMass, fSoft, fPot;
    FLOAT r[3];
    FLOAT v[3];

    double dTime;
    int nBodies, nDims, nGas, nDark, nStar, nPad;

    //! Parse command line arguments (flags).
    for(;;) {
	int c, option_index=0;
	
	static struct option long_options[] = {
	    { "help",         0, 0, OPT_HELP },
	    { "double",       0, 0, OPT_DOUBLE },
	    { 0,              0, 0, 0 }
	};
	c = getopt_long( argc, argv, "hd",
			 long_options, &option_index );
	if ( c == -1 ) break;

	switch(c) {
	case OPT_HELP:
	    bHelp = 1;
	    break;
	case OPT_DOUBLE:
	    bDoublePos = 1;
	    break;
	default:
	    usage(argv[0]);
	    break;
	}
    }
    if ( bHelp ) {
	usage(argv[0]);
	exit(0);
    }

    if ( optind>=argc ) usage(argv[0]);
    inName = argv[optind++];

    if ( optind>=argc ) usage(argv[0]);
    outName = argv[optind++];

    if ( optind!=argc ) usage(argv[0]);



    fin = fopen( inName, "rb" );
    if ( fin == NULL ) {
	fprintf( stderr, "Unable to open %s for reading\n", inName );
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

    fileID=H5Fcreate(outName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if ( fileID < 0 ) {
	fprintf( stderr, "Unable to create %s\n", outName );
	exit(2);
    }

    io = ioHDF5Initialize( fileID, 32768, IOHDF5_SINGLE );
    ioPot  = ioHDFF5NewVector( io, "potential",IOHDF5_SINGLE );
    ioHDF5WriteAttribute( io, "dTime", H5T_NATIVE_DOUBLE, &dTime );

    iOrder = 0;
    for( i=0; i<nDark; i++ ) {
	xdr_float(&xdr,&fTemp); fMass = fTemp;
	if ( bDoublePos ) {
	    xdr_double(&xdr,&dTemp); r[0] = dTemp;
	    xdr_double(&xdr,&dTemp); r[1] = dTemp;
	    xdr_double(&xdr,&dTemp); r[2] = dTemp;
	}
	else {
	    xdr_float(&xdr,&fTemp); r[0] = fTemp;
	    xdr_float(&xdr,&fTemp); r[1] = fTemp;
	    xdr_float(&xdr,&fTemp); r[2] = fTemp;
	}
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
