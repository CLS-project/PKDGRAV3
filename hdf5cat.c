#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <assert.h>
#include <getopt.h>
#include <stdlib.h>
#include "iohdf5.h"

#define OPT_HELP   'h'
#define OPT_DOUBLE 'd'
#define OPT_OUTPUT 'o'

static void usage(const char *argv0) {
    printf( "Usage: %s [-d,--double] [-o outfile] <infile...>\n", argv0 );
    exit(0);
}

int main( int argc, char *argv[] )
{
    int bHelp = 0;
    int bDoublePos = 0;
    const char *inName, *outName;

    int bFirst = 1;
    FILE *fin;
    hid_t fileID, inID;
    IOHDF5 ioOut, ioIn;
    //IOHDF5V ioPot;

    int i;
    PINDEX iOrder;
    float fTemp;
    double dTemp;
    FLOAT fMass, fSoft, fPot;
    FLOAT r[3];
    FLOAT v[3];

    double dTime;
    int nBodies, nDims, nGas, nDark, nStar, nPad;

    outName = "output.h5";

    //! Parse command line arguments (flags).
    for(;;) {
	int c, option_index=0;
	
	static struct option long_options[] = {
	    { "help",         0, 0, OPT_HELP },
	    { "double",       0, 0, OPT_DOUBLE },
	    { "output",       1, 0, OPT_OUTPUT },
	    { 0,              0, 0, 0 }
	};
	c = getopt_long( argc, argv, "hdo:",
			 long_options, &option_index );
	if ( c == -1 ) break;

	switch(c) {
	case OPT_HELP:
	    bHelp = 1;
	    break;
	case OPT_DOUBLE:
	    bDoublePos = 1;
	    break;
	case OPT_OUTPUT:
	    outName = optarg;
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

    fileID=H5Fcreate(outName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if ( fileID < 0 ) {
	fprintf( stderr, "Unable to create %s\n", outName );
	exit(2);
    }

    ioOut = ioHDF5Initialize( fileID, 32768, 
			      bDoublePos ? IOHDF5_DOUBLE : IOHDF5_SINGLE );
    //ioPot  = ioHDFF5NewVector( ioOut, "potential",IOHDF5_SINGLE );


    while( optind < argc ) {
	inName = argv[optind++];

	printf( "Opening %s\n", inName );


	inID=H5Fopen(inName, H5F_ACC_RDONLY, H5P_DEFAULT);
	if ( fileID < 0 ) {
	    fprintf( stderr, "Unable to open %s\n", inName );
	    exit(2);
	}
	ioIn = ioHDF5Initialize( inID, 32768, IOHDF5_SINGLE );

	assert(ioHDF5GasCount(ioIn)+ioHDF5DarkCount(ioIn)+ioHDF5StarCount(ioIn)<=UINT_MAX);

	if ( bFirst ) {
	    bFirst = 0;
	    assert(ioHDF5ReadAttribute(ioIn,"dTime",H5T_NATIVE_DOUBLE,&dTime));
	    ioHDF5WriteAttribute( ioOut, "dTime", H5T_NATIVE_DOUBLE, &dTime );
	}

	while( ioHDF5GetDark( ioIn, &iOrder, r, v, &fMass, &fSoft, &fPot ) ) {
	    ioHDF5AddDark( ioOut, iOrder, r, v, fMass, fSoft, fPot );
	    //ioHDF5AddVector( ioPot, iOrder, fPot );
	}


	ioHDF5Finish(ioIn);
	H5Fclose( inID );
    }

    ioHDF5Finish(ioOut);

    H5Fclose( fileID );

    return 0;
}
