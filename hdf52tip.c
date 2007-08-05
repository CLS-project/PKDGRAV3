#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <assert.h>
#include <rpc/xdr.h>
#include "iohdf5.h"

int main( int argc, char *argv[] )
{
    int bError = 0;
    const char *inName  = 0;
    const char *outName = 0;

    hid_t fileID;
    IOHDF5 io;
    FILE  *fout;
    XDR xdr;

    double dTime;
    float  fTemp;
    PINDEX iOrder;
    int nBodies, nDims, nGas, nDark, nStar, nPad;
    FLOAT fMass, fSoft, fPot;
    FLOAT r[3];
    FLOAT v[3];

    //! Parse command line arguments (flags).
    for(;;) {
        int c, option_index=0;

        static struct option long_options[] = {
            { NULL,   0, 0, 0 },
        };

        c = getopt_long( argc, argv, "",
                         long_options, &option_index );
        if ( c == -1 ) break;

        switch(c) {
	default:
	    bError = 1;
	    break;
	}
    }

    if ( optind < argc )
	inName = argv[optind++];
    else {
	fprintf(stderr, "Missing HDF5 input file\n" );
	bError = 1;
    }

    if ( optind < argc )
	outName = argv[optind++];
    else {
	fprintf(stderr, "Missing Tipsy output file\n" );
	bError = 1;
    }

    if ( bError ) {
	fprintf(stderr, "Usage: %s [-p] <inhdf5> <outtipsy>\n"
		"  -p,--pipe  Create a pipe instead of writing a regular file\n",
		argv[0] );
	exit(1);
    }


    fileID=H5Fopen(inName, H5F_ACC_RDONLY, H5P_DEFAULT);
    if ( fileID < 0 ) {
	fprintf( stderr, "Unable to open %s\n", inName );
	exit(2);
    }
    io = ioHDF5Initialize( fileID, 32768, IOHDF5_SINGLE );

    fout = fopen( outName, "wb" );
    if ( fout == NULL ) {
	fprintf( stderr, "Unable to open %s for writing\n", outName );
	exit(2);
    }
    xdrstdio_create(&xdr,fout,XDR_ENCODE);

    assert(ioHDF5ReadAttribute( io, "dTime", H5T_NATIVE_DOUBLE, &dTime ));

    nDims = 3;
    nGas  = ioHDF5GasCount(io);
    nDark = ioHDF5DarkCount(io);
    nStar = ioHDF5StarCount(io);
    nPad  = 0;
    nBodies = nGas + nDark + nStar;

    /* Tipsy Header */
    assert( xdr_double(&xdr,&dTime) );
    assert( xdr_int(&xdr,&nBodies) );
    assert( xdr_int(&xdr,&nDims) );
    assert( xdr_int(&xdr,&nGas) );
    assert( xdr_int(&xdr,&nDark) );
    assert( xdr_int(&xdr,&nStar) );
    assert( xdr_int(&xdr,&nPad) );

    /* Okay, we should really do this */
    assert( nGas == 0 );
    assert( nStar == 0 );

    while( ioHDF5GetDark( io, &iOrder, r, v, &fMass, &fSoft, &fPot ) ) {
	fTemp = fMass; xdr_float(&xdr,&fTemp);
	fTemp = r[0];  xdr_float(&xdr,&fTemp);
	fTemp = r[1];  xdr_float(&xdr,&fTemp);
	fTemp = r[2];  xdr_float(&xdr,&fTemp);
	fTemp = v[0];  xdr_float(&xdr,&fTemp);
	fTemp = v[1];  xdr_float(&xdr,&fTemp);
	fTemp = v[2];  xdr_float(&xdr,&fTemp);
	fTemp = fSoft; xdr_float(&xdr,&fTemp);
	fTemp = fPot;  xdr_float(&xdr,&fTemp);
    }

    xdr_destroy(&xdr);
    fclose(fout);
    H5Fclose(fileID);

    return 0;
}
