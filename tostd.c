#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <assert.h>
#include <errno.h>
#include <fio.h>

#define BUFFER_SIZE (8*1024*1024)

#define OPT_DOUBLE 'd'
#define OPT_NATIVE 'n'

int main( int argc, char *argv[] ) {
    int bError = 0;
    int bDouble = 0;
    int bNative = 0;
    uint64_t N, nSph, nDark, nStar, i;
    double dTime;

    uint64_t iOrder;
    double r[3], v[3];
    float fMass, fSoft, fPot, fRho, u, fMetals, fTimer;

    FIO fioIn, fioOut;
    FIO_SPECIES eSpecies;
    const char *inName  = 0;
    const char *outName = 0;

    //! Parse command line arguments (flags).
    for (;;) {
	int c, option_index=0;

	static struct option long_options[] = {
		{ "double",       0, 0, OPT_DOUBLE},
		{ "native",       0, 0, OPT_NATIVE},
		{ NULL,   0, 0, 0 },
	    };

	c = getopt_long( argc, argv, "dn",
			 long_options, &option_index );
	if ( c == -1 ) break;

	switch (c) {
	case OPT_DOUBLE:
	    bDouble = 1;
	    break;
	case OPT_NATIVE:
	    bNative = 1;
	    break;
	default:
	    bError = 1;
	    break;
	    }
	}

    if ( optind < argc )
	inName = argv[optind++];
    else {
	fprintf(stderr, "Missing input file\n" );
	bError = 1;
	}

    if ( optind < argc )
	outName = argv[optind++];
    else {
	fprintf(stderr, "Missing Tipsy output file\n" );
	bError = 1;
	}

    if ( bError ) {
	fprintf(stderr, "Usage: %s [-p] <input> <outtipsy>\n"
		"  -d,--double  Output double precision positions\n"
		"  -n,--native  Output a native tipsy binary\n",
		argv[0] );
	exit(1);
	}

    fioIn = fioOpen(inName,0.0,0.0);
    if (fioIn==NULL) {
	perror(inName);
	exit(errno);
	}
    N     = fioGetN(fioIn,FIO_SPECIES_ALL);
    nSph  = fioGetN(fioIn,FIO_SPECIES_SPH);
    nDark = fioGetN(fioIn,FIO_SPECIES_DARK);
    nStar = fioGetN(fioIn,FIO_SPECIES_STAR);
    if (!fioGetAttr(fioIn,"dTime",FIO_TYPE_DOUBLE,&dTime)) dTime = 0.0;

    fioOut = fioTipsyCreate(outName,bDouble,!bNative,dTime,nSph,nDark,nStar);
    if (fioOut==NULL) {
	perror(outName);
	exit(errno);
	}

    for( i=0; i<N; i++ ) {
        eSpecies = fioSpecies(fioIn);
        switch(eSpecies) {
        case FIO_SPECIES_SPH:
            fioReadSph(fioOut,&iOrder,r,v,&fMass,&fSoft,&fPot,&fRho,&u,&fMetals);
	    fioWriteSph(fioOut,iOrder,r,v,fMass,fSoft,fPot,fRho,u,fMetals);
            break;
        case FIO_SPECIES_DARK:
            fioReadDark(fioIn,&iOrder,r,v,&fMass,&fSoft,&fPot);
            fioWriteDark(fioOut,iOrder,r,v,fMass,fSoft,fPot);
            break;
        case FIO_SPECIES_STAR:
            fioReadStar(fioIn,&iOrder,r,v,&fMass,&fSoft,&fPot,&fMetals,&fTimer);
            fioWriteStar(fioOut,iOrder,r,v,fMass,fSoft,fPot,fMetals,fTimer);
            break;
        default:
            fprintf(stderr,"Unsupported particle type: %d\n",eSpecies);
            assert(0);
            }


	}

    fioClose(fioOut);
    fioClose(fioIn);



    return 0;
    }
