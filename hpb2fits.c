#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <getopt.h>
#include <assert.h>

#include "fitsio.h"

#define OPT_OUTPUT 'o'
#define OPT_NSIDE 'n'
#define OPT_MASS 'm'

int main(int argc, char *argv[]) {
    uint32_t *uSignal;
    const char *outName = "healpix.fits";
    float fMassUnit = 1.0;
    long nSideHealpix = 0;
    uint64_t nPix;
    uint64_t nParticleCount = 0;
    uint64_t nMaxPixelValue = 0;
    size_t i;

    //! Parse command line arguments (flags).
    for (;;) {
	int c, option_index=0;

	static struct option long_options[] = {
		{ "OUTPUT",       1, 0, OPT_OUTPUT},
		{ "nside",        1, 0, OPT_NSIDE},
		{ "mass",         1, 0, OPT_MASS},
		{ NULL,   0, 0, 0 },
	    };

	c = getopt_long( argc, argv, "o:n:m:",
			 long_options, &option_index );
	if ( c == -1 ) break;

	switch (c) {
	case OPT_OUTPUT:
	    outName = optarg;
	    break;
	case OPT_NSIDE:
	    nSideHealpix = atoi(optarg);
	    break;
	case OPT_MASS:
	    fMassUnit = atof(optarg);
	    break;
	default:
	    fprintf(stderr,"Usage: cat hpb.{0..2047} | %s -o healpix.fits -n 8192 -m 1.0\n",argv[0]);
	    fprintf(stderr,"       -o name of output file\n");
	    fprintf(stderr,"       -n nside Healpix parameter\n");
	    fprintf(stderr,"       -m particle mass unit\n");
	    return 9;
	    }
	}

    if (nSideHealpix==0) {
	fprintf(stderr,"Specify nside with -n\n");
	return 1;
	}
    nPix = 12ull * nSideHealpix*nSideHealpix;    

    uSignal = malloc(sizeof(uint32_t) * nPix);
    if (uSignal == NULL) {
	fprintf(stderr,"Unable to allocate healpix array!\n");
	return 2;
	}

    if (fread(uSignal,sizeof(uint32_t),nPix,stdin) != nPix) {
	fprintf(stderr,"Unable to read count array\n");
	return 3;
	}

    char *tform[1];
    long bzero = 0;

    for(i=0; i<nPix; ++i) {
	if (uSignal[i] > nMaxPixelValue) nMaxPixelValue = uSignal[i];
	nParticleCount += uSignal[i];
	}

    if (nMaxPixelValue <= 0xff) {
	tform[0] = "1B";    /* 8-bit integer */
	bzero = 0;          /* unsigned (yes) */
	}
    else if (nMaxPixelValue <= 0x7fff) {
	tform[0] = "1I";     /* 16-bit integer */
	bzero = 0;          /* signed */
	}
    else if (nMaxPixelValue <= 0xffff) {
	tform[0] = "1I";    /* 16-bit integer */
	bzero = 0x8000;     /* unsigned */
	}
    else if (nMaxPixelValue <= 0x7fffffff) {
	tform[0] = "1J";    /* 16-bit integer */
	bzero = 0;          /* signed */
	}
    else {
	tform[0] = "1J";    /* 32-bit integer */
	bzero = 0x80000000; /* unsigned */
	}

    if (bzero) {
	for(i=0; i<nPix; ++i) uSignal[i] -= bzero;
	}


    fprintf(stderr,"%llu particles in %llu pixels (max value %llu)\n",
	nParticleCount, nPix,nMaxPixelValue);

    fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
    int status=0, hdutype;

    long naxes[] = {0,0};

    char order[9];                 /* HEALPix ordering */
    char *ttype[] = { "SIGNAL" };
    char *tunit[] = { " " };
    char coordsys9[9];

    /* create new FITS file */
    fits_create_file(&fptr, outName, &status);
    fits_create_img(fptr, SHORT_IMG, 0, naxes, &status);
    fits_write_date(fptr, &status);
    fits_movabs_hdu(fptr, 1, &hdutype, &status);
    fits_create_tbl( fptr, BINARY_TBL, 12L*nSideHealpix*nSideHealpix, 1, ttype, tform,
	tunit, "BINTABLE", &status);
    fits_write_key(fptr, TSTRING, "PIXTYPE", "HEALPIX", "HEALPIX Pixelisation",
	&status);
    if (bzero) {
	fits_write_key(fptr, TLONG, "BZERO", &bzero,
	    "Zero point", &status);
	}

    strcpy(order, "RING    ");
    fits_write_key(fptr, TSTRING, "ORDERING", order,
	"Pixel ordering scheme, either RING or NESTED", &status);
    fits_write_key(fptr, TSTRING, "INDXSCHM", "IMPLICIT",
	"Pixel indexing scheme, either IMPLICIT or EXPLICIT", &status);
    fits_write_key(fptr, TLONG, "NSIDE", &nSideHealpix,
	"Resolution parameter for HEALPIX", &status);

    strcpy(coordsys9,"C       ");

    fits_write_key(fptr, TSTRING, "COORDSYS", coordsys9,
	"Pixelisation coordinate system", &status);

    fits_write_comment(fptr,
	"G = Galactic, E = ecliptic, C = celestial = equatorial", &status);

    fits_write_col(fptr, TUINT, 1, 1, 1, nPix, (void *)uSignal,&status);

    fits_close_file(fptr, &status);

    return 0;
    }
