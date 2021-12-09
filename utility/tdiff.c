/*  This file is part of PKDGRAV3 (http://www.pkdgrav.org/).
 *  Copyright (c) 2001-2018 Joachim Stadel & Douglas Potter
 *
 *  PKDGRAV3 is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  PKDGRAV3 is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with PKDGRAV3.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#else
#include "pkd_config.h"
#endif
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <errno.h>
#include <math.h>
#include "fio.h"

#define BUFFER_SIZE (8*1024*1024)

#define OPT_POSITION 'p'
#define OPT_VELOCITY 'v'

static double fix(double v) {
    if ( v >= 0.5 ) return v - 1.0;
    else if ( v < -0.5 ) return v + 1.0;
    else return v;
    }

double mag(double *r1, double *r2, int bFix) {
    double dx, dy, dz;

    dx = r1[0] - r2[0];
    dy = r1[1] - r2[1];
    dz = r1[2] - r2[2];
    if (bFix) {
	dx = fix(dx);
	dy = fix(dy);
	dz = fix(dz);
	}
    return sqrt(dx*dx + dy*dy + dz*dz);
    }



int main( int argc, char **argv ) {
    int bError = 0;
    int bVelocity = 0;
    int bPosition = 0;
    uint64_t N1, N2, nSph1, nSph2, nDark1, nDark2, nStar1, nStar2, i;
    double dTime1, dTime2;

    uint64_t iOrder1;
    double r1[3], v1[3];
    float fMass1, fSoft1, fPot1, fRho1, u1, fMetals1, fTimer1, afOtherData1[3];

    uint64_t iOrder2;
    double r2[3], v2[3];
    float fMass2, fSoft2, fPot2, fRho2, u2, fMetals2, fTimer2, afOtherData2[3];

    double v, vMin, vMax;

    FIO fioIn1, fioIn2;
    FIO_SPECIES eSpecies;
    const char *inName1, *inName2;

    for (;;) {
	int c, option_index=0;

	static struct option long_options[] = {
		{ "position",   0, 0, OPT_POSITION },
		{ "velocity",   0, 0, OPT_VELOCITY },
		{ NULL,   0, 0, 0 },
	    };

	c = getopt_long( argc, argv, "vp",
			 long_options, &option_index );
	if ( c == -1 ) break;

	switch (c) {
	case OPT_POSITION:
	    bPosition = 1;
	    break;
	case OPT_VELOCITY:
	    bVelocity = 1;
	    break;
	default:
	    bError = 1;
	    break;
	    }
	}

    if ( optind < argc ) {
	inName1 = argv[optind++];
	}
    else {
	fprintf(stderr, "Missing input file(s)\n" );
	bError = 1;
	}

    if ( optind < argc )
	inName2 = argv[optind++];
    else {
	fprintf(stderr, "Missing input file\n" );
	bError = 1;
	}

    if (bPosition && bVelocity) {
	fprintf(stderr,"Specify --position or --velocity, but not both\n");
	bError = 1;
	}
    else if ( !bPosition && !bVelocity) bPosition = 1;

    if ( bError ) {
	fprintf(stderr, "Usage: %s [-vp] <file1> <file2>\n",
		argv[0] );
	exit(1);
	}

    fioIn1 = fioOpen(inName1,0.0,0.0);
    if (fioIn1==NULL) {
	perror(inName1);
	exit(errno);
	}
    fioIn2 = fioOpen(inName2,0.0,0.0);
    if (fioIn2==NULL) {
	perror(inName2);
	exit(errno);
	}

    N1     = fioGetN(fioIn1,FIO_SPECIES_ALL);
    nSph1  = fioGetN(fioIn1,FIO_SPECIES_SPH);
    nDark1 = fioGetN(fioIn1,FIO_SPECIES_DARK);
    nStar1 = fioGetN(fioIn1,FIO_SPECIES_STAR);
    if (!fioGetAttr(fioIn1,"dTime",FIO_TYPE_DOUBLE,&dTime1)) dTime1 = 0.0;

    N2     = fioGetN(fioIn2,FIO_SPECIES_ALL);
    nSph2  = fioGetN(fioIn2,FIO_SPECIES_SPH);
    nDark2 = fioGetN(fioIn2,FIO_SPECIES_DARK);
    nStar2 = fioGetN(fioIn2,FIO_SPECIES_STAR);
    if (!fioGetAttr(fioIn2,"dTime",FIO_TYPE_DOUBLE,&dTime2)) dTime2 = 0.0;

    fprintf(stderr,"Comparing: %s N=%lu, nSph=%lu, nDark=%lu, nStar=%lu, dTime=%g\n",
	    inName1, N1, nSph1, nDark1, nStar1, dTime1);
    fprintf(stderr,"     with: %s N=%lu, nSph=%lu, nDark=%lu, nStar=%lu, dTime=%g\n",
	    inName2, N2, nSph2, nDark2, nStar2, dTime2);
    if ( N1!=N2 || nSph1!=nSph2 || nDark1!=nDark2 || nStar1!=nStar2) {
	fprintf(stderr,"File headers do not match!\n");
	exit(1);
	}

    printf("%lu\n",N1);

    vMin = HUGE_VAL;
    vMax = -HUGE_VAL;

    for( i=0; i<N1; i++ ) {
        eSpecies = fioSpecies(fioIn1);
        switch(eSpecies) {
        case FIO_SPECIES_SPH:
            fioReadSph(fioIn1,&iOrder1,r1,v1,&fMass1,&fSoft1,&fPot1,&fRho1,&u1,
		       &fMetals1,afOtherData1);
            fioReadSph(fioIn2,&iOrder2,r2,v2,&fMass2,&fSoft2,&fPot2,&fRho2,&u2,
		       &fMetals2,afOtherData2);
            break;
        case FIO_SPECIES_DARK:
            fioReadDark(fioIn1,&iOrder1,r1,v1,&fMass1,&fSoft1,&fPot1,&fRho1);
            fioReadDark(fioIn2,&iOrder2,r2,v2,&fMass2,&fSoft2,&fPot2,&fRho2);
            break;
        case FIO_SPECIES_STAR:
            fioReadStar(fioIn1,&iOrder1,r1,v1,&fMass1,&fSoft1,&fPot1,&fRho1,
			&fMetals1,&fTimer1,afOtherData1);
            fioReadStar(fioIn2,&iOrder2,r2,v2,&fMass2,&fSoft2,&fPot2,&fRho2,
			&fMetals2,&fTimer2,afOtherData2);
            break;
        default:
            fprintf(stderr,"Unsupported particle type: %d\n",eSpecies);
            abort();
            }
	if (bPosition) v = mag(r1,r2,1);
	else v = mag(v1,v2,0);

	if (v>vMax) vMax = v;
	if (v<vMin) vMin = v;

	printf("%g\n",v);
	}

    fprintf(stderr,"Min: %g, Max: %g\n", vMin, vMax);

    fioClose(fioIn2);
    fioClose(fioIn1);

    return 0;
    }
