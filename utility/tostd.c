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
#include "io/fio.h"

#define BUFFER_SIZE (8*1024*1024)

#define OPT_DOUBLE    'd'
#define OPT_NATIVE    'n'
#define OPT_HDF5      '5'
#define OPT_POTENTIAL 'p'
#define OPT_DENSITY   'r'
#define OPT_GADGET2   'g'

#define OPT_OMEGA0    'm'
#define OPT_LAMBDA    '^'
#define OPT_H0        'H'
#define OPT_LBOX      'L'

#define OPT_PERIOD    'P'

static inline void wrap(double *r,double dPeriod) {
    if (dPeriod > 0.0) {
        double dBound = dPeriod * 0.5;
        int i;
        for (i=0; i<3; ++i) {
            if (r[i] >= dBound) r[i] -= dPeriod;
            else if (r[i] < -dBound) r[i] += dPeriod;
        }
    }
}

int main( int argc, char *argv[] ) {
    int bError = 0;
    int bDouble = 0;
    int bNative = 0;
    int bHDF5 = 0;
    int bGADGET2 = 0;
    int bPotential = 0;
    int bDensity = 0;
    double Omega0 = 0.0, OmegaLambda = 0.0, HubbleParam = 0.0, Lbox=0.0;
    double dPeriod = 0.0;
    int bSetOmega0=0, bSetOmegaLambda=0, bSetHubbleParam=0;
    uint64_t N, nSph, nDark, nStar, i;
    double dTime;
    uint64_t iOrder;
    double r[3], v[3];
    float fMass, fSoft, fPot, fRho, u, fMetals, fTimer;
    uint64_t nPart[6];
    double dMass[6];

    FIO fioIn, fioOut;
    enum FIO_SPECIES eSpecies;
    int inNameIndex = 0;
    int inNameCount = 0;
    const char *outName = 0;

    //! Parse command line arguments (flags).
    for (;;) {
        int c, option_index=0;

        static struct option long_options[] = {
            { "double",       0, 0, OPT_DOUBLE},
            { "native",       0, 0, OPT_NATIVE},
            { "hdf5",         0, 0, OPT_HDF5},
            { "gadget2",      0, 0, OPT_GADGET2},
            { "potential",    0, 0, OPT_POTENTIAL},
            { "density",    0, 0, OPT_DENSITY},
            { "omega-matter", 1, 0, OPT_OMEGA0 },
            { "omega0",       1, 0, OPT_OMEGA0 },
            { "omega-lambda", 1, 0, OPT_LAMBDA },
            { "hubble0",      1, 0, OPT_H0 },
            { "H0",           1, 0, OPT_H0 },
            { "box-width",    1, 0, OPT_LBOX },
            { "period",       1, 0, OPT_PERIOD },
            { NULL,   0, 0, 0 },
        };

        c = getopt_long( argc, argv, "dn5gprm:^:H:L:P:",
                         long_options, &option_index );
        if ( c == -1 ) break;

        switch (c) {
        case OPT_DOUBLE:
            bDouble = 1;
            break;
        case OPT_NATIVE:
            bNative = 1;
            break;
        case OPT_HDF5:
            bHDF5 = 1;
            break;
        case OPT_GADGET2:
            bGADGET2 = 1;
            break;
        case OPT_POTENTIAL:
            bPotential = 1;
            break;
        case OPT_DENSITY:
            bDensity = 1;
            break;
        case OPT_OMEGA0:
            assert( optarg != 0 );
            bSetOmega0 = 1;
            Omega0 = atof(optarg);
            break;
        case OPT_LAMBDA:
            assert( optarg != 0 );
            bSetOmegaLambda=1;
            OmegaLambda = atof(optarg);
            break;
        case OPT_H0:
            assert( optarg != 0 );
            bSetHubbleParam = 1;
            HubbleParam = atof(optarg);
            break;
        case OPT_LBOX:
            assert( optarg != 0 );
            Lbox = atof(optarg);
            break;
        case OPT_PERIOD:
            assert( optarg != 0 );
            dPeriod = atof(optarg);
            break;
        default:
            bError = 1;
            break;
        }
    }

    if (bNative + bHDF5 + bGADGET2 > 1) {
        fprintf(stderr, "Specify only one of --hdf5, --gadget2 or --native\n" );
        bError = 1;
    }
    else if (bGADGET2 && !(bSetOmega0 && bSetOmegaLambda && bSetHubbleParam) ) {
        fprintf(stderr, "Specify Omega0, OmegaLambda and H0 with GADGET2 output format\n" );
        bError = 1;
    }
#ifndef USE_HDF5
    if (bHDF5) {
        fprintf(stderr, "HDF5 support was not compiled in.\n" );
        bError = 1;
    }
#endif

    if ( optind < argc ) {
        inNameIndex = optind++;
        inNameCount = argc - optind;
    }
    else {
        fprintf(stderr, "Missing input file(s)\n" );
        bError = 1;
    }

    if ( optind < argc )
        outName = argv[argc-1];
    else {
        fprintf(stderr, "Missing Tipsy output file\n" );
        bError = 1;
    }

    if ( bError ) {
        fprintf(stderr, "Usage: %s [-p] <input...> <outtipsy>\n"
                "  -d,--double    Output double precision positions\n"
                "  -n,--native    Output a native tipsy binary\n"
                "  -g,--gadget2   Output a GADGET2 binary\n"
                "  -m,--omega0=<omega-matter>\n"
                "  -v,--omega-lambda=<omega-lambda>\n"
                "  -H,--H0=<hubble-parameter>\n"
#ifdef USE_HDF5
                "  -5,--hdf5      Output in HDF5 format\n"
                "  -p,--potential Included potentials in HDF5 output\n"
#endif
                ,argv[0] );
        exit(1);
    }

    fioIn = fioOpenMany(inNameCount,(const char *const *)&argv[inNameIndex],0.0,0.0);
    if (fioIn==NULL) {
        perror(argv[inNameIndex]);
        exit(errno);
    }
    N     = fioGetN(fioIn,FIO_SPECIES_ALL);
    nSph  = fioGetN(fioIn,FIO_SPECIES_SPH);
    nDark = fioGetN(fioIn,FIO_SPECIES_DARK);
    nStar = fioGetN(fioIn,FIO_SPECIES_STAR);
    if (!fioGetAttr(fioIn,0,"Time",FIO_TYPE_DOUBLE,&dTime)) dTime = 0.0;

    printf("dTime=%g\n",dTime);

#ifdef USE_HDF5
    if (bHDF5) {
        int mFlag = FIO_FLAG_COMPRESS_MASS | FIO_FLAG_COMPRESS_SOFT;
        if (bDouble) mFlag |= FIO_FLAG_DOUBLE_POS | FIO_FLAG_DOUBLE_VEL;
        if (bPotential) mFlag |= FIO_FLAG_POTENTIAL;
        if (bDensity) mFlag |= FIO_FLAG_DENSITY;
        fioOut = fioHDF5Create(outName,mFlag);
    }
#else
    if (0) {}
#endif
    else if (bGADGET2) {
        int mFlag = 0;
        if (bDouble) mFlag |= FIO_FLAG_DOUBLE_POS | FIO_FLAG_DOUBLE_VEL;
        nPart[0] = nSph;
        nPart[1] = nDark;
        nPart[2] = 0;
        nPart[3] = 0;
        nPart[4] = nStar;
        nPart[5] = 0;
        dMass[0] = dMass[1] = dMass[2] = dMass[3] = dMass[4] = dMass[5] = 0.0;
        fioOut = fioGadgetCreate(outName,mFlag,dTime,Lbox,
                                 Omega0, OmegaLambda, HubbleParam,
                                 6, nPart, 1, nPart, dMass );
    }
    else {
        fioOut = fioTipsyCreate(outName,bDouble,!bNative,dTime,nSph,nDark,nStar);
    }

    if (fioOut==NULL) {
        perror(outName);
        exit(errno);
    }
    fioSetAttr(fioOut,0,"Time",FIO_TYPE_DOUBLE,1,&dTime);

    float otherData[10];
    for ( i=0; i<N; i++ ) {
        eSpecies = fioSpecies(fioIn);
        switch (eSpecies) {
        case FIO_SPECIES_SPH:
            fioReadSph(fioIn,&iOrder,r,v,&fMass,&fSoft,&fPot,&fRho,&u,&fMetals,otherData);
            wrap(r,dPeriod);
            fioWriteSph(fioOut,iOrder,r,v,fMass,fSoft,fPot,fRho,u,&fMetals,0.0,0.0,otherData);
            break;
        case FIO_SPECIES_DARK:
            fioReadDark(fioIn,&iOrder,r,v,&fMass,&fSoft,&fPot,&fRho);
            wrap(r,dPeriod);
            fioWriteDark(fioOut,iOrder,r,v,fMass,fSoft,fPot,fRho,otherData);
            break;
        case FIO_SPECIES_STAR:
            fioReadStar(fioIn,&iOrder,r,v,&fMass,&fSoft,&fPot,&fRho,&fMetals,otherData,otherData);
            wrap(r,dPeriod);
            fioWriteStar(fioOut,iOrder,r,v,fMass,fSoft,fPot,fRho,&fMetals,otherData);
            break;
        default:
            fprintf(stderr,"Unsupported particle type: %d\n",eSpecies);
            abort();
        }
    }

    fioClose(fioOut);
    fioClose(fioIn);

    return 0;
}
