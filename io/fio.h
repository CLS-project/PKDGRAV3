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

/******************************************************************************\
** File I/O (fio) Routines
**
** This module abstracts access to different format files such that a code can
** easily read all supported formats.
**
** To read a file, generally the generic "fioOpen" routine is called and the
** format is automatically detected by examining the file header.
**
** The following routines are most useful.
**
**   fioOpen     - Open the specified file (autodetects format)
**   fioClose    - Closes an open file.
**   fioGetN     - Total number of particles (or # of a specific species)
**   fioSeek     - Advance in the file to the specified particle
**   fioSpecies  - Returns the species of the next particle to read
**   fioReadDark - Reads a dark particle
**   fioReadSph  - Reads an SPH particle
**   fioReadStar - Reads a star particle
**
** Example (sequential read of all particles):
**
**   FIO fio;
**   uint64_t N, i;
**
**   fio = fioOpen("test.std");
**   if (fio==NULL) ...
**
**   N = fioGetN(fio,FIO_SPECIES_ALL);
**   for( i=0; i<N; i++ ) {
**       switch(fioSpecies(fio)) {
**       case FIO_SPECIES_DARK:
**           fioReadDark(fio,...);
**           break;
**       case FIO_SPECIES_SPH:
**           fioReadSph(fio,...);
**           break;
**       case FIO_SPECIES_STAR:
**           fioReadStar(fio,...);
**           break;
**           }
**       default:
**           perror("invalid/unknown particle type");
**           }
**       }
**
**  fioClose(fio);
**
** Example (read all dark matter particles)
**
**   FIO fio;
**   uint64_t nDark, i;
**
**   fio = fioOpen("test.std");
**   if (fio==NULL) ...
**
**   nDark = fioGetN(fio,FIO_SPECIES_DARK);
**   fioSeek(fio,0,FIO_SPECIES_DARK);
**   for( i=0; i<nDark; i++ ) {
**       fioReadDark(fio,...);
**       }
**
**   fioClose(fio);
**
\******************************************************************************/
#ifndef FIO_H
#define FIO_H

#include <stdint.h>
#include <assert.h>

/*
** These are the valid file formats
*/
typedef enum {
    FIO_FORMAT_MULTIPLE,
    FIO_FORMAT_TIPSY,
    FIO_FORMAT_HDF5,
    FIO_FORMAT_GRAFIC,
    FIO_FORMAT_GADGET2
} FIO_FORMAT;

/*
** Here are the valid flags for Create
*/
#define FIO_FLAG_DOUBLE_POS    1
#define FIO_FLAG_DOUBLE_VEL    2
#define FIO_FLAG_COMPRESS_MASS 4
#define FIO_FLAG_COMPRESS_SOFT 8
#define FIO_FLAG_CHECKPOINT    16 /* Restart - normally double */
#define FIO_FLAG_POTENTIAL     32 /* Include the potential */
#define FIO_FLAG_DENSITY       64 /* Include the density */
#define FIO_FLAG_ID           128

typedef enum {
    FIO_MODE_READING,
    FIO_MODE_WRITING
} FIO_MODE;

/*
** These are the valid data types for attributes.
*/
typedef enum {
    FIO_TYPE_FLOAT=0,
    FIO_TYPE_DOUBLE,
    FIO_TYPE_UINT32,
    FIO_TYPE_UINT64,
    FIO_TYPE_UINT8,
    FIO_TYPE_INT,
    FIO_TYPE_STRING,
} FIO_TYPE;

/*
** Particle species are a special thing.  Codes must know explicitly what data
** each particle will contain, so we can enumerate them here.
*/
enum FIO_SPECIES {
    FIO_SPECIES_DARK,
    FIO_SPECIES_SPH,
    FIO_SPECIES_STAR,
    FIO_SPECIES_BH,
    FIO_SPECIES_UNKNOWN,    /* Deleted for example */
    FIO_SPECIES_LAST,       /* This is the count of valid species (all is added on the end) */
    FIO_SPECIES_ALL=FIO_SPECIES_LAST /* Internally keep an ALL count */
};

typedef uint64_t fioSpeciesList[FIO_SPECIES_LAST+1];

typedef struct {
    uint64_t       iFirst;      /* Starting particle index */
    char          *pszFilename; /* Filename of this file */
    fioSpeciesList nSpecies;    /* # of each species in this file */
} fioFileInfo;

typedef struct {
    int iFile;                /* Current file */
    int nFiles;               /* Total number of files */
    fioFileInfo *fileInfo;    /* Array of information for each file */
} fioFileList;

/* This structure should be treated as PRIVATE.  Call the "fio" routines. */
typedef struct fioInfo {
    FIO_FORMAT eFormat;
    FIO_MODE   eMode;
    int        mFlags;
    fioSpeciesList nSpecies;

    /* This is for multi-file support */
    fioFileList fileList;

    void (*fcnClose)(struct fioInfo *fio);
    int  (*fcnSeek) (struct fioInfo *fio,uint64_t iPart,enum FIO_SPECIES eSpecies);
    enum FIO_SPECIES (*fcnSpecies) (struct fioInfo *fio);

    int  (*fcnReadDark) (struct fioInfo *fio,
                         uint64_t *piParticleID,double *pdPos,double *pdVel,
                         float *pfMass,float *pfSoft,float *pfPot,float *pfDen);
    int  (*fcnReadSph) (struct fioInfo *fio,
                        uint64_t *piParticleID,double *pdPos,double *pdVel,
                        float *pfMass,float *pfSoft,float *pfPot,float *pfDen,
                        float *pfTemp,float *pfMetals,float *pfOtherData);
    int  (*fcnReadStar) (struct fioInfo *fio,
                         uint64_t *piParticleID,double *pdPos,double *pdVel,
                         float *pfMass,float *pfSoft,float *pfPot,float *pfDen,
                         float *pfMetals,float *pfTform,float *pfOtherData);
    /* IA: We left some dummy arguments, that later can be renamed to actual IO variables interesting
     *  for black holes */
    int  (*fcnReadBH) (struct fioInfo *fio,
                       uint64_t *piParticleID,double *pdPos,double *pdVel,
                       float *pfMass,float *pfSoft,float *pfPot,float *pfDen,
                       float *pfOtherData,float *pfTform);

    int  (*fcnWriteDark) (struct fioInfo *fio,
                          uint64_t iParticleID,const double *pdPos,const double *pdVel,
                          float fMass,float fSoft,float fPot,float fDen,float *pfOtherData);
    int  (*fcnWriteSph) (struct fioInfo *fio,
                         uint64_t iParticleID,const double *pdPos,const double *pdVel,
                         float fMass,float fSoft,float fPot,float fDen,
                         float fTemp,float *pfMetals,float fBall,float fIntEnergy,float *pfOtherData);
    int  (*fcnWriteStar) (struct fioInfo *fio,
                          uint64_t iParticleID,const double *pdPos,const double *pdVel,
                          float fMass,float fSoft,float fPot,float fDen,
                          float *pfMetals,float *pfOtherData);
    int  (*fcnWriteBH) (struct fioInfo *fio,
                        uint64_t iParticleID,const double *pdPos,const double *pdVel,
                        float fMass,float fSoft,float fPot,float fDen,
                        float *pfOtherData,float fTform);

    int  (*fcnGetAttr)(struct fioInfo *fio, const int headerType,
                       const char *attr, FIO_TYPE dataType, void *data);
    int  (*fcnSetAttr)(struct fioInfo *fio, const int headerType,
                       const char *attr, FIO_TYPE dataType, int size, void *data);
} *FIO;

/******************************************************************************\
** Generic Routines
\******************************************************************************/
#ifdef __cplusplus
extern "C" {
#endif

/*
** Auto-detects the file format by looking at header information.
*/
FIO fioOpen(const char *fileName,double dOmega0,double dOmegab);
FIO fioOpenMany(int nFiles, const char *const *fileNames,double dOmega0,double dOmegab);
size_t fioDump(FIO fio, size_t nBytes, void *pBuffer);
FIO fioLoad(void *pBuffer,double dOmega0,double dOmegab);

/*
** Close an open file of any format.
*/
static inline void fioClose(struct fioInfo *fio) {
    (*fio->fcnClose)(fio);
}

/*
** Return the format of the currently open file.
*/
static inline FIO_FORMAT fioFormat(FIO fio) {
    return fio->eFormat;
}

/*
** Return the mode of the currently open file.
*/
static inline FIO_MODE fioMode(FIO fio) {
    return fio->eMode;
}

/*
** Returns the number of particles of a given species (or ALL).
*/
static inline uint64_t fioGetN(FIO fio,enum FIO_SPECIES eSpecies) {
    assert(eSpecies<=FIO_SPECIES_LAST);
    return fio->nSpecies[eSpecies];
}

/*
** Seek to the N'th particle or the N'th particle of a given species.
*/
static inline int fioSeek(struct fioInfo *fio,uint64_t iPart,enum FIO_SPECIES eSpecies) {
    return (*fio->fcnSeek)(fio,iPart,eSpecies);
}

/*
** Returns the species at the current file position (what will next be read)
*/
static inline enum FIO_SPECIES fioSpecies(struct fioInfo *fio) {
    return (*fio->fcnSpecies)(fio);
}

/*
** Read a particle.  Must already be positioned at the appropriate particle.
** Normally fioSpecies() is called to determine which type to read, or a
** specific seek is performed to start reading a particular type.
*/
static inline int fioReadDark(
    FIO fio,uint64_t *piParticleID,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft,float *pfPot,float *pfDen) {
    return (*fio->fcnReadDark)(fio,piParticleID,pdPos,pdVel,pfMass,pfSoft,pfPot,pfDen);
}
static inline int  fioReadSph(
    FIO fio,uint64_t *piParticleID,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft,float *pfPot,float *pfDen,
    float *pfTemp,float *pfMetals,float *pfOtherData) {
    return (*fio->fcnReadSph)(fio,piParticleID,pdPos,pdVel,pfMass,pfSoft,pfPot,pfDen,
                              pfTemp,pfMetals,pfOtherData);
}
static inline int fioReadStar(
    FIO fio,uint64_t *piParticleID,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft,float *pfPot,float *pfDen,
    float *pfMetals,float *pfTform,float *pfOtherData) {
    return (*fio->fcnReadStar)(fio,piParticleID,pdPos,pdVel,pfMass,pfSoft,pfPot,pfDen,
                               pfMetals,pfTform,pfOtherData);
}
static inline int fioReadBH(
    FIO fio,uint64_t *piParticleID,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft,float *pfPot,float *pfDen,
    float *pfOtherData,float *pfTform) {
    return (*fio->fcnReadBH)(fio,piParticleID,pdPos,pdVel,pfMass,pfSoft,pfPot,pfDen,
                             pfOtherData,pfTform);
}
/*
** Write a particle.  Must already be positioned at the appropriate particle.
*/
static inline int fioWriteDark(
    FIO fio,uint64_t iParticleID,const double *pdPos,const double *pdVel,
    float fMass,float fSoft,float fPot,float fDen,float *pfOtherData) {
    return (*fio->fcnWriteDark)(fio,iParticleID,pdPos,pdVel,fMass,fSoft,fPot,fDen,pfOtherData);
}
static inline int  fioWriteSph(
    FIO fio,uint64_t iParticleID,const double *pdPos,const double *pdVel,
    float fMass,float fSoft,float fPot,float fDen,
    float fTemp,float *pfMetals,float fBall,float fIntEnergy,float *pfOtherData) {
    return (*fio->fcnWriteSph)(fio,iParticleID,pdPos,pdVel,fMass,fSoft,fPot,fDen,
                               fTemp,pfMetals,fBall,fIntEnergy,pfOtherData);
}
static inline int fioWriteStar(
    FIO fio,uint64_t iParticleID,const double *pdPos,const double *pdVel,
    float fMass,float fSoft,float fPot,float fDen,
    float *pfMetals,float *pfOtherData) {
    return (*fio->fcnWriteStar)(fio,iParticleID,pdPos,pdVel,fMass,fSoft,fPot,fDen,
                                pfMetals,pfOtherData);
}
static inline int fioWriteBH(
    FIO fio,uint64_t iParticleID,const double *pdPos,const double *pdVel,
    float fMass,float fSoft,float fPot,float fDen,
    float *pfOtherData,float fTform) {
    return (*fio->fcnWriteBH)(fio,iParticleID,pdPos,pdVel,fMass,fSoft,fPot,fDen,
                              pfOtherData,fTform);
}
/*
** Returns the value of a given attribute.  Only "dTime" is available for
** Tipsy files, but HDF5 supports the inclusion of any arbitary attribute.
*/
static inline int fioGetAttr(FIO fio, int headerType,
                             const char *attr, FIO_TYPE dataType, void *data) {
    return (*fio->fcnGetAttr)(fio, headerType, attr,dataType,data);
}

/*
** Sets an arbitrary attribute.  Only supported for HDF5; other formats
** return 0 indicating that it was not successful.
*/
static inline int fioSetAttr(FIO fio, int headerType,
                             const char *attr, FIO_TYPE dataType, int size, void *data) {
    return (*fio->fcnSetAttr)(fio, headerType, attr,dataType,size,data);
}

static inline int fioGetFlags(FIO fio) {
    return fio->mFlags;
}

/******************************************************************************\
** TIPSY FORMAT
\******************************************************************************/

FIO fioTipsyCreate(const char *fileName,int mFlags,int bStandard,
                   double dTime,uint64_t nSph, uint64_t nDark, uint64_t nStar);
FIO fioTipsyAppend(const char *fileName,int mFlags,int bStandard);
FIO fioTipsyCreatePart(const char *fileName,int bAppend,int mFlags,int bStandard,
                       double dTime, uint64_t nSph, uint64_t nDark, uint64_t nStar,
                       uint64_t iStart);
int fioTipsyIsDouble(FIO fio);
int fioTipsyIsStandard(FIO fio);

/******************************************************************************\
** HDF5 FORMAT
\******************************************************************************/

/*
** Create an HDF5 file.
*/
#define HDF5_HEADER_G 0
#define HDF5_COSMO_G  1
#define HDF5_UNITS_G  2
#define HDF5_PARAM_G  3
FIO fioHDF5Create(const char *fileName,int mFlags);

/******************************************************************************\
** GADGET2 FORMAT
\******************************************************************************/

FIO fioGadgetCreate(
    const char *fileName,int mFlags, double dTime, double Lbox,
    double Omega0, double OmegaLambda, double HubbleParam,
    int nTypes, const uint64_t *nPart,
    int nFiles, const uint64_t *nAll,
    const double *dMass );

/******************************************************************************\
** GRAFIC FORMAT
\******************************************************************************/
#ifdef __cplusplus
}
#endif

#endif
