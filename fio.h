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
**   fio = fioOpen("test.std",0);
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
**   fio = fioOpen("test.std",0);
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

/*
** These are the valid file formats
*/
typedef enum {
    FIO_FORMAT_TIPSY,
    FIO_FORMAT_HDF5,
    FIO_FORMAT_GRAFIC
    } FIO_FORMAT;

typedef enum {
    FIO_MODE_READING,
    FIO_MODE_WRITING
    } FIO_MODE;

/*
** These are the valid data types for attributes.
*/
typedef enum {
    FIO_TYPE_FLOAT=0,
    FIO_TYPE_DOUBLE
    } FIO_TYPE;

/*
** Particle species are a special thing.  Codes must know explicitly what data
** each particle will contain, so we can enumerate them here.
*/
typedef enum {
    FIO_SPECIES_ALL=0,
    FIO_SPECIES_DARK,
    FIO_SPECIES_SPH,
    FIO_SPECIES_STAR,
    FIO_SPECIES_LAST /* Must be last */
    } FIO_SPECIES;

/* This structure should be treated as PRIVATE.  Call the "fio" routines. */
typedef struct fioInfo {
    FIO_FORMAT eFormat;
    FIO_MODE   eMode;
    uint64_t nSpecies[FIO_SPECIES_LAST];
    void (*fcnClose)(struct fioInfo *fio);
    int  (*fcnSeek) (struct fioInfo *fio,uint64_t iPart,FIO_SPECIES eSpecies);
    FIO_SPECIES (*fcnSpecies) (struct fioInfo *fio);
    int (*fcnOpenNext)(struct fioInfo *fio, const char *fileName);

    int  (*fcnReadDark) (struct fioInfo *fio,
	uint64_t *piOrder,double *pdPos,double *pdVel,
	float *pfMass,float *pfSoft,float *pfPot);
    int  (*fcnReadSph) (
	struct fioInfo *fio,uint64_t *piOrder,double *pdPos,double *pdVel,
	float *pfMass,float *pfSoft,float *pfPot,
	float *pfRho,float *pfTemp, float *pfMetals);
    int  (*fcnReadStar) (struct fioInfo *fio,
	uint64_t *piOrder,double *pdPos,double *pdVel,
	float *pfMass,float *pfSoft,float *pfPot,
	float *pfMetals, float *pfTform);

    int  (*fcnWriteDark) (struct fioInfo *fio,
	uint64_t iOrder,const double *pdPos,const double *pdVel,
	float fMass,float fSoft,float fPot);
    int  (*fcnWriteSph) (
	struct fioInfo *fio,uint64_t iOrder,const double *pdPos,const double *pdVel,
	float fMass,float fSoft,float fPot,
	float fRho,float fTemp,float fMetals);
    int  (*fcnWriteStar) (struct fioInfo *fio,
	uint64_t iOrder,const double *pdPos,const double *pdVel,
	float fMass,float fSoft,float fPot,
	float fMetals,float fTform);

    int  (*fcnGetAttr)(struct fioInfo *fio,
	const char *attr, FIO_TYPE dataType, void *data);
    } *FIO;

/******************************************************************************\
** Generic Routines
\******************************************************************************/

/*
** Auto-detects the file format by looking at header information.
** bDouble is required for Tipsy; there is no other way to reliably tell.
*/
FIO fioOpen(const char *fileName,int bDouble);

/*
** Open the next file in a sequence of files.
*/
static inline int fioOpenNext(struct fioInfo *fio,const char *fileName) {
    return (*fio->fcnOpenNext)(fio,fileName);
    }

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
static inline uint64_t fioGetN(FIO fio,FIO_SPECIES eSpecies) {
    assert(eSpecies>=FIO_SPECIES_ALL && eSpecies<FIO_SPECIES_LAST);
    return fio->nSpecies[eSpecies];
    }

/*
** Seek to the N'th particle or the N'th particle of a given species.
*/
static inline int fioSeek(struct fioInfo *fio,uint64_t iPart,FIO_SPECIES eSpecies) {
    return (*fio->fcnSeek)(fio,iPart,eSpecies);
    }

/*
** Returns the species at the current file position (what will next be read)
*/
static inline FIO_SPECIES fioSpecies(struct fioInfo *fio) {
    return (*fio->fcnSpecies)(fio);
    }

/*
** Read a particle.  Must already be positioned at the appropriate particle.
** Normally fioSpecies() is called to determine which type to read, or a
** specific seek is performed to start reading a particular type.
*/
static inline int fioReadDark(
    FIO fio,uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft,float *pfPot) {
    return (*fio->fcnReadDark)(fio,piOrder,pdPos,pdVel,pfMass,pfSoft,pfPot);
    }
static inline int  fioReadSph(
    FIO fio,uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft,float *pfPot,
    float *pfRho,float *pfTemp,float *pfMetals) {
    return (*fio->fcnReadSph)(fio,piOrder,pdPos,pdVel,pfMass,pfSoft,pfPot,
			      pfRho,pfTemp,pfMetals);
    }
static inline int fioReadStar(
    FIO fio,uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft,float *pfPot,
    float *pfMetals,float *pfTform) {
    return (*fio->fcnReadStar)(fio,piOrder,pdPos,pdVel,pfMass,pfSoft,pfPot,
			       pfMetals,pfTform);
    }
/*
** Write a particle.  Must already be positioned at the appropriate particle.
*/
static inline int fioWriteDark(
    FIO fio,uint64_t iOrder,const double *pdPos,const double *pdVel,
    float fMass,float fSoft,float fPot) {
    return (*fio->fcnWriteDark)(fio,iOrder,pdPos,pdVel,fMass,fSoft,fPot);
    }
static inline int  fioWriteSph(
    FIO fio,uint64_t iOrder,const double *pdPos,const double *pdVel,
    float fMass,float fSoft,float fPot,
    float fRho,float fTemp,float fMetals) {
    return (*fio->fcnWriteSph)(fio,iOrder,pdPos,pdVel,fMass,fSoft,fPot,
			      fRho,fTemp,fMetals);
    }
static inline int fioWriteStar(
    FIO fio,uint64_t iOrder,const double *pdPos,const double *pdVel,
    float fMass,float fSoft,float fPot,
    float fMetals,float fTform) {
    return (*fio->fcnWriteStar)(fio,iOrder,pdPos,pdVel,fMass,fSoft,fPot,
			       fMetals,fTform);
    }
/*
** Returns the value of a given attribute.  Only "dTime" is available for
** Tipsy files, but HDF5 supports the inclusion of any arbitary attribute.
*/
static inline int fioGetAttr(FIO fio,
    const char *attr, FIO_TYPE dataType, void *data) {
    return (*fio->fcnGetAttr)(fio,attr,dataType,data);
    }

/******************************************************************************\
** TIPSY FORMAT
\******************************************************************************/

/*
** Open a Tipsy file.  The "fioOpen" routine can also be called and Tipsy format
** will be auto-detected.  Note that Standard/Native is detected automatically.
*/
FIO fioTipsyOpen(const char *fileName,int bDouble);

/*
** Opens a Tipsy file fragment (a file that has been split ala pkdgrav2).  The
** TOTAL number of SPH, dark and star particles must be passed to this routine
** which can be determined by using fioTipsyOpen on the first fragment.  The
** iStart parameter is the starting particle number of this fragment.  Normally
** this is determined by opening the fragments in order and getting the number
** of particles in the file by calling fioGetN and keeping a running total.
*/
FIO fioTipsyOpenPart(const char *fileName,int bDouble,int bStandard,
		     uint64_t nSph, uint64_t nDark, uint64_t nStar,
		     uint64_t iStart);

FIO fioTipsyCreate(const char *fileName,int bDouble,int bStandard,
		   double dTime,uint64_t nSph, uint64_t nDark, uint64_t nStar);
FIO fioTipsyAppend(const char *fileName,int bDouble,int bStandard);
FIO fioTipsyCreatePart(const char *fileName,int bAppend,int bDouble,int bStandard,
		       double dTime, uint64_t nSph, uint64_t nDark, uint64_t nStar,
		       uint64_t iStart);


int fioTipsyIsDouble(FIO fio);
int fioTipsyIsStandard(FIO fio);

/******************************************************************************\
** HDF5 FORMAT
\******************************************************************************/

/*
** Open an HDF5 file.  The "fioOpen" routine can also be called and HDF5 format
** will be auto-detected.  Fragments need no special treatment as the starting
** iOrder is saved as part of the metadata.
*/
FIO fioHDF5Open(const char *fileName);

/*
** Create an HDF5 file.
*/
FIO fioHDF5Create(const char *fileName,int bDouble);

/******************************************************************************\
** GRAFIC FORMAT
\******************************************************************************/

/*
** Open a GRAFIC format initial condition file.  This can be detected automatically
** by fioOpen because the "file name" is the directory containing ic_vel[cb][xyz].
*/
FIO fioGraficOpen(const char *dirName,int bDouble);

#endif
