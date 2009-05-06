#ifndef FIO_H
#define FIO_H

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

typedef struct fioInfo {
    uint64_t nSpecies[FIO_SPECIES_LAST+1];
    uint64_t oSpecies[FIO_SPECIES_LAST+1];
    double   dTime;
    void (*fcnClose)(struct fioInfo *fio);
    int  (*fcnSeek) (struct fioInfo *fio,uint64_t iPart,FIO_SPECIES eSpecies);
    int  (*fcnReadDark) (struct fioInfo *fio,
	uint64_t *piOrder,double *pdPos,double *pdVel,
	float *pfMass,float *pfSoft,float *pfPot);
    int  (*fcnGetAttr)(struct fioInfo *fio,
	const char *attr, FIO_TYPE dataType, void *data);
    } *FIO;

/******************************************************************************\
** Generic Routines
\******************************************************************************/

FIO fioOpen(const char *fileName,int bDouble);

static inline int fioGetN(FIO fio,FIO_SPECIES eSpecies) {
    assert(eSpecies>=FIO_SPECIES_ALL && eSpecies<FIO_SPECIES_LAST);
    return fio->nSpecies[eSpecies];
    }

static inline int fioSeek(struct fioInfo *fio,uint64_t iPart,FIO_SPECIES eSpecies) {
    return (*fio->fcnSeek)(fio,iPart,eSpecies);
    }

static inline void fioClose(struct fioInfo *fio) {
    (*fio->fcnClose)(fio);
    }

static inline int fioReadDark(FIO fio,
	     uint64_t *piOrder,double *pdPos,double *pdVel,
	     float *pfMass,float *pfSoft,float *pfPot) {
    return (*fio->fcnReadDark)(fio,piOrder,pdPos,pdVel,pfMass,pfSoft,pfPot);
    }

static inline int fioGetAttr(FIO fio,
    const char *attr, FIO_TYPE dataType, void *data) {
    return (*fio->fcnGetAttr)(fio,attr,dataType,data);
    }

/******************************************************************************\
** TIPSY FORMAT
\******************************************************************************/

FIO fioTipsyOpen(const char *fileName,int bDouble);

/******************************************************************************\
** GRAFIC FORMAT
\******************************************************************************/

FIO fioGraficOpen(const char *dirName,int bDouble);

#endif
