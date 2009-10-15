#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <errno.h>
#include <sys/stat.h>
#include <rpc/xdr.h>
#if defined(HAVE_WORDEXP) && defined(HAVE_WORDFREE)
#include <wordexp.h>
#elif defined(HAVE_GLOB) && defined(HAVE_GLOBFREE)
#include <glob.h>
#endif

#include "fio.h"

const char *fio_module_id = "$Id$";
const char *fio_h_module_id = FIO_H_MODULE_ID;

/*
** This uses the best available seek routine to move to the specified
** 64-bit offset into a file.  The new "fseeko" function is preferred,
** but we fallback to multiple fseek calls if required.
*/
static int safe_fseek(FILE *fp,uint64_t lStart) {
#ifdef HAVE_FSEEKO
    uint64_t lSeeked;
    int iErr;
    errno = 0;
    iErr = fseeko(fp,lStart,SEEK_SET);
    if (iErr) {
	perror("safe_fseek failed");
	exit(errno);
	}
    lSeeked = ftello(fp);
    if(lSeeked != lStart) {
	fprintf(stderr,"fseeko() borked: %lu != %lu\n", lStart, lSeeked);
	exit(ESPIPE);
	}
    return 0;
#else
    /* fseek fails for offsets >= 2**31; this is an ugly workaround */
    static const off_t MAX_OFFSET = 2147483640;
    int iErr;

    if (lStart > MAX_OFFSET) {
	errno = 0;
	rewind(fp);
	if (errno) {
	    perror("safe_fseek failed");
	    exit(errno);
	    }
	while (lStart > MAX_OFFSET) {
	    fseek(fp,MAX_OFFSET,SEEK_CUR);
	    lStart -= MAX_OFFSET;
	    }
	iErr = fseek(fp,lStart,SEEK_CUR);
	if (iErr) {
	    perror("safe_fseek failed");
	    exit(errno);
	    }
	}
    else {
	iErr = fseek(fp,lStart,SEEK_SET);
	if (iErr) {
	    perror("safe_fseek failed");
	    exit(errno);
	    }
	}
    return 0;
#endif
    }

/*
** This is used for the HDF5 group name.
*/
const char *fioSpeciesName(FIO_SPECIES eSpecies) {
    switch(eSpecies) {
    case FIO_SPECIES_DARK: return "dark";
    case FIO_SPECIES_SPH:  return "sph";
    case FIO_SPECIES_STAR: return "star";
    default:
	break;
	}
    return NULL;
    }

static void fioFree(FIO fio) {
    if (fio->iFirst) free(fio->iFirst);
    if (fio->pszFiles) {
	free(fio->pszFiles[0]);
	free(fio->pszFiles);
	}
    }

/*
** Given a list of one or more files, this function will expand any wildcards
** present (if possible) and return a complete list of matching files.
*/
static int fileScan( FIO fio, int nFiles, const char * const *szFilenames) {
    int i, nScan, iIdx, nSize;
    char *pszFilename;
#if defined(HAVE_WORDEXP) && defined(HAVE_WORDFREE)
    wordexp_t files;
    int flag = 0;
#elif defined(HAVE_GLOB) && defined(HAVE_GLOBFREE)
    glob_t files;
    int flag = GLOB_ERR;
#endif

    /*
    ** Run through the list of files and expand them producing a total
    ** count of files, and the list in whatever format we can.
    */
    nSize = 0;
#if defined(HAVE_WORDEXP) && defined(HAVE_WORDFREE)
    for( iIdx = 0; iIdx<nFiles; iIdx++ ) {
	wordexp(szFilenames[iIdx], &files, flag);
	if ( files.we_wordc <= 0 ) {
	    return iIdx+1;
	    }
	flag |= WRDE_APPEND;
	}
    nScan = files.we_wordc;
    for(i=0; i<nScan; i++)
	nSize += strlen(files.we_wordv[i])+1;
#elif defined(HAVE_GLOB) && defined(HAVE_GLOBFREE)
    for( iIdx = 0; iIdx<nFiles; iIdx++ ) {
	if (glob(szFilenames[iIdx],flag,NULL,&files) || files.gl_pathc==0) {
	    return iIdx+1;
	    }
	flag |= GLOB_APPEND;
	}
    nScan = files.gl_pathc;
    for(i=0; i<nScan; i++)
	nSize += strlen(files.gl_pathv[i])+1;
#else
    nScan = nFiles;
    nSize += strlen(szFilesnames[iIdx]) + 1;
#endif

    /*
    ** Allocate space for the list of files.  We use a single buffer for the
    ** names -- we have already calculated how big this needs to be.
    */
    fio->nFiles = nScan;
    fio->iFirst = malloc(sizeof(uint64_t)*(nScan+1));
    assert(fio->iFirst);
    fio->pszFiles = malloc(sizeof(char *)*nScan);
    assert(fio->pszFiles);
    pszFilename = malloc(nSize);
    assert(pszFilename);

    /*
    ** Great, now just copy the filenames into the new buffer
    */
    for( i=0; i<nScan; i++ ) {
	fio->pszFiles[i] = pszFilename;
#if defined(HAVE_WORDEXP) && defined(HAVE_WORDFREE)
	strcpy( pszFilename, files.we_wordv[i] );
#elif defined(HAVE_GLOB) && defined(HAVE_GLOBFREE)
	strcpy( pszFilename, files.gl_pathv[i] );
#else
	strcpy( pszFilename, szFilenames[i] );
#endif
	pszFilename += strlen(pszFilename) + 1;
	}

    /* Finished: free memory used by wildcard routines */
#if defined(HAVE_WORDEXP) && defined(HAVE_WORDFREE)
    wordfree(&files);
#elif defined(HAVE_GLOB) && defined(HAVE_GLOBFREE)
    globfree(&files);
#endif
    return 0;
    }

/******************************************************************************\
** Stubs that are called when a feature is not supported.
\******************************************************************************/

static int  fioNoReadDark(
    struct fioInfo *fio,uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft,float *pfPot) {
    fprintf(stderr,"Reading dark particles is not supported\n");
    abort();
    }

static int fioNoReadSph(
    struct fioInfo *fio,uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft, float *pfPot,
    float *pfRho,float *pfTemp, float *pfMetals) {
    fprintf(stderr,"Reading SPH particles is not supported\n");
    abort();
    }

static int fioNoReadStar(
    struct fioInfo *fio,uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft,float *pfPot,float *pfMetals, float *pfTform) {
    fprintf(stderr,"Reading star particles is not supported\n");
    abort();
    }

static int  fioNoWriteDark(
    struct fioInfo *fio,uint64_t iOrder,const double *pdPos,const double *pdVel,
    float fMass,float fSoft,float fPot) {
    fprintf(stderr,"Writing dark particles is not supported\n");
    abort();
    }

static int fioNoWriteSph(
    struct fioInfo *fio,uint64_t iOrder,const double *pdPos,const double *pdVel,
    float fMass,float fSoft,float fPot,
    float fRho,float fTemp,float fMetals) {
    fprintf(stderr,"Writing SPH particles is not supported\n");
    abort();
    }

static int fioNoWriteStar(
    struct fioInfo *fio,uint64_t iOrder,const double *pdPos,const double *pdVel,
    float fMass,float fSoft,float fPot,float fMetals,float fTform) {
    fprintf(stderr,"Writing star particles is not supported\n");
    abort();
    }

/******************************************************************************\
** TIPSY FORMAT
\******************************************************************************/

#define TIO_BUFFER_SIZE (1024*1024)

typedef struct {
    double dTime;
    unsigned nBodies;
    unsigned nDim;
    unsigned nSph;
    unsigned nDark;
    unsigned nStar;
    unsigned nPad;
    } tipsyHdr;

typedef struct {
    float mass;
    float pos[3];
    float vel[3];
    float eps;
    float phi;
    } tipsyDark;

typedef struct {
    float mass;
    float pos[3];
    float vel[3];
    float rho;
    float temp;
    float hsmooth;
    float metals;
    float phi;
    } tipsySph;

typedef struct {
    float mass;
    float pos[3];
    float vel[3];
    float metals;
    float tform;
    float eps;
    float phi;
    } tipsyStar;

#define TIPSY_FD_SIZE(t,d) (sizeof(float)*((t)-(d))+sizeof(double)*(d))
#define TIPSY_DARK_SIZE(dark,dbl) ((dark)*TIPSY_FD_SIZE(9,dbl))
#define TIPSY_SPH_SIZE(sph,dbl)   ((sph)*TIPSY_FD_SIZE(12,dbl))
#define TIPSY_STAR_SIZE(star,dbl) ((star)*TIPSY_FD_SIZE(11,dbl))
#define TIPSY_TOTAL_SIZE(dark,gas,star,dbl) sizeof(tipsyHdr)+TIPSY_DARK_SIZE(dark,dbl)+TIPSY_SPH_SIZE(gas,dbl)+TIPSY_STAR_SIZE(star,dbl)

typedef struct {
    struct fioInfo fio; /* "base class" */
    FILE *fp;
    int iFile;
    double dTime;
    char *fpBuffer;
    XDR xdr;
    off_t nHdrSize;     /* Size of header; 0 if second or subsequent file */
    uint64_t iStart;    /* Index of first particle in the file fragment */
    uint64_t iOrder;    /* Current particle index */
    int bDoublePos;
    int bDoubleVel;
    } fioTipsy;

static int xdrHeader(XDR *pxdr,tipsyHdr *ph) {
    if (!xdr_double(pxdr,&ph->dTime)) return 0;
    if (!xdr_u_int(pxdr,&ph->nBodies)) return 0;
    if (!xdr_u_int(pxdr,&ph->nDim)) return 0;
    if (!xdr_u_int(pxdr,&ph->nSph)) return 0;
    if (!xdr_u_int(pxdr,&ph->nDark)) return 0;
    if (!xdr_u_int(pxdr,&ph->nStar)) return 0;
    if (!xdr_u_int(pxdr,&ph->nPad)) return 0;
    return 1;
    }

static int tipsyGetAttr(FIO fio,
    const char *attr, FIO_TYPE dataType, void *data) {
    fioTipsy *tio = (fioTipsy *)fio;
    assert(fio->eFormat == FIO_FORMAT_TIPSY && fio->eMode==FIO_MODE_READING);
    if ( strcmp(attr,"dTime")==0 ) {
	if (!tio->nHdrSize) return 0;
	switch(dataType) {
	case FIO_TYPE_FLOAT: *(float *)(data) = tio->dTime; break;
	case FIO_TYPE_DOUBLE:*(double *)(data) = tio->dTime; break;
	default: return 0;
	    }
	return 1;
	}

    return 0;
    }

static FIO_SPECIES tipsySpecies(FIO fio) {
    fioTipsy *tio = (fioTipsy *)fio;
    assert(fio->eFormat == FIO_FORMAT_TIPSY && fio->eMode==FIO_MODE_READING);
    if (tio->iOrder<tio->fio.nSpecies[FIO_SPECIES_SPH]) return FIO_SPECIES_SPH;
    else if (tio->iOrder<tio->fio.nSpecies[FIO_SPECIES_SPH]+tio->fio.nSpecies[FIO_SPECIES_DARK])
	return FIO_SPECIES_DARK;
    else if (tio->iOrder<tio->fio.nSpecies[FIO_SPECIES_SPH]+tio->fio.nSpecies[FIO_SPECIES_DARK]+tio->fio.nSpecies[FIO_SPECIES_STAR])
	return FIO_SPECIES_STAR;
    else return FIO_SPECIES_LAST;
    }

static int tipsySwitchFile(FIO fio) {
    fioTipsy *tio = (fioTipsy *)fio;
    int bStandard = fioTipsyIsStandard(fio);
    if (bStandard) xdr_destroy(&tio->xdr);
    fclose(tio->fp);

    tio->fp = fopen(fio->pszFiles[tio->iFile],"r");
    if (tio->fp == NULL) printf("Fail: %d %s\n",tio->iFile,fio->pszFiles[tio->iFile]);
    if (tio->fp == NULL) return 1;
    if (tio->fpBuffer != NULL) setvbuf(tio->fp,tio->fpBuffer,_IOFBF,TIO_BUFFER_SIZE);
    if (bStandard) xdrstdio_create(&tio->xdr,tio->fp,XDR_DECODE);
    return 0;
    }

static int tipsyNextFile(FIO fio) {
    fioTipsy *tio = (fioTipsy *)fio;
    if (tio->iFile+1 >= fio->nFiles) return 1;
    tio->iFile++;
    return tipsySwitchFile(fio);
    }

/* DARK PARTICLES */
static int tipsyReadNativeDark(FIO fio,
    uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft,float *pfPot) {
    fioTipsy *tio = (fioTipsy *)fio;
    int rc;
    int d;
    float fTmp[3];
    assert(fio->eFormat == FIO_FORMAT_TIPSY && fio->eMode==FIO_MODE_READING);

    if (tio->iOrder>=tio->fio.iFirst[tio->iFile+1])
	if (tipsyNextFile(fio)) return 0;
    *piOrder = tio->iOrder++ + tio->iStart;

    rc = fread(pfMass,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    if (tio->bDoublePos) {
	rc = fread(pdPos,sizeof(double),3,tio->fp); if (rc!=3) return 0;
	}
    else {
	rc = fread(fTmp,sizeof(float),3,tio->fp); if (rc!=3) return 0;
	for(d=0;d<3;d++) pdPos[d] = fTmp[d];
	}
    if (tio->bDoubleVel) {
	rc = fread(pdVel,sizeof(double),3,tio->fp); if (rc!=3) return 0;
	}
    else {
	rc = fread(fTmp,sizeof(float),3,tio->fp); if (rc!=3) return 0;
	for(d=0;d<3;d++) pdVel[d] = fTmp[d];
	}
    rc = fread(pfSoft,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fread(pfPot,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    return 1;
    }

static int tipsyWriteNativeDark(FIO fio,
    uint64_t iOrder,const double *pdPos,const double *pdVel,
    float fMass,float fSoft,float fPot) {
    fioTipsy *tio = (fioTipsy *)fio;
    int rc;
    int d;
    float fTmp[3];

    assert(fio->eFormat == FIO_FORMAT_TIPSY && fio->eMode==FIO_MODE_WRITING);
    assert(iOrder == tio->iOrder++ + tio->iStart);
    rc = fwrite(&fMass,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    if (tio->bDoublePos) {
	rc = fwrite(pdPos,sizeof(double),3,tio->fp); if (rc!=3) return 0;
	}
    else {
	for(d=0;d<3;d++) fTmp[d] = pdPos[d];
	rc = fwrite(fTmp,sizeof(float),3,tio->fp); if (rc!=3) return 0;
	}
    if (tio->bDoubleVel) {
	rc = fwrite(pdVel,sizeof(double),3,tio->fp); if (rc!=3) return 0;
	}
    else {
	for(d=0;d<3;d++) fTmp[d] = pdVel[d];
	rc = fwrite(fTmp,sizeof(float),3,tio->fp); if (rc!=3) return 0;
	}
    rc = fwrite(&fSoft,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fwrite(&fPot,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    return 1;
    }

static int tipsyReadStandardDark(FIO fio,
    uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft,float *pfPot) {
    fioTipsy *tio = (fioTipsy *)fio;
    int d;
    float fTmp;

    assert(fio->eFormat == FIO_FORMAT_TIPSY && fio->eMode==FIO_MODE_READING);
    if (tio->iOrder>=tio->fio.iFirst[tio->iFile+1])
	if (tipsyNextFile(fio)) return 0;
    *piOrder = tio->iOrder++ + tio->iStart;
    if (!xdr_float(&tio->xdr,pfMass)) return 0;
    for(d=0;d<3;d++) {
	if (tio->bDoublePos) {
	    if (!xdr_double(&tio->xdr,pdPos+d)) return 0;
	    }
	else {
	    if (!xdr_float(&tio->xdr,&fTmp)) return 0;
	    pdPos[d] = fTmp;
	    }
	}
    for(d=0;d<3;d++) {
	if (tio->bDoubleVel) {
	    if (!xdr_double(&tio->xdr,pdVel+d)) return 0;
	    }
	else {
	    if (!xdr_float(&tio->xdr,&fTmp)) return 0;
	    pdVel[d] = fTmp;
	    }
	}
    if (!xdr_float(&tio->xdr,pfSoft)) return 0;
    if (!xdr_float(&tio->xdr,pfPot)) return 0;
    return 1;
    }

static int tipsyWriteStandardDark(FIO fio,
    uint64_t iOrder,const double *pdPos,const double *pdVel,
    float fMass,float fSoft,float fPot) {
    fioTipsy *tio = (fioTipsy *)fio;
    int d;
    float fTmp;
    double dTmp;

    assert(fio->eFormat == FIO_FORMAT_TIPSY && fio->eMode==FIO_MODE_WRITING);
//    assert(iOrder == tio->iOrder++ + tio->iStart);  // JW: needs to handle non-contiguous iOrders
    if (!xdr_float(&tio->xdr,&fMass)) return 0;
    for(d=0;d<3;d++) {
	if (tio->bDoublePos) {
	    dTmp = pdPos[d]; if (!xdr_double(&tio->xdr,&dTmp)) return 0;
	    }
	else {
	    fTmp = pdPos[d]; if (!xdr_float(&tio->xdr,&fTmp)) return 0;
	    }
	}
    for(d=0;d<3;d++) {
	if (tio->bDoubleVel) {
	    dTmp = pdVel[d]; if (!xdr_double(&tio->xdr,&dTmp)) return 0;
	    }
	else {
	    fTmp = pdVel[d]; if (!xdr_float(&tio->xdr,&fTmp)) return 0;
	    }
	}
    if (!xdr_float(&tio->xdr,&fSoft)) return 0;
    if (!xdr_float(&tio->xdr,&fPot)) return 0;
    return 1;
    }

/* SPH PARTICLES */
static int tipsyReadNativeSph(
    FIO fio,uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft, float *pfPot,
    float *pfRho,float *pfTemp, float *pfMetals) {
    fioTipsy *tio = (fioTipsy *)fio;
    int rc;
    int d;
    float fTmp[3];

    assert(fio->eFormat == FIO_FORMAT_TIPSY && fio->eMode==FIO_MODE_READING);
    if (tio->iOrder>=tio->fio.iFirst[tio->iFile+1])
	if (tipsyNextFile(fio)) return 0;
    *piOrder = tio->iOrder++ + tio->iStart;
    rc = fread(pfMass,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    if (tio->bDoublePos) {
	rc = fread(pdPos,sizeof(double),3,tio->fp); if (rc!=3) return 0;
	}
    else {
	rc = fread(fTmp,sizeof(float),3,tio->fp); if (rc!=3) return 0;
	for(d=0;d<3;d++) pdPos[d] = fTmp[d];
	}
    if (tio->bDoubleVel) {
	rc = fread(pdVel,sizeof(double),3,tio->fp); if (rc!=3) return 0;
	}
    else {
	rc = fread(fTmp,sizeof(float),3,tio->fp); if (rc!=3) return 0;
	for(d=0;d<3;d++) pdVel[d] = fTmp[d];
	}
    rc = fread(pfRho,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fread(pfTemp,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fread(pfSoft,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fread(pfMetals,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fread(pfPot,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    return 1;
    }


static int tipsyWriteNativeSph(
    struct fioInfo *fio,uint64_t iOrder,const double *pdPos,const double *pdVel,
    float fMass,float fSoft,float fPot,
    float fRho,float fTemp,float fMetals) {
    fioTipsy *tio = (fioTipsy *)fio;
    int rc;
    int d;
    float fTmp[3];

    assert(fio->eFormat == FIO_FORMAT_TIPSY && fio->eMode==FIO_MODE_WRITING);
    assert(iOrder == tio->iOrder++ + tio->iStart);

    rc = fwrite(&fMass,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    if (tio->bDoublePos) {
	rc = fwrite(pdPos,sizeof(double),3,tio->fp); if (rc!=3) return 0;
	}
    else {
	for(d=0;d<3;d++) fTmp[d] = pdPos[d];
	rc = fwrite(fTmp,sizeof(float),3,tio->fp); if (rc!=3) return 0;
	}
    if (tio->bDoubleVel) {
	rc = fwrite(pdVel,sizeof(double),3,tio->fp); if (rc!=3) return 0;
	}
    else {
	for(d=0;d<3;d++) fTmp[d] = pdVel[d];
	rc = fwrite(fTmp,sizeof(float),3,tio->fp); if (rc!=3) return 0;
	}
    rc = fwrite(&fRho,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fwrite(&fTemp,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fwrite(&fSoft,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fwrite(&fMetals,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fwrite(&fPot,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    return 1;
    }

static int tipsyReadStandardSph(
    FIO fio,uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft, float *pfPot,
    float *pfRho, float *pfTemp, float *pfMetals) {
    fioTipsy *tio = (fioTipsy *)fio;
    int d;
    float fTmp;

    assert(fio->eFormat == FIO_FORMAT_TIPSY && fio->eMode==FIO_MODE_READING);
    if (tio->iOrder>=tio->fio.iFirst[tio->iFile+1])
	if (tipsyNextFile(fio)) return 0;
    *piOrder = tio->iOrder++ + tio->iStart;
    if (!xdr_float(&tio->xdr,pfMass)) return 0;
    for(d=0;d<3;d++) {
	if (tio->bDoublePos) {
	    if (!xdr_double(&tio->xdr,pdPos+d)) return 0;
	    }
	else {
	    if (!xdr_float(&tio->xdr,&fTmp)) return 0;
	    pdPos[d] = fTmp;
	    }
	}
    for(d=0;d<3;d++) {
	if (tio->bDoubleVel) {
	    if (!xdr_double(&tio->xdr,pdVel+d)) return 0;
	    }
	else {
	    if (!xdr_float(&tio->xdr,&fTmp)) return 0;
	    pdVel[d] = fTmp;
	    }
	}
    if (!xdr_float(&tio->xdr,pfRho)) return 0;
    if (!xdr_float(&tio->xdr,pfTemp)) return 0;
    if (!xdr_float(&tio->xdr,pfSoft)) return 0;
    if (!xdr_float(&tio->xdr,pfMetals)) return 0;
    if (!xdr_float(&tio->xdr,pfPot)) return 0;
    return 1;
    }

static int tipsyWriteStandardSph(
    struct fioInfo *fio,uint64_t iOrder,const double *pdPos,const double *pdVel,
    float fMass,float fSoft,float fPot,
    float fRho,float fTemp,float fMetals) {
    fioTipsy *tio = (fioTipsy *)fio;
    int d;
    float fTmp;
    double dTmp;

    assert(fio->eFormat == FIO_FORMAT_TIPSY && fio->eMode==FIO_MODE_WRITING);
//    assert(iOrder == tio->iOrder++ + tio->iStart); //JW -- non-contiguous iOrder fix needed

    if (!xdr_float(&tio->xdr,&fMass)) return 0;
    for(d=0;d<3;d++) {
	if (tio->bDoublePos) {
	    dTmp = pdPos[d];
	    if (!xdr_double(&tio->xdr,&dTmp)) return 0;
	    }
	else {
	    fTmp = pdPos[d];
	    if (!xdr_float(&tio->xdr,&fTmp)) return 0;
	    }
	}
    for(d=0;d<3;d++) {
	if (tio->bDoubleVel) {
	    dTmp = pdVel[d];
	    if (!xdr_double(&tio->xdr,&dTmp)) return 0;
	    }
	else {
	    fTmp = pdVel[d];
	    if (!xdr_float(&tio->xdr,&fTmp)) return 0;
	    }
	}
    if (!xdr_float(&tio->xdr,&fRho)) return 0;
    if (!xdr_float(&tio->xdr,&fTemp)) return 0;
    if (!xdr_float(&tio->xdr,&fSoft)) return 0;
    if (!xdr_float(&tio->xdr,&fMetals)) return 0;
    if (!xdr_float(&tio->xdr,&fPot)) return 0;
    return 1;
    }

/* STAR PARTICLES */
static int tipsyReadNativeStar(
    FIO fio,uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft,float *pfPot,
    float *pfMetals, float *pfTform) {
    fioTipsy *tio = (fioTipsy *)fio;
    int rc;
    int d;
    float fTmp[3];

    assert(fio->eFormat == FIO_FORMAT_TIPSY && fio->eMode==FIO_MODE_READING);
    if (tio->iOrder>=tio->fio.iFirst[tio->iFile+1])
	if (!tipsyNextFile(fio)) return 0;
    *piOrder = tio->iOrder++ + tio->iStart;
    rc = fread(pfMass,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    if (tio->bDoublePos) {
	rc = fread(pdPos,sizeof(double),3,tio->fp); if (rc!=3) return 0;
	}
    else {
	rc = fread(fTmp,sizeof(float),3,tio->fp); if (rc!=3) return 0;
	for(d=0;d<3;d++) pdPos[d] = fTmp[d];
	}
    if (tio->bDoubleVel) {
	rc = fread(pdVel,sizeof(double),3,tio->fp); if (rc!=3) return 0;
	}
    else {
	rc = fread(fTmp,sizeof(float),3,tio->fp); if (rc!=3) return 0;
	for(d=0;d<3;d++) pdVel[d] = fTmp[d];
	}
    rc = fread(pfMetals,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fread(pfTform,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fread(pfSoft,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fread(pfPot,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    return 1;
    }

static int tipsyWriteNativeStar(
    struct fioInfo *fio,uint64_t iOrder,const double *pdPos,const double *pdVel,
    float fMass,float fSoft,float fPot,float fMetals,float fTform) {
    fioTipsy *tio = (fioTipsy *)fio;
    int rc;
    int d;
    float fTmp[3];

    assert(fio->eFormat == FIO_FORMAT_TIPSY && fio->eMode==FIO_MODE_WRITING);
    assert(iOrder == tio->iOrder++ + tio->iStart);
    rc = fwrite(&fMass,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    if (tio->bDoublePos) {
	rc = fwrite(pdPos,sizeof(double),3,tio->fp); if (rc!=3) return 0;
	}
    else {
	for(d=0;d<3;d++) fTmp[d] = pdPos[d];
	rc = fwrite(fTmp,sizeof(float),3,tio->fp); if (rc!=3) return 0;
	}
    if (tio->bDoubleVel) {
	rc = fwrite(pdVel,sizeof(double),3,tio->fp); if (rc!=3) return 0;
	}
    else {
	for(d=0;d<3;d++) fTmp[d] = pdVel[d];
	rc = fwrite(fTmp,sizeof(float),3,tio->fp); if (rc!=3) return 0;
	}
    rc = fwrite(&fMetals,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fwrite(&fTform,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fwrite(&fSoft,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fwrite(&fPot,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    return 1;
    }

static int tipsyReadStandardStar(
    FIO fio,uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft,float *pfPot,
    float *pfMetals, float *pfTform) {
    fioTipsy *tio = (fioTipsy *)fio;
    int d;
    float fTmp;

    assert(fio->eFormat == FIO_FORMAT_TIPSY && fio->eMode==FIO_MODE_READING);
    if (tio->iOrder>=tio->fio.iFirst[tio->iFile+1])
	if (!tipsyNextFile(fio)) return 0;
    *piOrder = tio->iOrder++ + tio->iStart;
    if (!xdr_float(&tio->xdr,pfMass)) return 0;
    for(d=0;d<3;d++) {
	if (tio->bDoublePos) {
	    if (!xdr_double(&tio->xdr,pdPos+d)) return 0;
	    }
	else {
	    if (!xdr_float(&tio->xdr,&fTmp)) return 0;
	    pdPos[d] = fTmp;
	    }
	}
    for(d=0;d<3;d++) {
	if (tio->bDoubleVel) {
	    if (!xdr_double(&tio->xdr,pdVel+d)) return 0;
	    }
	else {
	    if (!xdr_float(&tio->xdr,&fTmp)) return 0;
	    pdVel[d] = fTmp;
	    }
	}
    if (!xdr_float(&tio->xdr,pfMetals)) return 0;
    if (!xdr_float(&tio->xdr,pfTform)) return 0;
    if (!xdr_float(&tio->xdr,pfSoft)) return 0;
    if (!xdr_float(&tio->xdr,pfPot)) return 0;
    return 1;
    }

static int tipsyWriteStandardStar(
    struct fioInfo *fio,uint64_t iOrder,const double *pdPos,const double *pdVel,
    float fMass,float fSoft,float fPot,float fMetals,float fTform) {
    fioTipsy *tio = (fioTipsy *)fio;
    int d;
    float fTmp;
    double dTmp;

    assert(fio->eFormat == FIO_FORMAT_TIPSY && fio->eMode==FIO_MODE_WRITING);
//    assert(iOrder == tio->iOrder++ + tio->iStart); // JW non-contiguous iOrder fix needed
    if (!xdr_float(&tio->xdr,&fMass)) return 0;
    for(d=0;d<3;d++) {
	if (tio->bDoublePos) {
	    dTmp = pdPos[d];
	    if (!xdr_double(&tio->xdr,&dTmp)) return 0;
	    }
	else {
	    fTmp = pdPos[d];
	    if (!xdr_float(&tio->xdr,&fTmp)) return 0;
	    }
	}
    for(d=0;d<3;d++) {
	if (tio->bDoubleVel) {
	    dTmp = pdVel[d];
	    if (!xdr_double(&tio->xdr,&dTmp)) return 0;
	    }
	else {
	    fTmp = pdVel[d];
	    if (!xdr_float(&tio->xdr,&fTmp)) return 0;
	    }
	}
    if (!xdr_float(&tio->xdr,&fMetals)) return 0;
    if (!xdr_float(&tio->xdr,&fTform)) return 0;
    if (!xdr_float(&tio->xdr,&fSoft)) return 0;
    if (!xdr_float(&tio->xdr,&fPot)) return 0;
    return 1;
    }

static void tipsyCloseNative(FIO fio) {
    fioTipsy *tio = (fioTipsy *)fio;
    assert(fio->eFormat == FIO_FORMAT_TIPSY);
    fclose(tio->fp);
    if (tio->fpBuffer) free(tio->fpBuffer);
    free(tio);
    }

static void tipsyCloseStandard(FIO fio) {
    fioTipsy *tio = (fioTipsy *)fio;
    assert(fio->eFormat == FIO_FORMAT_TIPSY);
    xdr_destroy(&tio->xdr);
    tipsyCloseNative(fio);
    }

int fioTipsyIsDouble(FIO fio) {
    fioTipsy *tio = (fioTipsy *)fio;
    assert(fio->eFormat == FIO_FORMAT_TIPSY);
    return tio->bDoublePos;
    }

int fioTipsyIsStandard(FIO fio) {
    fioTipsy *tio = (fioTipsy *)fio;
    assert(fio->eFormat == FIO_FORMAT_TIPSY);
    return tio->fio.fcnReadDark == tipsyReadStandardDark;
    }

/*
** Given an offset, calculate the particle number
*/
static uint64_t tipsyParticle(uint64_t iByte,uint64_t nHdrSize,
    uint64_t nSph, uint64_t nDark, uint64_t nStar,
    int nSphSize, int nDarkSize, int nStarSize ) {
    uint64_t iPart,n;

    iPart = 0;
    iByte -= nHdrSize;

    n = iByte / nSphSize;
    if ( n > nSph ) n = nSph;
    iPart += n;
    iByte -= nSphSize * n;

    n = iByte / nDarkSize;
    if ( n > nDark ) n = nDark;
    iPart += n;
    iByte -= nDarkSize * n;

    n = iByte / nStarSize;
    assert( n<=nStar );
    iPart += n;
    iByte -= nStarSize * n;

    assert(iByte==0);

    return iPart;
    }


/*
** Calculate the absolute offset of a given particle in a Tipsy file.
*/
static uint64_t tipsyOffset(
    uint64_t iPart,uint64_t nHdrSize,
    uint64_t nSph, uint64_t nDark, uint64_t nStar,
    int nSphSize, int nDarkSize, int nStarSize ) {
    uint64_t iByte;

    iByte = nHdrSize;
    if (iPart<nSph) {
	iByte += iPart * nSphSize;
	}
    else {
	iByte += nSph * nSphSize;
	iPart -= nSph;

	if (iPart<nDark) {
	    iByte += iPart * nDarkSize;
	    }
	else {
	    iByte += nDark * nDarkSize;
	    iPart -= nDark;

	    if (iPart<=nStar) {
		iByte += iPart * nStarSize;
		}
	    else {
		errno = ESPIPE;
		return 0;
		}
	    }
	}
    return iByte;
    }

static int tipsySeek(FIO fio,uint64_t iPart,FIO_SPECIES eSpecies) {
    fioTipsy *tio = (fioTipsy *)fio;
    uint64_t iByte;
    off_t iOffset;
    int rc;

    switch(eSpecies) {
    case FIO_SPECIES_STAR: iPart += tio->fio.nSpecies[FIO_SPECIES_DARK];
    case FIO_SPECIES_DARK: iPart += tio->fio.nSpecies[FIO_SPECIES_SPH];
    case FIO_SPECIES_SPH:  break;
    case FIO_SPECIES_ALL:  break;
    default:
	fprintf(stderr,"Invalid particle species\n");
	exit(1);
	}
    tio->iOrder = iPart;

    /* Seek to a different file if necessary */
    if (iPart < fio->iFirst[tio->iFile]) {
	while(tio->iFile!=0)
	    if (fio->iFirst[--tio->iFile]<iPart) break;
	rc = tipsySwitchFile(fio); assert(rc==0);
	if (rc) return rc;
	}
    else if (iPart >= fio->iFirst[tio->iFile+1]) {
	while(++tio->iFile<fio->nFiles)
	    if (fio->iFirst[tio->iFile]<=iPart) break;
	rc = tipsySwitchFile(fio); assert(rc==0);
	if (rc) return rc;
	}

    /* For file fragments, the first particle is not 0.  Adjust accordingly. */
    iByte = tipsyOffset(iPart,0,
			tio->fio.nSpecies[FIO_SPECIES_SPH],
			tio->fio.nSpecies[FIO_SPECIES_DARK],
			tio->fio.nSpecies[FIO_SPECIES_STAR],
			TIPSY_SPH_SIZE(1,3*tio->bDoublePos+3*tio->bDoubleVel),
			TIPSY_DARK_SIZE(1,3*tio->bDoublePos+3*tio->bDoubleVel),
			TIPSY_STAR_SIZE(1,3*tio->bDoublePos+3*tio->bDoubleVel))
	- tipsyOffset(fio->iFirst[tio->iFile],0,
		      tio->fio.nSpecies[FIO_SPECIES_SPH],
		      tio->fio.nSpecies[FIO_SPECIES_DARK],
		      tio->fio.nSpecies[FIO_SPECIES_STAR],
		      TIPSY_SPH_SIZE(1,3*tio->bDoublePos+3*tio->bDoubleVel),
		      TIPSY_DARK_SIZE(1,3*tio->bDoublePos+3*tio->bDoubleVel),
		      TIPSY_STAR_SIZE(1,3*tio->bDoublePos+3*tio->bDoubleVel))
	+ (tio->iFile?0:tio->nHdrSize);

    iOffset = iByte;
    assert(iOffset==iByte);
    rc = safe_fseek(tio->fp,iByte); assert(rc==0);
    return rc;
    }

static int tipsySeekStandard(FIO fio,uint64_t iPart,FIO_SPECIES eSpecies) {
    /*fioTipsy *tio = (fioTipsy *)fio;*/
    int rc;
    /*assert(fio->eFormat == FIO_FORMAT_TIPSY);*/
    /*xdr_destroy(&tio->xdr);*/
    rc = tipsySeek(fio,iPart,eSpecies);
    /*xdrstdio_create(&tio->xdr,tio->fp,fio->eMode==FIO_MODE_READING?XDR_DECODE:XDR_ENCODE);*/
    return rc;
    }

static uint64_t tipsyCountSpecies(uint64_t N, uint64_t *iStart, uint64_t *iLast) {
    uint64_t n;

    if ( *iStart >= N ) {
	n = 0;
	*iStart -= N; *iLast -= N;
	}
    else if ( *iLast > N) {
	n = N - *iStart;
	*iStart = 0; *iLast -= N;
	}
    else {
	n = *iLast - *iStart;
	*iStart = *iLast = 0;
	}

    return n;
    }

static void tipsySetFunctions(fioTipsy *tio, int bDouble, int bStandard) {
    tio->fio.fcnGetAttr = tipsyGetAttr;
    tio->fio.fcnSpecies = tipsySpecies;

    tio->bDoubleVel = 0;
    tio->bDoublePos = bDouble || tio->bDoubleVel;

    if ( bStandard ) {
	tio->fio.fcnClose    = tipsyCloseStandard;
	tio->fio.fcnSeek     = tipsySeekStandard;
	tio->fio.fcnReadDark = tipsyReadStandardDark;
	tio->fio.fcnReadSph  = tipsyReadStandardSph;
	tio->fio.fcnReadStar = tipsyReadStandardStar;
	tio->fio.fcnWriteDark= tipsyWriteStandardDark;
	tio->fio.fcnWriteSph = tipsyWriteStandardSph;
	tio->fio.fcnWriteStar= tipsyWriteStandardStar;
	}
    else {
	tio->fio.fcnClose    = tipsyCloseNative;
	tio->fio.fcnSeek     = tipsySeek;
	tio->fio.fcnReadDark = tipsyReadNativeDark;
	tio->fio.fcnReadSph  = tipsyReadNativeSph;
	tio->fio.fcnReadStar = tipsyReadNativeStar;
	tio->fio.fcnWriteDark= tipsyWriteNativeDark;
	tio->fio.fcnWriteSph = tipsyWriteNativeSph;
	tio->fio.fcnWriteStar= tipsyWriteNativeStar;
	}
    }

/*
** The pad field can extend the 32-bit integers to 40 bits, increasing
** the maximum number of particles from 4e9 to 1e12 (1 trillion).
** This checks if that is what has happened here.
*/
static void tipsySussHeader(tipsyHdr *h,uint64_t *pN, uint64_t *pDark,uint64_t *pSph, uint64_t *pStar) {
    *pN = h->nPad & 0x000000ff;
    *pN <<= 32;
    *pN += h->nBodies;

    *pSph = h->nPad & 0x0000ff00;
    *pSph <<= 32;
    *pSph += h->nSph;

    *pDark = h->nPad & 0x00ff0000;
    *pDark <<= 32;
    *pDark += h->nDark;

    *pStar = h->nPad & 0xff000000;
    *pStar <<= 32;
    *pStar += h->nStar;

    /* There may be junk in the pad field, in which case we ignore it */
    if ( *pN != *pDark + *pSph + *pStar ) {
	*pN = h->nBodies;
	*pDark = h->nDark;
	*pSph  = h->nSph;
	*pStar = h->nStar;
	}
    }


/*
** Read the header and detect native or standard
*/
static int tipsyDetectHeader(fioTipsy *tio, int bDouble) {
    tipsyHdr h;
    uint64_t N,nDark,nSph,nStar;
    int rc;

    assert(tio->fio.eFormat == FIO_FORMAT_TIPSY);
    rc = fread(&h,sizeof(h),1,tio->fp);
    if (rc!=1) {
	fprintf(stderr,"Error reading Tipsy header\n");
	return 0;
	}

    /* Check for native binary format */
    if (h.nDim>=1 && h.nDim<=3) {
	tipsySussHeader(&h,&N,&nDark,&nSph,&nStar);
	if ( nDark + nSph + nStar != N ) {
	    fprintf(stderr,"Tipsy Header mismatch:"
		    " nDim=%u nDark=%lu nSph=%lu nStar=%lu nBodies=%lu\n",
		    h.nDim, nDark, nSph, nStar, N);
	    return 0;
	    }
	tipsySetFunctions(tio,bDouble,0);
	}

    /* Okay, this must be "standard" format */
    else {
	rewind(tio->fp);
	xdrstdio_create(&tio->xdr,tio->fp,XDR_DECODE);
	rc = xdrHeader(&tio->xdr,&h);
	if (rc!=1) {
	    perror("Error reading tipsy header");
	    return 0;
	    }
	tipsySussHeader(&h,&N,&nDark,&nSph,&nStar);

	if (h.nDim<1 || h.nDim>3 || nDark + nSph + nStar != N ) {
	    fprintf(stderr,"Tipsy Header mismatch:"
		    " nDim=%u nDark=%lu nSph=%lu nStar=%lu nBodies=%lu\n",
		    h.nDim, nDark, nSph, nStar, N);
	    xdr_destroy(&tio->xdr);
	    return 0;
	    }
	tipsySetFunctions(tio,bDouble,1);
	}
    tio->nHdrSize = sizeof(h);
    tio->dTime = h.dTime;
    tio->fio.nSpecies[FIO_SPECIES_DARK] = nDark;
    tio->fio.nSpecies[FIO_SPECIES_SPH]  = nSph;
    tio->fio.nSpecies[FIO_SPECIES_STAR] = nStar;
    return 1;
    }

static int tipsyWriteHeader(fioTipsy *tio, int bDouble, int bStandard) {
    FIO fio = &tio->fio;
    tipsyHdr h;
    uint64_t N, nSph, nDark, nStar;
    int rc;

    nSph  = fioGetN(fio,FIO_SPECIES_SPH);
    nDark = fioGetN(fio,FIO_SPECIES_DARK);
    nStar = fioGetN(fio,FIO_SPECIES_STAR);
    N = nSph + nDark + nStar;

    h.dTime = tio->dTime;
    h.nBodies = N;
    h.nDim    = 3;
    h.nSph    = nSph;
    h.nDark   = nDark;
    h.nStar   = nStar;
    h.nPad    = ((N>>32)&0xff) + ((nSph>>24)&0xff00)
	+ ((nDark>>16)&0xff0000) + ((nStar>>8)&0xff000000);

    tio->nHdrSize = sizeof(h);
    if (bStandard) {
	rc = xdrHeader(&tio->xdr,&h);
	if (rc!=1) {
	    perror("Error writing tipsy header");
	    return 0;
	    }
	}
    else {
	rc = fwrite(&h,sizeof(h),1,tio->fp);
	if (rc!=1) {
	    fprintf(stderr,"Error writing Tipsy header\n");
	    return 0;
	    }
	}

    return 1;
    }

FIO fioTipsyOpenMany(int nFiles, const char * const *fileNames) {
    fioTipsy *tio;
    int i;
    off_t *nSizes;
    uint64_t nSize, nDark, nSph, nStar;

    tio = malloc(sizeof(fioTipsy));
    assert(tio!=NULL);
    tio->fio.eFormat = FIO_FORMAT_TIPSY;
    tio->fio.eMode   = FIO_MODE_READING;
    tio->fio.pszFiles= NULL;
    tio->fio.iFirst  = NULL;
    tio->fio.nFiles  = 0;
    tio->iOrder = tio->iStart = 0;

    fileScan(&tio->fio,nFiles,fileNames);

    tio->fp = fopen(fileNames[0],"r");
    tio->iFile = 0;
    if (tio->fp==NULL) {
	free(tio);
	return NULL;
	}
    tio->fpBuffer = malloc(TIO_BUFFER_SIZE);
    if (tio->fpBuffer != NULL) setvbuf(tio->fp,tio->fpBuffer,_IOFBF,TIO_BUFFER_SIZE);

    if ( !tipsyDetectHeader(tio,/*bDouble*/0) ) {
	fclose(tio->fp);
	if (tio->fpBuffer) free(tio->fpBuffer);
	free(tio);
	return NULL;
	}

    tio->fio.nSpecies[FIO_SPECIES_ALL] = 0;
    for( i=1; i<FIO_SPECIES_LAST; i++)
	tio->fio.nSpecies[FIO_SPECIES_ALL] += tio->fio.nSpecies[i];

    /*
    ** We know how many particles we are supposed to have at this point, and if
    ** this is a native or standard tipsy binary.  What we don't know is how
    ** many particles are in each file.  We have to stat the files to find this.
    */
    nSizes = malloc(sizeof(off_t)*tio->fio.nFiles);
    nSize = 0;
    for( i=0; i<tio->fio.nFiles; i++) {
	struct stat s;
	/* The file/directory needs to exist */
	if ( stat(tio->fio.pszFiles[i],&s) != 0 ) return NULL;
	if ( !S_ISREG(s.st_mode) ) return NULL;
	nSize += (nSizes[i] = s.st_size);
	}

    nDark = tio->fio.nSpecies[FIO_SPECIES_DARK];
    nSph  = tio->fio.nSpecies[FIO_SPECIES_SPH];
    nStar = tio->fio.nSpecies[FIO_SPECIES_STAR];

    if (TIPSY_TOTAL_SIZE(nDark,nSph,nStar,0) == nSize ) {
	tio->bDoublePos = 0;
	tio->bDoubleVel = 0;
	}
    else if (TIPSY_TOTAL_SIZE(nDark,nSph,nStar,3) == nSize ) {
	tio->bDoublePos = 1;
	tio->bDoubleVel = 0;
	}
    else if (TIPSY_TOTAL_SIZE(nDark,nSph,nStar,6) == nSize ) {
	tio->bDoublePos = 1;
	tio->bDoubleVel = 1;
	}
    else {
	fprintf(stderr,"Tipsy file size mismatch!\n");
	return NULL;
	}

    /*
    ** Okay, now we need to calculate the starting particle number for each
    ** file fragment.
    */
    tio->fio.iFirst[0] = 0;
    nSize = 0;
    for(i=1; i<=tio->fio.nFiles; i++ ) {
	nSize += nSizes[i-1];
	tio->fio.iFirst[i] = tipsyParticle(
	    nSize,tio->nHdrSize,nSph,nDark,nStar,
	    TIPSY_SPH_SIZE(1,3*tio->bDoublePos+3*tio->bDoubleVel),
	    TIPSY_DARK_SIZE(1,3*tio->bDoublePos+3*tio->bDoubleVel),
	    TIPSY_STAR_SIZE(1,3*tio->bDoublePos+3*tio->bDoubleVel));
	}
    assert(tio->fio.iFirst[tio->fio.nFiles]==tio->fio.nSpecies[FIO_SPECIES_ALL]);

    free(nSizes);

    return &tio->fio;
    }

FIO fioTipsyOpen(const char *fileName,int bDouble) {
    return fioTipsyOpenMany(1,&fileName);
    }

FIO fioTipsyCreate(const char *fileName,int bDouble,int bStandard,
		   double dTime,uint64_t nSph, uint64_t nDark, uint64_t nStar) {
    fioTipsy *tio;
    int i;

    tio = malloc(sizeof(fioTipsy));
    assert(tio!=NULL);
    tio->fio.eFormat = FIO_FORMAT_TIPSY;
    tio->fio.eMode   = FIO_MODE_WRITING;
    tio->fio.pszFiles= NULL;
    tio->fio.iFirst = malloc(sizeof(uint64_t)*2);
    assert(tio->fio.iFirst);
    tio->fio.nFiles  = 1;

    for( i=0; i<FIO_SPECIES_LAST; i++)
	tio->fio.nSpecies[i] = 0;
    tio->fio.nSpecies[FIO_SPECIES_SPH]  = nSph;
    tio->fio.nSpecies[FIO_SPECIES_DARK] = nDark;
    tio->fio.nSpecies[FIO_SPECIES_STAR] = nStar;
    for( i=1; i<FIO_SPECIES_LAST; i++)
	tio->fio.nSpecies[FIO_SPECIES_ALL] += tio->fio.nSpecies[i];

    tio->fio.iFirst[0] = 0;
    tio->fio.iFirst[1] = tio->fio.nSpecies[FIO_SPECIES_ALL];

    tio->iOrder = tio->iStart = 0;
    tio->dTime = dTime;

    tio->fp = fopen(fileName,"w");
    tio->iFile = 0;
    if (tio->fp==NULL) {
	free(tio);
	return NULL;
	}
    tio->fpBuffer = malloc(TIO_BUFFER_SIZE);
    if (tio->fpBuffer != NULL) setvbuf(tio->fp,tio->fpBuffer,_IOFBF,TIO_BUFFER_SIZE);
    if (bStandard) xdrstdio_create(&tio->xdr,tio->fp,XDR_ENCODE);

    tipsySetFunctions(tio,bDouble,bStandard);
    tipsyWriteHeader(tio,bDouble,bStandard);

    return &tio->fio;
    }

FIO fioTipsyAppend(const char *fileName,int bDouble,int bStandard) {
    fioTipsy *tio;
    int i;

    tio = malloc(sizeof(fioTipsy));
    assert(tio!=NULL);
    tio->fio.eFormat = FIO_FORMAT_TIPSY;
    tio->fio.eMode   = FIO_MODE_WRITING;
    tio->iOrder = tio->iStart = 0;
    tio->fio.pszFiles= NULL;
    tio->fio.iFirst = malloc(sizeof(uint64_t)*2);
    assert(tio->fio.iFirst);
    tio->fio.nFiles  = 1;

    tio->fp = fopen(fileName,"r+");
    tio->iFile = 0;
    if (tio->fp==NULL) {
	free(tio);
	return NULL;
	}
    tio->fpBuffer = malloc(TIO_BUFFER_SIZE);
    if (tio->fpBuffer != NULL) setvbuf(tio->fp,tio->fpBuffer,_IOFBF,TIO_BUFFER_SIZE);

    if ( !tipsyDetectHeader(tio,bDouble) ) {
	fclose(tio->fp);
	if (tio->fpBuffer) free(tio->fpBuffer);
	free(tio);
	return NULL;
	}
    tio->fio.nSpecies[FIO_SPECIES_ALL] = 0;
    for( i=1; i<FIO_SPECIES_LAST; i++)
	tio->fio.nSpecies[FIO_SPECIES_ALL] += tio->fio.nSpecies[i];

    tio->fio.iFirst[0] = 0;
    tio->fio.iFirst[1] = tio->fio.nSpecies[FIO_SPECIES_ALL];

    if (bStandard) xdrstdio_create(&tio->xdr,tio->fp,XDR_ENCODE);
    tipsySetFunctions(tio,bDouble,bStandard);

    return &tio->fio;
    }

/*
** Create a partial Tipsy file.  A partial Tipsy file is one that has been split
** into multiple files.  Only the first part has the header.  To rebuild the
** original file, the parts need only be concatenated together.
*/
FIO fioTipsyCreatePart(const char *fileName,int bAppend,int bDouble,int bStandard,
		       double dTime, uint64_t nSph, uint64_t nDark, uint64_t nStar,
		       uint64_t iStart) {
    fioTipsy *tio;
    uint64_t iLast;
    int i;

    tio = malloc(sizeof(fioTipsy));
    assert(tio!=NULL);
    tio->fio.eFormat = FIO_FORMAT_TIPSY;
    tio->fio.eMode   = FIO_MODE_WRITING;
    tio->fio.pszFiles= NULL;
    tio->fio.iFirst = malloc(sizeof(uint64_t)*2);
    assert(tio->fio.iFirst);
    tio->fio.nFiles  = 1;

    tio->iOrder = 0;
    tio->iStart = iStart;
    tio->dTime = dTime;

    for( i=0; i<FIO_SPECIES_LAST; i++)
	tio->fio.nSpecies[i] = 0;

    /* Adjust the local number of particles in this fragment */
    iLast = nSph + nDark + nStar;
    tio->fio.nSpecies[FIO_SPECIES_SPH]  = nSph;
    tio->fio.nSpecies[FIO_SPECIES_DARK] = nDark;
    tio->fio.nSpecies[FIO_SPECIES_STAR] = nStar;
    for( i=1; i<FIO_SPECIES_LAST; i++)
	tio->fio.nSpecies[FIO_SPECIES_ALL] += tio->fio.nSpecies[i];

    tio->fio.iFirst[0] = 0;
    tio->fio.iFirst[1] = tio->fio.nSpecies[FIO_SPECIES_ALL];

    tio->fp = fopen(fileName,bAppend?"r+":"w");
    tio->iFile = 0;
    if (tio->fp==NULL) {
	free(tio);
	return NULL;
	}
    tio->fpBuffer = malloc(TIO_BUFFER_SIZE);
    if (tio->fpBuffer != NULL) setvbuf(tio->fp,tio->fpBuffer,_IOFBF,TIO_BUFFER_SIZE);
    if (bStandard) xdrstdio_create(&tio->xdr,tio->fp,XDR_ENCODE);

    tipsySetFunctions(tio,bDouble,bStandard);
    if (tio->iStart) tio->nHdrSize = 0;
    else tipsyWriteHeader(tio,bDouble,bStandard);
    return &tio->fio;
    }

/******************************************************************************\
** General HDF5 support routines.  Not "fio" specific.
\******************************************************************************/
#ifdef USE_HDF5
#include <hdf5.h>

/* Returns the number of records in a dataset */
static hsize_t getSetSize(hid_t setID) {
    hid_t spaceID;
    hsize_t dims[2], maxs[2];
    hsize_t rc;

    spaceID = H5Dget_space(setID); assert(spaceID>=0);
    assert( H5Sis_simple(spaceID) > 0 );
    assert( H5Sget_simple_extent_ndims(spaceID) <= 3 );
    H5Sget_simple_extent_dims(spaceID,dims,maxs);
    H5Sclose(spaceID);
    return dims[0];
    }

/* Create a dataset inside a group (positions, velocities, etc.) */
static hid_t newSet(hid_t locID, const char *name, uint64_t chunk,
		    uint64_t count, int nDims, hid_t dataType ) {
    hid_t dataProperties, dataSpace, dataSet;
    hsize_t /*@alt uint64_t[]@*/ iDims[2];
    hsize_t /*@alt uint64_t[]@*/ iMax[2];
    hsize_t rc;

    /* Create a dataset property so we can set the chunk size */
    dataProperties = H5Pcreate( H5P_DATASET_CREATE );
    assert(dataProperties >= 0);
    iDims[0] = chunk;
    iDims[1] = 1;
    rc = H5Pset_chunk( dataProperties, nDims>1?2:1, iDims ); assert(rc>=0);

    /* Also request the FLETCHER checksum */
    rc = H5Pset_filter(dataProperties,H5Z_FILTER_FLETCHER32,0,0,NULL); assert(rc>=0);

    /* And the dataspace */
    iDims[0] = count;
    iDims[1] = nDims;
    iMax[0] = H5S_UNLIMITED;
    iMax[1] = nDims;

    dataSpace = H5Screate_simple( nDims>1?2:1, iDims, iMax );
    assert( dataSpace != H5I_INVALID_HID );

    /* Create the data set */
    dataSet = H5Dcreate( locID, name,
			 dataType, dataSpace, dataProperties );
    assert(dataSet != H5I_INVALID_HID);
    H5Pclose( dataProperties );
    H5Sclose( dataSpace );

    return dataSet;
    }

/* Read part of a set from disk into memory */
static void readSet(
    hid_t set_id,           /* set from which to read the data */
    void *pBuffer,          /* Buffer for the data */
    hid_t memType,          /* Data type in memory */
    hsize_t iOffset,        /* Offset into the set */
    hsize_t nRows,          /* Number of rows to read */
    hsize_t nCols ) {       /* Number of columns (normally 1 or 3) */
    hid_t memSpace, diskSpace;
    hsize_t dims[2], start[2];
    herr_t rc;

    dims[0] = nRows;
    dims[1] = nCols;
    memSpace = H5Screate_simple(nCols>1?2:1,dims,0);
    diskSpace = H5Dget_space(set_id);
    start[0] = iOffset;
    start[1] = 0;

    rc = H5Sselect_hyperslab(diskSpace,H5S_SELECT_SET,start,0,dims,0); assert(rc>=0);
    rc = H5Dread(set_id,memType,memSpace,diskSpace,H5P_DEFAULT,pBuffer); assert(rc>=0);
    rc = H5Sclose(memSpace); assert(rc>=0);
    rc = H5Sclose(diskSpace); assert(rc>=0);
    }


/* Add an attribute to a group */
static void writeAttribute( hid_t groupID, const char *name,
			    hid_t dataType, void *data ) {
    hid_t dataSpace, attrID;
    hsize_t dims = 1;
    herr_t rc;

    dataSpace = H5Screate_simple( 1, &dims, NULL ); assert(dataSpace!=H5I_INVALID_HID);
    attrID = H5Acreate( groupID,name,dataType,dataSpace,H5P_DEFAULT );
    assert(attrID!=H5I_INVALID_HID);
    rc = H5Awrite( attrID, dataType, data ); assert(rc>=0);
    rc = H5Aclose( attrID ); assert(rc>=0);
    rc = H5Sclose(dataSpace); assert(rc>=0);
    }

/* Read an attribute from a group */
static int readAttribute( hid_t groupID, const char *name,
			  hid_t dataType, void *data ) {
    hid_t attrID;
    herr_t rc;

    attrID = H5Aopen_name( groupID,name );
    if ( attrID == H5I_INVALID_HID ) return 0;
    rc = H5Aread( attrID, dataType, data ); assert(rc>=0);
    rc = H5Aclose( attrID ); assert(rc>=0);
    return 1;
    }

/* Write part of a set from memory */
static void writeSet(
    hid_t set_id,           /* set into which to write the data */
    const void *pBuffer,    /* Buffer containing the data */
    hid_t memType,          /* Data type in memory */
    hsize_t iOffset,        /* Offset into the set */
    hsize_t nRows,          /* Number of rows to write */
    hsize_t nCols ) {       /* Number of columns (normally 1 or 3) */
    hid_t memSpace, diskSpace;
    hsize_t dims[2], start[2], size[2];

    dims[0] = nRows;
    dims[1] = nCols;
    memSpace = H5Screate_simple(nCols>1?2:1,dims,0);
    size[0] = iOffset + nRows;
    size[1] = nCols;
    H5Dextend(set_id,size);
    diskSpace = H5Dget_space(set_id);
    start[0] = iOffset;
    start[1] = 0;
    H5Sselect_hyperslab(diskSpace,H5S_SELECT_SET,start,0,dims,0);
    H5Dwrite(set_id,memType,memSpace,diskSpace,H5P_DEFAULT,pBuffer);
    H5Sclose(memSpace);
    H5Sclose(diskSpace);
    }

/******************************************************************************\
** HDF5 FORMAT
\******************************************************************************/

#define CHUNK_SIZE 32768

#define ATTR_IORDER      "iOrder"

#define FIELD_POSITION   "position"
#define FIELD_VELOCITY   "velocity"
#define FIELD_POTENTIAL  "potential"
#define FIELD_ORDER      "order"
#define FIELD_CLASS      "class"
#define FIELD_CLASSES    "classes"

typedef struct {
    double v[3];
    } ioHDF5V3;

typedef uint64_t PINDEX;

typedef struct {
    uint8_t iClass;
    PINDEX  iOrderStart;       /* Start index of this class */
    double  fMass;             /* Particle mass */
    double  fSoft;             /* Softening */
    PINDEX nCount;            /* Number of particles with this class */
    } classEntry;

typedef struct {
    uint_fast32_t nClasses;   /* Number of different classes */
    classEntry Class[256];
    hid_t    setClass_id;
    uint8_t  *piClass;         /* Class index (or NULL) */
    double   *fMass;
    double   *fSoft;
    } IOCLASS;

typedef struct {
    PINDEX iStart;            /* Start iOrder number */
    PINDEX iNext;             /* Next iOrder number (starts at iStart) */
    PINDEX *iOrder;           /* Buffered iOrder numbers (or NULL) */
    hid_t  setOrder_id;
    } IOORDER;

typedef struct {
    PINDEX iOffset;           /* Particle offset into the file */
    PINDEX  nTotal;           /* Total number of particles in the file */
    uint_fast32_t iIndex;     /* Index of next particle in memory */
    uint_fast32_t nBuffered;  /* Number of buffered particles */
    hid_t group_id;
    hid_t setR_id;            /* set of Positions */
    ioHDF5V3 *R;              /* Positions (always present) */
    hid_t setV_id;            /* set of Velocities */
    ioHDF5V3 *V;              /* Velocities (always present) */
    hid_t setPot_id;
    float *pot;
    IOORDER order;
    IOCLASS class;
    } IOBASE;

typedef struct {
    struct fioInfo fio;
    hid_t fileID;
    FIO_SPECIES eCurrent;
    IOBASE base[FIO_SPECIES_LAST];
    } fioHDF5;

/* Create the type for the class table */
static hid_t makeClassType(hid_t floatType, int bStart) {
    hid_t tid;
    hid_t dataType = sizeof(PINDEX)==4 ? H5T_NATIVE_UINT32 : H5T_NATIVE_UINT64;

    tid = H5Tcreate (H5T_COMPOUND, sizeof(classEntry));
    assert(tid!=H5I_INVALID_HID);
    H5Tinsert(tid,"class",HOFFSET(classEntry,iClass), H5T_NATIVE_UINT8);
    H5Tinsert(tid,"mass", HOFFSET(classEntry,fMass), floatType);
    H5Tinsert(tid,"soft", HOFFSET(classEntry,fSoft), floatType);
    if ( bStart )
	H5Tinsert(tid,"start",HOFFSET(classEntry,iOrderStart), dataType);
    return tid;
    }

static void writeClassTable(IOBASE *Base ) {
    hid_t tid, set;

    if ( Base->nTotal > 0 ) {
	tid = makeClassType( H5T_NATIVE_DOUBLE, Base->class.piClass==NULL );

	set = newSet(Base->group_id, FIELD_CLASSES,
		     CHUNK_SIZE, Base->class.nClasses, 1, tid );

	writeSet( set, Base->class.Class, tid,
		  0, Base->class.nClasses, 1 );

	H5Dclose(set);
	H5Tclose(tid);
	}
    }

static void readClassTable( IOBASE *Base ) {
    hid_t tid, set;

    set = H5Dopen( Base->group_id, FIELD_CLASSES );
    if ( set != H5I_INVALID_HID ) {

	if ( Base->class.setClass_id != H5I_INVALID_HID
		&& Base->class.piClass == NULL ) {
	    Base->class.piClass = (uint8_t *)malloc( CHUNK_SIZE * sizeof(uint8_t) );
	    assert(Base->class.piClass != NULL );
	    }
	tid = makeClassType( H5T_NATIVE_DOUBLE, Base->class.piClass==NULL );
	Base->class.nClasses = getSetSize(set);
	readSet( set, Base->class.Class, tid,
		 0, Base->class.nClasses, 1 );
	H5Dclose(set);
	H5Tclose(tid);
	}

    }

static void addClass( IOBASE *base, PINDEX iOrder, double fMass, double fSoft ) {
    IOCLASS *Class = &base->class;
    uint_fast32_t i;

    /* See if we already have this class: Mass/Softening pair */
    for ( i=0; i<Class->nClasses; i++ ) {
	if ( Class->Class[i].fMass == fMass && Class->Class[i].fSoft == fSoft )
	    break;
	}

    /* Case 1: This is a new class */
    if ( i == Class->nClasses ) {
	assert( Class->nClasses < 256 ); /*TODO: handle this case */
	Class->Class[i].iClass = i;
	Class->Class[i].iOrderStart = iOrder;
	Class->Class[i].fMass = fMass;
	Class->Class[i].fSoft = fSoft;
	Class->Class[i].nCount= 0;
	Class->nClasses++;
	if ( Class->piClass != NULL )
	    Class->piClass[base->nBuffered] = i;
	}

    /* Case 2: This was the last class, and we might be compressing */
    else if ( i == Class->nClasses - 1 && Class->piClass==NULL ) {
	}

    /* Case 3: A match, but a prior class */
    else {
	//createClass(base);
	Class->piClass[base->nBuffered] = i;
	}
    Class->Class[i].nCount++;
    }

static int hdf5GetAttr(
    FIO fio,
    const char *attr, FIO_TYPE dataType, void *data) {
    fioHDF5 *hio = (fioHDF5 *)fio;
    assert(fio->eFormat == FIO_FORMAT_HDF5);


    return 0;
    }

static FIO_SPECIES hdf5Species (struct fioInfo *fio) {
    fioHDF5 *hio = (fioHDF5 *)fio;
    return hio->eCurrent;
    }

static int hdf5Seek(FIO fio,uint64_t iPart,FIO_SPECIES eSpecies) {
    fioHDF5 *hio = (fioHDF5 *)fio;
    IOBASE *base;
    int i;

    assert(fio->eFormat == FIO_FORMAT_HDF5);
    if (eSpecies==FIO_SPECIES_ALL) {
	for( i=1; i<FIO_SPECIES_LAST; i++) {
	    base = &hio->base[i];
	    if (iPart<base->nTotal) {
		eSpecies = i;
		break;
		}
	    iPart -= base->nTotal;
	    }
	assert(eSpecies!=FIO_SPECIES_ALL);
	}
    hio->eCurrent = eSpecies;
    base = &hio->base[eSpecies];
    base->iOffset = iPart;
    base->iIndex = base->nBuffered = 0;
    return 0;
    }

static int hdf5ReadDark(
    FIO fio,uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft,float *pfPot) {
    fioHDF5 *hio = (fioHDF5 *)(fio);
    IOBASE *base = &hio->base[hio->eCurrent];
    int d;
    int i, iClass;

    assert(fio->eFormat == FIO_FORMAT_HDF5);
    assert(hio->eCurrent == FIO_SPECIES_DARK);

    /* If we have exhausted our buffered data, read more */
    if (base->iIndex == base->nBuffered) {
	hsize_t N;

	base->iOffset += base->nBuffered;
	base->iIndex = base->nBuffered = 0;
	if ( base->iOffset >= base->nTotal ) return 0;

	N = base->nTotal - base->iOffset;
	base->nBuffered = N > CHUNK_SIZE ? CHUNK_SIZE : N;

	readSet( base->setR_id, base->R, H5T_NATIVE_DOUBLE,
		 base->iOffset, base->nBuffered, 3 );
	readSet( base->setV_id, base->V, H5T_NATIVE_DOUBLE,
		 base->iOffset, base->nBuffered, 3 );

	if (base->setPot_id != H5I_INVALID_HID) {
	    readSet( base->setPot_id, base->pot,
		     H5T_NATIVE_FLOAT, base->iOffset, base->nBuffered, 1 );
	    }
	if (base->order.setOrder_id != H5I_INVALID_HID) {
	    hid_t dataType = sizeof(PINDEX)==4
		? H5T_NATIVE_UINT32 : H5T_NATIVE_UINT64;
	    readSet( base->order.setOrder_id, base->order.iOrder,
		     dataType, base->iOffset, base->nBuffered, 1 );
	    }

	if (base->class.setClass_id!=H5I_INVALID_HID) {
	    readSet( base->class.setClass_id, base->class.piClass,
		     H5T_NATIVE_UINT8, base->iOffset, base->nBuffered, 1 );
	    }
	}

    /* Position and Velocity are always present */
    for(d=0;d<3;d++) {
	pdPos[d] = base->R[base->iIndex].v[d];
	pdVel[d] = base->V[base->iIndex].v[d];
	}

    /* Potential is optional */
    *pfPot = base->pot ? base->pot[base->iIndex] : 0.0;

    /* iOrder is either sequential, or is listed for each particle */
    if (base->order.setOrder_id != H5I_INVALID_HID) {
	*piOrder = base->order.iOrder[base->iIndex];
	}
    else {
	*piOrder = base->order.iStart + base->iOffset + base->iIndex;
	}

    /* If each particles has a unique class, use that */
    if ( base->class.piClass == NULL ) {
	assert(base->class.nClasses>=1);
	for ( iClass=0; iClass<base->class.nClasses; iClass++ )
	    if ( base->class.Class[iClass].iOrderStart > *piOrder )
		break;
	assert( iClass>0 );
	--iClass;
	}
    /* Otherwise, the particles were sorted by class to save space */
    else {
	iClass = base->class.piClass[base->iIndex];
	assert( iClass < base->class.nClasses );
	}
    *pfMass = base->class.Class[iClass].fMass;
    *pfSoft = base->class.Class[iClass].fSoft;

    /*
    ** Next particle.  If we are at the end of this species,
    ** advance to the start of the next species.
    */
    base->iIndex++;
    if (base->iIndex==base->nTotal) {
	for( i=hio->eCurrent+1; i<FIO_SPECIES_LAST; i++) {
	    base = &hio->base[i];
	    if ( base->nTotal ) {
		hio->eCurrent = i;
		base->iOffset = base->iIndex = base->nBuffered = 0;
		}
	    }
	}
    return 1;
    }

static int hdf5ReadSph(
    FIO fio,uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft, float *pfPot,
    float *pfRho, float *pfTemp, float *pfMetals) {
    return 0;
    }

static int hdf5ReadStar(
    FIO fio,uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft,float *pfPot,
    float *pfMetals, float *pfTform) {
    return 0;
    }

static void hdf5Close(FIO fio) {
    fioHDF5 *hio = (fioHDF5 *)(fio);
    int i;

    for( i=1; i<FIO_SPECIES_LAST; i++) {
	IOBASE *base = hio->base+i;
	if (base->group_id!=H5I_INVALID_HID) {
	    if (base->order.setOrder_id!=H5I_INVALID_HID) {
		free(base->order.iOrder);
		H5Dclose(base->group_id);
		}
	    free(base->V);
	    H5Dclose(base->setV_id);
	    free(base->R);
	    H5Dclose(base->setR_id);
	    H5Gclose(base->group_id);
	    }
	}
    H5Fclose(hio->fileID);
    }

FIO fioHDF5Open(const char *fileName) {
    fioHDF5 *hio;
    H5E_auto_t save_func;
    void *     save_data;
    int i;
    hid_t dataType = sizeof(PINDEX)==4 ? H5T_NATIVE_UINT32 : H5T_NATIVE_UINT64;

    hio = malloc(sizeof(fioHDF5));
    assert(hio!=NULL);
    hio->fio.eFormat = FIO_FORMAT_HDF5;

    for( i=0; i<FIO_SPECIES_LAST; i++)
	hio->fio.nSpecies[i] = 0;

    /* Open the HDF5 file. */
    hio->fileID = H5Fopen(fileName, H5F_ACC_RDONLY, H5P_DEFAULT);
    if ( hio->fileID < 0 ) {
	free(hio);
	return NULL;
	}

    /* Now open all of the available groups.  It's okay if some aren't there. */
    H5Eget_auto(&save_func,&save_data);
    H5Eset_auto(0,0);
    for( i=1; i<FIO_SPECIES_LAST; i++) {
	IOBASE *base = hio->base+i;
	base->group_id = H5Gopen(hio->fileID,fioSpeciesName(i));
	if (base->group_id!=H5I_INVALID_HID) {
	    base->iOffset = 0;
	    base->iIndex = base->nBuffered = 0;
	    base->setR_id = H5Dopen(base->group_id,FIELD_POSITION);
	    assert(base->setR_id!=H5I_INVALID_HID);
	    base->R = malloc(sizeof(*base->R)*CHUNK_SIZE);
	    base->setV_id = H5Dopen(base->group_id,FIELD_VELOCITY);
	    assert(base->setV_id!=H5I_INVALID_HID);
	    base->V = malloc(sizeof(*base->V)*CHUNK_SIZE);
	    base->nTotal = hio->fio.nSpecies[i] = getSetSize(base->setR_id);
	    assert(hio->fio.nSpecies[i] == getSetSize(base->setV_id));

	    base->setPot_id = H5Dopen(base->group_id,FIELD_POTENTIAL);
	    if (base->setPot_id!=H5I_INVALID_HID) {
		base->pot = malloc(sizeof(*base->pot)*CHUNK_SIZE);
		assert(base->pot!=NULL);
		}
	    else {
		base->pot = NULL;
		}

	    base->class.setClass_id = H5Dopen(base->group_id,FIELD_CLASS);
	    if (base->class.setClass_id!=H5I_INVALID_HID) {
		base->class.piClass = malloc(sizeof(*base->class.piClass)*CHUNK_SIZE);
		assert(base->class.piClass!=NULL);
		}
	    readClassTable( base );

	    /* iOrder can have a starting value if they are sequential, or a list */
	    base->order.setOrder_id = H5Dopen(base->group_id,FIELD_ORDER);
	    base->order.iStart = 0;
	    if (base->order.setOrder_id!=H5I_INVALID_HID) {
		base->order.iOrder = malloc(sizeof(*base->order.iOrder)*CHUNK_SIZE);
		assert(base->order.iOrder!=NULL);
		}
	    else {
		base->order.iOrder = NULL;
		readAttribute(base->group_id,ATTR_IORDER,dataType,&base->order.iStart);
		}
	    base->order.iNext = base->order.iStart;
	    }
	else base->nTotal = 0;
	}
    H5Eset_auto(save_func,save_data);

    hio->eCurrent = 0;
    for( i=1; i<FIO_SPECIES_LAST; i++) {
	if (hio->eCurrent==0 && hio->fio.nSpecies[i]) hio->eCurrent=0;
	hio->fio.nSpecies[FIO_SPECIES_ALL] += hio->fio.nSpecies[i];
	}

    hio->fio.fcnClose    = hdf5Close;
    hio->fio.fcnSeek     = hdf5Seek;
    hio->fio.fcnReadDark = hdf5ReadDark;
    hio->fio.fcnReadSph  = hdf5ReadSph;
    hio->fio.fcnReadStar = hdf5ReadStar;
    hio->fio.fcnGetAttr  = hdf5GetAttr;
    hio->fio.fcnSpecies  = hdf5Species;

    return &hio->fio;
    }

#endif

/******************************************************************************\
** GRAFIC FORMAT
\******************************************************************************/

typedef struct {
    int32_t n1, n2, n3;
    float dx, o1, o2, o3;
    float astart, omegam, omegav, H0;
    } GraficHdr4f;

typedef struct {
    int64_t n1, n2, n3;
    float dx, o1, o2, o3;
    float astart, omegam, omegav, H0;
    } GraficHdr8f;

typedef struct {
    int32_t n1, n2, n3;
    double dx, o1, o2, o3;
    double astart, omegam, omegav, H0;
    } GraficHdr4d;

typedef struct {
    int64_t n1, n2, n3;
    double dx, o1, o2, o3;
    double astart, omegam, omegav, H0;
    } GraficHdr8d;

typedef struct {
    FILE *fp;
    union {
	float *pFloat;
	double *pDouble;
	} data;
    int iIndex;
    int nPerSlab;
    int nSlabSize;
    off_t nHdrSize;
    int bDouble;
    GraficHdr8d hdr;
    } graficFile;

typedef struct {
    graficFile fp_velcx;
    graficFile fp_velcy;
    graficFile fp_velcz;
    } graficLevel;

typedef struct {
    struct fioInfo fio;
    uint64_t iOrder;
    double   dTime;
    double   vFactor;
    double   pFactor1;
    double   pFactor2;
    graficLevel level0;
    } fioGrafic;

/*
** Read the header information from a GRAFIC velocity file.
** The data types are inferred from the size of the record.
*/
static void graficReadHdr(graficFile *gf) {
    uint32_t w1,w2;
    GraficHdr4f hdr4f;
    GraficHdr8f hdr8f;
    GraficHdr4d hdr4d;
    GraficHdr8d hdr8d;
    int rc;

    rc = fread(&w1,sizeof(w1),1,gf->fp);
    assert(rc==1);
    switch(w1) {
    case sizeof(GraficHdr4f):
	rc = fread(&hdr4f,w1,1,gf->fp);
	assert(rc==1);
	gf->hdr.n1 = hdr4f.n1;
	gf->hdr.n2 = hdr4f.n2;
	gf->hdr.n3 = hdr4f.n3;
	gf->hdr.dx = hdr4f.dx;
	gf->hdr.o1 = hdr4f.o1;
	gf->hdr.o2 = hdr4f.o2;
	gf->hdr.o3 = hdr4f.o3;
	gf->hdr.astart = hdr4f.astart;
	gf->hdr.omegam = hdr4f.omegam;
	gf->hdr.omegav = hdr4f.omegav;
	gf->hdr.H0 = hdr4f.H0;
	gf->bDouble = 0;
	break;
    case sizeof(GraficHdr8f):
	rc = fread(&hdr8f,w1,1,gf->fp);
	assert(rc==1);
	gf->hdr.n1 = hdr8f.n1;
	gf->hdr.n2 = hdr8f.n2;
	gf->hdr.n3 = hdr8f.n3;
	gf->hdr.dx = hdr8f.dx;
	gf->hdr.o1 = hdr8f.o1;
	gf->hdr.o2 = hdr8f.o2;
	gf->hdr.o3 = hdr8f.o3;
	gf->hdr.astart = hdr8f.astart;
	gf->hdr.omegam = hdr8f.omegam;
	gf->hdr.omegav = hdr8f.omegav;
	gf->hdr.H0 = hdr8f.H0;
	gf->bDouble = 0;
	break;
    case sizeof(GraficHdr4d):
	rc = fread(&hdr4d,w1,1,gf->fp);
	assert(rc==1);
	gf->hdr.n1 = hdr4d.n1;
	gf->hdr.n2 = hdr4d.n2;
	gf->hdr.n3 = hdr4d.n3;
	gf->hdr.dx = hdr4d.dx;
	gf->hdr.o1 = hdr4d.o1;
	gf->hdr.o2 = hdr4d.o2;
	gf->hdr.o3 = hdr4d.o3;
	gf->hdr.astart = hdr4d.astart;
	gf->hdr.omegam = hdr4d.omegam;
	gf->hdr.omegav = hdr4d.omegav;
	gf->hdr.H0 = hdr4d.H0;
	gf->bDouble = 1;
	break;
    case sizeof(GraficHdr8d):
	rc = fread(&hdr8d,w1,1,gf->fp);
	assert(rc==1);
	gf->hdr.n1 = hdr8d.n1;
	gf->hdr.n2 = hdr8d.n2;
	gf->hdr.n3 = hdr8d.n3;
	gf->hdr.dx = hdr8d.dx;
	gf->hdr.o1 = hdr8d.o1;
	gf->hdr.o2 = hdr8d.o2;
	gf->hdr.o3 = hdr8d.o3;
	gf->hdr.astart = hdr8d.astart;
	gf->hdr.omegam = hdr8d.omegam;
	gf->hdr.omegav = hdr8d.omegav;
	gf->hdr.H0 = hdr8d.H0;
	gf->bDouble = 1;
	break;
    default:
	fprintf(stderr,"Invalid GRAFIC header size\n");
	abort();
	}
    rc = fread(&w2,sizeof(w2),1,gf->fp);
    assert(rc==1);
    assert(w1==w2);

    gf->nHdrSize = ftell(gf->fp);

    assert(sizeof(off_t)>=8);
    assert(gf->nHdrSize==w1+2*sizeof(w1));

    gf->nPerSlab = (uint64_t)gf->hdr.n1 * (uint64_t)gf->hdr.n2;
    gf->iIndex = gf->nPerSlab;
    if ( gf->bDouble ) {
	gf->nSlabSize = sizeof(double)*gf->nPerSlab;
	gf->data.pDouble = malloc(gf->nSlabSize);
	assert(gf->data.pDouble!=NULL);
	}
    else {
	gf->nSlabSize = sizeof(float)*gf->nPerSlab;
	gf->data.pFloat = malloc(gf->nSlabSize);
	assert(gf->data.pFloat!=NULL);
	}
    }

/*
** Open a single GRAFIC file and read in the header information
*/
static int graficOpen(graficFile *gf,const char *fileName) {
    gf->fp = fopen(fileName,"rb");
    if ( gf->fp != NULL ) graficReadHdr(gf);
    return gf->fp != NULL;
    }

/*
** Return the next velocity from the file, reading more if necessary.
*/
static double graficRead(graficFile *gf) {
    int rc;
    uint32_t w;
    assert( gf->iIndex <= gf->nPerSlab);
    if ( gf->iIndex == gf->nPerSlab) {
	gf->iIndex = 0;

	/* Read an entire slab */
	rc = fread(&w,sizeof(w),1,gf->fp);
	assert(rc==1 && w==gf->nSlabSize);

	if ( gf->bDouble )
	    rc = fread(gf->data.pDouble,sizeof(double),gf->nPerSlab,gf->fp);
	else
	    rc = fread(gf->data.pFloat,sizeof(float),gf->nPerSlab,gf->fp);
	assert(rc==gf->nPerSlab);

	rc = fread(&w,sizeof(w),1,gf->fp);
	assert(rc==1 && w==gf->nSlabSize);
	}
    if ( gf->bDouble ) return gf->data.pDouble[gf->iIndex++];
    else return gf->data.pFloat[gf->iIndex++];
    }

static void graficSeekFile(graficFile *gf,uint64_t iDark) {
    uint64_t
	sFloat,     /* Bytes in each float (could be double) */
	iSlab,      /* Index of the slab */
	iPart;      /* Index of the particle in the slab */
    uint64_t iByte; /* Byte offset into the file */
    off_t iOffset;
    uint32_t w;
    int rc;

    /* Calculate the slab, particle in slab and byte offset */
    sFloat = gf->bDouble ? sizeof(double) : sizeof(float);
    iSlab = iDark / gf->nPerSlab;
    assert( iSlab < gf->hdr.n3 );
    iPart = iDark - iSlab*gf->nPerSlab;

    iByte = gf->nHdrSize
	+ iSlab * (sFloat*gf->nPerSlab + 2*sizeof(uint32_t))
	+ iPart * sFloat + sizeof(uint32_t);
    iOffset = iByte;
    assert(iOffset==iByte);

    rc = safe_fseek(gf->fp,iByte); assert(rc==0);

    /* Read the remainder of this slab */
    if ( gf->bDouble )
	rc = fread(gf->data.pDouble+iPart,sizeof(double),gf->nPerSlab-iPart,gf->fp);
    else
	rc = fread(gf->data.pFloat+iPart,sizeof(float),gf->nPerSlab-iPart,gf->fp);
    assert(rc==gf->nPerSlab-iPart);

    /* Also verify that the FORTRAN record length is correct */
    rc = fread(&w,sizeof(w),1,gf->fp);
    assert(rc==1 && w==gf->nSlabSize);
    }

/*
** Compare two GRAFIC headers for equality
*/
static int graficCompare(graficFile *a,graficFile *b) {
    return a->hdr.n1 == b->hdr.n1
	&& a->hdr.n2 == b->hdr.n2
	&& a->hdr.n3 == b->hdr.n3
	&& a->hdr.dx == b->hdr.dx
	&& a->hdr.o1 == b->hdr.o1
	&& a->hdr.o2 == b->hdr.o2
	&& a->hdr.o3 == b->hdr.o3
	&& a->hdr.astart == b->hdr.astart
	&& a->hdr.omegam == b->hdr.omegam
	&& a->hdr.omegav == b->hdr.omegav
	&& a->hdr.H0 == b->hdr.H0;
    }

static int graficGetAttr(FIO fio,
    const char *attr, FIO_TYPE dataType, void *data) {
    fioGrafic *gio = (fioGrafic *)fio;
    assert(fio->eFormat == FIO_FORMAT_GRAFIC);
    if ( strcmp(attr,"dTime")==0 ) {
	switch(dataType) {
	case FIO_TYPE_FLOAT: *(float *)(data) = gio->dTime; break;
	case FIO_TYPE_DOUBLE:*(double *)(data) = gio->dTime; break;
	default: return 0;
	    }
	return 1;
	}

    return 0;
    }

static int graficSeek(FIO fio,uint64_t iDark,FIO_SPECIES eSpecies) {
    fioGrafic *gio = (fioGrafic *)fio;

    assert(fio->eFormat == FIO_FORMAT_GRAFIC);
    graficSeekFile(&gio->level0.fp_velcx,iDark);
    graficSeekFile(&gio->level0.fp_velcy,iDark);
    graficSeekFile(&gio->level0.fp_velcz,iDark);
    return 0;
    }

static void graficSetPV(FIO fio,double *r,double *v,double x,double y,double z) {
    fioGrafic *gio = (fioGrafic *)fio;
    r[0] = 0.0 + x * gio->pFactor1 * gio->pFactor2 - 0.5;
    r[1] = 0.0 + y * gio->pFactor1 * gio->pFactor2 - 0.5;
    r[2] = 0.0 + z * gio->pFactor1 * gio->pFactor2 - 0.5;

    v[0] = x * gio->vFactor;
    v[1] = y * gio->vFactor;
    v[2] = z * gio->vFactor;
    }

static int graficReadDark(FIO fio,
		   uint64_t *piOrder,double *pdPos,double *pdVel,
		   float *pfMass,float *pfSoft,float *pfPot) {
    fioGrafic *gio = (fioGrafic *)fio;
    assert(fio->eFormat == FIO_FORMAT_GRAFIC);
    *piOrder = gio->iOrder++;
    graficSetPV(fio,pdPos,pdVel,
		graficRead(&gio->level0.fp_velcx),
		graficRead(&gio->level0.fp_velcy),
		graficRead(&gio->level0.fp_velcz) );
    if ( pfPot) *pfPot = 0.0;
    return 1;
    }

static int graficReadSph(
    FIO fio,uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft,float *pfPot,
    float *pfRho,float *pfTemp, float *pfMetals) {
    fioGrafic *gio = (fioGrafic *)fio;
    assert(fio->eFormat == FIO_FORMAT_GRAFIC);
    *piOrder = gio->iOrder++;
#if 0
    graficSetPV(fio,pdPos,pdVel,
		graficRead(&gio->level0.fp_velbx),
		graficRead(&gio->level0.fp_velby),
		graficRead(&gio->level0.fp_velbz) );
    if ( pfPot) *pfPot = 0.0;
    *pfRho = *pfTemp = *pfSoft = *pfMetals = 0.0;
#endif
    fprintf(stderr,"Reading baryon particles is not supported\n");
    abort();
    }

static void graficClose(FIO fio) {
    fioGrafic *gio = (fioGrafic *)fio;
    if ( gio->level0.fp_velcx.fp!=NULL ) fclose(gio->level0.fp_velcx.fp);
    if ( gio->level0.fp_velcy.fp!=NULL ) fclose(gio->level0.fp_velcy.fp);
    if ( gio->level0.fp_velcz.fp!=NULL ) fclose(gio->level0.fp_velcz.fp);
    free(gio);
    }

typedef struct {
    double omegam;
    double omegav;
    } ddplus_ctx;

static double ddplus(void *ctx,double a) {
    ddplus_ctx *dc = ctx;
    double eta;
    if ( a == 0.0 ) return 0.0;
    eta = sqrt(dc->omegam/a + dc->omegav*a*a + 1.0 - dc->omegam - dc->omegav);
    return 2.5/(eta*eta*eta);
    }

double dRombergO(void *CTX, double (*func)(void *, double), double a,
                 double b, double eps);

static double dplus(double a,double omegam,double omegav) {
    double eta;
    ddplus_ctx ddctx;

    ddctx.omegam = omegam;
    ddctx.omegav = omegav;
    eta = sqrt(omegam/a + omegav*a*a + 1.0 - omegam - omegav);
    return eta/a * dRombergO(&ddctx,ddplus,0,a,1e-8);
}

static double fomega(double a,double omegam,double omegav) {
    double eta, omegak;
    if ( omegam == 1.0 && omegav == 0.0 ) return 1.0;
    omegak=1.0-omegam-omegav;
    eta=sqrt(omegam/a+omegav*a*a+omegak);
    return (2.5/dplus(a,omegam,omegav)-1.5*omegam/a-omegak)/(eta*eta);
}

static double dladt( double a, double omegam, double omegav ) {
    double eta;
    eta=sqrt(omegam/a+omegav*a*a+1.0-omegam-omegav);
    return a * eta;
    }

FIO fioGraficOpen(const char *dirName,int bDouble) {
    fioGrafic *gio;
    struct stat s;
    size_t n;
    int i;
    char *fileName;

    /*
    ** GRAFIC files are found in a specific directory, so verify
    ** that a directory was given as input.
    */
    if ( stat(dirName,&s) == 0 ) {
        if ( !S_ISDIR(s.st_mode) ) {
	    errno = ENOTDIR;
            return NULL;
	    }
	}

    gio = malloc(sizeof(fioGrafic));
    assert(gio!=NULL);
    gio->fio.eFormat = FIO_FORMAT_GRAFIC;
    gio->fio.eMode   = FIO_MODE_READING;
    gio->fio.pszFiles= NULL;
    gio->fio.iFirst = malloc(sizeof(uint64_t));
    assert(gio->fio.iFirst);
    gio->fio.nFiles  = 1;

    gio->fio.fcnClose    = graficClose;
    gio->fio.fcnSeek     = graficSeek;
    gio->fio.fcnReadDark = graficReadDark;
    gio->fio.fcnReadSph  = graficReadSph;
    gio->fio.fcnReadStar = fioNoReadStar;
    gio->fio.fcnGetAttr  = graficGetAttr;
    gio->fio.fcnSpecies  = NULL;

    gio->level0.fp_velcx.fp = gio->level0.fp_velcy.fp = gio->level0.fp_velcz.fp = NULL;

    n = strlen(dirName) + 1;
    fileName = malloc(n + 1 + strlen("ic_velcx"));
    assert(fileName!=NULL);
    strcpy(fileName,dirName);
    strcat(fileName,"/");

    strcpy(fileName+n,"ic_velcx");
    if ( !graficOpen(&gio->level0.fp_velcx,fileName) ) {
	free(fileName);
	graficClose(&gio->fio);
	return NULL;
	}
    strcpy(fileName+n,"ic_velcy");
    if ( !graficOpen(&gio->level0.fp_velcy,fileName) ) {
	free(fileName);
	graficClose(&gio->fio);
	return NULL;
	}
    strcpy(fileName+n,"ic_velcz");
    if ( !graficOpen(&gio->level0.fp_velcz,fileName) ) {
	free(fileName);
	graficClose(&gio->fio);
	return NULL;
	}
    gio->iOrder = 0L;

    assert(graficCompare(&gio->level0.fp_velcx,&gio->level0.fp_velcy));
    assert(graficCompare(&gio->level0.fp_velcx,&gio->level0.fp_velcz));

    gio->fio.nSpecies[FIO_SPECIES_DARK] = (uint64_t)gio->level0.fp_velcx.hdr.n1
	* (uint64_t)gio->level0.fp_velcx.hdr.n2
	* (uint64_t)gio->level0.fp_velcx.hdr.n3;
    gio->fio.nSpecies[FIO_SPECIES_SPH] = gio->fio.nSpecies[FIO_SPECIES_DARK];
    for( i=1; i<FIO_SPECIES_LAST; i++)
	gio->fio.nSpecies[FIO_SPECIES_ALL] += gio->fio.nSpecies[i];
    gio->dTime = gio->level0.fp_velcx.hdr.astart;

    assert(gio->level0.fp_velcx.hdr.n1==gio->level0.fp_velcx.hdr.n2&&gio->level0.fp_velcx.hdr.n2==gio->level0.fp_velcx.hdr.n3);

    /* Makes position dimensionless (i.e., be between 0 and 1) */
    gio->pFactor2 = 1.0 / (gio->level0.fp_velcx.hdr.n1*gio->level0.fp_velcx.hdr.dx);
    gio->pFactor1 = gio->level0.fp_velcx.hdr.astart / (
        fomega(gio->level0.fp_velcx.hdr.astart,gio->level0.fp_velcx.hdr.omegam,gio->level0.fp_velcx.hdr.omegav)
        * gio->level0.fp_velcx.hdr.H0
        * dladt(gio->level0.fp_velcx.hdr.astart,gio->level0.fp_velcx.hdr.omegam,gio->level0.fp_velcx.hdr.omegav) );
    gio->vFactor  = sqrt(8*M_PI/3) * gio->pFactor2 / (gio->level0.fp_velcx.hdr.H0*gio->level0.fp_velcx.hdr.astart);


    free(fileName);
    return &gio->fio;
    }

/******************************************************************************\
** Generic Routines
\******************************************************************************/

/* Attempt to determinate the file type by examining it */
FIO fioOpenMany(int nFiles, const char * const *fileNames) {
    struct stat s;
    const char *fileName = fileNames[0];

    /* The file/directory needs to exist */
    if ( stat(fileName,&s) != 0 ) return NULL;

    /* If given a directory, then it must be a GRAFIC file */
    if ( S_ISDIR(s.st_mode) ) {
	return fioGraficOpen(fileName,/*bDouble*/ 0);
	}

#ifdef USE_HDF5
    else if ( H5Fis_hdf5(fileName) ) {
	return fioHDF5Open(fileName);
	}
#endif

    /* Try tipsy as a last resort */
    else {
	return fioTipsyOpenMany(nFiles,fileNames);
	}

    return NULL;
    }

FIO fioOpen(const char *fileName) {
    return fioOpenMany(1,&fileName);
    }


#ifdef TEST_FIO
/*#include "tipsy.h"*/

#define ID(i,a) ((a)/32.0 + (i))


static int write1(int bStandard, uint64_t nDark, uint64_t nSph, uint64_t nStar) {
    FIO fio;
    int i;

    uint64_t iOrder;
    double dPos[3], dVel[3];
    float fMass, fSoft, fPot, fRho, fTemp, fMetals, fTform;

    printf("Creating %s tipsy test file\n",
	   bStandard ? "standard" : "native" );

    fio = fioTipsyCreate("test1.std",0,bStandard,0.0,nSph,nDark,nStar);
    iOrder = 0;
    for( i=0; i<nSph; i++) {
	dPos[0] = ID(iOrder,1);
	dPos[1] = ID(iOrder,2);
	dPos[2] = ID(iOrder,3);
	dVel[0] = ID(iOrder,4);
	dVel[1] = ID(iOrder,5);
	dVel[2] = ID(iOrder,6);
	fMass   = ID(iOrder,7);
	fSoft   = ID(iOrder,8);
	fPot    = ID(iOrder,9);
	fRho    = ID(iOrder,10);
	fTemp   = ID(iOrder,11);
	fMetals = ID(iOrder,12);
	fioWriteSph(fio,iOrder,dPos,dVel,fMass,fSoft,fPot,fRho,fTemp,fMetals);
	iOrder++;
	}
    for( i=0; i<nDark; i++) {
	dPos[0] = ID(iOrder,1);
	dPos[1] = ID(iOrder,2);
	dPos[2] = ID(iOrder,3);
	dVel[0] = ID(iOrder,4);
	dVel[1] = ID(iOrder,5);
	dVel[2] = ID(iOrder,6);
	fMass   = ID(iOrder,7);
	fSoft   = ID(iOrder,8);
	fPot    = ID(iOrder,9);
	fioWriteDark(fio,iOrder,dPos,dVel,fMass,fSoft,fPot);
	iOrder++;
	}
    for( i=0; i<nStar; i++) {
	dPos[0] = ID(iOrder,1);
	dPos[1] = ID(iOrder,2);
	dPos[2] = ID(iOrder,3);
	dVel[0] = ID(iOrder,4);
	dVel[1] = ID(iOrder,5);
	dVel[2] = ID(iOrder,6);
	fMass   = ID(iOrder,7);
	fSoft   = ID(iOrder,8);
	fPot    = ID(iOrder,9);
	fMetals = ID(iOrder,12);
	fTform  = ID(iOrder,13);
	fioWriteStar(fio,iOrder,dPos,dVel,fMass,fSoft,fPot,fMetals,fTform);
	iOrder++;
	}


    fioClose(fio);
    }


static int write2(int bStandard, uint64_t nDark, int nSph, int nStar) {
    FIO fio;
    int i,j;

    uint64_t iOrder, N;
    double dPos[3], dVel[3];
    float fMass, fSoft, fPot, fRho, fTemp, fMetals, fTform;
    int L[2];
    N = nDark + nSph + nStar;

    i = 0;
    L[0] = N/2;
    L[1] = N;

    for(j=1; j<=2; j++) {
	char fname[100];
	sprintf(fname,"test%d.std",j);

	printf("Creating %s tipsy test file : %s\n",
	       bStandard ? "standard" : "native", fname );

	fio = fioTipsyCreatePart(fname,0,0,bStandard, 0.0,
				 nSph, nDark, nStar, i);
	for( ; i<L[j-1]; i++) {
	    iOrder = i;
	    dPos[0] = ID(iOrder,1);
	    dPos[1] = ID(iOrder,2);
	    dPos[2] = ID(iOrder,3);
	    dVel[0] = ID(iOrder,4);
	    dVel[1] = ID(iOrder,5);
	    dVel[2] = ID(iOrder,6);
	    fMass   = ID(iOrder,7);
	    fSoft   = ID(iOrder,8);
	    fPot    = ID(iOrder,9);
	    fRho    = ID(iOrder,10);
	    fTemp   = ID(iOrder,11);
	    fMetals = ID(iOrder,12);
	    fTform  = ID(iOrder,13);

	    if ( i<nSph ) {
		fioWriteSph(fio,iOrder,dPos,dVel,fMass,fSoft,fPot,fRho,fTemp,fMetals);
		}
	    else if ( i<nSph+nDark ) {
		fioWriteDark(fio,iOrder,dPos,dVel,fMass,fSoft,fPot);
		}
	    else if ( i<nSph+nDark+nStar ) {
		fioWriteStar(fio,iOrder,dPos,dVel,fMass,fSoft,fPot,fMetals,fTform);
		}
	    else assert(0);
	    }
	fioClose(fio);
	}
    }

static int read1(FIO fio,uint64_t nDark) {
    int i;

    uint64_t iOrder, N;
    int rc;
    double dPos[3], dVel[3];
    float fMass, fSoft, fPot, fRho, fTemp, fMetals, fTform;

    printf("Reading sequentially\n");

    N = fioGetN(fio,FIO_SPECIES_ALL);
    assert(nDark == fioGetN(fio,FIO_SPECIES_DARK));
    for( i=0; i<N; i++) {
	switch(fioSpecies(fio)) {
	case FIO_SPECIES_DARK:
	    rc = fioReadDark(fio,&iOrder,dPos,dVel,&fMass,&fSoft,&fPot);
	    if ( rc!=1) {
		printf("%lu\n",i);
		}
	    assert(rc==1);
	    assert(iOrder == i);
	    assert(dPos[0] == ID(iOrder,1));
	    assert(dPos[1] == ID(iOrder,2));
	    assert(dPos[2] == ID(iOrder,3));
	    assert(dVel[0] == ID(iOrder,4));
	    assert(dVel[1] == ID(iOrder,5));
	    assert(dVel[2] == ID(iOrder,6));
	    assert(fMass == ID(iOrder,7));
	    assert(fSoft == ID(iOrder,8));
	    assert(fPot  == ID(iOrder,9));
	    break;
	case FIO_SPECIES_SPH:
	    rc = fioReadSph(fio,&iOrder,dPos,dVel,&fMass,&fSoft,&fPot,&fRho,&fTemp,&fMetals);
	    assert(rc==1);
	    assert(iOrder == i);
	    assert(dPos[0] == ID(iOrder,1));
	    assert(dPos[1] == ID(iOrder,2));
	    assert(dPos[2] == ID(iOrder,3));
	    assert(dVel[0] == ID(iOrder,4));
	    assert(dVel[1] == ID(iOrder,5));
	    assert(dVel[2] == ID(iOrder,6));
	    assert(fMass   == ID(iOrder,7));
	    assert(fSoft   == ID(iOrder,8));
	    assert(fPot    == ID(iOrder,9));
	    assert(fRho    == ID(iOrder,10));
	    assert(fTemp   == ID(iOrder,11));
	    assert(fMetals == ID(iOrder,12));
	    break;
	case FIO_SPECIES_STAR:
	    rc = fioReadStar(fio,&iOrder,dPos,dVel,&fMass,&fSoft,&fPot,&fMetals,&fTform);
	    assert(rc==1);
	    assert(iOrder == i);
	    assert(dPos[0] == ID(iOrder,1));
	    assert(dPos[1] == ID(iOrder,2));
	    assert(dPos[2] == ID(iOrder,3));
	    assert(dVel[0] == ID(iOrder,4));
	    assert(dVel[1] == ID(iOrder,5));
	    assert(dVel[2] == ID(iOrder,6));
	    assert(fMass   == ID(iOrder,7));
	    assert(fSoft   == ID(iOrder,8));
	    assert(fPot    == ID(iOrder,9));
	    assert(fMetals == ID(iOrder,12));
	    assert(fTform  == ID(iOrder,13));
	    break;
	    }
	}

    }

static int read2(FIO fio,uint64_t nDark) {
    int i;

    uint64_t iOrder, N;
    int rc;
    double dPos[3], dVel[3];
    float fMass, fSoft, fPot, fRho, fTemp, fMetals, fTform;

    printf("Reading in reverse order\n");

    N = fioGetN(fio,FIO_SPECIES_ALL);
    assert(nDark == fioGetN(fio,FIO_SPECIES_DARK));

    for( i=0; i<N; i++) {
	iOrder = N-i-1;
	fioSeek(fio,iOrder,FIO_SPECIES_ALL);

	switch(fioSpecies(fio)) {
	case FIO_SPECIES_DARK:
	    rc = fioReadDark(fio,&iOrder,dPos,dVel,&fMass,&fSoft,&fPot);
	    assert(rc==1);
	    assert(iOrder == N-i-1);
	    assert(dPos[0] == ID(iOrder,1));
	    assert(dPos[1] == ID(iOrder,2));
	    assert(dPos[2] == ID(iOrder,3));
	    assert(dVel[0] == ID(iOrder,4));
	    assert(dVel[1] == ID(iOrder,5));
	    assert(dVel[2] == ID(iOrder,6));
	    assert(fMass == ID(iOrder,7));
	    assert(fSoft == ID(iOrder,8));
	    assert(fPot  == ID(iOrder,9));
	    break;
	case FIO_SPECIES_SPH:
	    rc = fioReadSph(fio,&iOrder,dPos,dVel,&fMass,&fSoft,&fPot,&fRho,&fTemp,&fMetals);
	    assert(rc==1);
	    assert(iOrder == N-i-1);
	    assert(dPos[0] == ID(iOrder,1));
	    assert(dPos[1] == ID(iOrder,2));
	    assert(dPos[2] == ID(iOrder,3));
	    assert(dVel[0] == ID(iOrder,4));
	    assert(dVel[1] == ID(iOrder,5));
	    assert(dVel[2] == ID(iOrder,6));
	    assert(fMass   == ID(iOrder,7));
	    assert(fSoft   == ID(iOrder,8));
	    assert(fPot    == ID(iOrder,9));
	    assert(fRho    == ID(iOrder,10));
	    assert(fTemp   == ID(iOrder,11));
	    assert(fMetals == ID(iOrder,12));
	    break;
	case FIO_SPECIES_STAR:
	    rc = fioReadStar(fio,&iOrder,dPos,dVel,&fMass,&fSoft,&fPot,&fMetals,&fTform);
	    assert(rc==1);
	    assert(iOrder == N-i-1);
	    assert(dPos[0] == ID(iOrder,1));
	    assert(dPos[1] == ID(iOrder,2));
	    assert(dPos[2] == ID(iOrder,3));
	    assert(dVel[0] == ID(iOrder,4));
	    assert(dVel[1] == ID(iOrder,5));
	    assert(dVel[2] == ID(iOrder,6));
	    assert(fMass   == ID(iOrder,7));
	    assert(fSoft   == ID(iOrder,8));
	    assert(fPot    == ID(iOrder,9));
	    assert(fMetals == ID(iOrder,12));
	    assert(fTform  == ID(iOrder,13));
	    break;
	    }
	}
    }

static int read3(FIO fio,uint64_t iOrderBase,uint64_t n,FIO_SPECIES eSpecies) {
    int i;
    uint64_t iOrder1, iOrder2, N;
    int rc;
    double dPos[3], dVel[3];
    float fMass, fSoft, fPot, fRho, fTemp, fMetals, fTform;

    printf("Reading randomly\n");

    N = fioGetN(fio,eSpecies);
    if (n!=N ) {
	printf("n=%lu N=%lu species=%d\n", n, N, eSpecies);
	}
    assert(n==N);

    for( i=0; i<N; i++) {
	/* Randomly around the current position (faster as we stay in buffer) */
	iOrder1 = (i + (rand()%1000)) % N;
	fioSeek(fio,iOrder1,eSpecies);

	switch(fioSpecies(fio)) {
	case FIO_SPECIES_DARK:
	    rc = fioReadDark(fio,&iOrder2,dPos,dVel,&fMass,&fSoft,&fPot);
	    assert(rc==1);
	    assert(iOrder1+iOrderBase == iOrder2);
	    assert(dPos[0] == ID(iOrder2,1));
	    assert(dPos[1] == ID(iOrder2,2));
	    assert(dPos[2] == ID(iOrder2,3));
	    assert(dVel[0] == ID(iOrder2,4));
	    assert(dVel[1] == ID(iOrder2,5));
	    assert(dVel[2] == ID(iOrder2,6));
	    assert(fMass == ID(iOrder2,7));
	    assert(fSoft == ID(iOrder2,8));
	    assert(fPot  == ID(iOrder2,9));
	    break;
	case FIO_SPECIES_SPH:
	    rc = fioReadSph(fio,&iOrder2,dPos,dVel,&fMass,&fSoft,&fPot,&fRho,&fTemp,&fMetals);
	    assert(rc==1);
	    assert(iOrder1+iOrderBase == iOrder2);
	    assert(dPos[0] == ID(iOrder2,1));
	    assert(dPos[1] == ID(iOrder2,2));
	    assert(dPos[2] == ID(iOrder2,3));
	    assert(dVel[0] == ID(iOrder2,4));
	    assert(dVel[1] == ID(iOrder2,5));
	    assert(dVel[2] == ID(iOrder2,6));
	    assert(fMass   == ID(iOrder2,7));
	    assert(fSoft   == ID(iOrder2,8));
	    assert(fPot    == ID(iOrder2,9));
	    assert(fRho    == ID(iOrder2,10));
	    assert(fTemp   == ID(iOrder2,11));
	    assert(fMetals == ID(iOrder2,12));
	    break;
	case FIO_SPECIES_STAR:
	    rc = fioReadStar(fio,&iOrder2,dPos,dVel,&fMass,&fSoft,&fPot,&fMetals,&fTform);
	    assert(rc==1);
	    assert(iOrder1+iOrderBase == iOrder2);
	    assert(dPos[0] == ID(iOrder2,1));
	    assert(dPos[1] == ID(iOrder2,2));
	    assert(dPos[2] == ID(iOrder2,3));
	    assert(dVel[0] == ID(iOrder2,4));
	    assert(dVel[1] == ID(iOrder2,5));
	    assert(dVel[2] == ID(iOrder2,6));
	    assert(fMass   == ID(iOrder2,7));
	    assert(fSoft   == ID(iOrder2,8));
	    assert(fPot    == ID(iOrder2,9));
	    assert(fMetals == ID(iOrder2,12));
	    assert(fTform  == ID(iOrder2,13));
	    break;
	    }
	}
    }

int main(int argc, char *argv[]) {
    uint64_t nDark = 100000;
    uint64_t nSph  = 45001;
    uint64_t nStar = 100;
    FIO fio;
    const char *szFiles[] = {"test1.std","test2.std"};


#if 1
    write1(1,nDark,nSph,nStar);
    fio = fioOpen("test1.std");
    read1(fio,nDark);
    fioClose(fio);
    fio = fioOpen("test1.std");
    read2(fio,nDark);
    fioClose(fio);
    fio = fioOpen("test1.std");
    read3(fio,0,nDark+nSph+nStar,FIO_SPECIES_ALL);
    fioClose(fio);
    fio = fioOpen("test1.std");
    read3(fio,nSph,nDark,FIO_SPECIES_DARK);
    fioClose(fio);
    fio = fioOpen("test1.std");
    read3(fio,0,nSph,FIO_SPECIES_SPH);
    fioClose(fio);
    fio = fioOpen("test1.std");
    read3(fio,nSph+nDark,nStar,FIO_SPECIES_STAR);
    fioClose(fio);

    write1(0,nDark,nSph,nStar);
    fio = fioOpen("test1.std");
    read1(fio,nDark);
    fioClose(fio);
    fio = fioOpen("test1.std");
    read2(fio,nDark);
    fioClose(fio);
#endif
#if 1
    write2(1,nDark,nSph,nStar);
#endif
#if 1
    fio = fioOpenMany(2,szFiles);
    read1(fio,nDark);
    fioClose(fio);
#endif
#if 1
    fio = fioOpenMany(2,szFiles);
    read2(fio,nDark);
    fioClose(fio);

    fio = fioOpenMany(2,szFiles);
    read3(fio,nSph,nDark,FIO_SPECIES_DARK);
    fioClose(fio);
    fio = fioOpenMany(2,szFiles);
    read3(fio,0,nSph,FIO_SPECIES_SPH);
    fioClose(fio);
    fio = fioOpenMany(2,szFiles);
    read3(fio,nSph+nDark,nStar,FIO_SPECIES_STAR);
    fioClose(fio);
#endif
#if 0
    readP1(nDark);
#endif

    return 0;
    }
#endif
