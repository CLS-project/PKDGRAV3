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

#include "fio.h"

/*
** This uses the best available seek routine to move to the specified
** 64-bit offset into a file.  The new "fseeko" function is preferred,
** but we fallback to multiple fseek calls if required.
*/
static int safe_fseek(FILE *fp,uint64_t lStart) {
    int iErr;
#ifdef HAVE_FSEEKO
    uint64_t lSeeked;
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
    return iErr;
#else
    /* fseek fails for offsets >= 2**31; this is an ugly workaround */
    static const off_t MAX_OFFSET = 2147483640;

    if (lStart > MAX_OFFSET) {
	iErr = fseek(fp,0,SEEK_SET);
	if (iErr) {
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

/******************************************************************************\
** TIPSY FORMAT
\******************************************************************************/

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
    float phi;
    float eps;
    } tipsyDark;

typedef struct {
    float mass;
    float pos[3];
    float vel[3];
    float phi;
    float rho;
    float temp;
    float hsmooth;
    float metals;
    } tipsySph;

typedef struct {
    float mass;
    float pos[3];
    float vel[3];
    float phi;
    float eps;
    float metals;
    float tform;
    } tipsyStar;

typedef struct {
    struct fioInfo fio;
    FILE *fp;
    XDR xdr;
    off_t nHdrSize;
    uint64_t iOrder;
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

    if ( strcmp(attr,"dTime")==0 ) {
	switch(dataType) {
	case FIO_TYPE_FLOAT: *(float *)(data) = tio->fio.dTime; break;
	case FIO_TYPE_DOUBLE:*(double *)(data) = tio->fio.dTime; break;
	default: return 0;
	    }
	return 1;
	}

    return 0;
    }

static int tipsyReadNativeDark(FIO fio,
    uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft,float *pfPot) {
    fioTipsy *tio = (fioTipsy *)fio;
    tipsyDark dark;
    int rc;
    int d;

    rc = fread(&dark,sizeof(dark),1,tio->fp);
    if (rc!=1) return 0;
    *piOrder = tio->iOrder++;
    *pfMass = dark.mass;
    for(d=0;d<3;d++) {
	pdPos[d] = dark.pos[d];
	pdVel[d] = dark.vel[d];
	}
    *pfSoft = dark.eps;
    *pfPot = dark.phi;
    return 1;
    }

static int tipsyReadNativeDarkDbl(FIO fio,
    uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft,float *pfPot) {
    fioTipsy *tio = (fioTipsy *)fio;
    int rc;
    int d;
    int fTmp[3];

    *piOrder = tio->iOrder++;
    rc = fread(pfMass,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fread(pdPos,sizeof(double),3,tio->fp); if (rc!=3) return 0;
    rc = fread(fTmp,sizeof(float),3,tio->fp); if (rc!=3) return 0;
    for(d=0;d<3;d++) pdVel[d] = fTmp[d];
    rc = fread(pfPot,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    rc = fread(pfSoft,sizeof(float),1,tio->fp); if (rc!=1) return 0;
    return 1;
    }

static int tipsyReadStandardDark(FIO fio,
    uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft,float *pfPot) {
    fioTipsy *tio = (fioTipsy *)fio;
    int d;
    float fTmp;

    *piOrder = tio->iOrder++;
    if (!xdr_float(&tio->xdr,pfMass)) return 0;
    for(d=0;d<3;d++) {
	if (!xdr_float(&tio->xdr,&fTmp)) return 0;
	pdPos[d] = fTmp;
	}
    for(d=0;d<3;d++) {
	if (!xdr_float(&tio->xdr,&fTmp)) return 0;
	pdVel[d] = fTmp;
	}
    if (!xdr_float(&tio->xdr,pfSoft)) return 0;
    if (!xdr_float(&tio->xdr,pfPot)) return 0;
    return 1;
    }

static int tipsyReadStandardDarkDbl(FIO fio,
    uint64_t *piOrder,double *pdPos,double *pdVel,
    float *pfMass,float *pfSoft,float *pfPot) {
    fioTipsy *tio = (fioTipsy *)fio;
    int d;
    float fTmp;

    *piOrder = tio->iOrder++;
    if (!xdr_float(&tio->xdr,pfMass)) return 0;
    for(d=0;d<3;d++) {
	if (!xdr_double(&tio->xdr,pdPos+d)) return 0;
	}
    for(d=0;d<3;d++) {
	if (!xdr_float(&tio->xdr,&fTmp)) return 0;
	pdVel[d] = fTmp;
	}
    if (!xdr_float(&tio->xdr,pfSoft)) return 0;
    if (!xdr_float(&tio->xdr,pfPot)) return 0;
    return 1;
    }

static void tipsyCloseNative(FIO fio) {
    fioTipsy *tio = (fioTipsy *)fio;
    fclose(tio->fp);
    free(tio);
    }

static void tipsyCloseStandard(FIO fio) {
    fioTipsy *tio = (fioTipsy *)fio;
    xdr_destroy(&tio->xdr);
    tipsyCloseNative(fio);
    }

static int tipsySeek(fioTipsy *tio,uint64_t iPart,FIO_SPECIES eSpecies,
    int nSphSize, int nDarkSize, int nStarSize ) {
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
    if (iPart<=tio->fio.nSpecies[FIO_SPECIES_SPH]) {
	iByte += iPart * nSphSize;
	}
    else {
	iByte += tio->fio.nSpecies[FIO_SPECIES_SPH] * nSphSize;
	iPart -= tio->fio.nSpecies[FIO_SPECIES_SPH];

	if (iPart<=tio->fio.nSpecies[FIO_SPECIES_DARK]) {
	    iByte += iPart * nDarkSize;
	    }
	else {
	    iByte += tio->fio.nSpecies[FIO_SPECIES_DARK] * nDarkSize;
	    iPart -= tio->fio.nSpecies[FIO_SPECIES_DARK];

	    if (iPart<=tio->fio.nSpecies[FIO_SPECIES_STAR]) {
		iByte += iPart * nStarSize;
		}
	    else {
		errno = ESPIPE;
		return -1;
		}
	    }
	}
    iOffset = iByte;
    assert(iOffset==iByte);
    rc = safe_fseek(tio->fp,iByte); assert(rc==0);
    return rc;
    }

static int tipsySeekNative(FIO fio,uint64_t iPart,FIO_SPECIES eSpecies) {
    fioTipsy *tio = (fioTipsy *)fio;
    return tipsySeek(tio,iPart,eSpecies,
	sizeof(tipsySph),sizeof(tipsyDark),sizeof(tipsyStar));
    }

static int tipsySeekNativeDbl(FIO fio,uint64_t iPart,FIO_SPECIES eSpecies) {
    fioTipsy *tio = (fioTipsy *)fio;
    return tipsySeek(tio,iPart,eSpecies,
	sizeof(tipsySph)+3*sizeof(float),
	sizeof(tipsyDark)+3*sizeof(float),
	sizeof(tipsyStar)+3*sizeof(float));
    }

static int tipsySeekStandard(FIO fio,uint64_t iPart,FIO_SPECIES eSpecies) {
    fioTipsy *tio = (fioTipsy *)fio;
    int rc;
    xdr_destroy(&tio->xdr);
    rc = tipsySeekNative(fio,iPart,eSpecies);
    xdrstdio_create(&tio->xdr,tio->fp,XDR_DECODE);
    return rc;
    }

static int tipsySeekStandardDbl(FIO fio,uint64_t iPart,FIO_SPECIES eSpecies) {
    fioTipsy *tio = (fioTipsy *)fio;
    int rc;
    xdr_destroy(&tio->xdr);
    rc = tipsySeekNativeDbl(fio,iPart,eSpecies);
    xdrstdio_create(&tio->xdr,tio->fp,XDR_DECODE);
    return rc;
    }

FIO fioTipsyOpen(const char *fileName,int bDouble) {
    fioTipsy *tio;
    tipsyHdr h;
    int rc,i;

    tio = malloc(sizeof(fioTipsy));
    assert(tio!=NULL);
    tio->iOrder = 0;

    tio->fp = fopen(fileName,"r");
    if (tio->fp==NULL) {
	free(tio);
	return NULL;
	}
    rc = fread(&h,sizeof(h),1,tio->fp);
    if (rc!=1) {
	fprintf(stderr,"Error reading Tipsy header\n");
	fclose(tio->fp);
	free(tio);
	return NULL;
	}

    tio->fio.fcnGetAttr = tipsyGetAttr;

    /* Check for native binary format */
    if (h.nDim>=1 && h.nDim<=3) {
	if ( h.nDark + h.nSph + h.nStar != h.nBodies ) {
	    fprintf(stderr,"Tipsy Header mismatch:"
		    " nDim=%u nDark=%u nSph=%u nStar=%u nBodies=%u\n",
		    h.nDim, h.nDark, h.nSph, h.nStar, h.nBodies);
	    fclose(tio->fp);
	    free(tio);
	    return NULL;
	    }
	tio->fio.fcnClose    = tipsyCloseNative;
	tio->fio.fcnSeek     = bDouble?tipsySeekNativeDbl:tipsySeekNative;
	tio->fio.fcnReadDark = bDouble?tipsyReadNativeDarkDbl:tipsyReadNativeDark;
	}

    /* Okay, this must be "standard" format */
    else {
	rewind(tio->fp);
	xdrstdio_create(&tio->xdr,tio->fp,XDR_DECODE);
	rc = xdrHeader(&tio->xdr,&h);
	if (rc!=1) {
	    perror("Error reading tipsy header");
	    fclose(tio->fp);
	    free(tio);
	    return NULL;
	    }

	if (h.nDim<1 || h.nDim>3 || h.nDark + h.nSph + h.nStar != h.nBodies ) {
	    fprintf(stderr,"Tipsy Header mismatch:"
		    " nDim=%u nDark=%u nSph=%u nStar=%u nBodies=%u\n",
		    h.nDim, h.nDark, h.nSph, h.nStar, h.nBodies);
	    xdr_destroy(&tio->xdr);
	    fclose(tio->fp);
	    free(tio);
	    return NULL;
	    }
	tio->fio.fcnClose    = tipsyCloseStandard;
	tio->fio.fcnSeek     = bDouble ? tipsySeekStandardDbl : tipsySeekStandard;
	tio->fio.fcnReadDark = bDouble ? tipsyReadStandardDarkDbl : tipsyReadStandardDark;
	}

    tio->nHdrSize = sizeof(h);
    tio->fio.dTime = h.dTime;
    for( i=0; i<FIO_SPECIES_LAST; i++)
	tio->fio.nSpecies[i] = tio->fio.oSpecies[i] = 0;
    tio->fio.nSpecies[FIO_SPECIES_DARK] = h.nDark;
    tio->fio.nSpecies[FIO_SPECIES_SPH]  = h.nSph;
    tio->fio.nSpecies[FIO_SPECIES_STAR] = h.nStar;
    for( i=1; i<FIO_SPECIES_LAST; i++)
	tio->fio.nSpecies[FIO_SPECIES_ALL] += tio->fio.nSpecies[i];
    tio->fio.oSpecies[FIO_SPECIES_DARK] = tio->fio.oSpecies[FIO_SPECIES_SPH]
	+ tio->fio.nSpecies[FIO_SPECIES_SPH];
    tio->fio.oSpecies[FIO_SPECIES_STAR] = tio->fio.oSpecies[FIO_SPECIES_DARK]
	+ tio->fio.nSpecies[FIO_SPECIES_DARK];
    tio->fio.oSpecies[FIO_SPECIES_LAST] = tio->fio.oSpecies[FIO_SPECIES_STAR]
	+ tio->fio.nSpecies[FIO_SPECIES_STAR];
    return &tio->fio;
    }


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
    struct fioInfo fio;
    uint64_t iOrder;
    double   vFactor;
    double   pFactor1;
    double   pFactor2;
    graficFile fp_velcx;
    graficFile fp_velcy;
    graficFile fp_velcz;
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
	assert(0);
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
	nFloat,     /* Bytes in each float (could be double) */
	iSlab,      /* Index of the slab */
	iPart;      /* Index of the particle in the slab */
    uint64_t iByte; /* Byte offset into the file */
    off_t iOffset;
    uint32_t w;
    int rc;

    /* Calculate the slab, particle in slab and byte offset */
    nFloat = gf->bDouble ? sizeof(double) : sizeof(float);
    iSlab = iDark / gf->nPerSlab;
    assert( iSlab < gf->hdr.n3 );
    iPart = iDark - iSlab*gf->nPerSlab;

    iByte = gf->nHdrSize
	+ iSlab * (nFloat*gf->nPerSlab + 2*sizeof(uint32_t))
	+ iPart * nFloat + sizeof(uint32_t);
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
    fioGrafic *gfio = (fioGrafic *)fio;
    if ( strcmp(attr,"dTime")==0 ) {
	switch(dataType) {
	case FIO_TYPE_FLOAT: *(float *)(data) = gfio->fio.dTime; break;
	case FIO_TYPE_DOUBLE:*(double *)(data) = gfio->fio.dTime; break;
	default: return 0;
	    }
	return 1;
	}

    return 0;
    }

static int graficSeek(FIO fio,uint64_t iDark,FIO_SPECIES eSpecies) {
    fioGrafic *gfio = (fioGrafic *)fio;

    graficSeekFile(&gfio->fp_velcx,iDark);
    graficSeekFile(&gfio->fp_velcy,iDark);
    graficSeekFile(&gfio->fp_velcz,iDark);
    return 0;
    }

static int graficReadDark(FIO fio,
		   uint64_t *piOrder,double *pdPos,double *pdVel,
		   float *pfMass,float *pfSoft,float *pfPot) {
    double v[3];
    fioGrafic *gfio = (fioGrafic *)fio;
    *piOrder = gfio->iOrder++;

    v[0] = graficRead(&gfio->fp_velcx);
    v[1] = graficRead(&gfio->fp_velcy);
    v[2] = graficRead(&gfio->fp_velcz);

    pdVel[0] = v[0] * gfio->vFactor;
    pdVel[1] = v[1] * gfio->vFactor;
    pdVel[2] = v[2] * gfio->vFactor;


    if ( pfPot) *pfPot = 0.0;
    return 1;
    }

static void graficClose(FIO fio) {
    fioGrafic *gfio = (fioGrafic *)fio;
    if ( gfio->fp_velcx.fp!=NULL ) fclose(gfio->fp_velcx.fp);
    if ( gfio->fp_velcy.fp!=NULL ) fclose(gfio->fp_velcy.fp);
    if ( gfio->fp_velcz.fp!=NULL ) fclose(gfio->fp_velcz.fp);
    free(gfio);
    }

FIO fioGraficOpen(const char *dirName,int bDouble) {
    fioGrafic *gfio;
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

    gfio = malloc(sizeof(fioGrafic));
    assert(gfio!=NULL);

    gfio->fio.fcnClose    = graficClose;
    gfio->fio.fcnSeek     = graficSeek;
    gfio->fio.fcnReadDark = graficReadDark;
    gfio->fio.fcnGetAttr  = graficGetAttr;

    gfio->fp_velcx.fp = gfio->fp_velcy.fp = gfio->fp_velcz.fp = NULL;

    n = strlen(dirName) + 1;
    fileName = malloc(n + 1 + strlen("ic_velcx"));
    assert(fileName!=NULL);
    strcpy(fileName,dirName);
    strcat(fileName,"/");

    strcpy(fileName+n,"ic_velcx");
    if ( !graficOpen(&gfio->fp_velcx,fileName) ) {
	free(fileName);
	graficClose(&gfio->fio);
	return NULL;
	}
    strcpy(fileName+n,"ic_velcy");
    if ( !graficOpen(&gfio->fp_velcy,fileName) ) {
	free(fileName);
	graficClose(&gfio->fio);
	return NULL;
	}
    strcpy(fileName+n,"ic_velcz");
    if ( !graficOpen(&gfio->fp_velcz,fileName) ) {
	free(fileName);
	graficClose(&gfio->fio);
	return NULL;
	}
    gfio->iOrder = 0L;

    assert(graficCompare(&gfio->fp_velcx,&gfio->fp_velcy));
    assert(graficCompare(&gfio->fp_velcx,&gfio->fp_velcz));

    for( i=0; i<FIO_SPECIES_LAST; i++)
	gfio->fio.nSpecies[i] = gfio->fio.oSpecies[i] = 0;
    gfio->fio.nSpecies[FIO_SPECIES_DARK] = (uint64_t)gfio->fp_velcx.hdr.n1
	* (uint64_t)gfio->fp_velcx.hdr.n2
	* (uint64_t)gfio->fp_velcx.hdr.n3;
    gfio->fio.nSpecies[FIO_SPECIES_SPH] = gfio->fio.nSpecies[FIO_SPECIES_DARK];
    for( i=1; i<FIO_SPECIES_LAST; i++)
	gfio->fio.nSpecies[FIO_SPECIES_ALL] += gfio->fio.nSpecies[i];
    gfio->fio.oSpecies[FIO_SPECIES_SPH] = gfio->fio.oSpecies[FIO_SPECIES_DARK]
	+ gfio->fio.nSpecies[FIO_SPECIES_DARK];
    gfio->fio.dTime = gfio->fp_velcx.hdr.astart;

    assert(gfio->fp_velcx.hdr.n1==gfio->fp_velcx.hdr.n2&&gfio->fp_velcx.hdr.n2==gfio->fp_velcx.hdr.n3);

    /* Makes position dimensionless (i.e., be between 0 and 1) */
    gfio->pFactor2 = 1.0 / (gfio->fp_velcx.hdr.n1*gfio->fp_velcx.hdr.dx);
#if 0
    gfio->pFactor1 = gfio->fp_velcx.hdr.astart / (
        fomega(gfio->fp_velcx.hdr.astart,gfio->fp_velcx.hdr.omegam,gfio->fp_velcx.hdr.omegav)
        * gfio->fp_velcx.hdr.H0
        * dladt(gfio->fp_velcx.hdr.astart,gfio->fp_velcx.hdr.omegam,gfio->fp_velcx.hdr.omegav) );
#endif
    gfio->vFactor  = sqrt(8*M_PI/3) * gfio->pFactor2 / (gfio->fp_velcx.hdr.H0*gfio->fp_velcx.hdr.astart);


    free(fileName);
    return &gfio->fio;
    }

/******************************************************************************\
** Generic Routines
\******************************************************************************/

/* Attempt to determinate the file type by examining it */
FIO fioOpen(const char *fileName,int bDouble) {
    struct stat s;

    /* The file/directory needs to exist */
    if ( stat(fileName,&s) != 0 ) return NULL;

    /* If given a directory, then it must be a GRAFIC file */
    if ( S_ISDIR(s.st_mode) ) {
	return fioGraficOpen(fileName,bDouble);
	}

    /* Try tipsy as a last resort */
    else {
	return fioTipsyOpen(fileName,bDouble);
	}

    return NULL;
    }


#ifdef TEST_FIO
#include "tipsy.h"

int main(int argc, char *argv[]) {
    FIO fio;
    uint64_t iOrder;
    double pos[3], vel[3];
    double dTime1, dTime2;
    float mass,soft,pot;
    TCTX ctx;
    struct dark_particle *p;
    int iType;
    double dSoft;
    int i, d;
    uint64_t N, N2;
    int iFail;

    assert(argc==2);

    fio = fioOpen(argv[1],0);
    if ( fio == NULL ) {
	perror(argv[1]);
	exit(errno);
	}
    N = fioGetN(fio,FIO_SPECIES_ALL);
    fioGetAttr(fio,"dTime",FIO_TYPE_DOUBLE,&dTime1);

    TipsyInitialize(&ctx,0,argv[1]);
    N2 = iTipsyNumParticles(ctx);
    if ( N != N2 ) {
	printf("Particle counts do not match: %lu != %lu\n", N, N2);
	}
    dTime2 = dTipsyTime(ctx);
    if ( dTime1 != dTime2 ) {
	printf("Simulation times do not match: %g != %g\n",
	    dTime1, dTime2);
	}

    for(i=0;i<N;i++) {
	fioReadDark(fio,&iOrder,pos,vel,&mass,&soft,&pot);
	p = (struct dark_particle *)pTipsyRead(ctx,&iType,&dSoft);
	iFail = 0;
	for(d=0;d<3;d++) {
	    if (p->pos[d]!=pos[d]) iFail++;
	    if (p->vel[d]!=vel[d]) iFail++;
	    }
	if (iOrder!=i) iFail++;
	if (mass!=p->mass) iFail++;
	if (soft!=p->eps) iFail++;
	if (pot!=p->phi) iFail++;
	if ( iFail ) {
	    printf("NO:\n");
	    }
	}

    TipsyFinish(ctx);
    fioClose(fio);

    return 0;
    }
#endif
