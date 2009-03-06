#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <inttypes.h>
#include "pkd.h"
#include "outtype.h"
#ifndef HAVE_CONFIG_H
#include "floattype.h"
#endif
#include "tipsydefs.h"
#include <rpc/types.h>
#include <rpc/xdr.h>
#ifdef LARGEF
#include <fcntl.h>
#endif

const char *outtype_c_module_id = "$Id$";
const char *outtype_h_module_id = OUTTYPE_H_MODULE_ID;

/* Write an integer */
static uint64_t fetchInteger(PKD pkd,PARTICLE *p,int iType,int iDim) {
    uint64_t v;

    switch (iType) {
    case OUT_IORDER_ARRAY:
	v = p->iOrder;
	break;
    case OUT_GROUP_ARRAY:
	assert(pkd->oGroup);
	v = *pkdInt32(p,pkd->oGroup);
	break;
    default:
	v = 0;
	}
    return v;
    }

static void writeInteger(PKD pkd,void *ctxVoid,PARTICLE *p,int iType,int iDim) {
    PKDOUT ctx = (PKDOUT)ctxVoid;
    fprintf(ctx->fp,"%"PRIu64"\n",fetchInteger(pkd,p,iType,iDim));
    }

#ifdef HAVE_LIBBZ2
static void writeIntegerBZ2(PKD pkd, void *ctxVoid,PARTICLE *p,int iType,int iDim) {
    PKDOUT ctx = (PKDOUT)ctxVoid;
    char buffer[100];
    int bzerror;

    sprintf(buffer,"%"PRIu64"\n",fetchInteger(pkd,p,iType,iDim));
    BZ2_bzWrite(&bzerror,ctx->CTX.bz,buffer,strlen(buffer));
    }
#endif

#ifdef HAVE_LIBZ
static void writeIntegerZ(PKD pkd, void *ctxVoid,PARTICLE *p,int iType,int iDim) {
    PKDOUT ctx = (PKDOUT)ctxVoid;
    gzprintf(ctx->CTX.gz,"%"PRIu64"\n",fetchInteger(pkd,p,iType,iDim));
    }
#endif

static void writeHdr(PKD pkd,void *ctxVoid,uint64_t N) {
    PKDOUT ctx = (PKDOUT)ctxVoid;
    fprintf(ctx->fp,"%"PRIu64"\n",N);
    }


#ifdef HAVE_LIBBZ2
static void writeHdrBZ2(PKD pkd, void *ctxVoid,uint64_t N) {
    PKDOUT ctx = (PKDOUT)ctxVoid;
    char buffer[100];
    int bzerror;
    sprintf(buffer,"%"PRIu64"\n",N);
    BZ2_bzWrite(&bzerror,ctx->CTX.bz,buffer,strlen(buffer));
    }
#endif

#ifdef HAVE_LIBZ
static void writeHdrZ(PKD pkd, void *ctxVoid,uint64_t N) {
    PKDOUT ctx = (PKDOUT)ctxVoid;
    gzprintf(ctx->CTX.gz,"%"PRIu64"\n",N);
    }
#endif

static double fetchFloat(PKD pkd,PARTICLE *p,int iType,int iDim) {
    float *a;
    double v;
    VELSMOOTH *pvel;
    switch (iType) {
    case OUT_DENSITY_ARRAY:
	v = p->fDensity;
	break;
    case OUT_COLOR_ARRAY:
	assert(0);
    case OUT_POT_ARRAY:
	assert(pkd->oPotential);
	a = pkdPot(pkd,p);
	v = *a;
	break;
    case OUT_AMAG_ARRAY:
	assert(pkd->oAcceleration); /* Validate memory model */
	a = pkdAccel(pkd,p);
	v = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
	break;
    case OUT_RUNG_ARRAY:
	v = p->uRung;
	break;
    case OUT_SOFT_ARRAY:
	v = pkdSoft(pkd,p);
	break;
    case OUT_RELAX_ARRAY:
	assert(pkd->oRelaxation);
	a = pkdField(p,pkd->oRelaxation);
	v = *a;
	break;
    case OUT_DIVV_ARRAY:
	assert(pkd->oVelSmooth); /* Validate memory model */
	pvel = pkdField(p,pkd->oVelSmooth);
	v = pvel->divv;
	break;
    case OUT_VELDISP2_ARRAY:
	assert(pkd->oVelSmooth); /* Validate memory model */
	pvel = pkdField(p,pkd->oVelSmooth);
	v = pvel->veldisp2;
    case OUT_VELDISP_ARRAY:
	assert(pkd->oVelSmooth); /* Validate memory model */
	pvel = pkdField(p,pkd->oVelSmooth);
	v = sqrt(pvel->veldisp2);
	break;
    case OUT_PHASEDENS_ARRAY:
	assert(pkd->oVelSmooth); /* Validate memory model */
	pvel = pkdField(p,pkd->oVelSmooth);
	v = p->fDensity*pow(pvel->veldisp2,-1.5);
	break;
    case OUT_POS_VECTOR:
	v = p->r[iDim];
	break;
    case OUT_VEL_VECTOR:
	assert(pkd->oVelocity); /* Validate memory model */
	v = pkdVel(pkd,p)[iDim];
	break;
    case OUT_MEANVEL_VECTOR:
	assert(pkd->oVelSmooth); /* Validate memory model */
	pvel = pkdField(p,pkd->oVelSmooth);
	v = pvel->vmean[iDim];
	break;
    case OUT_ACCEL_VECTOR:
	assert(pkd->oAcceleration); /* Validate memory model */
	a = pkdAccel(pkd,p);
	v = a[iDim];
	break;
    default:
	v = 0.0;
	}
    return v;
    }
static void writeFloat(PKD pkd,void *ctxVoid,PARTICLE *p,int iType,int iDim) {
    PKDOUT ctx = (PKDOUT)ctxVoid;
    fprintf(ctx->fp,"%.8g\n",fetchFloat(pkd,p,iType,iDim));
    }
#ifdef HAVE_LIBBZ2
static void writeFloatBZ2(PKD pkd, void *ctxVoid,PARTICLE *p,int iType,int iDim) {
    PKDOUT ctx = (PKDOUT)ctxVoid;
    char buffer[100];
    int bzerror;

    sprintf(buffer,"%.8g\n",fetchFloat(pkd,p,iType,iDim));
    BZ2_bzWrite(&bzerror,ctx->CTX.bz,buffer,strlen(buffer));
    }
#endif

#ifdef HAVE_LIBZ
static void writeFloatZ(PKD pkd, void *ctxVoid,PARTICLE *p,int iType,int iDim) {
    PKDOUT ctx = (PKDOUT)ctxVoid;
    gzprintf(ctx->CTX.gz,"%.8g\n",fetchFloat(pkd,p,iType,iDim));
    }
#endif

static void closeFILE(PKD pkd, void *ctxVoid ) {
    PKDOUT ctx = (PKDOUT)ctxVoid;
    fclose(ctx->fp);
    }

static void closeBZ2(PKD pkd, void *ctxVoid ) {
    PKDOUT ctx = (PKDOUT)ctxVoid;
    int bzerror;
    unsigned int nbytes_in_lo32, nbytes_in_hi32,
	nbytes_out_lo32, nbytes_out_hi32;

    BZ2_bzWriteClose64( &bzerror, ctx->CTX.bz, 0, 
			&nbytes_in_lo32, &nbytes_in_hi32,
			&nbytes_out_lo32, &nbytes_out_hi32 );
    fclose(ctx->fp);
    }

static void closeZ(PKD pkd, void *ctxVoid ) {
    PKDOUT ctx = (PKDOUT)ctxVoid;

    if (gzclose(ctx->CTX.gz) != 0) {
	perror("pkdOutASCII: could not close file");
	exit(1);
	}
    }

PKDOUT pkdOpenOutASCII(PKD pkd,char *pszFileName,const char *mode,int iType) {
    PKDOUT ctx;
    void (*fnOut)(PKD pkd,void *fp,PARTICLE *p,int iType,int iDim);
#if defined(HAVE_LIBZ) || defined(HAVE_LIBBZ2)
    size_t n;
#ifdef HAVE_LIBBZ2
    void (*fnOutBZ2)(PKD pkd,BZFILE *fp,PARTICLE *p,int iType,int iDim);
#endif
#ifdef HAVE_LIBZ
    void (*fnOutZ)(PKD pkd,gzFile fp,PARTICLE *p,int iType,int iDim);
#endif
#endif

    ctx = malloc(sizeof(struct pkdout)); assert(ctx!=NULL);
    ctx->fp = NULL;

    switch(iType) {
    case OUT_IORDER_ARRAY:
    case OUT_GROUP_ARRAY:
	fnOut = writeInteger;
#ifdef HAVE_LIBBZ2
	fnOutBZ2 = writeIntegerBZ2;
#endif	
#ifdef HAVE_LIBZ
	fnOutZ = writeIntegerZ;
#endif	
	break;

    case OUT_DENSITY_ARRAY:
    case OUT_COLOR_ARRAY:
    case OUT_POT_ARRAY:
    case OUT_AMAG_ARRAY:
    case OUT_RUNG_ARRAY:
    case OUT_SOFT_ARRAY:
    case OUT_RELAX_ARRAY:
    case OUT_DIVV_ARRAY:
    case OUT_VELDISP2_ARRAY:
    case OUT_VELDISP_ARRAY:
    case OUT_PHASEDENS_ARRAY:

    case OUT_POS_VECTOR:
    case OUT_VEL_VECTOR:
    case OUT_MEANVEL_VECTOR:
    case OUT_ACCEL_VECTOR:
	fnOut = writeFloat;
#ifdef HAVE_LIBBZ2
	fnOutBZ2 = writeFloatBZ2;
#endif
#ifdef HAVE_LIBZ
	fnOutZ = writeFloatZ;
#endif
	break;

    default:
	assert(0);
	}

#if defined(HAVE_LIBZ) || defined(HAVE_LIBBZ2)
    n = strlen(pszFileName);
#ifdef HAVE_LIBBZ2
    if ( n>4 && strcmp(pszFileName+n-4,".bz2")==0 ) {
	int bzerror;
	ctx->fp = fopen (pszFileName,mode);
	assert(ctx->fp != NULL);
	ctx->CTX.bz = BZ2_bzWriteOpen( &bzerror, ctx->fp, 9, 0, 0 );
	ctx->fnOut = fnOutBZ2;
	ctx->fnHdr = writeHdrBZ2;
	ctx->fnClose = closeBZ2;
	return ctx;
	}
#endif
#ifdef HAVE_LIBZ
    if ( n>3 && strcmp(pszFileName+n-3,".gz")==0 ) {
	ctx->CTX.gz = gzopen(pszFileName,mode);
	assert(ctx->CTX.gz!=NULL);
	ctx->fnOut = fnOutZ;
	ctx->fnHdr = writeHdrZ;
	ctx->fnClose = closeZ;
	return ctx;
	}
#endif

#endif

    ctx->fp = fopen (pszFileName,mode);
    assert(ctx->fp != NULL);
    ctx->fnOut = fnOut;
    ctx->fnHdr = writeHdr;
    ctx->fnClose = closeFILE;
    return ctx;
    }

void pkdCloseOutASCII(PKD pkd,PKDOUT ctx) {
    ctx->fnClose(pkd,ctx);
    free(ctx);
    }

void pkdOutHdr(PKD pkd,PKDOUT ctx,uint64_t N) {
    (*ctx->fnHdr)(pkd,ctx,N);
    }

void pkdOutASCII(PKD pkd,PKDOUT ctx,int iType,int iDim) {
    int i;

    /*
    ** Write Elements!
    */
    for (i=0;i<pkd->nLocal;++i) {
	PARTICLE *p = pkdParticle(pkd,i);
	if ( pkdIsSrcActive(p,0,MAX_RUNG) )
	    (*ctx->fnOut)(pkd,ctx,p,iType,iDim);
	}
    }

#ifdef USE_HDF5
void pkdOutHDF5(PKD pkd,char *pszFileName,int iType,int iDim) {
    assert(0);
    }
#endif

void pkdOutGroup(PKD pkd,char *pszFileName,int iType, int nStart,double dvFac) {
    FILE *fp;
    int i,j,nout,lStart;

    /*
     ** Write Group Data!
     */
    if (iType == OUT_GROUP_STATS) {
	fp = fopen(pszFileName,"r+");
	assert(fp != NULL);
	for (i=0;i<pkd->nGroups;++i) {
	    fprintf(fp,"%d ",pkd->groupData[i].iGlobalId);
	    fprintf(fp,"%d ",pkd->groupData[i].nTotal);
	    fprintf(fp,"%.8g ",pkd->groupData[i].fMass);
            fprintf(fp,"%.8g ",pkd->groupData[i].fRMSRadius);
	    fprintf(fp,"%.8g ",pkd->groupData[i].r[0]);
	    fprintf(fp,"%.8g ",pkd->groupData[i].r[1]);
	    fprintf(fp,"%.8g ",pkd->groupData[i].r[2]);
	    fprintf(fp,"%.8g ",dvFac*pkd->groupData[i].v[0]);
	    fprintf(fp,"%.8g ",dvFac*pkd->groupData[i].v[1]);
	    fprintf(fp,"%.8g ",dvFac*pkd->groupData[i].v[2]);
	    fprintf(fp,"\n");
	    }
	}
    else if (iType == OUT_GROUP_TIPSY_NAT) {

	struct star_particle sp;
	fp = fopen(pszFileName,"r+");
	assert(fp != NULL);
	/*
	 ** Seek past the header
	 */
	lStart = sizeof(struct dump);
	fseek(fp,lStart,SEEK_SET);

	for (i=0;i<pkd->nGroups;++i) {
	    if (pkd->groupData[i].bMyGroup) {
		for (j=0;j<3;++j) {
		  sp.pos[j] = pkd->groupData[i].r[j];
		  sp.vel[j] = dvFac*pkd->groupData[i].v[j];
		}
		sp.mass = pkd->groupData[i].fMass;
		sp.eps = pkd->groupData[i].fRMSRadius;
		sp.tform = 0.0;
		sp.metals = 0.0;
		nout = fwrite(&sp,sizeof(struct star_particle),1,fp);
		mdlassert(pkd->mdl,nout == 1);
		}
	    }
	}
    else if (iType == OUT_GROUP_TIPSY_STD) {

	float fTmp;
	XDR xdrs;

	fp = fopen(pszFileName,"r+");
	assert(fp != NULL);
	/*
	 ** Seek past the header
	 */
	lStart = 32;
	lStart += nStart*44;

	fseek(fp,lStart,SEEK_SET);

	xdrstdio_create(&xdrs,fp,XDR_ENCODE);
	for (i=0;i<pkd->nGroups;++i) {
	    if (pkd->groupData[i].bMyGroup) {
		fTmp = pkd->groupData[i].fMass;
		xdr_float(&xdrs,&fTmp);
		for (j=0;j<3;++j) {
		    fTmp = pkd->groupData[i].r[j];
		    xdr_float(&xdrs,&fTmp);
		    }
		for (j=0;j<3;++j) {
		    fTmp = dvFac*pkd->groupData[i].v[j];
		    xdr_float(&xdrs,&fTmp);
		    }
		fTmp = 0.0;
		xdr_float(&xdrs,&fTmp); /* metals */
		xdr_float(&xdrs,&fTmp); /* t form */
		fTmp = pkd->groupData[i].fRMSRadius;
		xdr_float(&xdrs,&fTmp); /* softening eps*/
		fTmp = 0.0;
		xdr_float(&xdrs,&fTmp); /* grav. potential phi*/
		}
	    }
	xdr_destroy(&xdrs);

	}
    else if (iType == OUT_GROUP_PROFILES) {
	fp = fopen(pszFileName,"at");
	assert(fp != NULL);
	for (i=0;i< pkd->nBins;++i) {
	    fprintf(fp,"%.8g ",pkd->groupBin[i].fRadius);
	    fprintf(fp,"%d ",pkd->groupBin[i].nMembers);
	    fprintf(fp,"%.8g ",pkd->groupBin[i].fMassInBin);
	    fprintf(fp,"\n");
	    }
	}
    else assert(0);
    i = fclose(fp);
    if (i != 0) {
	perror("pkdOutGroup: could not close file");
	exit(1);
	}
    }
