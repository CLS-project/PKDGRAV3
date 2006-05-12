#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "pkd.h"
#include "outtype.h"
#include "floattype.h"
#include "tipsydefs.h"
#include <rpc/types.h>
#include <rpc/xdr.h>
#ifdef LARGEF
#include <fcntl.h>
#endif

FLOAT ArrType(PARTICLE *p,int iType)
{
	switch (iType) {
	case OUT_DENSITY_ARRAY:
		return(p->fDensity);
	case OUT_COLOR_ARRAY:
	case OUT_POT_ARRAY:
		return(p->fPot);
	case OUT_AMAG_ARRAY:
		return(sqrt(p->a[0]*p->a[0] + p->a[1]*p->a[1] + p->a[2]*p->a[2]));
	case OUT_RUNG_ARRAY:
		return(p->iRung);
	case OUT_DT_ARRAY:
		return(p->dt);
	case OUT_SOFT_ARRAY:
	        return(p->fSoft);
	case OUT_GROUP_ARRAY:
		return (p->pGroup);
#ifdef RELAXATION
	case OUT_RELAX_ARRAY:
		 return (p->fRelax);
#endif	
	    default:
		return(0.0);
		}
	}


FLOAT VecType(PARTICLE *p,int iDim,int iType)
{
	switch (iType) {
	case OUT_POS_VECTOR:
		return(p->r[iDim]);
	case OUT_VEL_VECTOR:
		return(p->v[iDim]);
	case OUT_ACCEL_VECTOR:
		return(p->a[iDim]);
	default:
		return(0.0);
		}
	}


void pkdOutArray(PKD pkd,char *pszFileName,int iArrType)
{
	FILE *fp;
	FLOAT fOut;
	int i;

	fp = fopen (pszFileName,"a");		
	assert(fp != NULL);
	/*
	 ** Write Array Elements!
	 */
	for (i=0;i<pkd->nLocal;++i) {
		fOut = ArrType(&pkd->pStore[i],iArrType);
		fprintf(fp,"%.8g\n",fOut);
		}
	i = fclose(fp);
	if (i != 0) {
		perror("pkdOutArray: could not close file");
		exit(1);
		}
	}


void pkdOutVector(PKD pkd,char *pszFileName,int iDim,int iVecType)
{
	FILE *fp;
	FLOAT fOut;
	int i;
	
	fp = fopen (pszFileName,"a");		
	assert(fp != NULL);
	/*
	 ** Write Vector Elements!
	 */
	for (i=0;i<pkd->nLocal;++i) {
		fOut = VecType(&pkd->pStore[i],iDim,iVecType);
		fprintf(fp,"%.8g\n",fOut);
		}
	i = fclose(fp);
	if (i != 0) {
		perror("pkdOutVector: could not close file");
		exit(1);
		}
	}

void pkdOutGroup(PKD pkd,char *pszFileName,int iType, int nStart,double dvFac)
{
	FILE *fp;
	int i,j,nout,lStart;
	
	/*
	 ** Write Group Data!
	 */
	if(iType == OUT_GROUP_STATS){
		fp = fopen(pszFileName,"r+");
		assert(fp != NULL);
		for (i=0;i<pkd->nGroups;++i) {
		  fprintf(fp,"%d ",pkd->groupData[i].iGlobalId);
		  fprintf(fp,"%d ",pkd->groupData[i].nTotal);
		  fprintf(fp,"%.4g ",pkd->groupData[i].fMass);
		  fprintf(fp,"%.4g ",pkd->groupData[i].fGasMass);
		  fprintf(fp,"%.4g ",pkd->groupData[i].fStarMass);			  
		  fprintf(fp,"%.4g ",pkd->groupData[i].fAvgDens);
		  fprintf(fp,"%.4g ",pkd->groupData[i].fRadius);
		  fprintf(fp,"%.4g ",pkd->groupData[i].fDeltaR2);
		  fprintf(fp,"%.4g ",dvFac*pkd->groupData[i].fVelDisp);
		  fprintf(fp,"%.4g ",dvFac*dvFac*pkd->groupData[i].fVelSigma2[0]);
		  fprintf(fp,"%.4g ",dvFac*dvFac*pkd->groupData[i].fVelSigma2[1]);
		  fprintf(fp,"%.4g ",dvFac*dvFac*pkd->groupData[i].fVelSigma2[2]);
		  fprintf(fp,"%.4g ",pkd->groupData[i].r[0]);
		  fprintf(fp,"%.4g ",pkd->groupData[i].r[1]);
		  fprintf(fp,"%.4g ",pkd->groupData[i].r[2]);
		  fprintf(fp,"%.4g ",pkd->groupData[i].rpotmin[0]);
		  fprintf(fp,"%.4g ",pkd->groupData[i].rpotmin[1]);
		  fprintf(fp,"%.4g ",pkd->groupData[i].rpotmin[2]);
		  fprintf(fp,"%.4g ",dvFac*pkd->groupData[i].v[0]);
		  fprintf(fp,"%.4g ",dvFac*pkd->groupData[i].v[1]);
		  fprintf(fp,"%.4g ",dvFac*pkd->groupData[i].v[2]);
		  fprintf(fp,"%.4g ",dvFac*pkd->groupData[i].vcircMax);
		  fprintf(fp,"%.4g ",pkd->groupData[i].rvcircMax);
		  fprintf(fp,"%.4g ",pkd->groupData[i].rvir);
		  fprintf(fp,"%.4g ",pkd->groupData[i].Mvir);
		  fprintf(fp,"%.4g ",pkd->groupData[i].lambda);
		  fprintf(fp,"\n");
		}
	} else if(iType == OUT_GROUP_TIPSY_NAT){ 

		struct star_particle sp;
		fp = fopen(pszFileName,"r+");
		assert(fp != NULL);
		/*
		 ** Seek past the header
		 */	
	    lStart = sizeof(struct dump);
		fseek(fp,lStart,SEEK_SET);

		for (i=0;i<pkd->nGroups;++i) {
			if(pkd->groupData[i].bMyGroup){
				for (j=0;j<3;++j) {
					sp.pos[j] = pkd->groupData[i].r[j];
					sp.vel[j] = dvFac*pkd->groupData[i].v[j];
					}
				sp.mass = pkd->groupData[i].fMass;
				sp.eps = pkd->groupData[i].fRadius;
				sp.tform = 0.0;
				sp.metals = 0.0;
				nout = fwrite(&sp,sizeof(struct star_particle),1,fp);
				mdlassert(pkd->mdl,nout == 1);
				}
			}
	} else if(iType == OUT_GROUP_TIPSY_STD){

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
			if(pkd->groupData[i].bMyGroup){
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
 				fTmp = pkd->groupData[i].fRadius; 
 				xdr_float(&xdrs,&fTmp); /* softening eps*/ 
 				fTmp = 0.0; 
 				xdr_float(&xdrs,&fTmp); /* grav. potential phi*/ 
				}
			}
		xdr_destroy(&xdrs); 
  
	} else if(iType == OUT_GROUP_PROFILES){
	 fp = fopen(pszFileName,"at");
		assert(fp != NULL);
		for (i=0;i< pkd->nBins;++i) {
 				fprintf(fp,"%d ",pkd->groupBin[i].iId);
				fprintf(fp,"%.4g ",pkd->groupBin[i].fRadius);
				fprintf(fp,"%d ",pkd->groupBin[i].nMembers);
				fprintf(fp,"%.4g ",pkd->groupBin[i].fDensity);	
				fprintf(fp,"%.4g ",pkd->groupBin[i].fMassEnclosed);		
				fprintf(fp,"%.4g ",
						pow(pkd->groupBin[i].fMassEnclosed/pkd->groupBin[i].fRadius,0.5) );
				fprintf(fp,"%.4g ",dvFac*dvFac*pkd->groupBin[i].v2[0]);
				fprintf(fp,"%.4g ",dvFac*dvFac*pkd->groupBin[i].v2[1]);
				fprintf(fp,"%.4g ",dvFac*dvFac*pkd->groupBin[i].v2[2]);
				fprintf(fp,"%.4g ",dvFac*dvFac*pkd->groupBin[i].L[0]);
				fprintf(fp,"%.4g ",dvFac*dvFac*pkd->groupBin[i].L[1]);
				fprintf(fp,"%.4g ",dvFac*dvFac*pkd->groupBin[i].L[2]);
				fprintf(fp,"%.4g ",pkd->groupBin[i].a);
				fprintf(fp,"%.4g ",pkd->groupBin[i].b);
				fprintf(fp,"%.4g ",pkd->groupBin[i].c);
				fprintf(fp,"%.4g ",pkd->groupBin[i].phi);
				fprintf(fp,"%.4g ",pkd->groupBin[i].theta);
				fprintf(fp,"%.4g ",pkd->groupBin[i].psi);
				fprintf(fp,"\n");
		}	
	}	
	i = fclose(fp);
	if (i != 0) {
		perror("pkdOutGroup: could not close file");
		exit(1);
		}
	}
