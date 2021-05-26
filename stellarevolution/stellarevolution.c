#ifdef STELLAR_EVOLUTION


#include <math.h>
#include <hdf5.h>

#include "stellarevolution.h"
#include "hydro/hydro.h"


void msrStellarEvolutionInit(MSR msr) {
   assert(STEV_N_ELEM == ELEMENT_COUNT);

   int i;
   char achPath[280];
   STEV_RAWDATA *CCSNdata, *AGBdata, *SNIadata, *LifetimesData;


   sprintf(achPath, "%s/CCSN.hdf5", msr->param.achStEvolPath);
   char *pszCCSNFields[] = {"Z_0.0004", "Z_0.004", "Z_0.008", "Z_0.02", "Z_0.05"};
   int nCCSNFields = 5;
   CCSNdata = stevReadTable(achPath, pszCCSNFields, nCCSNFields);
   assert(CCSNdata->nZs == STEV_CCSN_N_METALLICITY);
   assert(CCSNdata->nMasses == STEV_CCSN_N_MASS);
   assert(CCSNdata->nSpecs == STEV_N_ELEM);


   sprintf(achPath, "%s/AGB.hdf5", msr->param.achStEvolPath);
   char *pszAGBFields[] = {"Z_0.004", "Z_0.008", "Z_0.019"};
   int nAGBFields = 3;
   AGBdata = stevReadTable(achPath, pszAGBFields, nAGBFields);
   assert(AGBdata->nZs == STEV_AGB_N_METALLICITY);
   assert(AGBdata->nMasses == STEV_AGB_N_MASS);
   assert(AGBdata->nSpecs == STEV_N_ELEM);


   sprintf(achPath, "%s/SNIa.hdf5", msr->param.achStEvolPath);
   SNIadata = stevReadSNIaTable(achPath);
   assert(SNIadata->nSpecs == STEV_N_ELEM);


   sprintf(achPath, "%s/Lifetimes.hdf5", msr->param.achStEvolPath);
   LifetimesData = stevReadLifetimesTable(achPath);
   assert(LifetimesData->nZs == STEV_LIFETIMES_N_METALLICITY);
   assert(LifetimesData->nMasses == STEV_LIFETIMES_N_MASS);


   for (i = 0; i < CCSNdata->nMasses; i++)
      CCSNdata->pfMasses[i] = log10(CCSNdata->pfMasses[i]);
   for (i = 0; i < AGBdata->nMasses; i++)
      AGBdata->pfMasses[i] = log10(AGBdata->pfMasses[i]);


   STEV_DATA buffer;
   /* The lowest value of the initial mass array is chosen such that the corresponding 
      lifetime is greater than the age of the Universe, for all metallicities. The IMF
      is still normalized in the range [msr->param.dIMF_MinMass,msr->param.dIMF_MaxMass] */
   double dMinMass = 0.7;
   stevPopulateInitialMassArrays(dMinMass, msr->param.dIMF_MaxMass,
				 STEV_INTERP_N_MASS, buffer.afMasses, buffer.afLogMasses);

   if (strcmp(msr->param.achIMFtype, "chabrier") == 0) {
      stevChabrierIMF(buffer.afMasses, STEV_INTERP_N_MASS, msr->param.dIMF_MinMass,
		      msr->param.dIMF_MaxMass, buffer.afIMF);
   }
   else {
      printf("ERROR: Undefined IMF type has been given in achIMFtype parameter: %s\n",
	     msr->param.achIMFtype);
      assert(0);
   }

   if (dMinMass == msr->param.dIMF_MinMass) {
      float IMFnorm = 0.0f;
      for (i = 0; i < STEV_INTERP_N_MASS - 1; i++) {
	 IMFnorm += (buffer.afMasses[i + 1] * buffer.afIMF[i + 1] + buffer.afMasses[i] *
		     buffer.afIMF[i]) * (buffer.afMasses[i + 1] - buffer.afMasses[i]);
      }
      IMFnorm *= 0.5f;
      printf("IMF normalization is: %.4f\n", IMFnorm);
      if (fabs(IMFnorm - 1.0) > 1e-2) {
	 printf("ERROR: IMF normalization differs from unity by more than 1%%\n");
	 assert(0);
      }
   }
   else {
      printf("WARNING: IMF normalization has not been verified. Initial mass\n"
	     "         array's range is different from the IMF normalization range\n");
   }

   /* Yields array in the buffer is ordered as (Metallicities, Masses, Elements). MetalYield 
      and EjectedMass arrays, on the other hand, are ordered as (Metallicities, Masses). 
      To change this, the arguments of every call to stevRowMajorIndex that set the variable
      idxBuffer must be swapped appropriately. */
   int j, k, idxMass, idxTable, idxBuffer;
   float fDeltaMass, fMass, fMtrans;
   fMtrans = msr->param.fCCSN_MinMass;
   for (i = 0; i < STEV_INTERP_N_MASS; i++) {
      fMass = buffer.afMasses[i];
      if (fMass < fMtrans) {
	 stevGetIndex1D(AGBdata->pfMasses, STEV_AGB_N_MASS, log10(fMass), &idxMass, &fDeltaMass);

	 for (j = 0; j < STEV_N_ELEM; j++) {
	    for (k = 0; k < STEV_CCSN_N_METALLICITY; k++) {
	       idxBuffer = stevRowMajorIndex(k, i, j, STEV_CCSN_N_METALLICITY,
					     STEV_INTERP_N_MASS, STEV_N_ELEM);

	       buffer.afCCSN_Yields[idxBuffer] = 0.0f;
	    }

	    for (k = 0; k < STEV_AGB_N_METALLICITY; k++) {
	       idxTable  = stevRowMajorIndex(k, j, idxMass, STEV_AGB_N_METALLICITY, STEV_N_ELEM,
					     STEV_AGB_N_MASS);
	       idxBuffer = stevRowMajorIndex(k, i, j, STEV_AGB_N_METALLICITY,
					     STEV_INTERP_N_MASS, STEV_N_ELEM);

	       /* WARNING: The following formula is correct only when mass is the last index
		  in the tables */
	       buffer.afAGB_Yields[idxBuffer] =
		  AGBdata->pfYields[idxTable] * (1.0 - fDeltaMass) +
		  AGBdata->pfYields[idxTable + 1] * fDeltaMass;
	    }
	 }

	 for (k = 0; k < STEV_CCSN_N_METALLICITY; k++) {
	    idxBuffer = stevRowMajorIndex(k, i, 0, STEV_CCSN_N_METALLICITY,
					  STEV_INTERP_N_MASS, 1);

	    buffer.afCCSN_MetalYield[idxBuffer] = 0.0f;
	    buffer.afCCSN_EjectedMass[idxBuffer] = 0.0f;
	 }

	 for (k = 0; k < STEV_AGB_N_METALLICITY; k++) {
	    idxTable  = stevRowMajorIndex(k, idxMass, 0, STEV_AGB_N_METALLICITY,
					  STEV_AGB_N_MASS, 1);
	    idxBuffer = stevRowMajorIndex(k, i, 0, STEV_AGB_N_METALLICITY,
					  STEV_INTERP_N_MASS, 1);

	    /* WARNING: The following formulas are correct only when mass is the last index
	       in the tables */
	    buffer.afAGB_MetalYield[idxBuffer] =
	       AGBdata->pfMetalYield[idxTable] * (1.0 - fDeltaMass) +
	       AGBdata->pfMetalYield[idxTable + 1] * fDeltaMass;

	    buffer.afAGB_EjectedMass[idxBuffer] =
	       AGBdata->pfEjectedMass[idxTable] * (1.0 - fDeltaMass) +
	       AGBdata->pfEjectedMass[idxTable + 1] * fDeltaMass;
	 }
      }
      else {
	 stevGetIndex1D(CCSNdata->pfMasses, STEV_CCSN_N_MASS, log10(fMass),
			&idxMass, &fDeltaMass);

	 for (j = 0; j < STEV_N_ELEM; j++) {
	    for (k = 0; k < STEV_CCSN_N_METALLICITY; k++) {
	       idxTable  = stevRowMajorIndex(k, j, idxMass, STEV_CCSN_N_METALLICITY,
					     STEV_N_ELEM, STEV_CCSN_N_MASS);
	       idxBuffer = stevRowMajorIndex(k, i, j, STEV_CCSN_N_METALLICITY,
					     STEV_INTERP_N_MASS, STEV_N_ELEM);

	       /* WARNING: The following formula is correct only when mass is the last index
		  in the tables */
	       buffer.afCCSN_Yields[idxBuffer] =
		  CCSNdata->pfYields[idxTable] * (1.0 - fDeltaMass) +
		  CCSNdata->pfYields[idxTable + 1] * fDeltaMass;
	    }

	    for (k = 0; k < STEV_AGB_N_METALLICITY; k++) {
	       idxBuffer = stevRowMajorIndex(k, i, j, STEV_AGB_N_METALLICITY,
					     STEV_INTERP_N_MASS, STEV_N_ELEM);

	       buffer.afAGB_Yields[idxBuffer] = 0.0f;
	    }
	 }

	 for (k = 0; k < STEV_CCSN_N_METALLICITY; k++) {
	    idxTable  = stevRowMajorIndex(k, idxMass, 0, STEV_CCSN_N_METALLICITY,
					  STEV_CCSN_N_MASS, 1);
	    idxBuffer = stevRowMajorIndex(k, i, 0, STEV_CCSN_N_METALLICITY,
					  STEV_INTERP_N_MASS, 1);

	    /* WARNING: The following formulas are correct only when mass is the last index
	       in the tables */
	    buffer.afCCSN_MetalYield[idxBuffer] =
	       CCSNdata->pfMetalYield[idxTable] * (1.0 - fDeltaMass) +
	       CCSNdata->pfMetalYield[idxTable + 1] * fDeltaMass;

	    buffer.afCCSN_EjectedMass[idxBuffer] =
	       CCSNdata->pfEjectedMass[idxTable] * (1.0 - fDeltaMass) +
	       CCSNdata->pfEjectedMass[idxTable + 1] * fDeltaMass;
	 }

	 for (k = 0; k < STEV_AGB_N_METALLICITY; k++) {
	    idxBuffer = stevRowMajorIndex(k, i, 0, STEV_AGB_N_METALLICITY,
					  STEV_INTERP_N_MASS, 1);

	    buffer.afAGB_MetalYield[idxBuffer] = 0.0f;
	    buffer.afAGB_EjectedMass[idxBuffer] = 0.0f;
	 }
      }
   }


   for (i = 0; i < STEV_CCSN_N_METALLICITY; i++) {
      if (CCSNdata->pfZs[i] > 0.0f)
	 buffer.afCCSN_Zs[i] = log10(CCSNdata->pfZs[i]);
      else
	 buffer.afCCSN_Zs[i] = STEV_MIN_LOG_METALLICITY;
   }
   
   for (i = 0; i < STEV_AGB_N_METALLICITY; i++) {
      if (AGBdata->pfZs[i] > 0.0f)
	 buffer.afAGB_Zs[i] = log10(AGBdata->pfZs[i]);
      else
	 buffer.afAGB_Zs[i] = STEV_MIN_LOG_METALLICITY;
   }

   for (i = 0; i < STEV_N_ELEM; i++)
      buffer.afSNIa_EjectedMass[i] = SNIadata->pfEjectedMass[i];
   buffer.fSNIa_EjectedMetalMass = *SNIadata->pfMetalYield;

   for (i = 0; i < STEV_LIFETIMES_N_METALLICITY; i++) {
      if (LifetimesData->pfZs[i] > 0.0f)
	 buffer.afLifetimes_Zs[i] = log10(LifetimesData->pfZs[i]);
      else
	 buffer.afLifetimes_Zs[i] = STEV_MIN_LOG_METALLICITY;
   }
   for (i = 0; i < STEV_LIFETIMES_N_MASS; i++)
      buffer.afLifetimes_Masses[i] = log10(LifetimesData->pfMasses[i]);
   for (i = 0; i < STEV_LIFETIMES_N_METALLICITY * STEV_LIFETIMES_N_MASS; i++) {
      buffer.afLifetimes[i] = log10(LifetimesData->pfLifetimes[i] * SECPERYEAR /
				    msr->param.dSecUnit);
   }


   pstStellarEvolutionInit(msr->pst, (void *) &buffer, sizeof(buffer), NULL, 0);


   stevFreeTable(CCSNdata);
   stevFreeTable(AGBdata);
   stevFreeSNIaTable(SNIadata);
   stevFreeLifetimesTable(LifetimesData);
}


int pstStellarEvolutionInit(PST pst, void *vin, int nIn, void *vout, int nOut) {
   mdlassert(pst->mdl, nIn == sizeof(struct inStellarEvolution));
   if (pst->nLeaves > 1) {
      int rID = mdlReqService(pst->mdl, pst->idUpper, PST_STELLAREVOLUTIONINIT, vin, nIn);
      pstStellarEvolutionInit(pst->pstLower, vin, nIn, NULL, 0);
      mdlGetReply(pst->mdl, rID, NULL, NULL);
   }
   else {
      struct inStellarEvolution *in = vin;
      pkdStellarEvolutionInit(pst->plcl->pkd, in);
   }

   return 0;
}


int pkdStellarEvolutionInit(PKD pkd, STEV_DATA *data) {
   pkd->StelEvolData = (STEV_DATA *) malloc(sizeof(STEV_DATA));
   assert(pkd->StelEvolData != NULL);
   memcpy(pkd->StelEvolData, data, sizeof(STEV_DATA));

   if (strcmp(pkd->param.achSNIa_DTDtype, "exponential") == 0) {
      pkd->StelEvolData->fcnNumSNIa = stevExponentialNumSNIa;
   }
   else if (strcmp(pkd->param.achSNIa_DTDtype, "powerlaw") == 0) {
      pkd->StelEvolData->fcnNumSNIa = stevPowerlawNumSNIa;
   }

   if (pkd->param.nGrid > 0) {
      /* ICs have been generated */
      for (int i = 0; i < pkd->nLocal; ++i) {
	 PARTICLE *p = pkdParticle(pkd,i);
	 if (pkdIsGas(pkd, p)) {
	    SPHFIELDS *pSph = pkdSph(pkd,p);
	    float fMass = pkdMass(pkd, p);

            pSph->afElemMass[ELEMENT_H]  = pkd->param.dInitialH  * fMass;
            pSph->afElemMass[ELEMENT_He] = pkd->param.dInitialHe * fMass;
            pSph->afElemMass[ELEMENT_C]  = pkd->param.dInitialC  * fMass;
            pSph->afElemMass[ELEMENT_N]  = pkd->param.dInitialN  * fMass;
            pSph->afElemMass[ELEMENT_O]  = pkd->param.dInitialO  * fMass;
            pSph->afElemMass[ELEMENT_Ne] = pkd->param.dInitialNe * fMass;
            pSph->afElemMass[ELEMENT_Mg] = pkd->param.dInitialMg * fMass;
            pSph->afElemMass[ELEMENT_Si] = pkd->param.dInitialSi * fMass;
            pSph->afElemMass[ELEMENT_Fe] = pkd->param.dInitialFe * fMass;

	    pSph->fMetalMass = pkd->param.dInitialMetallicity * fMass;
	 }
      }
   }
   else if (pkd->param.achInFile[0]) {
      /* Input file has been read */
      for (int i = 0; i < pkd->nLocal; ++i) {
	 PARTICLE *p = pkdParticle(pkd,i);
	 if (pkdIsStar(pkd, p)) {
	    STARFIELDS *pStar = pkdStar(pkd, p);

	    for (int j = 0; j < ELEMENT_COUNT; j++)
	       pStar->afCumElemMassEj[j] = 0.0f;
	    pStar->fCumMetalMassEj = 0.0f;

	    int iIdxZ;
	    float fLogZ;
	    if (pStar->fMetalAbun > 0.0f)
	       fLogZ = log10(pStar->fMetalAbun);
	    else
	       fLogZ = STEV_MIN_LOG_METALLICITY;

	    stevGetIndex1D(pkd->StelEvolData->afCCSN_Zs, STEV_CCSN_N_METALLICITY,
			   fLogZ, &iIdxZ, &pStar->CCSN.fDeltaZ);
	    pStar->CCSN.oZ = iIdxZ * STEV_INTERP_N_MASS;
	    stevGetIndex1D(pkd->StelEvolData->afAGB_Zs, STEV_AGB_N_METALLICITY,
			   fLogZ, &iIdxZ, &pStar->AGB.fDeltaZ);
	    pStar->AGB.oZ = iIdxZ * STEV_INTERP_N_MASS;
	    stevGetIndex1D(pkd->StelEvolData->afLifetimes_Zs, STEV_LIFETIMES_N_METALLICITY,
			   fLogZ, &iIdxZ, &pStar->Lifetimes.fDeltaZ);
	    pStar->Lifetimes.oZ = iIdxZ * STEV_LIFETIMES_N_MASS;

	    pStar->fCCSNOnsetTime = stevLifetimeFunction(pkd, pStar, pkd->param.dIMF_MaxMass);
	    pStar->fSNIaOnsetTime = stevLifetimeFunction(pkd, pStar, pkd->param.fSNIa_MaxMass);

	    assert(pStar->fLastEnrichTime >= pStar->fCCSNOnsetTime);

	    pStar->fLastEnrichMass =
	       stevInverseLifetimeFunction(pkd, pStar, pStar->fLastEnrichTime);
	    pStar->iLastEnrichMassIdx =
	       stevGetIMFMassIndex(pkd->StelEvolData->afMasses, STEV_INTERP_N_MASS,
				   pStar->fLastEnrichMass, STEV_INTERP_N_MASS - 1) + 1;
	 }
      }
   }

   return 0;
}


void initChemEnrich(void *vpkd, void *vp) {
   PKD pkd = vpkd;
   PARTICLE *p = vp;

   if (pkdIsGas(pkd, p)) {
      *((float *) pkdField(p, pkd->oMass)) = 0.0f;

      SPHFIELDS *pSph = pkdSph(pkd,p);
      for (int j = 0; j < ELEMENT_COUNT; j++)
	 pSph->afElemMass[j]  = 0.0f;
      pSph->fMetalMass = 0.0f;
   }
}


void combChemEnrich(void *vpkd, void *vp1, void *vp2) {
   PKD pkd = vpkd;
   PARTICLE *p1 = vp1;
   PARTICLE *p2 = vp2;

   if (pkdIsGas(pkd, p1) && pkdIsGas(pkd, p2)) {
      *((float *) pkdField(p1, pkd->oMass)) += pkdMass(pkd, p2);

      SPHFIELDS *pSph1 = pkdSph(pkd, p1);
      SPHFIELDS *pSph2 = pkdSph(pkd, p2);
      for (int j = 0; j < ELEMENT_COUNT; j++)
	 pSph1->afElemMass[j] += pSph2->afElemMass[j];
      pSph1->fMetalMass += pSph2->fMetalMass;
   }
}

/* 
   One must be careful with the mass units when computing the ejected mass, such that
   the correct value is obtained in code mass units. Since the Stellar Evolution formalism
   allows one to compute the mass ejected per unit SSP mass, the only quantity that must be 
   in code mass units is the initial mass of the star particle. Units of ejecta and 
   initial mass from the tables, SNIa number and IMF normalization, and the code parameters
   to which they are compared can all be arbitrary, provided that they are all the same. 
   On the other hand, since time parameters/variables (e.g. stellar lifetimes) must be 
   compared and/or operated directly with the times in the simulation, they must all be in 
   code units. 
   Hence, here all the masses from the tables (pkd->StelEvolData) and the relevant parameters
   in pkd->param are in Msol. Conversely, pStar->[fInitialMass,afCumElemMassEj,fCumMetalMassEj]
   and pkdMass(pkd,p) are in code units. Furthermore, all times are in code units.
*/
void smChemEnrich(PARTICLE *p, float fBall, int nSmooth, NN *nnList, SMF *smf) {
   int i;
   PKD pkd = smf->pkd;

   assert(pkdIsStar(pkd, p));
   STARFIELDS *pStar = pkdStar(pkd, p);

   float ti, Mi;
   int idxMi;
   ti = pStar->fLastEnrichTime;
   Mi = pStar->fLastEnrichMass;
   idxMi = pStar->iLastEnrichMassIdx;

   const float tf = (float)smf->dTime - pStar->fTimer;
   const float Mf = stevInverseLifetimeFunction(pkd, pStar, tf);
   const int idxMf = stevGetIMFMassIndex(pkd->StelEvolData->afMasses, STEV_INTERP_N_MASS,
					 Mf, idxMi);

   pStar->fLastEnrichTime = tf;
   pStar->fLastEnrichMass = Mf;
   pStar->iLastEnrichMassIdx = idxMf + 1;


   /* Note: The parameter pStar->[AGB,CCSN,Lifetimes].oZ contains the index of the 
      interpolation's lower metallicity array multiplied by the number of mass bins.
      Since this multiplication must always be made, it is done once and for all in the
      the function pkdStarForm. */

   float afElemMass[STEV_N_ELEM];
   float fMetalMass;
   for (i = 0; i < STEV_N_ELEM; i++)
      afElemMass[i] = 0.0f;
   fMetalMass = 0.0f;

   const float Mtrans = pkd->param.fCCSN_MinMass;
   if (Mf > Mtrans) {
      /* CCSN only */
      stevComputeMassToEject(pkd->StelEvolData->afCCSN_Yields,
			     pkd->StelEvolData->afCCSN_MetalYield,
			     pkd->StelEvolData->afCCSN_EjectedMass,
			     pkd->StelEvolData->afMasses,
			     pkd->StelEvolData->afLogMasses,
			     pkd->StelEvolData->afIMF,
			     STEV_CCSN_N_METALLICITY, STEV_INTERP_N_MASS, STEV_N_ELEM,
			     pStar->afElemAbun, pStar->fMetalAbun, idxMf, idxMi, Mf, Mi,
			     pStar->CCSN.oZ, pStar->CCSN.fDeltaZ, afElemMass, &fMetalMass);
   }
   else if (Mi < Mtrans) {
      /* AGB only */
      stevComputeMassToEject(pkd->StelEvolData->afAGB_Yields,
			     pkd->StelEvolData->afAGB_MetalYield,
			     pkd->StelEvolData->afAGB_EjectedMass,
			     pkd->StelEvolData->afMasses,
			     pkd->StelEvolData->afLogMasses,
			     pkd->StelEvolData->afIMF,
			     STEV_AGB_N_METALLICITY, STEV_INTERP_N_MASS, STEV_N_ELEM,
			     pStar->afElemAbun, pStar->fMetalAbun, idxMf, idxMi, Mf, Mi,
			     pStar->AGB.oZ, pStar->AGB.fDeltaZ, afElemMass, &fMetalMass);
   }
   else {
      /* Mixed CCSN and AGB */
      int idxMtrans = stevGetIMFMassIndex(pkd->StelEvolData->afMasses, STEV_INTERP_N_MASS,
					  Mtrans, idxMi);

      stevComputeMassToEject(pkd->StelEvolData->afCCSN_Yields,
			     pkd->StelEvolData->afCCSN_MetalYield,
			     pkd->StelEvolData->afCCSN_EjectedMass,
			     pkd->StelEvolData->afMasses,
			     pkd->StelEvolData->afLogMasses,
			     pkd->StelEvolData->afIMF,
			     STEV_CCSN_N_METALLICITY, STEV_INTERP_N_MASS, STEV_N_ELEM,
			     pStar->afElemAbun, pStar->fMetalAbun, idxMtrans, idxMi, Mtrans, Mi,
			     pStar->CCSN.oZ, pStar->CCSN.fDeltaZ, afElemMass, &fMetalMass);

      idxMtrans++;
      stevComputeMassToEject(pkd->StelEvolData->afAGB_Yields,
			     pkd->StelEvolData->afAGB_MetalYield,
			     pkd->StelEvolData->afAGB_EjectedMass,
			     pkd->StelEvolData->afMasses,
			     pkd->StelEvolData->afLogMasses,
			     pkd->StelEvolData->afIMF,
			     STEV_AGB_N_METALLICITY, STEV_INTERP_N_MASS, STEV_N_ELEM,
			     pStar->afElemAbun, pStar->fMetalAbun, idxMf, idxMtrans, Mf, Mtrans,
			     pStar->AGB.oZ, pStar->AGB.fDeltaZ, afElemMass, &fMetalMass);
   }

   const float fNumSNIa = (*pkd->StelEvolData->fcnNumSNIa)(pkd, pStar, ti, tf);
   for (i = 0; i < STEV_N_ELEM; i++) {
      afElemMass[i] += fNumSNIa * pkd->StelEvolData->afSNIa_EjectedMass[i];
      afElemMass[i] *= pStar->fInitialMass;
      pStar->afCumElemMassEj[i] += afElemMass[i];
   }

   fMetalMass += fNumSNIa * pkd->StelEvolData->fSNIa_EjectedMetalMass;
   fMetalMass *= pStar->fInitialMass;
   pStar->fCumMetalMassEj += fMetalMass;


   const float fTotalMass = afElemMass[ELEMENT_H] + afElemMass[ELEMENT_He] + fMetalMass;
   *((float *) pkdField(p, pkd->oMass)) -= fTotalMass;
   assert(pkdMass(pkd, p) > 0.0f);


   PARTICLE *q;
   SPHFIELDS *qSph;
   float fWeightSum = 0.0f;
   for (i = 0; i < nSmooth; i++) {
      q = nnList[i].pPart;

      /* Three alternatives to ignore the star particle */
      if (!pkdIsGas(pkd, q)) continue;
      if (q == p) continue;
      if (nnList[i].dx == 0.0f && nnList[i].dy == 0.0f && nnList[i].dz == 0.0f) continue;

      const double dRpq = sqrt(nnList[i].fDist2);
      const float fWeight = (float)cubicSplineKernel(dRpq, fBall) / pkdDensity(pkd, q);

      qSph = pkdSph(pkd, q);
      for (int j = 0; j < STEV_N_ELEM; j++)
	 qSph->afElemMass[j] += fWeight * afElemMass[j];
      qSph->fMetalMass += fWeight * fMetalMass;

      *((float *) pkdField(q, pkd->oMass)) += fWeight * fTotalMass;

      fWeightSum += fWeight;
   }

   for (i = 0; i < nSmooth; i++) {
      q = nnList[i].pPart;

      /* Three alternatives to ignore the star particle */
      if (!pkdIsGas(pkd, q)) continue;
      if (q == p) continue;
      if (nnList[i].dx == 0.0f && nnList[i].dy == 0.0f && nnList[i].dz == 0.0f) continue;

      qSph = pkdSph(pkd, q);
      for (int j = 0; j < STEV_N_ELEM; j++)
	 qSph->afElemMass[j] /= fWeightSum;
      qSph->fMetalMass /= fWeightSum;

      *((float *) pkdField(q, pkd->oMass)) /= fWeightSum;
   }
}


STEV_RAWDATA *stevReadTable(char *pszPath, char **pszFields, int nFields) {
   STEV_RAWDATA *RawData = malloc(sizeof(STEV_RAWDATA));
   assert(RawData != NULL);


   hid_t file_id = H5Fopen(pszPath, H5F_ACC_RDONLY, H5P_DEFAULT);
   assert(file_id >= 0);


   hid_t dataset = H5Dopen(file_id, "/Number_of_metallicities", H5P_DEFAULT);
   hid_t datatype = H5Dget_type(dataset);
   assert(H5Tequal(datatype, H5T_NATIVE_INT));
   herr_t status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &RawData->nZs);
   assert(status >= 0);
   status = H5Tclose(datatype);
   assert(status >= 0);
   status = H5Dclose(dataset);
   assert(status >= 0);

   dataset = H5Dopen(file_id, "/Number_of_species", H5P_DEFAULT);
   datatype = H5Dget_type(dataset);
   assert(H5Tequal(datatype, H5T_NATIVE_INT));
   status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &RawData->nSpecs);
   assert(status >= 0);
   status = H5Tclose(datatype);
   assert(status >= 0);
   status = H5Dclose(dataset);
   assert(status >= 0);

   dataset = H5Dopen(file_id, "/Number_of_masses", H5P_DEFAULT);
   datatype = H5Dget_type(dataset);
   assert(H5Tequal(datatype, H5T_NATIVE_INT));
   status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &RawData->nMasses);
   assert(status >= 0);
   status = H5Tclose(datatype);
   assert(status >= 0);
   status = H5Dclose(dataset);
   assert(status >= 0);


   assert(RawData->nZs == nFields);


   RawData->pfZs = (float *) malloc(RawData->nZs * sizeof(float));
   assert(RawData->pfZs != NULL);
   RawData->pfMasses = (float *) malloc(RawData->nMasses * sizeof(float));
   assert(RawData->pfMasses != NULL);
   RawData->pfYields = (float *) malloc(RawData->nZs * RawData->nSpecs * RawData->nMasses *
					sizeof(float));
   assert(RawData->pfYields != NULL);
   RawData->pfMetalYield = (float *) malloc(RawData->nZs * RawData->nMasses * sizeof(float));
   assert(RawData->pfMetalYield != NULL);
   RawData->pfEjectedMass = (float *) malloc(RawData->nZs * RawData->nMasses * sizeof(float));
   assert(RawData->pfEjectedMass != NULL);


   dataset = H5Dopen(file_id, "/Metallicities", H5P_DEFAULT);
   datatype = H5Dget_type(dataset);
   assert(H5Tequal(datatype, H5T_NATIVE_FLOAT));
   status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->pfZs);
   assert(status >= 0);
   status = H5Tclose(datatype);
   assert(status >= 0);
   status = H5Dclose(dataset);
   assert(status >= 0);

   dataset = H5Dopen(file_id, "/Masses", H5P_DEFAULT);
   datatype = H5Dget_type(dataset);
   assert(H5Tequal(datatype, H5T_NATIVE_FLOAT));
   status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->pfMasses);
   assert(status >= 0);
   status = H5Tclose(datatype);
   assert(status >= 0);
   status = H5Dclose(dataset);
   assert(status >= 0);


   char achField[64];
   for (int i = 0; i < RawData->nZs; i++) {
      int n = i * RawData->nSpecs * RawData->nMasses;

      sprintf(achField, "/Yields/%s/Yield", pszFields[i]);
      dataset = H5Dopen(file_id, achField, H5P_DEFAULT);
      datatype = H5Dget_type(dataset);
      assert(H5Tequal(datatype, H5T_NATIVE_FLOAT));
      status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->pfYields + n);
      assert(status >= 0);
      status = H5Tclose(datatype);
      assert(status >= 0);
      status = H5Dclose(dataset);
      assert(status >= 0);


      n = i * RawData->nMasses;

      sprintf(achField, "/Yields/%s/Total_Metals", pszFields[i]);
      dataset = H5Dopen(file_id, achField, H5P_DEFAULT);
      datatype = H5Dget_type(dataset);
      assert(H5Tequal(datatype, H5T_NATIVE_FLOAT));
      status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->pfMetalYield + n);
      assert(status >= 0);
      status = H5Tclose(datatype);
      assert(status >= 0);
      status = H5Dclose(dataset);
      assert(status >= 0);

      sprintf(achField, "/Yields/%s/Ejected_mass", pszFields[i]);
      dataset = H5Dopen(file_id, achField, H5P_DEFAULT);
      datatype = H5Dget_type(dataset);
      assert(H5Tequal(datatype, H5T_NATIVE_FLOAT));
      status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->pfEjectedMass + n);
      assert(status >= 0);
      status = H5Tclose(datatype);
      assert(status >= 0);
      status = H5Dclose(dataset);
      assert(status >= 0);
   }

   status = H5Fclose(file_id);
   assert(status >= 0);

   return RawData;
}


/* This function assumes that the table of mass ejected from Type Ia SN is independent of
   the progenitor's initial mass and metallicity */
STEV_RAWDATA *stevReadSNIaTable(char *pszPath) {
   STEV_RAWDATA *RawData = malloc(sizeof(STEV_RAWDATA));
   assert(RawData != NULL);


   hid_t file_id = H5Fopen(pszPath, H5F_ACC_RDONLY, H5P_DEFAULT);
   assert(file_id >= 0);


   hid_t dataset = H5Dopen(file_id, "/Number_of_species", H5P_DEFAULT);
   hid_t datatype = H5Dget_type(dataset);
   assert(H5Tequal(datatype, H5T_NATIVE_INT));
   herr_t status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &RawData->nSpecs);
   assert(status >= 0);
   status = H5Tclose(datatype);
   assert(status >= 0);
   status = H5Dclose(dataset);
   assert(status >= 0);


   RawData->pfMetalYield = (float *) malloc(sizeof(float));
   assert(RawData->pfMetalYield != NULL);
   RawData->pfEjectedMass = (float *) malloc(RawData->nSpecs * sizeof(float));
   assert(RawData->pfEjectedMass != NULL);


   dataset = H5Dopen(file_id, "/Yield", H5P_DEFAULT);
   datatype = H5Dget_type(dataset);
   assert(H5Tequal(datatype, H5T_NATIVE_FLOAT));
   status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->pfEjectedMass);
   assert(status >= 0);
   status = H5Tclose(datatype);
   assert(status >= 0);
   status = H5Dclose(dataset);
   assert(status >= 0);

   dataset = H5Dopen(file_id, "/Total_Metals", H5P_DEFAULT);
   datatype = H5Dget_type(dataset);
   assert(H5Tequal(datatype, H5T_NATIVE_FLOAT));
   status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->pfMetalYield);
   assert(status >= 0);
   status = H5Tclose(datatype);
   assert(status >= 0);
   status = H5Dclose(dataset);
   assert(status >= 0);


   status = H5Fclose(file_id);
   assert(status >= 0);

   return RawData;
}


STEV_RAWDATA *stevReadLifetimesTable(char *pszPath) {
   STEV_RAWDATA *RawData = malloc(sizeof(STEV_RAWDATA));
   assert(RawData != NULL);


   hid_t file_id = H5Fopen(pszPath, H5F_ACC_RDONLY, H5P_DEFAULT);
   assert(file_id >= 0);


   hid_t dataset = H5Dopen(file_id, "/Number_of_metallicities", H5P_DEFAULT);
   hid_t datatype = H5Dget_type(dataset);
   assert(H5Tequal(datatype, H5T_NATIVE_INT));
   herr_t status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &RawData->nZs);
   assert(status >= 0);
   status = H5Tclose(datatype);
   assert(status >= 0);
   status = H5Dclose(dataset);
   assert(status >= 0);

   dataset = H5Dopen(file_id, "/Number_of_masses", H5P_DEFAULT);
   datatype = H5Dget_type(dataset);
   assert(H5Tequal(datatype, H5T_NATIVE_INT));
   status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &RawData->nMasses);
   assert(status >= 0);
   status = H5Tclose(datatype);
   assert(status >= 0);
   status = H5Dclose(dataset);
   assert(status >= 0);


   RawData->pfZs = (float *) malloc(RawData->nZs * sizeof(float));
   assert(RawData->pfZs != NULL);
   RawData->pfMasses = (float *) malloc(RawData->nMasses * sizeof(float));
   assert(RawData->pfMasses != NULL);
   RawData->pfLifetimes = (float *) malloc(RawData->nZs * RawData->nMasses * sizeof(float));
   assert(RawData->pfLifetimes != NULL);


   dataset = H5Dopen(file_id, "/Metallicities", H5P_DEFAULT);
   datatype = H5Dget_type(dataset);
   assert(H5Tequal(datatype, H5T_NATIVE_FLOAT));
   status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->pfZs);
   assert(status >= 0);
   status = H5Tclose(datatype);
   assert(status >= 0);
   status = H5Dclose(dataset);
   assert(status >= 0);

   dataset = H5Dopen(file_id, "/Masses", H5P_DEFAULT);
   datatype = H5Dget_type(dataset);
   assert(H5Tequal(datatype, H5T_NATIVE_FLOAT));
   status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->pfMasses);
   assert(status >= 0);
   status = H5Tclose(datatype);
   assert(status >= 0);
   status = H5Dclose(dataset);
   assert(status >= 0);

   dataset = H5Dopen(file_id, "/Lifetimes", H5P_DEFAULT);
   datatype = H5Dget_type(dataset);
   assert(H5Tequal(datatype, H5T_NATIVE_FLOAT));
   status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->pfLifetimes);
   assert(status >= 0);
   status = H5Tclose(datatype);
   assert(status >= 0);
   status = H5Dclose(dataset);
   assert(status >= 0);


   status = H5Fclose(file_id);
   assert(status >= 0);

   return RawData;
}


void stevFreeTable(STEV_RAWDATA *RawData) {
   free(RawData->pfZs);
   free(RawData->pfMasses);
   free(RawData->pfYields);
   free(RawData->pfMetalYield);
   free(RawData->pfEjectedMass);
   free(RawData);
}


void stevFreeSNIaTable(STEV_RAWDATA *RawData) {
   free(RawData->pfMetalYield);
   free(RawData->pfEjectedMass);
   free(RawData);
}


void stevFreeLifetimesTable(STEV_RAWDATA *RawData) {
   free(RawData->pfZs);
   free(RawData->pfMasses);
   free(RawData->pfLifetimes);
   free(RawData);
}


float stevExponentialNumSNIa(PKD pkd, STARFIELDS *pStar, float fInitialTime, float fFinalTime) {
   if (fFinalTime <= pStar->fSNIaOnsetTime) {
      return 0.0f;
   }
   else if (fInitialTime < pStar->fSNIaOnsetTime) {
      fInitialTime = pStar->fSNIaOnsetTime;
   }

   return pkd->param.dSNIa_Norm * (exp(-(fInitialTime / pkd->param.dSNIa_Scale)) -
				   exp(-(fFinalTime / pkd->param.dSNIa_Scale)));
}


float stevPowerlawNumSNIa(PKD pkd, STARFIELDS *pStar, float fInitialTime, float fFinalTime) {
   if (fFinalTime <= pStar->fSNIaOnsetTime) {
      return 0.0f;
   }
   else if (fInitialTime < pStar->fSNIaOnsetTime) {
      fInitialTime = pStar->fSNIaOnsetTime;
   }

   return pkd->param.dSNIa_Norm * (pow(fFinalTime, pkd->param.dSNIa_Scale + 1.0) -
				   pow(fInitialTime, pkd->param.dSNIa_Scale + 1.0));
}


#endif	/* STELLAR_EVOLUTION */
