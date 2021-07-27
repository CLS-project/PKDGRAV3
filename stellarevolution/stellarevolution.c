#ifdef STELLAR_EVOLUTION


#include <math.h>
#include <hdf5.h>

#include "stellarevolution.h"
#include "hydro/hydro.h"


void msrStellarEvolutionInit(MSR msr, double dTime) {
   assert(STEV_N_ELEM == ELEMENT_COUNT);

   int i;
   char achPath[280];
   STEV_RAWDATA *CCSNdata, *AGBdata, *SNIadata, *LifetimesData;

   /* Read the tables */
   sprintf(achPath, "%s/CCSN.hdf5", msr->param.achStEvolPath);
   CCSNdata = stevReadTable(achPath);
   assert(CCSNdata->nZs == STEV_CCSN_N_METALLICITY);
   assert(CCSNdata->nMasses == STEV_CCSN_N_MASS);
   assert(CCSNdata->nSpecs == STEV_N_ELEM);

   sprintf(achPath, "%s/AGB.hdf5", msr->param.achStEvolPath);
   AGBdata = stevReadTable(achPath);
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

   /* Allocate the chunk that will be sent to all workers. We append a double
      that will store the initial time of the simulation (dTime) */
   STEV_DATA *buffer = malloc(sizeof(STEV_DATA) + sizeof(double));

   /* NOTE: The lowest value of the initial mass array is set to the corresponding
      value of the Lifetimes table that is being used. The IMF is still normalized
      in the range [msr->param.dIMF_MinMass,msr->param.dIMF_MaxMass] */
   const double dMinMass = LifetimesData->pfMasses[0];
   double adMasses[STEV_INTERP_N_MASS], adIMF[STEV_INTERP_N_MASS];

   /* Setup the initial mass array */
   const double dMaxMass = msr->param.dIMF_MaxMass;
   const double dDeltaLog = (log10(dMaxMass) - log10(dMinMass)) / (STEV_INTERP_N_MASS - 1);
   const double dDelta = pow(10.0, dDeltaLog);
   double dMass = dMinMass;
   for (i = 0; i < STEV_INTERP_N_MASS; i++) {
      adMasses[i] = dMass;
      dMass *= dDelta;
   }

   /* Compute the IMF at the values of the initial mass array */
   if (strcmp(msr->param.achIMFtype, "chabrier") == 0) {
      stevChabrierIMF(adMasses, STEV_INTERP_N_MASS, msr->param.dIMF_MinMass,
		      msr->param.dIMF_MaxMass, adIMF);
   }
   else {
      printf("ERROR: Undefined IMF type has been given in achIMFtype parameter: %s\n",
	     msr->param.achIMFtype);
      assert(0);
   }

   /* Store data */
   for (i = 0; i < STEV_INTERP_N_MASS; i++) {
      buffer->afMasses[i] = adMasses[i];
      buffer->afIMFLogWeights[i] = adIMF[i] * adMasses[i];
   }
   buffer->fDeltaLogMass = dDeltaLog;

   /* Verify IMF normalization using the initial mass array, if possible */
   if (dMinMass == msr->param.dIMF_MinMass) {
      double IMFnorm = 0.0;
      for (i = 1; i < STEV_INTERP_N_MASS - 1; i++)
	 IMFnorm += adMasses[i] * buffer->afIMFLogWeights[i];
      IMFnorm += 0.5 * (adMasses[0] * buffer->afIMFLogWeights[0] +
			adMasses[i] * buffer->afIMFLogWeights[i]);
      IMFnorm *= M_LN10 * dDeltaLog;
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

   /* Convert CCSN/AGB initial mass arrays to log for interpolation */
   for (i = 0; i < CCSNdata->nMasses; i++)
      CCSNdata->pfMasses[i] = log10(CCSNdata->pfMasses[i]);
   for (i = 0; i < AGBdata->nMasses; i++)
      AGBdata->pfMasses[i] = log10(AGBdata->pfMasses[i]);

   /* Interpolate yields and ejected masses to the initial mass array */
   stevInterpToIMFSampling(buffer, CCSNdata, AGBdata, msr->param.dCCSN_MinMass);

   /* Store remaining data */
   for (i = 0; i < STEV_CCSN_N_METALLICITY; i++) {
      if (CCSNdata->pfZs[i] > 0.0f)
	 buffer->afCCSN_Zs[i] = log10(CCSNdata->pfZs[i]);
      else
	 buffer->afCCSN_Zs[i] = STEV_MIN_LOG_METALLICITY;
   }
   
   for (i = 0; i < STEV_AGB_N_METALLICITY; i++) {
      if (AGBdata->pfZs[i] > 0.0f)
	 buffer->afAGB_Zs[i] = log10(AGBdata->pfZs[i]);
      else
	 buffer->afAGB_Zs[i] = STEV_MIN_LOG_METALLICITY;
   }

   for (i = 0; i < STEV_N_ELEM; i++)
      buffer->afSNIa_EjectedMass[i] = SNIadata->pfEjectedMass[i];
   buffer->fSNIa_EjectedMetalMass = *SNIadata->pfMetalYield;

   for (i = 0; i < STEV_LIFETIMES_N_METALLICITY; i++) {
      if (LifetimesData->pfZs[i] > 0.0f)
	 buffer->afLifetimes_Zs[i] = log10(LifetimesData->pfZs[i]);
      else
	 buffer->afLifetimes_Zs[i] = STEV_MIN_LOG_METALLICITY;
   }
   for (i = 0; i < STEV_LIFETIMES_N_MASS; i++)
      buffer->afLifetimes_Masses[i] = log10(LifetimesData->pfMasses[i]);
   for (i = 0; i < STEV_LIFETIMES_N_METALLICITY * STEV_LIFETIMES_N_MASS; i++) {
      buffer->afLifetimes[i] = log10(LifetimesData->pfLifetimes[i] * SECPERYEAR /
				     msr->param.dSecUnit);
   }

   *((double *) (buffer + 1)) = dTime;

   /* Send data and finish */
   pstStellarEvolutionInit(msr->pst, buffer, sizeof(STEV_DATA) + sizeof(double), NULL, 0);

   free(buffer);
   stevFreeTable(CCSNdata);
   stevFreeTable(AGBdata);
   stevFreeSNIaTable(SNIadata);
   stevFreeLifetimesTable(LifetimesData);
}


int pstStellarEvolutionInit(PST pst, void *vin, int nIn, void *vout, int nOut) {
   mdlassert(pst->mdl, nIn == sizeof(struct inStellarEvolution) + sizeof(double));
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

   double dTime = *((double *) (data + 1));

   if (strcmp(pkd->param.achSNIa_DTDtype, "exponential") == 0) {
      pkd->StelEvolData->fcnNumSNIa = stevExponentialNumSNIa;
   }
   else if (strcmp(pkd->param.achSNIa_DTDtype, "powerlaw") == 0) {
      pkd->StelEvolData->fcnNumSNIa = stevPowerlawNumSNIa;
   }
   else {
      printf("ERROR: Undefined SNIa DTD type has been given in achSNIa_DTDtype\n"
	     "       parameter: %s\n", pkd->param.achSNIa_DTDtype);
      assert(0);
   }

   for (int i = 0; i < pkd->nLocal; ++i) {
      PARTICLE *p = pkdParticle(pkd, i);

      if (pkdIsGas(pkd, p)) {
	 SPHFIELDS *pSph = pkdSph(pkd,p);
	 float fMass = pkdMass(pkd, p);

	 if (pSph->afElemMass[ELEMENT_H] < 0.0f) {
            pSph->afElemMass[ELEMENT_H]  = pkd->param.dInitialH  * fMass;
            pSph->afElemMass[ELEMENT_He] = pkd->param.dInitialHe * fMass;
            pSph->afElemMass[ELEMENT_C]  = pkd->param.dInitialC  * fMass;
            pSph->afElemMass[ELEMENT_N]  = pkd->param.dInitialN  * fMass;
            pSph->afElemMass[ELEMENT_O]  = pkd->param.dInitialO  * fMass;
            pSph->afElemMass[ELEMENT_Ne] = pkd->param.dInitialNe * fMass;
            pSph->afElemMass[ELEMENT_Mg] = pkd->param.dInitialMg * fMass;
            pSph->afElemMass[ELEMENT_Si] = pkd->param.dInitialSi * fMass;
            pSph->afElemMass[ELEMENT_Fe] = pkd->param.dInitialFe * fMass;
	 }

	 if (pSph->fMetalMass < 0.0f)
	    pSph->fMetalMass = pkd->param.dInitialMetallicity * fMass;
      }
      else if (pkdIsStar(pkd, p)) {
	 STARFIELDS *pStar = pkdStar(pkd, p);

	 if (pStar->afElemAbun[ELEMENT_H] < 0.0f) {
	    pStar->afElemAbun[ELEMENT_H]  = pkd->param.dInitialH;
	    pStar->afElemAbun[ELEMENT_He] = pkd->param.dInitialHe;
	    pStar->afElemAbun[ELEMENT_C]  = pkd->param.dInitialC;
	    pStar->afElemAbun[ELEMENT_N]  = pkd->param.dInitialN;
	    pStar->afElemAbun[ELEMENT_O]  = pkd->param.dInitialO;
	    pStar->afElemAbun[ELEMENT_Ne] = pkd->param.dInitialNe;
	    pStar->afElemAbun[ELEMENT_Mg] = pkd->param.dInitialMg;
	    pStar->afElemAbun[ELEMENT_Si] = pkd->param.dInitialSi;
	    pStar->afElemAbun[ELEMENT_Fe] = pkd->param.dInitialFe;
	 }

	 if (pStar->fMetalAbun < 0.0f)
	    pStar->fMetalAbun = pkd->param.dInitialMetallicity;

	 if (pStar->fInitialMass < 0.0f)
	    pStar->fInitialMass = pkdMass(pkd, p);

	 if (pkd->param.bChemEnrich == 0) {
	    pStar->fNextEnrichTime = INFINITY;
	    continue;
	 }

	 if (pStar->fTimer <= 0.0f) {
	    pStar->fTimer = dTime;
	    pStar->fLastEnrichTime = 0.0f;
	 }
	 else {
	    pStar->fLastEnrichTime = (float)dTime - pStar->fTimer;
	 }

	 int iIdxZ;
	 float fLogZ;
	 if (pStar->fMetalAbun > 0.0f)
	    fLogZ = log10f(pStar->fMetalAbun);
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
	 pStar->fSNIaOnsetTime = stevLifetimeFunction(pkd, pStar, pkd->param.dSNIa_MaxMass);

	 if (pStar->fLastEnrichTime < pStar->fCCSNOnsetTime) {
	    pStar->fLastEnrichTime = pStar->fCCSNOnsetTime;
	    pStar->fLastEnrichMass = pkd->param.dIMF_MaxMass;
	    pStar->iLastEnrichMassIdx = STEV_INTERP_N_MASS - 1;
	    pStar->fNextEnrichTime = pStar->fCCSNOnsetTime + pStar->fTimer;
	 }
	 else {
	    pStar->fLastEnrichMass =
	       stevInverseLifetimeFunction(pkd, pStar, pStar->fLastEnrichTime);
	    pStar->iLastEnrichMassIdx =
	       stevGetIMFMassIndex(pkd->StelEvolData->afMasses, STEV_INTERP_N_MASS,
				   pStar->fLastEnrichMass, STEV_INTERP_N_MASS - 1) + 1;
	    pStar->fNextEnrichTime = dTime;
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

      pSph->mom[0] = pSph->mom[1] = pSph->mom[2] = 0.0;
      pSph->E = 0.0;

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
      const float fOldMass = pkdMass(pkd, p1);
      *((float *) pkdField(p1, pkd->oMass)) += pkdMass(pkd, p2);
      const float fNewMass = pkdMass(pkd, p1);

      SPHFIELDS *pSph1 = pkdSph(pkd, p1);
      SPHFIELDS *pSph2 = pkdSph(pkd, p2);

      const float fOldEkin = (pSph1->mom[0] * pSph1->mom[0] + pSph1->mom[1] * pSph1->mom[1] +
			      pSph1->mom[2] * pSph1->mom[2]) / (2.0f * fOldMass);

      pSph1->mom[0] += pSph2->mom[0];
      pSph1->mom[1] += pSph2->mom[1];
      pSph1->mom[2] += pSph2->mom[2];

      const float fNewEkin = (pSph1->mom[0] * pSph1->mom[0] + pSph1->mom[1] * pSph1->mom[1] +
			      pSph1->mom[2] * pSph1->mom[2]) / (2.0f * fNewMass);

      pSph1->E += pSph2->E;
      const float fDeltaUint = pSph2->E - (fNewEkin - fOldEkin);
      pSph1->Uint += fDeltaUint;
#ifdef ENTROPY_SWITCH
      pSph1->S += fDeltaUint * (pkd->param.dConstGamma - 1.0) *
	          pow(pkdDensity(pkd, p1), 1.0 - pkd->param.dConstGamma);
#endif

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
   in pkd->param are in Msol. Conversely, pStar->fInitialMass and pkdMass(pkd,p) are in code
   units. Furthermore, all times are in code units.
*/
void smChemEnrich(PARTICLE *p, float fBall, int nSmooth, NN *nnList, SMF *smf) {
   int i;
   PKD pkd = smf->pkd;
   STARFIELDS *pStar = pkdStar(pkd, p);

   const float ti = pStar->fLastEnrichTime;
   const float Mi = pStar->fLastEnrichMass;
   const int idxMi = pStar->iLastEnrichMassIdx;

   const float tf = (float)smf->dTime - pStar->fTimer;
   const float Mf = stevInverseLifetimeFunction(pkd, pStar, tf);
   const int idxMf = stevGetIMFMassIndex(pkd->StelEvolData->afMasses, STEV_INTERP_N_MASS,
					 Mf, idxMi);

   pStar->fLastEnrichTime = tf;
   pStar->fLastEnrichMass = Mf;
   pStar->iLastEnrichMassIdx = idxMf + 1;


   /* Note: The parameter pStar->[AGB,CCSN,Lifetimes].oZ contains the index of the 
      interpolation's lower metallicity array multiplied by the number of mass bins.
      Since this multiplication must always be made, it is done once and for all in
      pkdStarForm or pkdStellarEvolutionInit, depending on whether we are dealing with
      a star particle born during the run or one that was already in the ICs. */

   float afElemMass[STEV_N_ELEM];
   float fMetalMass;
   for (i = 0; i < STEV_N_ELEM; i++)
      afElemMass[i] = 0.0f;
   fMetalMass = 0.0f;

   const float Mtrans = pkd->param.dCCSN_MinMass;
   if (Mf > Mtrans) {
      /* CCSN only */
      stevComputeMassToEject(pkd->StelEvolData->afCCSN_Yields,
			     pkd->StelEvolData->afCCSN_MetalYield,
			     pkd->StelEvolData->afCCSN_EjectedMass,
			     pkd->StelEvolData->afMasses,
			     pkd->StelEvolData->afIMFLogWeights,
			     pkd->StelEvolData->fDeltaLogMass,
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
			     pkd->StelEvolData->afIMFLogWeights,
			     pkd->StelEvolData->fDeltaLogMass,
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
			     pkd->StelEvolData->afIMFLogWeights,
			     pkd->StelEvolData->fDeltaLogMass,
			     STEV_CCSN_N_METALLICITY, STEV_INTERP_N_MASS, STEV_N_ELEM,
			     pStar->afElemAbun, pStar->fMetalAbun, idxMtrans, idxMi, Mtrans, Mi,
			     pStar->CCSN.oZ, pStar->CCSN.fDeltaZ, afElemMass, &fMetalMass);

      idxMtrans++;
      stevComputeMassToEject(pkd->StelEvolData->afAGB_Yields,
			     pkd->StelEvolData->afAGB_MetalYield,
			     pkd->StelEvolData->afAGB_EjectedMass,
			     pkd->StelEvolData->afMasses,
			     pkd->StelEvolData->afIMFLogWeights,
			     pkd->StelEvolData->fDeltaLogMass,
			     STEV_AGB_N_METALLICITY, STEV_INTERP_N_MASS, STEV_N_ELEM,
			     pStar->afElemAbun, pStar->fMetalAbun, idxMf, idxMtrans, Mf, Mtrans,
			     pStar->AGB.oZ, pStar->AGB.fDeltaZ, afElemMass, &fMetalMass);
   }

   const float fNumSNIa = (*pkd->StelEvolData->fcnNumSNIa)(pkd, pStar, ti, tf);
   for (i = 0; i < STEV_N_ELEM; i++) {
      afElemMass[i] += fNumSNIa * pkd->StelEvolData->afSNIa_EjectedMass[i];
      afElemMass[i] *= pStar->fInitialMass;
   }
   fMetalMass += fNumSNIa * pkd->StelEvolData->fSNIa_EjectedMetalMass;
   fMetalMass *= pStar->fInitialMass;


   const float fTotalMass = afElemMass[ELEMENT_H] + afElemMass[ELEMENT_He] + fMetalMass;
   pStar->fNextEnrichTime = stevComputeNextEnrichTime(smf->dTime, pStar->fInitialMass,
						      fTotalMass, tf - ti);
   *((float *) pkdField(p, pkd->oMass)) -= fTotalMass;
   assert(pkdMass(pkd, p) > 0.0f);


   const float fScaleFactorInv = 1.0 / csmTime2Exp(pkd->csm, smf->dTime);
   const float fScaleFactorInvSq = fScaleFactorInv * fScaleFactorInv;
   vel_t *pStarVel = pkdVel(pkd, p);

   const float fStarDeltaEkin = 0.5f * fTotalMass * (pStarVel[0] * pStarVel[0] +
				pStarVel[1] * pStarVel[1] + pStarVel[2] * pStarVel[2]);
   const float fWindEkin = (float)pkd->param.dWindSpecificEkin * fTotalMass;
   const float fSNIaEjEnergy = (fNumSNIa * (float)pkd->param.dMsolUnit) * pStar->fInitialMass *
                               (float)pkd->param.dSNIaEnergy;

   const float fStarEjEnergy = (fStarDeltaEkin + fWindEkin) * fScaleFactorInvSq +
                               fSNIaEjEnergy;


   PARTICLE *q;
   float fWeights[nSmooth];
   float fNormFactor = 0.0f;
   for (i = 0; i < nSmooth; i++) {
      q = nnList[i].pPart;
      if (q == p) continue;

      const double dRpq = sqrt(nnList[i].fDist2);
      fWeights[i] = cubicSplineKernel(dRpq, fBall) / pkdDensity(pkd, q);
      fNormFactor += fWeights[i];
   }
   fNormFactor = 1.0f / fNormFactor;

   for (i = 0; i < nSmooth; i++) {
      q = nnList[i].pPart;
      if (q == p) continue;

      fWeights[i] *= fNormFactor;

      const float fDeltaMass = fWeights[i] * fTotalMass;
      const float fOldMass = pkdMass(pkd, q);
      *((float *) pkdField(q, pkd->oMass)) += fDeltaMass;

      SPHFIELDS *qSph = pkdSph(pkd, q);

      float fOldEkin;           /* Gas particles in cache have their mass set to zero */
      if (fOldMass > 0.0f) {	/* in initChemEnrich */
	 fOldEkin = (qSph->mom[0] * qSph->mom[0] + qSph->mom[1] * qSph->mom[1] +
		     qSph->mom[2] * qSph->mom[2]) / (2.0f * fOldMass);
      }

      qSph->mom[0] += fDeltaMass * pStarVel[0] * fScaleFactorInv;
      qSph->mom[1] += fDeltaMass * pStarVel[1] * fScaleFactorInv;
      qSph->mom[2] += fDeltaMass * pStarVel[2] * fScaleFactorInv;

      if (fOldMass > 0.0f) {
	 const float fNewMass = pkdMass(pkd, q);
	 const float fNewEkin = (qSph->mom[0] * qSph->mom[0] + qSph->mom[1] * qSph->mom[1] +
				 qSph->mom[2] * qSph->mom[2]) / (2.0f * fNewMass);
	 const float fDeltaUint = fWeights[i] * fStarEjEnergy - (fNewEkin - fOldEkin);
	 qSph->Uint += fDeltaUint;
#ifdef ENTROPY_SWITCH
	 qSph->S += fDeltaUint * (pkd->param.dConstGamma - 1.0) *
	            pow(pkdDensity(pkd, q), 1.0 - pkd->param.dConstGamma);
#endif
      }

      qSph->E += fWeights[i] * fStarEjEnergy;

      for (int j = 0; j < STEV_N_ELEM; j++)
	 qSph->afElemMass[j] += fWeights[i] * afElemMass[j];
      qSph->fMetalMass += fWeights[i] * fMetalMass;
   }
}


#define H5FIELD_MASS         "Masses"
#define H5FIELD_METALLICITY  "Metallicities"
#define H5FIELD_SPECNAME     "Species_names"
#define H5FIELD_TABLEZNAME   "Yield_names"
#define H5FIELD_TABLE        "Yields"
#define H5FIELD_YIELD        "Yield"
#define H5FIELD_METALYIELD   "Total_Metals"
#define H5FIELD_EJMASS       "Ejected_mass"
#define H5FIELD_LIFETIME     "Lifetimes"


STEV_RAWDATA *stevReadTable(char *pszPath) {
   hid_t fileID, datatype, dataspace;
   herr_t status;

   STEV_RAWDATA *RawData = malloc(sizeof(STEV_RAWDATA));
   assert(RawData != NULL);

   fileID = H5Fopen(pszPath, H5F_ACC_RDONLY, H5P_DEFAULT);
   assert(fileID >= 0);

   hid_t dsMetal = H5Dopen(fileID, H5FIELD_METALLICITY, H5P_DEFAULT);
   assert(dsMetal >= 0);
   dataspace = H5Dget_space(dsMetal);
   assert(dataspace >= 0);
   RawData->nZs = H5Sget_simple_extent_npoints(dataspace);
   status = H5Sclose(dataspace);
   assert(status >= 0);

   hid_t dsMass = H5Dopen(fileID, H5FIELD_MASS, H5P_DEFAULT);
   assert(dsMass >= 0);
   dataspace = H5Dget_space(dsMass);
   assert(dataspace >= 0);
   RawData->nMasses = H5Sget_simple_extent_npoints(dataspace);
   status = H5Sclose(dataspace);
   assert(status >= 0);

   hid_t dsSpecName = H5Dopen(fileID, H5FIELD_SPECNAME, H5P_DEFAULT);
   assert(dsSpecName >= 0);
   dataspace = H5Dget_space(dsSpecName);
   assert(dataspace >= 0);
   RawData->nSpecs = H5Sget_simple_extent_npoints(dataspace);
   status = H5Sclose(dataspace);
   assert(status >= 0);
   status = H5Dclose(dsSpecName);
   assert(status >= 0);


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


   status = H5Dread(dsMetal, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->pfZs);
   assert(status >= 0);
   status = H5Dclose(dsMetal);
   assert(status >= 0);

   status = H5Dread(dsMass, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->pfMasses);
   assert(status >= 0);
   status = H5Dclose(dsMass);
   assert(status >= 0);


   hid_t dsTableZName = H5Dopen(fileID, H5FIELD_TABLEZNAME, H5P_DEFAULT);
   assert(dsTableZName >= 0);
   dataspace = H5Dget_space(dsTableZName);
   assert(dataspace >= 0);
   int nTableZNames = H5Sget_simple_extent_npoints(dataspace);
   assert(nTableZNames == RawData->nZs);
   status = H5Sclose(dataspace);
   assert(status >= 0);
   char *apchTableZNames[RawData->nZs];
   datatype = H5Dget_type(dsTableZName);
   assert(datatype >= 0);
   status = H5Dread(dsTableZName, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, apchTableZNames);
   assert(status >= 0);
   status = H5Tclose(datatype);
   assert(status >= 0);
   status = H5Dclose(dsTableZName);
   assert(status >= 0);

   char achTable[64];
   int iOffset, nDims, nSize;
   hsize_t dimSize[2];
   for (int i = 0; i < RawData->nZs; i++) {
      iOffset = i * RawData->nSpecs * RawData->nMasses;

      sprintf(achTable, "%s/%s/%s", H5FIELD_TABLE, apchTableZNames[i], H5FIELD_YIELD);
      hid_t dsYield = H5Dopen(fileID, achTable, H5P_DEFAULT);
      assert(dsYield >= 0);
      dataspace = H5Dget_space(dsYield);
      assert(dataspace >= 0);
      nDims = H5Sget_simple_extent_dims(dataspace, dimSize, NULL);
      assert(nDims == 2);
      assert(dimSize[0] == RawData->nSpecs);
      assert(dimSize[1] == RawData->nMasses);
      status = H5Sclose(dataspace);
      assert(status >= 0);
      status = H5Dread(dsYield, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		       RawData->pfYields + iOffset);
      assert(status >= 0);
      status = H5Dclose(dsYield);
      assert(status >= 0);

      iOffset = i * RawData->nMasses;

      sprintf(achTable, "%s/%s/%s", H5FIELD_TABLE, apchTableZNames[i], H5FIELD_METALYIELD);
      hid_t dsMetalYield = H5Dopen(fileID, achTable, H5P_DEFAULT);
      assert(dsMetalYield >= 0);
      dataspace = H5Dget_space(dsMetalYield);
      assert(dataspace >= 0);
      nSize = H5Sget_simple_extent_npoints(dataspace);
      assert(nSize == RawData->nMasses);
      status = H5Sclose(dataspace);
      assert(status >= 0);
      status = H5Dread(dsMetalYield, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		       RawData->pfMetalYield + iOffset);
      assert(status >= 0);
      status = H5Dclose(dsMetalYield);
      assert(status >= 0);

      sprintf(achTable, "%s/%s/%s", H5FIELD_TABLE, apchTableZNames[i], H5FIELD_EJMASS);
      hid_t dsEjMass = H5Dopen(fileID, achTable, H5P_DEFAULT);
      assert(dsEjMass >= 0);
      dataspace = H5Dget_space(dsEjMass);
      assert(dataspace >= 0);
      nSize = H5Sget_simple_extent_npoints(dataspace);
      assert(nSize == RawData->nMasses);
      status = H5Sclose(dataspace);
      assert(status >= 0);
      status = H5Dread(dsEjMass, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		       RawData->pfEjectedMass + iOffset);
      assert(status >= 0);
      status = H5Dclose(dsEjMass);
      assert(status >= 0);
   }

   status = H5Fclose(fileID);
   assert(status >= 0);

   return RawData;
}


/* This function assumes that the table of mass ejected from Type Ia SN is independent of
   the progenitor's initial mass and metallicity */
STEV_RAWDATA *stevReadSNIaTable(char *pszPath) {
   hid_t fileID, dataspace;
   herr_t status;

   STEV_RAWDATA *RawData = malloc(sizeof(STEV_RAWDATA));
   assert(RawData != NULL);

   fileID = H5Fopen(pszPath, H5F_ACC_RDONLY, H5P_DEFAULT);
   assert(fileID >= 0);

   hid_t dsYield = H5Dopen(fileID, H5FIELD_YIELD, H5P_DEFAULT);
   assert(dsYield >= 0);
   dataspace = H5Dget_space(dsYield);
   assert(dataspace >= 0);
   RawData->nSpecs = H5Sget_simple_extent_npoints(dataspace);
   status = H5Sclose(dataspace);
   assert(status >= 0);

   hid_t dsMetalYield = H5Dopen(fileID, H5FIELD_METALYIELD, H5P_DEFAULT);
   assert(dsMetalYield >= 0);
   dataspace = H5Dget_space(dsMetalYield);
   assert(dataspace >= 0);
   int nSize = H5Sget_simple_extent_npoints(dataspace);
   assert(nSize == 1);
   status = H5Sclose(dataspace);
   assert(status >= 0);


   RawData->pfEjectedMass = (float *) malloc(RawData->nSpecs * sizeof(float));
   assert(RawData->pfEjectedMass != NULL);
   RawData->pfMetalYield = (float *) malloc(sizeof(float));
   assert(RawData->pfMetalYield != NULL);


   status = H5Dread(dsYield, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    RawData->pfEjectedMass);
   assert(status >= 0);
   status = H5Dclose(dsYield);
   assert(status >= 0);

   status = H5Dread(dsMetalYield, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    RawData->pfMetalYield);
   assert(status >= 0);
   status = H5Dclose(dsMetalYield);
   assert(status >= 0);

   status = H5Fclose(fileID);
   assert(status >= 0);

   return RawData;
}


STEV_RAWDATA *stevReadLifetimesTable(char *pszPath) {
   hid_t fileID, dataspace;
   herr_t status;

   STEV_RAWDATA *RawData = malloc(sizeof(STEV_RAWDATA));
   assert(RawData != NULL);

   fileID = H5Fopen(pszPath, H5F_ACC_RDONLY, H5P_DEFAULT);
   assert(fileID >= 0);

   hid_t dsMetal = H5Dopen(fileID, H5FIELD_METALLICITY, H5P_DEFAULT);
   assert(dsMetal >= 0);
   dataspace = H5Dget_space(dsMetal);
   assert(dataspace >= 0);
   RawData->nZs = H5Sget_simple_extent_npoints(dataspace);
   status = H5Sclose(dataspace);
   assert(status >= 0);

   hid_t dsMass = H5Dopen(fileID, H5FIELD_MASS, H5P_DEFAULT);
   assert(dsMass >= 0);
   dataspace = H5Dget_space(dsMass);
   assert(dataspace >= 0);
   RawData->nMasses = H5Sget_simple_extent_npoints(dataspace);
   status = H5Sclose(dataspace);
   assert(status >= 0);


   RawData->pfZs = (float *) malloc(RawData->nZs * sizeof(float));
   assert(RawData->pfZs != NULL);
   RawData->pfMasses = (float *) malloc(RawData->nMasses * sizeof(float));
   assert(RawData->pfMasses != NULL);
   RawData->pfLifetimes = (float *) malloc(RawData->nZs * RawData->nMasses * sizeof(float));
   assert(RawData->pfLifetimes != NULL);


   status = H5Dread(dsMetal, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->pfZs);
   assert(status >= 0);
   status = H5Dclose(dsMetal);
   assert(status >= 0);

   status = H5Dread(dsMass, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->pfMasses);
   assert(status >= 0);
   status = H5Dclose(dsMass);
   assert(status >= 0);

   hsize_t dimSize[2];
   hid_t dsLifetimes = H5Dopen(fileID, H5FIELD_LIFETIME, H5P_DEFAULT);
   assert(dsLifetimes >= 0);
   dataspace = H5Dget_space(dsLifetimes);
   assert(dataspace >= 0);
   int nDims = H5Sget_simple_extent_dims(dataspace, dimSize, NULL);
   assert(nDims == 2);
   assert(dimSize[0] == RawData->nZs);
   assert(dimSize[1] == RawData->nMasses);
   status = H5Sclose(dataspace);
   assert(status >= 0);
   status = H5Dread(dsLifetimes, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
		    RawData->pfLifetimes);
   assert(status >= 0);
   status = H5Dclose(dsLifetimes);
   assert(status >= 0);

   status = H5Fclose(fileID);
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

   return (float)pkd->param.dSNIa_Norm *
          (expf(-(fInitialTime / (float)pkd->param.dSNIa_Scale)) -
           expf(-(fFinalTime / (float)pkd->param.dSNIa_Scale)));
}


float stevPowerlawNumSNIa(PKD pkd, STARFIELDS *pStar, float fInitialTime, float fFinalTime) {
   if (fFinalTime <= pStar->fSNIaOnsetTime) {
      return 0.0f;
   }
   else if (fInitialTime < pStar->fSNIaOnsetTime) {
      fInitialTime = pStar->fSNIaOnsetTime;
   }

   return (float)pkd->param.dSNIa_Norm *
          (powf(fFinalTime, (float)pkd->param.dSNIa_Scale + 1.0f) -
	   powf(fInitialTime, (float)pkd->param.dSNIa_Scale + 1.0f));
}


#endif	/* STELLAR_EVOLUTION */
