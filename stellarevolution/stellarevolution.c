#ifdef STELLAR_EVOLUTION


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <hdf5.h>

#include "stellarevolution.h"
#include "interpolate.h"


void msrStellarEvolutionInit(MSR msr) {
   assert(STEV_N_ELEM == chemistry_element_count);

   int i;
   char achPath[256];
   STEV_RAWDATA *CCSNdata, *AGBdata, *SNIadata, *LifetimesData;


   sprintf(achPath, "%s/SNII.hdf5", msr->param.achStEvolPath);
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

   double MinMass =  msr->param.dIMF_MinMass;
   double MaxMass =  msr->param.dIMF_MaxMass;
   double Delta = pow(10.0, (log10(MaxMass) - log10(MinMass)) / (STEV_INTERP_N_MASS - 1));
   double MassVal;
   for (i = 0, MassVal = MinMass; i < STEV_INTERP_N_MASS; i++, MassVal *= Delta)
      buffer.afMasses[i] = MassVal;

   if (strcmp(msr->param.achIMFtype, "chabrier") == 0) {
      stevChabrierIMF(buffer.afMasses, STEV_INTERP_N_MASS, MinMass, MaxMass, buffer.afIMF);
   }
   else {
      printf("ERROR: Undefined IMF type has been given in achIMFtype parameter: %s\n",
	     msr->param.achIMFtype);
      abort();
   }

   float IMFnorm = 0.0f;
   for (i = 0; i < STEV_INTERP_N_MASS - 1; i++) {
      IMFnorm += (buffer.afMasses[i + 1] * buffer.afIMF[i + 1] + buffer.afMasses[i] *
		  buffer.afIMF[i]) * (buffer.afMasses[i + 1] - buffer.afMasses[i]);
   }
   IMFnorm *= 0.5f;
   printf("IMF normalization is: %.4f\n", IMFnorm);
   if (fabs(IMFnorm - 1.0) > 1e-2) {
      printf("ERROR: IMF normalization differs from unity by more than 1%%\n");
      abort();
   }

   /* WARNING: This assumes table masses in log and buffer masses not in log */
   int iM, j, k, idxTable, idxBuffer;
   float dM, Mass, Mtrans;
   Mtrans = msr->param.dCCSN_MinMass;
   for (i = 0; i < STEV_INTERP_N_MASS; i++) {
      Mass = buffer.afMasses[i];
      if (Mass < Mtrans) {
	 stevGetIndex1D(AGBdata->pfMasses, STEV_AGB_N_MASS, log10(Mass), &iM, &dM);

	 for (j = 0; j < STEV_N_ELEM; j++) {
	    for (k = 0; k < STEV_CCSN_N_METALLICITY; k++) {
	       idxBuffer = row_major_index_3d(j, k, i, STEV_N_ELEM, STEV_CCSN_N_METALLICITY,
					      STEV_INTERP_N_MASS);

	       buffer.afCCSN_Yields[idxBuffer] = 0.0f;
	    }

	    for (k = 0; k < STEV_AGB_N_METALLICITY; k++) {
	       idxTable  = row_major_index_3d(k, j, iM, STEV_AGB_N_METALLICITY, STEV_N_ELEM,
					      STEV_AGB_N_MASS);
	       idxBuffer = row_major_index_3d(j, k,  i, STEV_N_ELEM, STEV_AGB_N_METALLICITY,
					      STEV_INTERP_N_MASS);

	       /* WARNING: The following formula is correct only when mass is the last index
		  in the tables */
	       buffer.afAGB_Yields[idxBuffer] = AGBdata->pfYields[idxTable] * (1.0 - dM) + \
		                                AGBdata->pfYields[idxTable + 1] * dM;
	    }
	 }

	 for (k = 0; k < STEV_CCSN_N_METALLICITY; k++) {
	    idxBuffer = row_major_index_2d(k, i, STEV_CCSN_N_METALLICITY, STEV_INTERP_N_MASS);

	    buffer.afCCSN_MetalYield[idxBuffer] = 0.0f;
	    buffer.afCCSN_EjectedMass[idxBuffer] = 0.0f;
	 }

	 for (k = 0; k < STEV_AGB_N_METALLICITY; k++) {
	    idxTable  = row_major_index_2d(k, iM, STEV_AGB_N_METALLICITY, STEV_AGB_N_MASS);
	    idxBuffer = row_major_index_2d(k,  i, STEV_AGB_N_METALLICITY, STEV_INTERP_N_MASS);

	    /* WARNING: The following formulas are correct only when mass is the last index
	       in the tables */
	    buffer.afAGB_MetalYield[idxBuffer] = AGBdata->pfMetalYield[idxTable] * (1.0 - dM) + \
		                                 AGBdata->pfMetalYield[idxTable + 1] * dM;
	    buffer.afAGB_EjectedMass[idxBuffer] = AGBdata->pfEjectedMass[idxTable] * (1.0 - dM) + \
		                                  AGBdata->pfEjectedMass[idxTable + 1] * dM;
	 }
      }
      else {
	 stevGetIndex1D(CCSNdata->pfMasses, STEV_CCSN_N_MASS, log10(Mass), &iM, &dM);

	 for (j = 0; j < STEV_N_ELEM; j++) {
	    for (k = 0; k < STEV_CCSN_N_METALLICITY; k++) {
	       idxTable  = row_major_index_3d(k, j, iM, STEV_CCSN_N_METALLICITY, STEV_N_ELEM,
					      STEV_CCSN_N_MASS);
	       idxBuffer = row_major_index_3d(j, k,  i, STEV_N_ELEM, STEV_CCSN_N_METALLICITY,
					      STEV_INTERP_N_MASS);

	       /* WARNING: The following formula is correct only when mass is the last index
		  in the tables */
	       buffer.afCCSN_Yields[idxBuffer] = CCSNdata->pfYields[idxTable] * (1.0 - dM) + \
		                                 CCSNdata->pfYields[idxTable + 1] * dM;
	    }

	    for (k = 0; k < STEV_AGB_N_METALLICITY; k++) {
	       idxBuffer = row_major_index_3d(j, k, i, STEV_N_ELEM, STEV_AGB_N_METALLICITY,
					      STEV_INTERP_N_MASS);

	       buffer.afAGB_Yields[idxBuffer] = 0.0f;
	    }
	 }

	 for (k = 0; k < STEV_CCSN_N_METALLICITY; k++) {
	    idxTable  = row_major_index_2d(k, iM, STEV_CCSN_N_METALLICITY, STEV_CCSN_N_MASS);
	    idxBuffer = row_major_index_2d(k,  i, STEV_CCSN_N_METALLICITY, STEV_INTERP_N_MASS);

	    /* WARNING: The following formulas are correct only when mass is the last index
	       in the tables */
	    buffer.afCCSN_MetalYield[idxBuffer] = CCSNdata->pfMetalYield[idxTable] * (1.0 - dM) + \
		                                  CCSNdata->pfMetalYield[idxTable + 1] * dM;
	    buffer.afCCSN_EjectedMass[idxBuffer] = CCSNdata->pfEjectedMass[idxTable] * (1.0 - dM) + \
		                                   CCSNdata->pfEjectedMass[idxTable + 1] * dM;
	 }

	 for (k = 0; k < STEV_AGB_N_METALLICITY; k++) {
	    idxBuffer = row_major_index_2d(k, i, STEV_AGB_N_METALLICITY, STEV_INTERP_N_MASS);

	    buffer.afAGB_MetalYield[idxBuffer] = 0.0f;
	    buffer.afAGB_EjectedMass[idxBuffer] = 0.0f;
	 }
      }
   }


   for (i = 0; i < STEV_N_ELEM * STEV_CCSN_N_METALLICITY * STEV_INTERP_N_MASS; i++)
      buffer.afCCSN_Yields[i] /= msr->param.dMsolUnit;
   for (i = 0; i < STEV_CCSN_N_METALLICITY * STEV_INTERP_N_MASS; i++) {
      buffer.afCCSN_MetalYield[i]  /= msr->param.dMsolUnit;
      buffer.afCCSN_EjectedMass[i] /= msr->param.dMsolUnit;
   }
   for (i = 0; i < STEV_CCSN_N_METALLICITY; i++)
      buffer.afCCSN_Zs[i] = CCSNdata->pfZs[i];
   
   for (i = 0; i < STEV_N_ELEM * STEV_AGB_N_METALLICITY * STEV_INTERP_N_MASS; i++)
      buffer.afAGB_Yields[i] /= msr->param.dMsolUnit;
   for (i = 0; i < STEV_AGB_N_METALLICITY * STEV_INTERP_N_MASS; i++) {
      buffer.afAGB_MetalYield[i]  /= msr->param.dMsolUnit;
      buffer.afAGB_EjectedMass[i] /= msr->param.dMsolUnit;
   }
   for (i = 0; i < STEV_AGB_N_METALLICITY; i++)
      buffer.afAGB_Zs[i] = AGBdata->pfZs[i];

   for (i = 0; i < STEV_N_ELEM; i++)
      buffer.afSNIa_EjectedMass[i] = SNIadata->pfEjectedMass[i] / msr->param.dMsolUnit;
   buffer.fSNIa_EjectedMetalMass = *SNIadata->pfMetalYield / msr->param.dMsolUnit;

   for (i = 0; i < STEV_LIFETIMES_N_METALLICITY; i++)
      buffer.afLifetimes_Zs[i] = LifetimesData->pfZs[i];
   for (i = 0; i < STEV_LIFETIMES_N_MASS; i++)
      buffer.afLifetimes_Masses[i] = LifetimesData->pfMasses[i] / msr->param.dMsolUnit;
   for (i = 0; i < STEV_LIFETIMES_N_METALLICITY * STEV_LIFETIMES_N_MASS; i++)
      buffer.afLifetimes[i] = LifetimesData->pfLifetimes[i] * SECPERYEAR / msr->param.dSecUnit;


   pstStellarEvolutionInit(msr->pst, (void *) &buffer, sizeof(buffer), NULL, 0);


   stevFreeTable(CCSNdata);
   stevFreeTable(AGBdata);
   stevFreeSNIaTable(SNIadata);
   stevFreeLifetimesTable(LifetimesData);
}


int pstStellarEvolutionInit(PST pst, void *vin, int nIn, void *vout, int nOut) {
   mdlassert(pst->mdl, nIn == sizeof(struct inStellarEvolution));
   if (pst->nLeaves > 1) {
      int rID = mdlReqService(pst->mdl, pst->idUpper, STELLAREVOLUTIONINIT, vin, nIn);
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

#ifndef COOLING
   for (int i = 0; i < pkd->nLocal; ++i) {
      PARTICLE *p = pkdParticle(pkd,i);
      if (pkdIsGas(pkd, p)) {
         SPHFIELDS *psph = pkdSph(pkd,p);
	 float fMass = pkdMass(pkd, p);
         if (psph->chemistry[chemistry_element_H] == 0.0f) { /* TODO: Why this?????? */
            psph->chemistry[chemistry_element_H]  = pkd->param.dInitialH  * fMass;
            psph->chemistry[chemistry_element_He] = pkd->param.dInitialHe * fMass;
            psph->chemistry[chemistry_element_C]  = pkd->param.dInitialC  * fMass;
            psph->chemistry[chemistry_element_N]  = pkd->param.dInitialN  * fMass;
            psph->chemistry[chemistry_element_O]  = pkd->param.dInitialO  * fMass;
            psph->chemistry[chemistry_element_Ne] = pkd->param.dInitialNe * fMass;
            psph->chemistry[chemistry_element_Mg] = pkd->param.dInitialMg * fMass;
            psph->chemistry[chemistry_element_Si] = pkd->param.dInitialSi * fMass;
            psph->chemistry[chemistry_element_Fe] = pkd->param.dInitialFe * fMass;
         }
      }
   }
#endif

   return 0;
}


void initChemEnrich(void *vpkd, void *vp) {
   
}


void combChemEnrich(void *vpkd, void *vp1, void *vp2) {
   
}


void smChemEnrich(PARTICLE *p, float fBall, int nSmooth, NN *nnList, SMF *smf) {
   
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


/* This function assumes that the compiler can find the erf() function and the M_PI value
   of pi. In POSIX-compliant systems, math.h should declare them */
void stevChabrierIMF(float *pfMasses, int nSize, double fMinMass, double fMaxMass, float *pfIMF) {
   double Mc = 0.079;
   double Sigma = 0.69;
   double Slope = -2.3;

   double logMc = log10(Mc);
   double Sigma2 = pow(Sigma, 2.0);
   double C = exp(-0.5 * pow(-logMc, 2.0) / Sigma2);
   double ln10 = log(10.0);

   double I1 = sqrt(0.5 * M_PI) * exp(0.5 * Sigma2 * pow(ln10, 2.0)) * ln10 * Sigma * Mc * \
      (erf((-logMc / Sigma - ln10 * Sigma) / sqrt(2.0)) -
       erf(((log10(fMinMass) - logMc) / Sigma - ln10 * Sigma) / sqrt(2.0)));
   double I2 = C * (pow(fMaxMass, Slope + 2.0) - 1.0) / (Slope + 2.0);

   double A = 1.0 / (I1 + I2);
   double B = A * C;
   for (int i = 0; i < nSize; i++) {
      if (pfMasses[i] < 1.0)
	 pfIMF[i] = A * exp(-0.5 * pow(log10(pfMasses[i]) - logMc, 2.0) / Sigma2) / pfMasses[i];
      else
	 pfIMF[i] = B * pow(pfMasses[i], Slope);
   }
}


float stevExponentialNumSNIa(PKD pkd, STARFIELDS *pStar, float fTime) {
   float ti = pStar->fLastEnrichTime;
   float tf = fTime - pStar->fTimer;

   if (tf <= pStar->fTimeSNIa) {
      return 0.0f;
   }
   else if (ti < pStar->fTimeSNIa) {
      ti = pStar->fTimeSNIa;
   }

   return pStar->fInitialMass * pkd->param.dSNIa_Norm * (exp(-(ti / pkd->param.dSNIa_Scale)) -
							 exp(-(tf / pkd->param.dSNIa_Scale)));
}


float stevPowerlawNumSNIa(PKD pkd, STARFIELDS *pStar, float fTime) {
   float ti = pStar->fLastEnrichTime;
   float tf = fTime - pStar->fTimer;

   if (tf <= pStar->fTimeSNIa) {
      return 0.0f;
   }
   else if (ti < pStar->fTimeSNIa) {
      ti = pStar->fTimeSNIa;
   }

   return pStar->fInitialMass * pkd->param.dSNIa_Norm * (pow(tf, pkd->param.dSNIa_Scale + 1.0) -
							 pow(ti, pkd->param.dSNIa_Scale + 1.0));
}


#endif	/* STELLAR_EVOLUTION */
