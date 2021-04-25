#ifdef STELLAR_EVOLUTION

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <hdf5.h>

#include "stellarevolution.h"


void msrStellarEvolutionInit(MSR msr) {
   char achPath[256];


   STEV_RAWDATA *CCSNdata;
   CCSNdata = (STEV_RAWDATA *) malloc(sizeof(STEV_RAWDATA));
   assert(CCSNdata != NULL);

   sprintf(achPath, "%s/SNII.hdf5", msr->param.achStEvolPath);
   char *pszCCSNFields[] = {"Z_0.0004", "Z_0.004", "Z_0.008", "Z_0.02", "Z_0.05"};
   int nCCSNFields = 5;

   stevReadTable(CCSNdata, achPath, pszCCSNFields, nCCSNFields);

   assert(CCSNdata->nZs == STEV_CCSN_N_METALLICITY);
   assert(CCSNdata->nMasses == STEV_CCSN_N_MASS);
   assert(CCSNdata->nSpecs == STEV_N_ELEM);


   STEV_RAWDATA *AGBdata;
   AGBdata = (STEV_RAWDATA *) malloc(sizeof(STEV_RAWDATA));
   assert(AGBdata != NULL);

   sprintf(achPath, "%s/AGB.hdf5", msr->param.achStEvolPath);
   char *pszAGBFields[] = {"Z_0.004", "Z_0.008", "Z_0.019"};
   int nAGBFields = 3;

   stevReadTable(AGBdata, achPath, pszAGBFields, nAGBFields);

   assert(AGBdata->nZs == STEV_AGB_N_METALLICITY);
   assert(AGBdata->nMasses == STEV_AGB_N_MASS);
   assert(AGBdata->nSpecs == STEV_N_ELEM);


   STEV_RAWDATA *SNIadata;
   SNIadata = (STEV_RAWDATA *) malloc(sizeof(STEV_RAWDATA));
   assert(SNIadata != NULL);

   sprintf(achPath, "%s/SNIa.hdf5", msr->param.achStEvolPath);
   stevReadSNIaTable(SNIadata, achPath);

   assert(SNIadata->nSpecs == STEV_N_ELEM);


   STEV_RAWDATA *LifetimesData;
   LifetimesData = (STEV_RAWDATA *) malloc(sizeof(STEV_RAWDATA));
   assert(LifetimesData != NULL);

   sprintf(achPath, "%s/Lifetimes.hdf5", msr->param.achStEvolPath);
   stevReadLifetimesTable(LifetimesData, achPath);

   assert(LifetimesData->nZs == STEV_LIFETIMES_N_METALLICITY);
   assert(LifetimesData->nMasses == STEV_LIFETIMES_N_MASS);


   status = H5Fclose(file_id);
   assert(status >= 0);

   stevFreeTable(CCSNdata);
   stevFreeTable(AGBdata);
   stevFreeSNIaTable(SNIadata);
   stevFreeLifetimesTable(LifetimesData);
}


void stevReadTable(STEV_RAWDATA *RawData, char *pszPath, char **pszFields, int nFields) {
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


   RawData->Zs = (float *) malloc(RawData->nZs * sizeof(float));
   assert(RawData->Zs != NULL);
   RawData->Masses = (float *) malloc(RawData->nMasses * sizeof(float));
   assert(RawData->Masses != NULL);
   RawData->Yields = (float *) malloc(RawData->nZs * RawData->nSpecs * RawData->nMasses * sizeof(float));
   assert(RawData->Yields != NULL);
   RawData->MetalYield = (float *) malloc(RawData->nZs * RawData->nMasses * sizeof(float));
   assert(RawData->MetalYield != NULL);
   RawData->EjectedMass = (float *) malloc(RawData->nZs * RawData->nMasses * sizeof(float));
   assert(RawData->EjectedMass != NULL);


   dataset = H5Dopen(file_id, "/Metallicities", H5P_DEFAULT);
   datatype = H5Dget_type(dataset);
   assert(H5Tequal(datatype, H5T_NATIVE_FLOAT));
   status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->Zs);
   assert(status >= 0);
   status = H5Tclose(datatype);
   assert(status >= 0);
   status = H5Dclose(dataset);
   assert(status >= 0);

   dataset = H5Dopen(file_id, "/Masses", H5P_DEFAULT);
   datatype = H5Dget_type(dataset);
   assert(H5Tequal(datatype, H5T_NATIVE_FLOAT));
   status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->Masses);
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
      status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->Yields + n);
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
      status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->MetalYield + n);
      assert(status >= 0);
      status = H5Tclose(datatype);
      assert(status >= 0);
      status = H5Dclose(dataset);
      assert(status >= 0);

      sprintf(achField, "/Yields/%s/Ejected_mass", pszFields[i]);
      dataset = H5Dopen(file_id, achField, H5P_DEFAULT);
      datatype = H5Dget_type(dataset);
      assert(H5Tequal(datatype, H5T_NATIVE_FLOAT));
      status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->EjectedMass + n);
      assert(status >= 0);
      status = H5Tclose(datatype);
      assert(status >= 0);
      status = H5Dclose(dataset);
      assert(status >= 0);
   }

   status = H5Fclose(file_id);
   assert(status >= 0);
}


void stevFreeTable(STEV_RAWDATA *RawData) {
   free(RawData->Zs);
   free(RawData->Masses);
   free(RawData->Yields);
   free(RawData->MetalYield);
   free(RawData->EjectedMass);
   free(RawData);
}


/* This function assumes that the table of mass ejected from Type Ia SN is independent of
   the progenitor's initial mass and metallicity */
void stevReadSNIaTable(STEV_RAWDATA *RawData, char *pszPath) {
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


   RawData->MetalYield = (float *) malloc(sizeof(float));
   assert(RawData->MetalYield != NULL);
   RawData->EjectedMass = (float *) malloc(RawData->nSpecs * sizeof(float));
   assert(RawData->EjectedMass != NULL);


   dataset = H5Dopen(file_id, "/Yield", H5P_DEFAULT);
   datatype = H5Dget_type(dataset);
   assert(H5Tequal(datatype, H5T_NATIVE_FLOAT));
   status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->EjectedMass);
   assert(status >= 0);
   status = H5Tclose(datatype);
   assert(status >= 0);
   status = H5Dclose(dataset);
   assert(status >= 0);

   dataset = H5Dopen(file_id, "/Total_Metals", H5P_DEFAULT);
   datatype = H5Dget_type(dataset);
   assert(H5Tequal(datatype, H5T_NATIVE_FLOAT));
   status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->MetalYield);
   assert(status >= 0);
   status = H5Tclose(datatype);
   assert(status >= 0);
   status = H5Dclose(dataset);
   assert(status >= 0);


   status = H5Fclose(file_id);
   assert(status >= 0);
}


void stevFreeSNIaTable(STEV_RAWDATA *RawData) {
   free(RawData->MetalYield);
   free(RawData->EjectedMass);
   free(RawData);
}


void stevReadLifetimesTable(STEV_RAWDATA *RawData, char *pszPath) {
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


   RawData->Zs = (float *) malloc(RawData->nZs * sizeof(float));
   assert(RawData->Zs != NULL);
   RawData->Masses = (float *) malloc(RawData->nMasses * sizeof(float));
   assert(RawData->Masses != NULL);
   RawData->Lifetimes = (float *) malloc(RawData->nZs * RawData->nMasses * sizeof(float));
   assert(RawData->Lifetimes != NULL);


   dataset = H5Dopen(file_id, "/Metallicities", H5P_DEFAULT);
   datatype = H5Dget_type(dataset);
   assert(H5Tequal(datatype, H5T_NATIVE_FLOAT));
   status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->Zs);
   assert(status >= 0);
   status = H5Tclose(datatype);
   assert(status >= 0);
   status = H5Dclose(dataset);
   assert(status >= 0);

   dataset = H5Dopen(file_id, "/Masses", H5P_DEFAULT);
   datatype = H5Dget_type(dataset);
   assert(H5Tequal(datatype, H5T_NATIVE_FLOAT));
   status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->Masses);
   assert(status >= 0);
   status = H5Tclose(datatype);
   assert(status >= 0);
   status = H5Dclose(dataset);
   assert(status >= 0);

   dataset = H5Dopen(file_id, "/Lifetimes", H5P_DEFAULT);
   datatype = H5Dget_type(dataset);
   assert(H5Tequal(datatype, H5T_NATIVE_FLOAT));
   status = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->Lifetimes);
   assert(status >= 0);
   status = H5Tclose(datatype);
   assert(status >= 0);
   status = H5Dclose(dataset);
   assert(status >= 0);


   status = H5Fclose(file_id);
   assert(status >= 0);
}


void stevFreeLifetimesTable(STEV_RAWDATA *RawData) {
   free(RawData->Zs);
   free(RawData->Masses);
   free(RawData->Lifetimes);
   free(RawData);
}


#endif	/* STELLAR_EVOLUTION */
