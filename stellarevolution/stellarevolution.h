#ifdef STELLAR_EVOLUTION

#ifndef MASTER_HINCLUDED
#include "master.h"
#endif


/* NOTAS:
   - En EAGLE se impone conservacion de momento y energia al distribuir la masa (Schaye 2015).
   - Arrays del buffer que estan en logaritmo: 
      - CCSN y AGB: Metalicidades (array de masas tanto en lineal como en log)
      - Lifetimes: Todo
   - Masas se estan interpolando en "semilog", Lifetimes en "loglog".
   - El array de masas inicia en 0.7 Mo. La IMF se normaliza en su rango normal, [0.1,100].
   Esto significa que no se puede verificar la normalizacion usando el array de masas.
   - Se esta siendo poco preciso en la integracion alrededor de la masa de transicion CCSN-AGB.
*/


/* TODOs y DUDAS:
   - Implementar la inicializacion cuando se empieza de un archivo restart o de una snapshot
   - Se puede mejorar la lectura de las tablas. Los nombres de los datasets para cada 
   metalicidad estan especificados bajo Yield_names. Se puede leer el tamano de un dataset 
   directamente?
   - Considerar interpolar usando gsl

   - Llamada redundante a msrActiveRung(msr,uRung,1) en msrNewTopStepKDK despues de hacer
   el feedback
   - Por que no estoy utilizando las funciones de interpolacion de interpolate.h?
*/


#define STEV_CCSN_N_METALLICITY               5
#define STEV_CCSN_N_MASS                     11
#define STEV_AGB_N_METALLICITY                3
#define STEV_AGB_N_MASS                      23
#define STEV_LIFETIMES_N_METALLICITY          6
#define STEV_LIFETIMES_N_MASS                30
#define STEV_N_ELEM                           9

#define STEV_INTERP_N_MASS                  200

#define STEV_MIN_LOG_METALLICITY            -20


/*
 * ---------------------
 * STRUCTURE DEFINITIONS
 * --------------------- 
 */

typedef struct inStellarEvolution {
   /* Pointer to function that gives the number of SNIa per Msol in [fInitialTime, fFinalTime] */
   float (*fcnNumSNIa)(PKD pkd, STARFIELDS *pStar, float fInitialTime, float fFinalTime);

   /* Initial mass arrays for CCSN and AGB tables */
   float afMasses[STEV_INTERP_N_MASS];
   float afLogMasses[STEV_INTERP_N_MASS];

   /* Initial Mass Function values array */
   float afIMF[STEV_INTERP_N_MASS];

   /* Core Collapse SNe arrays */
   float afCCSN_Zs[STEV_CCSN_N_METALLICITY];
   float afCCSN_Yields[STEV_CCSN_N_METALLICITY * STEV_INTERP_N_MASS * STEV_N_ELEM];
   float afCCSN_MetalYield[STEV_CCSN_N_METALLICITY * STEV_INTERP_N_MASS];
   float afCCSN_EjectedMass[STEV_CCSN_N_METALLICITY * STEV_INTERP_N_MASS];

   /* AGB arrays */
   float afAGB_Zs[STEV_AGB_N_METALLICITY];
   float afAGB_Yields[STEV_AGB_N_METALLICITY * STEV_INTERP_N_MASS * STEV_N_ELEM];
   float afAGB_MetalYield[STEV_AGB_N_METALLICITY * STEV_INTERP_N_MASS];
   float afAGB_EjectedMass[STEV_AGB_N_METALLICITY * STEV_INTERP_N_MASS];

   /* Type Ia SNe */
   float afSNIa_EjectedMass[STEV_N_ELEM];
   float fSNIa_EjectedMetalMass;

   /* Stellar lifetimes arrays */
   float afLifetimes_Zs[STEV_LIFETIMES_N_METALLICITY];
   float afLifetimes_Masses[STEV_LIFETIMES_N_MASS];
   float afLifetimes[STEV_LIFETIMES_N_METALLICITY * STEV_LIFETIMES_N_MASS];
} STEV_DATA;


typedef struct StellarEvolutionRawData {
   int nZs;
   int nSpecs;
   int nMasses;
   float *pfZs;
   float *pfMasses;
   float *pfYields;
   float *pfMetalYield;
   float *pfEjectedMass;
   float *pfLifetimes;
} STEV_RAWDATA;



/*
 * --------------
 * MAIN FUNCTIONS
 * --------------
 */

void msrStellarEvolutionInit(MSR);
int pstStellarEvolutionInit(PST, void *, int, void *, int);
int pkdStellarEvolutionInit(PKD, STEV_DATA *);

void smChemEnrich(PARTICLE *p, float fBall, int nSmooth, NN *nnList, SMF *smf);
void initChemEnrich(void *vpkd, void *vp);
void combChemEnrich(void *vpkd, void *vp1, void *vp2);



/*
 * ----------------
 * HELPER FUNCTIONS
 * ----------------
 */

STEV_RAWDATA *stevReadTable(char *, char **, int);
STEV_RAWDATA *stevReadSNIaTable(char *);
STEV_RAWDATA *stevReadLifetimesTable(char *);
void stevFreeTable(STEV_RAWDATA *);
void stevFreeSNIaTable(STEV_RAWDATA *);
void stevFreeLifetimesTable(STEV_RAWDATA *);

float stevExponentialNumSNIa(PKD, STARFIELDS *, float, float);
float stevPowerlawNumSNIa(PKD, STARFIELDS *, float, float);



/*
 * -------
 * INLINES
 * -------
 */

static inline void stevPopulateInitialMassArrays(const double dMinMass, const double dMaxMass,
						 const int nSize, float *restrict pfMasses,
						 float *restrict pfLogMasses) {

   const double dDeltaLog = (log10(dMaxMass) - log10(dMinMass)) / (nSize - 1);
   const double dDelta = pow(10.0, dDeltaLog);

   const float *pfEnd = pfMasses + nSize;
   double dMass = dMinMass, dLogMass = log10(dMinMass);
   while (pfMasses < pfEnd) {
      *pfMasses++ = dMass;
      *pfLogMasses++ = dLogMass;
      dMass *= dDelta;
      dLogMass += dDeltaLog;
   }
}


static inline void stevChabrierIMF(const float *restrict pfMasses, const int nSize,
				   const double dMinMass, const double dMaxMass,
				   float *restrict pfIMF) {
   const double Mc = 0.079;
   const double Sigma = 0.69;
   const double Slope = -2.3;

   const double logMc = log10(Mc);
   const double Sigma2 = Sigma * Sigma;
   const double C = exp(-0.5 * logMc * logMc / Sigma2);

   const double I1 = sqrt(0.5 * M_PI) * exp(0.5 * Sigma2 * M_LN10 * M_LN10) * M_LN10 *
      Sigma * Mc * (erf((-logMc / Sigma - M_LN10 * Sigma) / sqrt(2.0)) -
		    erf(((log10(dMinMass) - logMc) / Sigma - M_LN10 * Sigma) / sqrt(2.0)));
   const double I2 = C * (pow(dMaxMass, Slope + 2.0) - 1.0) / (Slope + 2.0);

   const double A = 1.0 / (I1 + I2);
   const double B = A * C;
   for (int i = 0; i < nSize; i++) {
      if (pfMasses[i] < 1.0)
	 pfIMF[i] = A * exp(-0.5 * pow(log10(pfMasses[i]) - logMc, 2.0) / Sigma2) / pfMasses[i];
      else
	 pfIMF[i] = B * pow(pfMasses[i], Slope);
   }
}


static inline void stevGetIndex1D(const float *restrict pfTable, const int nSize,
				  const float fVal, int *restrict piIdx,
				  float *restrict pfDelta) {

   const float epsilon = 1e-4f;

   if (fVal < pfTable[0] + epsilon) {
      /* We are below the first element */
      *piIdx = 0;
      *pfDelta = 0.0f;
   }
   else if (fVal > pfTable[nSize - 1] - epsilon) {
      /* We are beyond the last element */
      *piIdx = nSize - 2;
      *pfDelta = 1.0f;
   }
   else {
      /* Normal case */
      int i;
      for (i = 1; (i < nSize - 1) && (fVal > pfTable[i]); i++)
	 ;

      *piIdx = --i;
      *pfDelta = (fVal - pfTable[i]) / (pfTable[i + 1] - pfTable[i]);
   }
}


static inline void stevGetIndex1DReversed(const float *restrict pfTable, const int nSize,
					  const float fVal, int *restrict piIdx,
					  float *restrict pfDelta) {

   const float epsilon = 1e-4f;

   if (fVal > pfTable[0] - epsilon) {
      /* We are above the first element */
      *piIdx = 0;
      *pfDelta = 0.0f;
   }
   else if (fVal < pfTable[nSize - 1] + epsilon) {
      /* We are below the last element */
      *piIdx = nSize - 2;
      *pfDelta = 1.0f;
   }
   else {
      /* Normal case */
      int i;
      for (i = nSize - 2; i > 0 && fVal > pfTable[i]; i--)
	 ;

      *piIdx = i;
      *pfDelta = (pfTable[i] - fVal) / (pfTable[i] - pfTable[i + 1]);
   }
}


static inline int stevRowMajorIndex(const int iX, const int iY, const int iZ,
				    const int nX, const int nY, const int nZ) {

   return iX * nY * nZ + iY * nZ + iZ;
}


static inline int stevGetIMFMassIndex(const float *restrict pfTable, const int nSize,
				      const float fVal, const int iStart) {

   assert(iStart > 0 && iStart < nSize);
   const float *pfTemp;
   for (pfTemp = pfTable + iStart - 1; pfTemp > pfTable && *pfTemp > fVal; pfTemp--)
      ;

   return pfTemp - pfTable;
}


/* Function to interpolate nSize numbers from the array pfTable along its rows (first or 
   X- axis). The argument iOffset is assumed to represent iX * nY + oY, where iX is the
   lower X-axis index for the interpolation and oY is the Y-axis' offset. When this is
   multiplied by nZ, it gives the displacement necessary to move from the beginning 
   of the array pfTable to the first value of the lower X-axis needed for the 
   interpolation. See the note in the function smChemEnrich. */
static inline void stevInterpolateXAxis(const float *restrict pfTable, const int nY,
					const int nZ, const int nSize, const int iOffset,
					const float fWeight, float *restrict pfResult) {

   const float *pfLower = pfTable + iOffset * nZ;
   const float *pfUpper = pfLower + nY * nZ;

   const float *pfLowerEnd = pfLower + nSize * nZ;
   while(pfLower < pfLowerEnd)
      *pfResult++ = *pfLower++ * (1.0f - fWeight) + *pfUpper++ * fWeight;
}


static inline void stevComputeAndCorrectSimulEjecta(
		       const float *restrict pfYields,
		       const float fMetalYield,
		       const float fEjectedMass,
		       const int nElems, const float *restrict pfElemAbun,
		       const float fMetalAbun, float *restrict pfElemEjMass,
		       float *restrict pfMetalEjMass) {

   int i;
   if (fEjectedMass > 0.0f) {
      *pfMetalEjMass = fMetalYield + fMetalAbun * fEjectedMass;

      for (i = 0; i < nElems; i++) {
	 pfElemEjMass[i] = pfYields[i] + pfElemAbun[i] * fEjectedMass;
	 if (pfElemEjMass[i] < 0.0f) {
	    if (i != ELEMENT_H && i != ELEMENT_He) *pfMetalEjMass -= pfElemEjMass[i];
	    pfElemEjMass[i] = 0.0f;
	 }
      }
      if (*pfMetalEjMass < 0.0f) *pfMetalEjMass = 0.0f;

      float fTotalMass = pfElemEjMass[ELEMENT_H] + pfElemEjMass[ELEMENT_He] + *pfMetalEjMass;
      assert(fTotalMass > 0.0f);
      float fNormFactor = fEjectedMass / fTotalMass;
      for (i = 0; i < nElems; i++) pfElemEjMass[i] *= fNormFactor;
      *pfMetalEjMass *= fNormFactor;
   }
   else {
      for (i = 0; i < nElems; i++) pfElemEjMass[i] = 0.0f;
      *pfMetalEjMass = 0.0f;
   }
}


static inline void stevComputeMassToEject(
		       const float *restrict pfYields,
		       const float *restrict pfMetalYield,
		       const float *restrict pfEjectedMass,
		       const float *restrict pfMasses,
		       const float *restrict pfLogMasses,
		       const float *restrict pfIMF,
		       const int nZs, const int nMasses, const int nElems,
		       const float *restrict pfElemAbun, const float fMetalAbun,
		       const int iStart, const int iEnd,
		       const float fMassStart, const float fMassEnd,
		       const int iOffset, const float fWeight,
		       float *pfElemMass, float *pfMetalMass) {

   int i, j, k;
   const int nSize = iEnd - iStart + 1;
   float afInterpYields[nSize * nElems];
   float afInterpMetalYield[nSize];
   float afInterpEjectedMass[nSize];

   stevInterpolateXAxis(pfYields, nMasses, nElems, nSize, iOffset + iStart, fWeight,
			afInterpYields);
   stevInterpolateXAxis(pfMetalYield, nMasses, 1, nSize, iOffset + iStart, fWeight,
			afInterpMetalYield);
   stevInterpolateXAxis(pfEjectedMass, nMasses, 1, nSize, iOffset + iStart, fWeight,
			afInterpEjectedMass);

   /* In a strict sense, this two new declarations are not needed. We could overwrite
      two of the previously declared arrays, although the code would become less readable.
      If this were to be done, the restrict qualifier in the first and last two arguments 
      of stevComputeAndCorrectSimulEjecta must be removed */
   float afElemSimulEjecta[nSize * nElems];
   float afMetalSimulEjecta[nSize];
   for (i = 0, j = 0; i < nSize; i++, j += nElems) {
      stevComputeAndCorrectSimulEjecta(afInterpYields + j, afInterpMetalYield[i],
				       afInterpEjectedMass[i], nElems, pfElemAbun, fMetalAbun,
				       afElemSimulEjecta + j, afMetalSimulEjecta + i);
   }


   pfIMF += iStart;
   pfMasses += iStart;
   pfLogMasses += iStart;


   for (i = 0, j = 0; i < nSize; i++) {
      float fIMFLogWeight = pfIMF[i] * pfMasses[i];
      for (k = 0; k < nElems; j++, k++)
	 afElemSimulEjecta[j] *= fIMFLogWeight;
      afMetalSimulEjecta[i] *= fIMFLogWeight;
   }


   float afElemMassTemp[nElems], fMetalMassTemp;
   for (i = 0; i < nElems; i++)
      afElemMassTemp[i] = 0.0f;
   fMetalMassTemp = 0.0f;


   /* Integration is done over the logarithms of the masses */
   const float fDeltaLogMass = pfLogMasses[1] - pfLogMasses[0];
   const float fWeightStart = ((float)log10(fMassStart) - pfLogMasses[0]) / fDeltaLogMass;
   const float fWeightEnd = (pfLogMasses[nSize - 1] - (float)log10(fMassEnd)) / fDeltaLogMass;

   const float fWeight0 = 0.5f * (1.0f - fWeightStart) * (1.0f - fWeightStart);
   const float fWeight1 = 0.5f * fWeightStart * fWeightStart;
   const float fWeightNm2 = 0.5f * fWeightEnd * fWeightEnd;
   const float fWeightNm1 = 0.5f * (1.0f - fWeightEnd) * (1.0f - fWeightEnd);


   for (i = 1, j = nElems; i < nSize - 1; i++) {
      for (k = 0; k < nElems; j++, k++)
	 afElemMassTemp[k] += afElemSimulEjecta[j];
      fMetalMassTemp += afMetalSimulEjecta[i];
   }

   for (i = 0; i < nElems; i++) {
      afElemMassTemp[i] += fWeight0 * afElemSimulEjecta[i] -
	                   fWeight1 * afElemSimulEjecta[i + nElems] -
	                   fWeightNm2 * afElemSimulEjecta[i + (nSize - 2) * nElems] +
	                   fWeightNm1 * afElemSimulEjecta[i + (nSize - 1) * nElems];
   }
   fMetalMassTemp += fWeight0 * afMetalSimulEjecta[0] -
                     fWeight1 * afMetalSimulEjecta[1] -
                     fWeightNm2 * afMetalSimulEjecta[nSize - 2] +
                     fWeightNm1 * afMetalSimulEjecta[nSize - 1];

   for (i = 0; i < nElems; i++)
      pfElemMass[i] += M_LN10 * fDeltaLogMass * afElemMassTemp[i];
   *pfMetalMass += M_LN10 * fDeltaLogMass * fMetalMassTemp;
}


static inline float stevLifetimeFunction(PKD pkd, STARFIELDS *pStar, const float fMass) {
   int iIdxMass;
   float fDeltaMass;
   stevGetIndex1D(pkd->StelEvolData->afLifetimes_Masses, STEV_LIFETIMES_N_MASS,
		  log10(fMass), &iIdxMass, &fDeltaMass);

   const float fDeltaZ = pStar->Lifetimes.fDeltaZ;

   const float *afTimesLowerZ = pkd->StelEvolData->afLifetimes + pStar->Lifetimes.oZ;
   const float *afTimesUpperZ = afTimesLowerZ + STEV_LIFETIMES_N_MASS;

   float fLogTime = afTimesLowerZ[iIdxMass] * (1.0f - fDeltaZ) * (1.0f - fDeltaMass);
   fLogTime += afTimesLowerZ[iIdxMass + 1] * (1.0f - fDeltaZ) * fDeltaMass;
   fLogTime += afTimesUpperZ[iIdxMass] * fDeltaZ * (1.0f - fDeltaMass);
   fLogTime += afTimesUpperZ[iIdxMass + 1] * fDeltaZ * fDeltaMass;

   return pow(10.0, fLogTime);
}


static inline float stevInverseLifetimeFunction(PKD pkd, STARFIELDS *pStar, const float fTime) {
   const float *afTimesLowerZ = pkd->StelEvolData->afLifetimes + pStar->Lifetimes.oZ;
   const float *afTimesUpperZ = afTimesLowerZ + STEV_LIFETIMES_N_MASS;

   int iIdxTimeLowerZ;
   float fDeltaTimeLowerZ;
   stevGetIndex1DReversed(afTimesLowerZ, STEV_LIFETIMES_N_MASS, log10(fTime),
			  &iIdxTimeLowerZ, &fDeltaTimeLowerZ);

   int iIdxTimeUpperZ;
   float fDeltaTimeUpperZ;
   stevGetIndex1DReversed(afTimesUpperZ, STEV_LIFETIMES_N_MASS, log10(fTime),
			  &iIdxTimeUpperZ, &fDeltaTimeUpperZ);

   const float fDeltaZ = pStar->Lifetimes.fDeltaZ;
   const float *afMasses = pkd->StelEvolData->afLifetimes_Masses;

   float fLogMass = afMasses[iIdxTimeLowerZ] * (1.0f - fDeltaZ) * (1.0f - fDeltaTimeLowerZ);
   fLogMass += afMasses[iIdxTimeLowerZ + 1] * (1.0f - fDeltaZ) * fDeltaTimeLowerZ;
   fLogMass += afMasses[iIdxTimeUpperZ] * fDeltaZ * (1.0f - fDeltaTimeUpperZ);
   fLogMass += afMasses[iIdxTimeUpperZ + 1] * fDeltaZ * fDeltaTimeUpperZ;

   return pow(10.0, fLogMass);
}


#endif	/* STELLAR_EVOLUTION */
