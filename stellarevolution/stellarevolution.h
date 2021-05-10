#ifdef STELLAR_EVOLUTION


#include "master.h"


/* NOTAS:
   - En EAGLE se impone conservacion de momento y energia al distribuir la masa (Schaye 2015).
   - Arrays del buffer que estan en logaritmo: 
      - CCSN y AGB: Metalicidades
      - Lifetimes: Todo
   - Masas se estan interpolando en "semilog", Lifetimes en "loglog".
*/


/* TODOs y DUDAS:
   - Armar archivos hdf5
   - Implementar la inicializacion cuando se empieza de un archivo restart o de una snapshot
   - Analizar pros y contras de hacer la integracion en logaritmos
   - Considerar almacenar valores interpolados del final del ultimo timestep en STARFIELDS
   - Llamada redundante a msrActiveRung(msr,uRung,1) en msrNewTopStepKDK despues de hacer
   el feedback.

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
   /* Pointer to the function that gives the number of SNIa in [fInitialTime, fFinalTime] */
   float (*fcnNumSNIa)(PKD pkd, STARFIELDS *pStar, float fInitialTime, float fFinalTime);

   /* Initial mass array for CCSN and AGB tables */
   float afMasses[STEV_INTERP_N_MASS];

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

void stevChabrierIMF(float *, int, double, double, float *);
float stevExponentialNumSNIa(PKD, STARFIELDS *, float, float);
float stevPowerlawNumSNIa(PKD, STARFIELDS *, float, float);



/*
 * -------
 * INLINES
 * -------
 */

static inline void stevGetIndex1D(const float *restrict pfTable, const int nSize,
				  const float fVal, int *piIdx, float *restrict pfDelta) {
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


static inline int stevGetIMFMassIndex(const float *restrict pfTable, const int nSize,
				      const float fVal, const int iStart) {

   assert(iStart > 0 && iStart < nSize);
   const float *pfTemp;
   for (pfTemp = pfTable + iStart - 1; pfTemp > pfTable && *pfTemp > fVal; pfTemp++)
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
					const float fWeight, float *pfResult) {

   const float *pfLower = pfTable + iOffset * nZ;
   const float *pfUpper = pfLower + nY * nZ;

   const float *pfLowerEnd = pfLower + nSize * nZ;
   while(pfLower < pfLowerEnd)
      *pfResult++ = *pfLower++ * (1.0f - fWeight) + *pfUpper++ * fWeight;
}


static inline void stevComputeAndCorrectEjecta(const float *restrict pfYields,
					       const float fMetalYield,
					       const float fEjectedMass,
					       const int nElems, const float *restrict pfElemAbun,
					       const float fMetalAbun, float *pfElemEjMass,
					       float *pfMetalEjMass) {

   *pfMetalEjMass = fMetalYield + fMetalAbun * fEjectedMass;

   int i;
   for (i = 0; i < nElems; i++) {
      pfElemEjMass[i] = pfYields[i] + pfElemAbun[i] * fEjectedMass;
      if (pfElemEjMass[i] < 0.0f) {
	 if (i != ELEMENT_H && i != ELEMENT_He) *pfMetalEjMass -= pfElemEjMass[i];
	 pfElemEjMass[i] = 0.0f;
      }
   }
   if (*pfMetalEjMass < 0.0f) *pfMetalEjMass = 0.0f;

   float fTotalMass = pfElemEjMass[ELEMENT_H] + pfElemEjMass[ELEMENT_He] + *pfMetalEjMass;
   float fNormFactor = fEjectedMass / fTotalMass;
   for (i = 0; i < nElems; i++) pfElemEjMass[i] *= fNormFactor;
   *pfMetalEjMass *= fNormFactor;
}


static inline void stevComputeMassToEject(const float *restrict pfYields,
					  const float *restrict pfMetalYield,
					  const float *restrict pfEjectedMass,
					  const float *restrict pfMasses,
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
      If this were to be done, the restrict qualifier in the first argument of
      stevComputeAndCorrectEjecta must be removed */
   float afElemEjMass[nSize * nElems];
   float afMetalEjMass[nSize];
   for (i = 0; i < nSize; i++) {
      stevComputeAndCorrectEjecta(afInterpYields + i * nElems, afInterpMetalYield[i],
				  afInterpEjectedMass[i], nElems, pfElemAbun, fMetalAbun,
				  afElemEjMass + i * nElems, afMetalEjMass + i);
   }


   pfIMF += iStart;
   pfMasses += iStart;


   for (i = 0, j = 0; i < nSize; i++) {
      for (k = 0; k < nElems; j++, k++)
	 afElemEjMass[j] *= pfIMF[i];
      afMetalEjMass[i] *= pfIMF[i];
   }

   /* For now, I am integrating in linear space. The pros and cons of switching to log space
     must be analyzed.
     This is being done over [pfMasses[iStart],pfMasses[iEnd]], which encloses [fMassStart,
     fMassEnd]. The excess is substracted later. */
   for (i = 0, j = 0; i < nSize - 1; i++) {
      for (k = 0; k < nElems; j++, k++) {
	 pfElemMass[k] += (afElemEjMass[j + nElems] + afElemEjMass[j]) *
	                  (pfMasses[i + 1] - pfMasses[i]);
      }
      *pfMetalMass += (afMetalEjMass[i + 1] + afMetalEjMass[i]) *
	              (pfMasses[i + 1] - pfMasses[i]);
   }

   /* Low-mass-end correction. The trapezoid formula looks funny because the ejecta value at
      fMassStart is being explicitly replaced by its linear interpolation. Same goes for the 
      high-mass-end correction. */
   float fDeltaM = (fMassStart - pfMasses[0]) / (pfMasses[1] - pfMasses[0]);
   assert(fDeltaM >= 0.0f);
   for (i = 0; i < nElems; i++) {
      pfElemMass[i] -= (afElemEjMass[i] * (2.0f - fDeltaM) +
			afElemEjMass[i + nElems] * fDeltaM) *
	               (fMassStart - pfMasses[0]);
   }
   *pfMetalMass -= (afMetalEjMass[0] * (2.0f - fDeltaM) +
		    afMetalEjMass[1] * fDeltaM) *
	           (fMassStart - pfMasses[0]);

   /* High-mass-end correction */
   fDeltaM = (fMassEnd - pfMasses[nSize - 2]) / (pfMasses[nSize - 1] - pfMasses[nSize - 2]);
   assert(fDeltaM >= 0.0f);
   for (i = 0; i < nElems; i++) {
      j = i + nElems * (nSize - 2);
      pfElemMass[i] -= (afElemEjMass[j] * (1.0f - fDeltaM) +
			afElemEjMass[j + nElems] * (1.0f + fDeltaM)) *
	               (pfMasses[nSize - 1] - fMassEnd);
   }
   *pfMetalMass -= (afMetalEjMass[nSize - 2] * (1.0f - fDeltaM) +
		    afMetalEjMass[nSize - 1] * (1.0f + fDeltaM)) *
	           (pfMasses[nSize - 1] - fMassEnd);

   /* Final 0.5 factor of the trapezoid rule */
   for (i = 0; i < nElems; i++)
      pfElemMass[i] *= 0.5f;
   *pfMetalMass *= 0.5f;
}


static inline float stevLifetimeFunction(PKD pkd, STARFIELDS *pStar, float fMass) {
   int iIdxMass;
   float fDeltaMass;
   stevGetIndex1D(pkd->StelEvolData->afLifetimes_Masses, STEV_LIFETIMES_N_MASS,
		  log10(fMass), &iIdxMass, &fDeltaMass);

   float fDeltaZ = pStar->Lifetimes.fDeltaZ;

   float *afTimesLowerZ = pkd->StelEvolData->afLifetimes + pStar->Lifetimes.oZ;
   float *afTimesUpperZ = afTimesLowerZ + STEV_LIFETIMES_N_MASS;

   float fLogTime = afTimesLowerZ[iIdxMass] * (1.0f - fDeltaZ) * (1.0f - fDeltaMass);
   fLogTime += afTimesLowerZ[iIdxMass + 1] * (1.0f - fDeltaZ) * fDeltaMass;
   fLogTime += afTimesUpperZ[iIdxMass] * fDeltaZ * (1.0f - fDeltaMass);
   fLogTime += afTimesUpperZ[iIdxMass + 1] * fDeltaZ * fDeltaMass;

   return pow(10.0, fLogTime);
}


static inline float stevInverseLifetimeFunction(PKD pkd, STARFIELDS *pStar, float fTime) {
   float *afTimesLowerZ = pkd->StelEvolData->afLifetimes + pStar->Lifetimes.oZ;
   float *afTimesUpperZ = afTimesLowerZ + STEV_LIFETIMES_N_MASS;

   int iIdxTimeLowerZ;
   float fDeltaTimeLowerZ;
   stevGetIndex1D(afTimesLowerZ, STEV_LIFETIMES_N_MASS, log10(fTime),
		  &iIdxTimeLowerZ, &fDeltaTimeLowerZ);

   int iIdxTimeUpperZ;
   float fDeltaTimeUpperZ;
   stevGetIndex1D(afTimesUpperZ, STEV_LIFETIMES_N_MASS, log10(fTime),
		  &iIdxTimeUpperZ, &fDeltaTimeUpperZ);

   float fDeltaZ = pStar->Lifetimes.fDeltaZ;
   float *afMasses = pkd->StelEvolData->afLifetimes_Masses;

   float fLogMass = afMasses[iIdxTimeLowerZ] * (1.0f - fDeltaZ) * (1.0f - fDeltaTimeLowerZ);
   fLogMass += afMasses[iIdxTimeLowerZ + 1] * (1.0f - fDeltaZ) * fDeltaTimeLowerZ;
   fLogMass += afMasses[iIdxTimeUpperZ] * fDeltaZ * (1.0f - fDeltaTimeUpperZ);
   fLogMass += afMasses[iIdxTimeUpperZ + 1] * fDeltaZ * fDeltaTimeUpperZ;

   return pow(10.0, fLogMass);
}


#endif	/* STELLAR_EVOLUTION */
