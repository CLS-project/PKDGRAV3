#ifdef STELLAR_EVOLUTION

#include "master.h"


#define STEV_CCSN_N_METALLICITY               5
#define STEV_CCSN_N_MASS                     11
#define STEV_AGB_N_METALLICITY                3
#define STEV_AGB_N_MASS                      23
#define STEV_LIFETIMES_N_METALLICITY          6
#define STEV_LIFETIMES_N_MASS                30

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

   /* Initial mass array for CCSN and AGB tables */
   float afMasses[STEV_INTERP_N_MASS];

   /* Logarithmic spacing of the values in the CCSN/AGB initial mass array */
   float fDeltaLogMass;

   /* Logarithmic weights (= IMF * Masses) for ejecta integration array */
   float afIMFLogWeights[STEV_INTERP_N_MASS];

   /* Core Collapse SNe arrays */
   float afCCSN_Zs[STEV_CCSN_N_METALLICITY];
   float afCCSN_Yields[STEV_CCSN_N_METALLICITY * STEV_INTERP_N_MASS * ELEMENT_COUNT];
   float afCCSN_MetalYield[STEV_CCSN_N_METALLICITY * STEV_INTERP_N_MASS];
   float afCCSN_EjectedMass[STEV_CCSN_N_METALLICITY * STEV_INTERP_N_MASS];

   /* AGB arrays */
   float afAGB_Zs[STEV_AGB_N_METALLICITY];
   float afAGB_Yields[STEV_AGB_N_METALLICITY * STEV_INTERP_N_MASS * ELEMENT_COUNT];
   float afAGB_MetalYield[STEV_AGB_N_METALLICITY * STEV_INTERP_N_MASS];
   float afAGB_EjectedMass[STEV_AGB_N_METALLICITY * STEV_INTERP_N_MASS];

   /* Type Ia SNe arrays */
   float afSNIa_EjectedMass[ELEMENT_COUNT];
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

void msrSetStellarEvolutionParam(MSR);
void msrStellarEvolutionInit(MSR, double);
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

STEV_RAWDATA *stevReadTable(char *);
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

static inline void stevChabrierIMF(const double *restrict pdMasses, const int nSize,
				   const double dMinMass, const double dMaxMass,
				   double *restrict pdIMF) {
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
      if (pdMasses[i] < 1.0)
	 pdIMF[i] = A * exp(-0.5 * pow(log10(pdMasses[i]) - logMc, 2.0) / Sigma2) / pdMasses[i];
      else
	 pdIMF[i] = B * pow(pdMasses[i], Slope);
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


static inline void stevInterpToIMFSampling(STEV_DATA *const Data, STEV_RAWDATA *const CCSN,
					   STEV_RAWDATA *const AGB, const float fCCSNMinMass) {

   /* Yields array in the buffer is ordered as (Metallicities, Masses, Elements). MetalYield
      and EjectedMass arrays, on the other hand, are ordered as (Metallicities, Masses).
      To change this, the arguments of every call to stevRowMajorIndex that set the variable
      idxData must be swapped appropriately. */
   int i, j, k, idxMass, idxTable, idxData;
   float fDeltaMass, fLogMass;

   const int idxCCSNMinMass = stevGetIMFMassIndex(Data->afMasses, STEV_INTERP_N_MASS,
						  fCCSNMinMass, STEV_INTERP_N_MASS - 1);

   for (i = 0; i < STEV_INTERP_N_MASS; i++) {
      fLogMass = log10(Data->afMasses[i]);
      if (i <= idxCCSNMinMass) {
	 stevGetIndex1D(AGB->pfMasses, STEV_AGB_N_MASS, fLogMass, &idxMass, &fDeltaMass);

	 for (j = 0; j < ELEMENT_COUNT; j++) {
	    for (k = 0; k < STEV_CCSN_N_METALLICITY; k++) {
	       idxData = stevRowMajorIndex(k, i, j, STEV_CCSN_N_METALLICITY,
					   STEV_INTERP_N_MASS, ELEMENT_COUNT);
	       if (i < idxCCSNMinMass) {
		  Data->afCCSN_Yields[idxData] = 0.0f;
	       }
	       else {
		  idxTable = stevRowMajorIndex(k, j, 0, STEV_CCSN_N_METALLICITY,
					       ELEMENT_COUNT, STEV_CCSN_N_MASS);
		  Data->afCCSN_Yields[idxData] = CCSN->pfYields[idxTable];
	       }
	    }

	    for (k = 0; k < STEV_AGB_N_METALLICITY; k++) {
	       idxTable = stevRowMajorIndex(k, j, idxMass, STEV_AGB_N_METALLICITY, ELEMENT_COUNT,
					    STEV_AGB_N_MASS);
	       idxData = stevRowMajorIndex(k, i, j, STEV_AGB_N_METALLICITY,
					   STEV_INTERP_N_MASS, ELEMENT_COUNT);

	       /* WARNING: The following formula is correct only when mass is the last index
		  in the tables */
	       Data->afAGB_Yields[idxData] =
		  AGB->pfYields[idxTable] * (1.0 - fDeltaMass) +
		  AGB->pfYields[idxTable + 1] * fDeltaMass;
	    }
	 }

	 for (k = 0; k < STEV_CCSN_N_METALLICITY; k++) {
	    idxData = stevRowMajorIndex(k, i, 0, STEV_CCSN_N_METALLICITY,
					STEV_INTERP_N_MASS, 1);
	    if (i < idxCCSNMinMass) {
	       Data->afCCSN_MetalYield[idxData]  = 0.0f;
	       Data->afCCSN_EjectedMass[idxData] = 0.0f;
	    }
	    else {
	       idxTable = stevRowMajorIndex(k, 0, 0, STEV_CCSN_N_METALLICITY,
					    STEV_CCSN_N_MASS, 1);
	       Data->afCCSN_MetalYield[idxData]  = CCSN->pfMetalYield[idxTable];
	       Data->afCCSN_EjectedMass[idxData] = CCSN->pfEjectedMass[idxTable];
	    }
	 }

	 for (k = 0; k < STEV_AGB_N_METALLICITY; k++) {
	    idxTable = stevRowMajorIndex(k, idxMass, 0, STEV_AGB_N_METALLICITY,
					 STEV_AGB_N_MASS, 1);
	    idxData = stevRowMajorIndex(k, i, 0, STEV_AGB_N_METALLICITY,
					STEV_INTERP_N_MASS, 1);

	    /* WARNING: The following formulas are correct only when mass is the last index
	       in the tables */
	    Data->afAGB_MetalYield[idxData] =
	       AGB->pfMetalYield[idxTable] * (1.0 - fDeltaMass) +
	       AGB->pfMetalYield[idxTable + 1] * fDeltaMass;

	    Data->afAGB_EjectedMass[idxData] =
	       AGB->pfEjectedMass[idxTable] * (1.0 - fDeltaMass) +
	       AGB->pfEjectedMass[idxTable + 1] * fDeltaMass;
	 }
      }
      else {
	 stevGetIndex1D(CCSN->pfMasses, STEV_CCSN_N_MASS, fLogMass, &idxMass, &fDeltaMass);

	 for (j = 0; j < ELEMENT_COUNT; j++) {
	    for (k = 0; k < STEV_CCSN_N_METALLICITY; k++) {
	       idxTable = stevRowMajorIndex(k, j, idxMass, STEV_CCSN_N_METALLICITY,
					    ELEMENT_COUNT, STEV_CCSN_N_MASS);
	       idxData = stevRowMajorIndex(k, i, j, STEV_CCSN_N_METALLICITY,
					   STEV_INTERP_N_MASS, ELEMENT_COUNT);

	       /* WARNING: The following formula is correct only when mass is the last index
		  in the tables */
	       Data->afCCSN_Yields[idxData] =
		  CCSN->pfYields[idxTable] * (1.0 - fDeltaMass) +
		  CCSN->pfYields[idxTable + 1] * fDeltaMass;
	    }

	    for (k = 0; k < STEV_AGB_N_METALLICITY; k++) {
	       idxData = stevRowMajorIndex(k, i, j, STEV_AGB_N_METALLICITY,
					   STEV_INTERP_N_MASS, ELEMENT_COUNT);
	       if (i > idxCCSNMinMass + 1) {
		  Data->afAGB_Yields[idxData] = 0.0f;
	       }
	       else {
		  idxTable = stevRowMajorIndex(k, j, STEV_AGB_N_MASS - 1,
					       STEV_AGB_N_METALLICITY, ELEMENT_COUNT,
					       STEV_AGB_N_MASS);
		  Data->afAGB_Yields[idxData] = AGB->pfYields[idxTable];
	       }
	    }
	 }

	 for (k = 0; k < STEV_CCSN_N_METALLICITY; k++) {
	    idxTable = stevRowMajorIndex(k, idxMass, 0, STEV_CCSN_N_METALLICITY,
					 STEV_CCSN_N_MASS, 1);
	    idxData = stevRowMajorIndex(k, i, 0, STEV_CCSN_N_METALLICITY,
					STEV_INTERP_N_MASS, 1);

	    /* WARNING: The following formulas are correct only when mass is the last index
	       in the tables */
	    Data->afCCSN_MetalYield[idxData] =
	       CCSN->pfMetalYield[idxTable] * (1.0 - fDeltaMass) +
	       CCSN->pfMetalYield[idxTable + 1] * fDeltaMass;

	    Data->afCCSN_EjectedMass[idxData] =
	       CCSN->pfEjectedMass[idxTable] * (1.0 - fDeltaMass) +
	       CCSN->pfEjectedMass[idxTable + 1] * fDeltaMass;
	 }

	 for (k = 0; k < STEV_AGB_N_METALLICITY; k++) {
	    idxData = stevRowMajorIndex(k, i, 0, STEV_AGB_N_METALLICITY,
					STEV_INTERP_N_MASS, 1);
	    if (i > idxCCSNMinMass + 1) {
	       Data->afAGB_MetalYield[idxData]  = 0.0f;
	       Data->afAGB_EjectedMass[idxData] = 0.0f;
	    }
	    else {
	       idxTable = stevRowMajorIndex(k, STEV_AGB_N_MASS - 1, 0,
					    STEV_AGB_N_METALLICITY, STEV_AGB_N_MASS, 1);
	       Data->afAGB_MetalYield[idxData]  = AGB->pfMetalYield[idxTable];
	       Data->afAGB_EjectedMass[idxData] = AGB->pfEjectedMass[idxTable];
	    }
	 }
      }
   }
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
   *pfMetalEjMass = fMetalYield + fMetalAbun * fEjectedMass;
   for (i = 0; i < nElems; i++) {
      pfElemEjMass[i] = pfYields[i] + pfElemAbun[i] * fEjectedMass;
      if (pfElemEjMass[i] < 0.0f) {
	 if (i != ELEMENT_H && i != ELEMENT_He) *pfMetalEjMass -= pfElemEjMass[i];
	 pfElemEjMass[i] = 0.0f;
      }
   }
   if (*pfMetalEjMass < 0.0f) *pfMetalEjMass = 0.0f;

   const float fTotalMass = pfElemEjMass[ELEMENT_H] + pfElemEjMass[ELEMENT_He] +
			    *pfMetalEjMass;
   const float fNormFactor = fEjectedMass / fTotalMass;
   for (i = 0; i < nElems; i++) pfElemEjMass[i] *= fNormFactor;
   *pfMetalEjMass *= fNormFactor;
}


static inline void stevComputeMassToEject(
		       const float *restrict pfYields,
		       const float *restrict pfMetalYield,
		       const float *restrict pfEjectedMass,
		       const float *restrict pfMasses,
		       const float *restrict pfIMFLogWeights,
		       const float fDeltaLogMass,
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

   pfIMFLogWeights += iStart;
   pfMasses += iStart;

   for (i = 0, j = 0; i < nSize; i++) {
      for (k = 0; k < nElems; j++, k++)
	 afElemSimulEjecta[j] *= pfIMFLogWeights[i];
      afMetalSimulEjecta[i] *= pfIMFLogWeights[i];
   }

   /* Integration is done over the logarithms of the masses. Here, the steps for this are
      followed in a rather cumbersome way. It is intended to avoid round-off error as much
      as possible. */
   float afElemMassTemp[nElems], fMetalMassTemp;

   /* Contribution from the first and last values */
   for (i = 0; i < nElems; i++) {
      afElemMassTemp[i] = 0.5f * (afElemSimulEjecta[i] +
				  afElemSimulEjecta[i + (nSize - 1) * nElems]);
   }
   fMetalMassTemp = 0.5f * (afMetalSimulEjecta[0] + afMetalSimulEjecta[nSize - 1]);

   /* Contribution from the middle values */
   for (i = 1, j = nElems; i < nSize - 1; i++) {
      for (k = 0; k < nElems; j++, k++)
	 afElemMassTemp[k] += afElemSimulEjecta[j];
      fMetalMassTemp += afMetalSimulEjecta[i];
   }

   /* Multiply by the logarithm of the spacing */
   for (i = 0; i < nElems; i++)
      afElemMassTemp[i] *= fDeltaLogMass;
   fMetalMassTemp *= fDeltaLogMass;

   /* Correction for initial and final values mismatch */
   const float fWeightStart = log10f(fMassStart / pfMasses[0]);
   const float fWeightEnd = log10f(pfMasses[nSize - 1] / fMassEnd);
   for (i = 0; i < nElems; i++) {
      afElemMassTemp[i] -=
	 0.5f * (fWeightStart * (afElemSimulEjecta[i] +
				 afElemSimulEjecta[i + nElems]) +
		 fWeightEnd * (afElemSimulEjecta[i + (nSize - 2) * nElems] +
			       afElemSimulEjecta[i + (nSize - 1) * nElems]));
   }
   fMetalMassTemp -=
      0.5f * (fWeightStart * (afMetalSimulEjecta[0] +
			      afMetalSimulEjecta[1]) +
	      fWeightEnd * (afMetalSimulEjecta[nSize - 2] +
			    afMetalSimulEjecta[nSize - 1]));

   /* Multiply by natural logarithm of 10 */
   for (i = 0; i < nElems; i++)
      pfElemMass[i] += M_LN10 * afElemMassTemp[i];
   *pfMetalMass += M_LN10 * fMetalMassTemp;
}


static inline float stevLifetimeFunction(PKD pkd, STARFIELDS *pStar, const double dMass) {
   int iIdxMass;
   float fDeltaMass;
   stevGetIndex1D(pkd->StelEvolData->afLifetimes_Masses, STEV_LIFETIMES_N_MASS,
		  log10(dMass), &iIdxMass, &fDeltaMass);

   const float fDeltaZ = pStar->Lifetimes.fDeltaZ;

   const float *afTimesLowerZ = pkd->StelEvolData->afLifetimes + pStar->Lifetimes.oZ;
   const float *afTimesUpperZ = afTimesLowerZ + STEV_LIFETIMES_N_MASS;

   float fLogTime = afTimesLowerZ[iIdxMass] * (1.0f - fDeltaZ) * (1.0f - fDeltaMass) +
                    afTimesLowerZ[iIdxMass + 1] * (1.0f - fDeltaZ) * fDeltaMass +
                    afTimesUpperZ[iIdxMass] * fDeltaZ * (1.0f - fDeltaMass) +
                    afTimesUpperZ[iIdxMass + 1] * fDeltaZ * fDeltaMass;

   return powf(10.0f, fLogTime);
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

   float fLogMass = afMasses[iIdxTimeLowerZ] * (1.0f - fDeltaZ) * (1.0f - fDeltaTimeLowerZ) +
                    afMasses[iIdxTimeLowerZ + 1] * (1.0f - fDeltaZ) * fDeltaTimeLowerZ +
                    afMasses[iIdxTimeUpperZ] * fDeltaZ * (1.0f - fDeltaTimeUpperZ) +
                    afMasses[iIdxTimeUpperZ + 1] * fDeltaZ * fDeltaTimeUpperZ;

   return powf(10.0f, fLogMass);
}


/* Function that estimates the time it will take a star particle to lose all
 * its initial mass by estimating its current ejecta rate. This is then used
 * to compute the time it needs to eject fMinMassFrac of its mass. */
static inline float stevComputeFirstEnrichTime(PKD pkd, STARFIELDS *pStar) {
   const float fMinMassFrac = 1e-3f;

   const int iEnd = pStar->iLastEnrichMassIdx;
   const int iStart = iEnd - 1;
   const float fMassStart = pkd->StelEvolData->afMasses[iStart];
   const float fMassEnd = pkd->StelEvolData->afMasses[iEnd];

   float afTblEjMass[2];
   if (fMassStart > pkd->param.dCCSN_MinMass) {
      stevInterpolateXAxis(pkd->StelEvolData->afCCSN_EjectedMass, STEV_INTERP_N_MASS,
            1, 2, pStar->CCSN.oZ + iStart, pStar->CCSN.fDeltaZ, afTblEjMass);
   }
   else {
      stevInterpolateXAxis(pkd->StelEvolData->afAGB_EjectedMass, STEV_INTERP_N_MASS,
            1, 2, pStar->AGB.oZ + iStart, pStar->AGB.fDeltaZ, afTblEjMass);
   }
   const float *const pfIMFLogWeights = pkd->StelEvolData->afIMFLogWeights + iStart;
   const float fEjMassFrac = 0.5f * M_LN10 * pkd->StelEvolData->fDeltaLogMass *
      (afTblEjMass[0] * pfIMFLogWeights[0] + afTblEjMass[1] * pfIMFLogWeights[1]);

   const float fTimeStart = stevLifetimeFunction(pkd, pStar, fMassEnd);
   const float fTimeEnd = stevLifetimeFunction(pkd, pStar, fMassStart);
   const float fDepletionTime = (fTimeEnd - fTimeStart) / fEjMassFrac;
   float fNextTime = pStar->fLastEnrichTime + fMinMassFrac * fDepletionTime;

   /* Here we make sure fNextTime is beyond the numerical artifacts of the
    * Lifetime and Inverse Lifetime functions */
   float fMass = stevInverseLifetimeFunction(pkd, pStar, fNextTime);
   const float fStep = 1.05f;
   while (fMass > pStar->fLastEnrichMass) {
      fNextTime *= fStep;
      fMass = stevInverseLifetimeFunction(pkd, pStar, fNextTime);
   }

   return pStar->fTimer + fNextTime;
}


/* Given the ejecta rate of a star particle, use it to compute the time it needs to
 * eject fMinMassFrac of its mass. */
static inline float stevComputeNextEnrichTime(const float fTime, const float fStarMass,
                                              const float fEjMass, const float fDt) {
   const float fMinMassFrac = 1e-3f;
   const float fEjMassRateInv = fDt / fEjMass;

   return fTime + fMinMassFrac * fStarMass * fEjMassRateInv;
}


static inline void stevStarParticleInit(PKD pkd, STARFIELDS *pStar) {
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

   pStar->fSNIaOnsetTime = stevLifetimeFunction(pkd, pStar, pkd->param.dSNIa_MaxMass);

   const float fCCSNOnsetTime = stevLifetimeFunction(pkd, pStar, pkd->param.dCCSN_MaxMass);
   if (pStar->fLastEnrichTime < fCCSNOnsetTime) {
      pStar->fLastEnrichTime = fCCSNOnsetTime;
      pStar->fLastEnrichMass = pkd->param.dCCSN_MaxMass;
      pStar->iLastEnrichMassIdx = STEV_INTERP_N_MASS - 1;
   }
   else {
      pStar->fLastEnrichMass = stevInverseLifetimeFunction(pkd, pStar, pStar->fLastEnrichTime);
      if (pStar->fLastEnrichMass < pkd->param.dCCSN_MaxMass) {
	 pStar->iLastEnrichMassIdx =
	    stevGetIMFMassIndex(pkd->StelEvolData->afMasses, STEV_INTERP_N_MASS,
				pStar->fLastEnrichMass, STEV_INTERP_N_MASS - 1) + 1;
      }
      else {
	 pStar->fLastEnrichTime = fCCSNOnsetTime;
	 pStar->fLastEnrichMass = pkd->param.dCCSN_MaxMass;
	 pStar->iLastEnrichMassIdx = STEV_INTERP_N_MASS - 1;
      }
   }

   pStar->fNextEnrichTime = stevComputeFirstEnrichTime(pkd, pStar);
}


#endif	/* STELLAR_EVOLUTION */
