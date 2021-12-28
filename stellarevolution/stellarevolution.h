#ifdef STELLAR_EVOLUTION

#include "pst.h"

#ifdef __cplusplus
extern "C" {
#endif

#define STEV_CCSN_N_METALLICITY               5
#define STEV_CCSN_N_MASS                     11
#define STEV_AGB_N_METALLICITY                3
#define STEV_AGB_N_MASS                      23
#define STEV_LIFETIME_N_METALLICITY           6
#define STEV_LIFETIME_N_MASS                 30

#define STEV_INTERP_N_MASS                  200

#define STEV_MIN_LOG_METALLICITY            -20


/*
 * ---------------------
 * STRUCTURE DEFINITIONS
 * ---------------------
 */

typedef struct StellarEvolutionData {
    /* Pointer to function that gives the number of SNIa per Msol in [fInitialTime, fFinalTime] */
    float (*fcnNumSNIa)(SMF *smf, STARFIELDS *pStar, float fInitialTime, float fFinalTime);

    /* Initial mass array for CCSN and AGB tables */
    float afInitialMass[STEV_INTERP_N_MASS];

    /* Logarithmic spacing of the values in the CCSN/AGB initial mass array */
    float fDeltaLogMass;

    /* Logarithmic weights (= IMF * Masses) for ejecta integration array */
    float afIMFLogWeight[STEV_INTERP_N_MASS];

    /* Core Collapse SNe arrays */
    float afCCSNMetallicity[STEV_CCSN_N_METALLICITY];
    float afCCSNYield[STEV_CCSN_N_METALLICITY * STEV_INTERP_N_MASS * ELEMENT_COUNT];
    float afCCSNMetalYield[STEV_CCSN_N_METALLICITY * STEV_INTERP_N_MASS];
    float afCCSNEjectedMass[STEV_CCSN_N_METALLICITY * STEV_INTERP_N_MASS];

    /* AGB arrays */
    float afAGBMetallicity[STEV_AGB_N_METALLICITY];
    float afAGBYield[STEV_AGB_N_METALLICITY * STEV_INTERP_N_MASS * ELEMENT_COUNT];
    float afAGBMetalYield[STEV_AGB_N_METALLICITY * STEV_INTERP_N_MASS];
    float afAGBEjectedMass[STEV_AGB_N_METALLICITY * STEV_INTERP_N_MASS];

    /* Type Ia SNe arrays */
    float afSNIaEjectedMass[ELEMENT_COUNT];
    float fSNIaEjectedMetalMass;

    /* Stellar lifetime arrays */
    float afLifetimeMetallicity[STEV_LIFETIME_N_METALLICITY];
    float afLifetimeInitialMass[STEV_LIFETIME_N_MASS];
    float afLifetime[STEV_LIFETIME_N_METALLICITY * STEV_LIFETIME_N_MASS];
} STEV_DATA;


typedef struct StellarEvolutionRawData {
    int nZs;
    int nElems;
    int nMasses;
    float *pfMetallicity;
    float *pfInitialMass;
    float *pfYield;
    float *pfMetalYield;
    float *pfEjectedMass;
    float *pfLifetime;
} STEV_RAWDATA;


struct inStellarEvolution {
    char achSNIaDTDType[32];
    int bChemEnrich;
    double dTime;
    double dSNIaMaxMass;
    double dCCSNMinMass;
    double dCCSNMaxMass;
    STEV_DATA StelEvolData;
};



/*
 * --------------
 * MAIN FUNCTIONS
 * --------------
 */

int pstStellarEvolutionInit(PST, void *, int, void *, int);
int pkdStellarEvolutionInit(PKD, struct inStellarEvolution *);

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
STEV_RAWDATA *stevReadLifetimeTable(char *);
void stevFreeTable(STEV_RAWDATA *);
void stevFreeSNIaTable(STEV_RAWDATA *);
void stevFreeLifetimeTable(STEV_RAWDATA *);

float stevExponentialNumSNIa(SMF *, STARFIELDS *, float, float);
float stevPowerlawNumSNIa(SMF *, STARFIELDS *, float, float);



/*
 * -------
 * INLINES
 * -------
 */

static inline void stevChabrierIMF(const double *restrict pdMass, const int nSize,
                                   const double dMinMass, const double dMaxMass,
                                   double *restrict pdIMF) {
    const double dMc = 0.079;
    const double dSigma = 0.69;
    const double dSlope = -2.3;

    const double dLogMc = log10(dMc);
    const double dSigmaSq = dSigma * dSigma;
    const double dConstFac = exp(-0.5 * dLogMc * dLogMc / dSigmaSq);

    const double dLogNormInt = sqrt(0.5 * M_PI) * exp(0.5 * dSigmaSq * M_LN10 * M_LN10) *
                               M_LN10 * dSigma * dMc *
                               (erf((-dLogMc / dSigma - M_LN10 * dSigma) / sqrt(2.0)) -
                                erf(((log10(dMinMass) - dLogMc) / dSigma - M_LN10 * dSigma) /
                                    sqrt(2.0)));
    const double dPowerLawInt = dConstFac * (pow(dMaxMass, dSlope + 2.0) - 1.0) / (dSlope + 2.0);

    const double dLogNormFac = 1.0 / (dLogNormInt + dPowerLawInt);
    const double dPowerLawFac = dLogNormFac * dConstFac;
    for (int i = 0; i < nSize; i++) {
        if (pdMass[i] < 1.0)
            pdIMF[i] = dLogNormFac * exp(-0.5 * pow(log10(pdMass[i]) - dLogMc, 2.0) / dSigmaSq) /
                       pdMass[i];
        else
            pdIMF[i] = dPowerLawFac * pow(pdMass[i], dSlope);
    }
}


static inline void stevGetIndex1D(const float *restrict pfTable, const int nSize,
                                  const float fVal, int *restrict piIdx,
                                  float *restrict pfDelta) {

    const float fEpsilon = 1e-4f;

    if (fVal < pfTable[0] + fEpsilon) {
        /* We are below the first element */
        *piIdx = 0;
        *pfDelta = 0.0f;
    }
    else if (fVal > pfTable[nSize - 1] - fEpsilon) {
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

    const float fEpsilon = 1e-4f;

    if (fVal > pfTable[0] - fEpsilon) {
        /* We are above the first element */
        *piIdx = 0;
        *pfDelta = 0.0f;
    }
    else if (fVal < pfTable[nSize - 1] + fEpsilon) {
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

    /* Yield arrays in Data are ordered as (Metallicity, Mass, Element). MetalYield
       and EjectedMass arrays, on the other hand, are ordered as (Metallicity, Mass).
       To change this, the arguments of every call to stevRowMajorIndex that set the
       variable iData must be swapped appropriately. */
    int i, j, k, iMass, iTable, iData;
    float fDeltaMass, fLogMass;

    const int iCCSNMinMass = stevGetIMFMassIndex(Data->afInitialMass, STEV_INTERP_N_MASS,
                             fCCSNMinMass, STEV_INTERP_N_MASS - 1);

    for (i = 0; i < STEV_INTERP_N_MASS; i++) {
        fLogMass = log10(Data->afInitialMass[i]);
        if (i <= iCCSNMinMass) {
            stevGetIndex1D(AGB->pfInitialMass, STEV_AGB_N_MASS, fLogMass, &iMass, &fDeltaMass);

            for (j = 0; j < ELEMENT_COUNT; j++) {
                for (k = 0; k < STEV_CCSN_N_METALLICITY; k++) {
                    iData = stevRowMajorIndex(k, i, j, STEV_CCSN_N_METALLICITY,
                                              STEV_INTERP_N_MASS, ELEMENT_COUNT);
                    if (i < iCCSNMinMass) {
                        Data->afCCSNYield[iData] = 0.0f;
                    }
                    else {
                        iTable = stevRowMajorIndex(k, j, 0, STEV_CCSN_N_METALLICITY,
                                                   ELEMENT_COUNT, STEV_CCSN_N_MASS);
                        Data->afCCSNYield[iData] = CCSN->pfYield[iTable];
                    }
                }

                for (k = 0; k < STEV_AGB_N_METALLICITY; k++) {
                    iTable = stevRowMajorIndex(k, j, iMass, STEV_AGB_N_METALLICITY,
                                               ELEMENT_COUNT, STEV_AGB_N_MASS);
                    iData = stevRowMajorIndex(k, i, j, STEV_AGB_N_METALLICITY,
                                              STEV_INTERP_N_MASS, ELEMENT_COUNT);

                    /* WARNING: The following formula is correct only when mass is the last index
                       in the tables */
                    Data->afAGBYield[iData] = AGB->pfYield[iTable] * (1.0 - fDeltaMass) +
                                              AGB->pfYield[iTable + 1] * fDeltaMass;
                }
            }

            for (k = 0; k < STEV_CCSN_N_METALLICITY; k++) {
                iData = stevRowMajorIndex(k, i, 0, STEV_CCSN_N_METALLICITY,
                                          STEV_INTERP_N_MASS, 1);
                if (i < iCCSNMinMass) {
                    Data->afCCSNMetalYield[iData]  = 0.0f;
                    Data->afCCSNEjectedMass[iData] = 0.0f;
                }
                else {
                    iTable = stevRowMajorIndex(k, 0, 0, STEV_CCSN_N_METALLICITY,
                                               STEV_CCSN_N_MASS, 1);
                    Data->afCCSNMetalYield[iData]  = CCSN->pfMetalYield[iTable];
                    Data->afCCSNEjectedMass[iData] = CCSN->pfEjectedMass[iTable];
                }
            }

            for (k = 0; k < STEV_AGB_N_METALLICITY; k++) {
                iTable = stevRowMajorIndex(k, iMass, 0, STEV_AGB_N_METALLICITY,
                                           STEV_AGB_N_MASS, 1);
                iData = stevRowMajorIndex(k, i, 0, STEV_AGB_N_METALLICITY,
                                          STEV_INTERP_N_MASS, 1);

                /* WARNING: The following formulas are correct only when mass is the last index
                   in the tables */
                Data->afAGBMetalYield[iData] = AGB->pfMetalYield[iTable] * (1.0 - fDeltaMass) +
                                               AGB->pfMetalYield[iTable + 1] * fDeltaMass;

                Data->afAGBEjectedMass[iData] = AGB->pfEjectedMass[iTable] * (1.0 - fDeltaMass) +
                                                AGB->pfEjectedMass[iTable + 1] * fDeltaMass;
            }
        }
        else {
            stevGetIndex1D(CCSN->pfInitialMass, STEV_CCSN_N_MASS, fLogMass, &iMass, &fDeltaMass);

            for (j = 0; j < ELEMENT_COUNT; j++) {
                for (k = 0; k < STEV_CCSN_N_METALLICITY; k++) {
                    iTable = stevRowMajorIndex(k, j, iMass, STEV_CCSN_N_METALLICITY,
                                               ELEMENT_COUNT, STEV_CCSN_N_MASS);
                    iData = stevRowMajorIndex(k, i, j, STEV_CCSN_N_METALLICITY,
                                              STEV_INTERP_N_MASS, ELEMENT_COUNT);

                    /* WARNING: The following formula is correct only when mass is the last index
                       in the tables */
                    Data->afCCSNYield[iData] = CCSN->pfYield[iTable] * (1.0 - fDeltaMass) +
                                               CCSN->pfYield[iTable + 1] * fDeltaMass;
                }

                for (k = 0; k < STEV_AGB_N_METALLICITY; k++) {
                    iData = stevRowMajorIndex(k, i, j, STEV_AGB_N_METALLICITY,
                                              STEV_INTERP_N_MASS, ELEMENT_COUNT);
                    if (i > iCCSNMinMass + 1) {
                        Data->afAGBYield[iData] = 0.0f;
                    }
                    else {
                        iTable = stevRowMajorIndex(k, j, STEV_AGB_N_MASS - 1,
                                                   STEV_AGB_N_METALLICITY, ELEMENT_COUNT, STEV_AGB_N_MASS);
                        Data->afAGBYield[iData] = AGB->pfYield[iTable];
                    }
                }
            }

            for (k = 0; k < STEV_CCSN_N_METALLICITY; k++) {
                iTable = stevRowMajorIndex(k, iMass, 0, STEV_CCSN_N_METALLICITY,
                                           STEV_CCSN_N_MASS, 1);
                iData = stevRowMajorIndex(k, i, 0, STEV_CCSN_N_METALLICITY,
                                          STEV_INTERP_N_MASS, 1);

                /* WARNING: The following formulas are correct only when mass is the last index
                   in the tables */
                Data->afCCSNMetalYield[iData] =
                    CCSN->pfMetalYield[iTable] * (1.0 - fDeltaMass) +
                    CCSN->pfMetalYield[iTable + 1] * fDeltaMass;

                Data->afCCSNEjectedMass[iData] =
                    CCSN->pfEjectedMass[iTable] * (1.0 - fDeltaMass) +
                    CCSN->pfEjectedMass[iTable + 1] * fDeltaMass;
            }

            for (k = 0; k < STEV_AGB_N_METALLICITY; k++) {
                iData = stevRowMajorIndex(k, i, 0, STEV_AGB_N_METALLICITY,
                                          STEV_INTERP_N_MASS, 1);
                if (i > iCCSNMinMass + 1) {
                    Data->afAGBMetalYield[iData]  = 0.0f;
                    Data->afAGBEjectedMass[iData] = 0.0f;
                }
                else {
                    iTable = stevRowMajorIndex(k, STEV_AGB_N_MASS - 1, 0,
                                               STEV_AGB_N_METALLICITY, STEV_AGB_N_MASS, 1);
                    Data->afAGBMetalYield[iData]  = AGB->pfMetalYield[iTable];
                    Data->afAGBEjectedMass[iData] = AGB->pfEjectedMass[iTable];
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
    while (pfLower < pfLowerEnd)
        *pfResult++ = *pfLower++ * (1.0f - fWeight) + *pfUpper++ * fWeight;
}


static inline void stevComputeAndCorrectSimulEjecta(
    const float *restrict pfYield,
    const float fMetalYield,
    const float fEjectedMass,
    const int nElems, const float *restrict pfElemAbun,
    const float fMetalAbun, float *restrict pfElemEjMass,
    float *restrict pfMetalEjMass) {

    int i;
    *pfMetalEjMass = fMetalYield + fMetalAbun * fEjectedMass;
    for (i = 0; i < nElems; i++) {
        pfElemEjMass[i] = pfYield[i] + pfElemAbun[i] * fEjectedMass;
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
    const float *restrict pfYield,
    const float *restrict pfMetalYield,
    const float *restrict pfEjectedMass,
    const float *restrict pfInitialMass,
    const float *restrict pfIMFLogWeight,
    const float fDeltaLogMass,
    const int nZs, const int nMasses, const int nElems,
    const float *restrict pfElemAbun, const float fMetalAbun,
    const int iStart, const int iEnd,
    const float fMassStart, const float fMassEnd,
    const int iOffset, const float fWeight,
    float *pfElemMass, float *pfMetalMass) {

    int i, j, k;
    const int nSize = iEnd - iStart + 1;
    float afInterpYield[nSize * nElems];
    float afInterpMetalYield[nSize];
    float afInterpEjectedMass[nSize];

    stevInterpolateXAxis(pfYield, nMasses, nElems, nSize, iOffset + iStart, fWeight,
                         afInterpYield);
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
        stevComputeAndCorrectSimulEjecta(afInterpYield + j, afInterpMetalYield[i],
                                         afInterpEjectedMass[i], nElems, pfElemAbun, fMetalAbun,
                                         afElemSimulEjecta + j, afMetalSimulEjecta + i);
    }

    pfIMFLogWeight += iStart;
    pfInitialMass += iStart;

    for (i = 0, j = 0; i < nSize; i++) {
        for (k = 0; k < nElems; j++, k++)
            afElemSimulEjecta[j] *= pfIMFLogWeight[i];
        afMetalSimulEjecta[i] *= pfIMFLogWeight[i];
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
    const float fWeightStart = log10f(fMassStart / pfInitialMass[0]);
    const float fWeightEnd = log10f(pfInitialMass[nSize - 1] / fMassEnd);
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
    int iMass;
    float fDeltaMass;
    stevGetIndex1D(pkd->StelEvolData->afLifetimeInitialMass, STEV_LIFETIME_N_MASS,
                   log10(dMass), &iMass, &fDeltaMass);

    const float fDeltaZ = pStar->Lifetime.fDeltaZ;

    const float *afTimeLowerZ = pkd->StelEvolData->afLifetime + pStar->Lifetime.oZ;
    const float *afTimeUpperZ = afTimeLowerZ + STEV_LIFETIME_N_MASS;

    float fLogTime = afTimeLowerZ[iMass] * (1.0f - fDeltaZ) * (1.0f - fDeltaMass) +
                     afTimeLowerZ[iMass + 1] * (1.0f - fDeltaZ) * fDeltaMass +
                     afTimeUpperZ[iMass] * fDeltaZ * (1.0f - fDeltaMass) +
                     afTimeUpperZ[iMass + 1] * fDeltaZ * fDeltaMass;

    return powf(10.0f, fLogTime);
}


static inline float stevInverseLifetimeFunction(PKD pkd, STARFIELDS *pStar, const float fTime) {
    const float *afTimeLowerZ = pkd->StelEvolData->afLifetime + pStar->Lifetime.oZ;
    const float *afTimeUpperZ = afTimeLowerZ + STEV_LIFETIME_N_MASS;

    int iTimeLowerZ;
    float fDeltaTimeLowerZ;
    stevGetIndex1DReversed(afTimeLowerZ, STEV_LIFETIME_N_MASS, log10(fTime),
                           &iTimeLowerZ, &fDeltaTimeLowerZ);

    int iTimeUpperZ;
    float fDeltaTimeUpperZ;
    stevGetIndex1DReversed(afTimeUpperZ, STEV_LIFETIME_N_MASS, log10(fTime),
                           &iTimeUpperZ, &fDeltaTimeUpperZ);

    const float fDeltaZ = pStar->Lifetime.fDeltaZ;
    const float *afMass = pkd->StelEvolData->afLifetimeInitialMass;

    float fLogMass = afMass[iTimeLowerZ] * (1.0f - fDeltaZ) * (1.0f - fDeltaTimeLowerZ) +
                     afMass[iTimeLowerZ + 1] * (1.0f - fDeltaZ) * fDeltaTimeLowerZ +
                     afMass[iTimeUpperZ] * fDeltaZ * (1.0f - fDeltaTimeUpperZ) +
                     afMass[iTimeUpperZ + 1] * fDeltaZ * fDeltaTimeUpperZ;

    return powf(10.0f, fLogMass);
}


/* Function that estimates the time it will take a star particle to lose all
 * its initial mass by estimating its current ejecta rate. This is then used
 * to compute the time it needs to eject fMinMassFrac of its mass. */
static inline float stevComputeFirstEnrichTime(PKD pkd, STARFIELDS *pStar,
        const double dCCSNMinMass) {
    const float fMinMassFrac = 1e-3f;

    const int iEnd = pStar->iLastEnrichMass;
    const int iStart = iEnd - 1;
    const float fMassStart = pkd->StelEvolData->afInitialMass[iStart];
    const float fMassEnd = pkd->StelEvolData->afInitialMass[iEnd];

    float afTblEjMass[2];
    if (fMassStart > dCCSNMinMass) {
        stevInterpolateXAxis(pkd->StelEvolData->afCCSNEjectedMass, STEV_INTERP_N_MASS,
                             1, 2, pStar->CCSN.oZ + iStart, pStar->CCSN.fDeltaZ, afTblEjMass);
    }
    else {
        stevInterpolateXAxis(pkd->StelEvolData->afAGBEjectedMass, STEV_INTERP_N_MASS,
                             1, 2, pStar->AGB.oZ + iStart, pStar->AGB.fDeltaZ, afTblEjMass);
    }
    const float *const pfIMFLogWeight = pkd->StelEvolData->afIMFLogWeight + iStart;
    const float fEjMassFrac = 0.5f * M_LN10 * pkd->StelEvolData->fDeltaLogMass *
                              (afTblEjMass[0] * pfIMFLogWeight[0] +
                               afTblEjMass[1] * pfIMFLogWeight[1]);

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
        const float fEjectedMass, const float fTimestep) {
    const float fMinMassFrac = 1e-3f;
    const float fEjMassRateInv = fTimestep / fEjectedMass;

    return fTime + fMinMassFrac * fStarMass * fEjMassRateInv;
}


static inline void stevStarParticleInit(PKD pkd, STARFIELDS *pStar, const double dSNIaMaxMass,
                                        const double dCCSNMinMass, const double dCCSNMaxMass) {
    int iZ;
    float fLogZ;
    if (pStar->fMetalAbun > 0.0f)
        fLogZ = log10f(pStar->fMetalAbun);
    else
        fLogZ = STEV_MIN_LOG_METALLICITY;

    stevGetIndex1D(pkd->StelEvolData->afCCSNMetallicity, STEV_CCSN_N_METALLICITY,
                   fLogZ, &iZ, &pStar->CCSN.fDeltaZ);
    pStar->CCSN.oZ = iZ * STEV_INTERP_N_MASS;
    stevGetIndex1D(pkd->StelEvolData->afAGBMetallicity, STEV_AGB_N_METALLICITY,
                   fLogZ, &iZ, &pStar->AGB.fDeltaZ);
    pStar->AGB.oZ = iZ * STEV_INTERP_N_MASS;
    stevGetIndex1D(pkd->StelEvolData->afLifetimeMetallicity, STEV_LIFETIME_N_METALLICITY,
                   fLogZ, &iZ, &pStar->Lifetime.fDeltaZ);
    pStar->Lifetime.oZ = iZ * STEV_LIFETIME_N_MASS;

    pStar->fSNIaOnsetTime = stevLifetimeFunction(pkd, pStar, dSNIaMaxMass);

    const float fCCSNOnsetTime = stevLifetimeFunction(pkd, pStar, dCCSNMaxMass);
    if (pStar->fLastEnrichTime < fCCSNOnsetTime) {
        pStar->fLastEnrichTime = fCCSNOnsetTime;
        pStar->fLastEnrichMass = dCCSNMaxMass;
        pStar->iLastEnrichMass = STEV_INTERP_N_MASS - 1;
    }
    else {
        pStar->fLastEnrichMass = stevInverseLifetimeFunction(pkd, pStar, pStar->fLastEnrichTime);
        if (pStar->fLastEnrichMass < dCCSNMaxMass) {
            pStar->iLastEnrichMass =
                stevGetIMFMassIndex(pkd->StelEvolData->afInitialMass, STEV_INTERP_N_MASS,
                                    pStar->fLastEnrichMass, STEV_INTERP_N_MASS - 1) + 1;
        }
        else {
            pStar->fLastEnrichTime = fCCSNOnsetTime;
            pStar->fLastEnrichMass = dCCSNMaxMass;
            pStar->iLastEnrichMass = STEV_INTERP_N_MASS - 1;
        }
    }

    pStar->fNextEnrichTime = stevComputeFirstEnrichTime(pkd, pStar, dCCSNMinMass);
}

#ifdef __cplusplus
}
#endif

#endif  /* STELLAR_EVOLUTION */
