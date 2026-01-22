
#include <hdf5.h>
#include <boost/container/small_vector.hpp>

#include "stellarevolution.h"
#include "imf.h"
#include "hydro/hydro.h"
#include "smooth/smooth.h"

struct stevPack {
    blitz::TinyVector<double, 3> position;
    float fDensity;
    uint8_t iClass;
};

struct stevFlush {
    blitz::TinyVector<float, 3> ReceivedMom;
    float fReceivedMass;
    float fReceivedE;
    blitz::TinyVector<float, ELEMENT_COUNT> ElemMass;
    float fMetalMass;
};

class SmoothHelperChemEnrich : public SmoothHelperPacked<stevPack, stevFlush> {
    double dTime, dWindSpecificEkin, dCCSNMinMass, dSNIaNorm, dSNIaScale;
    void Smooth(particleStore::ParticleReference& p, float fBall, std::span<NN> nnList) final;
    virtual void pack(void* dst, const void* src) override final;
    virtual void unpack(void* dst, const void* src, const void* key) override final;
    virtual void init(void* dst) override final;
    virtual void flush(void* dst, const void* src) override final;
    virtual void combine(void* dst, const void* src, const void* key) override final;
    virtual void init_shared(void* dst) override final { init_shared_gas_decorator(pkd, dst); }

  public:
    SmoothHelperChemEnrich(pkdContext* pkd, double dTime, double dWindSpecificEkin, double dCCSNMinMass, double dSNIaNorm, double dSNIaScale)
        : SmoothHelperPacked<stevPack, stevFlush>(pkd, true), dTime(dTime), dWindSpecificEkin(dWindSpecificEkin), dCCSNMinMass(dCCSNMinMass),
          dSNIaNorm(dSNIaNorm), dSNIaScale(dSNIaScale) {
        nCacheSize = pkd->decorators.DecoratorSize() + pkd->particles.ParticleSize();
    }
    auto range_neighbour(treeStore::NodePointer node) -> std::pair<int, int> override { return smooth::policy::OnlyGas::range(node); }
    auto accept_neighbour(particleStore::ParticleReference& p) -> bool override { return smooth::policy::OnlyGas::select(p); }
};

struct SmoothChemEnrich : smooth::policy::OnlyGas {
    double dTime;
    inline auto active(particleStore::ParticleReference& p) -> bool {
        if (p.is_star()) {
            auto& star = p.star();
            return dTime > star.fNextEnrichTime;
        }
        return false;
    }
};

auto ServiceChemEnrich::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto in = reinterpret_cast<const input*>(vin.data());
    auto out = reinterpret_cast<output*>(vout.data());
    auto pkd = pst->lcl()->pkd;
    auto helper = std::make_shared<SmoothHelperChemEnrich>(pkd, in->dTime, in->dWindSpecificEkin, in->dCCSNMinMass, in->dSNIaNorm, in->dSNIaScale);
    smContext smx(pkd, helper, in->nSmooth, in->bPeriodic);
    *out = smx.Smooth<SmoothChemEnrich>({
        .dTime = in->dTime,
    });
    return sizeof(output);
}

//**********************************************************************//

/*
 * -------
 * INLINES
 * -------
 */
namespace stev {
static inline void GetIndex1D(const float* restrict pfTable, const int nSize, const float fVal, int* restrict piIdx, float* restrict pfDelta) {

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
        for (i = 1; (i < nSize - 1) && (fVal > pfTable[i]); ++i);

        *piIdx = --i;
        *pfDelta = (fVal - pfTable[i]) / (pfTable[i + 1] - pfTable[i]);
    }
}

static inline void GetIndex1DReversed(const float* restrict pfTable, const int nSize, const float fVal, int* restrict piIdx, float* restrict pfDelta) {

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
        for (i = nSize - 2; i > 0 && fVal > pfTable[i]; i--);

        *piIdx = i;
        *pfDelta = (pfTable[i] - fVal) / (pfTable[i] - pfTable[i + 1]);
    }
}

static inline int RowMajorIndex(const int iX, const int iY, const int iZ, const int nX, const int nY, const int nZ) {

    return iX * nY * nZ + iY * nZ + iZ;
}

static inline int GetIMFMassIndex(const float* restrict pfTable, const int nSize, const float fVal, const int iStart) {

    assert(iStart > 0 && iStart < nSize);
    const float* pfTemp;
    for (pfTemp = pfTable + iStart - 1; pfTemp > pfTable && *pfTemp > fVal; pfTemp--);

    return pfTemp - pfTable;
}

void InterpToIMFSampling(STEV_DATA* const Data, STEV_RAWDATA* const CCSN, STEV_RAWDATA* const AGB, const float fCCSNMinMass) {

    /* Yield arrays in Data are ordered as (Metallicity, Mass, Element). MetalYield
       and EjectedMass arrays, on the other hand, are ordered as (Metallicity, Mass).
       To change this, the arguments of every call to stev::RowMajorIndex that set the
       variable iData must be swapped appropriately. */
    int iMass, iTable, iData;
    float fDeltaMass;

    const int iCCSNMinMass = GetIMFMassIndex(Data->afInitialMass, STEV_INTERP_N_MASS, fCCSNMinMass, STEV_INTERP_N_MASS - 1);

    for (auto i = 0; i < STEV_INTERP_N_MASS; ++i) {
        const float fLogMass = log10(Data->afInitialMass[i]);
        if (i <= iCCSNMinMass) {
            GetIndex1D(AGB->pfInitialMass, STEV_AGB_N_MASS, fLogMass, &iMass, &fDeltaMass);

            for (auto j = 0; j < ELEMENT_COUNT; ++j) {
                for (auto k = 0; k < STEV_CCSN_N_METALLICITY; ++k) {
                    iData = RowMajorIndex(k, i, j, STEV_CCSN_N_METALLICITY, STEV_INTERP_N_MASS, ELEMENT_COUNT);
                    if (i < iCCSNMinMass) {
                        Data->afCCSNYield[iData] = 0.0f;
                    }
                    else {
                        iTable = RowMajorIndex(k, j, 0, STEV_CCSN_N_METALLICITY, ELEMENT_COUNT, STEV_CCSN_N_MASS);
                        Data->afCCSNYield[iData] = CCSN->pfYield[iTable];
                    }
                }

                for (auto k = 0; k < STEV_AGB_N_METALLICITY; ++k) {
                    iTable = RowMajorIndex(k, j, iMass, STEV_AGB_N_METALLICITY, ELEMENT_COUNT, STEV_AGB_N_MASS);
                    iData = RowMajorIndex(k, i, j, STEV_AGB_N_METALLICITY, STEV_INTERP_N_MASS, ELEMENT_COUNT);

                    /* WARNING: The following formula is correct only when mass is the last index
                       in the tables */
                    Data->afAGBYield[iData] = AGB->pfYield[iTable] * (1.0 - fDeltaMass) + AGB->pfYield[iTable + 1] * fDeltaMass;
                }
            }

            for (auto k = 0; k < STEV_CCSN_N_METALLICITY; ++k) {
                iData = RowMajorIndex(k, i, 0, STEV_CCSN_N_METALLICITY, STEV_INTERP_N_MASS, 1);
                if (i < iCCSNMinMass) {
                    Data->afCCSNMetalYield[iData] = 0.0f;
                    Data->afCCSNEjectedMass[iData] = 0.0f;
                }
                else {
                    iTable = RowMajorIndex(k, 0, 0, STEV_CCSN_N_METALLICITY, STEV_CCSN_N_MASS, 1);
                    Data->afCCSNMetalYield[iData] = CCSN->pfMetalYield[iTable];
                    Data->afCCSNEjectedMass[iData] = CCSN->pfEjectedMass[iTable];
                }
            }

            for (auto k = 0; k < STEV_AGB_N_METALLICITY; ++k) {
                iTable = RowMajorIndex(k, iMass, 0, STEV_AGB_N_METALLICITY, STEV_AGB_N_MASS, 1);
                iData = RowMajorIndex(k, i, 0, STEV_AGB_N_METALLICITY, STEV_INTERP_N_MASS, 1);

                /* WARNING: The following formulas are correct only when mass is the last index
                   in the tables */
                Data->afAGBMetalYield[iData] = AGB->pfMetalYield[iTable] * (1.0 - fDeltaMass) + AGB->pfMetalYield[iTable + 1] * fDeltaMass;

                Data->afAGBEjectedMass[iData] = AGB->pfEjectedMass[iTable] * (1.0 - fDeltaMass) + AGB->pfEjectedMass[iTable + 1] * fDeltaMass;
            }
        }
        else {
            GetIndex1D(CCSN->pfInitialMass, STEV_CCSN_N_MASS, fLogMass, &iMass, &fDeltaMass);

            for (auto j = 0; j < ELEMENT_COUNT; ++j) {
                for (auto k = 0; k < STEV_CCSN_N_METALLICITY; ++k) {
                    iTable = RowMajorIndex(k, j, iMass, STEV_CCSN_N_METALLICITY, ELEMENT_COUNT, STEV_CCSN_N_MASS);
                    iData = RowMajorIndex(k, i, j, STEV_CCSN_N_METALLICITY, STEV_INTERP_N_MASS, ELEMENT_COUNT);

                    /* WARNING: The following formula is correct only when mass is the last index
                       in the tables */
                    Data->afCCSNYield[iData] = CCSN->pfYield[iTable] * (1.0 - fDeltaMass) + CCSN->pfYield[iTable + 1] * fDeltaMass;
                }

                for (auto k = 0; k < STEV_AGB_N_METALLICITY; ++k) {
                    iData = RowMajorIndex(k, i, j, STEV_AGB_N_METALLICITY, STEV_INTERP_N_MASS, ELEMENT_COUNT);
                    if (i > iCCSNMinMass + 1) {
                        Data->afAGBYield[iData] = 0.0f;
                    }
                    else {
                        iTable = RowMajorIndex(k, j, STEV_AGB_N_MASS - 1, STEV_AGB_N_METALLICITY, ELEMENT_COUNT, STEV_AGB_N_MASS);
                        Data->afAGBYield[iData] = AGB->pfYield[iTable];
                    }
                }
            }

            for (auto k = 0; k < STEV_CCSN_N_METALLICITY; ++k) {
                iTable = RowMajorIndex(k, iMass, 0, STEV_CCSN_N_METALLICITY, STEV_CCSN_N_MASS, 1);
                iData = RowMajorIndex(k, i, 0, STEV_CCSN_N_METALLICITY, STEV_INTERP_N_MASS, 1);

                /* WARNING: The following formulas are correct only when mass is the last index
                   in the tables */
                Data->afCCSNMetalYield[iData] = CCSN->pfMetalYield[iTable] * (1.0 - fDeltaMass) + CCSN->pfMetalYield[iTable + 1] * fDeltaMass;

                Data->afCCSNEjectedMass[iData] = CCSN->pfEjectedMass[iTable] * (1.0 - fDeltaMass) + CCSN->pfEjectedMass[iTable + 1] * fDeltaMass;
            }

            for (auto k = 0; k < STEV_AGB_N_METALLICITY; ++k) {
                iData = RowMajorIndex(k, i, 0, STEV_AGB_N_METALLICITY, STEV_INTERP_N_MASS, 1);
                if (i > iCCSNMinMass + 1) {
                    Data->afAGBMetalYield[iData] = 0.0f;
                    Data->afAGBEjectedMass[iData] = 0.0f;
                }
                else {
                    iTable = RowMajorIndex(k, STEV_AGB_N_MASS - 1, 0, STEV_AGB_N_METALLICITY, STEV_AGB_N_MASS, 1);
                    Data->afAGBMetalYield[iData] = AGB->pfMetalYield[iTable];
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
static inline void InterpolateXAxis(const float* restrict pfTable, const int nY, const int nZ, const int nSize, const int iOffset, const float fWeight,
                                    float* restrict pfResult) {

    const float* restrict pfLower = pfTable + iOffset * nZ;
    const float* restrict pfUpper = pfLower + nY * nZ;

    for (int i = 0; i < nSize * nZ; ++i) pfResult[i] = pfLower[i] * (1.0f - fWeight) + pfUpper[i] * fWeight;
}

static inline void ComputeAndCorrectSimulEjecta(const float* restrict pfYield, const float fMetalYield, const float fEjectedMass, const int nElems,
                                                const float* restrict pfElemAbun, const float fMetalAbun, float* restrict pfElemEjMass,
                                                float* restrict pfMetalEjMass) {

    int i;
    for (i = 0; i < nElems; ++i) pfElemEjMass[i] = pfYield[i] + pfElemAbun[i] * fEjectedMass;
    *pfMetalEjMass = fMetalYield + fMetalAbun * fEjectedMass;

    for (i = 0; i < nElems; ++i) {
        if (pfElemEjMass[i] < 0.0f) {
            if (i != ELEMENT_H && i != ELEMENT_He) *pfMetalEjMass -= pfElemEjMass[i];
            pfElemEjMass[i] = 0.0f;
        }
    }
    if (*pfMetalEjMass < 0.0f) *pfMetalEjMass = 0.0f;

    const float fTotalMass = pfElemEjMass[ELEMENT_H] + pfElemEjMass[ELEMENT_He] + *pfMetalEjMass;
    const float fNormFactor = fEjectedMass / fTotalMass;
    for (i = 0; i < nElems; ++i) pfElemEjMass[i] *= fNormFactor;
    *pfMetalEjMass *= fNormFactor;
}

template <int nZs, int nMasses, int nElems>
static inline void ComputeMassToEject(const float* restrict pfYield, const float* restrict pfMetalYield, const float* restrict pfEjectedMass,
                                      const float* restrict pfInitialMass, const float* restrict pfIMFLogWeight, const float fDeltaLogMass,
                                      const float* restrict pfElemAbun, const float fMetalAbun, const int iStart, const int iEnd, const float fMassStart,
                                      const float fMassEnd, const int iOffset, const float fWeight, float* pfElemMass, float* pfMetalMass) {

    int i, j, k;
    const int nSize = iEnd - iStart + 1;
    // NB: The number of pre-allocated elements may require some tuning
    // based on the typical sizes of these arrays
    boost::container::small_vector<float, 32> afInterpYield(nSize * nElems);
    boost::container::small_vector<float, 32> afInterpMetalYield(nSize);
    boost::container::small_vector<float, 32> afInterpEjectedMass(nSize);

    InterpolateXAxis(pfYield, nMasses, nElems, nSize, iOffset + iStart, fWeight, afInterpYield.data());
    InterpolateXAxis(pfMetalYield, nMasses, 1, nSize, iOffset + iStart, fWeight, afInterpMetalYield.data());
    InterpolateXAxis(pfEjectedMass, nMasses, 1, nSize, iOffset + iStart, fWeight, afInterpEjectedMass.data());

    /* In a strict sense, this two new declarations are not needed. We could overwrite
       two of the previously declared arrays, although the code would become less readable.
       If this were to be done, the restrict qualifier in the first and last two arguments
       of stev::ComputeAndCorrectSimulEjecta must be removed */
    boost::container::small_vector<float, 32> afElemSimulEjecta(nSize * nElems);
    boost::container::small_vector<float, 32> afMetalSimulEjecta(nSize);
    for (i = 0, j = 0; i < nSize; ++i, j += nElems) {
        ComputeAndCorrectSimulEjecta(afInterpYield.data() + j, afInterpMetalYield[i], afInterpEjectedMass[i], nElems, pfElemAbun, fMetalAbun,
                                     afElemSimulEjecta.data() + j, afMetalSimulEjecta.data() + i);
    }

    pfIMFLogWeight += iStart;
    pfInitialMass += iStart;

    for (i = 0, j = 0; i < nSize; ++i) {
        for (k = 0; k < nElems; ++j, ++k) afElemSimulEjecta[j] *= pfIMFLogWeight[i];
        afMetalSimulEjecta[i] *= pfIMFLogWeight[i];
    }

    /* Integration is done over the logarithms of the masses. Here, the steps for this are
       followed in a rather cumbersome way. It is intended to avoid round-off error as much
       as possible. */
    float afElemMassTemp[nElems], fMetalMassTemp;

    /* Contribution from the first and last values */
    for (i = 0; i < nElems; ++i) {
        afElemMassTemp[i] = 0.5f * (afElemSimulEjecta[i] + afElemSimulEjecta[i + (nSize - 1) * nElems]);
    }
    fMetalMassTemp = 0.5f * (afMetalSimulEjecta[0] + afMetalSimulEjecta[nSize - 1]);

    /* Contribution from the middle values */
    for (i = 1, j = nElems; i < nSize - 1; ++i) {
        for (k = 0; k < nElems; ++j, ++k) afElemMassTemp[k] += afElemSimulEjecta[j];
        fMetalMassTemp += afMetalSimulEjecta[i];
    }

    /* Multiply by the logarithm of the spacing */
    for (i = 0; i < nElems; ++i) afElemMassTemp[i] *= fDeltaLogMass;
    fMetalMassTemp *= fDeltaLogMass;

    /* Correction for initial and final values mismatch */
    const float fWeightStart = log10f(fMassStart / pfInitialMass[0]);
    const float fWeightEnd = log10f(pfInitialMass[nSize - 1] / fMassEnd);
    for (i = 0; i < nElems; ++i) {
        afElemMassTemp[i] -= 0.5f * (fWeightStart * (afElemSimulEjecta[i] + afElemSimulEjecta[i + nElems]) +
                                     fWeightEnd * (afElemSimulEjecta[i + (nSize - 2) * nElems] + afElemSimulEjecta[i + (nSize - 1) * nElems]));
    }
    fMetalMassTemp -=
        0.5f * (fWeightStart * (afMetalSimulEjecta[0] + afMetalSimulEjecta[1]) + fWeightEnd * (afMetalSimulEjecta[nSize - 2] + afMetalSimulEjecta[nSize - 1]));

    /* Multiply by natural logarithm of 10 */
    for (i = 0; i < nElems; ++i) pfElemMass[i] += M_LN10 * afElemMassTemp[i];
    *pfMetalMass += M_LN10 * fMetalMassTemp;
}

static inline float LifetimeFunction(PKD pkd, decorators::STAR& star, const double dMass) {
    int iMass;
    float fDeltaMass;
    GetIndex1D(pkd->StelEvolData->afLifetimeInitialMass, STEV_LIFETIME_N_MASS, log10(dMass), &iMass, &fDeltaMass);

    const float fDeltaZ = star.Lifetime.fDeltaZ;

    const float* afTimeLowerZ = pkd->StelEvolData->afLifetime + star.Lifetime.oZ;
    const float* afTimeUpperZ = afTimeLowerZ + STEV_LIFETIME_N_MASS;

    float fLogTime = afTimeLowerZ[iMass] * (1.0f - fDeltaZ) * (1.0f - fDeltaMass) + afTimeLowerZ[iMass + 1] * (1.0f - fDeltaZ) * fDeltaMass +
                     afTimeUpperZ[iMass] * fDeltaZ * (1.0f - fDeltaMass) + afTimeUpperZ[iMass + 1] * fDeltaZ * fDeltaMass;

    return powf(10.0f, fLogTime);
}

static inline float InverseLifetimeFunction(PKD pkd, decorators::STAR& star, const float fTime) {
    const float* afTimeLowerZ = pkd->StelEvolData->afLifetime + star.Lifetime.oZ;
    const float* afTimeUpperZ = afTimeLowerZ + STEV_LIFETIME_N_MASS;

    int iTimeLowerZ;
    float fDeltaTimeLowerZ;
    GetIndex1DReversed(afTimeLowerZ, STEV_LIFETIME_N_MASS, log10(fTime), &iTimeLowerZ, &fDeltaTimeLowerZ);

    int iTimeUpperZ;
    float fDeltaTimeUpperZ;
    GetIndex1DReversed(afTimeUpperZ, STEV_LIFETIME_N_MASS, log10(fTime), &iTimeUpperZ, &fDeltaTimeUpperZ);
    const float fDeltaZ = star.Lifetime.fDeltaZ;
    const float* afMass = pkd->StelEvolData->afLifetimeInitialMass;

    float fLogMass = afMass[iTimeLowerZ] * (1.0f - fDeltaZ) * (1.0f - fDeltaTimeLowerZ) + afMass[iTimeLowerZ + 1] * (1.0f - fDeltaZ) * fDeltaTimeLowerZ +
                     afMass[iTimeUpperZ] * fDeltaZ * (1.0f - fDeltaTimeUpperZ) + afMass[iTimeUpperZ + 1] * fDeltaZ * fDeltaTimeUpperZ;

    return powf(10.0f, fLogMass);
}

/* Function that estimates the time it will take a star particle to lose all
 * its initial mass by estimating its current ejecta rate. This is then used
 * to compute the time it needs to eject fMinMassFrac of its mass. */
static inline float ComputeFirstEnrichTime(PKD pkd, decorators::STAR& star, const double dCCSNMinMass) {
    const float fMinMassFrac = 1e-3f;

    const int iEnd = star.iLastEnrichMass;
    const int iStart = iEnd - 1;
    const float fMassStart = pkd->StelEvolData->afInitialMass[iStart];
    const float fMassEnd = pkd->StelEvolData->afInitialMass[iEnd];

    float afTblEjMass[2];
    if (fMassStart > dCCSNMinMass) {
        InterpolateXAxis(pkd->StelEvolData->afCCSNEjectedMass, STEV_INTERP_N_MASS, 1, 2, star.CCSN.oZ + iStart, star.CCSN.fDeltaZ, afTblEjMass);
    }
    else {
        InterpolateXAxis(pkd->StelEvolData->afAGBEjectedMass, STEV_INTERP_N_MASS, 1, 2, star.AGB.oZ + iStart, star.AGB.fDeltaZ, afTblEjMass);
    }
    const float* const pfIMFLogWeight = pkd->StelEvolData->afIMFLogWeight + iStart;
    const float fEjMassFrac = 0.5f * M_LN10 * pkd->StelEvolData->fDeltaLogMass * (afTblEjMass[0] * pfIMFLogWeight[0] + afTblEjMass[1] * pfIMFLogWeight[1]);

    const float fTimeStart = LifetimeFunction(pkd, star, fMassEnd);
    const float fTimeEnd = LifetimeFunction(pkd, star, fMassStart);
    const float fDepletionTime = (fTimeEnd - fTimeStart) / fEjMassFrac;
    float fNextTime = star.fLastEnrichTime + fMinMassFrac * fDepletionTime;

    /* Here we make sure fNextTime is beyond the numerical artifacts of the
     * Lifetime and Inverse Lifetime functions */
    float fMass = InverseLifetimeFunction(pkd, star, fNextTime);
    const float fStep = 1.05f;
    while (fMass > star.fLastEnrichMass) {
        fNextTime *= fStep;
        fMass = InverseLifetimeFunction(pkd, star, fNextTime);
    }

    return star.fTimer + fNextTime;
}

/* Given the ejecta rate of a star particle, use it to compute the time it needs to
 * eject fMinMassFrac of its mass. */
static inline float ComputeNextEnrichTime(const float fTime, const float fStarMass, const float fEjectedMass, const float fTimestep) {
    const float fMinMassFrac = 1e-3f;
    const float fEjMassRateInv = fTimestep / fEjectedMass;

    return fTime + fMinMassFrac * fStarMass * fEjMassRateInv;
}

void StarParticleInit(PKD pkd, decorators::STAR& star, const double dSNIaMaxMass, const double dCCSNMinMass, const double dCCSNMaxMass) {
    int iZ;
    float fLogZ;
    if (star.fMetalAbun > 0.0f) fLogZ = log10f(star.fMetalAbun);
    else fLogZ = STEV_MIN_LOG_METALLICITY;

    GetIndex1D(pkd->StelEvolData->afCCSNMetallicity, STEV_CCSN_N_METALLICITY, fLogZ, &iZ, &star.CCSN.fDeltaZ);
    star.CCSN.oZ = iZ * STEV_INTERP_N_MASS;
    GetIndex1D(pkd->StelEvolData->afAGBMetallicity, STEV_AGB_N_METALLICITY, fLogZ, &iZ, &star.AGB.fDeltaZ);
    star.AGB.oZ = iZ * STEV_INTERP_N_MASS;
    GetIndex1D(pkd->StelEvolData->afLifetimeMetallicity, STEV_LIFETIME_N_METALLICITY, fLogZ, &iZ, &star.Lifetime.fDeltaZ);
    star.Lifetime.oZ = iZ * STEV_LIFETIME_N_MASS;

    star.fSNIaOnsetTime = LifetimeFunction(pkd, star, dSNIaMaxMass);

    const float fCCSNOnsetTime = LifetimeFunction(pkd, star, dCCSNMaxMass);
    if (star.fLastEnrichTime < fCCSNOnsetTime) {
        star.fLastEnrichTime = fCCSNOnsetTime;
        star.fLastEnrichMass = dCCSNMaxMass;
        star.iLastEnrichMass = STEV_INTERP_N_MASS - 1;
    }
    else {
        star.fLastEnrichMass = InverseLifetimeFunction(pkd, star, star.fLastEnrichTime);
        if (star.fLastEnrichMass < dCCSNMaxMass) {
            star.iLastEnrichMass = GetIMFMassIndex(pkd->StelEvolData->afInitialMass, STEV_INTERP_N_MASS, star.fLastEnrichMass, STEV_INTERP_N_MASS - 1) + 1;
        }
        else {
            star.fLastEnrichTime = fCCSNOnsetTime;
            star.fLastEnrichMass = dCCSNMaxMass;
            star.iLastEnrichMass = STEV_INTERP_N_MASS - 1;
        }
    }

    star.fNextEnrichTime = ComputeFirstEnrichTime(pkd, star, dCCSNMinMass);
}
} // namespace stev
//**********************************************************************//

auto pkdStellarEvolutionInit(PKD pkd, const ServiceStellarEvolutionInit::input* in) -> int {
    pkd->StelEvolData = (STEV_DATA*)malloc(sizeof(STEV_DATA));
    assert(pkd->StelEvolData != NULL);
    memcpy(pkd->StelEvolData, &in->StelEvolData, sizeof(STEV_DATA));

    if (strcmp(in->achSNIaDTDType, "exponential") == 0) {
        pkd->StelEvolData->fcnNumSNIa = stev::ExponentialNumSNIa;
    }
    else if (strcmp(in->achSNIaDTDType, "powerlaw") == 0) {
        pkd->StelEvolData->fcnNumSNIa = stev::PowerlawNumSNIa;
    }
    else {
        std::cerr << "ERROR: Undefined SNIa DTD type has been given in " << "achSNIaDTDType parameter: " << in->achSNIaDTDType << std::endl;
        exit(1);
    }

    for (auto& p : pkd->particles) {

        if (p.is_star()) {
            auto& star = p.star();

            if (star.fInitialMass < 0.0f) star.fInitialMass = p.mass();

            if (star.fLastEnrichTime < 0.0f) star.fLastEnrichTime = 0.0f;

            if (star.fTimer <= 0.0f) star.fTimer = in->dTime;

            stev::StarParticleInit(pkd, star, in->dSNIaMaxMass, in->dCCSNMinMass, in->dCCSNMaxMass);
        }
        else if (p.is_gas()) {
            auto& sph = p.sph();
            sph.ReceivedMom = 0.0f;
            sph.fReceivedMass = 0.0f;
            sph.fReceivedE = 0.0f;
        }
    }

    return 0;
}

auto ServiceStellarEvolutionInit::Service(PST pst, std::span<std::byte> vin, std::span<std::byte> vout) -> int {
    auto pkd = pst->lcl()->pkd;
    auto in = reinterpret_cast<const struct input*>(vin.data());
    pkdStellarEvolutionInit(pkd, in);
    return 0;
}

void pkdAddStellarEjecta(PKD pkd, particleStore::ParticleReference& p, decorators::MFM& sph, const double dConstGamma) {
    const double dOldEkin = blitz::dot(sph.mom, sph.mom) / (2.0 * p.mass());

    sph.mom += sph.ReceivedMom;
    p.set_mass(p.mass() + sph.fReceivedMass);
    sph.E += sph.fReceivedE;

    const double dNewEkin = blitz::dot(sph.mom, sph.mom) / (2.0 * p.mass());

    const double dDeltaUint = sph.fReceivedE - (dNewEkin - dOldEkin);
    if (dDeltaUint > 0.0) {
        sph.Uint += dDeltaUint;
#ifdef ENTROPY_SWITCH
        sph.S += dDeltaUint * (dConstGamma - 1.0) * pow(p.density(), 1.0 - dConstGamma);
#endif
    }

    sph.ReceivedMom = 0.0f;
    sph.fReceivedMass = 0.0f;
    sph.fReceivedE = 0.0f;
}

void SmoothHelperChemEnrich::pack(void* dst, const void* src) {
    auto p1 = static_cast<pack_type*>(dst);
    auto p2 = pkd->particles[static_cast<const PARTICLE*>(src)];

    p1->iClass = p2.get_class();
    if (p2.is_gas()) {
        p1->position = p2.position();
        p1->fDensity = p2.density();
    }
}

void SmoothHelperChemEnrich::unpack(void* dst, const void* src, const void* key) {
    auto p1 = pkd->particles[static_cast<PARTICLE*>(dst)];
    auto p2 = static_cast<const pack_type*>(src);

    p1.set_class(p2->iClass);
    if (p1.is_gas()) {
        p1.set_adj_decorator();
        p1.set_position(p2->position);
        p1.set_density(p2->fDensity);
    }
}

void SmoothHelperChemEnrich::init(void* dst) {
    auto p = pkd->particles[static_cast<PARTICLE*>(dst)];

    if (p.is_gas()) {
        auto& sph = p.sph();

        sph.ReceivedMom = 0.0f;
        sph.fReceivedMass = 0.0f;
        sph.fReceivedE = 0.0f;

        sph.ElemMass = 0.0f;
        sph.fMetalMass = 0.0f;
    }
}

void SmoothHelperChemEnrich::flush(void* dst, const void* src) {
    auto p1 = static_cast<flush_type*>(dst);
    auto p2 = pkd->particles[static_cast<const PARTICLE*>(src)];

    if (p2.is_gas()) {
        const auto& sph = p2.sph();

        p1->ReceivedMom = sph.ReceivedMom;
        p1->fReceivedMass = sph.fReceivedMass;
        p1->fReceivedE = sph.fReceivedE;

        p1->ElemMass = sph.ElemMass;
        p1->fMetalMass = sph.fMetalMass;
    }
}

void SmoothHelperChemEnrich::combine(void* dst, const void* src, const void* key) {
    auto p1 = pkd->particles[static_cast<PARTICLE*>(dst)];
    auto p2 = static_cast<const flush_type*>(src);

    if (p1.is_gas()) {
        auto& sph = p1.sph();

        sph.ReceivedMom += p2->ReceivedMom;
        sph.fReceivedMass += p2->fReceivedMass;
        sph.fReceivedE += p2->fReceivedE;

        sph.ElemMass += p2->ElemMass;
        sph.fMetalMass += p2->fMetalMass;
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
   are in Msol. Conversely, star.fInitialMass and p.mass() are in code units.
   Furthermore, all times are in code units.
*/
void SmoothHelperChemEnrich::Smooth(particleStore::ParticleReference& p, float fBall, std::span<NN> nnList) {
    auto& star = p.star();

    const float fInitialTime = star.fLastEnrichTime;
    const float fInitialMass = star.fLastEnrichMass;
    const int iInitialMass = star.iLastEnrichMass;

    const float fFinalTime = (float)dTime - star.fTimer;
    const float fFinalMass = stev::InverseLifetimeFunction(pkd, star, fFinalTime);
    const int iFinalMass = stev::GetIMFMassIndex(pkd->StelEvolData->afInitialMass, STEV_INTERP_N_MASS, fFinalMass, iInitialMass);

    star.fLastEnrichTime = fFinalTime;
    star.fLastEnrichMass = fFinalMass;
    star.iLastEnrichMass = iFinalMass + 1;

    /* Note: The parameter star.[AGB,CCSN,Lifetime].oZ contains the index of the
       interpolation's lower metallicity array multiplied by the number of mass bins.
       Since this multiplication must always be made, it is done once and for all in
       the function stev::StarParticleInit. */

    blitz::TinyVector<float, ELEMENT_COUNT> ElemMass{0.0f};
    float fMetalMass{0.0f};

    const float fTransMass = dCCSNMinMass;
    if (fFinalMass > fTransMass) {
        /* CCSN only */
        stev::ComputeMassToEject<STEV_CCSN_N_METALLICITY, STEV_INTERP_N_MASS, ELEMENT_COUNT>(
            pkd->StelEvolData->afCCSNYield, pkd->StelEvolData->afCCSNMetalYield, pkd->StelEvolData->afCCSNEjectedMass, pkd->StelEvolData->afInitialMass,
            pkd->StelEvolData->afIMFLogWeight, pkd->StelEvolData->fDeltaLogMass, star.ElemAbun.data(), star.fMetalAbun, iFinalMass, iInitialMass, fFinalMass,
            fInitialMass, star.CCSN.oZ, star.CCSN.fDeltaZ, ElemMass.data(), &fMetalMass);
    }
    else if (fInitialMass < fTransMass) {
        /* AGB only */
        stev::ComputeMassToEject<STEV_AGB_N_METALLICITY, STEV_INTERP_N_MASS, ELEMENT_COUNT>(
            pkd->StelEvolData->afAGBYield, pkd->StelEvolData->afAGBMetalYield, pkd->StelEvolData->afAGBEjectedMass, pkd->StelEvolData->afInitialMass,
            pkd->StelEvolData->afIMFLogWeight, pkd->StelEvolData->fDeltaLogMass, star.ElemAbun.data(), star.fMetalAbun, iFinalMass, iInitialMass, fFinalMass,
            fInitialMass, star.AGB.oZ, star.AGB.fDeltaZ, ElemMass.data(), &fMetalMass);
    }
    else {
        /* Mixed CCSN and AGB */
        int iTransMass = stev::GetIMFMassIndex(pkd->StelEvolData->afInitialMass, STEV_INTERP_N_MASS, fTransMass, iInitialMass);

        stev::ComputeMassToEject<STEV_CCSN_N_METALLICITY, STEV_INTERP_N_MASS, ELEMENT_COUNT>(
            pkd->StelEvolData->afCCSNYield, pkd->StelEvolData->afCCSNMetalYield, pkd->StelEvolData->afCCSNEjectedMass, pkd->StelEvolData->afInitialMass,
            pkd->StelEvolData->afIMFLogWeight, pkd->StelEvolData->fDeltaLogMass, star.ElemAbun.data(), star.fMetalAbun, iTransMass, iInitialMass, fTransMass,
            fInitialMass, star.CCSN.oZ, star.CCSN.fDeltaZ, ElemMass.data(), &fMetalMass);

        ++iTransMass;
        stev::ComputeMassToEject<STEV_AGB_N_METALLICITY, STEV_INTERP_N_MASS, ELEMENT_COUNT>(
            pkd->StelEvolData->afAGBYield, pkd->StelEvolData->afAGBMetalYield, pkd->StelEvolData->afAGBEjectedMass, pkd->StelEvolData->afInitialMass,
            pkd->StelEvolData->afIMFLogWeight, pkd->StelEvolData->fDeltaLogMass, star.ElemAbun.data(), star.fMetalAbun, iFinalMass, iTransMass, fFinalMass,
            fTransMass, star.AGB.oZ, star.AGB.fDeltaZ, ElemMass.data(), &fMetalMass);
    }

    const float fNumSNIa = (*pkd->StelEvolData->fcnNumSNIa)(dSNIaNorm, dSNIaScale, star, fInitialTime, fFinalTime);
    for (auto i = 0; i < ELEMENT_COUNT; ++i) {
        ElemMass[i] += fNumSNIa * pkd->StelEvolData->afSNIaEjectedMass[i];
    }
    ElemMass *= star.fInitialMass;
    fMetalMass += fNumSNIa * pkd->StelEvolData->fSNIaEjectedMetalMass;
    fMetalMass *= star.fInitialMass;

    const double dTotalMass = (double)ElemMass[ELEMENT_H] + ElemMass[ELEMENT_He] + fMetalMass;
    star.fNextEnrichTime = stev::ComputeNextEnrichTime(dTime, star.fInitialMass, dTotalMass, fFinalTime - fInitialTime);
    p.set_mass(p.mass() - dTotalMass);
    assert(p.mass() > 0.0f);

    const double dScaleFactorInv = 1.0 / csmTime2Exp(pkd->csm, dTime);
    const double dScaleFactorInvSq = dScaleFactorInv * dScaleFactorInv;
    const auto& StarVel = p.velocity();

    const double dStarDeltaEkin = 0.5 * dTotalMass * blitz::dot(StarVel, StarVel);
    const double dWindEkin = dWindSpecificEkin * dTotalMass;
    const double dStarEjEnergy = dStarDeltaEkin * dScaleFactorInvSq + dWindEkin;

    const double dWeight = 1.0 / nnList.size();
    const double dDeltaMass = dWeight * dTotalMass;
    const double dDeltaE = dWeight * dStarEjEnergy;

    ElemMass *= dWeight;
    fMetalMass *= dWeight;

    for (auto& nn : nnList) {
        auto q = pkd->particles[nn.pPart];
        assert(q.is_gas());
        auto& qsph = q.sph();

        qsph.ReceivedMom += dDeltaMass * StarVel * dScaleFactorInv;
        qsph.fReceivedMass += dDeltaMass;
        qsph.fReceivedE += dDeltaE;

        qsph.ElemMass += ElemMass;
        qsph.fMetalMass += fMetalMass;
    }
}

#define H5FIELD_MASS "Masses"
#define H5FIELD_METALLICITY "Metallicities"
#define H5FIELD_SPECNAME "Species_names"
#define H5FIELD_TABLEZNAME "Yield_names"
#define H5FIELD_TABLE "Yields"
#define H5FIELD_YIELD "Yield"
#define H5FIELD_METALYIELD "Total_Metals"
#define H5FIELD_EJMASS "Ejected_mass"
#define H5FIELD_LIFETIME "Lifetimes"
namespace stev {
STEV_RAWDATA* ReadTable(char* pszPath) {
    hid_t fileID, datatype, dataspace;
    herr_t status;

    STEV_RAWDATA* RawData = (STEV_RAWDATA*)malloc(sizeof(STEV_RAWDATA));
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
    RawData->nElems = H5Sget_simple_extent_npoints(dataspace);
    status = H5Sclose(dataspace);
    assert(status >= 0);
    status = H5Dclose(dsSpecName);
    assert(status >= 0);

    RawData->pfMetallicity = (float*)malloc(RawData->nZs * sizeof(float));
    assert(RawData->pfMetallicity != NULL);
    RawData->pfInitialMass = (float*)malloc(RawData->nMasses * sizeof(float));
    assert(RawData->pfInitialMass != NULL);
    RawData->pfYield = (float*)malloc(RawData->nZs * RawData->nElems * RawData->nMasses * sizeof(float));
    assert(RawData->pfYield != NULL);
    RawData->pfMetalYield = (float*)malloc(RawData->nZs * RawData->nMasses * sizeof(float));
    assert(RawData->pfMetalYield != NULL);
    RawData->pfEjectedMass = (float*)malloc(RawData->nZs * RawData->nMasses * sizeof(float));
    assert(RawData->pfEjectedMass != NULL);

    status = H5Dread(dsMetal, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->pfMetallicity);
    assert(status >= 0);
    status = H5Dclose(dsMetal);
    assert(status >= 0);

    status = H5Dread(dsMass, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->pfInitialMass);
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
    std::vector<char*> apchTableZNames(RawData->nZs);
    datatype = H5Dget_type(dsTableZName);
    assert(datatype >= 0);
    status = H5Dread(dsTableZName, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, apchTableZNames.data());
    assert(status >= 0);
    status = H5Tclose(datatype);
    assert(status >= 0);
    status = H5Dclose(dsTableZName);
    assert(status >= 0);

    char achTable[64];
    int iOffset, nDims, nSize;
    hsize_t dimSize[2];
    for (int i = 0; i < RawData->nZs; ++i) {
        iOffset = i * RawData->nElems * RawData->nMasses;

        snprintf(achTable, sizeof(achTable), "%s/%s/%s", H5FIELD_TABLE, apchTableZNames[i], H5FIELD_YIELD);
        hid_t dsYield = H5Dopen(fileID, achTable, H5P_DEFAULT);
        assert(dsYield >= 0);
        dataspace = H5Dget_space(dsYield);
        assert(dataspace >= 0);
        nDims = H5Sget_simple_extent_dims(dataspace, dimSize, NULL);
        assert(nDims == 2);
        assert(dimSize[0] == RawData->nElems);
        assert(dimSize[1] == RawData->nMasses);
        status = H5Sclose(dataspace);
        assert(status >= 0);
        status = H5Dread(dsYield, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->pfYield + iOffset);
        assert(status >= 0);
        status = H5Dclose(dsYield);
        assert(status >= 0);

        iOffset = i * RawData->nMasses;

        snprintf(achTable, sizeof(achTable), "%s/%s/%s", H5FIELD_TABLE, apchTableZNames[i], H5FIELD_METALYIELD);
        hid_t dsMetalYield = H5Dopen(fileID, achTable, H5P_DEFAULT);
        assert(dsMetalYield >= 0);
        dataspace = H5Dget_space(dsMetalYield);
        assert(dataspace >= 0);
        nSize = H5Sget_simple_extent_npoints(dataspace);
        assert(nSize == RawData->nMasses);
        status = H5Sclose(dataspace);
        assert(status >= 0);
        status = H5Dread(dsMetalYield, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->pfMetalYield + iOffset);
        assert(status >= 0);
        status = H5Dclose(dsMetalYield);
        assert(status >= 0);

        snprintf(achTable, sizeof(achTable), "%s/%s/%s", H5FIELD_TABLE, apchTableZNames[i], H5FIELD_EJMASS);
        hid_t dsEjMass = H5Dopen(fileID, achTable, H5P_DEFAULT);
        assert(dsEjMass >= 0);
        dataspace = H5Dget_space(dsEjMass);
        assert(dataspace >= 0);
        nSize = H5Sget_simple_extent_npoints(dataspace);
        assert(nSize == RawData->nMasses);
        status = H5Sclose(dataspace);
        assert(status >= 0);
        status = H5Dread(dsEjMass, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->pfEjectedMass + iOffset);
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
STEV_RAWDATA* ReadSNIaTable(char* pszPath) {
    hid_t fileID, dataspace;
    herr_t status;

    STEV_RAWDATA* RawData = (STEV_RAWDATA*)malloc(sizeof(STEV_RAWDATA));
    assert(RawData != NULL);

    fileID = H5Fopen(pszPath, H5F_ACC_RDONLY, H5P_DEFAULT);
    assert(fileID >= 0);

    hid_t dsYield = H5Dopen(fileID, H5FIELD_YIELD, H5P_DEFAULT);
    assert(dsYield >= 0);
    dataspace = H5Dget_space(dsYield);
    assert(dataspace >= 0);
    RawData->nElems = H5Sget_simple_extent_npoints(dataspace);
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

    RawData->pfEjectedMass = (float*)malloc(RawData->nElems * sizeof(float));
    assert(RawData->pfEjectedMass != NULL);
    RawData->pfMetalYield = (float*)malloc(sizeof(float));
    assert(RawData->pfMetalYield != NULL);

    status = H5Dread(dsYield, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->pfEjectedMass);
    assert(status >= 0);
    status = H5Dclose(dsYield);
    assert(status >= 0);

    status = H5Dread(dsMetalYield, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->pfMetalYield);
    assert(status >= 0);
    status = H5Dclose(dsMetalYield);
    assert(status >= 0);

    status = H5Fclose(fileID);
    assert(status >= 0);

    return RawData;
}

STEV_RAWDATA* ReadLifetimeTable(char* pszPath) {
    hid_t fileID, dataspace;
    herr_t status;

    STEV_RAWDATA* RawData = (STEV_RAWDATA*)malloc(sizeof(STEV_RAWDATA));
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

    RawData->pfMetallicity = (float*)malloc(RawData->nZs * sizeof(float));
    assert(RawData->pfMetallicity != NULL);
    RawData->pfInitialMass = (float*)malloc(RawData->nMasses * sizeof(float));
    assert(RawData->pfInitialMass != NULL);
    RawData->pfLifetime = (float*)malloc(RawData->nZs * RawData->nMasses * sizeof(float));
    assert(RawData->pfLifetime != NULL);

    status = H5Dread(dsMetal, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->pfMetallicity);
    assert(status >= 0);
    status = H5Dclose(dsMetal);
    assert(status >= 0);

    status = H5Dread(dsMass, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->pfInitialMass);
    assert(status >= 0);
    status = H5Dclose(dsMass);
    assert(status >= 0);

    hsize_t dimSize[2];
    hid_t dsLifetime = H5Dopen(fileID, H5FIELD_LIFETIME, H5P_DEFAULT);
    assert(dsLifetime >= 0);
    dataspace = H5Dget_space(dsLifetime);
    assert(dataspace >= 0);
    int nDims = H5Sget_simple_extent_dims(dataspace, dimSize, NULL);
    assert(nDims == 2);
    assert(dimSize[0] == RawData->nZs);
    assert(dimSize[1] == RawData->nMasses);
    status = H5Sclose(dataspace);
    assert(status >= 0);
    status = H5Dread(dsLifetime, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, RawData->pfLifetime);
    assert(status >= 0);
    status = H5Dclose(dsLifetime);
    assert(status >= 0);

    status = H5Fclose(fileID);
    assert(status >= 0);

    return RawData;
}

void FreeTable(STEV_RAWDATA* RawData) {
    free(RawData->pfMetallicity);
    free(RawData->pfInitialMass);
    free(RawData->pfYield);
    free(RawData->pfMetalYield);
    free(RawData->pfEjectedMass);
    free(RawData);
}

void FreeSNIaTable(STEV_RAWDATA* RawData) {
    free(RawData->pfMetalYield);
    free(RawData->pfEjectedMass);
    free(RawData);
}

void FreeLifetimeTable(STEV_RAWDATA* RawData) {
    free(RawData->pfMetallicity);
    free(RawData->pfInitialMass);
    free(RawData->pfLifetime);
    free(RawData);
}
} // namespace stev
