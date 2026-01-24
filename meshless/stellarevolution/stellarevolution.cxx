
#include <hdf5.h>

#include "stellarevolution.h"
#include "master.h"
#include "imf.h"
#include "hydro/hydro.h"

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

int pstStellarEvolutionInit(PST pst, void* vin, int nIn, void* vout, int nOut) {
    mdlassert(pst->mdl(), nIn == sizeof(struct inStellarEvolutionInit));
    if (pst->NotCore()) {
        int rID = pst->mdl()->ReqService(pst->upper(), PST_STELLAREVOLUTIONINIT, vin, nIn);
        pstStellarEvolutionInit(pst->lower(), vin, nIn, NULL, 0);
        pst->mdl()->GetReply(rID);
    }
    else {
        struct inStellarEvolutionInit* in = (struct inStellarEvolutionInit*)vin;
        pkdStellarEvolutionInit(pst->lcl()->pkd, in);
    }

    return 0;
}

int pkdStellarEvolutionInit(PKD pkd, struct inStellarEvolutionInit* in) {
    pkd->StelEvolData = (STEV_DATA*)malloc(sizeof(STEV_DATA));
    assert(pkd->StelEvolData != NULL);
    memcpy(pkd->StelEvolData, &in->StelEvolData, sizeof(STEV_DATA));

    if (strcmp(in->achSNIaDTDType, "exponential") == 0) {
        pkd->StelEvolData->fcnNumSNIa = stevExponentialNumSNIa;
    }
    else if (strcmp(in->achSNIaDTDType, "powerlaw") == 0) {
        pkd->StelEvolData->fcnNumSNIa = stevPowerlawNumSNIa;
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

            stevStarParticleInit(pkd, star, in->dSNIaMaxMass, in->dCCSNMinMass, in->dCCSNMaxMass);
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
    const float fFinalMass = stevInverseLifetimeFunction(pkd, star, fFinalTime);
    const int iFinalMass = stevGetIMFMassIndex(pkd->StelEvolData->afInitialMass, STEV_INTERP_N_MASS, fFinalMass, iInitialMass);

    star.fLastEnrichTime = fFinalTime;
    star.fLastEnrichMass = fFinalMass;
    star.iLastEnrichMass = iFinalMass + 1;

    /* Note: The parameter star.[AGB,CCSN,Lifetime].oZ contains the index of the
       interpolation's lower metallicity array multiplied by the number of mass bins.
       Since this multiplication must always be made, it is done once and for all in
       the function stevStarParticleInit. */

    blitz::TinyVector<float, ELEMENT_COUNT> ElemMass{0.0f};
    float fMetalMass{0.0f};

    const float fTransMass = dCCSNMinMass;
    if (fFinalMass > fTransMass) {
        /* CCSN only */
        stevComputeMassToEject<STEV_CCSN_N_METALLICITY, STEV_INTERP_N_MASS, ELEMENT_COUNT>(
            pkd->StelEvolData->afCCSNYield, pkd->StelEvolData->afCCSNMetalYield, pkd->StelEvolData->afCCSNEjectedMass, pkd->StelEvolData->afInitialMass,
            pkd->StelEvolData->afIMFLogWeight, pkd->StelEvolData->fDeltaLogMass, star.ElemAbun.data(), star.fMetalAbun, iFinalMass, iInitialMass, fFinalMass,
            fInitialMass, star.CCSN.oZ, star.CCSN.fDeltaZ, ElemMass.data(), &fMetalMass);
    }
    else if (fInitialMass < fTransMass) {
        /* AGB only */
        stevComputeMassToEject<STEV_AGB_N_METALLICITY, STEV_INTERP_N_MASS, ELEMENT_COUNT>(
            pkd->StelEvolData->afAGBYield, pkd->StelEvolData->afAGBMetalYield, pkd->StelEvolData->afAGBEjectedMass, pkd->StelEvolData->afInitialMass,
            pkd->StelEvolData->afIMFLogWeight, pkd->StelEvolData->fDeltaLogMass, star.ElemAbun.data(), star.fMetalAbun, iFinalMass, iInitialMass, fFinalMass,
            fInitialMass, star.AGB.oZ, star.AGB.fDeltaZ, ElemMass.data(), &fMetalMass);
    }
    else {
        /* Mixed CCSN and AGB */
        int iTransMass = stevGetIMFMassIndex(pkd->StelEvolData->afInitialMass, STEV_INTERP_N_MASS, fTransMass, iInitialMass);

        stevComputeMassToEject<STEV_CCSN_N_METALLICITY, STEV_INTERP_N_MASS, ELEMENT_COUNT>(
            pkd->StelEvolData->afCCSNYield, pkd->StelEvolData->afCCSNMetalYield, pkd->StelEvolData->afCCSNEjectedMass, pkd->StelEvolData->afInitialMass,
            pkd->StelEvolData->afIMFLogWeight, pkd->StelEvolData->fDeltaLogMass, star.ElemAbun.data(), star.fMetalAbun, iTransMass, iInitialMass, fTransMass,
            fInitialMass, star.CCSN.oZ, star.CCSN.fDeltaZ, ElemMass.data(), &fMetalMass);

        ++iTransMass;
        stevComputeMassToEject<STEV_AGB_N_METALLICITY, STEV_INTERP_N_MASS, ELEMENT_COUNT>(
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
    star.fNextEnrichTime = stevComputeNextEnrichTime(dTime, star.fInitialMass, dTotalMass, fFinalTime - fInitialTime);
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

STEV_RAWDATA* stevReadTable(char* pszPath) {
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
STEV_RAWDATA* stevReadSNIaTable(char* pszPath) {
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

STEV_RAWDATA* stevReadLifetimeTable(char* pszPath) {
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

void stevFreeTable(STEV_RAWDATA* RawData) {
    free(RawData->pfMetallicity);
    free(RawData->pfInitialMass);
    free(RawData->pfYield);
    free(RawData->pfMetalYield);
    free(RawData->pfEjectedMass);
    free(RawData);
}

void stevFreeSNIaTable(STEV_RAWDATA* RawData) {
    free(RawData->pfMetalYield);
    free(RawData->pfEjectedMass);
    free(RawData);
}

void stevFreeLifetimeTable(STEV_RAWDATA* RawData) {
    free(RawData->pfMetallicity);
    free(RawData->pfInitialMass);
    free(RawData->pfLifetime);
    free(RawData);
}

float stevExponentialNumSNIa(double dSNIaNorm, double dSNIaScale, decorators::STAR& star, float fInitialTime, float fFinalTime) {
    if (fFinalTime <= star.fSNIaOnsetTime) {
        return 0.0f;
    }
    else if (fInitialTime < star.fSNIaOnsetTime) {
        fInitialTime = star.fSNIaOnsetTime;
    }

    return (float)dSNIaNorm * (expf(-(fInitialTime / (float)dSNIaScale)) - expf(-(fFinalTime / (float)dSNIaScale)));
}

float stevPowerlawNumSNIa(double dSNIaNorm, double dSNIaScale, decorators::STAR& star, float fInitialTime, float fFinalTime) {
    if (fFinalTime <= star.fSNIaOnsetTime) {
        return 0.0f;
    }
    else if (fInitialTime < star.fSNIaOnsetTime) {
        fInitialTime = star.fSNIaOnsetTime;
    }

    return (float)dSNIaNorm * (powf(fFinalTime, (float)dSNIaScale + 1.0f) - powf(fInitialTime, (float)dSNIaScale + 1.0f));
}
