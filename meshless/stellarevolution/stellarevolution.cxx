#ifdef STELLAR_EVOLUTION

#include <hdf5.h>

#include "stellarevolution.h"
#include "master.h"
#include "imf.h"
#include "hydro/hydro.h"

using blitz::TinyVector;

void MSR::SetStellarEvolutionParam() {
    const double dYrToTime = SECONDSPERYEAR / units.dSecUnit;
    auto achSNIaDTDType = parameters.get_achSNIaDTDType();

    if (achSNIaDTDType == "exponential") {
        calc.dSNIaNorm = parameters.get_dSNIaNumPerMass();
        calc.dSNIaScale = parameters.get_dSNIaExpScale() * dYrToTime;
    }
    else if (achSNIaDTDType == "powerlaw") {
        auto dSNIaPLScale = parameters.get_dSNIaPLScale();
        calc.dSNIaNorm = parameters.get_dSNIaNumPerMass() /
                         (pow(parameters.get_dSNIaPLFinalTime() * dYrToTime, dSNIaPLScale + 1.0) -
                          pow(parameters.get_dSNIaPLInitTime() * dYrToTime, dSNIaPLScale + 1.0));
        calc.dSNIaScale = dSNIaPLScale;
    }
    else {
        std::cerr << "ERROR: Undefined SNIa DTD type has been given in " <<
                  "achSNIaDTDType parameter: " << achSNIaDTDType << std::endl;
        Exit(1);
    }

    calc.dWindSpecificEkin = 0.5 * pow(parameters.get_dStellarWindSpeed() / units.dKmPerSecUnit, 2);
    /* The number of gas particles to enrich is set to the average number of
       neighbours within a smoothing length. The factor 0.5 comes from the cubic
       spline kernel used by the hydro */
    calc.nSmoothEnrich = 0.5 * parameters.get_nSmooth();
}

void MSR::StellarEvolutionInit(double dTime) {
    char achPath[280];
    struct inStellarEvolutionInit in;
    STEV_RAWDATA *CCSNData, *AGBData, *SNIaData, *LifetimeData;
    auto achStelEvolPath = parameters.get_achStelEvolPath();

    /* Read the tables */
    snprintf(achPath, sizeof(achPath), "%s/CCSN.hdf5", achStelEvolPath.data());
    CCSNData = stevReadTable(achPath);
    assert(CCSNData->nZs == STEV_CCSN_N_METALLICITY);
    assert(CCSNData->nMasses == STEV_CCSN_N_MASS);
    assert(CCSNData->nElems == ELEMENT_COUNT);

    snprintf(achPath, sizeof(achPath), "%s/AGB.hdf5", achStelEvolPath.data());
    AGBData = stevReadTable(achPath);
    assert(AGBData->nZs == STEV_AGB_N_METALLICITY);
    assert(AGBData->nMasses == STEV_AGB_N_MASS);
    assert(AGBData->nElems == ELEMENT_COUNT);

    snprintf(achPath, sizeof(achPath), "%s/SNIa.hdf5", achStelEvolPath.data());
    SNIaData = stevReadSNIaTable(achPath);
    assert(SNIaData->nElems == ELEMENT_COUNT);

    snprintf(achPath, sizeof(achPath), "%s/Lifetimes.hdf5", achStelEvolPath.data());
    LifetimeData = stevReadLifetimeTable(achPath);
    assert(LifetimeData->nZs == STEV_LIFETIME_N_METALLICITY);
    assert(LifetimeData->nMasses == STEV_LIFETIME_N_MASS);

    /* NOTE: The lowest value of the initial mass array is set to the corresponding value
       of the Lifetime table that is being used, while its highest value to the maximum
       stellar mass for a CCSN to occur, as set in dCCSNMaxMass. The IMF is still
       normalized in the range [dIMFMinMass,dIMFMaxMass] */
    const double dMinMass = LifetimeData->pfInitialMass[0];
    const double dMaxMass = parameters.get_dCCSNMaxMass();
    double adInitialMass[STEV_INTERP_N_MASS], adIMF[STEV_INTERP_N_MASS];

    auto IMF = ChooseIMF(parameters.get_achIMFType().data(), parameters.get_dIMFMinMass(), parameters.get_dIMFMaxMass());
    IMF->MassWeightedSample(dMinMass, dMaxMass, STEV_INTERP_N_MASS, adInitialMass, adIMF);

    /* Store the data */
    for (auto i = 0; i < STEV_INTERP_N_MASS; ++i) {
        in.StelEvolData.afInitialMass[i] = adInitialMass[i];
        in.StelEvolData.afIMFLogWeight[i] = adIMF[i];
    }
    const double dDeltaLog = (log10(dMaxMass) - log10(dMinMass)) / (STEV_INTERP_N_MASS - 1);
    in.StelEvolData.fDeltaLogMass = dDeltaLog;

    /* Convert CCSN/AGB initial mass arrays to log for interpolation */
    for (auto i = 0; i < CCSNData->nMasses; ++i)
        CCSNData->pfInitialMass[i] = log10(CCSNData->pfInitialMass[i]);
    for (auto i = 0; i < AGBData->nMasses; ++i)
        AGBData->pfInitialMass[i] = log10(AGBData->pfInitialMass[i]);

    /* Interpolate yields and ejected masses to the initial mass array */
    stevInterpToIMFSampling(&in.StelEvolData, CCSNData, AGBData, parameters.get_dCCSNMinMass());

    for (auto i = 0; i < STEV_CCSN_N_METALLICITY; ++i) {
        if (CCSNData->pfMetallicity[i] > 0.0f)
            in.StelEvolData.afCCSNMetallicity[i] = log10(CCSNData->pfMetallicity[i]);
        else
            in.StelEvolData.afCCSNMetallicity[i] = STEV_MIN_LOG_METALLICITY;
    }

    for (auto i = 0; i < STEV_AGB_N_METALLICITY; ++i) {
        if (AGBData->pfMetallicity[i] > 0.0f)
            in.StelEvolData.afAGBMetallicity[i] = log10(AGBData->pfMetallicity[i]);
        else
            in.StelEvolData.afAGBMetallicity[i] = STEV_MIN_LOG_METALLICITY;
    }

    for (auto i = 0; i < ELEMENT_COUNT; ++i)
        in.StelEvolData.afSNIaEjectedMass[i] = SNIaData->pfEjectedMass[i];
    in.StelEvolData.fSNIaEjectedMetalMass = *SNIaData->pfMetalYield;

    for (auto i = 0; i < STEV_LIFETIME_N_METALLICITY; ++i) {
        if (LifetimeData->pfMetallicity[i] > 0.0f)
            in.StelEvolData.afLifetimeMetallicity[i] = log10(LifetimeData->pfMetallicity[i]);
        else
            in.StelEvolData.afLifetimeMetallicity[i] = STEV_MIN_LOG_METALLICITY;
    }
    for (auto i = 0; i < STEV_LIFETIME_N_MASS; ++i)
        in.StelEvolData.afLifetimeInitialMass[i] = log10(LifetimeData->pfInitialMass[i]);
    for (auto i = 0; i < STEV_LIFETIME_N_METALLICITY * STEV_LIFETIME_N_MASS; ++i) {
        in.StelEvolData.afLifetime[i] = log10(LifetimeData->pfLifetime[i] * SECONDSPERYEAR /
                                              units.dSecUnit);
    }

    strcpy(in.achSNIaDTDType, parameters.get_achSNIaDTDType().data());
    in.dTime = dTime;
    in.dSNIaMaxMass = parameters.get_dSNIaMaxMass();
    in.dCCSNMinMass = parameters.get_dCCSNMinMass();
    in.dCCSNMaxMass = parameters.get_dCCSNMaxMass();

    /* Send data and finish */
    pstStellarEvolutionInit(pst, &in, sizeof(struct inStellarEvolutionInit), NULL, 0);

    stevFreeTable(CCSNData);
    stevFreeTable(AGBData);
    stevFreeSNIaTable(SNIaData);
    stevFreeLifetimeTable(LifetimeData);
}

int pstStellarEvolutionInit(PST pst, void *vin, int nIn, void *vout, int nOut) {
    mdlassert(pst->mdl, nIn == sizeof(struct inStellarEvolutionInit));
    if (pst->nLeaves > 1) {
        int rID = pst->mdl->ReqService(pst->idUpper, PST_STELLAREVOLUTIONINIT, vin, nIn);
        pstStellarEvolutionInit(pst->pstLower, vin, nIn, NULL, 0);
        pst->mdl->GetReply(rID);
    }
    else {
        struct inStellarEvolutionInit *in = (struct inStellarEvolutionInit *) vin;
        pkdStellarEvolutionInit(pst->plcl->pkd, in);
    }

    return 0;
}

int pkdStellarEvolutionInit(PKD pkd, struct inStellarEvolutionInit *in) {
    pkd->StelEvolData = (STEV_DATA *) malloc(sizeof(STEV_DATA));
    assert(pkd->StelEvolData != NULL);
    memcpy(pkd->StelEvolData, &in->StelEvolData, sizeof(STEV_DATA));

    if (strcmp(in->achSNIaDTDType, "exponential") == 0) {
        pkd->StelEvolData->fcnNumSNIa = stevExponentialNumSNIa;
    }
    else if (strcmp(in->achSNIaDTDType, "powerlaw") == 0) {
        pkd->StelEvolData->fcnNumSNIa = stevPowerlawNumSNIa;
    }
    else {
        std::cerr << "ERROR: Undefined SNIa DTD type has been given in " <<
                  "achSNIaDTDType parameter: " << in->achSNIaDTDType << std::endl;
        exit(1);
    }

    for (auto &p : pkd->particles) {

        if (p.is_star()) {
            auto &star = p.star();

            if (star.fInitialMass < 0.0f)
                star.fInitialMass = p.mass();

            if (star.fLastEnrichTime < 0.0f)
                star.fLastEnrichTime = 0.0f;

            if (star.fTimer <= 0.0f)
                star.fTimer = in->dTime;

            stevStarParticleInit(pkd, star, in->dSNIaMaxMass, in->dCCSNMinMass,
                                 in->dCCSNMaxMass);
        }
        else if (p.is_gas()) {
            auto &sph = p.sph();
            sph.ReceivedMom = 0.0f;
            sph.fReceivedMass = 0.0f;
            sph.fReceivedE = 0.0f;
        }
    }

    return 0;
}

void pkdAddStellarEjecta(PKD pkd, particleStore::ParticleReference &p, meshless::FIELDS &sph,
                         const double dConstGamma) {
    const double dOldEkin = blitz::dot(sph.mom, sph.mom) / (2.0 * p.mass());

    sph.mom += sph.ReceivedMom;
    p.set_mass(p.mass() + sph.fReceivedMass);
    sph.E += sph.fReceivedE;

    const double dNewEkin = blitz::dot(sph.mom, sph.mom) / (2.0 * p.mass());

    const double dDeltaUint = sph.fReceivedE - (dNewEkin - dOldEkin);
    if (dDeltaUint > 0.0) {
        sph.Uint += dDeltaUint;
#ifdef ENTROPY_SWITCH
        sph.S += dDeltaUint * (dConstGamma - 1.0) *
                 pow(p.density(), 1.0 - dConstGamma);
#endif
    }

    sph.ReceivedMom = 0.0f;
    sph.fReceivedMass = 0.0f;
    sph.fReceivedE = 0.0f;
}

void packChemEnrich(void *vpkd, void *dst, const void *src) {
    PKD pkd = (PKD) vpkd;
    auto p1 = static_cast<stevPack *>(dst);
    auto p2 = pkd->particles[static_cast<const PARTICLE *>(src)];

    p1->iClass = p2.get_class();
    if (p2.is_gas()) {
        p1->position = p2.position();
        p1->fDensity = p2.density();
    }
}

void unpackChemEnrich(void *vpkd, void *dst, const void *src) {
    PKD pkd = (PKD) vpkd;
    auto p1 = pkd->particles[static_cast<PARTICLE *>(dst)];
    auto p2 = static_cast<const stevPack *>(src);

    p1.set_class(p2->iClass);
    if (p1.is_gas()) {
        p1.set_position(p2->position);
        p1.set_density(p2->fDensity);
    }
}

void initChemEnrich(void *vpkd, void *dst) {
    PKD pkd = (PKD) vpkd;
    auto p = pkd->particles[static_cast<PARTICLE *>(dst)];

    if (p.is_gas()) {
        auto &sph = p.sph();

        sph.ReceivedMom = 0.0f;
        sph.fReceivedMass = 0.0f;
        sph.fReceivedE = 0.0f;

        sph.ElemMass = 0.0f;
        sph.fMetalMass = 0.0f;
    }
}

void flushChemEnrich(void *vpkd, void *dst, const void *src) {
    PKD pkd = (PKD) vpkd;
    auto p1 = static_cast<stevFlush *>(dst);
    auto p2 = pkd->particles[static_cast<const PARTICLE *>(src)];

    if (p2.is_gas()) {
        const auto &sph = p2.sph();

        p1->ReceivedMom = sph.ReceivedMom;
        p1->fReceivedMass = sph.fReceivedMass;
        p1->fReceivedE = sph.fReceivedE;

        p1->ElemMass = sph.ElemMass;
        p1->fMetalMass = sph.fMetalMass;
    }
}

void combChemEnrich(void *vpkd, void *dst, const void *src) {
    PKD pkd = (PKD) vpkd;
    auto p1 = pkd->particles[static_cast<PARTICLE *>(dst)];
    auto p2 = static_cast<const stevFlush *>(src);

    if (p1.is_gas()) {
        auto &sph = p1.sph();

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
   in smf are in Msol. Conversely, star.fInitialMass and p.mass() are in code units.
   Furthermore, all times are in code units.
*/
void smChemEnrich(PARTICLE *pIn, float fBall, int nSmooth, NN *nnList, SMF *smf) {
    PKD pkd = smf->pkd;
    auto p = pkd->particles[pIn];
    auto &star = p.star();

    const float fInitialTime = star.fLastEnrichTime;
    const float fInitialMass = star.fLastEnrichMass;
    const int iInitialMass = star.iLastEnrichMass;

    const float fFinalTime = (float)smf->dTime - star.fTimer;
    const float fFinalMass = stevInverseLifetimeFunction(pkd, star, fFinalTime);
    const int iFinalMass = stevGetIMFMassIndex(pkd->StelEvolData->afInitialMass,
                           STEV_INTERP_N_MASS, fFinalMass, iInitialMass);

    star.fLastEnrichTime = fFinalTime;
    star.fLastEnrichMass = fFinalMass;
    star.iLastEnrichMass = iFinalMass + 1;

    /* Note: The parameter star.[AGB,CCSN,Lifetime].oZ contains the index of the
       interpolation's lower metallicity array multiplied by the number of mass bins.
       Since this multiplication must always be made, it is done once and for all in
       the function stevStarParticleInit. */

    TinyVector<float,ELEMENT_COUNT> ElemMass{0.0f};
    float fMetalMass{0.0f};

    const float fTransMass = smf->dCCSNMinMass;
    if (fFinalMass > fTransMass) {
        /* CCSN only */
        stevComputeMassToEject(pkd->StelEvolData->afCCSNYield,
                               pkd->StelEvolData->afCCSNMetalYield,
                               pkd->StelEvolData->afCCSNEjectedMass,
                               pkd->StelEvolData->afInitialMass,
                               pkd->StelEvolData->afIMFLogWeight,
                               pkd->StelEvolData->fDeltaLogMass,
                               STEV_CCSN_N_METALLICITY, STEV_INTERP_N_MASS, ELEMENT_COUNT,
                               star.ElemAbun.data(), star.fMetalAbun, iFinalMass, iInitialMass,
                               fFinalMass, fInitialMass, star.CCSN.oZ, star.CCSN.fDeltaZ,
                               ElemMass.data(), &fMetalMass);
    }
    else if (fInitialMass < fTransMass) {
        /* AGB only */
        stevComputeMassToEject(pkd->StelEvolData->afAGBYield,
                               pkd->StelEvolData->afAGBMetalYield,
                               pkd->StelEvolData->afAGBEjectedMass,
                               pkd->StelEvolData->afInitialMass,
                               pkd->StelEvolData->afIMFLogWeight,
                               pkd->StelEvolData->fDeltaLogMass,
                               STEV_AGB_N_METALLICITY, STEV_INTERP_N_MASS, ELEMENT_COUNT,
                               star.ElemAbun.data(), star.fMetalAbun, iFinalMass, iInitialMass,
                               fFinalMass, fInitialMass, star.AGB.oZ, star.AGB.fDeltaZ,
                               ElemMass.data(), &fMetalMass);
    }
    else {
        /* Mixed CCSN and AGB */
        int iTransMass = stevGetIMFMassIndex(pkd->StelEvolData->afInitialMass,
                                             STEV_INTERP_N_MASS, fTransMass, iInitialMass);

        stevComputeMassToEject(pkd->StelEvolData->afCCSNYield,
                               pkd->StelEvolData->afCCSNMetalYield,
                               pkd->StelEvolData->afCCSNEjectedMass,
                               pkd->StelEvolData->afInitialMass,
                               pkd->StelEvolData->afIMFLogWeight,
                               pkd->StelEvolData->fDeltaLogMass,
                               STEV_CCSN_N_METALLICITY, STEV_INTERP_N_MASS, ELEMENT_COUNT,
                               star.ElemAbun.data(), star.fMetalAbun, iTransMass, iInitialMass,
                               fTransMass, fInitialMass, star.CCSN.oZ, star.CCSN.fDeltaZ,
                               ElemMass.data(), &fMetalMass);

        ++iTransMass;
        stevComputeMassToEject(pkd->StelEvolData->afAGBYield,
                               pkd->StelEvolData->afAGBMetalYield,
                               pkd->StelEvolData->afAGBEjectedMass,
                               pkd->StelEvolData->afInitialMass,
                               pkd->StelEvolData->afIMFLogWeight,
                               pkd->StelEvolData->fDeltaLogMass,
                               STEV_AGB_N_METALLICITY, STEV_INTERP_N_MASS, ELEMENT_COUNT,
                               star.ElemAbun.data(), star.fMetalAbun, iFinalMass, iTransMass,
                               fFinalMass, fTransMass, star.AGB.oZ, star.AGB.fDeltaZ,
                               ElemMass.data(), &fMetalMass);
    }

    const float fNumSNIa = (*pkd->StelEvolData->fcnNumSNIa)(smf, star, fInitialTime, fFinalTime);
    for (auto i = 0; i < ELEMENT_COUNT; ++i) {
        ElemMass[i] += fNumSNIa * pkd->StelEvolData->afSNIaEjectedMass[i];
    }
    ElemMass *= star.fInitialMass;
    fMetalMass += fNumSNIa * pkd->StelEvolData->fSNIaEjectedMetalMass;
    fMetalMass *= star.fInitialMass;

    const double dTotalMass = (double)ElemMass[ELEMENT_H] + ElemMass[ELEMENT_He] +
                              fMetalMass;
    star.fNextEnrichTime = stevComputeNextEnrichTime(smf->dTime, star.fInitialMass,
                           dTotalMass, fFinalTime - fInitialTime);
    p.set_mass(p.mass() - dTotalMass);
    assert(p.mass() > 0.0f);

    const double dScaleFactorInv = 1.0 / csmTime2Exp(pkd->csm, smf->dTime);
    const double dScaleFactorInvSq = dScaleFactorInv * dScaleFactorInv;
    const auto &StarVel = p.velocity();

    const double dStarDeltaEkin = 0.5 * dTotalMass * blitz::dot(StarVel,StarVel);
    const double dWindEkin = smf->dWindSpecificEkin * dTotalMass;
    const double dStarEjEnergy = dStarDeltaEkin * dScaleFactorInvSq + dWindEkin;

    const double dWeight = 1.0 / nSmooth;
    const double dDeltaMass = dWeight * dTotalMass;
    const double dDeltaE = dWeight * dStarEjEnergy;

    ElemMass *= dWeight;
    fMetalMass *= dWeight;

    for (auto i = 0; i < nSmooth; ++i) {
        auto q = pkd->particles[nnList[i].pPart];
        assert(q.is_gas());
        auto &qsph = q.sph();

        qsph.ReceivedMom += dDeltaMass * StarVel * dScaleFactorInv;
        qsph.fReceivedMass += dDeltaMass;
        qsph.fReceivedE += dDeltaE;

        qsph.ElemMass += ElemMass;
        qsph.fMetalMass += fMetalMass;
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

    STEV_RAWDATA *RawData = (STEV_RAWDATA *) malloc(sizeof(STEV_RAWDATA));
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

    RawData->pfMetallicity = (float *) malloc(RawData->nZs * sizeof(float));
    assert(RawData->pfMetallicity != NULL);
    RawData->pfInitialMass = (float *) malloc(RawData->nMasses * sizeof(float));
    assert(RawData->pfInitialMass != NULL);
    RawData->pfYield = (float *) malloc(RawData->nZs * RawData->nElems * RawData->nMasses *
                                        sizeof(float));
    assert(RawData->pfYield != NULL);
    RawData->pfMetalYield = (float *) malloc(RawData->nZs * RawData->nMasses * sizeof(float));
    assert(RawData->pfMetalYield != NULL);
    RawData->pfEjectedMass = (float *) malloc(RawData->nZs * RawData->nMasses * sizeof(float));
    assert(RawData->pfEjectedMass != NULL);

    status = H5Dread(dsMetal, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     RawData->pfMetallicity);
    assert(status >= 0);
    status = H5Dclose(dsMetal);
    assert(status >= 0);

    status = H5Dread(dsMass, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     RawData->pfInitialMass);
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
        status = H5Dread(dsYield, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                         RawData->pfYield + iOffset);
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
        status = H5Dread(dsMetalYield, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                         RawData->pfMetalYield + iOffset);
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

    STEV_RAWDATA *RawData = (STEV_RAWDATA *) malloc(sizeof(STEV_RAWDATA));
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

    RawData->pfEjectedMass = (float *) malloc(RawData->nElems * sizeof(float));
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

STEV_RAWDATA *stevReadLifetimeTable(char *pszPath) {
    hid_t fileID, dataspace;
    herr_t status;

    STEV_RAWDATA *RawData = (STEV_RAWDATA *) malloc(sizeof(STEV_RAWDATA));
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

    RawData->pfMetallicity = (float *) malloc(RawData->nZs * sizeof(float));
    assert(RawData->pfMetallicity != NULL);
    RawData->pfInitialMass = (float *) malloc(RawData->nMasses * sizeof(float));
    assert(RawData->pfInitialMass != NULL);
    RawData->pfLifetime = (float *) malloc(RawData->nZs * RawData->nMasses * sizeof(float));
    assert(RawData->pfLifetime != NULL);

    status = H5Dread(dsMetal, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     RawData->pfMetallicity);
    assert(status >= 0);
    status = H5Dclose(dsMetal);
    assert(status >= 0);

    status = H5Dread(dsMass, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     RawData->pfInitialMass);
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
    status = H5Dread(dsLifetime, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                     RawData->pfLifetime);
    assert(status >= 0);
    status = H5Dclose(dsLifetime);
    assert(status >= 0);

    status = H5Fclose(fileID);
    assert(status >= 0);

    return RawData;
}

void stevFreeTable(STEV_RAWDATA *RawData) {
    free(RawData->pfMetallicity);
    free(RawData->pfInitialMass);
    free(RawData->pfYield);
    free(RawData->pfMetalYield);
    free(RawData->pfEjectedMass);
    free(RawData);
}

void stevFreeSNIaTable(STEV_RAWDATA *RawData) {
    free(RawData->pfMetalYield);
    free(RawData->pfEjectedMass);
    free(RawData);
}

void stevFreeLifetimeTable(STEV_RAWDATA *RawData) {
    free(RawData->pfMetallicity);
    free(RawData->pfInitialMass);
    free(RawData->pfLifetime);
    free(RawData);
}

float stevExponentialNumSNIa(SMF *smf, meshless::STAR &star, float fInitialTime, float fFinalTime) {
    if (fFinalTime <= star.fSNIaOnsetTime) {
        return 0.0f;
    }
    else if (fInitialTime < star.fSNIaOnsetTime) {
        fInitialTime = star.fSNIaOnsetTime;
    }

    return (float)smf->dSNIaNorm *
           (expf(-(fInitialTime / (float)smf->dSNIaScale)) -
            expf(-(fFinalTime / (float)smf->dSNIaScale)));
}

float stevPowerlawNumSNIa(SMF *smf, meshless::STAR &star, float fInitialTime, float fFinalTime) {
    if (fFinalTime <= star.fSNIaOnsetTime) {
        return 0.0f;
    }
    else if (fInitialTime < star.fSNIaOnsetTime) {
        fInitialTime = star.fSNIaOnsetTime;
    }

    return (float)smf->dSNIaNorm *
           (powf(fFinalTime, (float)smf->dSNIaScale + 1.0f) -
            powf(fInitialTime, (float)smf->dSNIaScale + 1.0f));
}

#endif  /* STELLAR_EVOLUTION */
