#ifdef STELLAR_EVOLUTION


#include <hdf5.h>

#include "stellarevolution.h"
#include "master.h"
#include "imf.h"
#include "hydro/hydro.h"


void MSR::SetStellarEvolutionParam() {
    const double dYrToTime = SECONDSPERYEAR / param.units.dSecUnit;

    if (strcmp(param.achSNIaDTDType, "exponential") == 0) {
        param.dSNIaNorm = param.dSNIaNumPerMass;
        param.dSNIaScale = param.dSNIaExpScale * dYrToTime;
    }
    else if (strcmp(param.achSNIaDTDType, "powerlaw") == 0) {
        param.dSNIaNorm = param.dSNIaNumPerMass /
                          (pow(param.dSNIaPLFinalTime * dYrToTime, param.dSNIaPLScale + 1.0) -
                           pow(param.dSNIaPLInitTime * dYrToTime, param.dSNIaPLScale + 1.0));
        param.dSNIaScale = param.dSNIaPLScale;
    }
    else {
        std::cerr << "ERROR: Undefined SNIa DTD type has been given in " <<
                  "achSNIaDTDType parameter: " << param.achSNIaDTDType << std::endl;
        Exit(1);
    }

    param.dWindSpecificEkin = 0.5 * pow(param.dStellarWindSpeed / param.units.dKmPerSecUnit, 2);
}


void MSR::StellarEvolutionInit(double dTime) {
    int i;
    char achPath[280];
    struct inStellarEvolutionInit in;
    STEV_RAWDATA *CCSNData, *AGBData, *SNIaData, *LifetimeData;

    /* Read the tables */
    sprintf(achPath, "%s/CCSN.hdf5", param.achStelEvolPath);
    CCSNData = stevReadTable(achPath);
    assert(CCSNData->nZs == STEV_CCSN_N_METALLICITY);
    assert(CCSNData->nMasses == STEV_CCSN_N_MASS);
    assert(CCSNData->nElems == ELEMENT_COUNT);

    sprintf(achPath, "%s/AGB.hdf5", param.achStelEvolPath);
    AGBData = stevReadTable(achPath);
    assert(AGBData->nZs == STEV_AGB_N_METALLICITY);
    assert(AGBData->nMasses == STEV_AGB_N_MASS);
    assert(AGBData->nElems == ELEMENT_COUNT);

    sprintf(achPath, "%s/SNIa.hdf5", param.achStelEvolPath);
    SNIaData = stevReadSNIaTable(achPath);
    assert(SNIaData->nElems == ELEMENT_COUNT);

    sprintf(achPath, "%s/Lifetimes.hdf5", param.achStelEvolPath);
    LifetimeData = stevReadLifetimeTable(achPath);
    assert(LifetimeData->nZs == STEV_LIFETIME_N_METALLICITY);
    assert(LifetimeData->nMasses == STEV_LIFETIME_N_MASS);

    /* NOTE: The lowest value of the initial mass array is set to the corresponding value
       of the Lifetime table that is being used, while its highest value to the maximum
       stellar mass for a CCSN to occur, as set in param.dCCSNMaxMass. The IMF is still
       normalized in the range [param.dIMFMinMass,param.dIMFMaxMass] */
    const double dMinMass = LifetimeData->pfInitialMass[0];
    const double dMaxMass = param.dCCSNMaxMass;
    double adInitialMass[STEV_INTERP_N_MASS], adIMF[STEV_INTERP_N_MASS];

    auto IMF = ChooseIMF(param.achIMFType, param.dIMFMinMass, param.dIMFMaxMass);
    IMF->MassWeightedSample(dMinMass, dMaxMass, STEV_INTERP_N_MASS, adInitialMass, adIMF);

    /* Store the data */
    for (i = 0; i < STEV_INTERP_N_MASS; ++i) {
        in.StelEvolData.afInitialMass[i] = adInitialMass[i];
        in.StelEvolData.afIMFLogWeight[i] = adIMF[i];
    }
    const double dDeltaLog = (log10(dMaxMass) - log10(dMinMass)) / (STEV_INTERP_N_MASS - 1);
    in.StelEvolData.fDeltaLogMass = dDeltaLog;

    /* Convert CCSN/AGB initial mass arrays to log for interpolation */
    for (i = 0; i < CCSNData->nMasses; ++i)
        CCSNData->pfInitialMass[i] = log10(CCSNData->pfInitialMass[i]);
    for (i = 0; i < AGBData->nMasses; ++i)
        AGBData->pfInitialMass[i] = log10(AGBData->pfInitialMass[i]);

    /* Interpolate yields and ejected masses to the initial mass array */
    stevInterpToIMFSampling(&in.StelEvolData, CCSNData, AGBData, param.dCCSNMinMass);

    for (i = 0; i < STEV_CCSN_N_METALLICITY; ++i) {
        if (CCSNData->pfMetallicity[i] > 0.0f)
            in.StelEvolData.afCCSNMetallicity[i] = log10(CCSNData->pfMetallicity[i]);
        else
            in.StelEvolData.afCCSNMetallicity[i] = STEV_MIN_LOG_METALLICITY;
    }

    for (i = 0; i < STEV_AGB_N_METALLICITY; ++i) {
        if (AGBData->pfMetallicity[i] > 0.0f)
            in.StelEvolData.afAGBMetallicity[i] = log10(AGBData->pfMetallicity[i]);
        else
            in.StelEvolData.afAGBMetallicity[i] = STEV_MIN_LOG_METALLICITY;
    }

    for (i = 0; i < ELEMENT_COUNT; ++i)
        in.StelEvolData.afSNIaEjectedMass[i] = SNIaData->pfEjectedMass[i];
    in.StelEvolData.fSNIaEjectedMetalMass = *SNIaData->pfMetalYield;

    for (i = 0; i < STEV_LIFETIME_N_METALLICITY; ++i) {
        if (LifetimeData->pfMetallicity[i] > 0.0f)
            in.StelEvolData.afLifetimeMetallicity[i] = log10(LifetimeData->pfMetallicity[i]);
        else
            in.StelEvolData.afLifetimeMetallicity[i] = STEV_MIN_LOG_METALLICITY;
    }
    for (i = 0; i < STEV_LIFETIME_N_MASS; ++i)
        in.StelEvolData.afLifetimeInitialMass[i] = log10(LifetimeData->pfInitialMass[i]);
    for (i = 0; i < STEV_LIFETIME_N_METALLICITY * STEV_LIFETIME_N_MASS; ++i) {
        in.StelEvolData.afLifetime[i] = log10(LifetimeData->pfLifetime[i] * SECONDSPERYEAR /
                                              param.units.dSecUnit);
    }

    strcpy(in.achSNIaDTDType, param.achSNIaDTDType);
    in.dTime = dTime;
    in.dSNIaMaxMass = param.dSNIaMaxMass;
    in.dCCSNMinMass = param.dCCSNMinMass;
    in.dCCSNMaxMass = param.dCCSNMaxMass;

    /* Send data and finish */
    pstStellarEvolutionInit(pst, &in, sizeof(struct inStellarEvolutionInit), NULL, 0);

    stevFreeTable(CCSNData);
    stevFreeTable(AGBData);
    stevFreeSNIaTable(SNIaData);
    stevFreeLifetimeTable(LifetimeData);
}


#ifdef __cplusplus
extern "C" {
#endif


int pstStellarEvolutionInit(PST pst, void *vin, int nIn, void *vout, int nOut) {
    mdlassert(pst->mdl, nIn == sizeof(struct inStellarEvolutionInit));
    if (pst->nLeaves > 1) {
        int rID = mdlReqService(pst->mdl, pst->idUpper, PST_STELLAREVOLUTIONINIT, vin, nIn);
        pstStellarEvolutionInit(pst->pstLower, vin, nIn, NULL, 0);
        mdlGetReply(pst->mdl, rID, NULL, NULL);
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

    for (int i = 0; i < pkd->nLocal; ++i) {
        PARTICLE *p = pkdParticle(pkd, i);

        if (pkdIsStar(pkd, p)) {
            STARFIELDS *pStar = pkdStar(pkd, p);

            if (pStar->fInitialMass < 0.0f)
                pStar->fInitialMass = pkdMass(pkd, p);

            if (pStar->fLastEnrichTime < 0.0f)
                pStar->fLastEnrichTime = 0.0f;

            if (pStar->fTimer <= 0.0f)
                pStar->fTimer = in->dTime;

            stevStarParticleInit(pkd, pStar, in->dSNIaMaxMass, in->dCCSNMinMass,
                                 in->dCCSNMaxMass);
        }
        else if (pkdIsGas(pkd, p)) {
            SPHFIELDS *pSph = pkdSph(pkd, p);
            for (int j = 0; j < 3; ++j)
                pSph->afReceivedMom[j] = 0.0f;
            pSph->fReceivedMass = 0.0f;
            pSph->fReceivedE = 0.0f;
        }
    }

    return 0;
}


void pkdAddStellarEjecta(PKD pkd, PARTICLE *p, SPHFIELDS *pSph, const double dConstGamma) {
    const double dOldEkin = (pSph->mom[0] * pSph->mom[0] + pSph->mom[1] * pSph->mom[1] +
                             pSph->mom[2] * pSph->mom[2]) / (2.0 * pkdMass(pkd, p));

    for (int i = 0; i < 3; ++i)
        pSph->mom[i] += pSph->afReceivedMom[i];
    *((float *) pkdField(p, pkd->oFieldOffset[oMass])) += pSph->fReceivedMass;
    pSph->E += pSph->fReceivedE;

    const double dNewEkin = (pSph->mom[0] * pSph->mom[0] + pSph->mom[1] * pSph->mom[1] +
                             pSph->mom[2] * pSph->mom[2]) / (2.0 * pkdMass(pkd, p));

    const double dDeltaUint = pSph->fReceivedE - (dNewEkin - dOldEkin);
    if (dDeltaUint > 0.0) {
        pSph->Uint += dDeltaUint;
#ifdef ENTROPY_SWITCH
        pSph->S += dDeltaUint * (dConstGamma - 1.0) *
                   pow(pkdDensity(pkd, p), 1.0 - dConstGamma);
#endif
    }

    for (int i = 0; i < 3; ++i)
        pSph->afReceivedMom[i] = 0.0f;
    pSph->fReceivedMass = 0.0f;
    pSph->fReceivedE = 0.0f;
}


void initChemEnrich(void *vpkd, void *vp) {
    int i;
    PKD pkd = (PKD) vpkd;
    PARTICLE *p = (PARTICLE *) vp;

    if (pkdIsGas(pkd, p)) {
        SPHFIELDS *pSph = pkdSph(pkd,p);

        for (i = 0; i < 3; ++i)
            pSph->afReceivedMom[i] = 0.0f;
        pSph->fReceivedMass = 0.0f;
        pSph->fReceivedE = 0.0f;

        for (i = 0; i < ELEMENT_COUNT; ++i)
            pSph->afElemMass[i] = 0.0f;
        pSph->fMetalMass = 0.0f;
    }
}


void combChemEnrich(void *vpkd, void *vp1, const void *vp2) {
    int i;
    PKD pkd = (PKD) vpkd;
    PARTICLE *p1 = (PARTICLE *) vp1;
    PARTICLE *p2 = (PARTICLE *) vp2;

    if (pkdIsGas(pkd, p1) && pkdIsGas(pkd, p2)) {
        SPHFIELDS *pSph1 = pkdSph(pkd, p1);
        SPHFIELDS *pSph2 = pkdSph(pkd, p2);

        for (i = 0; i < 3; ++i)
            pSph1->afReceivedMom[i] += pSph2->afReceivedMom[i];
        pSph1->fReceivedMass += pSph2->fReceivedMass;
        pSph1->fReceivedE += pSph2->fReceivedE;

        for (i = 0; i < ELEMENT_COUNT; ++i)
            pSph1->afElemMass[i] += pSph2->afElemMass[i];
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
   in smf are in Msol. Conversely, pStar->fInitialMass and pkdMass(pkd,p) are in code units.
   Furthermore, all times are in code units.
*/
void smChemEnrich(PARTICLE *p, float fBall, int nSmooth, NN *nnList, SMF *smf) {
    int i, j;
    PKD pkd = smf->pkd;
    STARFIELDS *pStar = pkdStar(pkd, p);

    const float fInitialTime = pStar->fLastEnrichTime;
    const float fInitialMass = pStar->fLastEnrichMass;
    const int iInitialMass = pStar->iLastEnrichMass;

    const float fFinalTime = (float)smf->dTime - pStar->fTimer;
    const float fFinalMass = stevInverseLifetimeFunction(pkd, pStar, fFinalTime);
    const int iFinalMass = stevGetIMFMassIndex(pkd->StelEvolData->afInitialMass,
                           STEV_INTERP_N_MASS, fFinalMass, iInitialMass);

    pStar->fLastEnrichTime = fFinalTime;
    pStar->fLastEnrichMass = fFinalMass;
    pStar->iLastEnrichMass = iFinalMass + 1;


    /* Note: The parameter pStar->[AGB,CCSN,Lifetime].oZ contains the index of the
       interpolation's lower metallicity array multiplied by the number of mass bins.
       Since this multiplication must always be made, it is done once and for all in
       the function stevStarParticleInit. */

    float afElemMass[ELEMENT_COUNT];
    float fMetalMass;
    for (i = 0; i < ELEMENT_COUNT; ++i)
        afElemMass[i] = 0.0f;
    fMetalMass = 0.0f;

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
                               pStar->afElemAbun, pStar->fMetalAbun, iFinalMass, iInitialMass,
                               fFinalMass, fInitialMass, pStar->CCSN.oZ, pStar->CCSN.fDeltaZ,
                               afElemMass, &fMetalMass);
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
                               pStar->afElemAbun, pStar->fMetalAbun, iFinalMass, iInitialMass,
                               fFinalMass, fInitialMass, pStar->AGB.oZ, pStar->AGB.fDeltaZ,
                               afElemMass, &fMetalMass);
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
                               pStar->afElemAbun, pStar->fMetalAbun, iTransMass, iInitialMass,
                               fTransMass, fInitialMass, pStar->CCSN.oZ, pStar->CCSN.fDeltaZ,
                               afElemMass, &fMetalMass);

        ++iTransMass;
        stevComputeMassToEject(pkd->StelEvolData->afAGBYield,
                               pkd->StelEvolData->afAGBMetalYield,
                               pkd->StelEvolData->afAGBEjectedMass,
                               pkd->StelEvolData->afInitialMass,
                               pkd->StelEvolData->afIMFLogWeight,
                               pkd->StelEvolData->fDeltaLogMass,
                               STEV_AGB_N_METALLICITY, STEV_INTERP_N_MASS, ELEMENT_COUNT,
                               pStar->afElemAbun, pStar->fMetalAbun, iFinalMass, iTransMass,
                               fFinalMass, fTransMass, pStar->AGB.oZ, pStar->AGB.fDeltaZ,
                               afElemMass, &fMetalMass);
    }

    const float fNumSNIa = (*pkd->StelEvolData->fcnNumSNIa)(smf, pStar, fInitialTime, fFinalTime);
    for (i = 0; i < ELEMENT_COUNT; ++i) {
        afElemMass[i] += fNumSNIa * pkd->StelEvolData->afSNIaEjectedMass[i];
        afElemMass[i] *= pStar->fInitialMass;
    }
    fMetalMass += fNumSNIa * pkd->StelEvolData->fSNIaEjectedMetalMass;
    fMetalMass *= pStar->fInitialMass;


    const float fTotalMass = afElemMass[ELEMENT_H] + afElemMass[ELEMENT_He] + fMetalMass;
    pStar->fNextEnrichTime = stevComputeNextEnrichTime(smf->dTime, pStar->fInitialMass,
                             fTotalMass, fFinalTime - fInitialTime);
    *((float *) pkdField(p, pkd->oFieldOffset[oMass])) -= fTotalMass;
    assert(pkdMass(pkd, p) > 0.0f);


    const float fExpFacInv = 1.0 / csmTime2Exp(pkd->csm, smf->dTime);
    const float fExpFacInvSq = fExpFacInv * fExpFacInv;
    const vel_t *const pStarVel = pkdVel(pkd, p);

    const float fStarDeltaEkin = 0.5f * fTotalMass * (pStarVel[0] * pStarVel[0] +
                                 pStarVel[1] * pStarVel[1] + pStarVel[2] * pStarVel[2]);
    const float fWindEkin = (float)smf->dWindSpecificEkin * fTotalMass;
    const float fStarEjEnergy = fStarDeltaEkin * fExpFacInvSq + fWindEkin;


    PARTICLE *q;
    float fWeights[nSmooth];
    float fNormFactor = 0.0f;
    for (i = 0; i < nSmooth; ++i) {
        q = nnList[i].pPart;
        if (q == p) continue;

        const double dRpq = sqrt(nnList[i].fDist2);
        fWeights[i] = cubicSplineKernel(dRpq, 0.5*fBall) / pkdDensity(pkd, q);
        fNormFactor += fWeights[i];
    }
    fNormFactor = 1.0f / fNormFactor;

    for (i = 0; i < nSmooth; ++i) {
        q = nnList[i].pPart;
        if (q == p) continue;

        fWeights[i] *= fNormFactor;
        const float fDeltaMass = fWeights[i] * fTotalMass;
        SPHFIELDS *qSph = pkdSph(pkd, q);

        for (j = 0; j < 3; ++j)
            qSph->afReceivedMom[j] += fDeltaMass * pStarVel[j] * fExpFacInv;
        qSph->fReceivedMass += fDeltaMass;
        qSph->fReceivedE += fWeights[i] * fStarEjEnergy;

        for (j = 0; j < ELEMENT_COUNT; ++j)
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

        sprintf(achTable, "%s/%s/%s", H5FIELD_TABLE, apchTableZNames[i], H5FIELD_YIELD);
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


float stevExponentialNumSNIa(SMF *smf, STARFIELDS *pStar, float fInitialTime, float fFinalTime) {
    if (fFinalTime <= pStar->fSNIaOnsetTime) {
        return 0.0f;
    }
    else if (fInitialTime < pStar->fSNIaOnsetTime) {
        fInitialTime = pStar->fSNIaOnsetTime;
    }

    return (float)smf->dSNIaNorm *
           (expf(-(fInitialTime / (float)smf->dSNIaScale)) -
            expf(-(fFinalTime / (float)smf->dSNIaScale)));
}


float stevPowerlawNumSNIa(SMF *smf, STARFIELDS *pStar, float fInitialTime, float fFinalTime) {
    if (fFinalTime <= pStar->fSNIaOnsetTime) {
        return 0.0f;
    }
    else if (fInitialTime < pStar->fSNIaOnsetTime) {
        fInitialTime = pStar->fSNIaOnsetTime;
    }

    return (float)smf->dSNIaNorm *
           (powf(fFinalTime, (float)smf->dSNIaScale + 1.0f) -
            powf(fInitialTime, (float)smf->dSNIaScale + 1.0f));
}


#ifdef __cplusplus
}
#endif


#endif  /* STELLAR_EVOLUTION */
