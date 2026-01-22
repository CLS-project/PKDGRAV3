#ifndef STELLAREVOLUTION_H
#define STELLAREVOLUTION_H
#ifdef STELLAR_EVOLUTION
    #include "TraversePST.h"

    #define STEV_CCSN_N_METALLICITY 5
    #define STEV_CCSN_N_MASS 11
    #define STEV_AGB_N_METALLICITY 3
    #define STEV_AGB_N_MASS 23
    #define STEV_LIFETIME_N_METALLICITY 6
    #define STEV_LIFETIME_N_MASS 30

    #define STEV_INTERP_N_MASS 200

    #define STEV_MIN_LOG_METALLICITY -20

/*
 * ---------------------
 * STRUCTURE DEFINITIONS
 * ---------------------
 */

struct StellarEvolutionData {
    /* Pointer to function that gives the number of SNIa per Msol in [fInitialTime, fFinalTime] */
    float (*fcnNumSNIa)(double dSNIaNorm, double dSNIaScale, decorators::STAR& star, float fInitialTime, float fFinalTime);

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
};
using STEV_DATA = StellarEvolutionData;

struct StellarEvolutionRawData {
    int nZs;
    int nElems;
    int nMasses;
    float* pfMetallicity;
    float* pfInitialMass;
    float* pfYield;
    float* pfMetalYield;
    float* pfEjectedMass;
    float* pfLifetime;
};
using STEV_RAWDATA = struct StellarEvolutionRawData;

class ServiceStellarEvolutionInit : public TraversePST {
  public:
    static constexpr std::string_view name = "StellarEvolutionInit";
    struct input {
        char achSNIaDTDType[32];
        double dTime;
        double dSNIaMaxMass;
        double dCCSNMinMass;
        double dCCSNMaxMass;
        STEV_DATA StelEvolData;
    };
    using output = void;
    explicit ServiceStellarEvolutionInit(PST pst) : TraversePST(pst, name, sizeof(input)) {}

  protected:
    auto Service(PST pst, std::span<std::byte> in, std::span<std::byte> out) -> int override;
};

/*
 * --------------
 * MAIN FUNCTIONS
 * --------------
 */

void pkdAddStellarEjecta(PKD, particleStore::ParticleReference&, decorators::MFM&, const double);

/*
 * ----------------
 * HELPER FUNCTIONS
 * ----------------
 */
namespace stev {
auto ReadTable(char*) -> STEV_RAWDATA*;
auto ReadSNIaTable(char*) -> STEV_RAWDATA*;
auto ReadLifetimeTable(char*) -> STEV_RAWDATA*;
void FreeTable(STEV_RAWDATA*);
void FreeSNIaTable(STEV_RAWDATA*);
void FreeLifetimeTable(STEV_RAWDATA*);
auto ExponentialNumSNIa(double, double, decorators::STAR&, float, float) -> float;
auto PowerlawNumSNIa(double, double, decorators::STAR&, float, float) -> float;

void InterpToIMFSampling(STEV_DATA* const Data, STEV_RAWDATA* const CCSN, STEV_RAWDATA* const AGB, const float fCCSNMinMass);
void StarParticleInit(PKD pkd, decorators::STAR& star, const double dSNIaMaxMass, const double dCCSNMinMass, const double dCCSNMaxMass);
} // namespace stev

#endif /* STELLAR_EVOLUTION */
#endif /* STELLAREVOLUTION_H */