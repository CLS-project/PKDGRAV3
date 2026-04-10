#ifndef CORE_MESHLESS_H
#define CORE_MESHLESS_H
#include "core/datafields.h"
#include "blitz/array.h"
#include "chemistry.h"

namespace decorators {
typedef double myreal;

#define BH_ACCRETOR_FIELDS(X) X(iPid, int) X(iIndex, int)
MAKE_STRUCTURE(BH_ACCRETOR, BH_ACCRETOR_FIELDS)

// These fields are required for meshless finite volume (MFV),
// in addition to the meshless finate mass (MFM) fields (FIELDS).
#define MFV_FIELDS(X) X(lastDrDotFrho, blitz::TinyVector<myreal, 3>) X(drDotFrho, blitz::TinyVector<myreal, 3>)
MAKE_STRUCTURE(MFV, MFV_FIELDS)

#ifdef ENTROPY_SWITCH
    #define ENTROPY_FIELDS(X) X(S, double) X(lastS, double) X(maxEkin, double)
#else
    #define ENTROPY_FIELDS(X)
#endif
#ifdef STELLAR_EVOLUTION
    #define STELLAR_EVOLUTION_FIELDS(X) (ReceivedMom, blitz::TinyVector<float, 3>), (fReceivedMass, float), (fReceivedE, float)
#else
    #define STELLAR_EVOLUTION_FIELDS(X)
#endif

#define MFM_FIELDS(X)                                                                                                                                          \
    X(B, blitz::TinyVector<double, 6>)       /* IA: B matrix to 'easily' reconstruct faces'. Reminder: it is symmetric*/                                       \
    X(Ncond, myreal)                         /* IA: Condition number for pathological configurations */                                                        \
    X(gradRho, blitz::TinyVector<myreal, 3>) /* IA: Gradient of density */                                                                                     \
    X(gradVx, blitz::TinyVector<myreal, 3>)  /* IA: Gradient of x-velocity */                                                                                  \
    X(gradVy, blitz::TinyVector<myreal, 3>)  /* IA: Gradient of y-velocity */                                                                                  \
    X(gradVz, blitz::TinyVector<myreal, 3>)  /* IA: Gradient of z-velocity */                                                                                  \
    X(gradP, blitz::TinyVector<myreal, 3>)   /* IA: Gradient of pressure */                                                                                    \
    X(lastUpdateTime, myreal)                /* IA: last time this particle's primitve variables were updated */                                               \
    X(lastAcc, blitz::TinyVector<myreal, 3>) /**/                                                                                                              \
    X(lastMom, blitz::TinyVector<myreal, 3>) /**/                                                                                                              \
    X(lastE, myreal)                         /**/                                                                                                              \
    X(lastUint, myreal)                      /* TODO: Maybe there is a more intelligent way to avoid saving these... */                                        \
    X(lastHubble, myreal)                    /**/                                                                                                              \
    X(c, float)                              /* sound speed */                                                                                                 \
    X(lastMass, float)                       /**/                                                                                                              \
    X(omega, double)                         /* IA: normalization factor (Eq 7 Hopkins 2015) at the particle position */                                       \
    X(spheom, double)                        /**/                                                                                                              \
    X(Frho, myreal)                          /* IA: Fluxes */                                                                                                  \
    X(Fmom, blitz::TinyVector<myreal, 3>)    /**/                                                                                                              \
    X(Fene, myreal)                          /**/                                                                                                              \
    X(mom, blitz::TinyVector<double, 3>)     /* IA: Conserved variables */                                                                                     \
    X(E, double)                             /**/                                                                                                              \
    X(Uint, double)          /* IA: Internal energy, which is evolved in parallel and used for updating the pressure if we are in a cold flow  */              \
    ENTROPY_FIELDS(X)        /**/                                                                                                                              \
    X(P, double)             /* IA: Primitive variables */                                                                                                     \
    /*X(fLastBall, float)*/  /* IA: fBall from the last iteration. Used for the bisection algorithm */                                                         \
    /*X(nLastNeighs, int)*/  /* IA: Number of neighbors correspoding to that fBall */                                                                          \
    /*X(uNewRung, uint8_t)*/ /* IA: TODO temporarly */                                                                                                         \
    X(SFR, float)            /* IA: Star formation rate */                                                                                                     \
    X(ElemMass, blitz::TinyVector<float, ELEMENT_COUNT>) /*  */                                                                                                \
    X(fMetalMass, float)                                 /**/                                                                                                  \
    STELLAR_EVOLUTION_FIELDS(X)                          /**/                                                                                                  \
    X(fAccFBEnergy, float)                               /**/                                                                                                  \
    /* This could ideally be stored in a temporal buffer (pLite), but it has some limitations as it has to be shared when doing mdlAcquire.*/                  \
    X(BHAccretor, BH_ACCRETOR) /**/                                                                                                                            \
    X(uWake, uint8_t)
MAKE_STRUCTURE(MFM, MFM_FIELDS)

#define STAR_METALS_FIELDS(X) X(oZ, int) X(fDeltaZ, float)
MAKE_STRUCTURE(STAR_METALS, STAR_METALS_FIELDS)

#ifdef STELLAR_EVOLUTION
    #define STAR_STELLAR_EVOLUTION_FIELDS(X)                                                                                                                   \
        X(ElemAbun, blitz::TinyVector<float, ELEMENT_COUNT>) /* Formation abundances */                                                                        \
        X(fMetalAbun, float)                                 /* Formation metallicity */                                                                       \
        X(fInitialMass, float)                                                                                                                                 \
        X(fLastEnrichTime, float)                                                                                                                              \
        X(fLastEnrichMass, float)                                                                                                                              \
        X(iLastEnrichMass, int) X(fNextEnrichTime, float) X(CCSN, STAR_METALS) X(AGB, STAR_METALS) X(Lifetime, STAR_METALS) X(fSNIaOnsetTime, float)
#else
    #define STAR_STELLAR_EVOLUTION_FIELDS(X)
#endif
#define STAR_FIELDS(X)                                                                                                                                         \
    X(omega, double)                                                                                                                                           \
    STAR_STELLAR_EVOLUTION_FIELDS(X)                                                                                                                           \
    X(bCCSNFBDone, int)                                                                                                                                        \
    X(bSNIaFBDone, int)                                                                                                                                        \
    X(fSNEfficiency, float)                                                                                                                                    \
    X(fTimer, float)
MAKE_STRUCTURE(STAR, STAR_FIELDS)

#define GAS_PIN_FIELDS(X) X(iPid, int) X(iIndex, int)
MAKE_STRUCTURE(GAS_PIN, GAS_PIN_FIELDS)

#define BLACKHOLE_FIELDS(X)                                                                                                                                    \
    X(omega, double)                                                                                                                                           \
    X(dInternalMass, double)                                                                                                                                   \
    X(lastUpdateTime, double)                                                                                                                                  \
    X(dAccretionRate, double)                                                                                                                                  \
    X(dEddingtonRatio, double)                                                                                                                                 \
    X(dFeedbackRate, double)                                                                                                                                   \
    X(dAccEnergy, double)                                                                                                                                      \
    X(fTimer, float) /* Time of formation */                                                                                                                   \
    X(GasPin, GAS_PIN)                                                                                                                                         \
    X(bForceReposition, bool)
MAKE_STRUCTURE(BLACKHOLE, BLACKHOLE_FIELDS)

struct FIELDS {
    unsigned char zero[0];
};
static_assert(sizeof(FIELDS) == 0);

} // namespace decorators

#endif
