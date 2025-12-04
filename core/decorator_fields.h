#ifndef CORE_MESHLESS_H
#define CORE_MESHLESS_H
#include <boost/hana/define_struct.hpp>
#include "blitz/array.h"
#include "chemistry.h"

namespace decorators {
namespace hana = boost::hana;
typedef double myreal;

#ifdef BLACKHOLES
struct  BH_ACCRETOR {
    BOOST_HANA_DEFINE_STRUCT(
        BH_ACCRETOR,
        (int, iPid),
        (int, iIndex)
    );
};
#endif

// These fields are required for meshless finite volume (MFV),
// in addition to the meshless finate mass (MFM) fields (FIELDS).
struct MFV {
    BOOST_HANA_DEFINE_STRUCT(
        MFV,
        (blitz::TinyVector<myreal,3>, lastDrDotFrho),
        (blitz::TinyVector<myreal,3>, drDotFrho)
    );
};

struct MFM {
    BOOST_HANA_DEFINE_STRUCT(
        MFM,

        (blitz::TinyVector<double,6>, B), // IA: B matrix to 'easily' reconstruct faces'. Reminder: it is symmetric
        (myreal, Ncond),/* IA: Condition number for pathological configurations */

        /* IA: Gradients */
        (blitz::TinyVector<myreal,3>, gradRho),
        (blitz::TinyVector<myreal,3>, gradVx),
        (blitz::TinyVector<myreal,3>, gradVy),
        (blitz::TinyVector<myreal,3>, gradVz),
        (blitz::TinyVector<myreal,3>, gradP),

        /* IA: last time this particle's primitve variables were updated */
        (myreal, lastUpdateTime),
        (blitz::TinyVector<myreal,3>, lastAcc),
        (blitz::TinyVector<myreal,3>, lastMom),
        (myreal, lastE),
        (myreal, lastUint),
        (myreal, lastHubble), // TODO: Maybe there is a more intelligent way to avoid saving this...
        (float, c),        /* sound speed */

        (float, lastMass),

        /* IA: normalization factor (Eq 7 Hopkins 2015) at the particle position */
        (double, omega),
        (double, spheom),

        /* IA: Fluxes */
        (myreal, Frho),
        (blitz::TinyVector<myreal,3>,  Fmom),
        (myreal, Fene),

        /* IA: Conserved variables */
        (blitz::TinyVector<double,3>,  mom),
        (double, E),
        /* IA: Internal energy, which is evolved in parallel and used for updating the pressure if we are in a cold flow */
        (double, Uint),

#ifdef ENTROPY_SWITCH
        (double, S),
        (double, lastS),
        (double, maxEkin),
#endif

        /* IA: Primitive variables */
        (double, P),

        /* IA: fBall from the last iteration. Used for the bisection algorithm */
        //(float, fLastBall),
        /* IA: Number of neighbors correspoding to that fBall */
        //(int, nLastNeighs),

        /* IA: TODO temporarly */
        //(uint8_t, uNewRung),

#ifdef STAR_FORMATION
        (float, SFR),
#endif

        (blitz::TinyVector<float,ELEMENT_COUNT>, ElemMass),
#ifdef HAVE_METALLICITY
        (float, fMetalMass),
#endif

#ifdef STELLAR_EVOLUTION
        (blitz::TinyVector<float,3>, ReceivedMom),
        (float, fReceivedMass),
        (float, fReceivedE),
#endif

#if defined(FEEDBACK) || defined(BLACKHOLES)
        (float, fAccFBEnergy),
#endif

#ifdef BLACKHOLES
        // This could ideally be stored in a temporal buffer (pLite), but it has some
        // limitations as it has to be shared when doing mdlAcquire.
        (BH_ACCRETOR, BHAccretor),
#endif

        (uint8_t, uWake)
    );

}; // MFM

struct STAR_METALS {
    BOOST_HANA_DEFINE_STRUCT(
        STAR_METALS,
        (int, oZ),
        (float, fDeltaZ)
    );
}; // STAR_METALS

struct STAR {
    BOOST_HANA_DEFINE_STRUCT(
        STAR,
        (double, omega),
#ifdef STELLAR_EVOLUTION
        (blitz::TinyVector<float,ELEMENT_COUNT>, ElemAbun), /* Formation abundances */
        (float, fMetalAbun),            /* Formation metallicity */
        (float, fInitialMass),
        (float, fLastEnrichTime),
        (float, fLastEnrichMass),
        (int, iLastEnrichMass),
        (float, fNextEnrichTime),
        (STAR_METALS, CCSN),
        (STAR_METALS, AGB),
        (STAR_METALS, Lifetime),
        (float, fSNIaOnsetTime),
#endif
#ifdef FEEDBACK
        (int, bCCSNFBDone),
        (int, bSNIaFBDone),
        (float, fSNEfficiency),
#endif
        (float, fTimer)  /* Time of formation */
    );
}; // STAR

struct GAS_PIN {
    BOOST_HANA_DEFINE_STRUCT(
        GAS_PIN,
        (int, iPid),
        (int, iIndex)
    );
};

struct BLACKHOLE {
    BOOST_HANA_DEFINE_STRUCT(
        BLACKHOLE,
        (double, omega),
        (double, dInternalMass),
        (double, lastUpdateTime),
        (double, dAccretionRate),
        (double, dEddingtonRatio),
        (double, dFeedbackRate),
        (double, dAccEnergy),
        (float, fTimer),    /* Time of formation */
        (GAS_PIN, GasPin),
        (bool, bForceReposition)
    );
};

struct FIELDS { unsigned char zero[0]; };
static_assert(sizeof(FIELDS)==0);

} // meshless

#endif
