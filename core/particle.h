/*  This file is part of PKDGRAV3 (http://www.pkdgrav.org/).
 *  Copyright (c) 2001-2018 Joachim Stadel & Douglas Potter
 *
 *  PKDGRAV3 is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  PKDGRAV3 is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with PKDGRAV3.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CORE_PARTICLE_H
#define CORE_PARTICLE_H
#include "pkd_config.h"
#include <limits>
#include <vector>
#include <set>
#include <algorithm>
#include <functional>
#include <type_traits>
#include <iterator>
#include <boost/hana/define_struct.hpp>
#include "blitz/array.h"
#include "datastore.h"
#include "integerize.h"
#include "element.h"
#include "bound.h"
#include "io/fio.h"
#include "chemistry.h"

#define PKD_MAX_CLASSES 256

#define INTEGER_FACTOR 0x80000000u
#define pkdDblToIntPos(pkd,d) (int32_t)((d)*INTEGER_FACTOR)
#define pkdIntPosToDbl(pkd,pos) ((pos)*(1.0/INTEGER_FACTOR))

#ifdef NN_FLAG_IN_PARTICLE
    #define IORDERBITS 42
#else
    #define IORDERBITS 43
#endif
#define IORDERMAX ((((uint64_t) 1)<<IORDERBITS)-1)

#define IRUNGBITS 6
#define IRUNGMAX ((1<<IRUNGBITS)-1)

/* Regular particle with order and all the goodies */
struct PARTICLE {
    uint64_t  uRung      :  IRUNGBITS;
    uint64_t  bMarked    :  1;
    uint64_t  uNewRung   :  IRUNGBITS;  /* Optional with bNewKDK + bMemUnordered */
    uint64_t  iClass     :  8;          /* Optional with bMemUnordered */
    uint64_t  iOrder     :  IORDERBITS; /* Optional with bMemUnordered */
#ifdef NN_FLAG_IN_PARTICLE
    uint64_t bNNflag : 1;           /* Neighbor of Neighbor of active flag */
#endif
};
static_assert(sizeof(PARTICLE)==sizeof(uint64_t));

/* Abbreviated particle header with group id */
#define IGROUPBITS (32-IRUNGBITS-1)
#define IGROUPMAX ((1<<IGROUPBITS)-1)

typedef struct uparticle {
    uint32_t  uRung      :  IRUNGBITS;
    uint32_t  bMarked    :  1;
    uint32_t  iGroup     :  IGROUPBITS;
} UPARTICLE;
static_assert(sizeof(UPARTICLE)==sizeof(uint32_t));

struct PARTCLASS {
    float       fMass;    /* Particle mass */
    float       fSoft;    /* Current softening */
    int         iMat;     /* newSPH material id */
    FIO_SPECIES eSpecies; /* Species: dark, star, etc. */

    bool operator <(const PARTCLASS &b) const {
        if ( fMass < b.fMass ) return true;
        else if ( fMass > b.fMass ) return false;
        else if ( fSoft < b.fSoft ) return true;
        else if ( fSoft > b.fSoft ) return false;
        else if ( iMat  < b.iMat  ) return true;
        else if ( iMat  > b.iMat  ) return false;
        else return eSpecies < b.eSpecies;
    }
    bool operator==(const PARTCLASS &b) const {
        return fMass==b.fMass && fSoft==b.fSoft && iMat==b.iMat && eSpecies==b.eSpecies;
    }
    PARTCLASS() = default;
    PARTCLASS(FIO_SPECIES eSpecies,float fMass,float fSoft,int iMat)
        : fMass(fMass), fSoft(fSoft), iMat(iMat), eSpecies(eSpecies) {}
};
static_assert(std::is_trivial<PARTCLASS>());

namespace meshless {
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

struct FIELDS {
    BOOST_HANA_DEFINE_STRUCT(
        FIELDS,

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
#ifndef USE_MFM
        (blitz::TinyVector<myreal,3>, lastDrDotFrho),
#endif
        (float, c),        /* sound speed */

        (float, lastMass),

        /* IA: normalization factor (Eq 7 Hopkins 2015) at the particle position */
        (double, omega),

        /* IA: Fluxes */
        (myreal, Frho),
        (blitz::TinyVector<myreal,3>,  Fmom),
        (myreal, Fene),

#ifndef USE_MFM
        (blitz::TinyVector<double,3>, drDotFrho),
#endif
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

}; // FIELDS

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

#ifdef OPTIM_UNION_EXTRAFIELDS
union EXTRAFIELDS {
    FIELDS sph;
    STAR star;
    BLACKHOLE bh;
};
#endif

} // meshless

namespace sph {
namespace hana = boost::hana;

struct FIELDS {
    BOOST_HANA_DEFINE_STRUCT(
        FIELDS,
        (float, Omega),        // Correction factor
        (float, divv),         // Divergence of v
        (float, u),            // Thermodynamical variable, can be T, A(s) or u
        (float, uDot),         // Derivative of the thermodynamical variable
        (float, cs),           // Sound speed
        (float, P),            // Pressure
        (float, oldRho),       // Rho corresponding to where u is at
        (float, expImb2),      // exp(-imbalance^2)
        (float, T),            // Temperature
        (blitz::TinyVector<float,3>, vpred)  // predicted velocities
    );
}; // FIELDS

} // sph

struct VELSMOOTH {
    BOOST_HANA_DEFINE_STRUCT(
        VELSMOOTH,
        (blitz::TinyVector<float,3>, vmean),
        (float, divv),
        (float, veldisp2)
    );
};

//! \brief Enumerates all of the optional fields in a PARTICLE
//!
//! The PARTICLE structure can be minimally only four or eight bytes.
//! The fields listed below can be added based on the memory model.
enum PKD_FIELD {
    oPosition,
    oVelocity, /* Three vel_t */
    oAcceleration, /* Three float */
    oPotential, /* One float */
    oGroup, /* One int32 */
    oMass, /* One float */
    oSoft, /* One float */
    oDensity, /* One float */
    oBall, /* One float */
    oSph, /* Sph structure */
    oNewSph, /* NewSph structure */
    oStar, /* Star structure */
    oBH, /* BH structure */
    oVelSmooth,
    oRungDest, /* Destination processor for each rung */
    oParticleID,
    oGlobalGid, /* global group id, uint64 */

    MAX_FIELD
};

//! \brief Array of particles (main storage).
//!
//! The particleStore provides functions to access elements of the particle array,
//! and the individual fields of which there may be a varying number.
class particleStore : public dataStore<PARTICLE,PKD_FIELD>, protected Integerize {
protected:
    friend class Particle;
    bool bIntegerPosition = false;
    bool bNoParticleOrder = false;
    std::vector<PARTCLASS> ParticleClasses;
    float fSoftFix = -1.0;
    float fSoftFac = 1.0;
    float fSoftMax = HUGE_VALF;
    typedef uint8_t rung_type;
    rung_type uMinRungActive = 0;
    rung_type uMaxRungActive = std::numeric_limits<rung_type>::max();
public:
    using coord = blitz::TinyVector<double,3>;
    using icoord= blitz::TinyVector<int32_t,3>;
    void PhysicalSoft(double dSoftMax,double dFac,int bSoftMaxMul) {
        fSoftFac = dFac;
        fSoftMax = bSoftMaxMul ? HUGE_VALF : dSoftMax;
    }
    void SetSoft(double dSoft) {
        fSoftFix = dSoft;
    }
    void ActiveRung(int iRungMin, int iRungMax) {
        uMinRungActive = iRungMin;
        uMaxRungActive = iRungMax;
    }
    void ActiveRung(int iRungMin=0) {
        ActiveRung(iRungMin,std::numeric_limits<rung_type>::max());
    }

    // Naming convention. If a field is always part of the particle, then we can return
    // a reference and we use the name of the field (e.g. density). If the field has to
    // be converted or can come from different places then we use dont't return a reference,
    // and one needs to use the set version to set (e.g. set_position).
    coord position( PARTICLE *p ) const {
        if (bIntegerPosition) return convert(get<int32_t[3]>(p,PKD_FIELD::oPosition));
        else return get<double[3]>(p,PKD_FIELD::oPosition);
    }
    double position( const PARTICLE *p, int d ) const {
        if (bIntegerPosition) return convert(get<int32_t[3]>(p,PKD_FIELD::oPosition)[d]);
        else return get<double[3]>(p,PKD_FIELD::oPosition)[d];
    }
    void set_position( PARTICLE *p, coord r ) const {
        if (bIntegerPosition) get<int32_t[3]>(p,PKD_FIELD::oPosition) = convert(r);
        else get<double[3]>(p,PKD_FIELD::oPosition) = r;
    }
    float mass( const PARTICLE *p ) const {
        if (present(PKD_FIELD::oMass)) {
            return get<float>(p,PKD_FIELD::oMass);
        }
        else if (bNoParticleOrder) return ParticleClasses[0].fMass;
        else return ParticleClasses[p->iClass].fMass;
    }
    void set_mass( PARTICLE *p, float mass ) const {
        assert(present(PKD_FIELD::oMass));
        get<float>(p,PKD_FIELD::oMass) = mass;
    }
    uint32_t group(const PARTICLE *p ) const {
        if (bNoParticleOrder) return reinterpret_cast<const UPARTICLE *>(p)->iGroup;
        return get<int32_t>(p,PKD_FIELD::oGroup);
    }
    void set_group( PARTICLE *p, uint32_t gid ) const {
        if (bNoParticleOrder) reinterpret_cast<UPARTICLE *>(p)->iGroup = gid;
        else if (present(PKD_FIELD::oGroup)) get<int32_t>(p,PKD_FIELD::oGroup) = gid;
    }
    template<typename T>
    auto &raw_position( PARTICLE *p ) const {
        return get<T[3]>(p,PKD_FIELD::oPosition);
    }
    template<typename T>
    auto &raw_position( PARTICLE *p, int i ) const {
        return get<T[3]>(p,PKD_FIELD::oPosition)[i];
    }
    const auto &velocity(     const PARTICLE *p ) const {return get<float[3]>(p,PKD_FIELD::oVelocity);}
    auto       &velocity(           PARTICLE *p ) const {return get<float[3]>(p,PKD_FIELD::oVelocity);}
    const auto &acceleration( const PARTICLE *p ) const {return get<float[3]>(p,PKD_FIELD::oAcceleration);}
    auto       &acceleration(       PARTICLE *p ) const {return get<float[3]>(p,PKD_FIELD::oAcceleration);}
    auto       &potential(          PARTICLE *p ) const {return get<float>(p,PKD_FIELD::oPotential);}
    const auto &potential(    const PARTICLE *p ) const {return get<float>(p,PKD_FIELD::oPotential);}
    auto       &density(            PARTICLE *p ) const {return get<float>(p,PKD_FIELD::oDensity);}
    const auto &density(      const PARTICLE *p ) const {return get<float>(p,PKD_FIELD::oDensity);}
    auto       &ball(               PARTICLE *p ) const {return get<float>(p,PKD_FIELD::oBall);}
    const auto &ball(         const PARTICLE *p ) const {return get<float>(p,PKD_FIELD::oBall);}
    auto       &global_gid(         PARTICLE *p ) const {return get<int64_t>(p,PKD_FIELD::oGlobalGid);}
    const auto &global_gid(   const PARTICLE *p ) const {return get<int64_t>(p,PKD_FIELD::oGlobalGid);}
    auto       &ParticleID(         PARTICLE *p ) const {return get<uint64_t>(p,PKD_FIELD::oParticleID);}
    const auto &ParticleID(   const PARTICLE *p ) const {return get<uint64_t>(p,PKD_FIELD::oParticleID);}
    auto       &RungDest(           PARTICLE *p ) const {return get<uint16_t[8]>(p,PKD_FIELD::oRungDest);}
    const auto &RungDest(     const PARTICLE *p ) const {return get<uint16_t[8]>(p,PKD_FIELD::oRungDest);}
    auto       &VelSmooth(          PARTICLE *p ) const {return get<VELSMOOTH>(p,PKD_FIELD::oVelSmooth);}
    const auto &VelSmooth(    const PARTICLE *p ) const {return get<VELSMOOTH>(p,PKD_FIELD::oVelSmooth);}

    float soft0( const PARTICLE *p ) const {
        if (present(PKD_FIELD::oSoft)) {
            return get<float>(p,PKD_FIELD::oSoft);
        }
        else if (bNoParticleOrder) return ParticleClasses[0].fSoft;
        else return ParticleClasses[p->iClass].fSoft;
    }
    float fixedsoft() const { return fSoftFix; }
    float soft( const PARTICLE *p ) const {
        float fSoft = fSoftFac * (fSoftFix >= 0.0 ? fSoftFix : soft0(p));
        if ( fSoft > fSoftMax ) fSoft = fSoftMax;
        return fSoft;
    }
    FIO_SPECIES species( PARTICLE *p ) const {
        if (bNoParticleOrder) return ParticleClasses[0].eSpecies;
        else return ParticleClasses[p->iClass].eSpecies;
    }
    int iMat( PARTICLE *p) const {
        if (bNoParticleOrder) return ParticleClasses[0].iMat;
        else return ParticleClasses[p->iClass].iMat;
    }
    std::set<int> getMaterials() const {
        std::set<int> materials;
        for (auto p : ParticleClasses) {
            materials.insert(p.iMat);
        }
        return materials;
    }

    void set_density( PARTICLE *p, float fDensity ) const {
        if (present(PKD_FIELD::oDensity)) get<float>(p,PKD_FIELD::oDensity) = fDensity;
    }
    void set_ball(PARTICLE *p, float fBall) const {
        if (present(PKD_FIELD::oBall)) get<float>(p,PKD_FIELD::oBall) = fBall;
    }
    /* Sph variables */
    auto &sph( PARTICLE *p ) const {
#if defined(OPTIM_UNION_EXTRAFIELDS) && defined(DEBUG_UNION_EXTRAFIELDS)
        assert( species(p)==FIO_SPECIES_SPH);
#endif
        return get<meshless::FIELDS>(p,PKD_FIELD::oSph);
    }
    const auto &sph( const PARTICLE *p ) const {
#if defined(OPTIM_UNION_EXTRAFIELDS) && defined(DEBUG_UNION_EXTRAFIELDS)
        assert( species(p)==FIO_SPECIES_SPH);
#endif
        return get<meshless::FIELDS>(p,PKD_FIELD::oSph);
    }
    /* NewSph variables */
    auto &newsph( PARTICLE *p ) const {
        return get<sph::FIELDS>(p,PKD_FIELD::oNewSph);
    }
    const auto &newsph( const PARTICLE *p ) const {
        return get<sph::FIELDS>(p,PKD_FIELD::oNewSph);
    }

    auto &star( PARTICLE *p ) const {
#ifdef OPTIM_UNION_EXTRAFIELDS
#ifdef DEBUG_UNION_EXTRAFIELDS
        assert( species(p)==FIO_SPECIES_STAR);
#endif //DEBUG
        return get<meshless::STAR>(p,PKD_FIELD::oSph);
#else
        return get<meshless::STAR>(p,PKD_FIELD::oStar);
#endif
    }
    const auto &star( const PARTICLE *p ) const {
#ifdef OPTIM_UNION_EXTRAFIELDS
#ifdef DEBUG_UNION_EXTRAFIELDS
        assert( species(p)==FIO_SPECIES_STAR);
#endif //DEBUG
        return get<meshless::STAR>(p,PKD_FIELD::oSph);
#else
        return get<meshless::STAR>(p,PKD_FIELD::oStar);
#endif
    }
    auto &BH( PARTICLE *p ) const {
#ifdef OPTIM_UNION_EXTRAFIELDS
#ifdef DEBUG_UNION_EXTRAFIELDS
        assert( dpecies(p)==FIO_SPECIES_BH);
#endif //DEBUG
        return get<meshless::BLACKHOLE>(p,PKD_FIELD::oSph);
#else
        return get<meshless::BLACKHOLE>(p,PKD_FIELD::oBH);
#endif
    }
    const auto &BH( const PARTICLE *p ) const {
#ifdef OPTIM_UNION_EXTRAFIELDS
#ifdef DEBUG_UNION_EXTRAFIELDS
        assert( dpecies(p)==FIO_SPECIES_BH);
#endif //DEBUG
        return get<meshless::BLACKHOLE>(p,PKD_FIELD::oSph);
#else
        return get<meshless::BLACKHOLE>(p,PKD_FIELD::oBH);
#endif
    }
    auto Timer( PARTICLE *p ) const {
        return &get<meshless::STAR>(p,PKD_FIELD::oStar).fTimer;
    }
    auto is_deleted(PARTICLE *p) const {
        return (species(p) == FIO_SPECIES_UNKNOWN);
    }
    auto isNew(PARTICLE *p) const {
        return (p->iOrder == IORDERMAX);
    }
    bool isMarked(PARTICLE *p) const {
        return p->bMarked;
    }
    uint8_t rung(PARTICLE *p) const { return p->uRung; }
    bool is_rung_range(PARTICLE *p,uint8_t uRungLo,uint8_t uRungHi) const {
        return ((p->uRung >= uRungLo)&&(p->uRung <= uRungHi));
    }
    bool is_active(PARTICLE *p) const {
        return is_rung_range(p,uMinRungActive,uMaxRungActive);
    }
public:
    /// @brief Proxy object for a particle
    ///
    /// A particle in pkdgrav3 has a variable size which is tracked with a series of offsets.
    /// This class keeps track of the pointer to the particle, and a reference to the particle layout
    /// which is kept in the particle store.
    class Particle : public SingleElement<PARTICLE,particleStore> {
    public:
        using SingleElement::SingleElement;
    public:
        static double IntPosToDbl(int32_t pos) {return pos*1.0/INTEGER_FACTOR;}
    public:
        using coord = blitz::TinyVector<double,3>;
        using icoord= blitz::TinyVector<int32_t,3>;
        bool have_position()     const {return have(PKD_FIELD::oPosition);}
        bool have_velocity()     const {return have(PKD_FIELD::oVelocity);}
        bool have_acceleration() const {return have(PKD_FIELD::oAcceleration);}
        bool have_potential()    const {return have(PKD_FIELD::oPotential);}
        bool have_sph()          const {return have(PKD_FIELD::oSph);}
        bool have_newsph()       const {return have(PKD_FIELD::oNewSph);}
        bool have_star()         const {return have(PKD_FIELD::oStar);}
        bool have_ball()         const {return have(PKD_FIELD::oBall);}
        bool have_mass()         const {return have(PKD_FIELD::oMass);}
        bool have_group()        const {return have(PKD_FIELD::oGroup);}
        bool have_soft()         const {return have(PKD_FIELD::oSoft);}
        bool have_density()      const {return have(PKD_FIELD::oDensity);}
        bool have_bh()           const {return have(PKD_FIELD::oBH);}
        bool have_vel_smooth()   const {return have(PKD_FIELD::oVelSmooth);}
        bool have_rung_dest()    const {return have(PKD_FIELD::oRungDest);}
        bool have_particle_id()  const {return have(PKD_FIELD::oParticleID);}
        bool have_global_gid()   const {return have(PKD_FIELD::oGlobalGid);}

        auto position()          const {return store().position(p); }
        auto position(int d)     const {return store().position(p,d); }
        auto mass()              const {return store().mass(p);}
        auto group()             const {return store().group(p);}
        bool marked()            const {return p->bMarked;}
        uint8_t rung()           const {return p->uRung;}
        uint8_t new_rung()       const {return p->uNewRung;}
        uint8_t get_class()      const {return p->iClass;}
        uint64_t order()         const {return p->iOrder;}
#ifdef NN_FLAG_IN_PARTICLE
        bool NN_flag()           const {return p->bNNflag;}
        bool set_NN_flag(bool bFlag)    const {return (p->bNNflag = bFlag);}
#else
        bool NN_flag()           const {return false;}
        void set_NN_flag(bool bFlag)    const {}
#endif
        void set_position(coord r)      const {store().set_position(p,r);}
        void set_mass( float mass)      const {store().set_mass(p,mass);}
        void set_group( uint32_t gid )  const {store().set_group(p,gid); }
        bool set_marked(bool bMarked)   const {return (p->bMarked = bMarked);}
        auto set_rung(uint8_t uRung)    const {return (p->uRung = uRung);}
        auto set_new_rung(uint8_t uRung)const {return (p->uNewRung = uRung);}
        auto set_class(uint8_t uClass)  const {return (p->iClass = uClass);}
        void set_class( float fMass, float fSoft, int iMat, FIO_SPECIES eSpecies) {
            store().setClass(fMass,fSoft,iMat,eSpecies,p);
        }
        auto set_order(uint64_t uOrder) const {return (p->iOrder = uOrder);}

        template<typename T>
        auto &raw_position()        const {return store().raw_position<T>(p);}
        template<typename T>
        auto &raw_position(int i)   const {return store().raw_position<T>(p,i);}
        auto &velocity()            const {return store().velocity(p);}
        auto &acceleration()        const {return store().acceleration(p);}
        auto &potential()           const {return store().potential(p);}
        auto &density()             const {return store().density(p);}
        auto &ball()                const {return store().ball(p);}
        auto &global_gid()          const {return store().global_gid(p); }
        auto &ParticleID()          const {return store().ParticleID(p); }
        auto &RungDest()            const {return store().RungDest(p); }
        auto &VelSmooth()           const {return store().VelSmooth(p); }

        auto soft()         const {return store().soft(p);}
        auto soft0()        const {return store().soft0(p);}
        auto species()      const {return store().species(p);}
        auto imaterial()    const {return store().iMat(p);}
        auto &sph()         const { return store().sph(p); }
        auto &newsph()      const { return store().newsph(p); }
        auto &star()        const { return store().star(p); }
        auto &BH()          const { return store().BH(p); }
        auto Timer()        const { return store().Timer(p); }
        auto isNew()        const { return store().isNew(p); }
        auto isMarked()     const { return store().isMarked(p); }
        auto is_deleted()   const { return store().is_deleted(p); }
        auto is_new()       const { return store().isNew(p); }
        auto is_marked()    const { return store().isMarked(p); }
        auto is_active()    const { return store().is_active(p);}
        auto is_rung_range(uint8_t uRungLo,uint8_t uRungHi)  const { return store().is_rung_range(p,uRungLo,uRungHi); }
        void set_density(float fDensity) {store().set_density(p,fDensity);}
        void set_ball(float fBall) {store().set_ball(p,fBall);}

        auto is_dark()      const {return species() == FIO_SPECIES_DARK;}
        auto is_gas()       const {return species() == FIO_SPECIES_SPH;}
        auto is_star()      const {return species() == FIO_SPECIES_STAR;}
        auto is_bh()        const {return species() == FIO_SPECIES_BH;}
    };

    /// @brief A "pointer" reference to a "particle"
    ///
    /// This uses the default copy/assignment semantics (the internal pointers are copied).
    /// The Particle class uses direct assignment.
    using ParticlePointer = SinglePointer<Particle>;

    /// @brief A "reference" to a particle
    /// This allows a particle to be assigned and swapped
    using ParticleReference = SingleReference<Particle>;

    /// @brief Return a pointer proxy to a particle
    /// @param i Index of the particle
    /// @return pointer proxy
    auto operator[](int i) {
        return ParticleReference(*this,i);
    }
    auto operator[](PARTICLE *p) {
        return ParticleReference(*this,p);
    }
    const auto operator[](const PARTICLE *p) const {
        return ParticleReference(*this,p);
    }

    auto NewParticle() {
        assert(Local()<FreeStore());
        auto i = Local();
        SetLocal(Local()+1);
        auto p = ParticlePointer(*this,i);
        if (!bNoParticleOrder) p->set_order(IORDERMAX);
        return p;
    }

public:
    void initialize(bool bIntegerPosition,bool bNoParticleOrder) {
        this->bIntegerPosition = bIntegerPosition;
        this->bNoParticleOrder = bNoParticleOrder;
        dataStore<PARTICLE,PKD_FIELD>::initialize(bNoParticleOrder ? sizeof(UPARTICLE) : sizeof(PARTICLE));
        ParticleClasses.reserve(PKD_MAX_CLASSES);
    }
    auto integerized() const {return bIntegerPosition;}
    auto unordered() const {return bNoParticleOrder;}
    auto ParticleSize() const {return ElementSize(); }

public:
    typedef element_iterator<Particle> iterator;
    typedef element_iterator<Particle const> const_iterator;

    auto begin()        noexcept { return iterator(*this,Base()); }
    auto begin()  const noexcept { return const_iterator(*this,Base()); }
    auto cbegin() const noexcept { return const_iterator(*this,Base()); }
    auto end()          noexcept { return iterator(*this,Element(Local())); }
    auto end()    const noexcept { return const_iterator(*this,Element(Local())); }
    auto cend()   const noexcept { return const_iterator(*this,Element(Local())); }

public:
    void setClass( float fMass, float fSoft, int iMat, FIO_SPECIES eSpecies, PARTICLE *p ) {
        if ( present(PKD_FIELD::oMass) ) {
            get<float>(p,PKD_FIELD::oMass) = fMass;
            fMass = 0.0;
        }
        if ( present(PKD_FIELD::oSoft) ) {
            get<float>(p,PKD_FIELD::oSoft) = fSoft;
            fSoft = 0.0;
        }
        /* NOTE: The above can both be true, in which case a "zero" class is recorded */
        /* NOTE: Species is always part of the class table, so there will be at least one class per species */
        PARTCLASS newClass(eSpecies,fMass,fSoft,iMat);

        /* TODO: This is a linear search which is fine for a small number of classes */
        auto iclass = std::find(ParticleClasses.begin(),ParticleClasses.end(),newClass);
        if (iclass==ParticleClasses.end()) {
            assert( ParticleClasses.size() < PKD_MAX_CLASSES );
            p->iClass = ParticleClasses.size();
            ParticleClasses.emplace_back(newClass);
        }
        else p->iClass = std::distance(ParticleClasses.begin(),iclass);
        if (bNoParticleOrder) { assert(p->iClass==0); }
    }

    int getClasses( int nMax, PARTCLASS *pClass ) {
        std::copy(ParticleClasses.begin(),ParticleClasses.end(),pClass);
        return ParticleClasses.size();
    }

    void setClasses( int n, PARTCLASS *pClass, int bUpdate ) {
        uint8_t map[PKD_MAX_CLASSES];
        PARTICLE *p;

        if ( bUpdate && ParticleClasses.size() && !bNoParticleOrder) {
            /* Build a map from the old class to the new class */
            assert( n >= ParticleClasses.size() );
            for (auto i=0; i<ParticleClasses.size(); ++i) {
                auto jj = std::find(pClass,pClass+n,ParticleClasses[i]);
                map[i] = std::distance(pClass,jj);
            }

            /* Now update the class with the new value */
            for (auto i=0; i<Local(); ++i) {
                p = Element(i);
                assert( p->iClass < ParticleClasses.size() );
                p->iClass = map[p->iClass];
            }
        }

        /* Finally, set the new class table */
        ParticleClasses.clear();
        ParticleClasses.insert(ParticleClasses.end(),pClass,pClass+n);
    }

    void clearClasses() {ParticleClasses.clear();}

    /// @brief Calculate the bound for a range of particles
    /// @tparam T Either double or int32_t depending on the position type
    /// @param b beginning particle iterator
    /// @param e ending particle iterator
    /// @return bounding box of the particles
    template<typename T>
    typename bound_type<T>::type raw_bound(iterator b,iterator e) {
        blitz::TinyVector<T,3> vmin(std::numeric_limits<T>::max()),vmax(std::numeric_limits<T>::lowest());
        std::for_each(b,e,
        [&vmin,&vmax](auto &p) {
            auto v = p.template raw_position<T>();
            vmin = blitz::min(vmin,v);
            vmax = blitz::max(vmax,v);
        });
        return typename bound_type<T>::type(vmin,vmax);
    }

    Bound bound(iterator b,iterator e);
    Bound bound();

    int CountSelected();
    int SelActive(int setIfTrue=true, int clearIfFalse=true);
    int SelBlackholes(int setIfTrue=true, int clearIfFalse=true);
    int SelSpecies(uint64_t mSpecies, int setIfTrue=true, int clearIfFalse=true);
    int SelGroup(int iGroup, int setIfTrue=true, int clearIfFalse=true);
    int SelMass(double dMinMass, double dMaxMass, int setIfTrue=true, int clearIfFalse=true );
    int SelPhaseDensity(double dMinDensity, double dMaxDensity, int setIfTrue=true, int clearIfFalse=true );
    int SelById(uint64_t idStart, uint64_t idEnd, int setIfTrue=true, int clearIfFalse=true );
    int SelBox(blitz::TinyVector<double,3> dCenter, blitz::TinyVector<double,3> dSize, int setIfTrue=true, int clearIfFalse=true );
    int SelSphere(blitz::TinyVector<double,3> r, double dRadius, int setIfTrue=true, int clearIfFalse=true );
    int SelCylinder(blitz::TinyVector<double,3> dP1, blitz::TinyVector<double,3> dP2, double dRadius, int setIfTrue=true, int clearIfFalse=true );
};

#endif
