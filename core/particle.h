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

#include <vector>
#include <set>
#include <algorithm>
#include <type_traits>
#include <boost/iterator/iterator_facade.hpp>
#include "blitz/array.h"
#include "datastore.h"
#include "io/fio.h"
#include "chemistry.h"

#define PKD_MAX_CLASSES 256

#define INTEGER_FACTOR 0x80000000u
#define pkdDblToIntPos(pkd,d) (int32_t)((d)*INTEGER_FACTOR)
#define pkdIntPosToDbl(pkd,pos) ((pos)*(1.0/INTEGER_FACTOR))

#define IORDERBITS 43
#define IORDERMAX ((((uint64_t) 1)<<IORDERBITS)-1)

#define IRUNGBITS 6
#define IRUNGMAX ((1<<IRUNGBITS)-1)

/* Regular particle with order and all the goodies */
typedef struct particle {
uint64_t  uRung      :  IRUNGBITS;
    uint64_t  bMarked    :  1;
uint64_t  uNewRung   :  IRUNGBITS;  /* Optional with bNewKDK + bMemUnordered */
    uint64_t  iClass     :  8;          /* Optional with bMemUnordered */
uint64_t  iOrder     :  IORDERBITS; /* Optional with bMemUnordered */
} PARTICLE;
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
    PARTCLASS(float fMass,float fSoft,int iMat,FIO_SPECIES eSpecies)
        : fMass(fMass), fSoft(fSoft), iMat(iMat), eSpecies(eSpecies) {}
};
static_assert(std::is_trivial<PARTCLASS>());

typedef double myreal;

struct SPHFIELDS {
    char *pNeighborList; /* pointer to nearest neighbor list - compressed */
    double vPred[3];

    float c;        /* sound speed */
#ifndef OPTIM_REMOVE_UNUSED
    float u;            /* thermal energy */
    float uPred;    /* predicted thermal energy */
    float uDot;
    float divv;
    float BalsaraSwitch;    /* Balsara viscosity reduction */
    float fMetals;      /* mass fraction in metals, a.k.a, Z - tipsy output variable */

    /* diffusion */
    float diff;
    float fMetalsPred;
    float fMetalsDot;
#endif //OPTIM_REMOVE_UNUSED

    /* IA: B matrix to 'easily' reconstruct faces'. Reminder: it is symmetric */
    double B[6];

    /* IA: Condition number for pathological configurations */
    myreal Ncond;

    /* IA: Gradients */
    myreal gradRho[3];
    myreal gradVx[3];
    myreal gradVy[3];
    myreal gradVz[3];
    myreal gradP[3];

    /* IA: last time this particle's primitve variables were updated */
    myreal lastUpdateTime;
    myreal lastAcc[3];
    myreal lastMom[3];
    myreal lastE;
    myreal lastUint;
    myreal lastHubble; // TODO: Maybe there is a more intelligent way to avoid saving this...
#ifndef USE_MFM
    myreal lastDrDotFrho[3];
#endif
    float lastMass;

    /* IA: normalization factor (Eq 7 Hopkins 2015) at the particle position */
    double omega;

    /* IA: Fluxes */
    myreal Frho;
    myreal Fmom[3];
    myreal Fene;

#ifndef USE_MFM
    double drDotFrho[3];
#endif
    /* IA: Conserved variables */
    double mom[3];
    double E;
    /* IA: Internal energy, which is evolved in parallel and used for updating the pressure if we are in a cold flow */
    double Uint;

#ifdef ENTROPY_SWITCH
    double S;
    double lastS;
    double maxEkin;
#endif

    /* IA: Primitive variables */
    double P;

    /* IA: fBall from the last iteration. Used for the bisection algorithm */
    //float fLastBall;
    /* IA: Number of neighbors correspoding to that fBall */
    //int nLastNeighs;

    /* IA: TODO temporarly */
    //uint8_t uNewRung;

#ifdef STAR_FORMATION
    myreal SFR;
#endif

    float afElemMass[ELEMENT_COUNT];

#ifdef COOLING
    myreal lastCooling;
    float cooling_dudt;
#endif

#ifdef HAVE_METALLICITY
    float fMetalMass;
#endif

#ifdef STELLAR_EVOLUTION
    float afReceivedMom[3];
    float fReceivedMass;
    float fReceivedE;
#endif


#ifdef FEEDBACK
    float fAccFBEnergy;
#endif


    uint8_t uWake;

};

struct NEWSPHFIELDS {
    float Omega;        /* Correction factor */
    float divv;         /* Divergence of v */
    float u;            /* Thermodynamical variable, can be T, A(s) or u */
    float uDot;         /* Derivative of the thermodynamical variable */
    float cs;           /* Sound speed */
    float P;            /* Pressure */
    float oldRho;       /* Rho corresponding to where u is at */
    float expImb2;      /* exp(-imbalance^2) */
    float T;            /* Temperature */
};

struct STARFIELDS {
    double omega;
#ifdef STELLAR_EVOLUTION
    float afElemAbun[ELEMENT_COUNT]; /* Formation abundances */
    float fMetalAbun;            /* Formation metallicity */
    float fInitialMass;
    float fLastEnrichTime;
    float fLastEnrichMass;
    int iLastEnrichMass;
    float fNextEnrichTime;
    struct {
        int oZ;
        float fDeltaZ;
    } CCSN, AGB, Lifetime;
    float fSNIaOnsetTime;
#endif

    float fTimer;  /* Time of formation */
    float fSNEfficiency;
    int hasExploded; /* Has exploded as a supernova? */
};

struct BHFIELDS {
    PARTICLE *pLowPot;
    double omega;
    double dInternalMass;
    double newPos[3];
    double lastUpdateTime;
    double dAccretionRate;
    double dEddingtonRatio;
    double dFeedbackRate;
    double dAccEnergy;
    float fTimer;    /* Time of formation */
};

#ifdef OPTIM_UNION_EXTRAFIELDS
union EXTRAFIELDS {
    SPHFIELDS sph;
    STARFIELDS star;
    BHFIELDS bh;
};
#endif

//! \brief Enumerates all of the optional fields in a PARTICLE
//!
//! The PARTICLE structure can be minimally only four or eight bytes.
//! The fields listed below can be added based on the memory model.
enum class PKD_FIELD {
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
    oRelaxation,
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
class particleStore : public dataStore<PARTICLE,PKD_FIELD> {
protected:
    friend class Particle;
    bool bIntegerPosition = false;
    bool bNoParticleOrder = false;
    std::vector<PARTCLASS> ParticleClasses;
    float fSoftFix = -1.0;
    float fSoftFac = 1.0;
    float fSoftMax = HUGE_VALF;
public:
    void PhysicalSoft(double dSoftMax,double dFac,int bSoftMaxMul) {
        fSoftFac = dFac;
        fSoftMax = bSoftMaxMul ? HUGE_VALF : dSoftMax;
    }
    void SetSoft(double dSoft) {
        fSoftFix = dSoft;
    }
    blitz::TinyVector<double,3> position( PARTICLE *p ) const {
        if (bIntegerPosition) return get<int32_t,3>(p,PKD_FIELD::oPosition) * 1.0 / INTEGER_FACTOR;
        else return get<double,3>(p,PKD_FIELD::oPosition);
    }
    auto velocity( const PARTICLE *p )     const {return get<float,3>(p,PKD_FIELD::oVelocity);}
    auto acceleration( const PARTICLE *p ) const {return get<float,3>(p,PKD_FIELD::oAcceleration);}
    float mass( PARTICLE *p ) const {
        if (present(PKD_FIELD::oMass)) {
            return get<float>(p,PKD_FIELD::oMass)[0];
        }
        else if (bNoParticleOrder) return ParticleClasses[0].fMass;
        else return ParticleClasses[p->iClass].fMass;
    }
    float soft0( PARTICLE *p ) const {
        if (present(PKD_FIELD::oSoft)) {
            return get<float>(p,PKD_FIELD::oSoft)[0];
        }
        else if (bNoParticleOrder) return ParticleClasses[0].fSoft;
        else return ParticleClasses[p->iClass].fSoft;
    }
    float fixedsoft() const { return fSoftFix; }
    float soft( PARTICLE *p ) const {
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

    int32_t group(const PARTICLE *p ) const {
        if (bNoParticleOrder) return reinterpret_cast<const UPARTICLE *>(p)->iGroup;
        return get<int32_t>(p,PKD_FIELD::oGroup)[0];
    }
    void set_group( PARTICLE *p, uint32_t gid ) {
        if (bNoParticleOrder) reinterpret_cast<UPARTICLE *>(p)->iGroup = gid;
        else if (present(PKD_FIELD::oGroup)) get<int32_t>(p,PKD_FIELD::oGroup)[0] = gid;
    }

    int32_t global_gid(const PARTICLE *p ) const {
        return get<int64_t>(p,PKD_FIELD::oGlobalGid)[0];
    }
    void set_global_gid( PARTICLE *p, uint64_t gid ) {
        if (present(PKD_FIELD::oGlobalGid)) get<int64_t>(p,PKD_FIELD::oGlobalGid)[0] = gid;
    }

    auto density( const PARTICLE *p ) {
        return get<float>(p,PKD_FIELD::oDensity)[0];
    }
    void set_density( PARTICLE *p, float fDensity ) {
        if (present(PKD_FIELD::oDensity)) get<float>(p,PKD_FIELD::oDensity)[0] = fDensity;
    }
    auto ball( PARTICLE *p ) {
        return get<float>(p,PKD_FIELD::oBall)[0];
    }
    void set_ball(PARTICLE *p, float fBall) {
        if (present(PKD_FIELD::oBall)) get<float>(p,PKD_FIELD::oBall)[0] = fBall;
    }
    auto potential( PARTICLE *p ) const {
        return present(PKD_FIELD::oPotential) ? get<float>(p,PKD_FIELD::oPotential) : nullptr;
    }
    auto RungDest( PARTICLE *p ) const {
        return get<uint16_t>(p,PKD_FIELD::oRungDest);
    }
    auto ParticleID( PARTICLE *p ) const {
        return get<uint64_t>(p,PKD_FIELD::oParticleID);
    }
    /* Sph variables */
    auto sph( PARTICLE *p ) const {
#if defined(OPTIM_UNION_EXTRAFIELDS) && defined(DEBUG_UNION_EXTRAFIELDS)
        assert( species(p)==FIO_SPECIES_SPH);
#endif
        return get<SPHFIELDS>(p,PKD_FIELD::oSph);
    }
    auto sph( const PARTICLE *p ) const {
#if defined(OPTIM_UNION_EXTRAFIELDS) && defined(DEBUG_UNION_EXTRAFIELDS)
        assert( species(p)==FIO_SPECIES_SPH);
#endif
        return get<SPHFIELDS>(p,PKD_FIELD::oSph);
    }
    /* NewSph variables */
    auto newsph( PARTICLE *p ) const {
        return get<NEWSPHFIELDS>(p,PKD_FIELD::oNewSph);
    }
    auto newsph( const PARTICLE *p ) const {
        return get<NEWSPHFIELDS>(p,PKD_FIELD::oNewSph);
    }

    auto star( PARTICLE *p ) {
#ifdef OPTIM_UNION_EXTRAFIELDS
#ifdef DEBUG_UNION_EXTRAFIELDS
        assert( species(p)==FIO_SPECIES_STAR);
#endif //DEBUG
        return get<STARFIELDS>(p,PKD_FIELD::oSph);
#else
        return get<STARFIELDS>(p,PKD_FIELD::oStar);
#endif
    }
    auto star( const PARTICLE *p ) {
#ifdef OPTIM_UNION_EXTRAFIELDS
#ifdef DEBUG_UNION_EXTRAFIELDS
        assert( species(p)==FIO_SPECIES_STAR);
#endif //DEBUG
        return get<STARFIELDS>(p,PKD_FIELD::oSph);
#else
        return get<STARFIELDS>(p,PKD_FIELD::oStar);
#endif
    }
    auto BH( PARTICLE *p ) {
#ifdef OPTIM_UNION_EXTRAFIELDS
#ifdef DEBUG_UNION_EXTRAFIELDS
        assert( dpecies(p)==FIO_SPECIES_BH);
#endif //DEBUG
        return get<BHFIELDS>(p,PKD_FIELD::oSph);
#else
        return get<BHFIELDS>(p,PKD_FIELD::oBH);
#endif
    }
    auto BH( const PARTICLE *p ) {
#ifdef OPTIM_UNION_EXTRAFIELDS
#ifdef DEBUG_UNION_EXTRAFIELDS
        assert( dpecies(p)==FIO_SPECIES_BH);
#endif //DEBUG
        return get<BHFIELDS>(p,PKD_FIELD::oSph);
#else
        return get<BHFIELDS>(p,PKD_FIELD::oBH);
#endif
    }
    auto vPred( PARTICLE *p ) {
        return get<SPHFIELDS>(p,PKD_FIELD::oSph)->vPred;
    }
    auto pNeighborList( PARTICLE *p ) {
        return &get<SPHFIELDS>(p,PKD_FIELD::oSph)->pNeighborList;
    }
    auto Timer( PARTICLE *p ) {
        return &get<STARFIELDS>(p,PKD_FIELD::oStar)->fTimer;
    }
    auto isDeleted(PARTICLE *p) {
        return (species(p) == FIO_SPECIES_LAST);
    }
    auto isNew(PARTICLE *p) {
        return (p->iOrder == IORDERMAX);
    }
    bool isMarked(PARTICLE *p) {
        return p->bMarked;
    }
    uint8_t rung(PARTICLE *p) { return p->uRung; }
    bool isRungRange(PARTICLE *p,uint8_t uRungLo,uint8_t uRungHi) {
        return ((p->uRung >= uRungLo)&&(p->uRung <= uRungHi));
    }
public:
    class Particle {
    protected:
        PARTICLE *p;
        particleStore &store;
        template<typename T>
        T *get(PKD_FIELD f) const {return store.get<T>(p,f);}
        template<typename T,int n>
        blitz::TinyVector<T,n> get(PKD_FIELD f) const {
            return blitz::TinyVector<T,n>(store.get<T>(p,f));
        }
        bool have(PKD_FIELD f) const {return store.present(f);}
    public:
        Particle(particleStore &store,int i) : p(store.Element(i)), store(store) {}
        Particle(particleStore &store,void *p,int i=0) : p(store.Element(p,i)), store(store) {}
        operator PARTICLE *() {return p;}
        operator const PARTICLE *() const {return p;}
    public:
        static double IntPosToDbl(int32_t pos) {return pos*1.0/INTEGER_FACTOR;}
    public:
        using coord = blitz::TinyVector<double,3>;
        using icoord= blitz::TinyVector<int32_t,3>;
        bool have_position()     const {return have(PKD_FIELD::oPosition);}
        bool have_velocity()     const {return have(PKD_FIELD::oVelocity);}
        bool have_acceleration() const {return have(PKD_FIELD::oAcceleration);}
        bool have_sph()          const {return have(PKD_FIELD::oSph);}
        bool have_newsph()       const {return have(PKD_FIELD::oNewSph);}
        bool have_star()         const {return have(PKD_FIELD::oStar);}
        coord position() const {
            if (store.bIntegerPosition) return get<int32_t,3>(PKD_FIELD::oPosition) * 1.0 / INTEGER_FACTOR;
            else return get<double,3>(PKD_FIELD::oPosition);
        }
        double position(int d) const {
            if (store.bIntegerPosition) return IntPosToDbl(get<int32_t>(PKD_FIELD::oPosition)[d]);
            else return get<double>(PKD_FIELD::oPosition)[d];
        }
        auto velocity()     const {return get<float,3>(PKD_FIELD::oVelocity);}
        auto acceleration() const {return get<float,3>(PKD_FIELD::oAcceleration);}
        auto mass()         const {return store.mass(p);}
        auto soft0()        const {return store.soft0(p);}
        auto soft()         const {return store.soft(p);}
        auto species()      const {return store.species(p);}
        auto group()        const {return store.group(p);}
        auto density()      const {return store.density(p);}
        auto ball()         const {return store.ball(p);}
        auto potential()    const {return store.potential(p);}
        auto RungDest()     const { return store.RungDest(p); }
        auto ParticleID()   const { return store.ParticleID(p); }
        auto sph()          const { return store.sph(p); }
        auto newsph()       const { return store.newsph(p); }
        auto star()         const { return store.star(p); }
        auto BH()           const { return store.BH(p); }
        auto vPred()        const { return store.vPred(p); }
        auto pNeighborList()const { return store.pNeighborList(p); }
        auto Timer()        const { return store.Timer(p); }
        auto isDeleted()    const { return store.isDeleted(p); }
        auto isNew()        const { return store.isNew(p); }
        auto isMarked()     const { return store.isMarked(p); }
        auto rung()         const { return store.rung(p); }
        auto isRungRange(uint8_t uRungLo,uint8_t uRungHi)  const { return store.isRungRange(p,uRungLo,uRungHi); }

        void set_group( uint32_t gid ) { store.set_group(p,gid); }
        void set_density(float fDensity) {store.set_density(p,fDensity);}
        void set_ball(float fBall) {store.set_ball(p,fBall);}

    };
    Particle operator[](int i) {
        return Particle(*this,i);
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
    auto ParticleBase() { return Base(); }
    auto const ParticleBase() const { return Base(); }
    auto ParticleGet(void *pBase, int i) { return Element(pBase,i); }
    auto const ParticleGet(void *pBase, int i) const { return Element(pBase,i); }
    auto ParticleGet(int i) { return ParticleGet(Base(),i); }
    auto const ParticleGet(int i) const { return ParticleGet(Base(),i); }

protected:
    template <class Value,class Store>
    class particle_iter
        : public boost::iterator_facade<
          particle_iter<Value,Store>
        , Value
        , boost::iterators::random_access_traversal_tag
        , Particle
          > {
    public:
        particle_iter() : p(0), store(0) {}

        explicit particle_iter(Value *p,Store *store) : p(p),store(store) {}
        explicit particle_iter(int i,Store *store) : p(store->Element(i)),store(store) {}

        template <class OtherValue,class OtherStore>
        particle_iter(particle_iter<OtherValue,OtherStore> const &other)
            : p(other.p),store(other.store) {}

    private:
        friend class boost::iterator_core_access;
        template <class,class> friend class particle_iter;

        template <class OtherValue,class OtherStore>
        bool equal(particle_iter<OtherValue,OtherStore> const &other) const {
            return this->p == other.p;
        }

        void increment()
        { p = store->ParticleGet(p,1); }

        void decrement()
        { p = store->ParticleGet(p,-1); }

        void advance(int i)
        { p = store->ParticleGet(p,i); }

        template <class OtherValue,class OtherStore>
        int distance_to(particle_iter<OtherValue,OtherStore> other) const
        { return (reinterpret_cast<char *>(other.p) - reinterpret_cast<char *>(p)) / store->ParticleSize(); }

        auto dereference() const
        { return Particle(*store,p); }

        Value *p;
        Store *store;
    };
public:
    typedef particle_iter<PARTICLE,particleStore> particle_iterator;
    typedef particle_iter<PARTICLE const,particleStore const> particle_const_iterator;

    auto begin() { return particle_iterator(ParticleBase(),this); }
    auto begin() const { return particle_const_iterator(ParticleBase(),this); }
    auto end() { return particle_iterator(ParticleGet(Local()),this); }
    auto end() const { return particle_const_iterator(ParticleGet(Local()),this); }

public:
    void setClass( float fMass, float fSoft, int iMat, FIO_SPECIES eSpecies, PARTICLE *p ) {
        if ( present(PKD_FIELD::oMass) ) {
            auto pMass = get<float>(p,PKD_FIELD::oMass);
            *pMass = fMass;
            fMass = 0.0;
        }
        if ( present(PKD_FIELD::oSoft) ) {
            auto pSoft = get<float>(p,PKD_FIELD::oSoft);
            *pSoft = fSoft;
            fSoft = 0.0;
        }
        /* NOTE: The above can both be true, in which case a "zero" class is recorded */
        /* NOTE: Species is always part of the class table, so there will be at least one class per species */
        PARTCLASS newClass(fMass,fSoft,iMat,eSpecies);

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
                p = ParticleGet(i);
                assert( p->iClass < ParticleClasses.size() );
                p->iClass = map[p->iClass];
            }
        }

        /* Finally, set the new class table */
        ParticleClasses.clear();
        ParticleClasses.insert(ParticleClasses.end(),pClass,pClass+n);
    }

    void clearClasses() {ParticleClasses.clear();}

};

#endif
