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

#ifndef CORE_TREENODE_H
#define CORE_TREENODE_H

#include <vector>
#include "bound.h"
#include "splitstore.h"
#include "integerize.h"
#include "particle.h"
#include "gravity/moments.h"
#include "SPH/SPHOptions.h"

//! \brief Enumerates all of the optional fields in a TREE node
//!
//! The fields listed below can be added based on the memory model.
enum class KDN_FIELD {
    oNodePosition,      // Three double or int32_t (if bIntegerPosition)
    oNodeVelocity,      // Three vel_t
    oNodeAcceleration,  // Three doubles
    oNodeMass,          // One float (may point into oNodeMom!)
    oNodeBOB,           // Ball or Box of Balls
    oNodeMom,           // an FMOMR
    oNodeBnd,
    oNodeSphBounds,     // Three Bounds
    oNodeVBnd,          // Velocity bounds
    oNodeNgas,
    oNodeNstar,
    oNodeNbh,

    MAX_FIELD
};

struct KDN {
    int pLower;          /* also serves as thread id for the LTT */
    int pUpper;          /* pUpper < 0 indicates no particles in tree! */
    uint32_t iLower;         /* Local lower node (or remote processor w/bRemote=1) */
    uint32_t iDepth     : 13;
    uint32_t bGroup     : 1;
    uint32_t iSplitDim  : 2;
    uint32_t uMinRung   : 6;
    uint32_t uMaxRung   : 6;
    uint32_t bTopTree   : 1; /* This is a top tree node: pLower,pUpper are node indexes */
    uint32_t bRemote    : 1; /* children are remote */
    uint32_t bHasMarked : 1;         /* flag if node has a marked particle */
    uint32_t bHasNNflag : 1;         /* flag if node has a NNflagged particle */
    float bMax;
    float fSoft2;
};

typedef struct sphBounds {
    struct minmaxBound {
        blitz::TinyVector<double,3> min;
        blitz::TinyVector<double,3> max;
    } A,B,BI;
} SPHBNDS;

class treeStore : public splitStore<KDN,KDN_FIELD>, protected Integerize {
protected:
    bool bIntegerPosition = false;
    particleStore *pstore;
    typedef float vel_t;
public:
    void initialize(particleStore *pstore,bool bIntegerPosition,int Lo, int Hi, uint64_t nData=0, void *pData=nullptr) {
        this->bIntegerPosition = bIntegerPosition;
        this->pstore = pstore;
        setStore(Lo,Hi,nData,pData);
    }
    static_assert(std::is_standard_layout<blitz::TinyVector<double,3>>());
    static_assert(std::is_standard_layout<blitz::TinyVector<int32_t,3>>());

    auto is_integer_position() {return bIntegerPosition;}
    auto position(const KDN *n ) const {
        if (bIntegerPosition) return convert(get<int32_t[3]>(n,KDN_FIELD::oNodePosition));
        else return get<double[3]>(n,KDN_FIELD::oNodePosition);
    }
    auto set_position(KDN *n, blitz::TinyVector<double,3> r) const {
        if (bIntegerPosition) get<int32_t[3]>(n,KDN_FIELD::oNodePosition) = convert(r);
        else get<double[3]>(n,KDN_FIELD::oNodePosition) = r;
        return r;
    }
    template<typename T>
    auto &raw_position( KDN *n ) const {
        return get<T[3]>(n,KDN_FIELD::oNodePosition);
    }

    auto bound(const KDN *n ) const {
        if (bIntegerPosition) return Bound(get<IntegerBound>(n,KDN_FIELD::oNodeBnd),*this);
        else return get<Bound>(n,KDN_FIELD::oNodeBnd);
    }
    auto set_bound(KDN *n, const Bound &bnd) const {
        if (bIntegerPosition) get<IntegerBound>(n,KDN_FIELD::oNodeBnd) = IntegerBound(bnd,*this);
        else get<Bound>(n,KDN_FIELD::oNodeBnd) = bnd;
        return bnd;
    }
    template<typename T>
    auto &raw_bound( KDN *n ) const {
        return get<T>(n,KDN_FIELD::oNodeBnd);
    }

    const auto &velocity(     const KDN *n ) const {return get<vel_t[3]>(n,KDN_FIELD::oNodeVelocity);}
    auto       &velocity(           KDN *n ) const {return get<vel_t[3]>(n,KDN_FIELD::oNodeVelocity);}
    const auto &acceleration( const KDN *n ) const {return get<float[3]>(n,KDN_FIELD::oNodeAcceleration);}
    auto       &acceleration(       KDN *n ) const {return get<float[3]>(n,KDN_FIELD::oNodeAcceleration);}
    const auto &vbound(       const KDN *n ) const {return get<Bound>(n,KDN_FIELD::oNodeVBnd);}
    auto       &vbound(             KDN *n ) const {return get<Bound>(n,KDN_FIELD::oNodeVBnd);}
    const auto &sphbounds(    const KDN *n ) const {return get<SPHBNDS>(n,KDN_FIELD::oNodeSphBounds);}
    auto       &sphbounds(          KDN *n ) const {return get<SPHBNDS>(n,KDN_FIELD::oNodeSphBounds);}
    const auto &moment(       const KDN *n ) const {return get<FMOMR>(n,KDN_FIELD::oNodeMom);}
    auto       &moment(             KDN *n ) const {return get<FMOMR>(n,KDN_FIELD::oNodeMom);}
    const auto &mass(         const KDN *n ) const {return get<mass_t>(n,KDN_FIELD::oNodeMass);}
    auto       &mass(               KDN *n ) const {return get<mass_t>(n,KDN_FIELD::oNodeMass);}
    const auto &BOB(          const KDN *n ) const {return get<SPHBOB>(n,KDN_FIELD::oNodeBOB);}
    auto       &BOB(                KDN *n ) const {return get<SPHBOB>(n,KDN_FIELD::oNodeBOB);}
    const auto &Ngas(         const KDN *n ) const {return get<int32_t>(n,KDN_FIELD::oNodeNgas);}
    auto       &Ngas(               KDN *n ) const {return get<int32_t>(n,KDN_FIELD::oNodeNgas);}
    const auto &Nstar(        const KDN *n ) const {return get<int32_t>(n,KDN_FIELD::oNodeNstar);}
    auto       &Nstar(              KDN *n ) const {return get<int32_t>(n,KDN_FIELD::oNodeNstar);}
    const auto &Nbh(          const KDN *n ) const {return get<int32_t>(n,KDN_FIELD::oNodeNbh);}
    auto       &Nbh(                KDN *n ) const {return get<int32_t>(n,KDN_FIELD::oNodeNbh);}

    class Node : public SingleElement<KDN,treeStore> {
    public:
        using SingleElement::SingleElement; // Expose the constructors
        auto operator->() const {return p;}
        // When we iterate over a node, we are actually iterating over the particles in the node
        auto begin()              noexcept { return particleStore::iterator(      *store().pstore,store().pstore->Element(p->pLower  )); }
        auto begin()        const noexcept { return particleStore::const_iterator(*store().pstore,store().pstore->Element(p->pLower  )); }
        auto cbegin()       const noexcept { return particleStore::const_iterator(*store().pstore,store().pstore->Element(p->pLower  )); }
        auto end()                noexcept { return particleStore::iterator(      *store().pstore,store().pstore->Element(p->pUpper+1)); }
        auto end()          const noexcept { return particleStore::const_iterator(*store().pstore,store().pstore->Element(p->pUpper+1)); }
        auto cend()         const noexcept { return particleStore::const_iterator(*store().pstore,store().pstore->Element(p->pUpper+1)); }

        auto is_integer_position() {return store().is_integer_position();}
        auto position() const { return store().position(p); }
        auto set_position(blitz::TinyVector<double,3> r) const { return store().set_position(p,r); }
        template<typename T>
        auto &raw_position( KDN *n ) const { return store().raw_position<T>(p); }

        bool have_position()       const noexcept {return have(KDN_FIELD::oNodePosition);}
        bool have_velocity()       const noexcept {return have(KDN_FIELD::oNodeVelocity);}
        bool have_acceleration()   const noexcept {return have(KDN_FIELD::oNodeAcceleration);}
        bool have_BOB()            const noexcept {return have(KDN_FIELD::oNodeBOB);}
        bool have_moment()         const noexcept {return have(KDN_FIELD::oNodeMom);}
        bool have_bound()          const noexcept {return have(KDN_FIELD::oNodeBnd);}
        bool have_sph_bound()      const noexcept {return have(KDN_FIELD::oNodeSphBounds);}
        bool have_velocity_bound() const noexcept {return have(KDN_FIELD::oNodeVBnd);}
        bool have_Ngas()           const noexcept {return have(KDN_FIELD::oNodeNgas);}
        bool have_Nstar()          const noexcept {return have(KDN_FIELD::oNodeNstar);}
        bool have_Nbh()            const noexcept {return have(KDN_FIELD::oNodeNbh);}

        auto bound() const { return store().bound(p); }
        auto set_bound(const Bound &bnd) const { return store().set_bound(p,bnd); }
        template<typename T>
        auto &raw_bound()    const {return store().raw_bound<T>(p); }

        auto &velocity()     const {return store().velocity(p);}
        auto &acceleration() const {return store().acceleration(p);}
        auto &vbound()       const {return store().vbound(p);}
        auto &sphbounds()    const {return store().sphbounds(p);}
        auto &moment()       const {return store().moment(p);}
        auto &mass()         const {return store().mass(p);}
        auto &BOB()          const {return store().BOB(p);}
        auto &Ngas()         const {return store().Ngas(p);}
        auto &Nstar()        const {return store().Nstar(p);}
        auto &Nbh()          const {return store().Nbh(p);}

        auto bMax()         const noexcept { return p->bMax; }
        auto &bMax()              noexcept { return p->bMax; }
        auto fSoft2()       const noexcept { return p->fSoft2; }
        auto &fSoft2()            noexcept { return p->fSoft2; }
        auto lower()        const noexcept { return p->pLower; }
        auto &lower()             noexcept { return p->pLower; }
        auto upper()        const noexcept { return p->pUpper; }
        auto &upper()             noexcept { return p->pUpper; }
        auto count()        const noexcept { return upper() + 1 - lower(); }
        auto lchild()       const noexcept { return p->iLower; }
        auto rchild()       const noexcept { return p->iLower+1; }
        auto lnode()        const noexcept { return NodePointer(store(),lchild()); }
        auto rnode()        const noexcept { return NodePointer(store(),rchild()); }
        auto remote_id()    const noexcept { return p->pLower; }

        bool is_cell()      const noexcept { return p->iLower!=0; }
        bool is_bucket()    const noexcept { return p->iLower==0; }
        bool is_top_tree()  const noexcept { return p->bTopTree; }
        bool is_group()     const noexcept { return p->bGroup; }
        bool is_remote()    const noexcept { return p->bRemote; }
        bool is_local()     const noexcept { return !is_remote(); }
        bool is_marked()    const noexcept { return p->bHasMarked; }
        bool is_NN()        const noexcept { return p->bHasNNflag; }

        bool set_top_tree(bool b) noexcept { return (p->bTopTree=b); }
        bool set_group(bool b)    noexcept { return (p->bGroup=b); }
        bool set_marked(bool b)   noexcept { return (p->bHasMarked=b); }
        bool set_NN(bool b)       noexcept { return (p->bHasNNflag=b); }

        uint8_t min_rung() const noexcept { return p->uMinRung; }
        uint8_t max_rung() const noexcept { return p->uMaxRung; }
        void    set_min_rung(uint8_t uRung )  { p->uMinRung = uRung; }
        void    set_max_rung(uint8_t uRung )  { p->uMaxRung = uRung; }
        void    set_rung(uint8_t uRung) { p->uMinRung = p->uMaxRung = uRung;}
        void    set_rung(uint8_t uMinRung,uint8_t uMaxRung) {
            p->uMinRung = uMinRung;
            p->uMaxRung = uMaxRung;
        }

        int  depth() const noexcept { return p->iDepth; }
        void set_depth(int iDepth) noexcept { p->iDepth = iDepth; }

        int  split_dim() const noexcept { return p->iSplitDim; }
        void set_split_dim(int iSplitDim) noexcept { p->iSplitDim = iSplitDim; }

        void set_leaf() {p->iLower=0;}
        void set_local(int l=0, int u=0, int child=0) {
            p->pLower = l;
            p->pUpper = u;
            p->bRemote = 0;
            p->bTopTree = 0;
            p->iLower = child;
        }
        void set_remote(int id) {
            p->bTopTree = 1;
            p->bGroup = 0;
            p->bRemote = 1;
            p->pUpper = p->pLower = id;
        }
        void set_top_tree(int u) {
            p->bTopTree = 1;            /* The replicated top tree */
            p->bGroup = 0;              /* top tree can never be a group */
            p->bRemote = 0;             /* top tree is not remote */
            p->iLower = 1;              /* Relative index to lower cell */
            p->pUpper = u;              /* Relative index to upper cell */
            p->pLower = 0;              /* Not used */
        }

        /// @brief Split this cell into two children at the i'th particle
        /// @param i Absolute particle index of the split
        /// @return The left and right childs as a tuple of NodePointer
        /// The split must be between lower() and upper() inclusive.
        /// Particles below "i" will be in the left child and particles
        /// above and including "i" will be in the right child. The children
        /// are initialized as buckets (but can be later split).
        auto split(int i) {
            assert(i>lower() && i<=upper());
            p->iLower = store().AllocNode(2);
            NodePointer pLeft(store(),lchild());
            pLeft->set_local(lower(),i-1);
            pLeft->set_depth(depth()+1);
            pLeft->set_split_dim(3);
            NodePointer pRight(store(),rchild());
            pRight->set_local(i,upper());
            pRight->set_depth(depth()+1);
            pRight->set_split_dim(3);
            return std::make_tuple(pLeft,pRight);
        }

        auto get_child_cells(int id) {
            int idxLower = lchild();            // This is always true
            int idxUpper = rchild();            // This is usually true
            int idLower = id, idUpper = id;     // Default
            if (is_top_tree()) {
                if (is_remote()) {          // We may point off node, but then nodes are adjacent
                    idLower = lower();
                    idUpper = idLower;
                }
                else { idxUpper = upper(); }
            }
            return std::make_tuple(idxLower,idLower,idxUpper,idUpper);
        }

        void calc_open(double minside) {
            auto bnd = bound();
            auto d2 = bnd.maxdist(position());
            auto b = std::max(bnd.maxside(),minside);
            if (b*b < d2) b = sqrt(d2);
            bMax() = b;
        }

    protected:
        /// @brief update the node bound from the particles
        /// @tparam T Raw position type (double or int32_t)
        /// It is more efficient to scan the particles and min/max the coordinates if we do so
        /// using the native type (int32_t or double). This template allows us to do that.
        template<typename T>
        void update_bound() {
            raw_bound<typename bound_type<T>::type>() = store().pstore->raw_bound<T>(begin(),end());
        }
    public:
        /// @brief Update the bounds of the cell to shrink wrap the particles
        /// This will calculate the min and max coordinates of all of the particles and update
        /// the bounds of the cell to be the minimum box that encloses all of them.
        void update_bound() {
            if (is_integer_position()) update_bound<int32_t>();
            else update_bound<double>();
        }
        /// @brief returns true if there could be active particles in this subtree
        /// @param uRungLo Minimum active rung
        /// @param uRungHi Maximum active rung
        /// @return false if there are no active particles and true if there could be
        /// This is really means that the cell MIGHT have an active particle, but it
        /// is not for certain, since the contained DstActive particles may not lie in
        /// the rung range. To be certain that there are actually active particles
        /// one has to look at all particles of this cell (recursively walking the
        /// subcells).
        bool is_active(uint8_t uRungLo,uint8_t uRungHi) {
            return uRungLo <= p->uMaxRung && uRungHi >= p->uMinRung;
        }
    };

    /// @brief A "pointer" reference to a "Node"
    ///
    /// This uses the default copy/assignment semantics (the internal pointers are copied).
    using NodePointer = SinglePointer<Node>;

    /// @brief A "reference" reference to a "Node"
    ///
    /// Assignent of a reference copies the underlying data
    using NodeReference = SingleReference<Node>;

    auto operator[](int i) {
        return NodePointer(*this,i);
    }
    auto operator[](KDN *p) {
        return NodePointer(*this,p);
    }
    const auto operator[](const KDN *p) const {
        return NodePointer(*this,p);
    }
};
#endif