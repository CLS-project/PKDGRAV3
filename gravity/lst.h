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
#pragma once
#include "pkd_config.h"

#include <boost/intrusive/list.hpp>
#include <cstdint>
#include <cstring>
#include <type_traits>
#include <array>
#include <tuple>

#include "blitz/array.h"
#include "core/simd.h"

/*
** Generic List Management
** A list consists of a number of fields orgainized in way that they can be efficiently accessed and processed on a CPU or GPU.
** Each field has a type (e.g., float, int32_t, etc.) and a name (e.g., dx, dy, dz, etc.). The fields are stored in blocks of "n" elements (e.g., 32).
** Several blocks are stored contiguously in "tiles" of "m" blocks (e.g., 8), and the tiles are linked together in a list.
**
*/

/*
   FIELD-LIST MEMORY LAYOUT (CPU/GPU FRIENDLY)

   - The LIST is a boost::intrusive::list of TILE objects.
   - Each TILE is a std::array<Block, M>.
   - Each BLOCK stores fields as ilBlock<T, N> (SoA layout).


   LIST (boost::intrusive::list<Tile>)
   ─────────────────────────────────────────────────────────

   [TILE] <--> [TILE] <--> [TILE] <--> ...   (doubly-linked, intrusive)


   TILE STRUCTURE
   ─────────────────────────────────────────────────────────

   +------------------------------------------------------+
   | TILE                                                 |
   |------------------------------------------------------|
   | intrusive hooks (prev / next pointers)               |
   |                                                      |
   | std::array<Block, M> blocks                          |
   +------------------------------------------------------+


   BLOCK STRUCTURE
   ─────────────────────────────────────────────────────────

   +------------------------------------------------------+
   | BLOCK                                                |
   |------------------------------------------------------|
   | ilBlock<float, N> dx  : dx0 dx1 ... dx(N-1)          |
   | ilBlock<float, N> dy  : dy0 dy1 ... dy(N-1)          |
   | ilBlock<float, N> dz  : dz0 dz1 ... dz(N-1)          |
   | ... (other fields: int32_t, float, etc.)             |
   +------------------------------------------------------+


   ilBlock<T, N>
   ─────────────────────────────────────────────────────────

   +------------------------------------------------------+
   | ilBlock<T, N>                                        |
   |------------------------------------------------------|
   | T [N] : v0 v1 v2 ... v(N-1)                          |
   +------------------------------------------------------+


   Notes:
   - Intrusive list: no external node allocation (TILE contains links).
   - Excellent locality inside TILE (std::array of BLOCKS).
   - SoA inside BLOCK enables SIMD + GPU coalescing.
   - Scales dynamically by linking/unlinking TILEs.
*/

// We use __VA_ARGS__ in case the type of the field contains a comma (e.g., when using a template definition).
#define ILIST_PROTO_BLOCK_(name, ...) ilist::ilBlock<__VA_ARGS__, width> name; // Definition inside the block structure
#define ILIST_PROTO_ITEMS_(name, ...) __VA_ARGS__ name;                        // Definition inside the items structure
#define ILIST_ASSIGN_ITEM_(name, ...) b->name[i] = value.name;                 // Assignment from items to block when appending
#define ILIST_DEFINE(NAME, FIELD_LIST)                                                                                                                         \
    struct Items##NAME {                                                                                                                                       \
        FIELD_LIST(ILIST_PROTO_ITEMS_)                                                                                                                         \
    };                                                                                                                                                         \
    template <std::size_t N>                                                                                                                                   \
    struct Block##NAME {                                                                                                                                       \
        static constexpr std::size_t width = N;                                                                                                                \
        FIELD_LIST(ILIST_PROTO_BLOCK_)                                                                                                                         \
    };                                                                                                                                                         \
    template <int N, int M>                                                                                                                                    \
    class List##NAME : public ilist::List<Block##NAME<N>, M> {                                                                                                 \
      public:                                                                                                                                                  \
        using block_type = Block##NAME<N>;                                                                                                                     \
        using items_type = Items##NAME;                                                                                                                        \
        using base_type = ilist::List<block_type, M>;                                                                                                          \
        List##NAME() = default;                                                                                                                                \
        List##NAME(typename ilist::List<Block##NAME<N>, M>::list_type& freeList) : ilist::List<Block##NAME<N>, M>(freeList) {}                                 \
        void append(const items_type& value) {                                                                                                                 \
            auto [b, i] = this->create();                                                                                                                      \
            FIELD_LIST(ILIST_ASSIGN_ITEM_)                                                                                                                     \
        }                                                                                                                                                      \
    };

#if defined(__CUDACC__)
    #define LST_GPU_HD __host__ __device__
#else
    #define LST_GPU_HD
#endif

namespace ilist {
// This class maintains a center reference for the interaction list, which allows us to store relative positions as float in the list while avoiding precision
// issues. The center reference is typically set to the position of the first particle we encounter when we are building the interaction list.
class ilCenterReference {
    using double3 = blitz::TinyVector<double, 3>;
    using float3 = blitz::TinyVector<float, 3>;
    double3 dCenter = 0.0;

  public:
    [[nodiscard]] auto getReference() const { return dCenter; }
    [[nodiscard]] auto getReference(int i) const { return dCenter[i]; }
    void setReference(double3 r) { dCenter = r; }
    [[nodiscard]] auto get_dr(const double3& r) const { return float3(r - dCenter); }
};

// This defines a block of "n" elements (e.g., 32) for each field in the interaction list
template <typename SCALAR, typename VECTOR, int N>
struct alignas(fvec::width() * sizeof(float)) ilBlockBase {
    using scalar_t = SCALAR;
    using vector_t = VECTOR;
    static constexpr int scalar_count = N;
    static constexpr int vector_count = N * sizeof(scalar_t) / (sizeof(vector_t));
    union {
        scalar_t s[scalar_count]; // NOLINT(modernize-avoid-c-arrays)
        vector_t v[vector_count]; // NOLINT(modernize-avoid-c-arrays)
    };
    LST_GPU_HD auto operator[](std::size_t i) -> auto& { return s[i]; }
    LST_GPU_HD auto operator[](std::size_t i) const -> const auto& { return s[i]; }
};

// CAREFUL: the iBlockBase structure above must be properly aligned for this to work.
// The alignof(fvec) is very important or the GPU version will not work properly.
template <typename V, int N> struct ilBlock;
#if defined(__CUDACC__) || defined(__METAL_VERSION__)
template <int N> struct ilBlock<float, N> : public ilBlockBase<float, float, N> {};
template <int N> struct ilBlock<double, N> : public ilBlockBase<double, double, N> {};
template <int N> struct ilBlock<int32_t, N> : public ilBlockBase<int32_t, int32_t, N> {};
template <int N> struct ilBlock<int64_t, N> : public ilBlockBase<int64_t, int64_t, N> {};
#else
template <int N> struct ilBlock<float, N> : public ilBlockBase<float, fvec, N> {};
template <int N> struct ilBlock<double, N> : public ilBlockBase<double, dvec, N> {};
template <int N> struct ilBlock<int32_t, N> : public ilBlockBase<int32_t, i32v, N> {};
template <int N> struct ilBlock<int64_t, N> : public ilBlockBase<int64_t, i64v, N> {};
#endif

template <typename B, int N>
class Tile : public boost::intrusive::list_base_hook<>, public std::array<B, N> {
    template <typename B_, int N_> friend class List;

  protected:
    int nInTile = 0;

  public:
    int nRefs = 0;
    using block_type = B;
    static constexpr int num_blocks = N;
    static constexpr int width = B::width;
    static constexpr int max_interactions = num_blocks * width;
    Tile() { memset(this->data(), 0, N * sizeof(B)); }
    void clear() { nInTile = 0; }
    [[nodiscard]] auto count() const -> int { return nInTile; }
    [[nodiscard]] auto full() { return count() == max_interactions; }
    [[nodiscard]] auto empty() { return count() == 0; }
    [[nodiscard]] auto create() {
        auto iBlock = nInTile / width;
        auto iInter = nInTile % width;
        ++nInTile;
        return std::make_tuple(&(*this)[iBlock], iInter);
    }
};

// This will results tiles with "N" blocks of type "B".
template <typename B, int N = 8>
class List : public boost::intrusive::list<Tile<B, N>, boost::intrusive::constant_time_size<true>> {
  public:
    using tile_type = Tile<B, N>;
    using block_type = typename tile_type::block_type;
    using list_type = boost::intrusive::list<Tile<B, N>, boost::intrusive::constant_time_size<true>>;

  protected:
    list_type freeListDefault; // In case we have a private "free list"
    list_type* freeList;       // This can be the private free list, or a global shared free list
    struct delete_disposer {
        void operator()(tile_type* delete_this) { delete delete_this; }
    };
    struct move_disposer {
        list_type& other;
        explicit move_disposer(list_type& other) : other(other) {}
        void operator()(tile_type* tile) { other.push_back(*tile); }
    };
    // struct isref0 {
    //     explicit isref0() = default;
    //     auto operator()(const Tile<B, N>& tile) const -> bool { return tile.nRefs == 0; }
    // };
    void newTile() {
        if (freeList->size()) { // Pull a tile from the free list if we have any
            auto& t = freeList->front();
            freeList->pop_front();
            t.clear();
            this->push_back(t);
        }
        else this->push_back(*new tile_type);
        this->back().nRefs = 1;
    }

  public:
    static constexpr int max_blocks = N;
    static constexpr int max_interactions = N * B::width;
    [[nodiscard]] auto memory() -> std::size_t { return (freeListDefault.size() + this->size()) * sizeof(tile_type); }
    [[nodiscard]] auto count() -> std::size_t { return this->size() ? this->back().count() + (this->size() - 1) * max_interactions : 0; }
    void clear() { this->clear_and_dispose(move_disposer(*freeList)); }
    void shrink_to_fit() { freeListDefault.clear_and_dispose(delete_disposer()); }
    List() : freeList(&freeListDefault) {}
    List(list_type& freeList) : freeList(&freeList) {}
    ~List() {
        clear();
        shrink_to_fit(); // Delete our own free list. There will be no elements if using a shared list.
    }

    // Copy constructor
    List(const List& other) {
        if (other.freeList == &other.freeListDefault) freeList = &freeListDefault;
        else freeList = other.freeList; // keep shared/external association
        clone(other);
    }

    // Copy assignment
    auto operator=(const List& other) -> List& {
        if (this != &other) {
            clear();
            clone(other);
        }
        return *this;
    }

    // Move constructor. We need to:
    // - move elements in the active list
    // - move elements in our local free list (if any)
    // - update the free list pointer
    List(List&& other) noexcept(std::is_nothrow_move_constructible_v<list_type>)
        : list_type(std::move(static_cast<list_type&>(other))), freeListDefault(std::move(other.freeListDefault)), freeList(nullptr) {
        if (other.freeList == &other.freeListDefault) freeList = &freeListDefault;
        else freeList = other.freeList; // keep shared/external association
        // We leave the other freeList unchanged so it will continue to use the shared or its own free list.
        // Typically the destructor will be called and the local free list is cleaned up, and the shared ignored.
    }

    // Move assignment. Safest here is to perform a swap of the active list, and the default free list.
    auto operator=(List&& other) noexcept(std::is_nothrow_move_assignable_v<list_type>) -> List& {
        if (this != &other) this->swap(other);
        return *this;
    }

    void swap(List& other) noexcept(noexcept(std::declval<list_type&>().swap(std::declval<list_type&>()))) {
        using std::swap;

        // swap the base part
        list_type& self_base = *this;
        list_type& other_base = other;
        self_base.swap(other_base);

        // swap the member free list
        swap(freeListDefault, other.freeListDefault);

        // Now set the free list, preserving shared/external association
        list_type* new_freeList_self = (other.freeList == &other.freeListDefault) ? &freeListDefault : other.freeList;
        list_type* new_freeList_other = (freeList == &freeListDefault) ? &other.freeListDefault : freeList;
        freeList = new_freeList_self;
        other.freeList = new_freeList_other;
    }

    // ADL-enabled non-member swap
    friend void swap(List& a, List& b) noexcept(noexcept(a.swap(b))) { a.swap(b); }

    auto clone(const list_type& rhs) -> auto& {
        clear();
        for (auto const& tile : rhs) {
            newTile();
            memcpy(this->back().data(), tile.data(), sizeof(block_type) * tile_type::num_blocks);
            this->back().nInTile = tile.nInTile;
        }
        return *this;
    }
    auto create() {
        if (this->size() == 0 || this->back().full()) newTile();
        return this->back().create();
    }
};

template <typename RESULT, typename BLOCK> struct EvalBlock;
template <typename TILE, typename EVAL>
auto EvalTile(TILE& tile, EVAL& eval) -> typename EVAL::result_type {
    typename EVAL::result_type result;
    result.zero();

    // Do the full blocks, followed by the last (probably not full) block
    auto nBlocks = tile.count() / tile.width;
    for (auto iBlock = 0; iBlock < nBlocks; ++iBlock) {
        result += eval(tile.width, tile[iBlock]);
    }
    result += eval(tile.count() - nBlocks * tile.width, tile[nBlocks]);
    return result;
}
} // namespace ilist

#undef LST_GPU_HD