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

/*
** Generic List Management
**
*/

#ifndef LST_H
#define LST_H
#include "pkd_config.h"
#include <stdlib.h>
#include <stdint.h>

#include <cstdint>
#include <cstring>
//#include <memory>
//#include <cmath>
#include <type_traits>
#include <array>
#include <tuple>

#include <boost/intrusive/list.hpp>
#include <boost/fusion/include/algorithm.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/adapt_assoc_struct.hpp>
#include <boost/preprocessor/seq.hpp>
#include <boost/fusion/include/map.hpp>
#include <boost/fusion/include/at_key.hpp>

#include "core/simd.h"

// ILIST_FIELD_VALUES: comma separated list of field values, e.g., "dx, dy, dz, ..., u" or "float, float, ..."
#define ILIST_FIELD_VALUE_(r, data, elem) BOOST_PP_TUPLE_ELEM(2,data,elem),
#define ILIST_FIELD_VALUES(seq,j) BOOST_PP_SEQ_FOR_EACH(ILIST_FIELD_VALUE_, j, BOOST_PP_SEQ_POP_BACK(seq))BOOST_PP_TUPLE_ELEM(j,BOOST_PP_SEQ_HEAD(BOOST_PP_SEQ_REVERSE(seq)))

// ILIST_FIELD_PROTOS: function prototype parameters, e.g., "float dx, float dy, ..., float u"
#define ILIST_PROTO_(t) BOOST_PP_TUPLE_ELEM(2,0,t) BOOST_PP_TUPLE_ELEM(2,1,t)
#define ILIST_FIELD_PROTO_(r, data, elem) ILIST_PROTO_(elem),
#define ILIST_FIELD_PROTOS(seq) BOOST_PP_SEQ_FOR_EACH(ILIST_FIELD_PROTO_, _, BOOST_PP_SEQ_POP_BACK(seq))ILIST_PROTO_(BOOST_PP_SEQ_HEAD(BOOST_PP_SEQ_REVERSE(seq)))

#define ILIST_BLK_KEY_(r, data, elem) struct BOOST_PP_TUPLE_ELEM(2,1,elem);
#define ILIST_BLK_KEYS(seq) BOOST_PP_SEQ_FOR_EACH(ILIST_BLK_KEY_,_,seq)

// Structure with all fields as blocks, e.g.,
//   ilist::ilBlock<float,SIZE::value> dx;
//   ilist::ilBlock<float,SIZE::value> dy;
//   ...
#define ILIST_BLK_FIELD_(r, data, elem) ilist::ilBlock<BOOST_PP_TUPLE_ELEM(2,0,elem),SIZE::value> BOOST_PP_TUPLE_ELEM(2,1,elem);
#define ILIST_BLK_FIELDS(seq) BOOST_PP_SEQ_FOR_EACH(ILIST_BLK_FIELD_,_,seq)

#define ILIST_FIELD_KEY_(r, data, elem) (auto,BOOST_PP_TUPLE_ELEM(2,1,elem),Keys##data::BOOST_PP_TUPLE_ELEM(2,1,elem))
#define ILIST_FIELD_KEYS(seq,NAME) BOOST_PP_SEQ_FOR_EACH(ILIST_FIELD_KEY_,NAME,seq)

namespace ilist {
class ilCenterReference {
    double dCenter[3] {0.0,0.0,0.0};
public:
    double getReference(int i) const {return dCenter[i];}
    void setReference(double x,double y,double z) {
        dCenter[0] = x;
        dCenter[1] = y;
        dCenter[2] = z;
    }
    auto get_dr(double x,double y,double z) {
        return std::make_tuple<float,float,float>(x-dCenter[0],y-dCenter[1],z-dCenter[2]);
    }
};

// template<typename SCALAR,typename VECTOR,int N> struct ilBlockBase;
// template<typename SCALAR,typename VECTOR,int N>
// class ilBlockElement {
//     struct ilBlockBase<SCALAR,VECTOR,N> &base;
//     std::size_t i;
// public:
//     typedef SCALAR scalar_t;
//     typedef VECTOR vector_t;
//     operator SCALAR() const {return base.s[i];}
//     operator VECTOR() const {return base.v[i];}
//     struct ilBlockBase<SCALAR,VECTOR,N> &operator=(SCALAR s) {base.s[i]=s; return base;}
//     struct ilBlockBase<SCALAR,VECTOR,N> &operator=(VECTOR v) {base.v[i]=v; return base;}
//     ilBlockElement() = default;
//     explicit ilBlockElement(ilBlockBase<SCALAR,VECTOR,N> &base,std::size_t i) : base(base), i(i) {}
// };

template<typename SCALAR,typename VECTOR,int N>
struct alignas(SIMD_WIDTH *sizeof(float)) ilBlockBase {
    static constexpr int scalar_count = N;
#if !defined(__CUDACC__)
    typedef SCALAR scalar_t;
    typedef VECTOR vector_t;
    static constexpr int vector_count = N * sizeof(scalar_t) / (sizeof(vector_t));
    union {
        scalar_t s[scalar_count];
        vector_t v[vector_count];
    };
    auto &operator[](std::size_t i) {return s[i];}
    const auto &operator[](std::size_t i) const {return s[i];}
    // auto operator[](std::size_t i) {return ilBlockElement(*this,i);}
    // const auto operator[](std::size_t i) const {return ilBlockElement(*this,i);}
#else
    typedef SCALAR scalar_t;
    scalar_t s[N];
    __host__ __device__ auto &operator[](std::size_t i) {return s[i];}
    __host__ __device__ const auto &operator[](std::size_t i) const {return s[i];}
#endif
};

// CAREFUL: the iBlockBase structure above must be properly aligned for this to work.
// The alignof(fvec) is very important or the GPU version will not work properly.
template<typename V,int N> struct ilBlock;
#if !defined(__CUDACC__)
template<int N> struct ilBlock<float,N>  : public ilBlockBase<float,  fvec,N> {};
template<int N> struct ilBlock<double,N> : public ilBlockBase<double, dvec,N> {};
template<int N> struct ilBlock<int32_t,N>: public ilBlockBase<int32_t,i32v,N> {};
template<int N> struct ilBlock<int64_t,N>: public ilBlockBase<int64_t,i64v,N> {};
#else
template<int N> struct ilBlock<float,N>  : public ilBlockBase<float,  float,N> {};
template<int N> struct ilBlock<double,N> : public ilBlockBase<double, double,N> {};
template<int N> struct ilBlock<int32_t,N>: public ilBlockBase<int32_t,int32_t,N> {};
template<int N> struct ilBlock<int64_t,N>: public ilBlockBase<int64_t,int64_t,N> {};
#endif
namespace {
// Here we memcpy a single field from source to destination, e.g.,
//   dst.dx[ 0:31] = src.dx[0:31] (iteration 1)
//   dst.dx[31:63] = src.dx[0:31] (iteration 2, with src pointing to next block)
//   ... etc.
template<typename NGPU,typename NCPU>
struct copyone {
    int i;
    explicit copyone(int i) : i(i) {}
    template <typename ZipView>
    void operator() (const ZipView &t) const {
        std::memcpy(boost::fusion::at_c<0>(t).s+i,boost::fusion::at_c<1>(t).s,sizeof(boost::fusion::at_c<1>(t)));
    }
};
}

// This function iterates over the fields of the block (dx, dy, dz, ...) and calls copyone to memcpy the
// source (CPU) field into the appropriate destination (GPU) field. The destination size has to be at least
// as large as the source, and an even multiple of the size. If the destination block size is larger than
// the source then it iterates to fill up the destination.
template<template<typename> class BLOCK,typename NGPU,typename NCPU>
void copy(BLOCK<NGPU> *dst, BLOCK<NCPU> *src,std::enable_if_t<(NGPU::value>NCPU::value) &&NGPU::value%NCPU::value==0,bool> = true) {
    typedef boost::fusion::vector<BLOCK<NGPU> &, BLOCK<NCPU> &> sequences;
    for (auto i=0; i<NGPU::value; i+=NCPU::value) {
        boost::fusion::for_each(boost::fusion::zip_view<sequences>(sequences(*dst,*src)), copyone<NGPU,NCPU>(i));
        ++src;
    }
}

template<template<typename> class BLOCK,typename NGPU,typename NCPU>
void copy(BLOCK<NGPU> *dst, BLOCK<NCPU> *src,std::enable_if_t<NGPU::value==NCPU::value,bool> = true) {
    std::memcpy(dst,src,sizeof(*src));
}

template<typename B,int N>
class Tile : public boost::intrusive::list_base_hook<>, public std::array<B,N> {
    template<typename B_,int N_> friend class List;
protected:
    int nInTile = 0;
public:
    int nRefs = 0;
    typedef B block_type;
    static constexpr int num_blocks = N;
    static constexpr int width=  B::width;
    static constexpr int max_interactions = num_blocks * width;
    Tile() {memset(this->data(),0,N*sizeof(B));}
    void clear() {nInTile=0;}
    int count() const { return nInTile; }
    auto full() { return count()==max_interactions; }
    auto empty() { return count()==0; }
    auto create() {
        auto iBlock = nInTile / width;
        auto iInter = nInTile % width;
        ++nInTile;
        return std::make_tuple(&(*this)[iBlock],iInter);
    }
};

template<typename B,int N=32>
class List : public boost::intrusive::list<Tile<B,N>,boost::intrusive::constant_time_size<true>> {
public:
    typedef Tile<B,N> tile_type;
    typedef typename tile_type::block_type block_type;
protected:
    typedef boost::intrusive::list<Tile<B,N>,boost::intrusive::constant_time_size<true>> list_type;
    list_type freeListDefault;
    list_type *freeList;
    struct delete_disposer {
        void operator()(tile_type *delete_this) {  delete delete_this;  }
    };
    struct move_disposer {
        list_type &other;
        explicit move_disposer(list_type &other) : other(other) {}
        void operator() (tile_type *tile) {
            other.push_back(*tile);
        }
    };
    struct isref0 {
        explicit isref0() {}
        bool operator() (const Tile<B,N> &tile) const {
            return tile.nRefs==0;
        }
    };
    struct releaseone {
        list_type &free, &busy;
        explicit releaseone(list_type &free,list_type &busy) : free(free), busy(busy) {}
        void operator() (Tile<B,N> *tile) {
            if (--tile->nRefs) busy.push_back(*tile);
            else free.push_back(*tile);
        }
    };
    void newTile() {
        if (freeList->size()) {                  // Pull a tile from the free list if we have any
            auto &t = freeList->front();
            freeList->pop_front();
            t.clear();
            this->push_back(t);
        }
        else this->push_back(* new tile_type);
        this->back().nRefs = 1;
    }
    void setFreeList(list_type &newFreeList) { freeList = &newFreeList; }
public:
    static constexpr int max_blocks = N;
    static constexpr int max_interactions = N * B::width;
    std::size_t memory() { return (freeList->size() + this->size()) * sizeof(tile_type); }
    std::size_t count() { return this->size() ? this->back().count() + (this->size()-1) * max_interactions : 0; }
    void clear() {
        this->clear_and_dispose(move_disposer(*freeList));
        //this->clear_and_dispose(releaseone(freeList,busyList));
    }
    List() : freeList(&freeListDefault) {}
    ~List() {
        this->clear_and_dispose(delete_disposer());
        freeList->clear_and_dispose(delete_disposer());
    }
    auto &clone(const list_type &rhs) {
        this->clear_and_dispose(move_disposer(*freeList));
        for (auto tile : rhs) {
            newTile();
            memcpy(this->back().data(),tile.data(),sizeof(block_type) * tile_type::num_blocks);
            this->back().nInTile = tile.nInTile;
        }
        return *this;
    }
    auto create() {
        if (this->size()==0 || this->back().full()) newTile();
        return this->back().create();
    }

};
struct saveone {
    int i;
    explicit saveone(int i) : i(i) {}
    template <typename ZipView>
    void operator() (const ZipView &t) const {
        boost::fusion::at_c<0>(t).s[i] = boost::fusion::at_c<1>(t);
    }
};

template<typename RESULT,typename BLOCK> struct EvalBlock;
template<typename TILE,typename EVAL>
typename EVAL::result_type EvalTile(TILE &tile,EVAL &eval) {
    typename EVAL::result_type result;
    result.zero();

    // Do the full blocks, followed by the last (probably not full) block
    auto nBlocks = tile.count() / tile.width;
    for (auto iBlock=0; iBlock<nBlocks; ++iBlock) {
        result += eval(tile.width,tile[iBlock]);
    }
    result += eval(tile.count() - nBlocks*tile.width,tile[nBlocks]);
    return result;
}
} // namespace ilist

// Constructs: boost::fusion::at_key<KeysPC::dx>(*b)[i] = dx;
#define ILIST_ASSIGN_FIELD_(r, data, elem)\
 boost::fusion::at_key<BOOST_PP_TUPLE_ELEM(4,0,data)::BOOST_PP_TUPLE_ELEM(2,1,elem)>\
 (*BOOST_PP_TUPLE_ELEM(4,1,data))[BOOST_PP_TUPLE_ELEM(4,2,data)]\
  = BOOST_PP_TUPLE_ELEM(4,3,data)BOOST_PP_TUPLE_ELEM(2,1,elem);
#define ILIST_ASSIGN_FIELDS(seq,NAME,B,I,P) BOOST_PP_SEQ_FOR_EACH(ILIST_ASSIGN_FIELD_,(Keys##NAME,B,I,P),seq)

#define ILIST_DECLARE(NAME,FIELDS_SEQ)                                                                      \
template <typename SIZE> struct BLOCK##NAME {                                                               \
    static constexpr int width = SIZE::value;                                                               \
    ILIST_BLK_FIELDS(FIELDS_SEQ)                                                                            \
    };                                                                                                      \
namespace Keys##NAME {ILIST_BLK_KEYS(FIELDS_SEQ)}                                                           \
BOOST_FUSION_ADAPT_ASSOC_TPL_STRUCT(                                                                        \
    (SIZE),                                                                                                 \
    (BLOCK##NAME) (SIZE),                                                                                   \
    /*ILIST_FIELD_VALUES(FIELDS_SEQ,1)*/                                                                    \
    ILIST_FIELD_KEYS(FIELDS_SEQ,NAME)                                                                       \
)                                                                                                           \
template<int N> using Block##NAME = BLOCK##NAME<boost::mpl::int_<N>>;                                       \
template<int N,int M> using Tile##NAME = ilist::Tile<Block##NAME<N>,M>;                                     \
template<int N,int M>                                                                                       \
class List##NAME : public ilist::List<Block##NAME<N>,M> {                                                   \
public:                                                                                                     \
    void append(ILIST_FIELD_PROTOS(FIELDS_SEQ)) {                                                           \
        Block##NAME<N> *b; int i;                                                                           \
        std::tie(b,i) = this->create();                                                                     \
        ILIST_ASSIGN_FIELDS(FIELDS_SEQ,NAME,b,i,)                                                           \
    }                                                                                                       \
};                                                                                                          \
/**/
// typedef boost::fusion::vector<ILIST_FIELD_VALUES(FIELDS_SEQ,0)> argvec;
// typedef boost::fusion::vector<BLOCK##NAME<boost::mpl::int_<N>> &, argvec &> sequences;
// argvec args(ILIST_FIELD_VALUES(FIELDS_SEQ,1));
// boost::fusion::for_each(boost::fusion::zip_view<sequences>(sequences(*b,args)),ilist::saveone(i));
#endif
