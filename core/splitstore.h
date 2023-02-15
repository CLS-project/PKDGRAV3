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

#ifndef CORE_SPLITSTORE_H
#define CORE_SPLITSTORE_H

#include "datafields.h"

//! \brief A generic contiguous array of fixed-length elements.
//!
//! The length of the elements is choosen at run-time before the
//! container is created. This class provides a way to access the
//! individual fields
template<typename DATA,typename FIELD>
class splitStore : public dataFields<DATA,FIELD> {
public:
    using value_type = DATA;
    using dataFields<DATA,FIELD>::Element;
    using dataFields<DATA,FIELD>::ElementSize;
    splitStore() = default;
protected:
    std::vector<DATA *> tiles; //!< Array of "tiles"

    int nBitsLo = 14;
    int nBitsHi = 18;
    int iMask = (1<<nBitsLo) - 1;
    int nTilesReserved = 0;

    int nMaxNodes = 0;              //!< Maximum number of element that can be stored
    int nNodes = 0;                 //!< Current number of elements

    void extend() {
        auto tile_size = ElementSize() * (1<<nBitsLo);
        tiles.push_back(static_cast<DATA *>(malloc(tile_size)));
        nMaxNodes = tiles.size() * (1<<nBitsLo);
    }

public:
    auto Nodes() const { return nNodes; }

    void SetNodeCount(int n) {
        nNodes = n;
        while (nNodes > nMaxNodes) extend();
    }

    auto memory() const {
        return tiles.size() * ElementSize() * (1<<nBitsLo);
    }

    auto AlignNode() {
        if (nNodes&1) ++nNodes;
        return nNodes;
    }

    //! Allocate one or more elements. Extend if necessary
    //! \param n Number of elements to allocate
    auto AllocNode(int n=1) {
        int iNode = nNodes;
        nNodes += n;
        while (nNodes > nMaxNodes) extend();
        return iNode;
    }

    //! Set the size and location of the storage.
    //! \param p Pointer to the block of storage
    //! \param n Maximum number of element
    void setStore(int Lo, int Hi, uint64_t nData=0, void *pData=nullptr) {
        nBitsLo = Lo;
        nBitsHi = Hi;
        iMask = (1<<nBitsLo) - 1;
        tiles.reserve(1<<nBitsHi);
        tiles.clear();
        // If we were given a block of storage then try to use it for the first block of tiles
        auto tile_size = ElementSize() * (1<<nBitsLo);
        auto p = static_cast<char *>(pData);
        while (nData>tile_size) {
            tiles.push_back(reinterpret_cast<DATA *>(p));
            nData -= tile_size;
            p += tile_size;
        }
        nTilesReserved = tiles.size();
        nMaxNodes = tiles.size() * (1<<nBitsLo);
        if (tiles.empty()) extend(); // Make sure that we have at least one tile
    }

    DATA *Element(int iNode) const { return this->Element(tiles[iNode>>nBitsLo],iNode&iMask); }
    // int FreeStore() const { return nStore; }
    // int Local() const { return nLocal; }
    // int SetLocal(int n) { return (nLocal=n);}
    // int AddLocal(int n) { return (nLocal+=n);}
    DATA *operator[](int i) {return Element(i);}

    virtual ~splitStore() {
        // Free the tiles that we allocated
        while (tiles.size()>nTilesReserved) {
            free(tiles.back());
            tiles.pop_back();
        }
        tiles.clear();
    }

};

#endif
