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

#ifndef CORE_DATASTORE_H
#define CORE_DATASTORE_H

#include "datafields.h"

//! \brief A generic contiguous array of fixed-length elements.
//!
//! The length of the elements is choosen at run-time before the
//! container is created. This class provides a way to access the
//! individual fields
template<typename DATA,typename FIELD>
class dataStore : public datafields::Fields<DATA,FIELD> {
public:
    using datafields::Fields<DATA,FIELD>::ElementSize;
protected:
    DATA *pStore = nullptr;         //!< Pointer to the first element
    int nStore = 0;                 //!< Maximum number of element that can be stored
    int nLocal = 0;                 //!< Current number of elements
public:
    using value_type = DATA;
    //! Set the size and location of the storage.
    //! \param p Pointer to the block of storage
    //! \param n Maximum number of element
    void setStore(void *p,int n) {
        pStore = static_cast<DATA *>(p);
        nStore = n;
    }
    void setStore(int Lo, int Hi, uint64_t nData=0, void *pData=nullptr) {
        setStore(pData,nData);
    }
protected:
    DATA *Base() const { return pStore; }
public:
    using datafields::Fields<DATA,FIELD>::Element;
    size_t capacity() const { return nStore; }
    size_t size() const { return nLocal; }
    size_t append(int n=1) {
        int iNode = nLocal;
        nLocal += n;
        return iNode;
    }
    void resize(int n) {
        nLocal = n;
    }
    size_t align() {
        if (nLocal&1) ++nLocal;
        return nLocal;
    }

    size_t memory() const {
        return size() * ElementSize();
    }

    DATA *Element(int i) const { return this->Element(Base(),i); }

    DATA *operator[](int i) {return Element(i);}
    operator void *() {return Base();}
};

#endif
