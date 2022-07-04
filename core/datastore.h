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
class dataStore : public dataFields<DATA,FIELD> {
protected:
    DATA *pStore = nullptr;         //!< Pointer to the first element
    int nStore = 0;                 //!< Maximum number of element that can be stored
    int nLocal = 0;                 //!< Current number of elements
public:
    //! Set the size and location of the storage.
    //! \param p Pointer to the block of storage
    //! \param n Maximum number of element
    void setStore(void *p,int n) {
        pStore = static_cast<DATA *>(p);
        nStore = n;
    }
protected:
    DATA *Base() const { return pStore; }
public:
    using dataFields<DATA,FIELD>::Element;
    DATA *Element(int i) const { return this->Element(Base(),i); }
    int FreeStore() const { return nStore; }
    int Local() const { return nLocal; }
    int SetLocal(int n) { return (nLocal=n);}
    int AddLocal(int n) { return (nLocal+=n);}
    DATA *operator[](int i) {return Element(i);}
};

#endif
