/*  This file is part of PKDGRAV3 (http://www.pkdgrav.org/).
 *  Copyright (c) 2001-2022 Douglas Potter & Joachim Stadel
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

#ifndef CORE_ELEMENT_H
#define CORE_ELEMENT_H
#include "pkd_config.h"
#include <type_traits>
#include <functional>
#include <memory>

/// @brief Proxy object for a an ELEMENT (value)
/// @tparam ELEMENT The ELEMENT (Particle or Node for example)
template<typename DATA,typename STORE>
class SingleElement {
    template<typename ELEMENT> friend class SinglePointer;
    template<typename ELEMENT> friend class SingleReference;
    template <typename ELEMENT> friend class element_iterator;
public:
    using value_type = DATA;
    using store_type = STORE;
    // using store_type = typename std::conditional<std::is_const_v<DATA>,STORE const,STORE>::type;
protected:
    store_type *element_store;
    value_type *p;
    store_type &store() const {return *element_store;}
    template<typename T,std::enable_if_t<!std::is_array_v<T>,bool> = true>
    T &get(typename store_type::field f) const {return store().template get<T>(p,f);}
    template<typename T,std::enable_if_t<std::rank_v<T> == 1,bool> = true>
    auto get(typename store_type::field f) const {
        return store().template get<blitz::TinyVector<typename std::remove_extent<T>::type,std::extent_v<T>>>(p,f);
    }
    bool have(typename store_type::field f) const {return store().present(f);}
protected:
    SingleElement(const SingleElement &o) = delete;                 // Copy constructor (we dont' allow this)
    SingleElement(SingleElement &&o) = delete;                      // Move constructor (we dont' allow this)
    SingleElement &operator=(const SingleElement &rhs) {            // Copy assignment operator
        memcpy(p,rhs.p,store().ElementSize());                      // ... not that this is protected here
        return *this;                                               // ... because you cannot assign to a "value"
    }
public:
    explicit SingleElement(store_type &store) : element_store(&store), p(nullptr) {}
    SingleElement(store_type &store,int i) : element_store(&store), p(store.Element(i)) {}
    SingleElement(store_type &store,void *p,int i) : element_store(&store), p(store.Element(p,i)) {}
    SingleElement(store_type &store,value_type *p) : element_store(&store), p(p) {}
    SingleElement(const store_type &store,const value_type *p)
        : element_store(const_cast<store_type *>(&store)), p(const_cast<value_type *>(p)) {}
};

template<typename ELEMENT> class SinglePointer;

/// @brief Proxy object for a reference to an ELEMENT
/// @tparam ELEMENT The ELEMENT (Particle or Node for example)
template<typename ELEMENT>
class SingleReference : public ELEMENT {
protected:
    void swap(SingleReference &rhs) noexcept {
        auto n = this->store().ElementSize();
        char buffer[n];
        memcpy(buffer,this->p,n);
        memcpy(this->p,rhs.p,n);
        memcpy(rhs.p,buffer,n);
    }
public:
    using ELEMENT::ELEMENT; // Get all constructors (including the copy constructor)
    using ELEMENT::operator=;
    SinglePointer<ELEMENT> operator&() {
        return SinglePointer<ELEMENT>(this->store(),this->p);
    }
    friend void swap(SingleReference &lhs,SingleReference &rhs) noexcept {
        lhs.swap(rhs);
    }
};

/// @brief Proxy object for a pointer to an ELEMENT
/// @tparam ELEMENT The ELEMENT (Particle or Node for example)
///
/// Semantically this behaves as a pointer, so you can dereference it with "*"
/// to get a "reference" (see above), or use "->" to get the underlying element.
/// You can also convert this to a pointer to the base type (for compatibility).
template<typename ELEMENT>
class SinglePointer {
    template <typename T> friend class element_iterator;
protected:
    using value_type = typename ELEMENT::value_type;
    using store_type = typename ELEMENT::store_type;
    using reference         = SingleReference<ELEMENT>;
    SingleReference<ELEMENT> element;
    ptrdiff_t distance_to(SinglePointer const &other) const {
        return (reinterpret_cast<char *>(other.element.p)
                - reinterpret_cast<char *>(element.p))
               / element.store().ElementSize();
    }
public:
    SinglePointer() = delete; // We really need a particle store
    explicit SinglePointer(store_type &store) : element(store) {}
    SinglePointer(store_type &store,int i) : element(store,i) {}
    SinglePointer(store_type &store,void *p,int i) : element(store,p,i) {}
    SinglePointer(store_type &store,typename ELEMENT::value_type *p) : element(store,p) {}
    SinglePointer(const store_type &store,const typename ELEMENT::value_type *p) : element(store,p) {}
    // reference  operator*()  { return reference(element.store(),element.p); }
    reference &operator*()  { return element; }
    reference *operator->() { return std::addressof(element); }
    operator       typename ELEMENT::value_type *()       {return element.p;}
    operator const typename ELEMENT::value_type *() const {return element.p;}

    SinglePointer(const SinglePointer &o) : element(o.element.store(),o.element.p) {}   // Copy constructor
    SinglePointer(SinglePointer &&o)      : element(o.element.store(),o.element.p) {}   // Move constructor
    SinglePointer &operator=(const SinglePointer &rhs) {                                // Copy assignment operator
        element.p = rhs.element.p;
        element.element_store = rhs.element.element_store;
        return *this;
    }
    SinglePointer &operator=(SinglePointer &&rhs) {                                     // Move assignment operator
        element.p = rhs.element.p;
        element.element_store = rhs.element.element_store;
        return *this;
    }

    // Pointer atithmetic
    SinglePointer &operator+=(int i) {
        element.p = element.store().Element(element.p,i);
        return *this;
    }
    SinglePointer &operator-=(int i) { *this += -i; return *this; }
    SinglePointer &operator++()      { *this += 1; return *this; }
    SinglePointer  operator++(int)   { SinglePointer tmp = *this; ++(*this); return tmp; }
    SinglePointer &operator--()      { *this -= 1; return *this; }
    SinglePointer  operator--(int)   { SinglePointer tmp = *this; --(*this); return tmp; }

    SinglePointer operator+(std::ptrdiff_t i) const {
        SinglePointer tmp = *this; tmp += i; return tmp;
    }
    SinglePointer operator-(std::ptrdiff_t i) const {
        SinglePointer tmp = *this; tmp -= i; return tmp;
    }
    friend ptrdiff_t operator-(const SinglePointer &i,const SinglePointer &i2) {
        return i2.distance_to(i);
    }
    // Pointer comparison
    bool operator== (const SinglePointer &b) const { return element.p == b.element.p; };
    bool operator!= (const SinglePointer &b) const { return element.p != b.element.p; };
    bool operator>  (const SinglePointer &b) const { return element.p >  b.element.p; };
    bool operator>= (const SinglePointer &b) const { return element.p >= b.element.p; };
    bool operator<  (const SinglePointer &b) const { return element.p <  b.element.p; };
    bool operator<= (const SinglePointer &b) const { return element.p <= b.element.p; };
};

template <typename ELEMENT>
class element_iterator {
    typedef element_iterator<ELEMENT> Iterator;
    using store_type        = typename ELEMENT::store_type;
public:
    using iterator_category = std::random_access_iterator_tag;
    using difference_type   = std::ptrdiff_t;
    using value_type        = ELEMENT;
    using reference         = SingleReference<ELEMENT>;
    using pointer           = SinglePointer<ELEMENT>;

    explicit element_iterator(const store_type &store,typename ELEMENT::value_type *p) : element(store,p) {}

    reference &operator*()      { return element.element; }
    pointer   &operator->()     { return element; }

    Iterator &operator+=(int i) { element += i; return *this; }
    Iterator &operator-=(int i) { element -= i; return *this; }
    Iterator &operator++()      { ++element;    return *this; }
    Iterator &operator--()      { --element;    return *this; }
    Iterator  operator++(int)   { Iterator tmp = *this; ++(*this); return tmp; }
    Iterator  operator--(int)   { Iterator tmp = *this; --(*this); return tmp; }

    Iterator  operator+(difference_type i) const { Iterator tmp = *this; tmp += i; return tmp; }
    Iterator  operator-(difference_type i) const { Iterator tmp = *this; tmp -= i; return tmp; }
    friend ptrdiff_t operator-(const Iterator &i1,const Iterator &i2) {
        return i1.element - i2.element;
    }
    friend ptrdiff_t operator-(const SinglePointer<ELEMENT> &p1,const Iterator &i2) {
        return p1 - i2.element;
    }
    friend ptrdiff_t operator-(const Iterator &i1,const SinglePointer<ELEMENT> &p2) {
        return i1.element - p2;
    }

    friend bool operator== (const Iterator &a, const Iterator &b) { return a.element == b.element; };
    friend bool operator!= (const Iterator &a, const Iterator &b) { return a.element != b.element; };
    friend bool operator>  (const Iterator &a, const Iterator &b) { return a.element >  b.element; };
    friend bool operator>= (const Iterator &a, const Iterator &b) { return a.element >= b.element; };
    friend bool operator<  (const Iterator &a, const Iterator &b) { return a.element <  b.element; };
    friend bool operator<= (const Iterator &a, const Iterator &b) { return a.element <= b.element; };

private:
    SinglePointer<ELEMENT> element;
};

#endif