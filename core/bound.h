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

#ifndef BOUND_H
#define BOUND_H
#include "pkd_config.h"
#include <cstdint>
#include <ostream>
#include <utility>
#include <algorithm>
#include <type_traits>
#include "blitz/array.h"
#include "integerize.h"

namespace {
/// @brief Get the type needed to perform long operations (e.g., upper+lower)
template <typename S> struct get_long_type;
template <> struct get_long_type<double>        { using type = double; };
template <> struct get_long_type<std::int32_t>  { using type = std::uint64_t; };
template <typename S> struct get_width_type;
template <> struct get_width_type<double>       { using type = double; };
template <> struct get_width_type<std::int32_t> { using type = std::uint32_t; };
}

/// @brief Provides an implementation of a bounding box
/// @tparam T Coordinate type (double or int32_t)
/// This keeps track of box with a center coordinate and an apothem (fMax).
template<typename T>
class BoundBaseCenter {
    using this_type = BoundBaseCenter<T>;
public:
    // The scalar and vector types (double or int32_t)
    using value_type = T;
    using coord_type = blitz::TinyVector<value_type,3>;
    // Same but for sizes (double or uint32_t)
    using size_value_type = typename get_width_type<value_type>::type;
    using size_coord_type = blitz::TinyVector<size_value_type,3>;
    // Same but for squared values (double or uint64_t)
    using long_value_type = typename get_long_type<value_type>::type;
    using long_coord_type = blitz::TinyVector<long_value_type,3>;
protected:
    coord_type fCenter;
    coord_type fMax;
public:
    BoundBaseCenter() = default;
    BoundBaseCenter(const BoundBaseCenter &) = default;
    BoundBaseCenter(coord_type lower, coord_type upper) : fCenter(upper/2+lower/2),fMax((upper-lower)/2) {}
    template<typename O>
    BoundBaseCenter(const BoundBaseCenter<O> &ibnd,const Integerize &i)
        : fCenter(i.convert(ibnd.center())),fMax(i.convert(ibnd.apothem())) {}

    coord_type      apothem()     const {return fMax;}
    value_type      apothem(int d)const {return fMax[d];}
    size_coord_type width()       const {return 2*size_coord_type(fMax);}
    size_value_type width(int d)  const {return 2*size_value_type(fMax[d]);}
    coord_type      lower()       const {return fCenter-fMax;}
    value_type      lower(int d)  const {return fCenter[d]-fMax[d];}
    coord_type      upper()       const {return fCenter+fMax;}
    value_type      upper(int d)  const {return fCenter[d]+fMax[d];}
    coord_type      center()      const {return fCenter;}
    value_type      center(int d) const {return fCenter[d];}
    size_value_type minside()     const {return 2*size_value_type(blitz::min(apothem()));}
    size_value_type maxside()     const {return 2*size_value_type(blitz::max(apothem()));}
    /// @brief Returns the dimension (0, 1 or 2) in which the "width" of the bound is the largest.
    int             maxdim()      const {return (fMax[0]>fMax[2]) ? (fMax[0]>fMax[1]?0:1) : (fMax[2]>fMax[1]?2:1);}
    /// @brief Returns the dimension (0, 1 or 2) in which the "width" of the bound is the smallest.
    int             mindim()      const {return (fMax[0]<fMax[2]) ? (fMax[0]<fMax[1]?0:1) : (fMax[2]<fMax[1]?2:1);}

    long_value_type maxdist(coord_type r) const {
        long_coord_type x = blitz::abs(center()-r) + apothem();
        return blitz::dot(x,x);
    }
    long_value_type mindist(coord_type r) const {
        long_coord_type x = blitz::abs(center()-r) - apothem();
        x = blitz::where(x>0.0,x,0.0);
        return blitz::dot(x,x);
    }
    /// @brief Wrap the given coordinate so it lies in this bound
    /// @param r position
    /// @return new position
    /// Note that this assumes that the coordinate will be only slightly
    /// outside the bound, and will not work if it is more than width() away.
    coord_type wrap(coord_type r) {
        auto l=lower(), u=upper(), w=width();
        for (auto j=0; j<3; ++j) {
            if      (r[j] <  l[j]) r[j] += w[j];
            else if (r[j] >= u[j]) r[j] -= w[j];
        }
        return r;
    }
    /// @brief Adjust the apothem of the box
    /// @param a the new apothem
    void shrink(coord_type a) {fMax=a;}
    /// @brief Combine two bounds to form a new bound containing both
    /// @tparam D The derived class
    /// @param rhs bound to combine with
    /// @return a new bound containing both bounds
    template<typename D>
    D combine(const D &rhs) const {
        return D(blitz::min(lower(),rhs.lower()),blitz::max(upper(),rhs.upper()));
    }
    /// @brief Split bound at a specific point
    /// @param d dimension to split
    /// @param split where to split
    /// @return two new bounds
    std::pair<this_type,this_type> split(int d,double split) const {
        this_type l = *this, r = *this;
        l.fMax[d]    = (split - lower(d))/2;
        l.fCenter[d] = split/2 + lower(d)/2;
        r.fMax[d]    = (upper(d) - split)/2;
        r.fCenter[d] = upper(d)/2 + split/2;
        return std::make_pair(l,r);
    }
    /// @brief Split bound in the middle
    /// @param d dimension to split
    /// @return two new bounds
    std::pair<this_type,this_type> split(int d) const {
        this_type l = *this, r = *this;
        r.fMax[d] = l.fMax[d] = fMax[d]/2;
        l.fCenter[d] = fCenter[d] - l.fMax[d];
        r.fCenter[d] = fCenter[d] + r.fMax[d];
        return std::make_pair(l,r);
    }

    template<typename D>
    friend std::ostream &operator<<(std::ostream &os, const BoundBaseCenter<D> &b);
};

template<typename D>
inline std::ostream &operator<<(std::ostream &os, const BoundBaseCenter<D> &b) {
    os << "[" << b.lower() << "," << b.upper() << "]";
    return os;
}

//****************************************************************************************************

/// @brief Provides an implementation of a bounding box
/// @tparam T Coordinate type (double or int32_t)
/// This keeps track of box with a lower and upper coordinate.
template<typename T>
class BoundBaseMinMax {
    using this_type = BoundBaseMinMax<T>;
public:
    // The scalar and vector types (double or int32_t)
    using value_type = T;
    using coord_type = blitz::TinyVector<value_type,3>;
    // Same but for sizes (double or uint32_t)
    using size_value_type = typename get_width_type<value_type>::type;
    using size_coord_type = blitz::TinyVector<size_value_type,3>;
    // Same but for squared values (double or uint64_t)
    using long_value_type = typename get_long_type<value_type>::type;
    using long_coord_type = blitz::TinyVector<long_value_type,3>;
protected:
    coord_type fLower;
    coord_type fUpper;
public:
    BoundBaseMinMax() = default;
    BoundBaseMinMax(const BoundBaseMinMax &) = default;
    BoundBaseMinMax(coord_type lower, coord_type upper) : fLower(lower), fUpper(upper) {}
    template<typename O>
    BoundBaseMinMax(const BoundBaseMinMax<O> &ibnd,const Integerize &i)
        : fLower(i.convert(ibnd.lower())),fUpper(i.convert(ibnd.upper())) {}

    size_coord_type apothem()     const {return (upper()-lower())/2;}
    size_value_type apothem(int d)const {return (upper(d)-lower(d))/2;}
    size_coord_type width()       const {return size_coord_type(upper()-lower());}
    size_value_type width(int d)  const {return size_value_type(upper(d)-lower(d));}
    coord_type      lower()       const {return fLower;}
    value_type      lower(int d)  const {return fLower[d];}
    coord_type      upper()       const {return fUpper;}
    value_type      upper(int d)  const {return fUpper[d];}
    coord_type      center()      const {return lower()/2 + upper()/2;}
    value_type      center(int d) const {return lower(d)/2 + upper(d)/2;}
    size_value_type minside()     const {return blitz::min(width());}
    size_value_type maxside()     const {return blitz::max(width());}
    int             maxdim()      const {
        auto fMax = width();
        return (fMax[0]>fMax[2]) ? (fMax[0]>fMax[1]?0:1) : (fMax[2]>fMax[1]?2:1);
    }
    int             mindim()      const {
        auto fMax = width();
        return (fMax[0]<fMax[2]) ? (fMax[0]<fMax[1]?0:1) : (fMax[2]<fMax[1]?2:1);
    }
    long_value_type maxdist(coord_type r) const {
        long_coord_type x = blitz::abs(center()-r) + apothem();
        return blitz::dot(x,x);
    }
    long_value_type mindist(coord_type r) const {
        long_coord_type x = blitz::abs(center()-r) - apothem();
        x = blitz::where(x>0.0,x,0.0);
        return blitz::dot(x,x);
    }
    /// @brief Wrap the given coordinate so it lies in this bound
    /// @param r position
    /// @return new position
    /// Note that this assumes that the coordinate will be only slightly
    /// outside the bound, and will not work if it is more than width() away.
    coord_type wrap(coord_type r) {
        auto l=lower(), u=upper(), w=width();
        for (auto j=0; j<3; ++j) {
            if      (r[j] <  l[j]) r[j] += w[j];
            else if (r[j] >= u[j]) r[j] -= w[j];
        }
        return r;
    }
    /// @brief Adjust the apothem of the box
    /// @param a the new apothem
    void shrink(coord_type a) {
        auto c = center();
        fLower = c - a;
        fUpper = c + a;
    }
    /// @brief Combine two bounds to form a new bound containing both
    /// @tparam D The derived class
    /// @param rhs bound to combine with
    /// @return a new bound containing both bounds
    template<typename D>
    D combine(const D &rhs) const {
        return D(blitz::min(lower(),rhs.lower()),blitz::max(upper(),rhs.upper()));
    }
    /// @brief Split bound at a specific point
    /// @param d dimension to split
    /// @param split where to split
    /// @return two new bounds
    std::pair<this_type,this_type> split(int d,double split) const {
        this_type l = *this, r = *this;
        l.fUpper[d]  = split;
        r.fLower[d] = split;
        return std::make_pair(l,r);
    }
    /// @brief Split bound in the middle
    /// @param d dimension to split
    /// @return two new bounds
    std::pair<this_type,this_type> split(int d) const {
        this_type l = *this, r = *this;
        return split(d,center(d));
    }

    template<typename D>
    friend std::ostream &operator<<(std::ostream &os, const BoundBaseMinMax<D> &b);
};

template<typename D>
inline std::ostream &operator<<(std::ostream &os, const BoundBaseMinMax<D> &b) {
    os << "[" << b.lower() << "," << b.upper() << "]";
    return os;
}

//****************************************************************************************************

#ifdef BOUND_USES_MINMAX
    using Bound = BoundBaseMinMax<double>;
    using IntegerBound = BoundBaseMinMax<std::int32_t>;
#else
    using Bound = BoundBaseCenter<double>;
    using IntegerBound = BoundBaseCenter<std::int32_t>;
#endif
static_assert(std::is_standard_layout<Bound>());
static_assert(std::is_standard_layout<IntegerBound>());

template<typename T> struct bound_type {};
template<> struct bound_type<double>        {using type = Bound;};
template<> struct bound_type<std::int32_t>  {using type = IntegerBound;};

#endif
