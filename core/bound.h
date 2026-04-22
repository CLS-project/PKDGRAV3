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

namespace detail {
/// @brief Get the type needed to perform long operations (e.g., upper+lower)
template <typename S> struct get_long_type;
template <> struct get_long_type<double> {
    using type = double;
};
template <> struct get_long_type<std::int32_t> {
    using type = std::uint64_t;
};
template <typename S> struct get_width_type;
template <> struct get_width_type<double> {
    using type = double;
};
template <> struct get_width_type<std::int32_t> {
    using type = std::uint32_t;
};
} // namespace detail

/// @brief CRTP base providing shared bounding box operations
/// @tparam T Coordinate type (double or int32_t)
/// @tparam Derived The derived bounding box class
template <typename T, typename Derived>
class BoundBase {
  public:
    using value_type = T;
    using coord_type = blitz::TinyVector<value_type, 3>;
    using size_value_type = typename detail::get_width_type<value_type>::type;
    using size_coord_type = blitz::TinyVector<size_value_type, 3>;
    using long_value_type = typename detail::get_long_type<value_type>::type;
    using long_coord_type = blitz::TinyVector<long_value_type, 3>;

  private:
    BoundBase() = default;
    friend Derived;

    [[nodiscard]] auto derived() const -> const Derived& { return static_cast<const Derived&>(*this); }
    [[nodiscard]] auto derived() -> Derived& { return static_cast<Derived&>(*this); }

  public:
    [[nodiscard]] auto minside() const -> size_value_type { return blitz::min(derived().width()); }
    [[nodiscard]] auto maxside() const -> size_value_type { return blitz::max(derived().width()); }

    /// @brief Returns the dimension (0, 1, or 2) with the largest width.
    /// Ties are broken in favour of the lower index.
    [[nodiscard]] auto maxdim() const -> int {
        auto w = derived().width();
        if (w[0] >= w[1] && w[0] >= w[2]) {
            return 0;
        }
        if (w[1] >= w[2]) {
            return 1;
        }
        return 2;
    }
    /// @brief Returns the dimension (0, 1, or 2) with the smallest width.
    /// Ties are broken in favour of the lower index.
    [[nodiscard]] auto mindim() const -> int {
        auto w = derived().width();
        if (w[0] <= w[1] && w[0] <= w[2]) {
            return 0;
        }
        if (w[1] <= w[2]) {
            return 1;
        }
        return 2;
    }

    [[nodiscard]] auto maxdist(coord_type r) const -> long_value_type {
        long_coord_type x = blitz::abs(derived().center() - r) + derived().apothem();
        return blitz::dot(x, x);
    }
    [[nodiscard]] auto mindist(coord_type r) const -> long_value_type {
        long_coord_type x = blitz::abs(derived().center() - r) - derived().apothem();
        x = blitz::where(x > long_value_type{0}, x, long_value_type{0});
        return blitz::dot(x, x);
    }

    /// @brief Wrap the given coordinate so it lies in this bound
    /// @param r position
    /// @return new position
    /// Note that this assumes that the coordinate will be only slightly
    /// outside the bound, and will not work if it is more than width() away.
    [[nodiscard]] auto wrap(coord_type r) const -> coord_type {
        auto l = derived().lower(), u = derived().upper(), w = derived().width();
        for (auto j = 0; j < 3; ++j) {
            if (r[j] < l[j]) r[j] += w[j];
            else if (r[j] >= u[j]) r[j] -= w[j];
        }
        return r;
    }

    /// @brief Combine two bounds to form a new bound containing both
    /// @param rhs bound to combine with
    /// @return a new bound containing both bounds
    [[nodiscard]] auto combine(const Derived& rhs) const -> Derived {
        return Derived(blitz::min(derived().lower(), rhs.lower()), blitz::max(derived().upper(), rhs.upper()));
    }

    /// @brief Split bound in the middle
    /// @param d dimension to split
    /// @return two new bounds
    [[nodiscard]] auto split(int d) const -> std::pair<Derived, Derived> { return derived().split(d, derived().center(d)); }

    friend auto operator<<(std::ostream& os, const Derived& b) -> std::ostream& {
        os << "[" << b.lower() << "," << b.upper() << "]";
        return os;
    }
};

//****************************************************************************************************

/// @brief Provides an implementation of a bounding box
/// @tparam T Coordinate type (double or int32_t)
/// This keeps track of box with a center coordinate and an apothem (fMax).
template <typename T>
class BoundBaseCenter : public BoundBase<T, BoundBaseCenter<T>> {
    using Base = BoundBase<T, BoundBaseCenter<T>>;

  public:
    using typename Base::coord_type;
    using typename Base::long_coord_type;
    using typename Base::long_value_type;
    using typename Base::size_coord_type;
    using typename Base::size_value_type;
    using typename Base::value_type;

  protected:
    coord_type fCenter;
    coord_type fMax;

  public:
    BoundBaseCenter() = default;
    BoundBaseCenter(const BoundBaseCenter&) = default;
    BoundBaseCenter(coord_type lower, coord_type upper) : fCenter(upper / 2 + lower / 2), fMax((upper - lower) / 2) {}
    template <typename O>
    BoundBaseCenter(const BoundBaseCenter<O>& ibnd, const Integerize& i) : fCenter(i.convert(ibnd.center())), fMax(i.convert(ibnd.apothem())) {}

    [[nodiscard]] auto apothem() const -> coord_type { return fMax; }
    [[nodiscard]] auto apothem(int d) const -> value_type { return fMax[d]; }
    [[nodiscard]] auto width() const -> size_coord_type { return 2 * size_coord_type(fMax); }
    [[nodiscard]] auto width(int d) const -> size_value_type { return 2 * size_value_type(fMax[d]); }
    [[nodiscard]] auto lower() const -> coord_type { return fCenter - fMax; }
    [[nodiscard]] auto lower(int d) const -> value_type { return fCenter[d] - fMax[d]; }
    [[nodiscard]] auto upper() const -> coord_type { return fCenter + fMax; }
    [[nodiscard]] auto upper(int d) const -> value_type { return fCenter[d] + fMax[d]; }
    [[nodiscard]] auto center() const -> coord_type { return fCenter; }
    [[nodiscard]] auto center(int d) const -> value_type { return fCenter[d]; }

    /// @brief Adjust the apothem of the box
    /// @param a the new apothem
    void shrink(coord_type a) { fMax = a; }

    /// @brief Split bound at a specific point
    /// @param d dimension to split
    /// @param at where to split
    /// @return two new bounds
    [[nodiscard]] auto split(int d, value_type at) const -> std::pair<BoundBaseCenter, BoundBaseCenter> {
        BoundBaseCenter l = *this, r = *this;
        l.fMax[d] = (at - lower(d)) / 2;
        l.fCenter[d] = at / 2 + lower(d) / 2;
        r.fMax[d] = (upper(d) - at) / 2;
        r.fCenter[d] = upper(d) / 2 + at / 2;
        return std::make_pair(l, r);
    }
    using Base::split;
};

//****************************************************************************************************

/// @brief Provides an implementation of a bounding box
/// @tparam T Coordinate type (double or int32_t)
/// This keeps track of box with a lower and upper coordinate.
template <typename T>
class BoundBaseMinMax : public BoundBase<T, BoundBaseMinMax<T>> {
    using Base = BoundBase<T, BoundBaseMinMax<T>>;

  public:
    using typename Base::coord_type;
    using typename Base::long_coord_type;
    using typename Base::long_value_type;
    using typename Base::size_coord_type;
    using typename Base::size_value_type;
    using typename Base::value_type;

  protected:
    coord_type fLower;
    coord_type fUpper;

  public:
    BoundBaseMinMax() = default;
    BoundBaseMinMax(const BoundBaseMinMax&) = default;
    BoundBaseMinMax(coord_type lower, coord_type upper) : fLower(lower), fUpper(upper) {}
    template <typename O>
    BoundBaseMinMax(const BoundBaseMinMax<O>& ibnd, const Integerize& i) : fLower(i.convert(ibnd.lower())), fUpper(i.convert(ibnd.upper())) {}

    [[nodiscard]] auto apothem() const -> size_coord_type { return (upper() - lower()) / 2; }
    [[nodiscard]] auto apothem(int d) const -> size_value_type { return (upper(d) - lower(d)) / 2; }
    [[nodiscard]] auto width() const -> size_coord_type { return size_coord_type(upper() - lower()); }
    [[nodiscard]] auto width(int d) const -> size_value_type { return size_value_type(upper(d) - lower(d)); }
    [[nodiscard]] auto lower() const -> coord_type { return fLower; }
    [[nodiscard]] auto lower(int d) const -> value_type { return fLower[d]; }
    [[nodiscard]] auto upper() const -> coord_type { return fUpper; }
    [[nodiscard]] auto upper(int d) const -> value_type { return fUpper[d]; }
    [[nodiscard]] auto center() const -> coord_type { return lower() / 2 + upper() / 2; }
    [[nodiscard]] auto center(int d) const -> value_type { return lower(d) / 2 + upper(d) / 2; }

    /// @brief Adjust the apothem of the box
    /// @param a the new apothem
    void shrink(coord_type a) {
        auto c = center();
        fLower = c - a;
        fUpper = c + a;
    }

    /// @brief Split bound at a specific point
    /// @param d dimension to split
    /// @param at where to split
    /// @return two new bounds
    [[nodiscard]] auto split(int d, value_type at) const -> std::pair<BoundBaseMinMax, BoundBaseMinMax> {
        BoundBaseMinMax l = *this, r = *this;
        l.fUpper[d] = at;
        r.fLower[d] = at;
        return std::make_pair(l, r);
    }
    using Base::split;
};

//****************************************************************************************************

#ifdef BOUND_USES_MINMAX
using Bound = BoundBaseMinMax<double>;
using IntegerBound = BoundBaseMinMax<std::int32_t>;
#else
using Bound = BoundBaseCenter<double>;
using IntegerBound = BoundBaseCenter<std::int32_t>;
#endif
static_assert(std::is_standard_layout_v<Bound>);
static_assert(std::is_standard_layout_v<IntegerBound>);

template <typename T> struct bound_type {};
template <> struct bound_type<double> {
    using type = Bound;
};
template <> struct bound_type<std::int32_t> {
    using type = IntegerBound;
};

#endif
