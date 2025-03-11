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

#ifndef CORE_INTEGERIZE_H
#define CORE_INTEGERIZE_H
#include <cstdint>
#include "blitz/array.h"

class Integerize {
    double d_to_i = 0x80000000u;
    double i_to_d = 1.0 / d_to_i;
public:
    void set_factor(std::uint32_t factor) {
        d_to_i = factor;
        i_to_d = 1.0 / d_to_i;
    }
    double  convert(std::int32_t i) const {return i * i_to_d; }
    template<int n>
    blitz::TinyVector<double,n>  convert(blitz::TinyVector<std::int32_t,n> i) const {return i * i_to_d; }
    std::int32_t convert(double d) const {return d * d_to_i; }
    template<int n>
    blitz::TinyVector<std::int32_t,n> convert(blitz::TinyVector<double,n> d) const {return d * d_to_i; }
};

#endif
