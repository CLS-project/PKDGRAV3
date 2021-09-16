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
#ifndef WHITENOISE_HPP
#define WHITENOISE_HPP

#include "ic/RngStream.h"
#include "core/gridinfo.hpp"
#include "pkd.h"

class NoiseGenerator {
private:
    void pencilNoise(gridinfo::complex_vector_t &pencil,int nGrid,int j, int k);
protected:
    RngStream g;
    float fPhase;
    bool bFixed;
    virtual void update(gridinfo::complex_vector_t &pencil,gridinfo::complex_vector_t &noise,int j,int k);
public:
    explicit NoiseGenerator(unsigned long seed,bool bFixed=false,float fPhase=0);
    virtual ~NoiseGenerator();
    void FillNoise(gridinfo::complex_array_t &K,int nGrid,double *mean=0,double *csq=0);
    };
#endif
