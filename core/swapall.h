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
#include "pst.h"

class ServiceSwapAll : public mdl::BasicService {
    PST node_pst;
public:
    struct input {
    	int idSwap;
	input() = default;
	input(int idSwap) : idSwap(idSwap) {}
	};
    typedef void output;
    explicit ServiceSwapAll(PST node_pst)
	: BasicService(PST_SWAPALL, sizeof(input), "SwapAll"), node_pst(node_pst) {}
protected:
    virtual int operator()(int nIn, void *pIn, void *pOut) override;
    };
