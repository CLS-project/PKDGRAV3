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
#include "TraversePST.h"

class ServiceDistribTopTree : public TraversePST {
public:
    struct input {
	uint32_t uRoot; /* Which root node to use */
	uint32_t nTop;
	// KDN nodes[]; /* Array of tree nodes follow this structure */
	};
    typedef void output;
    explicit ServiceDistribTopTree(PST pst)
	: TraversePST(pst,PST_DISTRIBTOPTREE,
		sizeof(input)
		+ (mdlThreads(pst->mdl)==1?1:2*mdlThreads(pst->mdl)-1)*pkdMaxNodeSize(),
		"DistribTopTree") {}
protected:
    virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut);
    };
