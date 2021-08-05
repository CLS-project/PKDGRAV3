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
#ifndef TRAVERSEPST_H
#define TRAVERSEPST_H
#include "pst.h"
#include <string>

class TraversePST : public mdl::BasicService {
    PST node_pst;
public:
    explicit TraversePST(PST node_pst,int service_id,
	int nInBytes, int nOutBytes, const char *service_name="")
	: BasicService(service_id, nInBytes, nOutBytes, service_name), node_pst(node_pst) {}
    explicit TraversePST(PST node_pst,int service_id,
	int nInBytes, const char *service_name="")
	: BasicService(service_id, nInBytes, 0, service_name), node_pst(node_pst) {}
    explicit TraversePST(PST node_pst,int service_id,const char *service_name="")
	: BasicService(service_id, 0, 0, service_name), node_pst(node_pst) {}
    virtual ~TraversePST() = default;
protected:
    virtual int operator()(int nIn, void *pIn, void *pOut) final;
    virtual int Traverse(PST pst,void *vin,int nIn,void *vout,int nOut);

    virtual int OffNode(PST pst,void *vin,int nIn,void *vout,int nOut) {return Recurse(pst,vin,nIn,vout,nOut);}
    virtual int  AtNode(PST pst,void *vin,int nIn,void *vout,int nOut) {return Recurse(pst,vin,nIn,vout,nOut);}
    virtual int Recurse(PST pst,void *vin,int nIn,void *vout,int nOut);
    virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut) = 0;
    };

#endif
