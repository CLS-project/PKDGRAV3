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

// Base class for particle selection (setting bMark)
// The output is always the count of the particles selected
class ServiceSelect : public TraverseCombinePST {
public:
    typedef uint64_t output;
    explicit ServiceSelect(PST pst,int service_id,int nInBytes, const char *service_name="")
	: TraverseCombinePST(pst,service_id,nInBytes,sizeof(output),service_name) {}
    explicit ServiceSelect(PST pst,int service_id,const char *service_name="")
	: TraverseCombinePST(pst,service_id,0,sizeof(output),service_name) {}
protected:
    virtual int Combine(void *vout,void *vout2) final;
    virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut) override = 0;
    };

class ServiceCountSelected : public ServiceSelect {
public:
    typedef void input;
    explicit ServiceCountSelected(PST pst)
	: ServiceSelect(pst,PST_COUNTSELECTED,"CountSelected") {}
protected:
    virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut) override;
    };

class ServiceSelBlackholes : public ServiceSelect {
public:
    struct input {
	bool setIfTrue;
	bool clearIfFalse;
	input() = default;
	input(bool setIfTrue,bool clearIfFalse)
	    : setIfTrue(setIfTrue), clearIfFalse(clearIfFalse) {}
	};
    explicit ServiceSelBlackholes(PST pst)
	: ServiceSelect(pst,PST_SELBLACKHOLES,sizeof(input),"SelectBlackholes") {}
protected:
    virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut) override;
    };

class ServiceSelSpecies : public ServiceSelect {
public:
    struct input {
	uint64_t mSpecies;
	bool setIfTrue;
	bool clearIfFalse;
	input() = default;
	input(uint64_t mSpecies,bool setIfTrue=true,bool clearIfFalse=true)
	    : mSpecies(mSpecies), setIfTrue(setIfTrue), clearIfFalse(clearIfFalse) {}
	};
    explicit ServiceSelSpecies(PST pst)
	: ServiceSelect(pst,PST_SELSPECIES,sizeof(input),"SelSpecies") {}
protected:
    virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut) override;
    };

class ServiceSelGroup : public ServiceSelect {
public:
    struct input {
	uint64_t iGroup;
	bool setIfTrue;
	bool clearIfFalse;
	input() = default;
	input(uint64_t iGroup,bool setIfTrue=true,bool clearIfFalse=true)
	    : iGroup(iGroup), setIfTrue(setIfTrue), clearIfFalse(clearIfFalse) {}
	};
    explicit ServiceSelGroup(PST pst)
	: ServiceSelect(pst,PST_SELGROUP,sizeof(input),"SelGroup") {}
protected:
    virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut) override;
    };

class ServiceSelMass : public ServiceSelect {
public:
    struct input {
	double dMinMass;
	double dMaxMass;
	bool setIfTrue;
	bool clearIfFalse;
	input() = default;
	input(double dMinMass,double dMaxMass,bool setIfTrue=true,bool clearIfFalse=true)
	    : dMinMass(dMinMass), dMaxMass(dMaxMass), setIfTrue(setIfTrue), clearIfFalse(clearIfFalse) {}
	};
    explicit ServiceSelMass(PST pst)
	: ServiceSelect(pst,PST_SELMASS,sizeof(input),"SelMass") {}
protected:
    virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut) override;
    };

class ServiceSelPhaseDensity : public ServiceSelect {
public:
    struct input {
	double dMinDensity;
	double dMaxDensity;
	bool setIfTrue;
	bool clearIfFalse;
	input() = default;
	input(double dMinDensity,double dMaxDensity,bool setIfTrue=true,bool clearIfFalse=true)
	    : dMinDensity(dMinDensity), dMaxDensity(dMaxDensity), setIfTrue(setIfTrue), clearIfFalse(clearIfFalse) {}
	};
    explicit ServiceSelPhaseDensity(PST pst)
	: ServiceSelect(pst,PST_SELPHASEDENSITY,sizeof(input),"SelPhaseDensity") {}
protected:
    virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut) override;
    };

class ServiceSelById : public ServiceSelect {
public:
    struct input {
	uint64_t idStart;
	uint64_t idEnd;
	bool setIfTrue;
	bool clearIfFalse;
	input() = default;
	input(uint64_t idStart, uint64_t idEnd,bool setIfTrue=true,bool clearIfFalse=true)
	    : idStart(idStart), idEnd(idEnd), setIfTrue(setIfTrue), clearIfFalse(clearIfFalse) {}
	};
    explicit ServiceSelById(PST pst)
	: ServiceSelect(pst,PST_SELBYID,sizeof(input),"SelById") {}
protected:
    virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut) override;
    };

class ServiceSelBox : public ServiceSelect {
public:
    struct input {
	double dCenter[3];
	double dSize[3];
	bool setIfTrue;
	bool clearIfFalse;
	input() = default;
	input(const double *dCenter,const double *dSize,bool setIfTrue=true,bool clearIfFalse=true)
	    : dCenter{dCenter[0],dCenter[1],dCenter[2]}, dSize{dSize[0],dSize[1],dSize[2]},
	      setIfTrue(setIfTrue), clearIfFalse(clearIfFalse) {}
	input(double p0, double p1,double p2,double s0,double s1,double s2,bool setIfTrue=true,bool clearIfFalse=true)
	    : dCenter{p0,p1,p2}, dSize{s0,s1,s2}, setIfTrue(setIfTrue), clearIfFalse(clearIfFalse) {}
	};
    explicit ServiceSelBox(PST pst)
	: ServiceSelect(pst,PST_SELBOX,sizeof(input),"SelBox") {}
protected:
    virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut) override;
    };

class ServiceSelSphere : public ServiceSelect {
public:
    struct input {
	double r[3];
	double dRadius;
	bool setIfTrue;
	bool clearIfFalse;
	input() = default;
	input(const double *r,double dRadius,bool setIfTrue=true,bool clearIfFalse=true)
	    : r{r[0],r[1],r[2]}, dRadius(dRadius), setIfTrue(setIfTrue), clearIfFalse(clearIfFalse) {}
	input(double r0, double r1,double r2,double dRadius,bool setIfTrue=true,bool clearIfFalse=true)
	    : r{r0,r1,r2}, dRadius(dRadius), setIfTrue(setIfTrue), clearIfFalse(clearIfFalse) {}
	};
    explicit ServiceSelSphere(PST pst)
	: ServiceSelect(pst,PST_SELSPHERE,sizeof(input),"SelSphere") {}
protected:
    virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut) override;
    };

class ServiceSelCylinder : public ServiceSelect {
public:
    struct input {
	double dP1[3];
	double dP2[3];
	double dRadius;
	bool setIfTrue;
	bool clearIfFalse;
	input() = default;
	input(const double *dP1,const double *dP2,double dRadius,bool setIfTrue=true,bool clearIfFalse=true)
	    : dP1{dP1[0],dP1[1],dP1[2]}, dP2{dP2[0],dP2[1],dP2[2]}, dRadius(dRadius),
	      setIfTrue(setIfTrue), clearIfFalse(clearIfFalse) {}
	input(double ra0, double ra1,double ra2,double rb0,double rb1,double rb2,double dRadius,
	      bool setIfTrue=true,bool clearIfFalse=true)
	    : dP1{ra0,ra1,ra2}, dP2{rb0,rb1,rb2}, dRadius(dRadius), setIfTrue(setIfTrue), clearIfFalse(clearIfFalse) {}
	};
    explicit ServiceSelCylinder(PST pst)
	: ServiceSelect(pst,PST_SELCYLINDER,sizeof(input),"SelCylinder") {}
protected:
    virtual int Service(PST pst,void *vin,int nIn,void *vout,int nOut) override;
    };
