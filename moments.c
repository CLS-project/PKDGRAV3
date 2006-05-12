#include <stdio.h>
#include <math.h>
#include "moments.h"


/*
 ** This function adds the complete moment ma to the complete moment mc
 */
void momAddMomc(MOMC *mc,MOMC *ma)
{
	mc->m += ma->m;
	mc->xx += ma->xx;
	mc->yy += ma->yy;
	mc->xy += ma->xy;
	mc->xz += ma->xz;
	mc->yz += ma->yz;
	mc->xxx += ma->xxx;
	mc->xyy += ma->xyy;
	mc->xxy += ma->xxy;
	mc->yyy += ma->yyy;
	mc->xxz += ma->xxz;
	mc->yyz += ma->yyz;
	mc->xyz += ma->xyz;
	mc->xxxx += ma->xxxx;
	mc->xyyy += ma->xyyy;
	mc->xxxy += ma->xxxy;
	mc->yyyy += ma->yyyy;
	mc->xxxz += ma->xxxz;
	mc->yyyz += ma->yyyz;
	mc->xxyy += ma->xxyy;
	mc->xxyz += ma->xxyz;
	mc->xyyz += ma->xyyz;
	mc->zz += ma->zz;
	mc->xzz += ma->xzz;
	mc->yzz += ma->yzz;
	mc->zzz += ma->zzz;
	mc->xxzz += ma->xxzz;
	mc->xyzz += ma->xyzz;
	mc->xzzz += ma->xzzz;
	mc->yyzz += ma->yyzz;
	mc->yzzz += ma->yzzz;
	mc->zzzz += ma->zzzz;
	}


/*
 ** This function adds the reduced moment ma to the reduced moment mc
 */
void momAddMomr(MOMR *mr,MOMR *ma)
{
	mr->m += ma->m;
	mr->xx += ma->xx;
	mr->yy += ma->yy;
	mr->xy += ma->xy;
	mr->xz += ma->xz;
	mr->yz += ma->yz;
	mr->xxx += ma->xxx;
	mr->xyy += ma->xyy;
	mr->xxy += ma->xxy;
	mr->yyy += ma->yyy;
	mr->xxz += ma->xxz;
	mr->yyz += ma->yyz;
	mr->xyz += ma->xyz;
	mr->xxxx += ma->xxxx;
	mr->xyyy += ma->xyyy;
	mr->xxxy += ma->xxxy;
	mr->yyyy += ma->yyyy;
	mr->xxxz += ma->xxxz;
	mr->yyyz += ma->yyyz;
	mr->xxyy += ma->xxyy;
	mr->xxyz += ma->xxyz;
	mr->xyyz += ma->xyyz;
	}


/*
 ** This function multiply-adds the complete moment ma
 */
void momMulAddMomc(MOMC *mc,momFloat m,MOMC *ma)
{
	mc->m += m*ma->m;
	mc->xx += m*ma->xx;
	mc->yy += m*ma->yy;
	mc->xy += m*ma->xy;
	mc->xz += m*ma->xz;
	mc->yz += m*ma->yz;
	mc->xxx += m*ma->xxx;
	mc->xyy += m*ma->xyy;
	mc->xxy += m*ma->xxy;
	mc->yyy += m*ma->yyy;
	mc->xxz += m*ma->xxz;
	mc->yyz += m*ma->yyz;
	mc->xyz += m*ma->xyz;
	mc->xxxx += m*ma->xxxx;
	mc->xyyy += m*ma->xyyy;
	mc->xxxy += m*ma->xxxy;
	mc->yyyy += m*ma->yyyy;
	mc->xxxz += m*ma->xxxz;
	mc->yyyz += m*ma->yyyz;
	mc->xxyy += m*ma->xxyy;
	mc->xxyz += m*ma->xxyz;
	mc->xyyz += m*ma->xyyz;
	mc->zz += m*ma->zz;
	mc->xzz += m*ma->xzz;
	mc->yzz += m*ma->yzz;
	mc->zzz += m*ma->zzz;
	mc->xxzz += m*ma->xxzz;
	mc->xyzz += m*ma->xyzz;
	mc->xzzz += m*ma->xzzz;
	mc->yyzz += m*ma->yyzz;
	mc->yzzz += m*ma->yzzz;
	mc->zzzz += m*ma->zzzz;
	}


/*
 ** This function multiply-adds the reduced moment ma
 */
void momMulAddMomr(MOMR *mr,momFloat m,MOMR *ma)
{
	mr->m += m*ma->m;
	mr->xx += m*ma->xx;
	mr->yy += m*ma->yy;
	mr->xy += m*ma->xy;
	mr->xz += m*ma->xz;
	mr->yz += m*ma->yz;
	mr->xxx += m*ma->xxx;
	mr->xyy += m*ma->xyy;
	mr->xxy += m*ma->xxy;
	mr->yyy += m*ma->yyy;
	mr->xxz += m*ma->xxz;
	mr->yyz += m*ma->yyz;
	mr->xyz += m*ma->xyz;
	mr->xxxx += m*ma->xxxx;
	mr->xyyy += m*ma->xyyy;
	mr->xxxy += m*ma->xxxy;
	mr->yyyy += m*ma->yyyy;
	mr->xxxz += m*ma->xxxz;
	mr->yyyz += m*ma->yyyz;
	mr->xxyy += m*ma->xxyy;
	mr->xxyz += m*ma->xxyz;
	mr->xyyz += m*ma->xyyz;
	}


/*
 ** This function subtracts the complete moment ma from the complete moment mc
 */
void momSubMomc(MOMC *mc,MOMC *ma)
{
	mc->m -= ma->m;
	mc->xx -= ma->xx;
	mc->yy -= ma->yy;
	mc->xy -= ma->xy;
	mc->xz -= ma->xz;
	mc->yz -= ma->yz;
	mc->xxx -= ma->xxx;
	mc->xyy -= ma->xyy;
	mc->xxy -= ma->xxy;
	mc->yyy -= ma->yyy;
	mc->xxz -= ma->xxz;
	mc->yyz -= ma->yyz;
	mc->xyz -= ma->xyz;
	mc->xxxx -= ma->xxxx;
	mc->xyyy -= ma->xyyy;
	mc->xxxy -= ma->xxxy;
	mc->yyyy -= ma->yyyy;
	mc->xxxz -= ma->xxxz;
	mc->yyyz -= ma->yyyz;
	mc->xxyy -= ma->xxyy;
	mc->xxyz -= ma->xxyz;
	mc->xyyz -= ma->xyyz;
	mc->zz -= ma->zz;
	mc->xzz -= ma->xzz;
	mc->yzz -= ma->yzz;
	mc->zzz -= ma->zzz;
	mc->xxzz -= ma->xxzz;
	mc->xyzz -= ma->xyzz;
	mc->xzzz -= ma->xzzz;
	mc->yyzz -= ma->yyzz;
	mc->yzzz -= ma->yzzz;
	mc->zzzz -= ma->zzzz;
	}


/*
 ** This function subtracts the reduced moment ma from the reduced moment mc
 */
void momSubMomr(MOMR *mr,MOMR *ma)
{
	mr->m -= ma->m;
	mr->xx -= ma->xx;
	mr->yy -= ma->yy;
	mr->xy -= ma->xy;
	mr->xz -= ma->xz;
	mr->yz -= ma->yz;
	mr->xxx -= ma->xxx;
	mr->xyy -= ma->xyy;
	mr->xxy -= ma->xxy;
	mr->yyy -= ma->yyy;
	mr->xxz -= ma->xxz;
	mr->yyz -= ma->yyz;
	mr->xyz -= ma->xyz;
	mr->xxxx -= ma->xxxx;
	mr->xyyy -= ma->xyyy;
	mr->xxxy -= ma->xxxy;
	mr->yyyy -= ma->yyyy;
	mr->xxxz -= ma->xxxz;
	mr->yyyz -= ma->yyyz;
	mr->xxyy -= ma->xxyy;
	mr->xxyz -= ma->xxyz;
	mr->xyyz -= ma->xyyz;
	}


/*
 ** This function calculates a complete multipole from a single
 ** particle at position <x,y,z> from the center of mass.
 ** The strange order of evaluation reduces the number of 
 ** multiplications to a minimum.
 ** <x,y,z> := d := r(particle) - rcm.
 **
 ** OpCount (*,+) = (34,0)
 **
 */
void momMakeMomc(MOMC *mc,momFloat m,momFloat x,momFloat y,momFloat z)
{
	momFloat tx,ty,tz;

	mc->m = m;
	/*
	 ** Calculate the Quadrupole Moment.
	 */
	tx = m*x;
	ty = m*y;
	mc->xy = tx*y;
	mc->xz = tx*z;
	mc->yz = ty*z;
	tx *= x;
	ty *= y;
	tz = m*z*z;
	mc->xx = tx;
	mc->yy = ty;
	mc->zz = tz;
	/*
	 ** Calculate the Octopole Moment.
	 */
	mc->xxy = tx*y;
	mc->xxz = tx*z;
	mc->yyz = ty*z;
	mc->xyy = ty*x;
	mc->xzz = tz*x;
	mc->yzz = tz*y;
	mc->xyz = mc->xy*z;
	tx *= x;
	ty *= y;
	tz *= z;
 	mc->xxx = tx;
	mc->yyy = ty;
	mc->zzz = tz;
	/*
	 ** Calculate the Hexadecapole Moment.
	 */
	mc->xxxx = tx*x;
	mc->xxxy = tx*y;
	mc->xxxz = tx*z;
	mc->xyyy = ty*x;
	mc->yyyy = ty*y;
	mc->yyyz = ty*z;
	mc->xzzz = tz*x;
	mc->yzzz = tz*y;
	mc->zzzz = tz*z;
	mc->xxyy = mc->xxy*y;
	mc->xxyz = mc->xxy*z;
	mc->xyyz = mc->yyz*x;
	mc->yyzz = mc->yyz*z;
	mc->xxzz = mc->xzz*x;
	mc->xyzz = mc->xzz*y;
	}


/*
 ** This function calculates a reduced multipole from a single
 ** particle at position <x,y,z> from the center of mass.
 ** The strange order of evaluation reduces the number of 
 ** multiplications to a minimum.
 ** <x,y,z> := d := r(particle) - rcm.
 ** returns: d^2
 **
 ** OpCount (*,+) = (43,18) = ~60
 */
momFloat momMakeMomr(MOMR *mr,momFloat m,momFloat x,momFloat y,momFloat z)
{
	momFloat tx,ty,t,dx,dy;
	momFloat x2 = x*x;
	momFloat y2 = y*y;
	momFloat d2 = x2 + y2 + z*z;

	mr->m = m;
	/*
	 ** Calculate the Quadrupole Moment.
	 */
	tx = m*x;
	ty = m*y;
	mr->xy = tx*y;
	mr->xz = tx*z;
	mr->yz = ty*z;
	tx *= x;
	ty *= y;
	m *= d2;
	t = m/3;
	mr->xx = tx - t;
	mr->yy = ty - t;
	/*
	 ** Calculate the Octopole Moment.
	 */
	t = 0.2*m;
	dx = tx - t;
	dy = ty - t;
	mr->xxy = dx*y;
	mr->xxz = dx*z;
	mr->yyz = dy*z;
	mr->xyy = dy*x;
	mr->xyz = mr->xy*z;
	t *= 3;
 	mr->xxx = (tx - t)*x;
	mr->yyy = (ty - t)*y;
	/*
	 ** Calculate the Hexadecapole Moment.
	 */
	t = m/7;
	mr->xxyz = (tx - t)*y*z;
	mr->xyyz = (ty - t)*x*z;
	dx = (tx - 3*t)*x;
	dy = (ty - 3*t)*y;
	mr->xxxy = dx*y;
	mr->xxxz = dx*z;
	mr->xyyy = dy*x;
	mr->yyyz = dy*z;
	dx = t*(x2 - 0.1*d2);
	dy = t*(y2 - 0.1*d2);
	mr->xxxx = tx*x2 - 6*dx;
	mr->yyyy = ty*y2 - 6*dy;
	mr->xxyy = tx*y2 - dx - dy;

	return(d2);
	}


/*
 ** This function calculates a reduced multipole from a single
 ** particle at position <x,y,z> from the center of mass.
 ** This is the "straight forward" implementation which we 
 ** used in the original version of PKDGRAV. It remains a good
 ** test of more peculiar looking code.
 **
 ** <x,y,z> := d := r(particle) - rcm.
 ** 
 ** OpCount (*,+) = (115,20) = ~135
 */
void momOldMakeMomr(MOMR *mr,momFloat m,momFloat x,momFloat y,momFloat z)
{
	momFloat d2 = x*x + y*y + z*z;
	
	mr->xxxx = m*(x*x*x*x - 6.0/7.0*d2*(x*x - 0.1*d2));
	mr->xyyy = m*(x*y*y*y - 3.0/7.0*d2*x*y);
	mr->xxxy = m*(x*x*x*y - 3.0/7.0*d2*x*y);
	mr->yyyy = m*(y*y*y*y - 6.0/7.0*d2*(y*y - 0.1*d2));
	mr->xxxz = m*(x*x*x*z - 3.0/7.0*d2*x*z);
	mr->yyyz = m*(y*y*y*z - 3.0/7.0*d2*y*z);
	mr->xxyy = m*(x*x*y*y - 1.0/7.0*d2*(x*x + y*y - 0.2*d2));
	mr->xxyz = m*(x*x*y*z - 1.0/7.0*d2*y*z);
	mr->xyyz = m*(x*y*y*z - 1.0/7.0*d2*x*z);
	/*
	 ** Calculate reduced octopole moment...
	 */
	mr->xxx = m*(x*x*x - 0.6*d2*x);
	mr->xyy = m*(x*y*y - 0.2*d2*x);
	mr->xxy = m*(x*x*y - 0.2*d2*y);
	mr->yyy = m*(y*y*y - 0.6*d2*y);
	mr->xxz = m*(x*x*z - 0.2*d2*z);
	mr->yyz = m*(y*y*z - 0.2*d2*z);
	mr->xyz = m*x*y*z;
	/*
	 ** Calculate quadrupole moment...
	 */
	mr->xx = m*(x*x - 1.0/3.0*d2);
	mr->yy = m*(y*y - 1.0/3.0*d2);
	mr->xy = m*x*y;
	mr->xz = m*x*z;
	mr->yz = m*y*z;
	mr->m = m;
	}


/*
 ** This function shifts a complete multipole (MOMC) to a new center of mass.
 ** <x,y,z> := d := rcm(old) - rcm(new).
 **
 ** OpCount ShiftMomc   (*,+) = (111,84)
 **         MakeMomc    (*,+) = (34,0)
 **         MulAddMomc  (*,+) = (32,32)
 **         Total       (*,+) = (177,116) = 293
 */
void momShiftMomc(MOMC *m,momFloat x,momFloat y,momFloat z)
{
	MOMC f;

	momMakeMomc(&f,1,x,y,z);
	/*
	 ** Shift the Hexadecapole.
	 */
	m->xxxx += 4*m->xxx*x + 6*m->xx*f.xx;
	m->yyyy += 4*m->yyy*y + 6*m->yy*f.yy;
	m->zzzz += 4*m->zzz*z + 6*m->zz*f.zz;
	m->xyyy += m->yyy*x + 3*(m->xyy*y + m->yy*f.xy + m->xy*f.yy);
	m->xxxy += m->xxx*y + 3*(m->xxy*x + m->xx*f.xy + m->xy*f.xx);
	m->xxxz += m->xxx*z + 3*(m->xxz*x + m->xx*f.xz + m->xz*f.xx);
	m->yyyz += m->yyy*z + 3*(m->yyz*y + m->yy*f.yz + m->yz*f.yy);
	m->xzzz += m->zzz*x + 3*(m->xzz*z + m->zz*f.xz + m->xz*f.zz);
	m->yzzz += m->zzz*y + 3*(m->yzz*z + m->zz*f.yz + m->yz*f.zz);
	m->xxyy += 2*(m->xxy*y + m->xyy*x) + m->xx*f.yy + m->yy*f.xx + 4*m->xy*f.xy;
	m->xxzz += 2*(m->xxz*z + m->xzz*x) + m->xx*f.zz + m->zz*f.xx + 4*m->xz*f.xz;
	m->yyzz += 2*(m->yyz*z + m->yzz*y) + m->yy*f.zz + m->zz*f.yy + 4*m->yz*f.yz;
	m->xxyz += m->xxy*z + m->xxz*y + m->xx*f.yz + m->yz*f.xx + 2*(m->xyz*x + m->xy*f.xz + m->xz*f.xy);
	m->xyyz += m->xyy*z + m->yyz*x + m->yy*f.xz + m->xz*f.yy + 2*(m->xyz*y + m->xy*f.yz + m->yz*f.xy);
	m->xyzz += m->yzz*x + m->xzz*y + m->zz*f.xy + m->xy*f.zz + 2*(m->xyz*z + m->yz*f.xz + m->xz*f.yz);
	/*
	 ** Now shift the Octopole.
	 */
	m->xxx += 3*m->xx*x;
	m->yyy += 3*m->yy*y;
	m->zzz += 3*m->zz*z;
	m->xyy += 2*m->xy*y + m->yy*x;
	m->xzz += 2*m->xz*z + m->zz*x;
	m->yzz += 2*m->yz*z + m->zz*y;
	m->xxy += 2*m->xy*x + m->xx*y; 
	m->xxz += 2*m->xz*x + m->xx*z;
	m->yyz += 2*m->yz*y + m->yy*z;
	m->xyz += m->xy*z + m->xz*y + m->yz*x;
	/*
	 ** Now deal with the monopole terms.
	 */
	f.m = 0;
	momMulAddMomc(m,m->m,&f);
	}


/*
 ** This function shifts a reduced multipole (MOMR) to a new center of mass.
 ** <x,y,z> := d := rcm(old) - rcm(new).
 **
 ** OpCount ShiftMomr  (*,+) = (128,111)
 **         MakeMomr   (*,+) = (43,18)
 **         MulAddMomr (*,+) = (22,22)
 **         Total      (*,+) = (193,151) = 344
 */
void momShiftMomr(MOMR *m,momFloat x,momFloat y,momFloat z)
{
	MOMR f;
	momFloat t,tx,ty,tz,txx,tyy,txy,tyz,txz;
	const momFloat twosevenths = 2.0/7.0;

	momMakeMomr(&f,1,x,y,z);
	/*
	 ** Calculate the correction terms.
	 */
	tx = 0.4*(m->xx*x + m->xy*y + m->xz*z);
	ty = 0.4*(m->xy*x + m->yy*y + m->yz*z);
	tz = 0.4*(m->xz*x + m->yz*y - (m->xx + m->yy)*z);
	t = tx*x + ty*y + tz*z;
	txx = twosevenths*(m->xxx*x + m->xxy*y + m->xxz*z + 2*(m->xx*f.xx + m->xy*f.xy + m->xz*f.xz) - 0.5*t);
	tyy = twosevenths*(m->xyy*x + m->yyy*y + m->yyz*z + 2*(m->xy*f.xy + m->yy*f.yy + m->yz*f.yz) - 0.5*t);
	txy = twosevenths*(m->xxy*x + m->xyy*y + m->xyz*z + m->xy*(f.xx + f.yy) + (m->xx + m->yy)*f.xy + m->yz*f.xz + m->xz*f.yz);
	tyz = twosevenths*(m->xyz*x + m->yyz*y - (m->xxy + m->yyy)*z - m->yz*f.xx - m->xx*f.yz + m->xz*f.xy + m->xy*f.xz);
	txz = twosevenths*(m->xxz*x + m->xyz*y - (m->xxx + m->xyy)*z - m->xz*f.yy - m->yy*f.xz + m->yz*f.xy + m->xy*f.yz);
	/*
	 ** Shift the Hexadecapole.
	 */
	m->xxxx += 4*m->xxx*x + 6*(m->xx*f.xx - txx);
	m->yyyy += 4*m->yyy*y + 6*(m->yy*f.yy - tyy);
	m->xyyy += m->yyy*x + 3*(m->xyy*y + m->yy*f.xy + m->xy*f.yy - txy);
	m->xxxy += m->xxx*y + 3*(m->xxy*x + m->xx*f.xy + m->xy*f.xx - txy);
	m->xxxz += m->xxx*z + 3*(m->xxz*x + m->xx*f.xz + m->xz*f.xx - txz);
	m->yyyz += m->yyy*z + 3*(m->yyz*y + m->yy*f.yz + m->yz*f.yy - tyz);
	m->xxyy += 2*(m->xxy*y + m->xyy*x) + m->xx*f.yy + m->yy*f.xx + 4*m->xy*f.xy - txx - tyy;
	m->xxyz += m->xxy*z + m->xxz*y + m->xx*f.yz + m->yz*f.xx + 2*(m->xyz*x + m->xy*f.xz + m->xz*f.xy) - tyz;
	m->xyyz += m->xyy*z + m->yyz*x + m->yy*f.xz + m->xz*f.yy + 2*(m->xyz*y + m->xy*f.yz + m->yz*f.xy) - txz;
	/*
	 ** Now shift the Octopole.
	 */
	m->xxx += 3*(m->xx*x - tx);
	m->xyy += 2*m->xy*y + m->yy*x - tx;
	m->yyy += 3*(m->yy*y - ty);
	m->xxy += 2*m->xy*x + m->xx*y - ty;
	m->xxz += 2*m->xz*x + m->xx*z - tz;
	m->yyz += 2*m->yz*y + m->yy*z - tz;
	m->xyz += m->xy*z + m->xz*y + m->yz*x;
	/*
	 ** Now deal with the monopole terms.
	 */
	f.m = 0;
	momMulAddMomr(m,m->m,&f);
	}


/*
 ** This function converts a complete multipole (MOMC) to a reduced one (MOMR).
 */
void momReduceMomc(MOMC *mc,MOMR *mr)
{
	momFloat  t,tx,ty,tz,txx,txy,txz,tyy,tyz,tzz;

	/*
	 ** First reduce Hexadecapole.
	 */
	txx = (mc->xxxx + mc->xxyy + mc->xxzz)/7;
	txy = (mc->xxxy + mc->xyyy + mc->xyzz)/7;
	txz = (mc->xxxz + mc->xyyz + mc->xzzz)/7;
	tyy = (mc->xxyy + mc->yyyy + mc->yyzz)/7;
	tyz = (mc->xxyz + mc->yyyz + mc->yzzz)/7;
	tzz = (mc->xxzz + mc->yyzz + mc->zzzz)/7;
	t = 0.1*(txx + tyy + tzz);
	mr->xxxx = mc->xxxx - 6*(txx - t);
	mr->xyyy = mc->xyyy - 3*txy;
	mr->xxxy = mc->xxxy - 3*txy;
	mr->yyyy = mc->yyyy - 6*(tyy - t);
	mr->xxxz = mc->xxxz - 3*txz;
	mr->yyyz = mc->yyyz - 3*tyz;
	mr->xxyy = mc->xxyy - (txx + tyy - 2*t);
	mr->xxyz = mc->xxyz - tyz;
	mr->xyyz = mc->xyyz - txz;	
	/*
	 ** Now reduce the Octopole.
	 */
	tx = (mc->xxx + mc->xyy + mc->xzz)/5;
	ty = (mc->xxy + mc->yyy + mc->yzz)/5;
	tz = (mc->xxz + mc->yyz + mc->zzz)/5;
	mr->xxx = mc->xxx - 3*tx;
	mr->xyy = mc->xyy - tx;
	mr->xxy = mc->xxy - ty;
	mr->yyy = mc->yyy - 3*ty;
	mr->xxz = mc->xxz - tz;
	mr->yyz = mc->yyz - tz;
	mr->xyz = mc->xyz;
	/*
	 ** Now reduce the Quadrupole.
	 */
	t = (mc->xx + mc->yy + mc->zz)/3;
	mr->xx = mc->xx - t;
	mr->yy = mc->yy - t;
	mr->xy = mc->xy;
	mr->xz = mc->xz;
	mr->yz = mc->yz;
	/*
	 ** Finally the mass remains the same.
	 */
	mr->m = mc->m;
	}


/*
 ** This function converts a complete local moment (LOCC) to a reduced one (LOCR).
 ** This is just test code, since Locc is trace-free anyway this function is not 
 ** needed.
 */
void momReduceLocc(LOCC *mc,LOCR *mr)
{
	momFloat  t,tx,ty,tz,txx,txy,txz,tyy,tyz,tzz;

	/*
	 ** First reduce Hexadecapole.
	 */
	txx = (mc->xxxx + mc->xxyy + mc->xxzz)/7;
	txy = (mc->xxxy + mc->xyyy + mc->xyzz)/7;
	txz = (mc->xxxz + mc->xyyz + mc->xzzz)/7;
	tyy = (mc->xxyy + mc->yyyy + mc->yyzz)/7;
	tyz = (mc->xxyz + mc->yyyz + mc->yzzz)/7;
	tzz = (mc->xxzz + mc->yyzz + mc->zzzz)/7;
	t = 0.1*(txx + tyy + tzz);
	mr->xxxx = mc->xxxx - 6*(txx - t);
	mr->xyyy = mc->xyyy - 3*txy;
	mr->xxxy = mc->xxxy - 3*txy;
	mr->yyyy = mc->yyyy - 6*(tyy - t);
	mr->xxxz = mc->xxxz - 3*txz;
	mr->yyyz = mc->yyyz - 3*tyz;
	mr->xxyy = mc->xxyy - (txx + tyy - 2*t);
	mr->xxyz = mc->xxyz - tyz;
	mr->xyyz = mc->xyyz - txz;	
	/*
	 ** Now reduce the Octopole.
	 */
	tx = (mc->xxx + mc->xyy + mc->xzz)/5;
	ty = (mc->xxy + mc->yyy + mc->yzz)/5;
	tz = (mc->xxz + mc->yyz + mc->zzz)/5;
	mr->xxx = mc->xxx - 3*tx;
	mr->xyy = mc->xyy - tx;
	mr->xxy = mc->xxy - ty;
	mr->yyy = mc->yyy - 3*ty;
	mr->xxz = mc->xxz - tz;
	mr->yyz = mc->yyz - tz;
	mr->xyz = mc->xyz;
	/*
	 ** Now reduce the Quadrupole.
	 */
	t = (mc->xx + mc->yy + mc->zz)/3;
	mr->xx = mc->xx - t;
	mr->yy = mc->yy - t;
	mr->xy = mc->xy;
	mr->xz = mc->xz;
	mr->yz = mc->yz;
	mr->x = mc->x;
	mr->y = mc->y;
	mr->z = mc->z;
	/*
	 ** Finally the mass remains the same.
	 */
	mr->m = mc->m;
	}


/*
 ** This is a new fast version of QEVAL which evaluates
 ** the interaction due to the reduced moment 'm'.
 ** This version is nearly two times as fast as a naive
 ** implementation.
 **
 ** OpCount = (*,+) = (103,72) = 175 - 8 = 167
 */
void momEvalMomr(MOMR *m,momFloat dir,momFloat x,momFloat y,momFloat z,
				 momFloat *fPot,momFloat *ax,momFloat *ay,momFloat *az)
{
	const momFloat onethird = 1.0/3.0;
	momFloat xx,xy,xz,yy,yz,zz;
	momFloat xxx,xxy,xxz,xyy,yyy,yyz,xyz;
	momFloat tx,ty,tz,dir2,g2,g3,g4;

	dir = -dir;
	dir2 = dir*dir;
	g2 = 3*dir*dir2*dir2;
	g3 = -5*g2*dir2;
	g4 = -7*g3*dir2;
	/*
	 ** Calculate the funky distance terms.
	 */
	xx = 0.5*x*x;
	xy = x*y;
	xz = x*z;
    yy = 0.5*y*y;
	yz = y*z;
	zz = 0.5*z*z;
	xxx = x*(onethird*xx - zz);
	xxz = z*(xx - onethird*zz);
	yyy = y*(onethird*yy - zz);
	yyz = z*(yy - onethird*zz);
	xx -= zz;
	yy -= zz;
	xxy = y*xx;
	xyy = x*yy;
	xyz = xy*z;
	/*
	 ** Now calculate the interaction up to Hexadecapole order.
	 */
	tx = g4*(m->xxxx*xxx + m->xyyy*yyy + m->xxxy*xxy + m->xxxz*xxz + m->xxyy*xyy + m->xxyz*xyz + m->xyyz*yyz);
	ty = g4*(m->xyyy*xyy + m->xxxy*xxx + m->yyyy*yyy + m->yyyz*yyz + m->xxyy*xxy + m->xxyz*xxz + m->xyyz*xyz);
	tz = g4*(-m->xxxx*xxz - (m->xyyy + m->xxxy)*xyz - m->yyyy*yyz + m->xxxz*xxx + m->yyyz*yyy - m->xxyy*(xxz + yyz) + m->xxyz*xxy + m->xyyz*xyy);
	g4 = 0.25*(tx*x + ty*y + tz*z);
	xxx = g3*(m->xxx*xx + m->xyy*yy + m->xxy*xy + m->xxz*xz + m->xyz*yz);
	xxy = g3*(m->xyy*xy + m->xxy*xx + m->yyy*yy + m->yyz*yz + m->xyz*xz);
	xxz = g3*(-(m->xxx + m->xyy)*xz - (m->xxy + m->yyy)*yz + m->xxz*xx + m->yyz*yy + m->xyz*xy);
	g3 = onethird*(xxx*x + xxy*y + xxz*z);
	xx = g2*(m->xx*x + m->xy*y + m->xz*z);
	xy = g2*(m->yy*y + m->xy*x + m->yz*z);
	xz = g2*(-(m->xx + m->yy)*z + m->xz*x + m->yz*y);
	g2 = 0.5*(xx*x + xy*y + xz*z);
	dir *= m->m;
	dir2 *= -(dir + 5*g2 + 7*g3 + 9*g4);
	*fPot += dir + g2 + g3 + g4;
	*ax += xx + xxx + tx + x*dir2;
	*ay += xy + xxy + ty + y*dir2;
	*az += xz + xxz + tz + z*dir2;
	}


void momClearLocc(LOCC *l)
{
	l->m = 0;
	l->x = 0;
	l->y = 0;
	l->z = 0;
	l->xx = 0;
	l->xy = 0;
	l->yy = 0;
	l->xz = 0;
	l->yz = 0;
	l->zz = 0;
	l->xxx = 0;
	l->xxy = 0;
	l->xyy = 0;
	l->yyy = 0;
	l->xxz = 0;
	l->xyz = 0;
	l->yyz = 0;
	l->xzz = 0;
	l->yzz = 0;
	l->zzz = 0;
	l->xxxx = 0;
	l->xxxy = 0;
	l->xxyy = 0;
	l->xyyy = 0;
	l->yyyy = 0;
	l->xxxz = 0;
	l->xxyz = 0;
	l->xyyz = 0;
	l->yyyz = 0;
	l->xxzz = 0;
	l->xyzz = 0;
	l->yyzz = 0;
	l->xzzz = 0;
	l->yzzz = 0;
	l->zzzz = 0;
	}


void momClearLocr(LOCR *l)
{
	l->m = 0;
	l->x = 0;
	l->y = 0;
	l->z = 0;
	l->xx = 0;
	l->xy = 0;
	l->yy = 0;
	l->xz = 0;
	l->yz = 0;
	l->xxx = 0;
	l->xxy = 0;
	l->xyy = 0;
	l->yyy = 0;
	l->xxz = 0;
	l->xyz = 0;
	l->yyz = 0;
	l->xxxx = 0;
	l->xxxy = 0;
	l->xxyy = 0;
	l->xyyy = 0;
	l->yyyy = 0;
	l->xxxz = 0;
	l->xxyz = 0;
	l->xyyz = 0;
	l->yyyz = 0;
	}


void momMomr2Momc(MOMR *ma,MOMC *mc)
{
	mc->m = ma->m;
	mc->xx = ma->xx;
	mc->yy = ma->yy;
	mc->xy = ma->xy;
	mc->xz = ma->xz;
	mc->yz = ma->yz;
	mc->xxx = ma->xxx;
	mc->xyy = ma->xyy;
	mc->xxy = ma->xxy;
	mc->yyy = ma->yyy;
	mc->xxz = ma->xxz;
	mc->yyz = ma->yyz;
	mc->xyz = ma->xyz;
	mc->xxxx = ma->xxxx;
	mc->xyyy = ma->xyyy;
	mc->xxxy = ma->xxxy;
	mc->yyyy = ma->yyyy;
	mc->xxxz = ma->xxxz;
	mc->yyyz = ma->yyyz;
	mc->xxyy = ma->xxyy;
	mc->xxyz = ma->xxyz;
	mc->xyyz = ma->xyyz;
	mc->zz = -(ma->xx + ma->yy);
	mc->xzz = -(ma->xxx + ma->xyy);
	mc->yzz = -(ma->xxy + ma->yyy);
	mc->zzz = -(ma->xxz + ma->yyz);
	mc->xxzz = -(ma->xxxx + ma->xxyy);
	mc->xyzz = -(ma->xxxy + ma->xyyy);
	mc->xzzz = -(ma->xxxz + ma->xyyz);
	mc->yyzz = -(ma->xxyy + ma->yyyy);
	mc->yzzz = -(ma->xxyz + ma->yyyz);
	mc->zzzz = -(mc->xxzz + mc->yyzz);
	}


void momLoccAddMomr(LOCC *l,MOMC *m,momFloat g0,momFloat x,momFloat y,momFloat z)
{
	const momFloat onethird = 1.0/3.0;
	momFloat dir2,g1,g2,g3,g4,g5,g6,g7,g8;
	momFloat R,R2,Rx,Ry,Rz,Rxx,Rxy,Ryy,Rxz,Ryz,Rzz;
	momFloat Rxxx,Rxxy,Rxyy,Ryyy,Rxxz,Rxyz,Ryyz,Rxzz,Ryzz,Rzzz;
	momFloat Rxxxx,Rxxxy,Rxxyy,Rxyyy,Ryyyy,Rxxxz,Rxxyz,Rxyyz,Ryyyz,Rxxzz,Rxyzz,Ryyzz,Rxzzz,Ryzzz,Rzzzz;
	momFloat xx,xy,yy,xz,yz,zz,xxx,xxy,xyy,yyy,xxz,xyz,yyz,xzz,yzz,zzz;
	momFloat Q2x,Q2y,Q2z,Q2,Q3xx,Q3xy,Q3yy,Q3xz,Q3yz,Q3zz,Q3x,Q3y,Q3z,Q3;
	momFloat Q4xxx,Q4xxy,Q4xyy,Q4yyy,Q4xxz,Q4xyz,Q4yyz,Q4xzz,Q4yzz,Q4zzz;
	momFloat Q4xx,Q4xy,Q4yy,Q4xz,Q4yz,Q4zz,Q4x,Q4y,Q4z,Q4;
	momFloat tx,ty;

	g0 = -g0;
	dir2 = g0*g0;
    g1 = -g0*dir2;
	g2 = -3*g1*dir2;
	g3 = -5*g2*dir2;
	g4 = -7*g3*dir2;
	g5 = -9*g4*dir2;
	g6 = -11*g5*dir2;
	g7 = -13*g6*dir2;
	g8 = -15*g7*dir2;
	xx = x*x;
	xy = x*y;
	yy = y*y;
	xz = x*z;
	yz = y*z;
	zz = z*z;
	xxx = x*xx;
	xxy = y*xx;
	xyy = x*yy;
	yyy = y*yy;
	xxz = x*xz;
	xyz = x*yz;
	yyz = y*yz;
	xzz = x*zz;
	yzz = y*zz;
	zzz = z*zz;
	Q2x = x*m->xx + y*m->xy + z*m->xz;
	Q2y = x*m->xy + y*m->yy + z*m->yz;
	Q2z = x*m->xz + y*m->yz + z*m->zz;
	Q2 = 0.5*(x*Q2x + y*Q2y + z*Q2z);
	Q3xx = x*m->xxx + y*m->xxy + z*m->xxz;
	Q3xy = x*m->xxy + y*m->xyy + z*m->xyz;
	Q3yy = x*m->xyy + y*m->yyy + z*m->yyz;
	Q3xz = x*m->xxz + y*m->xyz + z*m->xzz;
	Q3yz = x*m->xyz + y*m->yyz + z*m->yzz;
	Q3zz = x*m->xzz + y*m->yzz + z*m->zzz;
	Q3x = 0.5*(x*Q3xx + y*Q3xy + z*Q3xz);
	Q3y = 0.5*(x*Q3xy + y*Q3yy + z*Q3yz);
	Q3z = 0.5*(x*Q3xz + y*Q3yz + z*Q3zz);
	Q3 = onethird*(x*Q3x + y*Q3y + z*Q3z);
	Q4xxx = x*m->xxxx + y*m->xxxy + z*m->xxxz;
	Q4xxy = x*m->xxxy + y*m->xxyy + z*m->xxyz;
	Q4xyy = x*m->xxyy + y*m->xyyy + z*m->xyyz;
	Q4yyy = x*m->xyyy + y*m->yyyy + z*m->yyyz;
	Q4xxz = x*m->xxxz + y*m->xxyz + z*m->xxzz;
	Q4xyz = x*m->xxyz + y*m->xyyz + z*m->xyzz;
	Q4yyz = x*m->xyyz + y*m->yyyz + z*m->yyzz;
	Q4xzz = x*m->xxzz + y*m->xyzz + z*m->xzzz;
	Q4yzz = x*m->xyzz + y*m->yyzz + z*m->yzzz;
	Q4zzz = x*m->xzzz + y*m->yzzz + z*m->zzzz;
	Q4xx = 0.5*(x*Q4xxx + y*Q4xxy + z*Q4xxz);
	Q4xy = 0.5*(x*Q4xxy + y*Q4xyy + z*Q4xyz);
	Q4yy = 0.5*(x*Q4xyy + y*Q4yyy + z*Q4yyz);
	Q4xz = 0.5*(x*Q4xxz + y*Q4xyz + z*Q4xzz);
	Q4yz = 0.5*(x*Q4xyz + y*Q4yyz + z*Q4yzz);
	Q4zz = 0.5*(x*Q4xzz + y*Q4yzz + z*Q4zzz);
	Q4x = onethird*(x*Q4xx + y*Q4xy + z*Q4xz);
	Q4y = onethird*(x*Q4xy + y*Q4yy + z*Q4yz);
	Q4z = onethird*(x*Q4xz + y*Q4yz + z*Q4zz);
	Q4 = 0.25*(x*Q4x + y*Q4y + z*Q4z);
	R = g0*m->m + g2*Q2 + g3*Q3 + g4*Q4;
	l->m += R;
	Rx = g2*Q2x + g3*Q3x + g4*Q4x;
	Ry = g2*Q2y + g3*Q3y + g4*Q4y;
	Rz = g2*Q2z + g3*Q3z + g4*Q4z;
	R = g1*m->m + g3*Q2 + g4*Q3 + g5*Q4;
	l->x += x*R + Rx;
	l->y += y*R + Ry;
	l->z += z*R + Rz;
	Rxx = g2*m->xx + g3*Q3xx + g4*Q4xx + R;
	Rxy = g2*m->xy + g3*Q3xy + g4*Q4xy;
	Ryy = g2*m->yy + g3*Q3yy + g4*Q4yy + R;
	Rxz = g2*m->xz + g3*Q3xz + g4*Q4xz;
	Ryz = g2*m->yz + g3*Q3yz + g4*Q4yz;
	Rzz = g2*m->zz + g3*Q3zz + g4*Q4zz + R;
	Rx = g3*Q2x + g4*Q3x + g5*Q4x;
	Ry = g3*Q2y + g4*Q3y + g5*Q4y;
	Rz = g3*Q2z + g4*Q3z + g5*Q4z;
	R = g2*m->m + g4*Q2 + g5*Q3 + g6*Q4;
	R2 = R;									/* need to save this one! */
	l->xx += xx*R + 2*x*Rx + Rxx;
	l->xy += xy*R + y*Rx + x*Ry + Rxy;
	l->yy += yy*R + 2*y*Ry + Ryy;
	l->xz += xz*R + z*Rx + x*Rz + Rxz;
	l->yz += yz*R + z*Ry + y*Rz + Ryz;
	l->zz += zz*R + 2*z*Rz + Rzz;
	Rxxx = g3*m->xxx + g4*Q4xxx + 3*Rx;
	Rxxy = g3*m->xxy + g4*Q4xxy + Ry;
	Rxyy = g3*m->xyy + g4*Q4xyy + Rx;
	Ryyy = g3*m->yyy + g4*Q4yyy + 3*Ry;
	Rxxz = g3*m->xxz + g4*Q4xxz + Rz;
	Rxyz = g3*m->xyz + g4*Q4xyz;
	Ryyz = g3*m->yyz + g4*Q4yyz + Rz;
	Rxzz = g3*m->xzz + g4*Q4xzz + Rx;
	Ryzz = g3*m->yzz + g4*Q4yzz + Ry;
	Rzzz = g3*m->zzz + g4*Q4zzz + 3*Rz;
	Rxx = g3*m->xx + g4*Q3xx + g5*Q4xx + R;
	Rxy = g3*m->xy + g4*Q3xy + g5*Q4xy;
	Ryy = g3*m->yy + g4*Q3yy + g5*Q4yy + R;
	Rxz = g3*m->xz + g4*Q3xz + g5*Q4xz;
	Ryz = g3*m->yz + g4*Q3yz + g5*Q4yz;
	Rzz = g3*m->zz + g4*Q3zz + g5*Q4zz + R;
	Rx = g4*Q2x + g5*Q3x + g6*Q4x;
	Ry = g4*Q2y + g5*Q3y + g6*Q4y;
	Rz = g4*Q2z + g5*Q3z + g6*Q4z;
	R = g3*m->m + g5*Q2 + g6*Q3 + g7*Q4;
	l->xxx += xxx*R + 3*xx*Rx + 3*x*Rxx + Rxxx;
	l->xxy += xxy*R + 2*xy*Rx + xx*Ry + y*Rxx + 2*x*Rxy + Rxxy;
	l->xyy += xyy*R + yy*Rx + 2*xy*Ry + 2*y*Rxy + x*Ryy + Rxyy;
	l->yyy += yyy*R + 3*yy*Ry + 3*y*Ryy + Ryyy;
	l->xxz += xxz*R + 2*xz*Rx + xx*Rz + z*Rxx + 2*x*Rxz + Rxxz;
	l->xyz += xyz*R + yz*Rx + xz*Ry + xy*Rz + z*Rxy + y*Rxz + x*Ryz + Rxyz;
	l->yyz += yyz*R + 2*yz*Ry + yy*Rz + z*Ryy + 2*y*Ryz + Ryyz;
	l->xzz += xzz*R + zz*Rx + 2*xz*Rz + 2*z*Rxz + x*Rzz + Rxzz;
	l->yzz += yzz*R + zz*Ry + 2*yz*Rz + 2*z*Ryz + y*Rzz + Ryzz;
	l->zzz += zzz*R + 3*zz*Rz + 3*z*Rzz + Rzzz;
	Rxxxx = g4*m->xxxx + 6*Rxx - 3*R2;
	Rxxxy = g4*m->xxxy + 3*Rxy;
	Rxxyy = g4*m->xxyy + Rxx + Ryy - R2;
	Rxyyy = g4*m->xyyy + 3*Rxy;
	Ryyyy = g4*m->yyyy + 6*Ryy - 3*R2;
	Rxxxz = g4*m->xxxz + 3*Rxz;
	Rxxyz = g4*m->xxyz + Ryz;
	Rxyyz = g4*m->xyyz + Rxz;
	Ryyyz = g4*m->yyyz + 3*Ryz;
	Rxxzz = g4*m->xxzz + Rxx + Rzz - R2;
	Rxyzz = g4*m->xyzz + Rxy;
	Ryyzz = g4*m->yyzz + Ryy + Rzz - R2;
	Rxzzz = g4*m->xzzz + 3*Rxz;
	Ryzzz = g4*m->yzzz + 3*Ryz;
	Rzzzz = g4*m->zzzz + 6*Rzz - 3*R2;
	Rxxx = g4*m->xxx + g5*Q4xxx + 3*Rx;
	Rxxy = g4*m->xxy + g5*Q4xxy + Ry;
	Rxyy = g4*m->xyy + g5*Q4xyy + Rx;
	Ryyy = g4*m->yyy + g5*Q4yyy + 3*Ry;
	Rxxz = g4*m->xxz + g5*Q4xxz + Rz;
	Rxyz = g4*m->xyz + g5*Q4xyz;
	Ryyz = g4*m->yyz + g5*Q4yyz + Rz;
	Rxzz = g4*m->xzz + g5*Q4xzz + Rx;
	Ryzz = g4*m->yzz + g5*Q4yzz + Ry;
	Rzzz = g4*m->zzz + g5*Q4zzz + 3*Rz;
	Rxx = g4*m->xx + g5*Q3xx + g6*Q4xx + R;
	Rxy = g4*m->xy + g5*Q3xy + g6*Q4xy;
	Ryy = g4*m->yy + g5*Q3yy + g6*Q4yy + R;
	Rxz = g4*m->xz + g5*Q3xz + g6*Q4xz;
	Ryz = g4*m->yz + g5*Q3yz + g6*Q4yz;
	Rzz = g4*m->zz + g5*Q3zz + g6*Q4zz + R;
	Rx = g5*Q2x + g6*Q3x + g7*Q4x;
	Ry = g5*Q2y + g6*Q3y + g7*Q4y;
	Rz = g5*Q2z + g6*Q3z + g7*Q4z;
	R = g4*m->m + g6*Q2 + g7*Q3 + g8*Q4;
	tx = x*R;
	ty = y*R;
	l->xxxx += xxx*(tx + 4*Rx) + 6*xx*Rxx + 4*x*Rxxx + Rxxxx;
	l->xxxy += xxy*(tx + 3*Rx) + xxx*Ry + 3*xy*Rxx + 3*xx*Rxy + y*Rxxx + 3*x*Rxxy + Rxxxy;
	l->xxyy += xyy*(tx + 2*Rx) + 2*xxy*Ry + yy*Rxx + 4*xy*Rxy + xx*Ryy + 2*y*Rxxy + 2*x*Rxyy + Rxxyy;
	l->xyyy += yyy*(tx + Rx) + 3*xyy*Ry + 3*yy*Rxy + 3*xy*Ryy + 3*y*Rxyy + x*Ryyy + Rxyyy;
	l->yyyy += yyy*(ty + 4*Ry) + 6*yy*Ryy + 4*y*Ryyy + Ryyyy;
	l->xxxz += xxz*(tx + 3*Rx) + xxx*Rz + 3*xz*Rxx + 3*xx*Rxz + z*Rxxx + 3*x*Rxxz + Rxxxz;
	l->xxyz += xyz*(tx + 2*Rx) + xxz*Ry + xxy*Rz + yz*Rxx + 2*xz*Rxy + 2*xy*Rxz + xx*Ryz + z*Rxxy + y*Rxxz + 2*x*Rxyz + Rxxyz;
	l->xyyz += yyz*(tx + Rx) + 2*xyz*Ry + xyy*Rz + 2*yz*Rxy + yy*Rxz + xz*Ryy + 2*xy*Ryz + z*Rxyy + 2*y*Rxyz + x*Ryyz + Rxyyz;
	l->yyyz += yyz*(ty + 3*Ry) + yyy*Rz + 3*yz*Ryy + 3*yy*Ryz + z*Ryyy + 3*y*Ryyz + Ryyyz;
	l->xxzz += xzz*(tx + 2*Rx) + 2*xxz*Rz + zz*Rxx + 4*xz*Rxz + xx*Rzz + 2*z*Rxxz + 2*x*Rxzz + Rxxzz;
	l->xyzz += yzz*(tx + Rx) + xzz*Ry + 2*xyz*Rz + zz*Rxy + 2*yz*Rxz + 2*xz*Ryz + xy*Rzz + 2*z*Rxyz + y*Rxzz + x*Ryzz + Rxyzz;
	l->yyzz += yzz*(ty + 2*Ry) + 2*yyz*Rz + zz*Ryy + 4*yz*Ryz + yy*Rzz + 2*z*Ryyz + 2*y*Ryzz + Ryyzz;
	l->xzzz += zzz*(tx + Rx) + 3*xzz*Rz + 3*zz*Rxz + 3*xz*Rzz + 3*z*Rxzz + x*Rzzz + Rxzzz;
	l->yzzz += zzz*(ty + Ry) + 3*yzz*Rz + 3*zz*Ryz + 3*yz*Rzz + 3*z*Ryzz + y*Rzzz + Ryzzz;
	l->zzzz += zzz*(z*R + 4*Rz) + 6*zz*Rzz + 4*z*Rzzz + Rzzzz;
	}


void momLocrAddMomr(LOCR *l,MOMR *m,momFloat g0,momFloat x,momFloat y,momFloat z)
{
	const momFloat onethird = 1.0/3.0;
	momFloat dir2,g1,g2,g3,g4,g5,g6,g7,g8;
	momFloat R,Rx,Ry,Rz,Rxx,Rxy,Ryy,Rxz,Ryz;
	momFloat Rxxx,Rxxy,Rxyy,Ryyy,Rxxz,Rxyz,Ryyz;
	momFloat Rxxxx,Rxxxy,Rxxyy,Rxyyy,Ryyyy,Rxxxz,Rxxyz,Rxyyz,Ryyyz;
	momFloat xx,xy,yy,xz,yz,xxx,xxy,xyy,yyy,xxz,xyz,yyz;
	momFloat Q2x,Q2y,Q2z,Q2,Q3xx,Q3xy,Q3yy,Q3xz,Q3yz,Q3x,Q3y,Q3z,Q3;
	momFloat Q4xxx,Q4xxy,Q4xyy,Q4yyy,Q4xxz,Q4xyz,Q4yyz;
	momFloat Q4xx,Q4xy,Q4yy,Q4xz,Q4yz,Q4x,Q4y,Q4z,Q4;
	momFloat tx,ty;

	g0 = -g0;
	dir2 = g0*g0;
    g1 = -g0*dir2;
	g2 = -3*g1*dir2;
	g3 = -5*g2*dir2;
	g4 = -7*g3*dir2;
	g5 = -9*g4*dir2;
	g6 = -11*g5*dir2;
	g7 = -13*g6*dir2;
	g8 = -15*g7*dir2;
	xx = x*x;
	xy = x*y;
	yy = y*y;
	xz = x*z;
	yz = y*z;
	xxx = x*xx;
	xxy = y*xx;
	xyy = x*yy;
	yyy = y*yy;
	xxz = x*xz;
	xyz = x*yz;
	yyz = y*yz;
	Q2x = x*m->xx + y*m->xy + z*m->xz;
	Q2y = x*m->xy + y*m->yy + z*m->yz;
	Q2z = x*m->xz + y*m->yz - z*(m->xx + m->yy);
	Q2 = 0.5*(x*Q2x + y*Q2y + z*Q2z);
	Q3xx = x*m->xxx + y*m->xxy + z*m->xxz;
	Q3xy = x*m->xxy + y*m->xyy + z*m->xyz;
	Q3yy = x*m->xyy + y*m->yyy + z*m->yyz;
	Q3xz = x*m->xxz + y*m->xyz - z*(m->xxx + m->xyy);
	Q3yz = x*m->xyz + y*m->yyz - z*(m->xxy + m->yyy);
	Q3x = 0.5*(x*Q3xx + y*Q3xy + z*Q3xz);
	Q3y = 0.5*(x*Q3xy + y*Q3yy + z*Q3yz);
	Q3z = 0.5*(x*Q3xz + y*Q3yz - z*(Q3xx + Q3yy));
	Q3 = onethird*(x*Q3x + y*Q3y + z*Q3z);
	Q4xxx = x*m->xxxx + y*m->xxxy + z*m->xxxz;
	Q4xxy = x*m->xxxy + y*m->xxyy + z*m->xxyz;
	Q4xyy = x*m->xxyy + y*m->xyyy + z*m->xyyz;
	Q4yyy = x*m->xyyy + y*m->yyyy + z*m->yyyz;
	Q4xxz = x*m->xxxz + y*m->xxyz - z*(m->xxxx + m->xxyy);
	Q4xyz = x*m->xxyz + y*m->xyyz - z*(m->xxxy + m->xyyy);
	Q4yyz = x*m->xyyz + y*m->yyyz - z*(m->xxyy + m->yyyy);
	Q4xx = 0.5*(x*Q4xxx + y*Q4xxy + z*Q4xxz);
	Q4xy = 0.5*(x*Q4xxy + y*Q4xyy + z*Q4xyz);
	Q4yy = 0.5*(x*Q4xyy + y*Q4yyy + z*Q4yyz);
	Q4xz = 0.5*(x*Q4xxz + y*Q4xyz - z*(Q4xxx + Q4xyy));
	Q4yz = 0.5*(x*Q4xyz + y*Q4yyz - z*(Q4xxy + Q4yyy));
	Q4x = onethird*(x*Q4xx + y*Q4xy + z*Q4xz);
	Q4y = onethird*(x*Q4xy + y*Q4yy + z*Q4yz);
	Q4z = onethird*(x*Q4xz + y*Q4yz - z*(Q4xx + Q4yy));
	Q4 = 0.25*(x*Q4x + y*Q4y + z*Q4z);
	R = g0*m->m + g2*Q2 + g3*Q3 + g4*Q4;
	l->m += R;
	Rx = g2*Q2x + g3*Q3x + g4*Q4x;
	Ry = g2*Q2y + g3*Q3y + g4*Q4y;
	Rz = g2*Q2z + g3*Q3z + g4*Q4z;
	R = g1*m->m + g3*Q2 + g4*Q3 + g5*Q4;
	l->x += x*R + Rx;
	l->y += y*R + Ry;
	l->z += z*R + Rz;
	Rxx = g2*m->xx + g3*Q3xx + g4*Q4xx + R;
	Rxy = g2*m->xy + g3*Q3xy + g4*Q4xy;
	Ryy = g2*m->yy + g3*Q3yy + g4*Q4yy + R;
	Rxz = g2*m->xz + g3*Q3xz + g4*Q4xz;
	Ryz = g2*m->yz + g3*Q3yz + g4*Q4yz;
	Rx = g3*Q2x + g4*Q3x + g5*Q4x;
	Ry = g3*Q2y + g4*Q3y + g5*Q4y;
	Rz = g3*Q2z + g4*Q3z + g5*Q4z;
	R = g2*m->m + g4*Q2 + g5*Q3 + g6*Q4;
	l->xx += xx*R + 2*x*Rx + Rxx;
	l->xy += xy*R + y*Rx + x*Ry + Rxy;
	l->yy += yy*R + 2*y*Ry + Ryy;
	l->xz += xz*R + z*Rx + x*Rz + Rxz;
	l->yz += yz*R + z*Ry + y*Rz + Ryz;
	Rxxx = g3*m->xxx + g4*Q4xxx + 3*Rx;
	Rxxy = g3*m->xxy + g4*Q4xxy + Ry;
	Rxyy = g3*m->xyy + g4*Q4xyy + Rx;
	Ryyy = g3*m->yyy + g4*Q4yyy + 3*Ry;
	Rxxz = g3*m->xxz + g4*Q4xxz + Rz;
	Rxyz = g3*m->xyz + g4*Q4xyz;
	Ryyz = g3*m->yyz + g4*Q4yyz + Rz;
	Rxx = g3*m->xx + g4*Q3xx + g5*Q4xx + R;
	Rxy = g3*m->xy + g4*Q3xy + g5*Q4xy;
	Ryy = g3*m->yy + g4*Q3yy + g5*Q4yy + R;
	Rxz = g3*m->xz + g4*Q3xz + g5*Q4xz;
	Ryz = g3*m->yz + g4*Q3yz + g5*Q4yz;
	Rx = g4*Q2x + g5*Q3x + g6*Q4x;
	Ry = g4*Q2y + g5*Q3y + g6*Q4y;
	Rz = g4*Q2z + g5*Q3z + g6*Q4z;
	Rxxxx = g4*m->xxxx + 6*Rxx - 3*R;
	Rxxxy = g4*m->xxxy + 3*Rxy;
	Rxxyy = g4*m->xxyy + Rxx + Ryy - R;
	Rxyyy = g4*m->xyyy + 3*Rxy;
	Ryyyy = g4*m->yyyy + 6*Ryy - 3*R;
	Rxxxz = g4*m->xxxz + 3*Rxz;
	Rxxyz = g4*m->xxyz + Ryz;
	Rxyyz = g4*m->xyyz + Rxz;
	Ryyyz = g4*m->yyyz + 3*Ryz;
	R = g3*m->m + g5*Q2 + g6*Q3 + g7*Q4;
	l->xxx += xxx*R + 3*xx*Rx + 3*x*Rxx + Rxxx;
	l->xxy += xxy*R + 2*xy*Rx + xx*Ry + y*Rxx + 2*x*Rxy + Rxxy;
	l->xyy += xyy*R + yy*Rx + 2*xy*Ry + 2*y*Rxy + x*Ryy + Rxyy;
	l->yyy += yyy*R + 3*yy*Ry + 3*y*Ryy + Ryyy;
	l->xxz += xxz*R + 2*xz*Rx + xx*Rz + z*Rxx + 2*x*Rxz + Rxxz;
	l->xyz += xyz*R + yz*Rx + xz*Ry + xy*Rz + z*Rxy + y*Rxz + x*Ryz + Rxyz;
	l->yyz += yyz*R + 2*yz*Ry + yy*Rz + z*Ryy + 2*y*Ryz + Ryyz;
	Rxxx = g4*m->xxx + g5*Q4xxx + 3*Rx;
	Rxxy = g4*m->xxy + g5*Q4xxy + Ry;
	Rxyy = g4*m->xyy + g5*Q4xyy + Rx;
	Ryyy = g4*m->yyy + g5*Q4yyy + 3*Ry;
	Rxxz = g4*m->xxz + g5*Q4xxz + Rz;
	Rxyz = g4*m->xyz + g5*Q4xyz;
	Ryyz = g4*m->yyz + g5*Q4yyz + Rz;
	Rxx = g4*m->xx + g5*Q3xx + g6*Q4xx + R;
	Rxy = g4*m->xy + g5*Q3xy + g6*Q4xy;
	Ryy = g4*m->yy + g5*Q3yy + g6*Q4yy + R;
	Rxz = g4*m->xz + g5*Q3xz + g6*Q4xz;
	Ryz = g4*m->yz + g5*Q3yz + g6*Q4yz;
	Rx = g5*Q2x + g6*Q3x + g7*Q4x;
	Ry = g5*Q2y + g6*Q3y + g7*Q4y;
	Rz = g5*Q2z + g6*Q3z + g7*Q4z;
	R = g4*m->m + g6*Q2 + g7*Q3 + g8*Q4;
	tx = x*R;
	ty = y*R;
	l->xxxx += xxx*(tx + 4*Rx) + 6*xx*Rxx + 4*x*Rxxx + Rxxxx;
	l->xxxy += xxy*(tx + 3*Rx) + xxx*Ry + 3*xy*Rxx + 3*xx*Rxy + y*Rxxx + 3*x*Rxxy + Rxxxy;
	l->xxyy += xyy*(tx + 2*Rx) + 2*xxy*Ry + yy*Rxx + 4*xy*Rxy + xx*Ryy + 2*y*Rxxy + 2*x*Rxyy + Rxxyy;
	l->xyyy += yyy*(tx + Rx) + 3*xyy*Ry + 3*yy*Rxy + 3*xy*Ryy + 3*y*Rxyy + x*Ryyy + Rxyyy;
	l->yyyy += yyy*(ty + 4*Ry) + 6*yy*Ryy + 4*y*Ryyy + Ryyyy;
	l->xxxz += xxz*(tx + 3*Rx) + xxx*Rz + 3*xz*Rxx + 3*xx*Rxz + z*Rxxx + 3*x*Rxxz + Rxxxz;
	l->xxyz += xyz*(tx + 2*Rx) + xxz*Ry + xxy*Rz + yz*Rxx + 2*xz*Rxy + 2*xy*Rxz + xx*Ryz + z*Rxxy + y*Rxxz + 2*x*Rxyz + Rxxyz;
	l->xyyz += yyz*(tx + Rx) + 2*xyz*Ry + xyy*Rz + 2*yz*Rxy + yy*Rxz + xz*Ryy + 2*xy*Ryz + z*Rxyy + 2*y*Rxyz + x*Ryyz + Rxyyz;
	l->yyyz += yyz*(ty + 3*Ry) + yyy*Rz + 3*yz*Ryy + 3*yy*Ryz + z*Ryyy + 3*y*Ryyz + Ryyyz;
	}


void momEvalLocc(LOCC *l,momFloat x,momFloat y,momFloat z,
				 momFloat *fPot,momFloat *ax,momFloat *ay,momFloat *az)
{
	const momFloat sixth = 1.0/6.0;
	const momFloat third = 1.0/3.0;
	momFloat xx,xy,yy,xz,yz,zz,xxx,xxy,xyy,yyy,xxz,xyz,yyz,xzz,yzz,zzz;
	momFloat L1,L2,L3,L4,L2x,L2y,L2z,L3x,L3y,L3z,L4x,L4y,L4z;

	L1 = x*l->x + y*l->y + z*l->z;
	L2x = x*l->xx + y*l->xy + z*l->xz;
	L2y = x*l->xy + y*l->yy + z*l->yz;
	L2z = x*l->xz + y*l->yz + z*l->zz;
	L2 = 0.5*(x*L2x + y*L2y + z*L2z);
	xx = x*x;
	xy = x*y;
	yy = y*y;
	xz = x*z;
	yz = y*z;
	zz = z*z;
	L3x = 0.5*(l->xxx*xx + l->xyy*yy + l->xzz*zz + 2*(l->xxy*xy + l->xxz*xz + l->xyz*yz));
	L3y = 0.5*(l->xxy*xx + l->yyy*yy + l->yzz*zz + 2*(l->xyy*xy + l->xyz*xz + l->yyz*yz));
	L3z = 0.5*(l->xxz*xx + l->yyz*yy + l->zzz*zz + 2*(l->xyz*xy + l->xzz*xz + l->yzz*yz));
	L3 = third*(x*L3x + y*L3y + z*L3z);
	xxx = x*xx;
	xxy = y*xx;
	xyy = x*yy;
	yyy = y*yy;
	yyz = z*yy;
	xxz = x*xz;
	xyz = y*xz;
	xzz = x*zz;
	yzz = y*zz;
	zzz = z*zz;
	L4x = sixth*(l->xxxx*xxx + l->xyyy*yyy + l->xzzz*zzz + 3*(l->xxxy*xxy + l->xxyy*xyy + l->xxxz*xxz + l->xyyz*yyz + l->xyzz*yzz + 2*l->xxyz*xyz));
	L4y = sixth*(l->xxxy*xxx + l->yyyy*yyy + l->yzzz*zzz + 3*(l->xxyy*xxy + l->xyyy*xyy + l->xxyz*xxz + l->yyyz*yyz + l->yyzz*yzz + 2*l->xyyz*xyz));
	L4z = sixth*(l->xxxz*xxx + l->yyyz*yyy + l->zzzz*zzz + 3*(l->xxyz*xxy + l->xyyz*xyy + l->xxzz*xxz + l->yyzz*yyz + l->yzzz*yzz + 2*l->xyzz*xyz));
	L4 = 0.25*(x*L4x + y*L4y + z*L4z);
	*ax += l->x + L2x + L3x + L4x;
	*ay += l->y + L2y + L3y + L4y;
	*az += l->z + L2z + L3z + L4z;
	*fPot += l->m + L1 + L2 + L3 + L4;
	}


/*
 ** This function subtracts the complete local moment ma from the complete 
 ** local moment mc
 */
void momSubLocc(LOCC *mc,LOCC *ma)
{
	mc->m -= ma->m;
	mc->x -= ma->x;
	mc->y -= ma->y;
	mc->z -= ma->z;	
	mc->xx -= ma->xx;
	mc->yy -= ma->yy;
	mc->xy -= ma->xy;
	mc->xz -= ma->xz;
	mc->yz -= ma->yz;
	mc->xxx -= ma->xxx;
	mc->xyy -= ma->xyy;
	mc->xxy -= ma->xxy;
	mc->yyy -= ma->yyy;
	mc->xxz -= ma->xxz;
	mc->yyz -= ma->yyz;
	mc->xyz -= ma->xyz;
	mc->xxxx -= ma->xxxx;
	mc->xyyy -= ma->xyyy;
	mc->xxxy -= ma->xxxy;
	mc->yyyy -= ma->yyyy;
	mc->xxxz -= ma->xxxz;
	mc->yyyz -= ma->yyyz;
	mc->xxyy -= ma->xxyy;
	mc->xxyz -= ma->xxyz;
	mc->xyyz -= ma->xyyz;
	mc->zz -= ma->zz;
	mc->xzz -= ma->xzz;
	mc->yzz -= ma->yzz;
	mc->zzz -= ma->zzz;
	mc->xxzz -= ma->xxzz;
	mc->xyzz -= ma->xyzz;
	mc->xzzz -= ma->xzzz;
	mc->yyzz -= ma->yyzz;
	mc->yzzz -= ma->yzzz;
	mc->zzzz -= ma->zzzz;
	}


void momPrintMomc(MOMC *m) {
	printf("MOMC:%20.15g\n",(double)m->m);
	printf("  xx:%20.15g   yy:%20.15g   zz:%20.15g\n",(double)m->xx,(double)m->yy,(double)m->zz);
    printf("  xy:%20.15g   yz:%20.15g   xz:%20.15g\n",(double)m->xy,(double)m->yz,(double)m->xz);
	printf(" xxx:%20.15g  xyy:%20.15g  xzz:%20.15g\n",(double)m->xxx,(double)m->xyy,(double)m->xzz);
    printf(" xxy:%20.15g  yyy:%20.15g  yzz:%20.15g\n",(double)m->xxy,(double)m->yyy,(double)m->yzz);
	printf(" xxz:%20.15g  yyz:%20.15g  zzz:%20.15g\n",(double)m->xxz,(double)m->yyz,(double)m->zzz);
	printf(" xyz:%20.15g\n",(double)m->xyz);
	printf("xxxx:%20.15g xxxy:%20.15g xxxz:%20.15g\n",(double)m->xxxx,(double)m->xxxy,(double)m->xxxz);
	printf("xyyy:%20.15g yyyy:%20.15g yyyz:%20.15g\n",(double)m->xyyy,(double)m->yyyy,(double)m->yyyz);
	printf("xzzz:%20.15g yzzz:%20.15g zzzz:%20.15g\n",(double)m->xzzz,(double)m->yzzz,(double)m->zzzz);
	printf("xxyy:%20.15g xxyz:%20.15g xyyz:%20.15g\n",(double)m->xxyy,(double)m->xxyz,(double)m->xyyz);
	printf("yyzz:%20.15g xxzz:%20.15g xyzz:%20.15g\n",(double)m->yyzz,(double)m->xxzz,(double)m->xyzz);
	}


void momPrintMomr(MOMR *m) {
	printf("MOMR:%20.15g\n",(double)m->m);
	printf("  xx:%20.15g   yy:%20.15g\n",(double)m->xx,(double)m->yy);
    printf("  xy:%20.15g   yz:%20.15g   xz:%20.15g\n",(double)m->xy,(double)m->yz,(double)m->xz);
	printf(" xxx:%20.15g  xyy:%20.15g\n",(double)m->xxx,(double)m->xyy);
    printf(" xxy:%20.15g  yyy:%20.15g\n",(double)m->xxy,(double)m->yyy);
	printf(" xxz:%20.15g  yyz:%20.15g\n",(double)m->xxz,(double)m->yyz);
	printf(" xyz:%20.15g\n",(double)m->xyz);
	printf("xxxx:%20.15g xxxy:%20.15g xxxz:%20.15g\n",(double)m->xxxx,(double)m->xxxy,(double)m->xxxz);
	printf("xyyy:%20.15g yyyy:%20.15g yyyz:%20.15g\n",(double)m->xyyy,(double)m->yyyy,(double)m->yyyz);
	printf("xxyy:%20.15g xxyz:%20.15g xyyz:%20.15g\n",(double)m->xxyy,(double)m->xxyz,(double)m->xyyz);
	}


void momPrintLocc(LOCC *m) {
	printf("LOCC:%20.15g\n",(double)m->m);
	printf("   x:%20.15g    y:%20.15g    z:%20.15g\n",(double)m->x,(double)m->y,(double)m->z);
	printf("  xx:%20.15g   yy:%20.15g   zz:%20.15g\n",(double)m->xx,(double)m->yy,(double)m->zz);
    printf("  xy:%20.15g   yz:%20.15g   xz:%20.15g\n",(double)m->xy,(double)m->yz,(double)m->xz);
	printf(" xxx:%20.15g  xyy:%20.15g  xzz:%20.15g\n",(double)m->xxx,(double)m->xyy,(double)m->xzz);
    printf(" xxy:%20.15g  yyy:%20.15g  yzz:%20.15g\n",(double)m->xxy,(double)m->yyy,(double)m->yzz);
	printf(" xxz:%20.15g  yyz:%20.15g  zzz:%20.15g\n",(double)m->xxz,(double)m->yyz,(double)m->zzz);
	printf(" xyz:%20.15g\n",(double)m->xyz);
	printf("xxxx:%20.15g xxxy:%20.15g xxxz:%20.15g\n",(double)m->xxxx,(double)m->xxxy,(double)m->xxxz);
	printf("xyyy:%20.15g yyyy:%20.15g yyyz:%20.15g\n",(double)m->xyyy,(double)m->yyyy,(double)m->yyyz);
	printf("xzzz:%20.15g yzzz:%20.15g zzzz:%20.15g\n",(double)m->xzzz,(double)m->yzzz,(double)m->zzzz);
	printf("xxyy:%20.15g xxyz:%20.15g xyyz:%20.15g\n",(double)m->xxyy,(double)m->xxyz,(double)m->xyyz);
	printf("yyzz:%20.15g xxzz:%20.15g xyzz:%20.15g\n",(double)m->yyzz,(double)m->xxzz,(double)m->xyzz);
	}


void momPrintLocr(LOCR *m) {
	printf("LOCR:%20.15g\n",(double)m->m);
	printf("   x:%20.15g    y:%20.15g    z:%20.15g\n",(double)m->x,(double)m->y,(double)m->z);
	printf("  xx:%20.15g   yy:%20.15g\n",(double)m->xx,(double)m->yy);
    printf("  xy:%20.15g   yz:%20.15g   xz:%20.15g\n",(double)m->xy,(double)m->yz,(double)m->xz);
	printf(" xxx:%20.15g  xyy:%20.15g\n",(double)m->xxx,(double)m->xyy);
    printf(" xxy:%20.15g  yyy:%20.15g\n",(double)m->xxy,(double)m->yyy);
	printf(" xxz:%20.15g  yyz:%20.15g\n",(double)m->xxz,(double)m->yyz);
	printf(" xyz:%20.15g\n",(double)m->xyz);
	printf("xxxx:%20.15g xxxy:%20.15g xxxz:%20.15g\n",(double)m->xxxx,(double)m->xxxy,(double)m->xxxz);
	printf("xyyy:%20.15g yyyy:%20.15g yyyz:%20.15g\n",(double)m->xyyy,(double)m->yyyy,(double)m->yyyz);
	printf("xxyy:%20.15g xxyz:%20.15g xyyz:%20.15g\n",(double)m->xxyy,(double)m->xxyz,(double)m->xyyz);
	}


