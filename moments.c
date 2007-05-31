#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <math.h>
#include <assert.h>
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
 ** This function shifts a reduced local expansion (LOCR) to a new center of expansion.
 ** <x,y,z> := d := rexp(new) - rexp(old).
 **
 ** Op Count (*,+) = (97,100) = 197
 */
void momShiftLocr(LOCR *l,momFloat x,momFloat y,momFloat z)
{
  const momFloat onethird = 1.0/3.0;
  momFloat L,Lx,Ly,Lz,Lxx,Lxy,Lxz,Lyy,Lyz,Lxxx,Lxxy,Lxxz,Lxyy,Lxyz,Lyyy,Lyyz,hx,hy,hz;

  L = l->x*x + l->y*y + l->z*z; 
  l->m += L;

  hx = 0.5*x;
  hy = 0.5*x;
  hz = 0.5*x;

  Lx = l->xx*x + l->xy*y + l->xz*z;
  Ly = l->xy*x + l->yy*y + l->yz*z;
  Lz = l->xz*x + l->yz*y - (l->xx + l->yy)*z;
  L = Lx*hx + Ly*hy + Lz*hz;
  l->x += Lx;
  l->y += Ly;
  l->z += Lz;
  l->m += L;

  Lxx = l->xxx*x + l->xxy*y + l->xxz*z;
  Lxy = l->xxy*x + l->xyy*y + l->xyz*z;
  Lxz = l->xxz*x + l->xyz*y - (l->xxx + l->xyy)*z;
  Lyy = l->xyy*x + l->yyy*y + l->yyz*z;
  Lyz = l->xyz*x + l->yyz*y - (l->xxy + l->yyy)*z;
  Lx = Lxx*hx + Lxy*hy + Lxz*hz;
  Ly = Lxy*hx + Lyy*hy + Lyz*hz;
  Lz = Lxz*hx + Lyz*hy - (Lxx + Lyy)*hz;

  l->xx += Lxx;
  l->xy += Lxy;
  l->xz += Lxz;
  l->yy += Lyy;
  l->yz += Lyz;
  l->x += Lx;
  l->y += Ly;
  l->z += Lz;

  Lxxx = l->xxxx*x + l->xxxy*y + l->xxxz*z;
  Lxxy = l->xxxy*x + l->xxyy*y + l->xxyz*z;
  Lxxz = l->xxxz*x + l->xxyz*y - (l->xxxx + l->xxyy)*z;
  Lxyy = l->xxyy*x + l->xyyy*y + l->xyyz*z;
  Lxyz = l->xxyz*x + l->xyyz*y - (l->xxxy + l->xyyy)*z;
  Lyyy = l->xyyy*x + l->yyyy*y + l->yyyz*z;
  Lyyz = l->xyyz*x + l->yyyz*y - (l->xxyy + l->yyyy)*z;
  Lxx = Lxxx*hx + Lxxy*hy + Lxxz*hz;
  Lxy = Lxxy*hx + Lxyy*hy + Lxyz*hz;
  Lxz = Lxxz*hx + Lxyz*hy - (Lxxx + Lxyy)*hz;
  Lyy = Lxyy*hx + Lyyy*hy + Lyyz*hz;
  Lyz = Lxyz*hx + Lyyz*hy - (Lxxy + Lyyy)*hz;
  x *= onethird;
  y *= onethird;
  z *= onethird;
  L = Lx*x + Ly*y + Lz*z;  /* note here we use the old Lx,Ly,Lz */
  l->m += L;
  Lx = Lxx*x + Lxy*y + Lxz*z;
  Ly = Lxy*x + Lyy*y + Lyz*z;
  Lz = Lxz*x + Lyz*y - (Lxx + Lyy)*z;
  L = Lx*hx + Ly*hy + Lz*hz;

  l->xxx += Lxxx;
  l->xxy += Lxxy;
  l->xxz += Lxxz;
  l->xyy += Lxyy;
  l->xyz += Lxyz;
  l->yyy += Lyyy;
  l->yyz += Lyyz;
  l->xx += Lxx;
  l->xy += Lxy;
  l->xz += Lxz;
  l->yy += Lyy;
  l->yz += Lyz;
  l->x += Lx;
  l->y += Ly;
  l->z += Lz;
  l->m += 0.5*L;
}


/*
 ** This function shifts a reduced local expansion of 2ND ORDER (quadrupole) to a new center of expansion.
 ** <x,y,z> := d := rexp(new) - rexp(old).
 */
void momShiftLocr2(LOCR *l,momFloat x,momFloat y,momFloat z)
{
  momFloat L,Lx,Ly,Lz;

  L = l->x*x + l->y*y + l->z*z; 
  l->m += L;
  Lx = l->xx*x + l->xy*y + l->xz*z;
  Ly = l->xy*x + l->yy*y + l->yz*z;
  Lz = l->xz*x + l->yz*y - (l->xx + l->yy)*z;
  L = Lx*x + Ly*y + Lz*z;
  l->x += Lx;
  l->y += Ly;
  l->z += Lz;
  l->m += 0.5*L;
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
 ** March 23, 2007: This function now uses unit vectors 
 ** which reduces the required precision in the exponent
 ** since the highest power of r is now 5 (g4 ~ r^(-5)).
 **
 ** OpCount = (*,+) = (105,72) = 177 - 8 = 169
 **
 ** CAREFUL: this function no longer accumulates on fPot,ax,ay,az!
 */
void momEvalMomr(MOMR *m,momFloat dir,momFloat x,momFloat y,momFloat z,
		 momFloat *fPot,momFloat *ax,momFloat *ay,momFloat *az,momFloat *magai)
{
	const momFloat onethird = 1.0/3.0;
	momFloat xx,xy,xz,yy,yz,zz;
	momFloat xxx,xxy,xxz,xyy,yyy,yyz,xyz;
	momFloat tx,ty,tz,g0,g2,g3,g4;

	g0 = dir;
	g2 = 3*dir*dir*dir;
	g3 = 5*g2*dir;
	g4 = 7*g3*dir;
	/*
	 ** Calculate the funky distance terms.
	 */
	x *= dir;
	y *= dir;
	z *= dir;
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
	g0 *= m->m;
	*fPot = -(g0 + g2 + g3 + g4);
	g0 += 5*g2 + 7*g3 + 9*g4;
	*ax = dir*(xx + xxx + tx - x*g0);
	*ay = dir*(xy + xxy + ty - y*g0);
	*az = dir*(xz + xxz + tz - z*g0);
	*magai = g0*dir;
	}


/*
** The generalized version of the above.
**
** CAREFUL: this function no longer accumulates on fPot,ax,ay,az!
*/
void momGenEvalMomr(MOMR *m,momFloat dir,momFloat g0,momFloat t1,momFloat t2,
		    momFloat t3r,momFloat t4r,momFloat x,momFloat y,momFloat z,
		    momFloat *fPot,momFloat *ax,momFloat *ay,momFloat *az,momFloat *magai)
{
	const momFloat onethird = 1.0/3.0;
	momFloat xx,xy,xz,yy,yz,zz;
	momFloat xxx,xxy,xxz,xyy,yyy,yyz,xyz;
	momFloat tx,ty,tz,g2,g3,g4;

	g2 = t2*t1*g0;
	g3 = t3r*dir*g2;
	g4 = t4r*dir*g3;
	/*
	 ** Calculate the funky distance terms.
	 */
	x *= dir;
	y *= dir;
	z *= dir;
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
	g0 *= m->m;
	*fPot = g0 + g2 + g3 + g4;
	g0 += t3r*g2 + t4r*g3;      /* for exact deriv + t5r*g4 */
	*ax = dir*(xx + xxx + tx + x*g0);
	*ay = dir*(xy + xxy + ty + y*g0);
	*az = dir*(xz + xxz + tz + z*g0);
	*magai = -g0*dir;
	}


/*
 ** This is a new version of EvalMomr which is designed to be 
 ** less sensitive to the size of the exponent in single precision.
 ** Highest direct power of r is 2 (as it should be for a force)
 ** when using scaled moments (providing a scale factor of s, some size
 ** of the cell).
 **
 ** OpCount = (*,+) = (108,73) = 181 - 8 = 173
 */
void momEvalMomrs(MOMR *m,momFloat s,momFloat dir,
		  momFloat x,momFloat y,momFloat z,
		  momFloat *fPot,momFloat *ax,momFloat *ay,momFloat *az)
{
	const momFloat onethird = 1.0/3.0;
	momFloat xx,xy,xz,yy,yz,zz;
	momFloat xxx,xxy,xxz,xyy,yyy,yyz,xyz;
	momFloat tx,ty,tz,sir,g0,g2,g3,g4;

	sir = s*dir;
	g2 = -3*dir*sir*sir;
	g3 = -5*g2*sir;
	g4 = -7*g3*sir;
	/*
	 ** Calculate the funky distance terms.
	 */
	x *= dir;
	y *= dir;
	z *= dir;
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
	g0 = m->m*dir;
	*fPot += g2 + g3 + g4 - g0;
	g0 -= 5*g2 + 7*g3 + 9*g4;
	*ax += dir*(xx + xxx + tx + x*g0);
	*ay += dir*(xy + xxy + ty + y*g0);
	*az += dir*(xz + xxz + tz + z*g0);
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

void momClearMomr(MOMR *l)
{
	l->m = 0;
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

/*
** Op Count (*,+) = (240,188) = 428
** Op Count per Multipole = 214
*/
void momSymLocrAddMomr(LOCR *l1,LOCR *l2,MOMR *q1,MOMR *q2,momFloat dir,momFloat x,momFloat y,momFloat z) {
    const momFloat onethird = 1.0/3.0;
    momFloat xx,xy,xz,yy,yz,zz;
    momFloat m1,m2,Ax,Ay,Az,A,Bx,By,Bz,B,C1,C2,R,T,t,f,t3xx,t3yy,t4xx,t4yy;
    momFloat g0,g1,t2,g2,g2t,t3,g3,g3t,t4,g4;

    g0 = -dir;
    g1 = -g0*dir;
    t2 = -3*dir;
    g2 = t2*g1;
    g2t = g2*dir;
    t3 = -5*dir;
    g3 = t3*g2;
    g3t = g3*dir;
    t4 = -7*dir;

    x *= dir;
    y *= dir;
    z *= dir;
    xx = 0.5*x*x;
    xy = x*y;
    yy = 0.5*y*y;
    xz = x*z;
    yz = y*z;
    zz = 0.5*z*z;
    xx -= zz;
    yy -= zz;

    Ax = x*q2->xx + y*q2->xy + z*q2->xz;
    Ay = x*q2->xy + y*q2->yy + z*q2->yz;
    Az = x*q2->xz + y*q2->yz - z*(q2->xx + q2->yy);
    A = 0.5*g2*(x*Ax + y*Ay + z*Az);
    Ax *= g2t;
    Ay *= g2t;
    Az *= g2t;

    Bx = xx*q2->xxx + xy*q2->xxy + xz*q2->xxz + yy*q2->xyy + yz*q2->xyz;
    By = xx*q2->xxy + xy*q2->xyy + xz*q2->xyz + yy*q2->yyy + yz*q2->yyz;
    Bz = xx*q2->xxz + xy*q2->xyz - xz*(q2->xxx + q2->xyy) + yy*q2->yyz - yz*(q2->xxy + q2->yyy);
    B = onethird*g3*(x*Bx + y*By + z*Bz);
    Bx *= g3t;
    By *= g3t;
    Bz *= g3t;

    l1->m += g0*q2->m + A + B;

    A *= t3;

    m2 = g1*q2->m;
    R = m2 + A + t4*B;
    l1->x += Ax + Bx + x*R;
    l1->y += Ay + By + y*R;
    l1->z += Az + Bz + z*R;

    Ax *= t3;
    Ay *= t3;
    Az *= t3;

    T = dir*(m2 + A);
    m2 *= t2;
    R = m2 + t4*A;
    l1->xx += T + (R*x + 2*Ax)*x;
    l1->yy += T + (R*y + 2*Ay)*y;
    l1->xy += R*xy + (Ax*y + Ay*x);
    l1->xz += R*xz + (Ax*z + Az*x);
    l1->yz += R*yz + (Ay*z + Az*y);

    Ax = x*q1->xx + y*q1->xy + z*q1->xz;
    Ay = x*q1->xy + y*q1->yy + z*q1->yz;
    Az = x*q1->xz + y*q1->yz - z*(q1->xx + q1->yy);
    A = 0.5*g2*(x*Ax + y*Ay + z*Az);


    Bx = xx*q1->xxx + xy*q1->xxy + xz*q1->xxz + yy*q1->xyy + yz*q1->xyz;
    By = xx*q1->xxy + xy*q1->xyy + xz*q1->xyz + yy*q1->yyy + yz*q1->yyz;
    Bz = xx*q1->xxz + xy*q1->xyz - xz*(q1->xxx + q1->xyy) + yy*q1->yyz - yz*(q1->xxy + q1->yyy);
    B = onethird*g3*(x*Bx + y*By + z*Bz);

    l2->m += g0*q1->m + A - B;

    A *= t3;

    m1 = g1*q1->m;
    R = m1 + A - t4*B;
    l2->x += -Ax + Bx - x*R;
    l2->y += -Ay + By - y*R;
    l2->z += -Az + Bz - z*R;

    Ax *= t3;
    Ay *= t3;
    Az *= t3;

    T = dir*(m1 + A);
    m1 *= t2;
    R = m1 + t4*A;
    l2->xx += T + (R*x + 2*Ax)*x;
    l2->yy += T + (R*y + 2*Ay)*y;
    l2->xy += R*xy + (Ax*y + Ay*x);
    l2->xz += R*xz + (Ax*z + Az*x);
    l2->yz += R*yz + (Ay*z + Az*y);

    /*
    ** No more A's and no more B's used here!
    */
    t = xx*xx - xz*xz;     C1 =  t*q2->xxxx; C2 =  t*q1->xxxx;
    t = yy*yy - yz*yz;     C1 += t*q2->yyyy; C2 += t*q1->yyyy; 
    t = 2*xx*xz;           C1 += t*q2->xxxz; C2 += t*q1->xxxz; 
    t = 2*yy*yz;           C1 += t*q2->yyyz; C2 += t*q1->yyyz; 
    t = 6*yy*xx - 4*zz*zz; C1 += t*q2->xxyy; C2 += t*q1->xxyy;
    xx = x*x;
    yy = y*y;
    zz *= 2;
    t = (3*xx - zz)*yz;    C1 += t*q2->xxyz; C2 += t*q1->xxyz;
    t = (3*yy - zz)*xz;    C1 += t*q2->xyyz; C2 += t*q1->xyyz;
    zz *= 3;
    t = xy*(xx - zz);      C1 += t*q2->xxxy; C2 += t*q1->xxxy;
    t = xy*(yy - zz);      C1 += t*q2->xyyy; C2 += t*q1->xyyy; 

    g4 = (7.0/6.0)*g3t;
    l1->m += g4*C1;
    l2->m += g4*C2;

    m1 *= dir;
    m2 *= dir;
    t3xx = -5*xx;
    t3yy = -5*yy;
    t = (3 + t3xx)*x; l1->xxx += t*m2; l2->xxx -= t*m1;
    t = (3 + t3yy)*y; l1->yyy += t*m2; l2->yyy -= t*m1;
    f = (1 + t3xx);
    t = f*y;          l1->xxy += t*m2; l2->xxy -= t*m1;
    t = f*z;          l1->xxz += t*m2; l2->xxz -= t*m1;
    f = (1 + t3yy);
    t = f*x;          l1->xyy += t*m2; l2->xyy -= t*m1;
    t = f*z;          l1->yyz += t*m2; l2->yyz -= t*m1;
    t = -5*xy*z;      l1->xyz += t*m2; l2->xyz -= t*m1;

    m1 *= dir;
    m2 *= dir;
    t4xx = t4*xx;
    t4yy = t4*yy;
    t = (3 + (6 + t4xx)*t3xx);        l1->xxxx += t*m2; l2->xxxx += t*m1;
    t = (1 + t3xx + t3yy*(1 + t4xx)); l1->xxyy += t*m2; l2->xxyy += t*m1;
    t = (3 + (6 + t4yy)*t3yy);        l1->yyyy += t*m2; l2->yyyy += t*m1;
    m1 *= -5;
    m2 *= -5;
    f = (3 + t4xx);
    t = f*xy;                         l1->xxxy += t*m2; l2->xxxy += t*m1;
    t = f*xz;                         l1->xxxz += t*m2; l2->xxxz += t*m1;
    f = (3 + t4yy);
    t = f*xy;                         l1->xyyy += t*m2; l2->xyyy += t*m1;
    t = f*yz;                         l1->yyyz += t*m2; l2->yyyz += t*m1;
    t = (1 + t4xx)*yz;                l1->xxyz += t*m2; l2->xxyz += t*m1;
    t = (1 + t4yy)*xz;                l1->xyyz += t*m2; l2->xyyz += t*m1;
    }

/*
** Op Count (*,+) = (154,105) = 259 
** Op Count (*,+) = (151,92) = 243
** via: gcc -c -S -03 t.c
**      grep mul t.s|wc
**      grep 'add.*xmm' t.s|wc
*/
double momLocrAddMomr(LOCR *l,MOMR *q,momFloat dir,momFloat x,momFloat y,momFloat z) {
    const momFloat onethird = 1.0/3.0;
    momFloat xx,xy,xz,yy,yz,zz;
    momFloat Ax,Ay,Az,A,Bx,By,Bz,B,C,R,T;
    momFloat g0,g1,g2,g2t,g3,g3t,m,t2,t3,t4,t3xx,t3yy,t4xx,t4yy,f;

    g0 = -dir;
    g1 = -g0*dir;
    t2 = -3*dir;
    g2 = t2*g1;
    g2t = g2*dir;
    t3 = -5*dir;
    g3 = t3*g2;
    g3t = g3*dir;
    t4 = -7*dir;

    x *= dir;
    y *= dir;
    z *= dir;
    xx = 0.5*x*x;
    xy = x*y;
    yy = 0.5*y*y;
    xz = x*z;
    yz = y*z;
    zz = 0.5*z*z;
    xx -= zz;
    yy -= zz;

    Ax = x*q->xx + y*q->xy + z*q->xz;
    Ay = x*q->xy + y*q->yy + z*q->yz;
    Az = x*q->xz + y*q->yz - z*(q->xx + q->yy);
    A = 0.5*g2*(x*Ax + y*Ay + z*Az);
    Ax *= g2t;
    Ay *= g2t;
    Az *= g2t;

    Bx = xx*q->xxx + xy*q->xxy + xz*q->xxz + yy*q->xyy + yz*q->xyz;
    By = xx*q->xxy + xy*q->xyy + xz*q->xyz + yy*q->yyy + yz*q->yyz;
    Bz = xx*q->xxz + xy*q->xyz - xz*(q->xxx + q->xyy) + yy*q->yyz - yz*(q->xxy + q->yyy);
    B = onethird*g3*(x*Bx + y*By + z*Bz);
    Bx *= g3t;
    By *= g3t;
    Bz *= g3t;

    l->m += g0*q->m + A + B;

    A *= t3;

    m = g1*q->m;
    R = m + A + t4*B;
    l->x += Ax + Bx + x*R;
    l->y += Ay + By + y*R;
    l->z += Az + Bz + z*R;
    /*
    ** No more B's used here!
    */
    Ax *= t3;
    Ay *= t3;
    Az *= t3;

    T = dir*(m + A);
    m *= t2;
    R = m + t4*A;
    l->xx += T + (R*x + 2*Ax)*x;
    l->yy += T + (R*y + 2*Ay)*y;
    l->xy += R*xy + (Ax*y + Ay*x);
    l->xz += R*xz + (Ax*z + Az*x);
    l->yz += R*yz + (Ay*z + Az*y);
    /*
    ** No more A's and no more B's used here!
    */
    C = (xx*xx - xz*xz)*q->xxxx + (yy*yy - yz*yz)*q->yyyy + 
      2*(xx*xz*q->xxxz + yy*yz*q->yyyz) + (6*yy*xx - 4*zz*zz)*q->xxyy;
    xx = x*x;
    yy = y*y;
    zz *= 2;
    C += (3*xx - zz)*yz*q->xxyz  + (3*yy - zz)*xz*q->xyyz;
    zz *= 3;
    C += xy*((xx - zz)*q->xxxy + (yy - zz)*q->xyyy); 

    l->m -= (7.0/6.0)*g3t*C;

    m *= dir;
    t3xx = -5*xx;
    t3yy = -5*yy;
    l->xxx += (3 + t3xx)*x*m;
    l->yyy += (3 + t3yy)*y*m;
    f = m*(1 + t3xx);
    l->xxy += f*y;
    l->xxz += f*z;
    f = m*(1 + t3yy);
    l->xyy += f*x;
    l->yyz += f*z;
    l->xyz += -5*xy*z*m;

    m *= dir;
    t4xx = -7*xx;
    t4yy = -7*yy;
    l->xxxx += (3 + (6 + t4xx)*t3xx)*m;
    l->xxyy += (1 + t3xx + t3yy*(1 + t4xx))*m;
    l->yyyy += (3 + (6 + t4yy)*t3yy)*m;
    m *= -5;
    f = m*(3 + t4xx);
    l->xxxy += f*xy;
    l->xxxz += f*xz;
    f = m*(3 + t4yy);
    l->xyyy += f*xy;
    l->yyyz += f*yz;
    l->xxyz += (1 + t4xx)*yz*m;
    l->xyyz += (1 + t4yy)*xz*m;

    return 243.0;
    }

/*
** Op Count (*,+) = (154,105) = 259 
**  
**    Below an explanation to the factors that are passed to this function and their definitions
**    for the Newtonian Green's function. They can be set to perform the Ewald summation terms
**    as well.
**
**    g0 = -dir;
**    g1 = -g0*dir;
**    g2 = -3*g1*dir;
**    g3 = -5*g2*dir;
**    g4 = -7*g3*dir;
**
**    t1 = g1/g0 = -dir;
**    t2 = g2/g1 = -3*dir;
**    t3 = g3/g2 = -5*dir;
**    t4 = g4/g3 = -7*dir;
**
**    t3r = t3*r = -5;
**    t4r = t4*r = -7;
**
**    This is how we would calculate these factors for part of the Ewald
**    summation. This gets a little complicated, but it is equivalent to 
**    the existing code.
** 
**	L = 1.0;
**	alpha = 2.0/L;
**	alpha2 = alpha*alpha;
**	ka = 2.0*alpha/sqrt(M_PI);
**	a = exp(-r2*alpha2);
**	a *= ka*dir;
**	alpha /= dir;
**	alpha2 /= dir;
**	if (bInHole) g0 = -erf(alpha);
**	else g0 = erfc(alpha);
**	g0 *= dir;
**	g1 = g0*dir + a;
**	alphan = 2*alpha2;
**	g2 = 3*g1*dir + alphan*a;
**	alphan *= 2*alpha2/dir;
**	g3r = 5*g2 + alphan*a;
**	g3 = g3r*dir;
**	alphan *= 2*alpha2;
**	g4r = 7*g3 + alphan*a;
**	g4 = g4r*dir;
**
**	t1 = g1/g0;
**	t2 = g2/g1;
**	t3r = g3r/g2;
**	t4r = g4r/g3;
*/
void momGenLocrAddMomr(LOCR *l,MOMR *q,momFloat dir,
		       momFloat g0,momFloat t1,momFloat t2,momFloat t3r,momFloat t4r,
		       momFloat x,momFloat y,momFloat z) {
    const momFloat onethird = 1.0/3.0;
    momFloat xx,xy,xz,yy,yz,zz;
    momFloat Ax,Ay,Az,A,Bx,By,Bz,B,C,R,T;
    momFloat g1,g2,g3,g2t,g3t,m,t3,t4,t3xx,t3yy,t4xx,t4yy,f;

    t3 = t3r*dir;
    t4 = t4r*dir;
    g1 = t1*g0;
    g2 = t2*g1;
    g3 = t3*g2;
    g2t = g2*dir;
    g3t = g3*dir;

    x *= dir;
    y *= dir;
    z *= dir;

    zz = z*z;
    yy = 0.5*(y*y - zz);
    xx = 0.5*(x*x - zz);
    xy = x*y;
    xz = x*z;
    yz = y*z;

    Ax = x*q->xx + y*q->xy + z*q->xz;
    Ay = x*q->xy + y*q->yy + z*q->yz;
    Az = x*q->xz + y*q->yz - z*(q->xx + q->yy);
    A = 0.5*g2*(x*Ax + y*Ay + z*Az);
    Ax *= g2t;
    Ay *= g2t;
    Az *= g2t;

    Bx = xx*q->xxx + xy*q->xxy + xz*q->xxz + yy*q->xyy + yz*q->xyz;
    By = xx*q->xxy + xy*q->xyy + xz*q->xyz + yy*q->yyy + yz*q->yyz;
    Bz = xx*q->xxz + xy*q->xyz - xz*(q->xxx + q->xyy) + yy*q->yyz - yz*(q->xxy + q->yyy);
    B = onethird*g3*(x*Bx + y*By + z*Bz);
    Bx *= g3t;
    By *= g3t;
    Bz *= g3t;
    l->m += g0*q->m + A + B;

    A *= t3;

    m = g1*q->m;
    R = m + A + t4*B;
    l->x += Ax + Bx + x*R;
    l->y += Ay + By + y*R;
    l->z += Az + Bz + z*R;
    /*
    ** No more B's used here!
    */
    Ax *= t3;
    Ay *= t3;
    Az *= t3;

    T = dir*(m + A);
    m *= t2;
    R = m + t4*A;
    l->xx += T + (R*x + 2*Ax)*x;
    l->yy += T + (R*y + 2*Ay)*y;
    l->xy += R*xy + (Ax*y + Ay*x);
    l->xz += R*xz + (Ax*z + Az*x);
    l->yz += R*yz + (Ay*z + Az*y);
    /*
    ** No more A's and no more B's used here!
    */
    C = (xx*xx - xz*xz)*q->xxxx + (yy*yy - yz*yz)*q->yyyy + 
      2*(xx*xz*q->xxxz + yy*yz*q->yyyz) + (6*yy*xx - zz*zz)*q->xxyy;
    xx = x*x;
    yy = y*y;
    C += (3*xx - zz)*yz*q->xxyz  + (3*yy - zz)*xz*q->xyyz;
    zz *= 3;
    C += xy*((xx - zz)*q->xxxy + (yy - zz)*q->xyyy); 

    l->m += (1.0/6.0)*t4*g3*C;

    m *= dir;
    t3xx = t3r*xx;
    t3yy = t3r*yy;
    l->xxx += (3 + t3xx)*x*m;
    l->yyy += (3 + t3yy)*y*m;
    f = m*(1 + t3xx);
    l->xxy += f*y;
    l->xxz += f*z;
    f = m*(1 + t3yy);
    l->xyy += f*x;
    l->yyz += f*z;
    l->xyz += t3r*xy*z*m;

    m *= dir;
    t4xx = t4r*xx;
    t4yy = t4r*yy;
    l->xxxx += (3 + (6 + t4xx)*t3xx)*m;
    l->xxyy += (1 + t3xx + t3yy*(1 + t4xx))*m;
    l->yyyy += (3 + (6 + t4yy)*t3yy)*m;
    m *= t3r;
    f = m*(3 + t4xx);
    l->xxxy += f*xy;
    l->xxxz += f*xz;
    f = m*(3 + t4yy);
    l->xyyy += f*xy;
    l->yyyz += f*yz;
    l->xxyz += (1 + t4xx)*yz*m;
    l->xyyz += (1 + t4yy)*xz*m;
    }

void momEwaldLocrAddMomr(LOCR *l,MOMR *m,momFloat r2,int bInHole,momFloat x,momFloat y,momFloat z) {
    momFloat xx,xy,xz,yy,yz,zz;
    momFloat Ax,Ay,Az,A,Bx,By,Bz,B,C,R,T;
    momFloat g0,g1,g2,g3,g4,dir,dir2,a,alphan,alpha,alpha2,ka,L;

    L = 1.0;
    alpha = 2.0/L;
    alpha2 = alpha*alpha;
    ka = 2.0*alpha/sqrt(M_PI);

    dir = 1.0/sqrt(r2);
    dir2 = dir*dir;
    a = exp(-r2*alpha2);
    a *= ka*dir2;
    if (bInHole) g0 = -erf(alpha/dir);
    else g0 = erfc(alpha/dir);
    g0 *= dir;
    g1 = g0*dir2 + a;
    alphan = 2*alpha2;
    g2 = 3*g1*dir2 + alphan*a;
    alphan *= 2*alpha2;
    g3 = 5*g2*dir2 + alphan*a;
    alphan *= 2*alpha2;
    g4 = 7*g3*dir2 + alphan*a;

    printf("g0:%.20g g1:%.20g g2:%.20g g3:%.20g g4:%.20g\n",
	   g0,g1,g2,g3,g4);

    xx = 0.5*x*x;
    xy = x*y;
    yy = 0.5*y*y;
    xz = x*z;
    yz = y*z;
    zz = 0.5*z*z;
    xx -= zz;
    yy -= zz;

    Ax = x*m->xx + y*m->xy + z*m->xz;
    Ay = x*m->xy + y*m->yy + z*m->yz;
    Az = x*m->xz + y*m->yz - z*(m->xx + m->yy);
    A = x*Ax + y*Ay + z*Az;

    Bx = xx*m->xxx + xy*m->xxy + xz*m->xxz + yy*m->xyy + yz*m->xyz;
    By = xx*m->xxy + xy*m->xyy + xz*m->xyz + yy*m->yyy + yz*m->yyz;
    Bz = xx*m->xxz + xy*m->xyz - xz*(m->xxx + m->xyy) + yy*m->yyz - yz*(m->xxy + m->yyy);
    B = x*Bx + y*By + z*Bz;

    l->m += g0*m->m + 0.5*g2*A + (1/3.0)*g3*B;

    R = g1*m->m + 0.5*g3*A + (1/3.0)*g4*B;
    l->x += g2*Ax + g3*Bx + x*R;
    l->y += g2*Ay + g3*By + y*R;
    l->z += g2*Az + g3*Bz + z*R;
    /*
    ** No more B's used here!
    */
    T = g1*m->m + 0.5*g3*A;
    R = g2*m->m + 0.5*g4*A;
    l->xx += T + (R*x + 2*g3*Ax)*x;
    l->xy += R*xy + g3*(Ax*y + Ay*x);
    l->yy += T + (R*y + 2*g3*Ay)*y;
    l->xz += R*xz + g3*(Ax*z + Az*x);
    l->yz += R*yz + g3*(Ay*z + Az*y);
    /*
    ** No more A's and no more B's used here!
    */
    C = (xx*xx - xz*xz)*m->xxxx + (yy*yy - yz*yz)*m->yyyy + 2*(xx*xz*m->xxxz + yy*yz*m->yyyz) + (6*yy*xx - 4*zz*zz)*m->xxyy;
    xx = x*x;
    yy = y*y;
    zz *= 2;
    C += (3*xx - zz)*yz*m->xxyz  + (3*yy - zz)*xz*m->xyyz;
    zz *= 3;
    C += xy*((xx - zz)*m->xxxy + (yy - zz)*m->xyyy); 

    l->m += (1/6.0)*g4*C;

    l->xxx += (3*g2 + g3*x*x)*x*m->m;
    l->xxy += (g2 + g3*x*x)*y*m->m;
    l->xyy += (g2 + g3*y*y)*x*m->m;
    l->yyy += (3*g2 + g3*y*y)*y*m->m;
    l->xxz += (g2 + g3*x*x)*z*m->m;
    l->xyz += g3*x*y*z*m->m;
    l->yyz += (g2 + g3*y*y)*z*m->m;

    l->xxxx += (3*g2 + (6*g3 + g4*x*x)*x*x)*m->m;
    l->xxxy += (3*g3 + g4*x*x)*x*y*m->m;
    l->xxyy += (g2 + g3*(x*x + y*y) + g4*x*x*y*y)*m->m;
    l->xyyy += (3*g3 + g4*y*y)*x*y*m->m;
    l->yyyy += (3*g2 + (6*g3 + g4*y*y)*y*y)*m->m;
    l->xxxz += (3*g3 + g4*x*x)*x*z*m->m;
    l->xxyz += (g3 + g4*x*x)*y*z*m->m;
    l->xyyz += (g3 + g4*y*y)*x*z*m->m;
    l->yyyz += (3*g3 + g4*y*y)*y*z*m->m;
    }



void momNooptLocrAddMomr(LOCR *l,MOMR *m,momFloat dir,momFloat x,momFloat y,momFloat z) {
    momFloat xx,xy,xz,yy,yz,zz;
    momFloat Ax,Ay,Az,A,Bx,By,Bz,B,C,R,T;
    momFloat g0,g1,g2,g3,g4,dir2;

    g0 = -dir;
    dir2 = dir*dir;
    g1 = -g0*dir2;
    g2 = -3*g1*dir2;
    g3 = -5*g2*dir2;
    g4 = -7*g3*dir2;

    xx = 0.5*x*x;
    xy = x*y;
    yy = 0.5*y*y;
    xz = x*z;
    yz = y*z;
    zz = 0.5*z*z;
    xx -= zz;
    yy -= zz;

    Ax = x*m->xx + y*m->xy + z*m->xz;
    Ay = x*m->xy + y*m->yy + z*m->yz;
    Az = x*m->xz + y*m->yz - z*(m->xx + m->yy);
    A = x*Ax + y*Ay + z*Az;

    Bx = xx*m->xxx + xy*m->xxy + xz*m->xxz + yy*m->xyy + yz*m->xyz;
    By = xx*m->xxy + xy*m->xyy + xz*m->xyz + yy*m->yyy + yz*m->yyz;
    Bz = xx*m->xxz + xy*m->xyz - xz*(m->xxx + m->xyy) + yy*m->yyz - yz*(m->xxy + m->yyy);
    B = x*Bx + y*By + z*Bz;

    l->m += g0*m->m + 0.5*g2*A + (1/3.0)*g3*B;

    R = g1*m->m + 0.5*g3*A + (1/3.0)*g4*B;
    l->x += g2*Ax + g3*Bx + x*R;
    l->y += g2*Ay + g3*By + y*R;
    l->z += g2*Az + g3*Bz + z*R;
    /*
    ** No more B's used here!
    */
    T = g1*m->m + 0.5*g3*A;
    R = g2*m->m + 0.5*g4*A;
    l->xx += T + (R*x + 2*g3*Ax)*x;
    l->xy += R*xy + g3*(Ax*y + Ay*x);
    l->yy += T + (R*y + 2*g3*Ay)*y;
    l->xz += R*xz + g3*(Ax*z + Az*x);
    l->yz += R*yz + g3*(Ay*z + Az*y);
    /*
    ** No more A's and no more B's used here!
    */
    C = (xx*xx - xz*xz)*m->xxxx + (yy*yy - yz*yz)*m->yyyy + 2*(xx*xz*m->xxxz + yy*yz*m->yyyz) + (6*yy*xx - 4*zz*zz)*m->xxyy;
    xx = x*x;
    yy = y*y;
    zz *= 2;
    C += (3*xx - zz)*yz*m->xxyz  + (3*yy - zz)*xz*m->xyyz;
    zz *= 3;
    C += xy*((xx - zz)*m->xxxy + (yy - zz)*m->xyyy); 

    l->m += (1/6.0)*g4*C;

    l->xxx += (3*g2 + g3*x*x)*x*m->m;
    l->xxy += (g2 + g3*x*x)*y*m->m;
    l->xyy += (g2 + g3*y*y)*x*m->m;
    l->yyy += (3*g2 + g3*y*y)*y*m->m;
    l->xxz += (g2 + g3*x*x)*z*m->m;
    l->xyz += g3*x*y*z*m->m;
    l->yyz += (g2 + g3*y*y)*z*m->m;

    l->xxxx += (3*g2 + (6*g3 + g4*x*x)*x*x)*m->m;
    l->xxxy += (3*g3 + g4*x*x)*x*y*m->m;
    l->xxyy += (g2 + g3*(x*x + y*y) + g4*x*x*y*y)*m->m;
    l->xyyy += (3*g3 + g4*y*y)*x*y*m->m;
    l->yyyy += (3*g2 + (6*g3 + g4*y*y)*y*y)*m->m;
    l->xxxz += (3*g3 + g4*x*x)*x*z*m->m;
    l->xxyz += (g3 + g4*x*x)*y*z*m->m;
    l->xyyz += (g3 + g4*y*y)*x*z*m->m;
    l->yyyz += (3*g3 + g4*y*y)*y*z*m->m;
    }



void momLoccAddMomrAccurate(LOCC *l,MOMC *m,momFloat g0,momFloat x,momFloat y,momFloat z)
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


void momLocrAddMomrAccurate(LOCR *l,MOMR *m,momFloat g0,momFloat x,momFloat y,momFloat z)
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


void momEvalLocr(LOCR *l,momFloat x,momFloat y,momFloat z,
				 momFloat *fPot,momFloat *ax,momFloat *ay,momFloat *az)
{
	const momFloat onethird = 1.0/3.0;
	momFloat xx,xy,xz,yy,yz,zz,xxx,xxz,yyy,yyz,xxy,xyy,xyz,tx,ty,tz,g1,g2,g3,g4;

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
	tx = l->xxxx*xxx + l->xyyy*yyy + l->xxxy*xxy + l->xxxz*xxz + l->xxyy*xyy + l->xxyz*xyz + l->xyyz*yyz;
	ty = l->xyyy*xyy + l->xxxy*xxx + l->yyyy*yyy + l->yyyz*yyz + l->xxyy*xxy + l->xxyz*xxz + l->xyyz*xyz;
	tz = -l->xxxx*xxz - (l->xyyy + l->xxxy)*xyz - l->yyyy*yyz + l->xxxz*xxx + l->yyyz*yyy - l->xxyy*(xxz + yyz) + l->xxyz*xxy + l->xyyz*xyy;
	g4 = 0.25*(tx*x + ty*y + tz*z);
	xxx = l->xxx*xx + l->xyy*yy + l->xxy*xy + l->xxz*xz + l->xyz*yz;
	xxy = l->xyy*xy + l->xxy*xx + l->yyy*yy + l->yyz*yz + l->xyz*xz;
	xxz = -(l->xxx + l->xyy)*xz - (l->xxy + l->yyy)*yz + l->xxz*xx + l->yyz*yy + l->xyz*xy;
	g3 = onethird*(xxx*x + xxy*y + xxz*z);
	xx = l->xx*x + l->xy*y + l->xz*z;
	xy = l->yy*y + l->xy*x + l->yz*z;
	xz = -(l->xx + l->yy)*z + l->xz*x + l->yz*y;
	g2 = 0.5*(xx*x + xy*y + xz*z);
	g1 = x*l->x + y*l->y + z*l->z;
	*ax += l->x + xx + xxx + tx;
	*ay += l->y + xy + xxy + ty;
	*az += l->z + xz + xxz + tz;
	*fPot += l->m + g1 + g2 + g3 + g4;
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
	printf("  xx:%20.15g   xy:%20.15g   xy:%20.15g\n",(double)m->xx,(double)m->xy,(double)m->xz);
	printf("  yy:%20.15g   yz:%20.15g  xxx:%20.15g\n",(double)m->yy,(double)m->yz,(double)m->xxx);
	printf(" xxy:%20.15g  xxz:%20.15g  xyy:%20.15g\n",(double)m->xxy,(double)m->xxz,(double)m->xyy);
	printf(" xyz:%20.15g  yyy:%20.15g  yyz:%20.15g\n",(double)m->xyz,(double)m->yyy,(double)m->yyz);
	printf("xxxx:%20.15g xxxy:%20.15g xxxz:%20.15g\n",(double)m->xxxx,(double)m->xxxy,(double)m->xxxz);
	printf("xxyy:%20.15g xxyz:%20.15g xyyy:%20.15g\n",(double)m->xxyy,(double)m->xxyz,(double)m->xyyy);
	printf("xyyz:%20.15g yyyy:%20.15g yyyz:%20.15g\n",(double)m->xyyz,(double)m->yyyy,(double)m->yyyz);
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
	printf("  xx:%20.15g   xy:%20.15g   xz:%20.15g\n",(double)m->xx,(double)m->xy,(double)m->xz);
	printf("  yy:%20.15g   yz:%20.15g  xxx:%20.15g\n",(double)m->yy,(double)m->yz,(double)m->xxx);
	printf(" xxy:%20.15g  xxz:%20.15g  xyy:%20.15g\n",(double)m->xxy,(double)m->xxz,(double)m->xyy);
	printf(" xyz:%20.15g  yyy:%20.15g  yyz:%20.15g\n",(double)m->xyz,(double)m->yyy,(double)m->yyz);
	printf("xxxx:%20.15g xxxy:%20.15g xxxz:%20.15g\n",(double)m->xxxx,(double)m->xxxy,(double)m->xxxz);
	printf("xxyy:%20.15g xxyz:%20.15g xyyy:%20.15g\n",(double)m->xxyy,(double)m->xxyz,(double)m->xyyy);
	printf("xyyz:%20.15g yyyy:%20.15g yyyz:%20.15g\n",(double)m->xyyz,(double)m->yyyy,(double)m->yyyz);
	}


