#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#define MPICH_SKIP_MPICXX
#include "simd.h"
#include "pkd.h"

extern "C"
void pkdGravEvalPC(PINFOIN *pPart, int nBlocks, int nInLast, ILC_BLK *blk,  PINFOOUT *pOut ) {

    const fvec onethird = 1.0f/3.0f;
    fvec u,g0,g1,g2,g3,g4;
    fvec x,y,z;
    fvec adotai;
    fvec tx,ty,tz;
    fvec xx,xy,xz,yy,yz,zz;
    fvec xxx,xxz,yyy,yyz,xxy,xyy,xyz;
    fvec dir;
    fvec dimaga, tax, tay, taz;

    int j, n, nLeft, nIntr;

    fvec fx = pPart->r[0];
    fvec fy = pPart->r[1];
    fvec fz = pPart->r[2];
    float *a = pPart->a;

    fvec ax,ay,az,fPot,dirsum,normsum;

    float a2 = a[0]*a[0] + a[1]*a[1] + a[2]*a[2];
    dimaga = a2 > 0.0f ? 1.0f / sqrtf(a2) : 0.0f;

    /* Pad the last value if necessary */
    for( j = nInLast; j&SIMD_MASK; j++) {
	blk[nBlocks].dx.f[j] = blk[nBlocks].dy.f[j] = blk[nBlocks].dz.f[j] = 1e18f;
	blk[nBlocks].m.f[j] = 0.0f;
	blk[nBlocks].u.f[j] = 0.0f;
	}

    ax = 0.0;
    ay = 0.0;
    az = 0.0;
    fPot= 0.0;
    dirsum = 0.0;
    normsum = 0.0;
    nIntr = nBlocks * ILP_PART_PER_BLK + nInLast;
    for( nLeft=nBlocks; nLeft >= 0; --nLeft,++blk ) {
	int n = ((nLeft ? ILP_PART_PER_BLK : nInLast) + SIMD_MASK) >> SIMD_BITS;
	for (j=0; j<n; ++j) {
	    fvec dx = blk->dx.p[j] + fx;
	    fvec dy = blk->dy.p[j] + fy;
	    fvec dz = blk->dz.p[j] + fz;
	    fvec d2 = dx*dx + dy*dy + dz*dz;
	    dir = rsqrt(d2);
	    u = blk->u.p[j]*dir;
	    g0 = dir;
	    g1 = g0*u;
	    g2 = 3*g1*u;
	    g3 = 5*g2*u;
	    g4 = 7*g3*u;
	    /*
	    ** Calculate the funky distance terms.
	    */
	    x = dx*dir;
	    y = dy*dir;
	    z = dz*dir;
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
	    tx = g4*(blk->xxxx.p[j]*xxx + blk->xyyy.p[j]*yyy + fvec(blk->xxxy.p[j]*xxy) + blk->xxxz.p[j]*xxz + blk->xxyy.p[j]*xyy + blk->xxyz.p[j]*xyz + blk->xyyz.p[j]*yyz);
	    ty = g4*(blk->xyyy.p[j]*xyy + blk->xxxy.p[j]*xxx + blk->yyyy.p[j]*yyy + blk->yyyz.p[j]*yyz + blk->xxyy.p[j]*xxy + blk->xxyz.p[j]*xxz + blk->xyyz.p[j]*xyz);
	    tz = g4*(-fvec(blk->xxxx.p[j])*xxz - (fvec(blk->xyyy.p[j]) + blk->xxxy.p[j])*xyz - blk->yyyy.p[j]*yyz + blk->xxxz.p[j]*xxx + blk->yyyz.p[j]*yyy - blk->xxyy.p[j]*(xxz + yyz) + blk->xxyz.p[j]*xxy + blk->xyyz.p[j]*xyy);
	    g4 = 0.25*(tx*x + ty*y + tz*z);
	    xxx = g3*(blk->xxx.p[j]*xx + blk->xyy.p[j]*yy + blk->xxy.p[j]*xy + blk->xxz.p[j]*xz + blk->xyz.p[j]*yz);
	    xxy = g3*(blk->xyy.p[j]*xy + blk->xxy.p[j]*xx + blk->yyy.p[j]*yy + blk->yyz.p[j]*yz + blk->xyz.p[j]*xz);
	    xxz = g3*(-(fvec(blk->xxx.p[j]) + blk->xyy.p[j])*xz - (fvec(blk->xxy.p[j]) + blk->yyy.p[j])*yz + blk->xxz.p[j]*xx + blk->yyz.p[j]*yy + blk->xyz.p[j]*xy);
	    g3 = onethird*(xxx*x + xxy*y + xxz*z);
	    xx = g2*(blk->xx.p[j]*x + blk->xy.p[j]*y + blk->xz.p[j]*z);
	    xy = g2*(blk->yy.p[j]*y + blk->xy.p[j]*x + blk->yz.p[j]*z);
	    xz = g2*(-(fvec(blk->xx.p[j]) + blk->yy.p[j])*z + blk->xz.p[j]*x + blk->yz.p[j]*y);
	    g2 = 0.5*(xx*x + xy*y + xz*z);
	    g0 *= fvec(blk->m.p[j]);
	    fPot -= g0 + g2 + g3 + g4;
	    g0 += 5*g2 + 7*g3 + 9*g4;
#ifdef USE_DIAPOLE
	    yy = g1*blk->x.p[j];
	    yz = g1*blk->y.p[j];
	    zz = g1*blk->z.p[j];
	    g1 = (yy*x + yz*y + zz*z);
	    fPot -= g1;
	    g0 += 3*g1;
#else
	    yy = 0.0f;
	    yz = 0.0f;
	    zz = 0.0f;
#endif
	    tax = dir*(yy + xx + xxx + tx - x*g0);
	    tay = dir*(yz + xy + xxy + ty - y*g0);
	    taz = dir*(zz + xz + xxz + tz - z*g0);
	    /*
	    ** Calculations for determining the timestep.
	    */
#if 0
	    adotai = pPart->a[0]*tax + pPart->a[1]*tay + pPart->a[2]*taz;
	    if (adotai > 0.0f) {
		adotai *= dimaga;
		dirsum += dir*adotai*adotai;
		normsum += adotai*adotai;
		}
#endif
	    ax += tax;
	    ay += tay;
	    az += taz;
	    }
	}

    pOut->a[0] += hadd(ax);
    pOut->a[1] += hadd(ay);
    pOut->a[2] += hadd(az);
    pOut->fPot += hadd(fPot);
    pOut->dirsum += hadd(dirsum);
    pOut->normsum += hadd(normsum);
    }
