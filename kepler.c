#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef SYMBA /* also needs PLANET */ 
/* The routines in kepler.c solve kepler's equation for various ecc. and dm */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "kepler.h"

int drift_dan(double mu, double *r, double *v, double dt){
  /*This subroutine does the Danby and decides which vbles to use
  ** 
  **  Input:                
  *     mu          ==>  solar mass 
  **    r[3],v[3]   ==>  initial position and velocity in Democratic coordinate
  **    dt          ==>  time step
  **  Output:
  **    r[3],v[3]   ==>  final position and velocity in Democratic coordinate
  **    iflg        ==>  integer flag (zero if satisfactory)
  **					      (non-zero if nonconvergence)
  ** Authors:  Hal Levison & Martin Duncan  
  ** Date:    2/10/93
  ** Last revision: April 6/93 - MD adds dt and keeps dt0 unchanged
  ** Transferred to pkdgrav2: RM 06/07
  */
  
  int iflg =0;
  double d, v2, u, alpha;
  double a, asq, en, ec, es, esq, dm;
  double fchk, fp, c1, c2, c3, xkep, s, c;
  double temp, f, g, fdot, gdot; 
  int i;
  
  d = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
  v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
  u = r[0]*v[0] + r[1]*v[1] + r[2]*v[2];
  alpha = 2.0*mu/d - v2; /* -2.0*energy */
  
  if (alpha > 0.0){ 
    a = mu/alpha;
    asq = a*a;
    en = sqrt(mu/(a*asq)); /* mean motion */
    ec = 1.0 - d/a;      /* ecosE */
    es = u/(en*asq);       /* esinE */
    esq = ec*ec + es*es;   /* e^2 */
    dm = dt*en - floor(0.5*dt*en*M_1_PI)*2*M_PI;  /* delta (mean anormaly) */
    dt = dm/en;

    /* printf("check1 a = %e esq = %e \n",a,esq);*/
    /* printf("a=%e, ecc=%e, dm=%e, dt =%e \n",a,sqrt(esq),dm,dt);*/
    /*printf("check1 a = %e esq = %e \n",a,esq); */

    if((dm*dm > 0.16) || (esq > 0.36)) goto universal;
    
    if(esq*dm*dm < 0.0016){
	/* printf("enter drift_kepmd \n");*/
	drift_kepmd(dm,es,ec,&xkep,&s,&c); /* xkep: increment in eccentric anomaly */
	fchk = (xkep - ec*s +es*(1.0-c) - dm); /* residure */      
      
	/*printf("xkep = %e, s = %e, c = %e\n", xkep,s,c);*/
	
	if(fchk*fchk > DANBYB){
	    /*printf("drift_kepmd failed: fchk = %e, danby05 = %e\n", 
	  fchk, sqrt(DANBYB));*/
	    return(1);
	}      
	fp = 1.0 - ec*c + es*s;
	f = (a/d)*(c-1.0) + 1.0;
	g = dt + (s-xkep)/en;
	fdot = - (a/(d*fp))*en*s;
	gdot = (c-1.0)/fp + 1.0;
      for (i=0;i<3;i++) {
	  temp = f*r[i]+g*v[i];
	  v[i] = fdot*r[i]+gdot*v[i];
	  r[i] = temp;
      }
      return(0);
    }
    
  } /* alpha>0 */
  
 universal:
  /* printf("enter drift_kepu \n");*/
  iflg = drift_kepu(dt,d,mu,alpha,u,&fp,&c1,&c2,&c3);
  /* printf("exit drift_kepu, iflg = %d \n", iflg);*/
  if(iflg == 0){
      f = 1.0 - (mu/d)*c2;
      g = dt - mu*c3;
      fdot = -(mu/(fp*d))*c1;
      gdot = 1.0 - (mu/fp)*c2;
      for (i=0;i<3;i++) {
	temp = f*r[i]+g*v[i];
	v[i] = fdot*r[i]+gdot*v[i];
	r[i] = temp;
      }
  } 
 
  return(iflg);
}

void drift_kepmd(double dm, double es,double ec,double *px,double *ps,double *pc)
{
/*  Subroutine for solving kepler's equation in difference form for an
 **  ellipse, given SMALL dm and SMALL eccentricity.  See DRIFT_DAN.F
 ** for the criteria.
 ** WARNING - BUILT FOR SPEED : DOES NOT CHECK HOW WELL THE ORIGINAL
 ** EQUATION IS SOLVED! (CAN DO THAT IN THE CALLING ROUTINE BY
 ** CHECKING HOW CLOSE (x - ec*s +es*(1.-c) - dm) IS TO ZERO.

 **	Input:
 **	    dm		==> increment in mean anomaly M
 **	    es,ec       ==> ecc. times sin and cos of E_0 

 **      Output:
 **           x          ==> increment in eccentric anomaly M
 **           s,c        ==> sin and cosine of x 
 */
  double x, s, c, fac1, fac2, q; 
  double f, fp, fpp, fppp, dx, y;

  double A0 = 39916800.0;
  double A1 = 6652800.0;
  double A2 = 332640.0;
  double A3 = 7920.0;
  double A4 = 110.0;

  /*    calc initial guess for root */
  fac1 = 1.0/(1.0 - ec);
  q = fac1*dm;
  fac2 = es*es*fac1 - ec/3.0;
  x = q*(1.0 -0.5*fac1*q*(es -q*fac2));
  /*  excellent approx. to sin and cos of x for small x */
  y = x*x;
  s = x*(A0-y*(A1-y*(A2-y*(A3-y*(A4-y)))))/A0;
  c = sqrt(1.0 - s*s);

  /* Compute better value for the root using quartic Newton method */
  f = x - ec*s + es*(1.0-c) - dm;
  fp = 1.0 - ec*c + es*s;
  fpp = ec*s + es*c;
  fppp = ec*c - es*s;
  dx = -f/fp;
  dx = -f/(fp + 0.5*dx*fpp);
  dx = -f/(fp + 0.5*dx*fpp + 0.16666666666666666*dx*dx*fppp);
  x = x + dx;

  /*  excellent approx. to sin and cos of x for small x */
  y = x*x;
  s = x*(A0-y*(A1-y*(A2-y*(A3-y*(A4-y)))))/A0;
  c = sqrt(1.0 - s*s);

  *px = x;
  *ps = s;
  *pc = c;
}

int drift_kepu(double dt,double d,double mu,double alpha,double u,double *pfp,
	   double *pc1,double *pc2,double *pc3){
  /* subroutine for solving kepler's equation using universal variables.
  **
  **             Input:
  **                 dt            ==>  time step
  **                 d             ==>  Distance between `Sun' and paritcle
  **                 mu            ==>  solar mass 
  **                 alpha         ==>  -2.0*energy 
  **                 u             ==>  angular momentun  
  **             Output:
  **                 fp            ==>  f' from p170                                     
  **                 c1,c2,c3      ==>  c's from p171-172  
  **                 iflg          ==>  =0 if converged; !=0 if not
  ** Author:  Hal Levison  
  ** revision: 2/3/93
  */
  double s, st, fo, fn;
  double fp, c1,c2,c3;
  int iflg;
    
  
  /* initial guess */
  s = drift_kepu_guess(dt,d,mu,alpha,u);

  /* store initial guess for possible use later in 
     Laguerre's method, in case newton's method fails. */
  st = s;
  /* printf("enter drift_kepu_new \n");*/
  iflg = drift_kepu_new(&s,dt,d,mu,alpha,u,&fp,&c1,&c2,&c3);
  /* printf("exit drift_kepu_new, iflg = %d \n",iflg);*/
  if(iflg != 0){
    fo = drift_kepu_fchk(dt,d,mu,alpha,u,st);
    fn = drift_kepu_fchk(dt,d,mu,alpha,u,s);
    if(fabs(fo) < fabs(fn)) s = st;
    /* printf("enter drift_kepu_lag \n");*/
    iflg = drift_kepu_lag(s,dt,d,mu,alpha,u,&fp,&c1,&c2,&c3);
    /* printf("exit drift_kepu_lag, iflg = %d \n",iflg); */
  }
  *pfp = fp;
  *pc1 = c1;
  *pc2 = c2;
  *pc3 = c3;

  return(iflg);
}

double drift_kepu_guess(double dt,double d,double mu,double alpha,double u){
  /* Initial guess for solving kepler's equation using universal variables
  **           Input:
  **               dt            ==>  time step 
  **               d             ==>  Distance between `Sun' and paritcle
  **               mu            ==>  solar mass 
  **               alpha         ==>  energy
  **               u             ==>  angular momentun 
  **           Output:
  **               s             ==>  initial guess for the value of 
  **                                  universal variable
  ** Author:  Hal Levison & Martin Duncan 
  ** Date:    3/12/93
  ** Last revision: April 6/93
  ** Modified by JEC: 8/6/98
  */

  double s, a, en, ec, es, e, y, x;
  double sy, cy;
  int iflg;

  if (alpha > 0.0){ /* find initial guess for elliptic motion */
    if( dt/d <= 0.4){
      s = (dt/d)*(1.0 - (dt*u)/(2.0*d*d));
      return(s);
    }else{
      a = mu/alpha;
      en = sqrt(mu/(a*a*a));
      ec = 1.0 - d/a;
      es = u/(en*a*a);
      e = sqrt(ec*ec + es*es);
      y = en*dt - es;       
      y = fmod(y,2*M_PI);
      if (y < 0) y += 2*M_PI; 
      cy = cos(y);
      sy = sqrt(1.0 - cy*cy);
      if (y > M_PI) sy *= -1.0;
   
      if((es*cy + ec*sy) > 0) x = y + 0.85*e;
      else x = y - 0.85*e;
      s = x/sqrt(alpha);
    }
  }else{ /* find initial guess for hyperbolic and parabolic motion */
    iflg = drift_kepu_p3solve(dt,d,mu,alpha,u,&s);
    if(iflg != 0) s = dt/d;
  }
  return(s);
}

int drift_kepu_p3solve(double dt,double d,double mu,double alpha,double u,double *ps){
  /* Returns the real root of cubic often found in solving kepler
  ** problem in universal variables.
  **
  **           Input:
  **               dt            ==>  time step 
  **               d             ==>  Distance between `Sun' and paritcle
  **               mu            ==>  solar
  **               alpha         ==>  -2.0*energy 
  **               u             ==>  Vel. dot radial vector 
  **           Output:
  **               s             ==>  solution of cubic eqn for the  
  **                                  universal variable
  **               iflg          ==>  success flag ( = 0 if O.K.) (integer)

  ** Author:  Martin Duncan  
  ** Date:    March 12/93
  ** Last revision: March 12/93
  */
    double denom, a0, a1, a2;
    double q, r, sq2, sq, p1, p2;
    int iflg;

  denom = (mu - alpha*d)/6.0;
  a2 = 0.5*u/denom;
  a1 = d/denom;
  a0 =-dt/denom;
 
  q = (a1 - a2*a2/3.0)/3.0;
  r = (a1*a2 -3.0*a0)/6.0 - (a2*a2*a2)/27.0;
  sq2 = q*q*q + r*r;

  if(sq2 > 0.0){
    sq = sqrt(sq2);

    if ((r+sq) <= 0.0){
      p1 =  -pow(-(r + sq), 1.0/3.0);
    }else{
      p1 = pow(r + sq, 1.0/3.0);
    }
	  
    if ((r-sq) <= 0.0){
      p2 =  -pow(-(r - sq), 1.0/3.0);
    }else{
      p2 = pow(r - sq, 1.0/3.0);
    }
    
    iflg = 0;
    *ps = p1 + p2 - a2/3.0;
    
  }else{
    iflg = 1;
    *ps = 0;
  }
  return(iflg);
}
   
int drift_kepu_new(double *ps,double dt,double d,double mu,double alpha,double u,
		   double *pfp,double *pc1,double *pc2,double *pc3){
/* subroutine for solving kepler's equation in universal variables.
** using NEWTON'S METHOD
**
**             Input:
**                 s             ==>  inital value of universal variable
**                 dt            ==>  time step 
**                 d             ==>  Distance between `Sun' and paritcle
**                 mu            ==>  solar mass
**                 alpha         ==>  -2.0*energy 
**                 u             ==>  angular momentun  
**             Output:
**                 s             ==>  final value of universal variable
**                 fp            ==>  f' from p170  
**                 c1,c2,c3      ==>  c's from p171-172
**                 iflgn         ==>  =0 if converged; !=0 if not
**
** Author:  Hal Levison  
** Date:    2/3/93
** Last revision: 4/21/93
** Modified by JEC: 31/3/98
*/
  int nc;
  double s = *ps;
  double s2, x, c0, c1, c2, c3;
  double f, fp, fpp, fppp, ds, fdt;

  /* printf("in drift_kepu_new \n");*/


     for(nc=0;nc<6;nc++){
       s2 = s*s;
       x = s2*alpha;
       drift_kepu_stumpff(x,&c0,&c1,&c2,&c3);
       c1 *= s; 
       c2 *= s2; 
       c3 *= s*s2;
       f = d*c1 + u*c2 + mu*c3 - dt;
       fp = d*c0 + u*c1 + mu*c2;
       fpp = (mu - d*alpha)*c1 + u*c0;
       fppp = (mu - d*alpha)*c0 - u*alpha*c1;
       ds = - f/fp;
       ds = - f/(fp + 0.5*ds*fpp);
       ds = -f/(fp + 0.5*ds*fpp + ds*ds*fppp*0.1666666666666667);
       s += ds;
       fdt = f/dt;

       /* quartic convergence */
       if(fdt*fdt < DANBYB*DANBYB){ /* newton's method succeeded */
	 *ps = s;
	 *pfp = fp;
	 *pc1 = c1;
	 *pc2 = c2;
	 *pc3 = c3;
	 return(0);
       }       
     }
     
     /* newton's method failed */
     *ps = s;
     *pfp = fp;
     *pc1 = c1;
     *pc2 = c2;
     *pc3 = c3;   
     return(1);
}

double drift_kepu_fchk(double dt,double d,double mu,double alpha,double u,double s){
  /* Returns the value of the function f of which we are trying to find the root
  ** in universal variables.
  **
  **             Input:
  **                 dt            ==>  time step
  **                 d             ==>  Distance between `Sun' and particle
  **                 mu            ==>  solar mass  
  **                 alpha         ==>  Twice the binding energy 
  **                 u             ==>  Vel. dot radial vector
  **                 s             ==>  Approx. root of f 
  **             Output:
  **                 f             ==>  function value ( = 0 if O.K.) (integer)
  ** Author:  Martin Duncan  
  ** Date:    March 12/93
  ** Last revision: March 12/93
  */
  double c0, c1, c2, c3, f, x; 
  x=s*s*alpha;
  drift_kepu_stumpff(x,&c0,&c1,&c2,&c3);
  c1 *= s;
  c2 *= s*s;
  c3 *= s*s*s;
  f = d*c1 + u*c2 + mu*c3 - dt;
  return(f);
}

int drift_kepu_lag(double s,double dt,double d,double mu,double alpha,double u,
		   double *pfp,double *pc1,double *pc2,double *pc3){
  /* subroutine for solving kepler's equation in universal variables.
  ** using LAGUERRE'S METHOD
  **
  **           Input:
  **               s             ==>  inital value of universal variable
  **               dt            ==>  time step
  **               d             ==>  Distance between `Sun' and paritcle
  **               mu            ==>  solar mass
  **               alpha         ==>  -2.0*energy 
  **               u             ==>  angular momentun 
  **           Output:
  **               fp            ==>  f' from p170  
  **               c1,c2,c3      ==>  c's from p171-172
  **               iflgn         ==>  =0 if converged; !=0 if not

  ** Author:  Hal Levison  
  ** Date:    2/3/93
  ** Last revision: 4/21/93
  */

  int NTMP = NLAG2+1;
  int nc, ncmax;
  double ln, x, c0, c1, c2, c3;
  double f, fp, fpp, sigma, ds, fdt;

  /* To get close approch needed to take lots of iterations if alpha<0 */
  ncmax = NLAG2;
   
  ln = 5.0;
  /* start laguere's method */
  for (nc =0;nc<ncmax;nc ++){
    x = s*s*alpha;
    drift_kepu_stumpff(x,&c0,&c1,&c2,&c3);
    c1 *= s; 
    c2 *= s*s; 
    c3 *= s*s*s;
    f = d*c1 + u*c2 + mu*c3 - dt;
    fp = d*c0 + u*c1 + mu*c2;
    fpp = (-40.0*alpha + mu)*c1 + u*c0;
    if(fp > 0) sigma = 1.0;
    else sigma = -1.0;
    ds = - ln*f/(fp + sigma*sqrt(fabs((ln - 1.0)*(ln - 1.0)*fp*fp 
					      - (ln - 1.0)*ln*f*fpp)));
    s += ds;
    fdt = f/dt;

    /* quartic convergence */
    if( fdt*fdt < DANBYB*DANBYB){
      *pfp = fp;
      *pc1 = c1;
      *pc2 = c2;
      *pc3 = c3;	
      return(0);
    } 
    /* Laguerre's method succeeded */
  }
  *pfp = fp;
  *pc1 = c1;
  *pc2 = c2;
  *pc3 = c3;	
  return(2);
}
    

void drift_kepu_stumpff(double x,double *pc0,double *pc1,double *pc2,double *pc3){
  /* subroutine for the calculation of stumpff functions
  ** see Danby p.172  equations 6.9.15
  **
  **          Input:
  **               x             ==>  argument
  **           Output:
  **               c0,c1,c2,c3   ==>  c's from p171-172
                                
  ** Author:  Hal Levison  
  ** Date:    2/3/93
  ** Last revision: 2/3/93
  ** Modified by JEC: 31/3/98
  */

  int n,i;
  double xm,x2,x3,x4,x5,x6;
  double c0,c1,c2,c3;
 
  n = 0;
  xm = 0.1;
  while(fabs(x) >= xm){
    n ++;
    x *= 0.25;
  }

  x2 = x*x;
  x3 = x*x2;
  x4 = x2*x2;
  x5 = x2*x3;
  x6 = x3*x3;

  c2 = 1.147074559772972e-11*x6 - 2.087675698786810e-9*x5
    + 2.755731922398589e-7*x4  - 2.480158730158730e-5*x3
    + 1.388888888888889e-3*x2  - 4.166666666666667e-2*x + 0.5;

  c3 = 7.647163731819816e-13*x6 - 1.605904383682161e-10*x5
    + 2.505210838544172e-8*x4  - 2.755731922398589e-6*x3
    + 1.984126984126984e-4*x2  - 8.333333333333333e-3*x
    + 1.666666666666667e-1;

  c1 = 1.0 - x*c3;
  c0 = 1.0 - x*c2;

  if(n != 0){
    for (i=n;n<1;n--){
      c3 = (c2 + c0*c3)*0.25;
      c2 = c1*c1*0.5;
      c1 = c0*c1;
      c0 = 2.*c0*c0 - 1.0;
      x = x * 4.0;
    }
  }   
  *pc0 = c0; 
  *pc1 = c1; 
  *pc2 = c2; 
  *pc3 = c3; 
}

/* void mco_sine (double x,double *psx,double *pcx){
   x %=  2*M_PI;
   if (x < 0) x += 2*M_PI;
   *pcx = cos(x);
   *psx = sqrt(1.d0 - (*pcx)*(*pcx));
   if (x > M_PI) *psx *= -1.0;
   }*/

#endif /* SYMBA */
