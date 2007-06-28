#ifdef SYMBA
#ifndef KEPLER_HINCLUDED
#define KEPLER_HINCLUDED

#define DANBYB  1.0e-13
#define NLAG2   400;

void drift_kepmd(double dm, double es,double ec,double *px,double *ps,
		 double *pc);
int drift_kepu(double dt,double d,double mu,double alpha,double u,double fp,
	       double *pc1,double *pc2,double *pc3);
double drift_kepu_guess(double dt,double d,double mu,double alpha,double u);
int drift_kepu_p3solve(double dt,double d,double mu,double alpha,
		       double u,double *ps);
int drift_kepu_new(double *ps,double dt,double d,double mu,double alpha,
		   double u,double *pfp,double *pc1,double *pc2,double *pc3);
double drift_kepu_fchk(double dt,double d,double mu,double alpha,
		       double u,double s);
int drift_kepu_lag(double s,double dt,double d,double mu,double alpha,double u,
		   double *pfp,double *pc1,double *pc2,double *pc3);
void drift_kepu_stumpff(double x,double *pc0,double *pc1,
			double *pc2,double *pc3);

#endif
#endif /* SYMBA */
