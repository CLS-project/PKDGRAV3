
#include "pkd.h"
#include "smooth/smoothfcn.h"


#define XX 0
#define YY 3
#define ZZ 5
#define XY 1
#define XZ 2
#define YZ 4

typedef double my_real;

/* -----------------
 * MAIN FUNCTIONS
 * -----------------
 */

/* Density loop */
void hydroDensity(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf);
void hydroDensity_node(PKD pkd, SMF *smf, BND bnd_node, PARTICLE **sinks, NN *nnList,
                       int nCnt_own, int nCnt);

/* Gradient loop */
void hydroGradients(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf);

/* Flux loop */
void initHydroFluxes(void *vpkd, void *vp);
void initHydroFluxesCached(void *vpkd, void *vp);
void hydroRiemann(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf);
void hydroRiemann_vec(PARTICLE *p,float fBall,int nSmooth,
                      my_real **restrict input_buffer,
                      my_real **restrict output_buffer,
                      SMF *smf);
void pkdResetFluxes(PKD pkd,double dTime,double dDelta,double,double);

void combThirdHydroLoop(void *vpkd, void *p1,const void *p2);
void hydroFluxFillBuffer(my_real **buffer, PARTICLE *q, int i,
                         double dr2, double dx, double dy, double dz, SMF *);
void hydroFluxUpdateFromBuffer(my_real **out_buffer, my_real **in_buffer,
                               PARTICLE *p, PARTICLE *q, int i, SMF *);
void hydroFluxGetNvars(int *in, int *out);

/* Time step loop */
void initHydroStep(void *vpkd, void *vp);
void hydroStep(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf);
void combHydroStep(void *vpkd, void *p1,const void *p2);
void pkdWakeParticles(PKD pkd,int iRoot, double dTime, double dDelta);


/* Source terms */
void hydroSourceGravity(PKD pkd, PARTICLE *p, SPHFIELDS *psph,
                        double pDelta, double *pa, double dScaleFactor,
                        int bComove);
void hydroSourceExpansion(PKD pkd, PARTICLE *p, SPHFIELDS *psph,
                          double pDelta, double dScaleFactor, double dHubble,
                          int bComove, double dConstGamma);
void hydroSyncEnergies(PKD pkd, PARTICLE *p, SPHFIELDS *psph, double pa[3],
                       double dConstGamma);

void hydroSetPrimitives(PKD pkd, PARTICLE *p, SPHFIELDS *psph, double dTuFac, double dConstGamma);

void hydroSetLastVars(PKD pkd, PARTICLE *p, SPHFIELDS *psph, double *pa,
                      double dScaleFactor, double dTime, double dDelta,
                      double dConstGamma);

/* -----------------
 * HELPERS
 * -----------------
 */
void inverseMatrix(double *E, double *B);
double conditionNumber(double *E, double *B);
double cubicSplineKernel(double r, double h);
void BarthJespersenLimiter(double *limVar, double *gradVar,
                           double var_max, double var_min,
                           double dx, double dy, double dz);
void ConditionedBarthJespersenLimiter(double *limVar, myreal *gradVar,
                                      double var_max, double var_min,
                                      double dx, double dy, double dz,
                                      double Ncrit, double Ncond);
void genericPairwiseLimiter(double Lstate, double Rstate,
                            double *Lstate_face, double *Rstate_face);
void compute_Ustar(double rho_K, double S_K, double v_K,
                   double p_K, double h_K, double S_s,
                   double *rho_sK, double *rhov_sK, double *e_sK);


