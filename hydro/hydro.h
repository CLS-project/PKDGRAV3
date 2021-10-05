#ifndef HYDRO_HINCLUDED
#define HYDRO_HINCLUDED

#include "master.h"
#include "pkd.h"
#include "master.h"
#include "smoothfcn.h"


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
void msrComputeSmoothing(MSR msr,double dTime);
void hydroDensity(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf);
void hydroDensity_node(PKD pkd, BND bnd_node, PARTICLE **sinks, NN *nnList,
                       int nCnt_own, int nCnt);

/* Gradient loop */
void msrMeshlessGradients(MSR msr,double dTime);
void hydroGradients(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf);

/* Flux loop */
void msrMeshlessFluxes(MSR msr,double dTime,double dDelta,int iRoot);
void initHydroFluxes(void *vpkd, void *vp);
void initHydroFluxesCached(void *vpkd, void *vp);
void hydroRiemann(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf);
void hydroRiemann_vec(PARTICLE *p,float fBall,int nSmooth,
                      my_real** restrict input_buffer,
                      my_real** restrict output_buffer,
                      SMF *smf);
void combThirdHydroLoop(void *vpkd, void *p1,void *p2);
void hydroFluxFillBuffer(PKD pkd, my_real **buffer, PARTICLE* q, int i,
                         double dDelta, double dr2, double dx, double dy, double dz);
void hydroFluxUpdateFromBuffer(PKD pkd, my_real **out_buffer, my_real **in_buffer,
                               PARTICLE* p, PARTICLE* q, int i, double dDelta);
void hydroFluxAllocateBuffer(my_real *input_buffer, my_real **input_pointers,
                             my_real *output_buffer, my_real**output_pointers,
                             int N);

/* Time step loop */
void msrHydroStep(MSR msr,uint8_t uRungLo,uint8_t uRungHi,double dTime);
void initHydroStep(void *vpkd, void *vp);
void hydroStep(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf);
void combHydroStep(void *vpkd, void *p1,void *p2);
void pkdWakeParticles(PKD pkd,int iRoot, double dTime, double dDelta);


/* Source terms */
void hydroSourceGravity(PKD pkd, PARTICLE* p, SPHFIELDS* psph,
                        double pDelta, double *pa, double dScaleFactor,
                        int bComove);
void hydroSourceExpansion(PKD pkd, PARTICLE* p, SPHFIELDS* psph,
                          double pDelta, double dScaleFactor, double dHubble,
                          int bComove);
void hydroSyncEnergies(PKD pkd, PARTICLE* p, SPHFIELDS* psph, double pa[3]);

void hydroSetPrimitives(PKD pkd, PARTICLE* p, SPHFIELDS* psph);

void hydroSetLastVars(PKD pkd, PARTICLE *p, SPHFIELDS *psph, double *pa,
                      double dScaleFactor, double dTime, double dDelta);

/* -----------------
 * HELPERS
 * -----------------
 */
#define SIGN(x) (((x) > 0) ? 1 : (((x) < 0) ? -1 : 0) )
#define MIN(X, Y)  ((X) < (Y) ? (X) : (Y))
#define MAX(X, Y)  ((X) > (Y) ? (X) : (Y))
void inverseMatrix(double* E, double* B);
double conditionNumber(double *E, double* B);
double cubicSplineKernel(double r, double h);
void BarthJespersenLimiter(double* limVar, double* gradVar,
                           double var_max, double var_min,
                           double dx, double dy, double dz);
void ConditionedBarthJespersenLimiter(double* limVar, myreal* gradVar,
                                      double var_max, double var_min,
                                      double dx, double dy, double dz,
                                      double Ncrit, double Ncond);
void genericPairwiseLimiter(double Lstate, double Rstate,
                            double *Lstate_face, double *Rstate_face);
void compute_Ustar(double rho_K, double S_K, double v_K,
                   double p_K, double h_K, double S_s,
                   double *rho_sK, double *rhov_sK, double *e_sK);


#endif
