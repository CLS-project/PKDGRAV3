#ifndef HYDRO_HINCLUDED
#define HYDRO_HINCLUDED

#include "master.h"
#include "pkd.h"
#include "master.h"
#include "smoothfcn.h"

#ifdef OPTIM_FLUX_VEC
/* When doing SIMD, a structure of arrays is used,
 * rather than an array of structures.
 *
 * To simplify the indexing of the elements, these enum should be
 * always used
 */
enum q_data{
      q_mass,
      q_ball,
      q_dx, q_dy, q_dz, q_dr,
      q_rung,
      q_rho,
      q_P,
#ifdef ENTROPY_SWITCH
      q_S,
#endif
      q_vx, q_vy, q_vz,
      q_gradRhoX, q_gradRhoY, q_gradRhoZ,
      q_gradPX, q_gradPY, q_gradPZ,
      q_gradVxX, q_gradVxY, q_gradVxZ,
      q_gradVyX, q_gradVyY, q_gradVyZ,
      q_gradVzX, q_gradVzY, q_gradVzZ,
      q_lastAccX, q_lastAccY, q_lastAccZ,
      q_lastUpdateTime,
      q_B_XX, q_B_YY, q_B_ZZ, q_B_XY, q_B_XZ, q_B_YZ,
      q_omega,
      q_last};

enum FLUX_OUT{
      out_minDt,
      out_Frho,
      out_FmomX,out_FmomY,out_FmomZ,
      out_Fene,
#ifdef ENTROPY_SWITCH
      out_FS,
#endif
      out_last};
#endif // OPTIM_FLUX_VEC

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

/* Time step loop */
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
