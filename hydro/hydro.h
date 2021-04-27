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


typedef double my_real;


/*-----------------
* SMOOTH OPERATIONS
* -----------------
*/


/* Density loop */
void hydroDensity(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf);

/* Gradient loop */
void hydroGradients(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf);

/* Flux loop */
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


/*-----------------
* MAIN FUNCTIONS
* -----------------
*/


void pkdWakeParticles(PKD pkd,int iRoot, double dTime, double dDelta);
void pkdComputePrimVars(PKD pkd,int iRoot, double dTime, double dDelta);
void msrUpdatePrimVars(MSR msr,double dTime,double dDelta,int iRoot);



/*-----------------
* HELPERS
* -----------------
*/
void inverseMatrix(double* E, double* B);
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


#endif
