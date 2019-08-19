/* File added by Isaac Alonso for computing
 * the hydrodynamical part using a mesh-free
 * method, following the work of Hopkins 2015
 */




double cubicSplineKernel(double r, double h);
void inverseMatrix(double* E, double* B);

void initHydroLoop(void *vpkd, void *vp);
void initHydroLoopCached(void *vpkd, void *vp);
void initHydroFluxes(void *vpkd, void *vp);
void initHydroFluxesCached(void *vpkd, void *vp);
void initHydroGradients(void *vpkd, void *vp);
void combFirstHydroLoop(void *vpkd, void *p1,void *p2);
void combSecondHydroLoop(void *vpkd, void *p1,void *p2);
void combThirdHydroLoop(void *vpkd, void *p1,void *p2);
void hydroDensity(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf);
void hydroGradients(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf);
void hydroRiemann(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf);
void hydroStep(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf);
void BarthJespersenLimiter(double* limVar, double* gradVar, double var_max, double var_min, double dx, double dy, double dz);
void ConditionedBarthJespersenLimiter(double* limVar, double* gradVar, double var_max, double var_min, double dx, double dy, double dz, double Ncrit, double Ncond);
