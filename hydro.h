/* File added by Isaac Alonso for computing
 * the hydrodynamical part using a mesh-free
 * method, following the work of Hopkins 2015
 */




double cubicSplineKernel(double r, double h);
void inverseMatrix(double* E, double* B);

void initHydroForces(void *vpkd, void *vp);
void initHydroForcesCached(void *vpkd, void *vp);
void combHydroForces(void *vpkd, void *p1,void *p2);
void hydroGradients(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf);
void hydroRiemann(PARTICLE *p,float fBall,int nSmooth,NN *nnList,SMF *smf);
