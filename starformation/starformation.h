#include "master.h"
#ifdef STELLAR_EVOLUTION
#include "stellarevolution/stellarevolution.h"
#endif

/*
 * ---------------------
 * MAIN FUNCTIONS
 * ---------------------
 */
void msrStarFormInit(MSR msr, double dTime);
int pstStarFormInit(PST,void *,int,void *,int);
void pkdStarFormInit(PKD pkd, double dTime, int *nFormed);

void msrStarForm( MSR, double, double, int);
int pstStarForm(PST,void *,int,void *,int);
void pkdStarForm(PKD pkd, double dTime, double dDelta, double dScaleFactor,
      double dDenMin, int *nFormed, double *dMassFormed, int *nDeleted);


/*
 * ---------------------
 * STRUCTURE DEFINITIONS
 * ---------------------
 */
struct inStarForm
    {
    double dTime;
    double dScaleFactor;
    double dDenMin;
    double dDelta;
    };

struct outStarForm
    {
    int nFormed;
    int nDeleted;
    double dMassFormed;
    };

/*
 * ---------------------
 * HELPER FUNCTIONS
 * ---------------------
 */
static inline double pressure_SFR(PKD pkd, double a_m3, double dDenMin,
      PARTICLE *p, SPHFIELDS *psph);
static inline double density_SFR(PKD pkd, double a_m3, double dDenMin,
      PARTICLE *p, SPHFIELDS *psph);
