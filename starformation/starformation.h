#include "master.h"

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

