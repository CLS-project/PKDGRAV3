#include "pst.h"
#include "units.h"
#ifdef __cplusplus
extern "C" {
#endif
/*
 * ---------------------
 * MAIN FUNCTIONS
 * ---------------------
 */
int pstStarFormInit(PST,void *,int,void *,int);
void pkdStarFormInit(PKD pkd, double dTime, double dSNFBDelay, int *nFormed);

int pstStarForm(PST,void *,int,void *,int);
void pkdStarForm(PKD pkd, struct inStarForm in, 
                 int *nFormed, double *dMassFormed, int *nDeleted);


/*
 * ---------------------
 * STRUCTURE DEFINITIONS
 * ---------------------
 */
struct outStarForm
    {
    int nFormed;
    int nDeleted;
    double dMassFormed;
    };
struct inStarFormInit
    {
    double dTime;
    double dSNFBDelay;
    };

#ifdef __cplusplus
}
#endif
