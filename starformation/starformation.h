#include "pst.h"
#include "units.h"
#ifdef __cplusplus
extern "C" {
#endif

/*
 * ---------------------
 * STRUCTURE DEFINITIONS
 * ---------------------
 */
struct outStarForm {
    int nFormed;
    int nDeleted;
    double dMassFormed;
};
struct inStarFormInit {
#ifdef FEEDBACK
    double dTime;
    int bCCSNFeedback;
    int bSNIaFeedback;
    double dCCSNFBDelay;
    double dSNIaFBDelay;
#endif
};

/*
 * ---------------------
 * MAIN FUNCTIONS
 * ---------------------
 */
int pstStarFormInit(PST,void *,int,void *,int);
void pkdStarFormInit(PKD pkd, struct inStarFormInit in, int *nFormed);

int pstStarForm(PST,void *,int,void *,int);
void pkdStarForm(PKD pkd, struct inStarForm in,
                 int *nFormed, double *dMassFormed, int *nDeleted);

#ifdef __cplusplus
}
#endif
