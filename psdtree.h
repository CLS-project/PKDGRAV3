#ifndef PSDTREE_HINCLUDED
#define PSDTREE_HINCLUDED

#include <stdint.h>
#include <string.h>

#include "mdl.h"
#include "pst.h"
#ifndef HAVE_CONFIG_H
#include "floattype.h"
#endif

#ifdef NEW_TREE
#define MAX_NBUCKET 5

typedef struct psdNew {
    float s;  /* scale of the cell, if -ve then it indicates a bucket! */
    union celltype {
    struct cell {
	double dSplit;
	uint32_t iLower;
	uint32_t iUpper;
	uint16_t idLower;
	uint16_t idUpper;
	} c;
    struct bucket {
	uint32_t iPart[MAX_NBUCKET];
	} b;
    };
    double r[6];
    } PSDNEW;
#endif

static inline int max_bnd(const double *fMax)
{
    int i;
    int d = 0;
    double max = fMax[0];
    for (i=1; i < 3; i++)
	if (fMax[i] > max) { max = fMax[i]; d = i; }
    return d;
}

static inline int max_bnd_range(const double *fMax, const int m, const int M)
{
    int i;
    int d = m;
    double max = fMax[m];
    for (i=m+1; i < M; i++)
	if (fMax[i] > max) { max = fMax[i]; d = i; }
    if (max == 0) 
	return -1;
    return d;
}

static inline int min_bnd(const double *fMax)
{
    int i;
    int d = 0;
    double min = fMax[0];
    for (i=1; i < 3; i++)
	if (fMax[i] < min) { min = fMax[i]; d = i; }
    return d;
}

static inline int max_side(const double *fMax)
{
    return 2 * max_bnd(fMax);
}

static inline int min_side(const double *fMax)
{
    return 2 * min_bnd(fMax);
}

void psdBuildTree(PKD pkd, PSX psx, struct inPSD *in,KDN *pkdn);
int psdTreeDepth(PKD pkd);

#endif
