#ifndef MDLFFT_H
#define MDLFFT_H
#include <stdint.h>
/*
** Grid Operations
*/

typedef struct mdlGridContext {
    uint32_t n1,n2,n3;     /* Real dimensions */
    uint32_t a1;           /* Actual size of dimension 1 */
    uint32_t sSlab, nSlab; /* Start and number of slabs */
    uint64_t nLocal;       /* Number of local elements */
    uint32_t *rs;  /* Starting slab for each processor */
    uint32_t *rn;  /* Number of slabs on each processor */
    uint32_t *id;  /* Which processor has this slab */
    } * MDLGRID;

typedef struct {
    MDLGRID grid;
    uint64_t II;
    int x, y, z; /* Real coordinate */
    int i;       /* Index into local array */
    } mdlGridCoord;


#ifdef MDL_FFTW
#include <fftw3-mpi.h>
#ifdef USE_SINGLE
#define FFTW3(name) fftwf_ ## name
typedef float fftwf_real;
#else
#define FFTW3(name) fftw_ ## name
typedef double fftw_real;
#endif
typedef struct mdlFFTContext {
    MDLGRID rgrid;
    MDLGRID kgrid;
    FFTW3(plan) fplan, iplan;
    } * MDLFFT;
#endif
#endif
