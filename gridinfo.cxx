#include "gridinfo.hpp"
using namespace gridinfo;

GridInfo::GridInfo(MDL mdl,MDLFFT fft) {
    m_grid = fft->rgrid->n1, fft->rgrid->n2, fft->rgrid->n3;
    m_nlocal = fft->rgrid->nLocal;
    m_sz = fft->rgrid->sSlab;
    m_nz = fft->rgrid->nSlab;
    m_sy = fft->kgrid->sSlab;
    m_ny = fft->kgrid->nSlab;
    m_iCore = mdlCore(mdl);
    m_nCore = mdlCores(mdl);
    }


/*
** This will change the given array to be a reference to the data
** buffer with the correct geometry.
*/
void GridInfo::setupArray(real_t *dataFirst,real_array_t &rspace) {
    // Make a view of the raw data with the proper extents

    // The raw real array has the following properties:
    // - The x dimension is larger than the grid size: a1r()
    // - The z dimension is a slab on this processor: sz() -> ez()
    //   but it has dimensions 0 -> nz() to avoid overflowing index calculations
    // Note that we may end up with zero elements on some processors! nz()==0
    if (nz()) {
        real_array_t rawr(
            dataFirst,
            blitz::shape(a1r(),n2(),nz()), blitz::neverDeleteData,
            RegularArray());
        real_array_t rspace2=
            rawr(blitz::Range(0,n1r()-1),
                blitz::Range(0,n2()-1),
                blitz::Range(0,nz()-1));
        rspace.reference(rspace2);
        }
    // Create an empty array; note that the data pointer is NULL
    else {
        real_array_t rspace2;
        rspace.reference(rspace2);
        }
    }

//pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
//pthread_mutex_lock(&mutex);
//pthread_mutex_unlock(&mutex);

void GridInfo::setupArray(real_t *dataFirst,complex_array_t &kspace) {
    // The complex array shares the storage of the real array, but:
    // - The x dimension is 1 greater than half the grid size
    // - The z and y dimensions are transposed.
    // - The y dimension is a slab: sy() -> ey()
    //   but has dimensions 0 -> ny()
    auto nz = (n3()+m_nCore-1) / m_nCore;
    auto sz = m_iCore * nz;
    auto ez = sz + nz;
    if (sz >= n3()) sz = ez = 0;
    else if (ez > n3()) ez = n3();
    nz = ez - sz;
    if (ny() && nz) {
        complex_array_t rawk( // Raw dimensions. May include holes
            reinterpret_cast<complex_t *>(dataFirst),
            blitz::shape(a1k(),ny(),n3()), blitz::neverDeleteData,
            TransposedArray());
        complex_array_t kspace2 = // Just the valid parts for our core
            rawk(blitz::Range(0,n1k()-1),
                blitz::Range(0,ny()-1),
                blitz::Range(sz,ez-1));
        kspace.reference(kspace2);
        kspace.reindexSelf(dimension_t(0,0,sz)); // Correct "z" dimension
        }
    // Create an empty array; note that the data pointer is NULL
    else {
        complex_array_t kspace2;
        kspace.reference(kspace2);
        }
    }
