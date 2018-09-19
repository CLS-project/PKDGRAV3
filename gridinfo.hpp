#include "blitz/array.h"
#include "mdl.h"

namespace gridinfo {
typedef blitz::TinyVector<int,3>    dimension_t;

typedef float                       real_t;
typedef blitz::Array<real_t,3>      real_array_t;
typedef blitz::Array<real_t,2>      real_slice_t;
typedef blitz::Array<real_t,1>      real_vector_t;
typedef std::complex<real_t>        complex_t;
typedef blitz::Array<complex_t,3>   complex_array_t;
typedef blitz::Array<complex_t,2>   complex_slice_t;
typedef blitz::Array<complex_t,1>   complex_vector_t;

// This lays out the r-space array in normal order.
class RegularArray : public blitz::GeneralArrayStorage<3> {
public:
    RegularArray(int kstart=0)
        : blitz::GeneralArrayStorage<3>(noInitializeFlag()) {
        ordering_(0) = 0;
        ordering_(1) = 1;
        ordering_(2) = 2;
        ascendingFlag_ = true;
        base_[0] = 0;
        base_[1] = 0;
        base_[2] = kstart;
	}
    };

// This lays out the k-space array in transposed order.
class TransposedArray : public blitz::GeneralArrayStorage<3> {
public:
    TransposedArray(int jstart=0)
        : blitz::GeneralArrayStorage<3>(noInitializeFlag()) {
        ordering_(0) = 0;
        ordering_(1) = 2;
        ordering_(2) = 1;
        ascendingFlag_ = true;
        base_[0] = 0;
        base_[1] = jstart;
        base_[2] = 0;
	}
    };

class GridInfo {
    dimension_t m_grid; // full dimensions of the grid
    int m_sz, m_nz;           // Start and number of of "z" slabs
    int m_sy, m_ny;           // Ditto for y after transpose
    int m_nlocal;             // Size of local array in elements    
    int m_iCore, m_nCore;
public:
    const dimension_t &grid()   const { return m_grid; }

    inline int sz()             const { return m_sz; }
    inline int nz()             const { return m_nz; }
    inline int ez()             const { return m_sz+m_nz; }

    inline int sy()             const { return m_sy; }
    inline int ny()             const { return m_ny; }
    inline int ey()             const { return m_sy+m_ny; }

    inline int n1r()            const { return grid()[0]; }
    inline int n1k()            const { return n1r()/2 + 1; }
    inline int n2()             const { return grid()[1]; }
    inline int n3()             const { return grid()[2]; }
    inline uint64_t nt()        const { return (uint64_t)grid()[0]*grid()[1]*grid()[2]; }

    inline int a1k()            const { return n1k(); }
    inline int a1r()            const { return 2*a1k(); }

    void setupArray(real_t *dataFirst,real_array_t &rspace);
    void setupArray(real_t *dataFirst,complex_array_t &kspace);

    explicit GridInfo(MDL mdl,MDLFFT fft);
    };
}