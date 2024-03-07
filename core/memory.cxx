#include "core/memory.h"
#include "core/fftsizes.h"
#include "pst.h"

// Constructor with MSR, grid, and count
EphemeralMemory::EphemeralMemory(mdl::mdlClass *mdl, int grid, int count) {
    if (grid>0 && count>0) {
        ServiceFftSizes::input inFFTSizes;
        ServiceFftSizes::output outFFTSizes;
        inFFTSizes.nx = inFFTSizes.ny = inFFTSizes.nz = grid;
        mdl->RunService(PST_GETFFTMAXSIZES,sizeof(inFFTSizes),&inFFTSizes,&outFFTSizes);
        this->per_process += outFFTSizes.nMaxLocal*sizeof(FFTW3(real)) * count;
    }
}
