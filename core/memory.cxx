#include "core/memory.h"
#include "master.h"

// Constructor with MSR, grid, and count
EphemeralMemory::EphemeralMemory(class MSR *msr, int grid, int count) {
    if (grid>0 && count>0)
        this->per_process += msr->getLocalGridMemory(grid) * count;
}
