#include "mdlgpu.h"

namespace mdl {
namespace gpu {

int Client::flushCompleted() {
    while (!done.empty()) {
        auto &M = done.dequeue();
        M.finish();
        --nGPU;
    }
    return nGPU;
}

} // namespace gpu
} // namespace mdl
