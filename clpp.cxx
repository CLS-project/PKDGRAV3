#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "clutil.h"

static const char *PP_kernel_src =
    "";

/* Compile/create ewald kernel */
void CLkernelPP(CLCTX cl) {
    cl_int rc;
    cl->context->programEwald = CL_compile(cl,PP_kernel_src);
//    cl->kernelEwald = clCreateKernel(cl->context->programEwald, "PP", &rc);
    assert(rc == CL_SUCCESS);
    }
