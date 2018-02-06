cmake_minimum_required(VERSION 2.4)
project(openpa)

INCLUDE(CheckTypeSize)
check_type_size ("int" OPA_SIZEOF_INT)

INCLUDE (CheckIncludeFiles)
CHECK_INCLUDE_FILES (pthread.h OPA_HAVE_PTHREAD_H)
CHECK_INCLUDE_FILES (stddef.h OPA_HAVE_STDDEF_H)

INCLUDE (CheckFunctionExists)
CHECK_FUNCTION_EXISTS(sched_yield OPA_HAVE_SCHED_YIELD)

INCLUDE (CheckCSourceCompiles)

CHECK_C_SOURCE_COMPILES("
int foo(char *,...) __attribute__ ((format(printf,1,2)));
int main() {return 0;}
" OPA_HAVE_GCC_ATTRIBUTE)

macro(CompilePrimitive HEADER NAME)
CHECK_C_SOURCE_COMPILES("
#define OPA_SIZEOF_INT ${OPA_SIZEOF_INT}
#define OPA_SIZEOF_VOID_P  ${CMAKE_SIZEOF_VOID_P}
#ifndef _opa_inline
#define _opa_inline inline
#endif
#ifndef _opa_restrict
#define _opa_restrict restrict
#endif
#ifndef _opa_const
#define _opa_const const
#endif
#ifdef HAVE_GCC_ATTRIBUTE
#define OPA_HAVE_GCC_ATTRIBUTE 1
#endif
#include \"${CMAKE_CURRENT_SOURCE_DIR}/opa_util.h\"
#include \"${CMAKE_CURRENT_SOURCE_DIR}/primitives/${HEADER}\"
int main() {
    OPA_int_t a, b;
    int c;

    OPA_store_int(&a, 0);
    OPA_store_int(&b, 1);
    c = OPA_load_int(&a);

    OPA_add_int(&a, 10);
    OPA_incr_int(&a);
    OPA_decr_int(&a);

    c = OPA_decr_and_test_int(&a);
    c = OPA_fetch_and_add_int(&a, 10);
    c = OPA_fetch_and_incr_int(&a);
    c = OPA_fetch_and_decr_int(&a);

    c = OPA_cas_int(&a, 10, 11);
    c = OPA_swap_int(&a, OPA_load_int(&b));

    OPA_write_barrier();
    OPA_read_barrier();
    OPA_read_write_barrier();
    return 0;
}" ${NAME})
endmacro(CompilePrimitive)

CompilePrimitive("opa_gcc_intel_32_64.h" OPA_HAVE_GCC_X86_32_64)
CompilePrimitive("opa_gcc_intel_32_64_p3.h" OPA_HAVE_GCC_X86_32_64_P3)
CompilePrimitive("opa_gcc_ia64.h" OPA_HAVE_GCC_AND_IA64_ASM)
CompilePrimitive("opa_gcc_ppc.h" OPA_HAVE_GCC_AND_POWERPC_ASM)
CompilePrimitive("opa_gcc_arm.h" OPA_HAVE_GCC_AND_ARM_ASM)
CompilePrimitive("opa_gcc_sicortex.h" OPA_HAVE_GCC_AND_SICORTEX_ASM)
CompilePrimitive("opa_gcc_intrinsics.h" OPA_HAVE_GCC_INTRINSIC_ATOMICS)
CompilePrimitive("opa_nt_intrinsics.h" OPA_HAVE_NT_INTRINSICS)
CompilePrimitive("opa_sun_atomic_ops.h" OPA_HAVE_SUN_ATOMIC_OPS)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/opa_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/opa_config.h)

add_library(${PROJECT_NAME} STATIC opa_primitives.c opa_queue.c)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

install(TARGETS ${PROJECT_NAME} DESTINATION "lib")
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/opa_config.h opa_queue.h opa_primitives.h opa_util.h DESTINATION "include")
