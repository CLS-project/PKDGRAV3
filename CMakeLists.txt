cmake_minimum_required(VERSION 3.1)
project(pkdgrav3 VERSION 3.0.4)
if(CMAKE_VERSION VERSION_GREATER "3.12.0")
   cmake_policy(SET CMP0074 NEW)
endif(CMAKE_VERSION VERSION_GREATER "3.12.0")

if (CMAKE_VERSION VERSION_GREATER "3.19")
cmake_policy(SET CMP0074 NEW)
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/mdl2)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})

include(CheckFunctionExists)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Release build with debug info selected")
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -Wall ")
  set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -Wall ")
endif()

include(CheckCXXCompilerFlag)

set(TARGET_ARCHITECTURE "auto" CACHE STRING "CPU architecture to optimize for.")
if(TARGET_ARCHITECTURE STREQUAL "auto")
  CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_OPT_ARCH_NATIVE_SUPPORTED)
  if (COMPILER_OPT_ARCH_NATIVE_SUPPORTED)
    add_compile_options(-march=native)
  endif()
  CHECK_CXX_COMPILER_FLAG("/arch:AVX" COMPILER_OPT_ARCH_AVX_SUPPORTED)
  if (COMPILER_OPT_ARCH_AVX_SUPPORTED)
    add_compile_options(/arch:AVX)
  endif()
else()
    add_compile_options(-march=${TARGET_ARCHITECTURE})
endif()
CHECK_CXX_COMPILER_FLAG("-Wall" COMPILER_OPT_WARN_ALL_SUPPORTED)
if (COMPILER_OPT_WARN_ALL_SUPPORTED)
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
  set(CMAKE_C_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
endif()
if(CMAKE_VERSION VERSION_GREATER "3.17.0")
   find_package(CUDAToolkit)
else()
   message(WARNING "Using old FindCUDA, consider updating cmake.")
   find_package(CUDA)
endif()
find_package(GSL REQUIRED)      # GNU Scientific Library
find_package(HDF5 COMPONENTS C HL)
find_package(FFTW REQUIRED)
#find_package(PythonLibs)
# _GNU_SOURCE gives us more options
INCLUDE(CheckCSourceCompiles)
check_c_source_compiles("
#include <features.h>
#ifdef __GNU_LIBRARY__
  int main() {return 0;} 
#endif
" _GNU_SOURCE)
if (_GNU_SOURCE)
  set(CMAKE_REQUIRED_DEFINITIONS -D_GNU_SOURCE)
endif()

# Check for restrict keyword
# Builds the macro A_C_RESTRICT form automake
foreach(ac_kw __restrict __restrict__ _Restrict restrict)
  check_c_source_compiles(
  "
  typedef int * int_ptr;
  int foo (int_ptr ${ac_kw} ip) {
    return ip[0];
  }
  int main(){
    int s[1];
    int * ${ac_kw} t = s;
    t[0] = 0;
    return foo(t);
  }
  "
  RESTRICT)
  if(RESTRICT)
    set(ac_cv_c_restrict ${ac_kw})
    break()
  endif()
endforeach()
if(RESTRICT)
  add_definitions("-Drestrict=${ac_cv_c_restrict}")
else()
  add_definitions("-Drestrict=")
endif()


INCLUDE (CheckIncludeFiles)
INCLUDE (CheckLibraryExists)
CHECK_INCLUDE_FILES (malloc.h HAVE_MALLOC_H)
CHECK_INCLUDE_FILES (signal.h HAVE_SIGNAL_H)
CHECK_INCLUDE_FILES (sys/time.h HAVE_SYS_TIME_H)
CHECK_INCLUDE_FILES (sys/stat.h HAVE_SYS_STAT_H)
CHECK_INCLUDE_FILES (sys/param.h HAVE_SYS_PARAM_H)
CHECK_INCLUDE_FILES (unistd.h HAVE_UNISTD_H)
CHECK_INCLUDE_FILES (inttypes.h HAVE_INTTYPES_H)
CHECK_INCLUDE_FILES(rpc/types.h HAVE_RPC_TYPES_H)
CHECK_INCLUDE_FILES(rpc/xdr.h HAVE_RPC_XDR_H)

INCLUDE(CheckSymbolExists)
check_symbol_exists(floor math.h HAVE_FLOOR)
check_symbol_exists(pow math.h HAVE_POW)
check_symbol_exists(sqrt math.h HAVE_SQRT)
check_symbol_exists(strchr string.h HAVE_STRCHR)
check_symbol_exists(strrchr string.h HAVE_STRRCHR)
check_symbol_exists(strdup string.h HAVE_STRDUP)
check_symbol_exists(strstr string.h HAVE_STRSTR)
check_symbol_exists(memmove string.h HAVE_MEMMOVE)
check_symbol_exists(memset string.h HAVE_MEMSET)
check_symbol_exists(gettimeofday sys/time.h HAVE_GETTIMEOFDAY)
check_symbol_exists(posix_memalign stdlib.h HAVE_POSIX_MEMALIGN)

check_symbol_exists(wordexp wordexp.h HAVE_WORDEXP)
check_symbol_exists(wordfree wordexp.h HAVE_WORDFREE)
check_symbol_exists(glob glob.h HAVE_GLOB)
check_symbol_exists(globfree glob.h HAVE_GLOBFREE)
check_symbol_exists(gethostname unistd.h HAVE_GETHOSTNAME)
check_symbol_exists(getpagesize unistd.h HAVE_GETPAGESIZE)
check_symbol_exists(mkdir sys/stat.h HAVE_MKDIR)
check_symbol_exists(strverscmp string.h HAVE_STRVERSCMP)

check_symbol_exists(backtrace execinfo.h USE_BT)

#AC_CHECK_FUNCS([gethrtime read_real_time time_base_to_time clock_gettime mach_absolute_time])
check_symbol_exists(atexit stdlib.h HAVE_ATEXIT)

add_subdirectory(blitz)
add_subdirectory(mdl2)

add_executable(${PROJECT_NAME} "")
target_compile_features(${PROJECT_NAME} PRIVATE cxx_auto_type cxx_lambdas)
set_property(TARGET ${PROJECT_NAME} APPEND PROPERTY COMPILE_DEFINITIONS _LARGEFILE_SOURCE)
target_sources(${PROJECT_NAME} PRIVATE
	main.c cosmo.c master.c pst.c fio.c illinois.c param.c
	pkd.c analysis.c smooth.c smoothfcn.c outtype.c output.c
	walk2.c grav2.c ewald.cxx ic.cxx tree.cxx opening.cxx pp.cxx pc.cxx cl.c
	lst.c moments.c ilp.c ilc.c iomodule.c
	fof.c hop.c group.c groupstats.c RngStream.c listcomp.c healpix.c
	tinypy.c pkdtinypy.c
	gridinfo.cxx assignmass.cxx measurepk.cxx whitenoise.cxx pmforces.cxx
        hydro/hydro.c eEOS/eEOS.c
)
add_executable(tostd tostd.c fio.c)
add_executable(psout psout.c cosmo.c)

target_link_libraries(${PROJECT_NAME} m)
target_link_libraries(tostd m)
target_link_libraries(psout m)

if(USE_BT)
target_sources(${PROJECT_NAME} PRIVATE bt.c)
endif()

CHECK_INCLUDE_FILES(libaio.h HAVE_LIBAIO_H)
if (HAVE_LIBAIO_H)
  CHECK_LIBRARY_EXISTS(aio io_setup "" HAVE_LIBAIO)
  if(HAVE_LIBAIO)
    find_library(LIBAIO_LIBRARY aio)
    target_link_libraries(${PROJECT_NAME} ${LIBAIO_LIBRARY})
  endif()
endif()
CHECK_INCLUDE_FILES(aio.h HAVE_AIO_H)
if (HAVE_AIO_H)
  CHECK_LIBRARY_EXISTS(rt aio_read "" HAVE_RT)
  if(HAVE_RT)
    target_link_libraries(${PROJECT_NAME} rt)
  endif()
endif()

if (CUDAToolkit_FOUND)
  enable_language(CUDA)
  set(USE_CUDA TRUE)
  set_target_properties(${PROJECT_NAME} PROPERTIES CUDA_ARCHITECTURES "35;60")
  target_sources(${PROJECT_NAME} PRIVATE cudaewald.cu cudapppc.cu cudautil.cu)
endif(CUDAToolkit_FOUND)
if (CUDA_FOUND)
  set(USE_CUDA TRUE)
  get_property(MDL_INCLUDES TARGET mdl2 PROPERTY INCLUDE_DIRECTORIES)
  get_property(OPENPA_INCLUDES TARGET openpa PROPERTY INCLUDE_DIRECTORIES)
#  target_link_libraries(${PROJECT_NAME} ${CUDA_LIBRARIES})
  CUDA_INCLUDE_DIRECTORIES(${MDL_INCLUDES} ${OPENPA_INCLUDES})
  CUDA_COMPILE(cuda_files cudaewald.cu cudapppc.cu cudautil.cu
	OPTIONS -arch compute_35
		-I${CMAKE_CURRENT_SOURCE_DIR}/mdl2/openpa
		-I${CMAKE_CURRENT_BINARY_DIR})
  target_sources(${PROJECT_NAME} PRIVATE ${cuda_files})
endif(CUDA_FOUND)
if (HDF5_FOUND)
  set(USE_HDF5 TRUE)
  target_include_directories(${PROJECT_NAME} PUBLIC ${HDF5_INCLUDE_DIRS})
  target_link_libraries(${PROJECT_NAME} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
  #set_property(TARGET ${PROJECT_NAME} APPEND PROPERTY COMPILE_DEFINITIONS ${HDF5_DEFINITIONS})
  #set_property(TARGET ${PROJECT_NAME} APPEND PROPERTY COMPILE_DEFINITIONS H5_USE_16_API)
  target_include_directories(tostd PUBLIC ${HDF5_INCLUDE_DIRS})
  target_link_libraries(tostd ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
  target_include_directories(psout PUBLIC ${HDF5_INCLUDE_DIRS})
  target_link_libraries(psout ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
endif(HDF5_FOUND)
if (PYTHONLIBS_FOUND)
  set(USE_PYTHON TRUE)
  target_include_directories(${PROJECT_NAME} PUBLIC ${PYTHON_INCLUDE_DIRS})
  target_link_libraries(${PROJECT_NAME} ${PYTHON_LIBRARIES})
  target_sources(${PROJECT_NAME} PRIVATE pkdpython.c)
endif(PYTHONLIBS_FOUND)

#if USE_SIMD
#if USE_SIMD_FMM
target_sources(${PROJECT_NAME} PRIVATE vmoments.cxx)
#endif
#if USE_SIMD_LC
target_sources(${PROJECT_NAME} PRIVATE lightcone.cxx)
#endif
#endif




foreach(flag_var
        CMAKE_C_FLAGS CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE
        CMAKE_C_FLAGS_MINSIZEREL CMAKE_C_FLAGS_RELWITHDEBINFO
        CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
        CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
   if(${flag_var} MATCHES "-DNDEBUG")
      string(REPLACE "-DNDEBUG" "" ${flag_var} "${${flag_var}}")
   endif(${flag_var} MATCHES "-DNDEBUG")
endforeach(flag_var)

option(HYDRO "Compile with hydrodynamics support" ON)
if (HYDRO)
   message("-- Compiling with hydrodynamics support")
   target_sources(${PROJECT_NAME} PRIVATE  hydro/hydro.c hydro/density.c
                  hydro/gradient.c hydro/flux.c hydro/step.c)
endif(HYDRO)

option(USE_MFM "Use the MFM hydro scheme instead of MFV" ON)
if (USE_MFM)
   message("-- Compiling with MFM hydrodynamics")
   target_compile_options(${PROJECT_NAME} PRIVATE -DUSE_MFM)
else(USE_MFM)
   message("-- Compiling with MFV hydrodynamics")
endif(USE_MFM)

option(FORCE_2D "Flag for running 2D cases" OFF)
if (FORCE_2D)
   message("-- NOT RECOMMENDED: Ignoring fluxes in the z-direction and interactions among misaligned particles")
   target_compile_options(${PROJECT_NAME} PRIVATE -DFORCE_2D)
endif(FORCE_2D)

option(FORCE_1D "Flag for running 1D cases" OFF)
if (FORCE_1D)
   message("-- NOT RECOMMENDED: Ignoring fluxes in the z and y directions and interactions among misaligned particles")
   target_compile_options(${PROJECT_NAME} PRIVATE -DFORCE_1D)
endif(FORCE_1D)

option(LIMITER_BARTH "Use the Barth-Jespersen limiter" OFF)
if (LIMITER_BARTH)
   message("-- Using the Barth-Jespersen limiter")
   target_compile_options(${PROJECT_NAME} PRIVATE -DLIMITER_BARTH)
endif(LIMITER_BARTH)

option(LIMITER_CONDBARTH "Use the Conditioned Barth-Jespersen limiter" ON)
if (LIMITER_CONDBARTH)
   message("-- Using the Conditioned Barth-Jespersen limiter")
   target_compile_options(${PROJECT_NAME} PRIVATE -DLIMITER_CONDBARTH)
endif(LIMITER_CONDBARTH)

option(ENTROPY_SWITCH "Assume entropy conservation when a cold flow is detected" OFF)
if (ENTROPY_SWITCH)
   message("-- Assuming entropy conservation for smooth cold flows")
   target_compile_options(${PROJECT_NAME} PRIVATE -DENTROPY_SWITCH)
endif(ENTROPY_SWITCH)

option(MAKE_GLASS "Create a glass distribution using an approximated Lloyd's algorithm" OFF)
if (MAKE_GLASS)
   message("-- Compiled to generate GLASS")
   target_compile_options(${PROJECT_NAME} PRIVATE -DMAKE_GLASS)
endif(MAKE_GLASS)

option(REGULARIZE_MESH "Try to maintain a regular mesh - only valid for MFV!" OFF)
if (REGULARIZE_MESH )
   message("-- Mesh regularization active. Only valid with MFV!")
   target_compile_options(${PROJECT_NAME} PRIVATE -DREGULARIZE_MESH)
endif(REGULARIZE_MESH)

option(FIXED_NSMOOTH_STRICT "Maintain a constant number of neighbors for computing the smoothing length")
option(FIXED_NSMOOTH_RELAXED "Maintain a approximately constant number of neighbors for computing the smoothing length")
if (FIXED_NSMOOTH_STRICT)
   if (FIXED_NSMOOTH_RELAXED)
      message("-- CAUTION! Only one mode for computing the smoothing length may be active")
   else(FIXED_NSMOOTH_RELAXED)
      message("-- Mantaining a strictly constant number of neighbours. Recommended for simulations with high density contrast (e.g., cosmological simulations)")
      target_compile_options(${PROJECT_NAME} PRIVATE -DFIXED_NSMOOTH_STRICT)
   endif(FIXED_NSMOOTH_RELAXED)
else(FIXED_SMOOTH_STRICT)
   if (FIXED_NSMOOTH_RELAXED)
      message("-- Mainting an approximately constant number of neighbours. This is not done with the density, thus a high dNeighborsStd must be given (~ a few)")
      target_compile_options(${PROJECT_NAME} PRIVATE -DFIXED_NSMOOTH_RELAXED)
   else(FIXED_NSMOOTH_RELAXED)
      message("-- Computing the smoothing length assuming a constant number density within a sphere (as in SPH)")
   endif(FIXED_NSMOOTH_RELAXED)
endif(FIXED_NSMOOTH_STRICT)

option(COOLING "Activate the chemistry and the cooling for the gas particles. Also makes effective a internal energy floor" OFF)
if (COOLING)
   target_sources(${PROJECT_NAME} PRIVATE  cooling/cooling.c cooling/cooling_tables.c )
   target_compile_options(${PROJECT_NAME} PRIVATE -DCOOLING)
   target_compile_options(${PROJECT_NAME} PRIVATE -DHAVE_HELIUM -DHAVE_CARBON -DHAVE_NITROGEN -DHAVE_OXYGEN -DHAVE_NEON -DHAVE_MAGNESIUM -DHAVE_SILICON -DHAVE_IRON)
   message("-- Cooling is ACTIVE")
endif(COOLING)

option(GRACKLE "Activate the cooling using grackle" OFF)
if (GRACKLE)
   find_package(GRACKLE REQUIRED)
   target_include_directories(${PROJECT_NAME} PUBLIC ${GRACKLE_INCLUDE})
   target_link_libraries(${PROJECT_NAME} ${GRACKLE_LIBRARY})
   target_compile_options(${PROJECT_NAME} PRIVATE -DGRACKLE -DCONFIG_BFLOAT_8)
   target_compile_options(${PROJECT_NAME} PRIVATE -DHAVE_METALLICITY)
   target_sources(${PROJECT_NAME} PRIVATE  cooling_grackle/cooling_grackle.c )
   message("-- GRACKLE is ACTIVE")
endif(GRACKLE)

option(STAR_FORMATION "Activate the formation of star particles from cold, dense gas" OFF)
if (STAR_FORMATION)
   target_sources(${PROJECT_NAME} PRIVATE  starformation/starformation.c starformation/init.c)
   target_compile_options(${PROJECT_NAME} PRIVATE -DSTAR_FORMATION)
   message("-- Star formation is ACTIVE")
endif(STAR_FORMATION)

option(EEOS_POLYTROPE "Apply a polytropic effective EOS" OFF)
if (EEOS_POLYTROPE)
   target_compile_options(${PROJECT_NAME} PRIVATE -DEEOS_POLYTROPE)
   message("-- Polytropic effective EOS")
endif(EEOS_POLYTROPE)

option(EEOS_JEANS "Apply an effective EOS to enforce that the Jeans length is resolved" OFF)
if (EEOS_JEANS)
   target_compile_options(${PROJECT_NAME} PRIVATE -DEEOS_JEANS)
   message("-- Jeans effective EOS")
endif(EEOS_JEANS)


option(FEEDBACK "Allow for spawned star particles to explode as type II supernova" OFF)
if (FEEDBACK)
   target_sources(${PROJECT_NAME} PRIVATE  starformation/feedback.c starformation/init.c) 
   target_compile_options(${PROJECT_NAME} PRIVATE -DFEEDBACK)
   message("-- Feedback from SNII is ACTIVE")
endif(FEEDBACK)

option(STELLAR_EVOLUTION "Activate time-resolved mass return from star particles into their surrounding media" OFF)
if (STELLAR_EVOLUTION)
   target_sources(${PROJECT_NAME} PRIVATE  stellarevolution/stellarevolution.c) 
   target_compile_options(${PROJECT_NAME} PRIVATE -DSTELLAR_EVOLUTION)
   target_compile_options(${PROJECT_NAME} PRIVATE -DHAVE_HELIUM -DHAVE_CARBON -DHAVE_NITROGEN -DHAVE_OXYGEN -DHAVE_NEON -DHAVE_MAGNESIUM -DHAVE_SILICON -DHAVE_IRON -DHAVE_METALLICITY)
   message("-- Stellar evolution is ACTIVE")
endif(STELLAR_EVOLUTION)

option(HERNQUIST_POTENTIAL "Add the contribution of DM using a Hernquist potential" OFF)
if (HERNQUIST_POTENTIAL)
  target_compile_options(${PROJECT_NAME} PRIVATE -DHERNQUIST_POTENTIAL)
  message("-- Analytical gravity ACTIVE -> Hernquist potential")
endif(HERNQUIST_POTENTIAL)

option(BLACKHOLES "Activate black hole model" OFF)
if (BLACKHOLES)
  target_compile_options(${PROJECT_NAME} PRIVATE -DBLACKHOLES)
  target_sources(${PROJECT_NAME} PRIVATE  blackhole/merger.c)
  target_sources(${PROJECT_NAME} PRIVATE  blackhole/seed.c)
  target_sources(${PROJECT_NAME} PRIVATE  blackhole/evolve.c)
  target_sources(${PROJECT_NAME} PRIVATE  blackhole/init.c)
  message("-- Black holes are ACTIVE")
endif(BLACKHOLES)
option(DEBUG_BH_NODRIFT "Deactivate BH drift by the gas" OFF)
if (DEBUG_BH_NODRIFT)
   target_compile_options(${PROJECT_NAME} PRIVATE -DDEBUG_BH_NODRIFT)
   message("-- Deactivated BH drift by the gas")
endif(DEBUG_BH_NODRIFT)
option(DEBUG_BH_ONLY "Allow for BH-only simulations for debugging purposes" OFF)
if (DEBUG_BH_ONLY)
   target_compile_options(${PROJECT_NAME} PRIVATE -DDEBUG_BH_ONLY)
   message("-- Allow for BH-only simulations for debugging purposes")
endif(DEBUG_BH_ONLY)

#### OPTIMIZATION OPTIONS
option(OPTIM_NO_REDUNDANT_FLUXES "Avoid any redundant flux computation using the particle IDs" ON)
if (OPTIM_NO_REDUNDANT_FLUXES)
   target_compile_options(${PROJECT_NAME} PRIVATE -DOPTIM_NO_REDUNDANT_FLUXES)
   message("-- OPTIM_NO_REDUNDANT_FLUXES: \t Avoid all redundant flux computations. Requires unique particleIDs and OPTIM_SMOOTH_NODE.")
endif(OPTIM_NO_REDUNDANT_FLUXES)

option(OPTIM_REMOVE_UNUSED "Remove unused variables in SPHFIELDS" ON)
if(OPTIM_REMOVE_UNUSED)
   target_compile_options(${PROJECT_NAME} PRIVATE -DOPTIM_REMOVE_UNUSED)
   message("-- OPTIM_REMOVE_UNUSED: \t Removed unused variables in particle structs")
endif(OPTIM_REMOVE_UNUSED)

option(OPTIM_REDUCE_PRECISION "Reduce precision of certain variable types to save memory" OFF)
if (OPTIM_REDUCE_PRECISION)
   target_compile_options(${PROJECT_NAME} PRIVATE -DOPTIM_REDUCE_PRECISION)
   message("-- OPTIM_REDUCE_PRECISION: \t Reduced precision to decrease particle size")
endif(OPTIM_REDUCE_PRECISION)

option(OPTIM_DENSITY_REITER "Reduce the number of density iterations when we have neighbors than expected" OFF)
if (OPTIM_DENSITY_REITER)
   target_compile_options(${PROJECT_NAME} PRIVATE -DOPTIM_DENSITY_REITER)
   message("-- OPTIM_DENSITY_REITER: \t Reduced number of density iterations when we have a larger fBall than expected")
endif(OPTIM_DENSITY_REITER)

option(OPTIM_INVERSE_WALK "When searching neighbours, start walking the tree from the leaves" OFF)
if (OPTIM_INVERSE_WALK)
   target_compile_options(${PROJECT_NAME} PRIVATE -DOPTIM_INVERSE_WALK)
   message("-- OPTIM_INVERSE_WALK: \t When searching neighbours, start walking the tree from the leaves (*NOTE*: broken for parallel runs)")
endif (OPTIM_INVERSE_WALK)

option(OPTIM_SMOOTH_NODE "Do the smooth operator for a whole node at once, rather than particle by particle" ON)
if (OPTIM_SMOOTH_NODE)
   target_compile_options(${PROJECT_NAME} PRIVATE -DOPTIM_SMOOTH_NODE)
   message("-- OPTIM_SMOOTH_NODE: \t\t Do the smooth operator for a whole node at once, rather than particle by particle. This overrrides other density loop optimizations!")
endif (OPTIM_SMOOTH_NODE)

option(OPTIM_REORDER_IN_NODES "Reorder the particle within the each cell memory space to decrease branching in hydrodynamic loops" ON)
if (OPTIM_REORDER_IN_NODES)
   target_compile_options(${PROJECT_NAME} PRIVATE -DOPTIM_REORDER_IN_NODES)
   message("-- OPTIM_REORDER_IN_NODES: \t Reorder the particle within the each cell memory space to decrease branching in hydrodynamic loops (only useful if OPTIM_SMOOTH_NODE is set).")
endif (OPTIM_REORDER_IN_NODES)

option(OPTIM_FLUX_VEC "Reorganize the flux interactions such that they can be vectorized by the compiler" ON)
if (OPTIM_FLUX_VEC)
   target_compile_options(${PROJECT_NAME} PRIVATE -DOPTIM_FLUX_VEC)
   message("-- OPTIM_FLUX_VEC: \t\t Reorganize the flux interactions such that they can be vectorized by the compiler (requires OPTIM_SMOOTH_NODE=on, USE_MFM=on)")
endif (OPTIM_FLUX_VEC)

option(OPTIM_UNION_EXTRAFIELDS "Use a union of structures to decrease memory footprint in simulations with different particle types" ON)
option(DEBUG_UNION_EXTRAFIELDS "Add extra asserts when accessing SPH,STAR or BH FIELDS for debugging" OFF)
if (OPTIM_UNION_EXTRAFIELDS)
   target_compile_options(${PROJECT_NAME} PRIVATE -DOPTIM_UNION_EXTRAFIELDS)
   message("-- OPTIM_UNION_EXTRAFIELDS: \t Use a union of structures to decrease memory footprint in simulations with different particle types. Requires OPTIM_SMOOTH_NODE.")
   if (DEBUG_UNION_EXTRAFIELDS)
      target_compile_options(${PROJECT_NAME} PRIVATE -DDEBUG_UNION_EXTRAFIELDS)
      message("-- DEBUG_UNION_EXTRAFIELDS: \t Added asserts when accessing EXTRAFIELDS. Disable for performance!!")
   endif (DEBUG_UNION_EXTRAFIELDS)
endif (OPTIM_UNION_EXTRAFIELDS)

option(DEBUG_FLUX_INFO "Show cache debug information for the flux computation" OFF)
if (DEBUG_FLUX_INFO)
   target_compile_options(${PROJECT_NAME} PRIVATE -DDEBUG_FLUX_INFO)
   message("-- DEBUG_FLUX_INFO: \t\t Show cache debug information for the flux computation")
endif (DEBUG_FLUX_INFO)

option(OPTIM_AVOID_IS_ACTIVE "Select particles rather than repeatedly using pkdIsActive" ON)
if (OPTIM_AVOID_IS_ACTIVE)
   target_compile_options(${PROJECT_NAME} PRIVATE -DOPTIM_AVOID_IS_ACTIVE)
   message("-- OPTIM_AVOID_IS_ACTIVE: \t Select particles rather than repeatedly using pkdIsActive")
endif (OPTIM_AVOID_IS_ACTIVE)

option(OPTIM_EXTRA "Activate extra, small, optimizations" OFF)
if (OPTIM_EXTRA)
   target_compile_options(${PROJECT_NAME} PRIVATE -DOPTIM_EXTRA)
   message("-- OPTIM_EXTRA: \t\t Activate extra, small, optimizations (WRONG - do not use!!)")
endif (OPTIM_EXTRA)
####


CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/pkd_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/pkd_config.h)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(${PROJECT_NAME} mdl2 openpa)
target_include_directories(${PROJECT_NAME} PRIVATE ${GSL_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${GSL_LIBRARIES})
target_link_libraries(${PROJECT_NAME} blitz)

target_include_directories(tostd PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(psout PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(psout PRIVATE ${GSL_INCLUDE_DIRS})
target_link_libraries(psout ${GSL_LIBRARIES})

install(TARGETS ${PROJECT_NAME} tostd DESTINATION "bin")
