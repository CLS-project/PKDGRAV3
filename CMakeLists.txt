cmake_minimum_required(VERSION 3.14)
cmake_policy(VERSION 3.14...3.22)
include(CheckLanguage)
include(CheckFunctionExists)
include(CheckCXXCompilerFlag)

# cmake uses -isysroot to specify the SDK, but this removes /usr/local/lib
# from consideration. We need this for the dependant libraries of Python3,
# most specifically libintl.
if (APPLE)
  link_directories(/usr/local/lib)
endif(APPLE)

project(pkdgrav3 VERSION 3.1.2)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/mdl2)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Release build with debug info selected")
  set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()


set(TARGET_ARCHITECTURE "auto" CACHE STRING "CPU architecture to optimize for.")
if(TARGET_ARCHITECTURE STREQUAL "auto")
  CHECK_CXX_COMPILER_FLAG("-march=native" COMPILER_OPT_ARCH_NATIVE_SUPPORTED)
  if (COMPILER_OPT_ARCH_NATIVE_SUPPORTED)
    add_compile_options(-march=native)
  endif()
  if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    CHECK_CXX_COMPILER_FLAG("/arch:AVX" COMPILER_OPT_ARCH_AVX_SUPPORTED)
    if (COMPILER_OPT_ARCH_AVX_SUPPORTED)
      add_compile_options(/arch:AVX)
    endif(COMPILER_OPT_ARCH_AVX_SUPPORTED)
  endif()
else()
    add_compile_options(-march=${TARGET_ARCHITECTURE})
endif()
CHECK_CXX_COMPILER_FLAG("-Wall" COMPILER_OPT_WARN_ALL_SUPPORTED)
if (COMPILER_OPT_WARN_ALL_SUPPORTED)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
    CHECK_CXX_COMPILER_FLAG("-Wno-sign-compare" COMPILER_OPT_WARN_NO_SIGN_COMPARE_SUPPORTED)
    if(COMPILER_OPT_WARN_NO_SIGN_COMPARE_SUPPORTED)
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wno-sign-compare")
        set(CMAKE_C_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wno-sign-compare")
    endif(COMPILER_OPT_WARN_NO_SIGN_COMPARE_SUPPORTED)
endif(COMPILER_OPT_WARN_ALL_SUPPORTED)
find_package(Git)
if(Git_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --git-dir
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_DIRECTORY
    ERROR_VARIABLE GIT_ERROR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if(NOT "${GIT_DIRECTORY}" STREQUAL "")
    string(REGEX REPLACE "[^/]+" ".." HOOK_RELATIVE_PATH "${GIT_DIRECTORY}")
    if(NOT EXISTS "${CMAKE_SOURCE_DIR}/${GIT_DIRECTORY}/hooks/pre-commit")
      ADD_CUSTOM_TARGET(link_target ALL COMMAND ${CMAKE_COMMAND} -E create_symlink
        "../${HOOK_RELATIVE_PATH}/git-pre-commit-hook"
        "${CMAKE_SOURCE_DIR}/${GIT_DIRECTORY}/hooks/pre-commit" )
    endif()
  endif()
endif()
find_package(GSL REQUIRED)      # GNU Scientific Library
find_package(HDF5 COMPONENTS C HL)
find_package(FFTW REQUIRED)
find_package(Python3 REQUIRED COMPONENTS Interpreter Development NumPy)
set(PYTHON ${Python3_EXECUTABLE})
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

add_subdirectory(fmt)
add_subdirectory(blitz)
add_subdirectory(mdl2)
add_subdirectory(cuda)
add_subdirectory(metal)
add_subdirectory(tests)

add_custom_command(
    OUTPUT ${CMAKE_BINARY_DIR}/parse.pyx
    MAIN_DEPENDENCY ${CMAKE_SOURCE_DIR}/parameters.toml
    DEPENDS ${CMAKE_SOURCE_DIR}/make_parse.py
    COMMAND ${PYTHON} ${CMAKE_SOURCE_DIR}/make_parse.py ${CMAKE_SOURCE_DIR}/parameters.toml ${CMAKE_BINARY_DIR}/parse.pyx
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    VERBATIM
)
add_custom_command(
    OUTPUT ${CMAKE_BINARY_DIR}/parse.cpp
    MAIN_DEPENDENCY ${CMAKE_BINARY_DIR}/parse.pyx
    COMMAND ${PYTHON} -m cython --3str --cplus --no-docstrings
        -o ${CMAKE_BINARY_DIR}/parse.cpp
        ${CMAKE_BINARY_DIR}/parse.pyx
    VERBATIM
)

add_custom_command(
    OUTPUT ${CMAKE_BINARY_DIR}/pkd_parameters.h ${CMAKE_BINARY_DIR}/pkd_parameters.pxi
    MAIN_DEPENDENCY ${CMAKE_SOURCE_DIR}/parameters.toml
    DEPENDS ${CMAKE_SOURCE_DIR}/make_parameters.py
    COMMAND ${PYTHON} ${CMAKE_SOURCE_DIR}/make_parameters.py ${CMAKE_SOURCE_DIR}/parameters.toml ${CMAKE_BINARY_DIR}/pkd_parameters.h
    ${CMAKE_BINARY_DIR}/pkd_parameters.pxi
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    VERBATIM
)

add_executable(${PROJECT_NAME} "")
target_compile_features(${PROJECT_NAME} PRIVATE cxx_auto_type cxx_range_for cxx_lambdas cxx_strong_enums)
set_target_properties(${PROJECT_NAME} PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)
set_property(TARGET ${PROJECT_NAME} APPEND PROPERTY COMPILE_DEFINITIONS _LARGEFILE_SOURCE)
target_sources(${PROJECT_NAME} PRIVATE
	main.cxx cosmo.c master.cxx simulate.cxx pst.cxx TraversePST.cxx io/fio.c core/illinois.c param.c pyrameters.cxx ${CMAKE_BINARY_DIR}/pkd_parameters.h
	pkd.cxx analysis/analysis.cxx smooth/smooth.cxx smooth/smoothfcn.cxx io/outtype.cxx io/output.cxx io/service.cxx
	gravity/walk2.cxx gravity/grav2.cxx gravity/ewald.cxx ic/ic.cxx domains/tree.cxx gravity/opening.cxx gravity/pp.cxx gravity/pc.cxx gravity/cl.cxx
	gravity/lst.cxx gravity/moments.c gravity/ilp.cxx gravity/ilc.cxx io/iomodule.cxx io/iochunk.cxx io/restore.cxx
	group/fof.cxx group/hop.cxx group/group.cxx group/groupstats.cxx ic/RngStream.c smooth/listcomp.c core/healpix.c
	core/gridinfo.cxx analysis/interlace.cxx analysis/contrast.cxx analysis/assignmass.cxx analysis/measurepk.cxx bispectrum.cxx ic/whitenoise.cxx gravity/pmforces.cxx
	core/setadd.cxx core/hostname.cxx core/initcosmology.cxx core/calcroot.cxx core/swapall.cxx core/select.cxx core/particle.cxx
	domains/calcbound.cxx domains/combinebound.cxx domains/distribtoptree.cxx domains/distribroot.cxx domains/dumptrees.cxx
	domains/enforceperiodic.cxx domains/freestore.cxx domains/olddd.cxx domains/getordsplits.cxx
	gravity/setsoft.cxx gravity/activerung.cxx gravity/countrungs.cxx gravity/zeronewrung.cxx
	analysis/rsloadids.cxx analysis/rssaveids.cxx analysis/rsextract.cxx analysis/rsreorder.cxx
	eEOS/eEOS.cxx
	SPH/SPHOptions.cxx SPH/SPHEOS.cxx SPH/SPHpredict.cxx initlightcone.cxx
)
set_property(SOURCE io/fio.c APPEND PROPERTY COMPILE_DEFINITIONS "USE_PTHREAD")
#add_executable(tostd utility/tostd.c io/fio.c)
#add_executable(psout utility/psout.c cosmo.c)

target_link_libraries(${PROJECT_NAME} m)
#target_link_libraries(tostd m)
#target_link_libraries(psout m)

make_directory(${CMAKE_BINARY_DIR}/modules)
add_custom_command(
    OUTPUT ${CMAKE_BINARY_DIR}/modules/accuracy.cpp
    MAIN_DEPENDENCY ${CMAKE_SOURCE_DIR}/modules/accuracy.py
    COMMAND ${PYTHON} -m cython --3str --cplus --no-docstrings
        -o ${CMAKE_BINARY_DIR}/modules/accuracy.cpp
        ${CMAKE_SOURCE_DIR}/modules/accuracy.py
    VERBATIM
)
target_sources(${PROJECT_NAME} PRIVATE ${CMAKE_BINARY_DIR}/modules/accuracy.cpp)
add_custom_command(
    OUTPUT ${CMAKE_BINARY_DIR}/modules/checkpoint.cpp ${CMAKE_BINARY_DIR}/modules/checkpoint.h
    MAIN_DEPENDENCY ${CMAKE_SOURCE_DIR}/modules/checkpoint.pyx
    COMMAND ${PYTHON} -m cython --3str --cplus --no-docstrings
        -o ${CMAKE_BINARY_DIR}/modules/checkpoint.cpp
        ${CMAKE_SOURCE_DIR}/modules/checkpoint.pyx
    VERBATIM
)
target_sources(${PROJECT_NAME} PRIVATE ${CMAKE_BINARY_DIR}/modules/checkpoint.cpp)
add_custom_command(
    OUTPUT ${CMAKE_BINARY_DIR}/modules/PKDGRAV.cpp ${CMAKE_BINARY_DIR}/modules/PKDGRAV.h
    MAIN_DEPENDENCY ${CMAKE_SOURCE_DIR}/modules/PKDGRAV.py
    DEPENDS ${CMAKE_SOURCE_DIR}/modules/PKDGRAV.pxd ${CMAKE_BINARY_DIR}/pkd_parameters.pxi
    COMMAND ${PYTHON} -m cython --3str --cplus --no-docstrings
        -I ${CMAKE_BINARY_DIR}
        -o ${CMAKE_BINARY_DIR}/modules/PKDGRAV.cpp
        ${CMAKE_SOURCE_DIR}/modules/PKDGRAV.py
    VERBATIM
)
target_sources(${PROJECT_NAME} PRIVATE ${CMAKE_BINARY_DIR}/modules/PKDGRAV.cpp)
add_custom_command(
    OUTPUT ${CMAKE_BINARY_DIR}/modules/cosmology.cpp
    MAIN_DEPENDENCY ${CMAKE_SOURCE_DIR}/modules/cosmology.py
    DEPENDS ${CMAKE_SOURCE_DIR}/modules/cosmology.pxd
    COMMAND ${PYTHON} -m cython --3str --cplus --no-docstrings
        -o ${CMAKE_BINARY_DIR}/modules/cosmology.cpp
        ${CMAKE_SOURCE_DIR}/modules/cosmology.py
    VERBATIM
)
target_sources(${PROJECT_NAME} PRIVATE ${CMAKE_BINARY_DIR}/modules/cosmology.cpp)

add_library(cosmology MODULE ${CMAKE_CURRENT_BINARY_DIR}/modules/cosmology.cpp ${CMAKE_CURRENT_SOURCE_DIR}/cosmo.c)
target_include_directories(cosmology PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${GSL_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})
target_link_libraries(cosmology Python3::Python Python3::NumPy ${GSL_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} mdl2)
set_target_properties(cosmology PROPERTIES PREFIX "")  

if (NOT APPLE)
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
endif()
if (HDF5_FOUND)
  set(USE_HDF5 TRUE)
  target_include_directories(${PROJECT_NAME} PUBLIC ${HDF5_INCLUDE_DIRS})
  target_link_libraries(${PROJECT_NAME} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
  #set_property(TARGET ${PROJECT_NAME} APPEND PROPERTY COMPILE_DEFINITIONS ${HDF5_DEFINITIONS})
  #set_property(TARGET ${PROJECT_NAME} APPEND PROPERTY COMPILE_DEFINITIONS H5_USE_16_API)
  #target_include_directories(tostd PUBLIC ${HDF5_INCLUDE_DIRS})
  #target_link_libraries(tostd ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
  #target_include_directories(psout PUBLIC ${HDF5_INCLUDE_DIRS})
  #target_link_libraries(psout ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})
endif(HDF5_FOUND)
if (Python3_Development_FOUND)
  set(USE_PYTHON TRUE)
  target_include_directories(${PROJECT_NAME} PUBLIC ${Python3_INCLUDE_DIRS})
  target_link_libraries(${PROJECT_NAME} Python3::Python Python3::NumPy)
  target_sources(${PROJECT_NAME} PRIVATE msrpython.cxx ${CMAKE_BINARY_DIR}/parse.cpp csmpython.c)
  if (Python3_NumPy_FOUND)
    target_include_directories(${PROJECT_NAME} PUBLIC ${Python3_NumPy_INCLUDE_DIRS})
  endif()
else()
  target_sources(${PROJECT_NAME} PRIVATE tinypy.c pkdtinypy.c)
endif(Python3_Development_FOUND)

set(USE_SIMD "yes" CACHE STRING "Should SIMD be used.")
if(USE_SIMD)
#if USE_SIMD_FMM
target_sources(${PROJECT_NAME} PRIVATE gravity/vmoments.cxx)
#endif
target_sources(${PROJECT_NAME} PRIVATE lightcone.cxx)
endif()

foreach(flag_var
        CMAKE_C_FLAGS CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE
        CMAKE_C_FLAGS_MINSIZEREL CMAKE_C_FLAGS_RELWITHDEBINFO
        CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
        CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
   if(${flag_var} MATCHES "-DNDEBUG")
      string(REPLACE "-DNDEBUG" "" ${flag_var} "${${flag_var}}")
   endif(${flag_var} MATCHES "-DNDEBUG")
endforeach(flag_var)

# The header files for rockstar conflict with C++ so make symbolic links
set(ROCKSTAR_PATH CACHE STRING "Path to the Rockstar source code")
if(ROCKSTAR_PATH)
    set(HAVE_ROCKSTAR TRUE)
    file(REAL_PATH ${ROCKSTAR_PATH} rdir EXPAND_TILDE)
    make_directory(${CMAKE_BINARY_DIR}/rockstar)
    make_directory(${CMAKE_BINARY_DIR}/rockstar/io)
    foreach(file particle.h halo.h io/io_internal.h)
        #find_file(rs_src ${file} PATHS ${rdir} NO_DEFAULT_PATH NO_CACHE REQUIRED)
        add_custom_command(
            OUTPUT ${CMAKE_BINARY_DIR}/rockstar/${file}
            MAIN_DEPENDENCY ${rdir}/${file}
            COMMAND ${CMAKE_COMMAND} -E create_symlink
                ${rdir}/${file}
                ${CMAKE_BINARY_DIR}/rockstar/${file}
            VERBATIM
        )
        target_sources(${PROJECT_NAME} PRIVATE ${CMAKE_BINARY_DIR}/rockstar/${file})
    endforeach()
    target_sources(${PROJECT_NAME} PRIVATE analysis/rshaloloadids.cxx analysis/rshalocount.cxx)
endif()

option(HYDRO "Compile with hydrodynamics support" ON)
if (HYDRO)
   message("-- Compiling with hydrodynamics support")
   target_sources(${PROJECT_NAME} PRIVATE  hydro/hydro.cxx hydro/density.cxx
                  hydro/gradient.cxx hydro/flux.cxx hydro/step.cxx)
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
   target_sources(${PROJECT_NAME} PRIVATE  cooling/cooling.cxx cooling/cooling_tables.c )
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
   target_sources(${PROJECT_NAME} PRIVATE  cooling_grackle/cooling_grackle.cxx )
   message("-- GRACKLE is ACTIVE")
endif(GRACKLE)

option(STAR_FORMATION "Activate the formation of star particles from cold, dense gas" OFF)
if (STAR_FORMATION)
   target_sources(${PROJECT_NAME} PRIVATE  starformation/starformation.cxx starformation/init.cxx)
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
   target_sources(${PROJECT_NAME} PRIVATE  starformation/feedback.cxx starformation/init.cxx)
   target_compile_options(${PROJECT_NAME} PRIVATE -DFEEDBACK)
   message("-- Feedback from SNII is ACTIVE")
endif(FEEDBACK)

option(STELLAR_EVOLUTION "Activate time-resolved mass return from star particles into their surrounding media" OFF)
if (STELLAR_EVOLUTION)
   target_sources(${PROJECT_NAME} PRIVATE  stellarevolution/stellarevolution.cxx)
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
  target_sources(${PROJECT_NAME} PRIVATE  blackhole/merger.cxx)
  target_sources(${PROJECT_NAME} PRIVATE  blackhole/seed.cxx)
  target_sources(${PROJECT_NAME} PRIVATE  blackhole/evolve.cxx)
  target_sources(${PROJECT_NAME} PRIVATE  blackhole/init.cxx)
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

option(OPTIM_DENSITY_REITER "Reduce the number of density iterations when we have neighbors than expected" OFF)
if (OPTIM_DENSITY_REITER)
   target_compile_options(${PROJECT_NAME} PRIVATE -DOPTIM_DENSITY_REITER)
   message("-- OPTIM_DENSITY_REITER: \t Reduced number of density iterations when we have a larger fBall than expected")
endif(OPTIM_DENSITY_REITER)

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

option(DEBUG_FLUX_NOLIMITER "Avoid the use of the limiters to achieve second-order accuracy at the cost of worse stability" OFF)
if (DEBUG_FLUX_NOLIMITER)
   target_compile_options(${PROJECT_NAME} PRIVATE -DDEBUG_FLUX_NOLIMITER)
   message("-- DEBUG_FLUX_NOLIMITER: \t\t Avoid the use of limiters, increasing accuracy at the cost of stability")
endif (DEBUG_FLUX_NOLIMITER)

option(OPTIM_AVOID_IS_ACTIVE "Select particles rather than repeatedly using pkdIsActive" ON)
if (OPTIM_AVOID_IS_ACTIVE)
   target_compile_options(${PROJECT_NAME} PRIVATE -DOPTIM_AVOID_IS_ACTIVE)
   message("-- OPTIM_AVOID_IS_ACTIVE: \t Select particles rather than repeatedly using pkdIsActive")
endif (OPTIM_AVOID_IS_ACTIVE)
####

find_path(EOSLIB_PATH EOSlib.h HINTS "${CMAKE_SOURCE_DIR}/../EOSlib" NO_DEFAULT_PATH)
find_path(ANEOSMATERIAL_PATH ANEOSmaterial.h HINTS "${CMAKE_SOURCE_DIR}/../ANEOSmaterial" NO_DEFAULT_PATH)
find_path(TILLOTSON_PATH tillotson.h HINTS "${CMAKE_SOURCE_DIR}/../tillotson" NO_DEFAULT_PATH)
find_path(REOS3_PATH reos3.h HINTS "${CMAKE_SOURCE_DIR}/../reos3" NO_DEFAULT_PATH)
find_path(SCVHEOS_PATH scvheos.h HINTS "${CMAKE_SOURCE_DIR}/../scvh" NO_DEFAULT_PATH)
if(EOSLIB_PATH)
    message(STATUS "found EOSlib.h in ${EOSLIB_PATH}")
    target_sources(${PROJECT_NAME} PRIVATE ${EOSLIB_PATH}/EOSlib.c ${EOSLIB_PATH}/igeos.c)
    target_include_directories(${PROJECT_NAME} PUBLIC ${EOSLIB_PATH})
    target_compile_options(${PROJECT_NAME} PRIVATE -DHAVE_EOSLIB_H)
endif()
if(EOSLIB_PATH AND ANEOSMATERIAL_PATH)
    message(STATUS "found ANEOSmaterial.h in ${ANEOSMATERIAL_PATH}")
    target_sources(${PROJECT_NAME} PRIVATE ${ANEOSMATERIAL_PATH}/ANEOSmaterial.c ${ANEOSMATERIAL_PATH}/interpBilinear.c)
    target_include_directories(${PROJECT_NAME} PUBLIC ${ANEOSMATERIAL_PATH})
    target_compile_options(${PROJECT_NAME} PRIVATE -DHAVE_ANEOSMATERIAL_H)
endif()
if(EOSLIB_PATH AND TILLOTSON_PATH)
    message(STATUS "found tillotson.h in ${TILLOTSON_PATH}")
    target_sources(${PROJECT_NAME} PRIVATE ${TILLOTSON_PATH}/tillotson.c ${TILLOTSON_PATH}/tillinitlookup.c ${TILLOTSON_PATH}/tillsplint.c ${TILLOTSON_PATH}/interpol/brent.c)
    target_include_directories(${PROJECT_NAME} PUBLIC ${TILLOTSON_PATH})
    target_compile_options(${PROJECT_NAME} PRIVATE -DHAVE_TILLOTSON_H)
endif()
if(EOSLIB_PATH AND REOS3_PATH)
    message(STATUS "found reos3.h in ${REOS3_PATH}")
    target_sources(${PROJECT_NAME} PRIVATE ${REOS3_PATH}/reos3.c)
    target_include_directories(${PROJECT_NAME} PUBLIC ${REOS3_PATH})
    target_compile_options(${PROJECT_NAME} PRIVATE -DHAVE_REOS3_H)
endif()
if(EOSLIB_PATH AND SCVHEOS_PATH)
    message(STATUS "found scvheos.h in ${SCVHEOS_PATH}")
    target_sources(${PROJECT_NAME} PRIVATE ${SCVHEOS_PATH}/scvheos.c)
    target_include_directories(${PROJECT_NAME} PUBLIC ${SCVHEOS_PATH})
    target_compile_options(${PROJECT_NAME} PRIVATE -DHAVE_SCVHEOS_H)
endif()

option(NN_FLAG_IN_PARTICLE "Add a NN flag to the particle." ON)
if (NN_FLAG_IN_PARTICLE)
   target_compile_options(${PROJECT_NAME} PRIVATE -DNN_FLAG_IN_PARTICLE)
   message("-- NN_FLAG_IN_PARTICLE: A flag is added to the base particle to allow FastGas with interface correction. This reduces the maximum number of particles from 2^43 to 2^42. Set this option to OFF if more than 2^42 particles are needed in non-NewSPH runs.")
endif (NN_FLAG_IN_PARTICLE)

option(BOUND_USES_MINMAX "Bounds use lower and upper coordinate if ON, otherwise center and apothem" OFF)
option(NEW_REORDER "Used the new MDL swap based particle reorder method" OFF)
if(NEW_REORDER)
   target_sources(${PROJECT_NAME} PRIVATE domains/reorder.cxx)
endif(NEW_REORDER)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/pkd_config.h.in ${CMAKE_CURRENT_BINARY_DIR}/pkd_config.h)
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(${PROJECT_NAME} mdl2 pkdcuda pkdmetal)
target_include_directories(${PROJECT_NAME} PRIVATE ${GSL_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} ${GSL_LIBRARIES})
target_link_libraries(${PROJECT_NAME} blitz fmt)

enable_testing()
add_test(NAME cosmology COMMAND ${PROJECT_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/tests/cosmology.py WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(NAME gravity   COMMAND ${PROJECT_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/tests/gravity.py   WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(NAME select    COMMAND ${PROJECT_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/tests/select.py    WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(NAME fof   COMMAND ${PROJECT_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/tests/foftest.py   WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
if(EOSLIB_PATH AND ANEOSMATERIAL_PATH)
add_test(NAME NewSPH   COMMAND ${PROJECT_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/tests/NewSPH.py   WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
endif()

#target_include_directories(tostd PUBLIC ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
#target_include_directories(psout PUBLIC ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
#target_include_directories(psout PRIVATE ${GSL_INCLUDE_DIRS})
#target_link_libraries(psout ${GSL_LIBRARIES})

install(TARGETS ${PROJECT_NAME} ${PROJECT_NAME} DESTINATION "bin")
#install(TARGETS ${PROJECT_NAME} tostd DESTINATION "bin")
