/*  This file is part of PKDGRAV3 (http://www.pkdgrav.org/).
 *  Copyright (c) 2001-2018 Joachim Stadel & Douglas Potter
 *
 *  PKDGRAV3 is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  PKDGRAV3 is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with PKDGRAV3.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "analysis/orbits.h"
#include "blitz/tinyvec2.h"
#include "core/partclass.h"
#include "io/service.h"
#include "pkd_config.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdio.h>
#include <stdlib.h>
#ifdef HAVE_UNISTD_H
    #include <unistd.h> /* for unlink() */
#endif
#ifdef _MSC_VER
    #define mkdir _mkdir
    #define unlink _unlink
#endif
#include <stddef.h>
#include <string.h>
#include <ctype.h>
#include <cinttypes>
#include <limits.h>
#include <stdarg.h>
#include <assert.h>
#include "core/abort_error.h"
#include <time.h>
#ifdef HAVE_SYS_TIME_H
    #include <sys/time.h>
#endif
#ifdef HAVE_SYS_STAT_H
    #include <sys/stat.h>
#endif
#include <math.h>
#if defined(HAVE_WORDEXP) && defined(HAVE_WORDFREE)
    #include <wordexp.h>
#elif defined(HAVE_GLOB) && defined(HAVE_GLOBFREE)
    #include <glob.h>
#endif
#include <sys/stat.h>
#include <algorithm>
#include <numeric>
#include <ranges>
#include <functional>
#include <fstream>
#include <filesystem>
#ifdef HAVE_SYS_PARAM_H
    #include <sys/param.h> /* for MAXHOSTNAMELEN, if available */
#endif

#include "master.h"
using namespace fmt::literals;        // Gives us ""_a and ""_format literals
using namespace pyrameters::literals; // Gives us ""_key literals

struct DeferFormat {
    std::string name;
};

// Formatter that reproduces "{name[:specs]}" literally.
template <>
struct fmt::formatter<DeferFormat> {
    // Capture any format spec after ':' up to the closing '}'
    template <typename ParseContext>
    constexpr auto parse(ParseContext& ctx) {
        auto it = ctx.begin();
        auto end = ctx.end();

        for (int depth = 0; it != end; ++it) {
            char c = *it;
            if (c == '{') ++depth;
            else if (c == '}') {
                if (depth == 0) break;
                --depth;
            }
        }
        specs_.assign(ctx.begin(), it); // whatever was written after ':'
        // Advance to the end (fmt expects us to return the iterator at '}' or end)
        return it;
    }

    template <typename FormatContext>
    auto format(const DeferFormat& d, FormatContext& ctx) const {
        auto out = ctx.out();
        // Re-emit the original placeholder, preserving any specs
        if (specs_.empty()) {
            return fmt::format_to(out, "{{{}}}", d.name);
        }
        else {
            return fmt::format_to(out, "{{{}:{}}}", d.name, specs_);
        }
    }

  private:
    std::string specs_;
};

#include "core/illinois.h"
#include "io/outtype.h"
#include "smooth/smoothfcn.h"
#include "io/fio.h"
#include "SPH/SPHOptions.h"

#include "services/smooth.h"
#include "services/reorderinnodes.h"
#include "services/chemcompinit.h"
#include "services/sendarray.h"
#include "services/checkpoint.h"
#include "services/calceandl.h"
#include "core/setadd.h"
#include "core/swapall.h"
#include "core/hostname.h"
#include "core/memstatus.h"
#include "core/ignoresigbus.h"
#include "core/initcosmology.h"
#include "core/calcroot.h"
#include "core/select.h"
#include "core/fftgrid.h"
#include "core/countspecies.h"
#include "core/removedeleted.h"
#include "core/decorator_reorder.h"
#include "io/restore.h"
#include "io/output.h"
#include "io/classes.h"
#include "io/ascii.h"
#include "ic/ic.h"
#include "analysis/lightcone.h"
#ifdef INSTRUMENT_MPI
    #include "core/cachestats.h"
#endif

#ifdef HAVE_ROCKSTAR
    #include "analysis/rshalocount.h"
    #include "analysis/rshaloloadids.h"
namespace rockstar {
    #include "rockstar/io/io_internal.h"
    #include "rockstar/halo.h"
} // namespace rockstar
#endif

#include "analysis/analysis.h"
#include "analysis/measurepk.h"
#include "analysis/assignmass.h"
#include "analysis/bispectrum.h"
#include "analysis/rsloadids.h"
#include "analysis/rsreorder.h"
#include "analysis/rssaveids.h"
#include "analysis/rsextract.h"
#include "group/fof.h"
#include "group/hop.h"
#include "group/groupstats.h"

#include "domains/buildtree.h"
#include "domains/calcbound.h"
#include "domains/distribtoptree.h"
#include "domains/distribroot.h"
#include "domains/dumptrees.h"
#include "domains/getordsplits.h"
#include "domains/olddd.h"
#include "domains/reorder.h"
#include "domains/enforceperiodic.h"
#include "domains/combinebound.h"
#include "gravity/setsoft.h"
#include "gravity/activerung.h"
#include "gravity/countrungs.h"
#include "gravity/pmforces.h"
#include "gravity/updaterung.h"
#include "gravity/zeronewrung.h"
#include "gravity/timestep.h"

#include "SPH/gasvalues.h"
#include "SPH/resetcom.h"
#include "imf.h"
#include "meshless/hydro/hydro.h"
#include "meshless/hydro/endstep.h"
#include "meshless/hydro/wakeup.h"
#include "meshless/blackhole/init.h"
#include "meshless/blackhole/seed.h"
#include "meshless/blackhole/drift.h"
#include "meshless/blackhole/accretion.h"
#ifdef STELLAR_EVOLUTION
    #include "stellarevolution/stellarevolution.h"
#endif
#include "meshless/starformation/init.h"
#include "meshless/starformation/starformation.h"
#include "meshless/cooling/cooling.h"
#include "meshless/cooling/cooling_tables.h"
#include "eEOS/eEOS.h"

#define LOCKFILE ".lockfile"   /* for safety lock */
#define STOPFILE "STOP"        /* for user interrupt */
#define CHECKFILE "CHECKPOINT" /* for user interrupt */

// This function will take a buffer (e.g., char filename[256]) and safely copy
// the contents of a std::string_view into it, ensuring null-termination.
template <std::size_t N>
void marshal_string(char (&dst)[N], std::string_view src) {
    if (src.length() >= N) throw std::runtime_error("string is too long");
    memcpy(dst, src.data(), src.length());
    dst[src.length()] = '\0';
}

#ifdef _MSC_VER
double MSR::Time() {
    FILETIME ft;
    uint64_t clock;
    GetSystemTimeAsFileTime(&ft);
    clock = ft.dwHighDateTime;
    clock <<= 32;
    clock |= ft.dwLowDateTime;
    /* clock is in 100 nano-second units */
    return clock / 10000000.0;
}
#else
double MSR::Time() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (tv.tv_sec + (tv.tv_usec * 1e-6));
}
#endif

void MSR::TimerStart(int iTimer) {
    ti[iTimer].sec = Time();
}

void MSR::TimerStop(int iTimer) {
    ti[iTimer].sec = Time() - ti[iTimer].sec;
    ti[iTimer].acc += ti[iTimer].sec;
}

// Query the timing of the last call, given that TimerStop was called
// for that iTimer
double MSR::TimerGet(int iTimer) {
    return ti[iTimer].sec;
}

double MSR::TimerGetAcc(int iTimer) {
    return ti[iTimer].acc;
}

// The order should be the same than in the enumerate above!
static const char* timer_names[TOTAL_TIMERS] = {"Gravity",  "IO",   "Tree",     "DomainDecom", "KickOpen", "KickClose", "Density",  "EndTimeStep",
                                                "Gradient", "Flux", "TimeStep", "Drift",       "FoF",      "Feedback",  "StarForm", "BHs",
#ifdef STELLAR_EVOLUTION
                                                "Stev",
#endif
                                                "Others"};
static_assert(sizeof(timer_names) / sizeof(timer_names[0]) == TOTAL_TIMERS);

void MSR::TimerHeader() {
    auto achFile = MakePath(OutName(), ".timing");
    auto fpLog = fopen(achFile.c_str(), "a");
    fprintf(fpLog, "# Step");
    for (int iTimer = 0; iTimer < TOTAL_TIMERS; iTimer++) {
        fprintf(fpLog, " %s", timer_names[iTimer]);
    }
    fprintf(fpLog, "\n");
    fclose(fpLog);
}

void MSR::TimerDump(int iStep) {
    auto achFile = MakePath(OutName(), ".timing");
    auto fpLog = fopen(achFile.c_str(), "a");

    fprintf(fpLog, "%d", iStep);
    for (int iTimer = 0; iTimer < TOTAL_TIMERS; iTimer++) {
        fprintf(fpLog, " %f", TimerGetAcc(iTimer));
    }
    fprintf(fpLog, "\n");
    fclose(fpLog);
}

void MSR::TimerRestart() {
    for (int iTimer = 0; iTimer < TOTAL_TIMERS; iTimer++) {
        ti[iTimer].acc = 0.0;
        ti[iTimer].sec = 0.0;
    }
}

void MSR::Leader(void) {
    puts("pkdgrav" PACKAGE_VERSION " Joachim Stadel & Doug Potter Sept 2015");
    puts("USAGE: pkdgrav3 [SETTINGS | FLAGS] [SIM_FILE]");
    puts("SIM_FILE: Configuration file of a particular simulation, which");
    puts("          includes desired settings and relevant input and");
    puts("          output files. Settings specified in this file override");
    puts("          the default settings.");
    puts("SETTINGS");
    puts("or FLAGS: Command line settings or flags for a simulation which");
    puts("          will override any defaults and any settings or flags");
    puts("          specified in the SIM_FILE.");
}
void MSR::Trailer(void) {
    puts("(see man page for more information)");
}

void MSR::Exit(int status) {
    exit(status);
}

static void make_directories(std::string name) {
    auto i = name.rfind('/');
    if (i > 0) {
        if (i != std::string::npos) {
            name = name.substr(0, i);
            make_directories(name);
            mkdir(name.c_str(), 0755);
        }
    }
}

std::string MSR::BuildName(std::string_view path, int iStep, std::string_view type) {
    if (type.data() == nullptr) type = "";
    std::string_view kind = type.length() ? type.substr(1) : "output";
    if (!path.length()) path = "{name}.{step:05d}{type}";
    auto name = fmt::format(path, "name"_a = OutName(), "step"_a = iStep, "type"_a = type, "kind"_a = kind,
                            "octant"_a = DeferFormat{"octant"}, // Defers the evaluation of octant (if present)
                            "i"_a = DeferFormat{"i"});          // Defers the evaluation of i (if present)
    // If the octant field is present then we may need to create multiple directories
    if (fmt::format(name, "octant"_a = 0, "i"_a = DeferFormat{"i"}) != name) {
        int mask = 0;
        if (parameters.has_lstLightConeOctants()) {
            auto octants = parameters.get_lstLightConeOctants();
            for (auto o : octants) mask |= (1 << o);
        }
        else mask = 1; // weird, write to "octant 0"
        for (auto octant = 0; octant < 8; ++octant) {
            if (mask & (1 << octant)) {
                make_directories(fmt::format(name, "octant"_a = octant, "i"_a = DeferFormat{"i"}));
            }
        }
    }
    else {
        make_directories(name);
    }
    return name;
}
std::string MSR::BuildNameWithI(std::string_view path, int iStep, std::string_view type) {
    auto name = BuildName(path, iStep, type);
    // If the string does not contain an 'i' field, then we append ".{i}"
    if (fmt::format(name, "i"_a = 0) == name) name += ".{i}";
    return name;
}

std::string MSR::BuildName(int iStep, std::string_view type) {
    std::string_view kind = type.length() ? type.substr(1) : "output";
    return BuildName(parameters.get_achOutPath("type"_key = type, "kind"_key = kind), iStep, type);
}
std::string MSR::BuildNameWithI(int iStep, std::string_view type) {
    auto name = BuildName(iStep, type);
    // If the string does not contain an 'i' field, then we append ".{i}"
    if (fmt::format(name, "i"_a = 0) == name) name += ".{i}";
    return name;
}
std::string MSR::BuildNameWithOI(int iStep, std::string_view type) {
    auto name = BuildName(iStep, type);
    // If the string does not contain an 'octant' field, then we append ".{octant}"
    if (fmt::format(name, "octant"_a = 0, "i"_a = DeferFormat{"i"}) == name) name += ".{octant}";
    // If the string does not contain an 'i' field, then we append ".{i}"
    if (fmt::format(name, "i"_a = 0, "octant"_a = DeferFormat{"octant"}) == name) name += ".{i}";
    return name;
}

std::string MSR::BuildIoName(int iStep, std::string_view type) {
    std::string_view kind = type.length() ? type.substr(1) : "output";
    auto achIoPath = parameters.get_achIoPath("type"_key = type, "kind"_key = kind);
    if (achIoPath.length()) return BuildName(achIoPath, iStep, type);
    else return BuildName(iStep, type);
}
std::string MSR::BuildCpName(int iStep) {
    static constexpr std::string_view type = ".chk";
    static constexpr std::string_view kind = "chk";
    auto achCheckpointPath = parameters.get_achCheckpointPath("type"_key = type, "kind"_key = kind);
    if (achCheckpointPath.length() == 0) achCheckpointPath = parameters.get_achOutPath("type"_key = type, "kind"_key = kind);
    if (achCheckpointPath.length() == 0) achCheckpointPath = "{name}.{step:05d}{type}";
    auto name = fmt::format(achCheckpointPath, "name"_a = OutName(), "step"_a = iStep, "type"_a = type, "kind"_a = kind);
    make_directories(name);
    return name;
}

// prepends "dir" to "base" and returns the result.
std::string MSR::MakePath(std::string_view base, std::string_view ext) {
    auto dir = parameters.get_achDataSubPath();
    std::string path;
    path.reserve(dir.length() + 1 + base.length() + ext.length());
    if (dir.length()) {
        path = dir;
        path += '/';
    }
    path += base;
    path += ext;
    make_directories(path);
    return path;
}

uint64_t MSR::getMemoryModel() {
    uint64_t mMemoryModel = 0;
    /*
    ** Figure out what memory models are in effect.  Configuration flags
    ** can be used to request a specific model, but certain operations
    ** will force these flags to be on.
    */
    if (parameters.get_bFindGroups()) {
        mMemoryModel |= PKD_MODEL_GROUPS | PKD_MODEL_VELOCITY;
        if (parameters.get_bMemGlobalGid()) {
            mMemoryModel |= PKD_MODEL_GLOBALGID;
        }
    }
    if (DoGravity()) {
        mMemoryModel |= PKD_MODEL_VELOCITY | PKD_MODEL_NODE_MOMENT;
        if (!parameters.get_bNewKDK()) mMemoryModel |= PKD_MODEL_ACCELERATION;
    }
    if (parameters.get_bDoDensity()) mMemoryModel |= PKD_MODEL_DENSITY;
    if (parameters.get_bMemIntegerPosition()) mMemoryModel |= PKD_MODEL_INTEGER_POS;
    if (parameters.get_bMemUnordered() && parameters.get_bNewKDK()) mMemoryModel |= PKD_MODEL_UNORDERED;
    if (parameters.get_bMemParticleID()) mMemoryModel |= PKD_MODEL_PARTICLE_ID;
    if (parameters.get_bMemAcceleration() || parameters.get_bDoAccOutput()) mMemoryModel |= PKD_MODEL_ACCELERATION;
    if (parameters.get_bMemVelocity()) mMemoryModel |= PKD_MODEL_VELOCITY;
    if (parameters.get_bMemPotential() || parameters.get_bDoPotOutput()) mMemoryModel |= PKD_MODEL_POTENTIAL;
    if (parameters.get_bFindHopGroups()) mMemoryModel |= PKD_MODEL_GROUPS | PKD_MODEL_DENSITY | PKD_MODEL_BALL;
    if (parameters.get_bMemGroups()) mMemoryModel |= PKD_MODEL_GROUPS;
    if (parameters.get_bMemMass()) mMemoryModel |= PKD_MODEL_MASS;
    if (parameters.get_bMemSoft()) mMemoryModel |= PKD_MODEL_SOFTENING;
    if (parameters.get_bMemVelSmooth()) mMemoryModel |= PKD_MODEL_VELSMOOTH;
    if (parameters.get_bMemNewRung()) mMemoryModel |= PKD_MODEL_NEWRUNG;

    if (parameters.get_bMemNodeAcceleration()) mMemoryModel |= PKD_MODEL_NODE_ACCEL;
    if (parameters.get_bMemNodeVelocity()) mMemoryModel |= PKD_MODEL_NODE_VEL;
    if (parameters.get_bMemNodeMoment()) mMemoryModel |= PKD_MODEL_NODE_MOMENT;
    if (parameters.get_bMemNodeSphBounds()) mMemoryModel |= PKD_MODEL_NODE_SPHBNDS;

    if (parameters.get_bMemNodeBnd()) mMemoryModel |= PKD_MODEL_NODE_BND;
    if (parameters.get_bMemNodeVBnd()) mMemoryModel |= PKD_MODEL_NODE_VBND;

    if (MeshlessHydro()) mMemoryModel |= (PKD_MODEL_MFM | PKD_MODEL_ACCELERATION);
    if (MeshlessFiniteVolume()) mMemoryModel |= PKD_MODEL_MFV;

    if (parameters.get_bStars()) {
        mMemoryModel |= PKD_MODEL_STAR;
    }
    if (parameters.get_bBlackholes()) {
        mMemoryModel |= PKD_MODEL_BALL;
        mMemoryModel |= PKD_MODEL_BH;
    }

    if (parameters.get_bMemBall()) mMemoryModel |= PKD_MODEL_BALL;

    return mMemoryModel;
}

EphemeralMemory MSR::EphemeralMemoryGrid(int nGrid, int nCount) {
    return EphemeralMemory(mdl, nGrid, nCount);
}

std::tuple<int, int, int> MSR::InitializePStore(uint64_t* nSpecies, uint64_t mMemoryModel, uint64_t nEphemeral) {
    struct inInitializePStore ps;
    double dStorageAmount = (1.0 + parameters.get_dExtraStore());
    double dDecoratorStorageAmount = (1.0 + parameters.get_dExtraDecoratorStore());
    int i;
    for (i = 0; i <= FIO_SPECIES_LAST; ++i) ps.nSpecies[i] = nSpecies[i];
    ps.nStore = ceil(dStorageAmount * ps.nSpecies[FIO_SPECIES_ALL] / mdl->Threads());
    ps.nDecoratorStore =
        ceil(dDecoratorStorageAmount * (ps.nSpecies[FIO_SPECIES_SPH] + ps.nSpecies[FIO_SPECIES_BH] + ps.nSpecies[FIO_SPECIES_STAR]) / mdl->Threads());
    ps.nTreeBitsLo = parameters.get_nTreeBitsLo();
    ps.nTreeBitsHi = parameters.get_nTreeBitsHi();
    ps.iCacheSize = parameters.get_iCacheSize();
    ps.iCacheMaxInflight = parameters.get_iCacheMaxInflight();
    ps.iWorkQueueSize = parameters.get_iWorkQueueSize();
    ps.mpi_recv_mode = mdl::MPI_RECV_MODE(parameters.get_iMpiRecvMode());
    ps.mpi_recv_count = parameters.get_nMpiRecvCount();
    ps.mpi_RDMA = parameters.get_bMpiRDMA();
    ps.fPeriod = parameters.get_dPeriod();
    ps.mMemoryModel = mMemoryModel | PKD_MODEL_VELOCITY;
    ps.nIntegerFactor = parameters.get_nIntegerFactor();
    ps.nBucket = parameters.get_nBucket();

#define SHOW(m) ((ps.mMemoryModel & PKD_MODEL_##m) ? " " #m : "")
    print(
        "Memory Models:{position}{unordered}{velocity}{acceleration}{potential}{groups}{mass}{density}{ball}{new_rung}{softening}{velsmooth}{mfm}{mfv}{new_sph}"
        "{star}{particle_id}{bh}{globalgid}{node_moment}{node_accel}{node_vel}{node_sphbnds}{node_bnd}{node_vbnd}{node_bob}\n",
        "position"_a = parameters.get_bMemIntegerPosition() ? " INTEGER_POSITION" : " DOUBLE_POSITION", "unordered"_a = SHOW(UNORDERED),
        "velocity"_a = SHOW(VELOCITY), "acceleration"_a = SHOW(ACCELERATION), "potential"_a = SHOW(POTENTIAL), "groups"_a = SHOW(GROUPS), "mass"_a = SHOW(MASS),
        "density"_a = SHOW(DENSITY), "new_rung"_a = SHOW(NEWRUNG), "ball"_a = SHOW(BALL), "softening"_a = SHOW(SOFTENING), "velsmooth"_a = SHOW(VELSMOOTH),
        "mfm"_a = SHOW(MFM), "mfv"_a = SHOW(MFV), "new_sph"_a = SHOW(NEW_SPH), "star"_a = SHOW(STAR), "particle_id"_a = SHOW(PARTICLE_ID), "bh"_a = SHOW(BH),
        "globalgid"_a = SHOW(GLOBALGID), "node_moment"_a = SHOW(NODE_MOMENT), "node_accel"_a = SHOW(NODE_ACCEL), "node_vel"_a = SHOW(NODE_VEL),
        "node_sphbnds"_a = SHOW(NODE_SPHBNDS), "node_bnd"_a = SHOW(NODE_BND), "node_vbnd"_a = SHOW(NODE_VBND), "node_bob"_a = SHOW(NODE_BOB));
#undef SHOW

    // Calculate the Ephemeris memory requirements
    EphemeralMemory e(nEphemeral);
    if (parameters.get_bFindGroups()) e |= EphemeralMemory(4);
    if (parameters.get_bFindHopGroups()) e |= EphemeralMemory(8);
    if (parameters.get_iPkInterval()) e |= EphemeralMemory(4);
    if (parameters.get_bGravStep()) e |= EphemeralMemory(8);
    if (DoGas()) e |= EphemeralMemory(8);
    if (parameters.get_bMemBall()) e |= EphemeralMemory(8);
    if (parameters.get_bDoDensity()) e |= EphemeralMemory(12);
    if (parameters.get_bBlackholes()) e |= EphemeralMemory(8);

    // We need one grid to measure P(k); two if we are interlacing
    e |= EphemeralMemory(mdl, parameters.get_nGridPk(), parameters.get_bPkInterlace() ? 2 : 1);
    // Add some ephemeral memory (if needed) for the linGrid. 3 grids are stored : forceX, forceY, forceZ
    e |= EphemeralMemory(mdl, parameters.get_nGridLin(), 3);

    // Check all registered Python analysis routines and account for their memory requirements
    for (msr_analysis_callback& i : analysis_callbacks) {
        e |= EphemeralMemory(i.memory);
    }
    ps.nEphemeralBytes = e.per_particle;
    ps.nMinEphemeral = e.per_process;

    // Calculate constraint for generating initial conditions
    EphemeralMemory ic_memory(mdl, parameters.get_nGrid(), getGridInfoLPT(parameters.get_iLPT()).nGrids);
    ps.nMinTotalStore = ic_memory.per_process;

    outInitializePStore pout;
    pstInitializePStore(pst, &ps, sizeof(ps), &pout, sizeof(pout));
    PKD pkd = pst->lcl()->pkd;
    print("Allocated {} MB for {} particles on each processor.\n", pkd->ParticleMemory() / (1024 * 1024), pkd->FreeStore());
    print("Allocated {} MB for {} decorators on each processor.\n", pkd->DecoratorMemory() / (1024 * 1024), pkd->FreeDecoratorStore());
    print("Particles: {} bytes (persistent) + {} bytes (decorator) + {} bytes (ephemeral) , Nodes: {} bytes\n", pkd->particles.ParticleSize(),
          pkd->decorators.DecoratorSize(), ps.nEphemeralBytes, pkd->NodeSize());
    if (pkd->particles.ParticleSize() > MDL_CACHE_DATA_SIZE) {
        print_error("ERROR! MDL_CACHE_DATA_SIZE ({} bytes) is too small for the given particle size, please increasing it\n", MDL_CACHE_DATA_SIZE);
        abort();
    }
    if (ps.nMinEphemeral) print("Ephemeral will be at least {} MB per node.\n", ps.nMinEphemeral / (1024 * 1024));
    return std::make_tuple(pout.nSizeParticle, pout.nSizeDecorator, pout.nSizeNode);
}

/// @brief Return the number of threads to use for parallel reading/writing
/// @param bParallel True if parallel reading/writing is enabled
/// @param nParallel number of threads to use for parallel reading/writing
/// @return Actual number of threads to use for parallel reading/writing
int64_t MSR::parallel_count(int64_t nParallel, bool bParallel) {
    if (!bParallel) return 1;                           // Parallel reading/writing is disabled, use 1 thread
    else if (nParallel < 1) return nThreads;            // Use maximum number of threads
    else return std::min(nParallel, int64_t(nThreads)); // Use specified number of threads (if less than maximum)
}

/// @brief Return the number of threads to use for parallel reading
/// @return Number of threads to use for parallel reading
int64_t MSR::parallel_read_count() {
    return parallel_count(parameters.get_nParaRead(), parameters.get_bParaRead());
}

/// @brief Return the number of threads to use for parallel writing
/// @return Number of threads to use for parallel writing
int64_t MSR::parallel_write_count() {
    return parallel_count(parameters.get_nParaWrite(), parameters.get_bParaWrite());
}

void MSR::stat_files(std::vector<uint64_t>& counts, const std::string_view& filename_template, uint64_t element_size) {
    ServiceFileSizes::input hdr;
    marshal_string(hdr.filename, filename_template);
    hdr.nSimultaneous = hdr.nTotalActive = parallel_read_count();
    hdr.iReaderWriter = 0;
    hdr.nElementSize = 1;

    auto out = new ServiceFileSizes::output[ServiceFileSizes::max_files];
    auto out_span = std::span(out, ServiceFileSizes::max_files);
    auto n = mdl->RunService<ServiceFileSizes>(hdr, out_span);
    // auto n = mdl->RunService<ServiceFileSizes>(hdr, std::span(out, ServiceFileSizes::max_files));
    n /= sizeof(ServiceFileSizes::output);
    counts.resize(n);

    for (auto i = 0; i < n; ++i) {
        assert(out[i].iFileIndex < n);
        counts[out[i].iFileIndex] = out[i].nFileBytes / element_size;
    }
    delete[] out;
}

void MSR::Restore(const std::string& baseName, int nSizeParticle, int nSizeDecorator) {
    std::vector<uint64_t> counts;
    std::vector<uint64_t> counts_d;
    std::string filename_template = baseName + ".{i}";
    std::string filename_d_template = baseName + ".{i}.dec";

    print("Scanning Checkpoint files...\n");

    TimerStart(TIMER_NONE);
    stat_files(counts, filename_template, nSizeParticle);

    using mdl::ServiceBuffer;
    ServiceBuffer msg{ServiceBuffer::Field<ServiceInput::input>(), ServiceBuffer::Field<ServiceInput::io_elements>(counts.size())};
    auto hdr = static_cast<ServiceInput::input*>(msg.data(0));
    auto elements = static_cast<ServiceInput::io_elements*>(msg.data(1));
    hdr->nFiles = counts.size();
    hdr->nElements = std::accumulate(counts.begin(), counts.end(), uint64_t(0));
    std::copy(counts.begin(), counts.end(), elements);
    hdr->iBeg = 0;
    hdr->iEnd = hdr->nElements;
    marshal_string(hdr->io.filename, filename_template);
    hdr->io.nSimultaneous = parallel_read_count();
    hdr->io.nSegment = hdr->io.iThread = 0; // setup later
    hdr->io.iReaderWriter = 0;
    mdl->RunService<ServiceRestore>(msg);
    TimerStop(TIMER_NONE);
    auto dsec = TimerGet(TIMER_NONE);
    print("... identified {} particles in {} files, Wallclock: {:.5f} secs.\n", std::accumulate(counts.begin(), counts.end(), uint64_t(0)), int(counts.size()),
          dsec);

    // Repeat for decorators
    TimerStart(TIMER_NONE);
    stat_files(counts_d, filename_d_template, nSizeDecorator + nSizeParticle);
    ServiceBuffer msg_d{ServiceBuffer::Field<ServiceInput::input>(), ServiceBuffer::Field<ServiceInput::io_elements>(counts_d.size())};
    auto hdr_d = static_cast<ServiceInput::input*>(msg_d.data(0));
    auto elements_d = static_cast<ServiceInput::io_elements*>(msg_d.data(1));
    hdr_d->nFiles = counts_d.size();
    hdr_d->nElements = std::accumulate(counts_d.begin(), counts_d.end(), uint64_t(0));
    std::copy(counts_d.begin(), counts_d.end(), elements_d);
    hdr_d->iBeg = 0;
    hdr_d->iEnd = hdr_d->nElements;
    strncpy(hdr_d->io.filename, filename_d_template.c_str(), sizeof(hdr_d->io.filename));
    hdr_d->io.nSimultaneous = parallel_read_count();
    assert(hdr_d->io.nSimultaneous > 0);
    hdr_d->io.nSegment = hdr_d->io.iThread = 0; // setup later
    hdr_d->io.iReaderWriter = 0;
    mdl->RunService<ServiceRestoreDecorator>(msg_d);

    // mdl->RunService(PST_RESET_DECORATOR);
    TimerStop(TIMER_NONE);
    dsec = TimerGet(TIMER_NONE);
    print("... identified {} decorators in {} files, Wallclock: {:.5f} secs.\n", std::accumulate(counts_d.begin(), counts_d.end(), uint64_t(0)),
          int(counts_d.size()), dsec);
}

template <>
PyObject* MSR::restore(PyObject* file, PyObject* replace) {
    PyObject* args = PyTuple_Pack(1, file);
    auto result = PyObject_CallObject(pDill_load, args);
    Py_DECREF(args);
    if (!result) {
        PyErr_Print();
        abort();
    }
    // If we are replacing the parameter, throw away the old one
    // and return the new one
    if (replace != Py_None) {
        Py_DECREF(result);
        Py_INCREF(replace);
        return replace;
    }
    return result;
}

// Specialization for int
template <>
int MSR::restore(PyObject* file, PyObject* replace) {
    auto result = restore<PyObject*>(file, replace);
    auto value = PyLong_AsLong(result);
    Py_DECREF(result);
    return value;
}

// Specialization for double
template <>
double MSR::restore(PyObject* file, PyObject* replace) {
    auto result = restore<PyObject*>(file, replace);
    auto value = PyFloat_AsDouble(result);
    Py_DECREF(result);
    return value;
}

std::tuple<double, double, int64_t, int64_t, int64_t, int64_t> // dTime,dDelta,iStep,nSteps,nSizeParticle,nSizeNode
MSR::ReadCheckpoint(const char* filename, PyObject* kwargs, PyObject* species, PyObject* classes, PyObject* step, PyObject* steps, PyObject* time,
                    PyObject* delta, PyObject* E, PyObject* U, PyObject* Utime) {
    std::string pkl_filename = filename;
    pkl_filename += ".pkl";

    auto pFile = PyObject_CallFunction(PyDict_GetItemString(PyEval_GetBuiltins(), "open"), "ss", pkl_filename.c_str(), "rb");
    if (!pFile) {
        PyErr_Print();
        abort();
    }

    auto version = restore<int>(pFile);
    if (version < 1 || version > 2) {
        PyErr_SetString(PyExc_ValueError, "Invalid checkpoint file version");
        PyErr_Print();
        abort();
    }

    auto species_list = restore<PyObject*>(pFile, species);
    fioSpeciesList nSpecies;
    for (auto i = 0; i < FIO_SPECIES_LAST; ++i) {
        nSpecies[i] = PyLong_AsUnsignedLongLong(PyList_GetItem(species_list, i));
    }
    this->nDark = nSpecies[FIO_SPECIES_DARK];
    this->nGas = nSpecies[FIO_SPECIES_SPH];
    this->nStar = nSpecies[FIO_SPECIES_STAR];
    this->nBH = nSpecies[FIO_SPECIES_BH];
    this->N = nSpecies[FIO_SPECIES_ALL] = nDark + nGas + nStar + nBH;
    nMaxOrder = N - 1; // iOrder goes from 0 to N-1

    auto classes_list = restore<PyObject*>(pFile);
    static_assert(PKD_MAX_CLASSES <= 256); // Hopefully nobody will be mean to us (we use the stack)
    std::vector<PARTCLASS> CheckpointClasses;
    CheckpointClasses.reserve(PKD_MAX_CLASSES);
    auto nCheckpointClasses = PyList_Size(classes_list);
    for (int i = 0; i < nCheckpointClasses; ++i) {
        auto class_list = PyList_GetItem(classes_list, i); // Borrowed reference, no need to Py_DECREF
        auto eSpeciesObj = PyList_GetItem(class_list, 0);
        auto fMassObj = PyList_GetItem(class_list, 1);
        auto fSoftObj = PyList_GetItem(class_list, 2);
        auto iMatObj = PyList_GetItem(class_list, 3);
        auto accFac = PyList_Size(class_list) > 4 ? PyFloat_AsDouble(PyList_GetItem(class_list, 4)) : 1.0;
        CheckpointClasses.push_back(PARTCLASS{.fMass = static_cast<float>(PyFloat_AsDouble(fMassObj)),
                                              .fSoft = static_cast<float>(PyFloat_AsDouble(fSoftObj)),
                                              .iMat = static_cast<int>(PyLong_AsLong(iMatObj)),
                                              .accFac = static_cast<float>(accFac),
                                              .eSpecies = FIO_SPECIES(PyLong_AsLong(eSpeciesObj))});
    }

    auto iStep = restore<int>(pFile, step);
    auto nSteps = restore<int>(pFile, steps);
    auto dTime = restore<double>(pFile, time);
    auto dDelta = restore<double>(pFile, delta);
    this->dEcosmo = restore<double>(pFile, E);
    this->dUOld = restore<double>(pFile, U);
    this->dTimeOld = restore<double>(pFile, Utime);

    auto arguments = restore<PyObject*>(pFile);
    auto specified = restore<PyObject*>(pFile);
    parameters.merge(pkd_parameters(arguments, specified));
    parameters.update(kwargs, false);

    if (version >= 2) {
        auto analysis_list = restore<PyObject*>(pFile);
        analysis_callbacks.clear();
        for (auto i = 0; i < PyList_Size(analysis_list); ++i) {
            auto item = PyList_GetItem(analysis_list, i);
            auto callback = PyList_GetItem(item, 0);
            auto per_particle = PyLong_AsLong(PyList_GetItem(item, 1));
            auto per_process = PyLong_AsLong(PyList_GetItem(item, 2));
            analysis_callbacks.emplace_back(callback, per_particle, per_process);
        }
    }

    // The load_module function relies on peek(256) which is not guaranteed to return 256 bytes.
    // We seek to the start of the file and back to the current position as a workaround (best effort)
    auto where = PyObject_CallMethod(pFile, "tell", nullptr);
    if (!where) {
        PyErr_Print();
        abort();
    }
    Py_XDECREF(PyObject_CallMethod(pFile, "seek", "n", (Py_ssize_t)0));
    Py_XDECREF(PyObject_CallMethod(pFile, "seek", "O", where));
    Py_DECREF(where);

    // Restore the interpreter state so that functions that reference variables will still work
    PyObject* args = PyTuple_Pack(1, pFile);
    PyObject* result = PyObject_CallObject(pDill_load_module, args);
    Py_DECREF(args);
    if (!result) {
        PyErr_Print();
        abort();
    }
    Py_DECREF(result);

    PyObject_CallMethod(pFile, "close", NULL);
    Py_DECREF(pFile);

    ValidateParameters(); // Should be okay, but other stuff happens here (cosmo is setup for example)

    uint64_t mMemoryModel = 0;
    mMemoryModel = getMemoryModel();
    abort_assert(1001, nGas == 0 || DoGas());
    if (NewSPH()) mMemoryModel |= (PKD_MODEL_NEW_SPH | PKD_MODEL_ACCELERATION | PKD_MODEL_VELOCITY | PKD_MODEL_DENSITY | PKD_MODEL_BALL | PKD_MODEL_NODE_BOB);
    if (NewSPH() && parameters.get_bShearStrengthModel()) mMemoryModel |= PKD_MODEL_STRENGTH;
    if (parameters.get_bTrackPeakPressure()) mMemoryModel |= PKD_MODEL_PEAKPRESSURE;
    auto [nSizeParticle, nSizeDecorator, nSizeNode] = InitializePStore(nSpecies, mMemoryModel, parameters.get_nMemEphemeral());

    mdl->RunService<ServiceSetClasses>(CheckpointClasses);
    Restore(filename, nSizeParticle, nSizeDecorator);
    CalcBound();
    CountRungs(NULL);

    /* We can indicate that the DD was already done at rung 0 */
    iLastRungRT = 0;
    iLastRungDD = 0;

    InitCosmology(csm);

    SetDerivedParameters(true);

    if (parameters.has_dSoft()) SetSoft(Soft());

    if (NewSPH()) {
        /*
        ** Initialize kernel target with nSmooth
        */
        SetSPHoptions();
        InitializeEOS();
    }
    if (parameters.get_bAddDelete()) CountSpecies();
    return {dTime, dDelta, iStep, nSteps, nSizeParticle, nSizeNode};
}

void MSR::Restart(const char* filename, PyObject* kwargs, PyObject* species, PyObject* classes, PyObject* step, PyObject* steps, PyObject* time,
                  PyObject* delta, PyObject* E, PyObject* U, PyObject* Utime) {
    if (parameters.get_bVStart()) fmt::print("Restoring from checkpoint\n");
    TimerStart(TIMER_IO);

    auto [dTime, dDelta, iStep, nSteps, nSizeParticle, nSizeNode] = ReadCheckpoint(filename, kwargs, species, classes, step, steps, time, delta, E, U, Utime);

    TimerStop(TIMER_IO);
    auto dsec = TimerGet(TIMER_IO);
    double dExp = csmTime2Exp(csm, dTime);
    if (dsec > 0.0) {
        double rate = N * nSizeParticle / dsec;
        const char* units = "B";
        if (rate > 10000) {
            rate /= 1024;
            units = "KB";
        }
        if (rate > 10000) {
            rate /= 1024;
            units = "MB";
        }
        if (rate > 10000) {
            rate /= 1024;
            units = "GB";
        }
        print_detail("Checkpoint Restart Complete @ a={a:g}, Wallclock: {seconds:.5} secs ({rate:.2} {units}/s)\n", "a"_a = dExp, "seconds"_a = dsec,
                     "rate"_a = rate, "units"_a = units);
    }
    else print_detail("Checkpoint Restart Complete @ a={a:g}, Wallclock: {seconds:.5} secs\n\n", "a"_a = dExp, "seconds"_a = dsec);

    Simulate(dTime, iStep, true);
}

void MSR::persist(PyObject* file, PyObject* obj) {
    auto args = PyTuple_Pack(3, obj, file, Py_True);
    auto result = PyObject_CallObject(pDill_dump, args);
    Py_DECREF(args);
    if (!result) {
        PyErr_Print();
        abort();
    }
    Py_DECREF(result);
}
void MSR::persist(PyObject* file, int n) {
    auto number = PyLong_FromLongLong(n);
    persist(file, number);
    Py_DECREF(number);
}
void MSR::persist(PyObject* file, double d) {
    auto number = PyFloat_FromDouble(d);
    persist(file, number);
    Py_DECREF(number);
}

void MSR::writeParameters(const std::string& baseName, int iStep, int nSteps, double dTime, double dDelta) {
    fioSpeciesList nSpecies;
    int i;
    int nBytes;

    // ******************************************************************
    // Collect the information to checkpoint
    // ******************************************************************

    for (i = 0; i <= FIO_SPECIES_LAST; ++i) nSpecies[i] = 0;
    nSpecies[FIO_SPECIES_ALL] = N;
    nSpecies[FIO_SPECIES_SPH] = nGas;
    nSpecies[FIO_SPECIES_DARK] = nDark;
    nSpecies[FIO_SPECIES_STAR] = nStar;
    nSpecies[FIO_SPECIES_BH] = nBH;

    // Create a list with a count of each species
    auto species_list = PyList_New(FIO_SPECIES_LAST);
    for (int i = 0; i < FIO_SPECIES_LAST; ++i) {
        PyList_SetItem(species_list, i, PyLong_FromUnsignedLongLong(nSpecies[i])); // PyList_SetItem steals a reference to num
    }

    // Create a list with the checkpoint classes
    auto CheckpointClasses = GetClasses();
    auto classes_list = PyList_New(CheckpointClasses.size());
    for (int i = 0; i < CheckpointClasses.size(); ++i) {
        auto class_list = PyList_New(5); // Each inner list has 5 elements

        // Convert structure members to Python objects and add them to the class_list
        PyObject* eSpecies = PyLong_FromLong(CheckpointClasses[i].eSpecies);
        PyObject* fMass = PyFloat_FromDouble(CheckpointClasses[i].fMass);
        PyObject* fSoft = PyFloat_FromDouble(CheckpointClasses[i].fSoft);
        PyObject* iMat = PyLong_FromLong(CheckpointClasses[i].iMat);
        PyObject* accFac = PyFloat_FromDouble(CheckpointClasses[i].accFac);
        // Note: PyList_SetItem steals a reference to the item
        PyList_SetItem(class_list, 0, eSpecies);
        PyList_SetItem(class_list, 1, fMass);
        PyList_SetItem(class_list, 2, fSoft);
        PyList_SetItem(class_list, 3, iMat);
        PyList_SetItem(class_list, 4, accFac);

        // Add the inner list to the outer list
        PyList_SetItem(classes_list, i, class_list); // This also steals a reference
    }

    auto analysis_list = PyList_New(0);
    for (auto& cb : analysis_callbacks) {
        auto list = PyList_New(3); // Each inner list has 3 elements
        // Convert structure members to Python objects and add them to the class_list
        PyObject* callback = cb.callback;
        PyObject* per_particle = PyLong_FromLong(cb.memory.per_particle);
        PyObject* per_process = PyLong_FromLong(cb.memory.per_process);

        // Note: PyList_SetItem steals a reference to the item
        PyList_SetItem(list, 0, callback);
        PyList_SetItem(list, 1, per_particle);
        PyList_SetItem(list, 2, per_process);

        // Add the inner list to the outer list
        PyList_Append(analysis_list, list); // This also steals a reference
    }

    auto a = parameters.arguments();
    auto s = parameters.specified();

    // ******************************************************************
    // Write the interpreter state and the checkpoint variables to a file
    // ******************************************************************
    // Write the interpreter and the checkpoint variables to a file
    auto achOutName = baseName + ".pkl";
    auto pFile = PyObject_CallFunction(PyDict_GetItemString(PyEval_GetBuiltins(), "open"), "ss", achOutName.c_str(), "wb");
    if (!pFile) {
        PyErr_Print();
        abort();
    }

    // Checkpoint the important variables
    persist(pFile, 2);             // checkpoint version
    persist(pFile, species_list);  // 1: species
    persist(pFile, classes_list);  // 1: classes
    persist(pFile, iStep);         // 1: step
    persist(pFile, nSteps);        // 1: steps
    persist(pFile, dTime);         // 1: time
    persist(pFile, dDelta);        // 1: delta
    persist(pFile, dEcosmo);       // 1: E
    persist(pFile, dUOld);         // 1: U
    persist(pFile, dTimeOld);      // 1: Utime
    persist(pFile, a);             // 1: arguments
    persist(pFile, s);             // 1: specified
    persist(pFile, analysis_list); // 2: analysis callbacks

    // 1: Persist the interpreter state
    PyObject* args = PyTuple_Pack(3, pFile, Py_None, Py_True);
    PyObject* result = PyObject_CallObject(pDill_dump_module, args);
    Py_DECREF(args);
    if (!result) {
        PyErr_Print();
        abort();
    }
    Py_DECREF(result);

    PyObject_CallMethod(pFile, "close", NULL);
    Py_DECREF(pFile);

    Py_DECREF(species_list);
    Py_DECREF(classes_list);
    Py_DECREF(analysis_list);

    // ******************************************************************
    // Write the restart file
    // ******************************************************************
    std::ofstream restart_file(baseName);
    if (!restart_file) {
        perror(baseName.c_str());
        abort();
    }
    restart_file << "import PKDGRAV as msr\n";
    restart_file << "msr.restore(__file__)\n";
    restart_file.close();
}

void MSR::Checkpoint(int iStep, int nSteps, double dTime, double dDelta) {
    double dsec;

    auto filename = BuildCpName(iStep);
    auto pattern = filename + ".{i}";

    if (csm->val.bComove) {
        double dExp = csmTime2Exp(csm, dTime);
        print_detail("Writing checkpoint for Step: {step} Time:{time:g} Redshift:{z:g}\n", "step"_a = iStep, "time"_a = dTime, "z"_a = (1.0 / dExp - 1.0));
    }
    else {
        print_detail("Writing checkpoint for Step: {step} Time:{time:g}\n", "step"_a = iStep, "time"_a = dTime);
    }

    TimerStart(TIMER_IO);

    ServiceCheckpoint::input in{
        .nSimultaneous = static_cast<uint32_t>(parallel_write_count()),
        .iReaderWriter = 0,
        .iThread = 0,
    };
    abort_assert(1003, pattern.size() < sizeof(in.filename));
    marshal_string(in.filename, pattern);
    mdl->RunService<ServiceCheckpoint>(in);

    writeParameters(filename, iStep, nSteps, dTime, dDelta);

    /* This is not necessary, but it means the bounds will be identical upon restore */
    CalcBound();

    TimerStop(TIMER_IO);
    dsec = TimerGet(TIMER_IO);
    print_detail("Checkpoint has been successfully written, Wallclock: {seconds:.5} secs.\n", "seconds"_a = dsec);
}

void MSR::SetDerivedParameters(bool bRestart) {
    /**********************************************************************\
    * The following "parameters" are derived from real parameters.
    \**********************************************************************/

    /* Temperature-to-specific-internal-energy conversion factors for different values
     * of mean molecular weight (mu) */
    units = UNITS(parameters, csm->val.h);
    const double dTuPrefac = units.dGasConst / (parameters.get_dConstGamma() - 1.);
    dTuFac = dTuPrefac / parameters.get_dMeanMolWeight();

    const double dInvPrimNeutralMu = 0.25 + 0.75 * parameters.get_dInitialH();
    const double dInvPrimIonisedMu = 0.75 + 1.25 * parameters.get_dInitialH();
    dTuFacPrimNeutral = dTuPrefac * dInvPrimNeutralMu;
    dTuFacPrimIonised = dTuPrefac * dInvPrimIonisedMu;

    if (parameters.get_bDoCooling()) {
        SetCoolingParam();
    }
    if (parameters.get_bDoStarFormation()) SetStarFormationParam();
    if (DoStellarFeedback()) SetFeedbackParam();
    SetEOSParam();
    if (parameters.get_bBlackholes()) SetBlackholeParam();
#ifdef STELLAR_EVOLUTION
    SetStellarEvolutionParam();
#endif
}

void MSR::Initialize() {
    lcl.pkd = NULL;
    nThreads = mdl->Threads();
    lStart = time(0);
    fCenter = 0; // Center is at (0,0,0)
    /* Storage for output times*/

    iCurrMaxRung = 0;
    iRungDD = 0;
    iRungDT = 0;
    iLastRungRT = -1;
    iLastRungDD = -1; /* Domain decomposition is not done */
    nRung.resize(MAX_RUNG + 1, 0);
    csmInitialize(&csm);

    /*
    ** Create the processor subset tree.
    */
    if (nThreads > 1) {
        ServiceSetAdd::input inAdd(nThreads);
        mdl->RunService<ServiceSetAdd>(inAdd);
    }
}

void MSR::IgnoreSIGBUS() {
    mdl->RunService<ServiceIgnoreSIGBUS>();
}

int MSR::GetLock() {
    /*
    ** Attempts to lock run directory to prevent overwriting. If an old lock
    ** is detected with the same achOutName, an abort is signaled. Otherwise
    ** a new lock is created. The bOverwrite parameter flag can be used to
    ** suppress lock checking.
    */
    auto bOverwrite = parameters.get_bOverwrite();

    FILE* fp = NULL;
    char achTmp[PATH_MAX];

    auto achFile = MakePath(LOCKFILE);
    if (!bOverwrite && (fp = fopen(achFile.c_str(), "r"))) {
        if (fscanf(fp, "%s", achTmp) != 1) achTmp[0] = '\0';
        (void)fclose(fp);
        if (parameters.get_achOutName() == achTmp) {
            (void)print_error("ABORT: {} detected.\nPlease ensure data is safe to "
                              "overwrite. Delete lockfile and try again.\n",
                              achFile);
            return 0;
        }
    }
    if (!(fp = fopen(achFile.c_str(), "w"))) {
        if (bOverwrite && parameters.get_bVWarnings()) {
            (void)print_warning("WARNING: Unable to create {}...ignored.\n", achFile);
            return 1;
        }
        else {
            (void)print_error("Unable to create {}\n", achFile);
            return 0;
        }
    }
    (void)fprintf(fp, "%s", parameters.get_achOutName().data());
    (void)fclose(fp);
    return 1;
}

int MSR::CheckForStop(const char* achStopFile) {
    /*
    ** Checks for existence of STOPFILE in run directory. If found, the file
    ** is removed and the return status is set to 1, otherwise 0.
    */
    FILE* fp = NULL;
    auto achFile = MakePath(achStopFile);
    if ((fp = fopen(achFile.c_str(), "r"))) {
        (void)fmt::print("User interrupt detected.\n");
        (void)fclose(fp);
        (void)unlink(achFile.c_str());
        return 1;
    }
    return 0;
}

MSR::MSR(MDL mdl, PST pst) : pst(pst), mdl(static_cast<mdl::mdlClass*>(mdl)) {
}

MSR::~MSR() {
    csmFinish(csm);
    if (Py_IsInitialized()) Py_Finalize();
}

auto MSR::GetClasses() -> std::vector<PARTCLASS> {
    std::vector<PARTCLASS> classes(PKD_MAX_CLASSES);
    // Fetch a sorted list of particle classes
    auto nClassBytes = mdl->RunService<ServiceGetClasses>({}, classes);
    auto n = nClassBytes / sizeof(PARTCLASS);
    abort_assert(1002, n * sizeof(PARTCLASS) == nClassBytes);
    classes.resize(n);
    return classes;
}

void MSR::SetClasses() {
    auto classes = GetClasses();
    mdl->RunService<ServiceSetClasses>(classes);
}

void MSR::SwapClasses(int id) {
    auto pkd = pst->lcl()->pkd;
    pkd->swap_classes(id);
}

void MSR::OneNodeRead(struct inReadFile* in, FIO fio) {
    int id;
    LCL* plcl = pst->lcl();
    ServiceSwapAll::input inswap;
    int rID;
    fioSpeciesList nRead, nStart, nBasePart;
    for (auto s = FIO_SPECIES(0); s < FIO_SPECIES_LAST; s = FIO_SPECIES(s + 1)) {
        // NB: The first nBasePart particles are read by this thread at the end
        // of this function
        nBasePart[s] = (in->nSpeciesEnd[s] - in->nSpeciesStart[s]) / nThreads;
        nStart[s] = in->nSpeciesStart[s] + nBasePart[s];
        nRead[s] = 0;
    }

    for (id = 1; id < nThreads; ++id) {
        for (auto s = FIO_SPECIES(0); s < FIO_SPECIES_LAST; s = FIO_SPECIES(s + 1)) {
            if (nBasePart[s] == 0) continue;
            auto mod = (in->nSpeciesEnd[s] - in->nSpeciesStart[s]) % nThreads;

            auto nPart = nBasePart[s];
            if (id >= nThreads - mod) nPart++;
            /*
             * Read particles into the local storage.
             */
            assert(plcl->pkd->FreeStore() >= nPart);
            pkdReadFIO(plcl->pkd, fio, nStart[s], nPart, in->dvFac, in->dTuFac, FIO_SPECIES(s));
            nStart[s] += nPart;
            nRead[s] += nPart;
        }
        /*
         * Now shove them over to the remote processor.
         */
        SwapClasses(id);
        inswap.idSwap = 0;
        rID = mdl->ReqService<ServiceSwapAll>(id, inswap);
        pkdSwapAll(plcl->pkd, id);
        pst->mdl()->GetReply(rID);
    }
    for (auto s = FIO_SPECIES(0); s < FIO_SPECIES_LAST; s = FIO_SPECIES(s + 1)) {
        if (nBasePart[s] == 0) {
            // It may happen that we have so little particles that we can
            // not divide them by threads, so we just load it here
            auto mod = (in->nSpeciesEnd[s] - in->nSpeciesStart[s]) % nThreads;
            if (mod) {
                nBasePart[s] = mod;
            }
            else {
                continue;
            }
        }
        /*
         * Now read our own particles.
         */
        pkdReadFIO(plcl->pkd, fio, 0, nBasePart[s], in->dvFac, in->dTuFac, FIO_SPECIES(s));
        nRead[s] += nBasePart[s];
    }

    assert(nRead[FIO_SPECIES_DARK] == nDark);
    assert(nRead[FIO_SPECIES_SPH] == nGas);
    assert(nRead[FIO_SPECIES_STAR] == nStar);
    assert(nRead[FIO_SPECIES_BH] == nBH);
}

double MSR::SwitchDelta(double dTime, int iStep) {
    while (iStep_list.size() && iStep_list.front() <= iStep) {
        iStep_list.erase(iStep_list.begin());
        dDelta_list.erase(dDelta_list.begin());
    }
    if (dDelta_list.size() == 0) return 0.0;
    return dDelta_list.front();
}

double MSR::getTime(double dExpansion) {
    if (csm->val.bComove) {
        if (csm->val.dHubble0 == 0.0) {
            print_error("No hubble constant specified\n");
            Exit(1);
        }
        return csmExp2Time(csm, dExpansion);
    }
    else return dExpansion;
}

double MSR::getVfactor(double dExpansion) {
    return csm->val.bComove ? dExpansion * dExpansion : 1.0;
}

void MSR::RecvArray(void* vBuffer, PKD_FIELD field, int iUnitSize, double dTime, bool bMarked) {
    PKD pkd = pst->lcl()->pkd;
    ServiceSendArray::input in{
        .dvFac =
            [&] {
                if (csm->val.bComove) {
                    auto dExp = csmTime2Exp(csm, dTime);
                    return 1.0 / (dExp * dExp);
                }
                return 1.0;
            }(),
        .field = field,
        .iTo = 0, // Send to thread 0 (this thread)
        .iUnitSize = iUnitSize,
        .bMarked = bMarked,
    };
    int iIndex = 0;
    vBuffer = pkd->PackArray(iUnitSize * pkd->Local(), vBuffer, &iIndex, pkd->Local(), field, iUnitSize, in.dvFac, in.bMarked);
    assert(iIndex == pkd->Local());
    for (auto i = 1; i < nThreads; ++i) {
        auto rID = pkd->mdl->ReqService<ServiceSendArray>(i, in);
        vBuffer = pkd->RecvArray(i, vBuffer, iUnitSize);
        pkd->mdl->GetReply(rID);
    }
}

/*
** This function makes some potentially problematic assumptions!!!
** Main problem is that it calls pkd level routines, bypassing the
** pst level. It uses plcl pointer which is not desirable.
*/
void MSR::AllNodeWrite(std::string_view FileName, double dTime, double dvFac, int bDouble) {
    int nProcessors;
    struct inWrite in;

    /*
    ** Add Data Subpath for local and non-local names.
    */
    marshal_string(in.achOutFile, FileName);
    in.bStandard = parameters.get_bStandard();
    nProcessors = parallel_write_count();
    in.iIndex = 0;

    in.dTime = dTime;
    if (csm->val.bComove) {
        in.dExp = csmTime2Exp(csm, dTime);
        in.dvFac = 1.0 / (in.dExp * in.dExp);
    }
    else {
        in.dExp = 1.0;
        in.dvFac = 1.0;
    }

    /* We need to enforce periodic boundaries (when applicable) */
    auto period = parameters.get_dPeriod();
    if (parameters.get_bPeriodic() && blitz::all(period < FLOAT_MAXVAL)) {
        Bound::coord_type offset(0.5 * period);
        in.bnd = Bound(fCenter - offset, fCenter + offset);
    }
    else {
        in.bnd = Bound(Bound::coord_type(-std::numeric_limits<double>::max()), Bound::coord_type(std::numeric_limits<double>::max()));
    }

    in.dEcosmo = dEcosmo;
    in.dTimeOld = dTimeOld;
    in.dUOld = dUOld;
    in.dTuFac = dTuFac;
    in.units = units;
    in.dBoxSize = in.units.dKpcUnit * 1e-3 * csm->val.h;
    in.Omega0 = csm->val.dOmega0;
    in.OmegaLambda = csm->val.dLambda;
    in.HubbleParam = csm->val.h;

    in.nDark = nDark;
    in.nGas = nGas;
    in.nStar = nStar;
    in.nBH = nBH;

    in.bHDF5 = parameters.get_bHDF5();
    in.mFlags = FIO_FLAG_POTENTIAL | FIO_FLAG_DENSITY | (bDouble ? FIO_FLAG_CHECKPOINT : 0) | (parameters.get_bDoublePos() ? FIO_FLAG_DOUBLE_POS : 0) |
                (parameters.get_bDoubleVel() ? FIO_FLAG_DOUBLE_VEL : 0) | ((parameters.get_bMemParticleID() || MeshlessHydro()) ? FIO_FLAG_ID : 0) |
                (parameters.get_bMemMass() ? 0 : FIO_FLAG_COMPRESS_MASS) | (parameters.get_bMemSoft() ? 0 : FIO_FLAG_COMPRESS_SOFT) |
                (parameters.has_dSoft() ? FIO_FLAG_GLOBAL_SOFT : 0) | (parameters.get_bMemUnordered() ? FIO_FLAG_UNORDERED : 0);

    if (!in.bHDF5 && fmt::format(in.achOutFile, fmt::arg("i", 0)) == in.achOutFile) {
        FIO fio;
        fio = fioTipsyCreate(in.achOutFile, in.mFlags & FIO_FLAG_CHECKPOINT, in.bStandard, NewSPH() ? in.dTime : in.dExp, in.nGas, in.nDark, in.nStar);
        fioClose(fio);
    }
    in.iLower = 0;
    in.iUpper = nThreads;
    in.iIndex = 0;
    in.nProcessorsTotal = in.nProcessors = nProcessors;
    pstWrite(pst, &in, sizeof(in), NULL, 0);
}

auto MSR::SetTotal() -> uint64_t {
    ServiceSetTotal::output out;
    mdl->RunService<ServiceSetTotal>({}, out);
    return out;
}

auto MSR::CalcWriteStart() -> uint64_t {
    auto nTotal = SetTotal();
    assert(nTotal <= N);
    ServiceSetWriteStart::input nWriteStart = 0;
    mdl->RunService<ServiceSetWriteStart>(nWriteStart);
    return nTotal;
}

void MSR::Write(std::string_view FileName, double dTime, int bCheckpoint) {
    int nProcessors;
    double dvFac, dExp;
    double dsec;

#ifdef NAIVE_DOMAIN_DECOMP
    Reorder();
    ReorderDecorator();
#else
    if (iLastRungRT >= 0) Reorder();
#endif

    /*
    ** Calculate where to start writing.
    ** This sets plcl->nWriteStart.
    */
    /*uint64_t N =*/CalcWriteStart();

    nProcessors = parallel_write_count();
    auto bHDF5 = parameters.get_bHDF5();

    if (csm->val.bComove) {
        dExp = csmTime2Exp(csm, dTime);
        dvFac = 1.0 / (dExp * dExp);
    }
    else {
        dExp = dTime;
        dvFac = 1.0;
    }
    if (nProcessors == 1) {
        print_detail("Writing {name} in {format} format serially ...\n", "name"_a = FileName, "format"_a = (bHDF5 ? "HDF5" : "Tipsy"));
    }
    else {
        print_detail("Writing {name} in {format} format in parallel (but limited to {n}) ...\n", "name"_a = FileName, "format"_a = (bHDF5 ? "HDF5" : "Tipsy"),
                     "n"_a = nProcessors);
    }

    if (csm->val.bComove) print_detail("Time:{time:g} Redshift:{z:g}\n", "time"_a = dTime, "z"_a = (1.0 / dExp - 1.0));
    else print_detail("Time:{time:g}\n", "time"_a = dTime);

    TimerStart(TIMER_IO);
    AllNodeWrite(FileName, dTime, dvFac, bCheckpoint);
    TimerStop(TIMER_IO);
    dsec = TimerGet(TIMER_IO);

    print_detail("Output file has been successfully written, Wallclock: {seconds:.5} secs.\n", "seconds"_a = dsec);
}

void MSR::SetSoft(double dSoft) {
    print_detail("Set Softening...\n");
    ServiceSetSoft::input in(dSoft);
    mdl->RunService<ServiceSetSoft>(in);
}

// Use the Smooth routine with an empty smooth function to initialize fBall
void MSR::InitBall() {
    print("Computing a first guess for the smoothing length\n");

    bUpdateBall = 1;
    Smooth(1., 0.0, SMX_NULL, 0, parameters.get_nSmooth());
    bUpdateBall = 0;
}

void MSR::DomainDecompOld(int iRung) {
    OldDD::ServiceDomainDecomp::input in;
    uint64_t nActive;
    const uint64_t nDT = (N * parameters.get_dFracDualTree());
    const uint64_t nDD = (N * parameters.get_dFracNoDomainDecomp());
    const uint64_t nRT = (N * parameters.get_dFracNoDomainRootFind());
    const uint64_t nSD = (N * parameters.get_dFracNoDomainDimChoice());
    double dsec;
    int iRungDT, iRungDD = 0, iRungRT, iRungSD;
    int i;
    int bRestoreActive = 0;

    in.bDoRootFind = 1;
    in.bDoSplitDimFind = 1;
    if (iRung > 0) {
        /*
        ** All of this could be calculated once for the case that the number
        ** of particles don't change. Or calculated every time the number of
        ** particles does change.
        */
        nActive = 0;
        iRungDT = 0;
        iRungDD = 0;
        iRungRT = 0;
        iRungSD = 0;
        for (i = iCurrMaxRung; i >= 0; --i) {
            nActive += nRung[i];
            if (nActive > nDT && !iRungDT) iRungDT = i;
            if (nActive > nDD && !iRungDD) iRungDD = i;
            if (nActive > nRT && !iRungRT) iRungRT = i;
            if (nActive > nSD && !iRungSD) iRungSD = i;
        }
        assert(iRungDD >= iRungRT);
        assert(iRungRT >= iRungSD);
#ifdef NAIVE_DOMAIN_DECOMP
        if (iLastRungRT < 0) {
            /*
            ** We need to do a full domain decompotition with iRungRT particles being active.
            ** However, since I am not sure what the exact state of the domains can be at this point
            ** I had better do a split dim find as well.
            */
            iLastRungRT = 0;
            ActiveRung(iLastRungRT, 1);
            bRestoreActive = 1;
            in.bDoRootFind = 1;
            in.bDoSplitDimFind = 1;
            if (parameters.get_bVRungStat()) {
                print("Doing Domain Decomposition (nActive = {}/{} iRung:{} iRungRT:{})\n", nActive, N, iRung, iRungRT);
            }
        }
        else if (iRung <= iRungRT) {
            /*
            ** We need to do a full domain decomposition with *ALL* particles being active.
            */
            in.bDoRootFind = 1;
            if (iRung <= iRungSD) {
                if (parameters.get_bVRungStat()) {
                    print("Doing Domain Decomposition (nActive = {}/{}, iRung:{} iRungRT:{})\n", nActive, N, iRung, iRungRT);
                }
                in.bDoSplitDimFind = 1;
            }
            else {
                if (parameters.get_bVRungStat()) {
                    print("Skipping Domain Dim Choice (nActive = {}/{}, iRung:{} iRungSD:{})\n", nActive, N, iRung, iRungSD);
                }
                in.bDoSplitDimFind = 0;
            }
            ActiveRung(0, 1); /* Here we activate all particles. */
            bRestoreActive = 1;
        }
        else if (iRung <= iRungDD) {
            if (parameters.get_bVRungStat()) {
                print("Skipping Root Finder (nActive = {}/{}, iRung:{} iRungRT:{} iRungDD:{})\n", nActive, N, iRung, iRungRT, iRungDD);
            }
            in.bDoRootFind = 0;
            in.bDoSplitDimFind = 0;
            bRestoreActive = 0;
        }
        else {
            if (parameters.get_bVRungStat()) {
                print("Skipping Domain Decomposition (nActive = {}/{}, iRung:{} iRungDD:{})\n", nActive, N, iRung, iRungDD);
            }
            return; /* do absolutely nothing! */
        }
#else
        if (iLastRungRT < 0) {
            /*
            ** We need to do a full domain decompotition with iRungRT particles being active.
            ** However, since I am not sure what the exact state of the domains can be at this point
            ** I had better do a split dim find as well.
            */
            iLastRungRT = iRungRT;
            ActiveRung(iRungRT, 1);
            bRestoreActive = 1;
            in.bDoRootFind = 1;
            in.bDoSplitDimFind = 1;
        }
        else if (iRung == iLastRungDD) {
            if (parameters.get_bVRungStat()) {
                print("Skipping Domain Decomposition (nActive = {}/{}, iRung:{} iRungDD:{} iLastRungRT:{})\n", nActive, N, iRung, iRungDD, iLastRungRT);
            }
            return; /* do absolutely nothing! */
        }
        else if (iRung >= iRungDD && !bSplitVA) {
            if (iLastRungRT < iRungRT) {
                iLastRungRT = iRungRT;
                ActiveRung(iRungRT, 1);
                bRestoreActive = 1;
                in.bDoRootFind = 1;
                in.bDoSplitDimFind = 0;
            }
            else {
                if (parameters.get_bVRungStat()) {
                    print("Skipping Domain Decomposition (nActive = {}/{}, iRung:{} iRungDD:{} iLastRungRT:{})\n", nActive, N, iRung, iRungDD, iLastRungRT);
                }
                return; /* do absolutely nothing! */
            }
        }
        else if (iRung > iRungRT) {
            if (iLastRungRT < iRungRT) {
                iLastRungRT = iRungRT;
                ActiveRung(iRungRT, 1);
                bRestoreActive = 1;
                in.bDoRootFind = 1;
                in.bDoSplitDimFind = 0;
            }
            else {
                if (parameters.get_bVRungStat()) {
                    print("Skipping Root Finder (nActive = {}/{}, iRung:{} iRungRT:{} iRungDD:{} iLastRungRT:{})\n", nActive, N, iRung, iRungRT, iRungDD,
                          iLastRungRT);
                }
                in.bDoRootFind = 0;
                in.bDoSplitDimFind = 0;
            }
        }
        else if (iRung > iRungSD) {
            if (iLastRungRT == iRung) {
                if (parameters.get_bVRungStat()) {
                    print("Skipping Root Finder (nActive = {}/{}, iRung:{} iRungRT:{} iRungDD:{} iLastRungRT:{})\n", nActive, N, iRung, iRungRT, iRungDD,
                          iLastRungRT);
                }
                in.bDoRootFind = 0;
                in.bDoSplitDimFind = 0;
            }
            else {
                if (parameters.get_bVRungStat()) {
                    print("Skipping Domain Dim Choice (nActive = {}/{}, iRung:{} iRungSD:{} iLastRungRT:{})\n", nActive, N, iRung, iRungSD, iLastRungRT);
                }
                iLastRungRT = iRung;
                in.bDoRootFind = 1;
                in.bDoSplitDimFind = 0;
            }
        }
        else {
            if (iLastRungRT == iRung) {
                in.bDoRootFind = 0;
                in.bDoSplitDimFind = 0;
            }
            else {
                iLastRungRT = iRung;
                in.bDoRootFind = 1;
                in.bDoSplitDimFind = 1;
            }
        }
#endif
    }
    else nActive = N;
    iLastRungDD = iLastRungRT;
    in.nActive = nActive;
    in.nTotal = N;

    in.nBndWrap[0] = 0;
    in.nBndWrap[1] = 0;
    in.nBndWrap[2] = 0;
    /*
    ** If we are dealing with a nice periodic volume in all
    ** three dimensions then we can set the initial bounds
    ** instead of calculating them.
    */
    auto period = parameters.get_dPeriod();
    if (parameters.get_bPeriodic() && blitz::all(period < FLOAT_MAXVAL)) {
        Bound::coord_type offset(0.5 * period);
        in.bnd = Bound(fCenter - offset, fCenter + offset);
        mdl->RunService<ServiceEnforcePeriodic>(in.bnd);
    }
    else {
        mdl->RunService<ServiceCombineBound>({}, in.bnd);
    }
    /* We make sure that the classes are synchronized among all the domains,
     * otherwise a new class type being moved to another DD region could cause
     * very nasty bugs!
     */
    SetClasses();

    print_detail("Domain Decomposition: nActive (Rung {rung}) {active}\n", "rung"_a = iLastRungRT, "active"_a = nActive);
    print_detail("Domain Decomposition... \n");
    TimerStart(TIMER_DOMAIN);

    mdl->RunService<OldDD::ServiceDomainDecomp>(in);
    ReorderDecorator();
    TimerStop(TIMER_DOMAIN);
    dsec = TimerGet(TIMER_DOMAIN);
    print("Domain Decomposition complete, Wallclock: {:.5} secs\n\n", dsec);
    if (bRestoreActive) {
        /* Restore Active data */
        ActiveRung(iRung, 1);
    }
}

void MSR::DomainDecomp(int iRung) {
    if (iRung > parameters.get_iMaxRungDomainDecomp()) {
        print_detail("Skipping Domain Decomposition, as rung {rung} is above the maximum defined by iMaxRungDomainDecomp ({maxrung}).\n", "rung"_a = iRung,
                     "maxrung"_a = parameters.get_iMaxRungDomainDecomp());
        return;
    }
    DomainDecompOld(iRung);
}

/*
** This the meat of the tree build, but will be called by differently named
** functions in order to implement special features without recoding...
*/
void MSR::BuildTree(int bNeedEwald, uint32_t uRoot, uint32_t utRoot) {
    ServiceBuildTree::input in;
    const double ddHonHLimit = parameters.get_ddHonHLimit();
    LCL* plcl = pst->lcl();
    PKD pkd = plcl->pkd;
    double dsec;

    auto nTopTree = pkd->NodeSize() * (2 * nThreads - 1);
    auto nMsgSize = sizeof(ServiceDistribTopTree::input) + nTopTree;

    auto buffer = std::make_unique<std::byte[]>(nMsgSize);
    auto inDistribTop = std::span(buffer.get(), nMsgSize);
    auto pDistribTop = new (buffer.get()) ServiceDistribTopTree::input;
    auto topTree = std::span(reinterpret_cast<std::byte*>(pDistribTop + 1), nTopTree);
    auto pkdn = pkd->tree[reinterpret_cast<KDN*>(pDistribTop + 1)];

    pDistribTop->uRoot = uRoot;
    pDistribTop->allocateMemory = 1;

    in.nBucket = parameters.get_nBucket();
    in.nGroup = parameters.get_nGroup();
    in.uRoot = uRoot;
    in.utRoot = utRoot;
    in.ddHonHLimit = ddHonHLimit;
    TimerStart(TIMER_TREE);
    nTopTree = mdl->RunService<ServiceBuildTree>(in, topTree);
    // nTopTree = pstBuildTree(pst, &in, sizeof(in), pkdn, nTopTree);
    pDistribTop->nTop = nTopTree / pkd->NodeSize();
    assert(pDistribTop->nTop == (2 * nThreads - 1));
    mdl->RunService<ServiceDistribTopTree>(inDistribTop);
    TimerStop(TIMER_TREE);
    dsec = TimerGet(TIMER_TREE);
    print("Tree built, Wallclock: {:.5} secs\n\n", dsec);

    if (bNeedEwald) {
        /*
        ** For simplicity we will skip calculating the Root for all particles
        ** with exclude very active since there are missing particles which
        ** could add to the mass and because it probably is not important to
        ** update the root so frequently.
        */
        ServiceCalcRoot::input calc;
        ServiceCalcRoot::output root;
        calc.com = pkdn->position();
        calc.uRoot = uRoot;

        mdl->RunService<ServiceCalcRoot>(calc, root);
        momTreeRoot[uRoot] = root.momc;
        momTreeCom[uRoot][0] = calc.com[0];
        momTreeCom[uRoot][1] = calc.com[1];
        momTreeCom[uRoot][2] = calc.com[2];
    }
}

void MSR::BuildTree(int bNeedEwald) {
    print_detail("Building local trees...\n\n");
    ServiceDumpTrees::input dump(IRUNGMAX);
    mdl->RunService<ServiceDumpTrees>(dump);
    BuildTree(bNeedEwald, ROOT, 0);

    if (bNeedEwald) {
        ServiceDistribRoot::input droot;
        droot.momc = momTreeRoot[ROOT];
        droot.r[0] = momTreeCom[ROOT][0];
        droot.r[1] = momTreeCom[ROOT][1];
        droot.r[2] = momTreeCom[ROOT][2];
        mdl->RunService<ServiceDistribRoot>(droot);
    }
#ifdef OPTIM_REORDER_IN_NODES
    if (MeshlessHydro()) {
        ReorderInNodes();
    }
#endif
}

/*
** Separates the particles into two trees, and builds the "fixed" tree.
*/
void MSR::BuildTreeFixed(int bNeedEwald, uint8_t uRungDD) {
    print_detail("Building fixed local trees...\n\n");
    BuildTree(bNeedEwald, FIXROOT, 0);
}

void MSR::BuildTreeActive(int bNeedEwald, uint8_t uRungDD) {
    /*
     ** The trees reset/removed. This does the following:
     **   1. Closes any open cell cache (it will be subsequently invalid)
     **   2. Resets the number of used nodes to zero (or more if we keep the main tree)
     **   3. Sets up the ROOT and FIXROOT node (either of which may have zero particles).
     */

    print_detail("Building active local trees...\n\n");

    ServiceDumpTrees::input dump(uRungDD, true);
    mdl->RunService<ServiceDumpTrees>(dump);

    /* New build the very active tree */
    BuildTree(bNeedEwald, ROOT, FIXROOT);

    /* For ewald we have to shift and combine the individual tree moments */
    if (bNeedEwald) {
        ServiceDistribRoot::input droot;
        MOMC momc;
        double* com1 = momTreeCom[FIXROOT];
        double m1 = momTreeRoot[FIXROOT].m;
        double* com2 = momTreeCom[ROOT];
        double m2 = momTreeRoot[ROOT].m;
        double ifMass = 1.0 / (m1 + m2);
        double x, y, z;
        int j;

        /* New Center of Mass, then shift and scale the moments */
        for (j = 0; j < 3; ++j) droot.r[j] = ifMass * (m1 * com1[j] + m2 * com2[j]);

        droot.momc = momTreeRoot[FIXROOT];
        x = com1[0] - droot.r[0];
        y = com1[1] - droot.r[1];
        z = com1[2] - droot.r[2];
        momShiftMomc(&droot.momc, x, y, z);

        momc = momTreeRoot[ROOT];
        x = com2[0] - droot.r[0];
        y = com2[1] - droot.r[1];
        z = com2[2] - droot.r[2];
        momShiftMomc(&momc, x, y, z);

        momAddMomc(&droot.momc, &momc);

        mdl->RunService<ServiceDistribRoot>(droot);
    }
}

void MSR::BuildTreeMarked(int bNeedEwald) {
    ServiceDumpTrees::input dump(IRUNGMAX);
    mdl->RunService<ServiceDumpTrees>(dump);

    mdl->RunService<ServiceTreeInitMarked>();
    BuildTree(bNeedEwald, FIXROOT, 0);
    BuildTree(bNeedEwald, ROOT, FIXROOT);

    /* For ewald we have to shift and combine the individual tree moments */
    if (bNeedEwald) {
        ServiceDistribRoot::input droot;
        MOMC momc;
        double* com1 = momTreeCom[FIXROOT];
        double m1 = momTreeRoot[FIXROOT].m;
        double* com2 = momTreeCom[ROOT];
        double m2 = momTreeRoot[ROOT].m;
        double ifMass = 1.0 / (m1 + m2);
        double x, y, z;
        int j;

        /* New Center of Mass, then shift and scale the moments */
        for (j = 0; j < 3; ++j) droot.r[j] = ifMass * (m1 * com1[j] + m2 * com2[j]);

        droot.momc = momTreeRoot[FIXROOT];
        x = com1[0] - droot.r[0];
        y = com1[1] - droot.r[1];
        z = com1[2] - droot.r[2];
        momShiftMomc(&droot.momc, x, y, z);

        momc = momTreeRoot[ROOT];
        x = com2[0] - droot.r[0];
        y = com2[1] - droot.r[1];
        z = com2[2] - droot.r[2];
        momShiftMomc(&momc, x, y, z);

        momAddMomc(&droot.momc, &momc);

        mdl->RunService<ServiceDistribRoot>(droot);
    }
}

void MSR::Reorder() {
    if (!parameters.get_bMemUnordered()) {
        double sec, dsec;

        print_detail("Ordering...\n");
        sec = Time();
#ifdef NEW_REORDER
        // Start by dividing the particles by processor; cores will follow
        auto nPerProc = (N + mdl->Procs() - 1) / mdl->Procs();
        print("Divided {} particles into {} domains ({})\n", N, mdl->Procs(), nPerProc);
        NewDD::ServiceReorder::input indomain(nPerProc, MaxOrder());
        mdl->RunService<NewDD::ServiceReorder>(indomain);
#else
        OldDD::ServiceDomainOrder::input indomain(MaxOrder());
        mdl->RunService<OldDD::ServiceDomainOrder>(indomain);
        OldDD::ServiceLocalOrder::input inlocal(MaxOrder());
        mdl->RunService<OldDD::ServiceLocalOrder>(inlocal);
#endif
        dsec = Time() - sec;
        print_detail("Order established, Wallclock {seconds} secs\n\n", "seconds"_a = dsec);

        /*
        ** Mark domain decomp as not done.
        */
        iLastRungRT = -1;
        iLastRungDD = -1;
    }
}

void MSR::OutASCII(std::string_view File, int iType, int nDims, int iFileType) {
    LCL* plcl = pst->lcl();
    int id, iDim;
    ServiceSwapAll::input inswap;
    PKDOUT pkdout;
    const char* arrayOrVector;
    int rID;

    switch (nDims) {
    case 1:
        arrayOrVector = "array";
        break;
    case 3:
        arrayOrVector = "vector";
        break;
    default:
        arrayOrVector = NULL;
        assert(nDims == 1 || nDims == 3);
    }

    auto total = SetTotal();

    /*
    ** Add Data Subpath for local and non-local names.
    */
    std::string achOutFile{File};

    switch (iFileType) {
#ifdef HAVE_LIBBZ2
    case PKDOUT_TYPE_BZIP2:
        achOutFile += ".bz2";
        break;
#endif
#ifdef HAVE_LIBZ
    case PKDOUT_TYPE_ZLIB:
        achOutFile += ".gz";
        break;
#endif
    default:
        break;
    }
    print_detail("Writing {name} to {file}\n", "name"_a = arrayOrVector, "file"_a = achOutFile);

    if (parallel_write_count() > 1 && iFileType > 1) {
        FILE* fp;

        fp = fopen(achOutFile.c_str(), "wb");
        if (fp == NULL) {
            print_error("Could not create {} Output File:{}\n", arrayOrVector, achOutFile);
            Exit(1);
        }
        fclose(fp);

        for (iDim = 0; iDim < nDims; iDim++) {
            ServiceCompressASCII::input in{.nTotal = total, .iFile = iFileType, .iType = iType, .iDim = iDim};
            ServiceCompressASCII::output out;
            mdl->RunService<ServiceCompressASCII>(in, out);
            ServiceWriteASCII::input inWrite;
            marshal_string(inWrite.achOutFile, achOutFile);
            mdl->RunService<ServiceWriteASCII>(inWrite);
        }
    }
    else {
        pkdout = pkdOpenOutASCII(plcl->pkd, achOutFile.c_str(), "wb", iFileType, iType);
        if (!pkdout) {
            print_error("Could not open {} Output File:{}\n", arrayOrVector, achOutFile);
            Exit(1);
        }

        pkdOutHdr(plcl->pkd, pkdout, total);

        /*
         * First write our own particles.
         */
        for (iDim = 0; iDim < nDims; ++iDim) {
            pkdOutASCII(plcl->pkd, pkdout, iType, iDim);
            for (id = 1; id < nThreads; ++id) {
                /*
                 * Swap particles with the remote processor.
                 */
                inswap.idSwap = 0;
                rID = mdl->ReqService<ServiceSwapAll>(id, inswap);
                pkdSwapAll(plcl->pkd, id);
                pst->mdl()->GetReply(rID);
                /*
                 * Write the swapped particles.
                 */
                pkdOutASCII(plcl->pkd, pkdout, iType, iDim);
                /*
                 * Swap them back again.
                 */
                inswap.idSwap = 0;
                rID = mdl->ReqService<ServiceSwapAll>(id, inswap);
                pkdSwapAll(plcl->pkd, id);
                pst->mdl()->GetReply(rID);
            }
        }
        pkdCloseOutASCII(plcl->pkd, pkdout);
    }
}

void MSR::OutArray(std::string_view File, int iType, int iFileType) {
    OutASCII(File, iType, 1, iFileType);
}
void MSR::OutArray(std::string_view File, int iType) {
    OutArray(File, iType, parameters.get_iCompress());
}

void MSR::OutVector(std::string_view File, int iType, int iFileType) {
    OutASCII(File, iType, 3, iFileType);
}
void MSR::OutVector(std::string_view File, int iType) {
    OutVector(File, iType, parameters.get_iCompress());
}

auto MSR::get_eEOS() -> eEOSparam {
    return eEOSparam{.dPolyFloorMinOD = calc.dEOSPolyFloorMinBaryonOD,
                     .dPolyFloorExponent = calc.dEOSPolyFloorExponent,
                     .dPolyFloorDen = calc.dEOSPolyFloorDen,
                     .dPolyFlooru = calc.dEOSPolyFlooru,
#ifdef EEOS_JEANS
                     .dNJeans = parameters.get_dEOSNJeans(),
#endif
                     .dFlooru = calc.dEOSFlooru,
                     .dFloorDen = calc.dEOSFloorDen,
                     .dFloorMinOD = calc.dEOSFloorMinBaryonOD,
                     .bDoEOS = parameters.get_bDoEOS()};
}

void MSR::Smooth(double dTime, double dDelta, int iSmoothType, int bSymmetric, int nSmooth) {
    ServiceSmooth::input in{.nSmooth = nSmooth, .bPeriodic = parameters.get_bPeriodic(), .bSymmetric = bSymmetric, .iSmoothType = iSmoothType};
    uint64_t nSmoothed;
    double sec, dsec;
    if (parameters.get_bVStep()) print("Smoothing...\n");
    sec = MSR::Time();
    mdl->RunService<ServiceSmooth>(in, nSmoothed);
    dsec = MSR::Time() - sec;
    if (parameters.get_bVStep()) print("Smooth Calculated, Wallclock: {:.5} secs\n\n", dsec);
}

uint64_t MSR::ReSmooth(double dTime, double dDelta, int iSmoothType, int bSymmetric) {
    ServiceReSmooth::input in{
        .nSmooth = static_cast<int>(parameters.get_nSmooth()), .bPeriodic = parameters.get_bPeriodic(), .bSymmetric = bSymmetric, .iSmoothType = iSmoothType};
    ServiceReSmooth::output out;
    double sec, dsec;
    if (parameters.get_bVStep()) print("ReSmoothing...\n");
    sec = MSR::Time();
    mdl->RunService<ServiceReSmooth>(in, out);
    dsec = MSR::Time() - sec;
    if (parameters.get_bVStep()) print("ReSmooth Calculated, Wallclock: {:.5} secs\n\n", dsec);
    return out;
}

void MSR::ReorderDecorator() {
    printf("Reordering decorator...\n");
    mdl->RunService<ServiceReorderDecorator>();
}

#ifdef OPTIM_REORDER_IN_NODES
void MSR::ReorderInNodes() {
    double dsec;
    TimerStart(TIMER_TREE);
    mdl->RunService<ServiceReorderInNodes>();

    TimerStop(TIMER_TREE);
    dsec = TimerGet(TIMER_TREE);
    print("Reordering nodes took {:.5} secs \n", dsec);
}
#endif

void MSR::UpdateSoft(double dTime) {
    const auto dSoft = parameters.get_dSoft();
    const auto dSoftMax = parameters.get_dSoftMax();
    const auto dMaxPhysicalSoft = parameters.get_dMaxPhysicalSoft();
    if (parameters.get_bPhysicalSoft()) {
        struct inPhysicalSoft in;

        in.dFac = 1. / csmTime2Exp(csm, dTime);
        in.bSoftMaxMul = parameters.get_bSoftMaxMul();
        in.dSoftMax = dSoftMax;

        if (in.bSoftMaxMul && in.dFac > in.dSoftMax) in.dFac = in.dSoftMax;

        pstPhysicalSoft(pst, &in, sizeof(in), NULL, 0);
    }
    if (dMaxPhysicalSoft > 0) {
        double dFac = csmTime2Exp(csm, dTime);
        if (dSoft > 0.0) {
            // Change the global softening
            if (dSoft * dFac > dMaxPhysicalSoft) SetSoft(dMaxPhysicalSoft / dFac);
        }
        else {
            // Individual (either particle or classes) softening are used
            struct inPhysicalSoft in;

            in.bSoftMaxMul = parameters.get_bSoftMaxMul();
            in.dSoftMax = dSoftMax;

            if (in.bSoftMaxMul) {
                const float fPhysFac = dFac; // Factor to convert to physical
                if (fPhysFac < dMaxPhysicalSoft) {
                    // nothing happens, still in the comoving softening regime
                    in.dFac = 1.;
                }
                else {
                    // late-times, softening limited by physical
                    in.dFac = dMaxPhysicalSoft / dFac;
                }
            }

            pstPhysicalSoft(pst, &in, sizeof(in), NULL, 0);
        }
    }
}

template <typename T, typename M>
void print_grid(int width, int nThreads, const std::unique_ptr<T[]>& data, M T::* member) {
    fmt::print(fmt::emphasis::bold, "      {:>{w}} {:>{w}} {:>{w}} {:>{w}} {:>{w}} {:>{w}} {:>{w}} {:>{w}} {:>{w}} {:>{w}}\n", 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
               "w"_a = width);
    for (auto i = 0; i < nThreads; i += 10) {
        std::string line = fmt::format(fmt::emphasis::bold, "{:4d}:", i);
        auto end = std::min(i + 10, nThreads);
        for (auto j = i; j < end; ++j) {
            line += fmt::format(" {:>{w}}", data[j].*member, "w"_a = width);
        }
        fmt::print("{}\n", line);
    }
}

void MSR::Hostname() {
    std::unique_ptr<ServiceHostname::output[]> out{new ServiceHostname::output[nThreads]};
    std::span out_span(out.get(), nThreads);
    mdl->RunService<ServiceHostname>({}, out_span);
    print("Host Names:\n");
    print_grid(12, nThreads, out, &ServiceHostname::output::szHostname);
    print("MPI Rank:\n");
    print_grid(8, nThreads, out, &ServiceHostname::output::iMpiID);
}

void MSR::MemStatus() {
    if (parameters.get_bVDetails()) {
        auto out = std::make_unique<ServiceMemStatus::output[]>(nThreads);
        mdl->RunService<ServiceMemStatus>({}, std::span(out.get(), nThreads));
        // Normalize to an appropriate unit
        for (auto i = 0; i < nThreads; ++i) {
            out[i].nBytesTree = out[i].nBytesTree / 1024 / 1024;
            out[i].nBytesCl = out[i].nBytesCl / 1024;
            out[i].nBytesIlp = out[i].nBytesIlp / 1024;
            out[i].nBytesIlc = out[i].nBytesIlc / 1024;
        }
#ifdef __linux__
        print("Resident (MB):\n");
        print_grid(8, nThreads, out, &ServiceMemStatus::output::rss);
        print("Free Memory (MB):\n");
        print_grid(8, nThreads, out, &ServiceMemStatus::output::freeMemory);
#endif
        print("Tree size (MB):\n");
        print_grid(8, nThreads, out, &ServiceMemStatus::output::nBytesTree);
        print("Checklist size (KB):\n");
        print_grid(8, nThreads, out, &ServiceMemStatus::output::nBytesCl);
        print("Particle List size (KB):\n");
        print_grid(8, nThreads, out, &ServiceMemStatus::output::nBytesIlp);
        print("Cell List size (KB):\n");
        print_grid(8, nThreads, out, &ServiceMemStatus::output::nBytesIlc);
    }
}

void MSR::PrintStat(const STAT& ps, char const* pszPrefix, int p) {
    double dSum = ps.dSum;
    double dMax = ps.dMax;
    const char* minmax = "max";

    if (dSum < 0) {
        dSum = -dSum;
        dMax = -dMax;
        minmax = "min";
    }
    if (ps.n > 1) {
        print("{} {}={:>8.{prec}f} @{:5d} avg={:>8.{prec}f} of {:5d} std-dev={:8.{prec}f}\n", pszPrefix, minmax, dMax, ps.idMax, dSum / ps.n, ps.n,
              sqrt((ps.dSum2 - ps.dSum * ps.dSum / ps.n) / (ps.n - 1)), "prec"_a = p);
    }
    else if (ps.n == 1) {
        print("{} {}={:>8.{prec}f} @{:5d}\n", pszPrefix, minmax, dMax, ps.idMax, "prec"_a = p);
    }
    else {
        print("{} no data\n", pszPrefix);
    }
}

#ifdef INSTRUMENT_MPI
void MSR::PrintStat(const WelfordStats<double>& in_seconds, char const* pszPrefix, int p) {
    auto us = in_seconds * 1e6;
    const char* minmax = "max";
    print("{} {}={:>8.{prec}f} avg={:>8.{prec}f} std-dev={:8.{prec}f}\n", pszPrefix, minmax, us.max(), us.mean_value(), us.stddev(), "prec"_a = p);

    static constexpr size_t bin_width = 1; // Bin the histogram entries by this width
    // Find the index of the start of the all-zero tail
    auto rit = std::find_if(std::rbegin(us.histogram), std::rend(us.histogram), [](int x) { return x != 0; });
    auto idx = (rit == std::rend(us.histogram)) ? 0 : std::distance(rit, std::rend(us.histogram));
    // Output the histogram and truncate the all-zero tail. Output log2 scale (plus 1)
    fmt::print("      ", idx);
    for (int i = 0; i < idx; i += bin_width) {
        uint64_t s = std::accumulate(us.histogram + i, us.histogram + std::min(i + bin_width, size_t(idx)), uint64_t{0});
        s = s ? 64u - __builtin_clzll(s) : 0u; // log2(s)+1 or 0
        fmt::print(" {:>2}", s);
    }
    fmt::print("\n");
}
#endif

uint8_t MSR::Gravity(uint8_t uRungLo, uint8_t uRungHi, int iRoot1, int iRoot2, double dTime, double dDelta, double dStep, double dTheta, int bKickClose,
                     int bKickOpen, int bEwald, int bGravStep, int nPartRhoLoc, int iTimeStepCrit) {
    SPHOptions SPHoptions = initializeSPHOptions(parameters, csm, dTime, dDelta);
    SPHoptions.doGravity = parameters.get_bDoGravity();
    SPHoptions.nPredictRung = uRungLo;
    return Gravity(uRungLo, uRungHi, iRoot1, iRoot2, dTime, dDelta, dStep, dTheta, bKickClose, bKickOpen, bEwald, bGravStep, nPartRhoLoc, iTimeStepCrit,
                   SPHoptions);
}

uint8_t MSR::Gravity(uint8_t uRungLo, uint8_t uRungHi, int iRoot1, int iRoot2, double dTime, double dDelta, double dStep, double dTheta, int bKickClose,
                     int bKickOpen, int bEwald, int bGravStep, int nPartRhoLoc, int iTimeStepCrit, SPHOptions SPHoptions) {
    struct inGravity in;
    uint64_t nRungSum[IRUNGMAX + 1];
    int i;
    double dsec, dTotFlop, dt, a;
    double dTimeLCP;
    uint8_t uRungMax = 0;
    char c;

    a = csmTime2Exp(csm, dTime);

    if (parameters.get_bVStep()) {
        if (SPHoptions.doDensity && SPHoptions.useDensityFlags) print("Calculating Density using FastGas, Step:{:f} (rung {})", dStep, uRungLo);
        if (SPHoptions.doDensity && !SPHoptions.useDensityFlags) print("Calculating Density without FastGas, Step:{:f} (rung {})", dStep, uRungLo);
        if (SPHoptions.doDensityCorrection && SPHoptions.useDensityFlags)
            print("Calculating Density Correction using FastGas, Step:{:f} (rung {})", dStep, uRungLo);
        if (SPHoptions.doDensityCorrection && !SPHoptions.useDensityFlags)
            print("Calculating Density Correction without FastGas, Step:{:f} (rung {})", dStep, uRungLo);
        if (SPHoptions.doGravity && SPHoptions.doSPHForces) print("Calculating Gravity and SPH forces, Step:{:f} (rung {})", dStep, uRungLo);
        if (SPHoptions.doGravity && !SPHoptions.doSPHForces) print("Calculating Gravity, Step:{:f} (rung {})", dStep, uRungLo);
        if (!SPHoptions.doGravity && SPHoptions.doSPHForces) print("Calculating SPH forces, Step:{:f} (rung {})", dStep, uRungLo);
        if (SPHoptions.doSetDensityFlags) print("Marking Neighbors for FastGas, Step:{:f} (rung {})", dStep, uRungLo);
        if (SPHoptions.doSetNNflags) print("Marking Neighbors of Neighbors for FastGas, Step:{:f} (rung {})", dStep, uRungLo);
        if (csm->val.bComove) print(", Redshift:{:.6f}\n", 1. / a - 1.);
        else print(", Time:{:.6f}\n", dTime);
    }
    in.dTime = dTime;
    in.iRoot1 = iRoot1;
    in.iRoot2 = iRoot2;

    in.dTheta = dTheta;

    in.bPeriodic = parameters.get_bPeriodic();
    in.bEwald = bEwald;
    in.bGPU = parameters.get_bGPU();
    in.dEwCut = parameters.get_dEwCut();
    in.dEwhCut = parameters.get_dEwhCut();
    in.nReps = in.bPeriodic ? parameters.get_nReplicas() : 0;

    // Parameters related to timestepping
    in.ts.iTimeStepCrit = iTimeStepCrit;
    in.ts.nPartRhoLoc = nPartRhoLoc;
    in.ts.nPartColl = parameters.get_nPartColl();
    in.ts.dEccFacMax = parameters.get_dEccFacMax();
    if (bGravStep) in.ts.dRhoFac = 1.0 / (a * a * a);
    else in.ts.dRhoFac = 0.0;
    in.ts.dTime = dTime;
    in.ts.dDelta = dDelta;
    in.ts.dEta = Eta();
    in.ts.dPreFacRhoLoc = parameters.get_dPreFacRhoLoc();
    in.ts.bGravStep = bGravStep;
    in.ts.uRungLo = uRungLo;
    in.ts.uRungHi = uRungHi;
    in.ts.uMaxRung = parameters.get_iMaxRung();
    if (csm->val.bComove) in.ts.dAccFac = 1.0 / (a * a * a);
    else in.ts.dAccFac = 1.0;

    CalculateKickParameters(&in.kick, uRungLo, dTime, dDelta, dStep, bKickClose, bKickOpen, SPHoptions);

    in.lc.bLightConeParticles = parameters.get_bLightConeParticles();
    in.lc.dBoxSize = parameters.get_dBoxSize();
    if (parameters.get_bLightCone()) {
        in.lc.dLookbackFac = csmComoveKickFac(csm, dTime, (csmExp2Time(csm, 1.0) - dTime));
        dTimeLCP = csmExp2Time(csm, 1.0 / (1.0 + parameters.get_dRedshiftLCP()));
        in.lc.dLookbackFacLCP = csmComoveKickFac(csm, dTimeLCP, (csmExp2Time(csm, 1.0) - dTimeLCP));
        auto sqdegLCP = parameters.get_sqdegLCP();
        if (sqdegLCP <= 0 || sqdegLCP >= 4 * M_1_PI * 180.0 * 180.0) {
            in.lc.tanalpha_2 = -1; // indicates we want an all sky lightcone, a bit weird but this is the flag.
        }
        else {
            double alpha = sqrt(sqdegLCP * M_1_PI) * (M_PI / 180.0);
            in.lc.tanalpha_2 = tan(0.5 * alpha); // it is tangent of the half angle that we actually need!
        }
    }
    else {
        in.lc.dLookbackFac = 0.0;
        in.lc.dLookbackFacLCP = 0.0;
    }
    in.lc.hLCP = parameters.get_hLCP();

    /*
    ** Note that in this loop we initialize dt with the full step, not a half step!
    */
    for (i = 0, dt = dDelta; i <= parameters.get_iMaxRung(); ++i, dt *= 0.5) {
        in.lc.dtLCDrift[i] = 0.0;
        in.lc.dtLCKick[i] = 0.0;
        if (i >= uRungLo) {
            if (csm->val.bComove) {
                in.lc.dtLCDrift[i] = csmComoveDriftFac(csm, dTime, dt);
                in.lc.dtLCKick[i] = csmComoveKickFac(csm, dTime, dt);
            }
            else {
                in.lc.dtLCDrift[i] = dt;
                in.lc.dtLCKick[i] = dt;
            }
        }
    }

    in.SPHoptions = SPHoptions;

    outGravityReduct outr;

    TimerStart(TIMER_GRAVITY);

    pstGravity(pst, &in, sizeof(in), &outr, sizeof(outr));

    TimerStop(TIMER_GRAVITY);
    dsec = TimerGet(TIMER_GRAVITY);

    if (bKickOpen) {
        for (i = IRUNGMAX; i >= uRungLo; --i) {
            if (outr.nRung[i]) break;
        }
        assert(i >= uRungLo);
        uRungMax = i;
        iCurrMaxRung = uRungMax; /* this assignment shouldn't be needed */
        /*
        ** Update only the active rung counts in the master rung counts.
        ** We need to go all the way to IRUNGMAX to clear any prior counts at rungs
        ** deeper than the current uRungMax!
        */
        for (i = uRungLo; i <= IRUNGMAX; ++i) nRung[i] = outr.nRung[i];

        const uint64_t nDT = (N * parameters.get_dFracDualTree());
        const uint64_t nDD = (N * parameters.get_dFracNoDomainDecomp());
        uint64_t nActive = 0;
        iRungDD = 0;
        iRungDT = 0;
        for (i = iCurrMaxRung; i >= 0; --i) {
            nActive += nRung[i];
            if (nActive > nDT && !iRungDT) iRungDT = i;
            if (nActive > nDD && !iRungDD) iRungDD = i;
        }
    }
    if (parameters.get_bVStep()) {
#ifdef INSTRUMENT_MPI
        ServiceCacheStats::output stats;
        auto n = mdl->RunService<ServiceCacheStats>(stats);
#endif
        /*
        ** Output some info...
        */
        dTotFlop = outr.sFlop.dSum;
        if (dsec > 0.0) {
            double dGFlops = dTotFlop / dsec;
            print("Gravity Calculated, Wallclock: {:.5} secs, Gflops:{:.1f}, Total Gflop:{:.3g}\n", dsec, dGFlops, dTotFlop);
            print("  Gflops: CPU:{:.1f},{:.1f} GPU:{:.1f},{:.1f}", outr.dFlopSingleCPU / dsec, outr.dFlopDoubleCPU / dsec, outr.dFlopSingleGPU / dsec,
                  outr.dFlopDoubleGPU / dsec);
        }
        else {
            print("Gravity Calculated, Wallclock: {:.5} secs, Gflops:unknown, Total Gflop:{:.3g}\n", dsec, dTotFlop);
        }
        print("  Gflop: CPU:{:.3g},{:.3g} GPU:{:.3g},{:.3g}\n", outr.dFlopSingleCPU, outr.dFlopDoubleCPU, outr.dFlopSingleGPU, outr.dFlopDoubleGPU);
        PrintStat(outr.sLocal, "  particle  load:", 0);
        PrintStat(outr.sActive, "  actives   load:", 0);
        PrintStat(outr.sFlop, "  Gflop     load:", 1);
        PrintStat(outr.sPart, "  P-P per active:", 2);
        PrintStat(outr.sCell, "  P-C per active:", 2);
#ifdef INSTRUMENT
        PrintStat(outr.sComputing, "     % computing:", 3);
        PrintStat(outr.sWaiting, "     %   waiting:", 3);
        PrintStat(outr.sSynchronizing, "     %   syncing:", 3);
#endif
#ifdef __linux__
        PrintStat(outr.sFreeMemory, "free memory (GB):", 3);
        PrintStat(outr.sRSS, "   resident size:", 3);
#endif
        print("  (cache access statistics are given per active particle)\n");
        PrintStat(outr.sPartNumAccess, "  P-cache access:", 1);
        PrintStat(outr.sCellNumAccess, "  C-cache access:", 1);
        PrintStat(outr.sPartMissRatio, "  P-cache miss %:", 2);
        PrintStat(outr.sCellMissRatio, "  C-cache miss %:", 2);
#ifdef INSTRUMENT_MPI
        print("  (cache timeing statistics are given in microseconds)\n");
        PrintStat(stats.mpi_loop, "   MPI Loop time:", 3);
        PrintStat(stats.mpi_test, "   MPI Test time:", 3);
        PrintStat(stats.gpu_test, "   GPU Test time:", 3);
        PrintStat(stats.msg_test, "   MSG Test time:", 3);
#endif
    }
    if (outr.nTilesTotal > 0) {
        print("Total tiles processed: {:.5e}, on the GPU: {:.5e}, ratio: {:.2f} %\n", (double)outr.nTilesTotal, (double)(outr.nTilesTotal - outr.nTilesCPU),
              100.0 - ((double)outr.nTilesCPU) / ((double)outr.nTilesTotal) * 100.0);
    }
    if (parameters.get_bVRungStat() && bKickOpen) {
        print("Rung distribution:\n\n");
        nRungSum[uRungMax] = nRung[uRungMax];
        for (i = uRungMax - 1; i >= 0; --i) {
            nRungSum[i] = nRungSum[i + 1] + nRung[i];
        }
        for (i = 0; i <= uRungMax; ++i) {
            if (nRung[i]) break;
        }
        if (nRungSum[0] > 0)
            for (; i <= uRungMax; ++i) {
                c = ' ';
                print(" {:c} rung:{} {:14} {:14} {:3.0f} %\n", c, i, nRung[i], nRungSum[i], std::ceil(100.0 * nRungSum[i] / nRungSum[0]));
            }
        print("\n");
    }
    return (uRungMax);
}

void MSR::CalcEandL(int bFirst, double dTime, double& E, double& T, double& U, double& Eth, blitz::TinyVector<double, 3>& L, blitz::TinyVector<double, 3>& F,
                    double& W) {
    ServiceCalcEandL::output out;

    mdl->RunService<ServiceCalcEandL>({}, out);
    T = out.T;
    U = out.U;
    Eth = out.Eth;
    L = out.L;
    F = out.F;
    W = out.W;
    /*
    ** Do the comoving coordinates stuff.
    ** Currently L is not adjusted for this. Should it be?
    */
    auto a = csmTime2Exp(csm, dTime);
    if (!csm->val.bComove) T *= pow(a, 4.0);
    /*
     * Estimate integral (\dot a*U*dt) over the interval.
     * Note that this is equal to integral (W*da) and the latter
     * is more accurate when a is changing rapidly.
     */
    if (csm->val.bComove && !bFirst) {
        dEcosmo += 0.5 * (a - csmTime2Exp(csm, dTimeOld)) * (U + dUOld);
    }
    else {
        dEcosmo = 0.0;
    }
    dTimeOld = dTime;
    dUOld = U;
    U *= a;
    E = T + U - dEcosmo + a * a * Eth;
}

void MSR::Drift(double dTime, double dDelta, int iRoot) {
    double dsec;

    if (parameters.get_bBlackholes() && parameters.get_bBHPin()) {
        TimerStart(TIMER_BHS);
        BHGasPin(dTime, dDelta);
        TimerStop(TIMER_BHS);
    }

    TimerStart(TIMER_DRIFT);

    ServiceDrift::input in{.dTime = dTime,
                           .dDelta = csm->val.bComove ? csmComoveDriftFac(csm, dTime, dDelta) : dDelta,
                           .dDeltaVPred = csm->val.bComove ? csmComoveKickFac(csm, dTime, dDelta) : dDelta,
                           .dDeltaUPred = dDelta,
                           .iRoot = iRoot,
                           .bDoGas = DoGas(),
                           .bGasEvolveDensity = parameters.get_bGasEvolveDensity()};
    mdl->RunService<ServiceDrift>(in);

    TimerStop(TIMER_DRIFT);
    dsec = TimerGet(TIMER_DRIFT);
    print("Drift took {:.5f} seconds \n", dsec);

    if (parameters.get_bBlackholes() && parameters.get_bBHPin()) {
        TimerStart(TIMER_BHS);
        BHReposition();
        TimerStop(TIMER_BHS);
    }
}

void MSR::OutputFineStatistics(double dStep, double dTime) {
    if (!parameters.get_bOutFineStatistics()) return;
    if (dTime == -1) {
        std::string achFile(OutName());
        achFile += ".finelog";
        fpFineLog = fopen(achFile.c_str(), "a");
        assert(fpFineLog != NULL);
        setvbuf(fpFineLog, NULL, _IONBF, 0); // Disable buffering

        /* Write the header */
        fprintf(fpFineLog, "#First line:\n#dStep dTime ");
        fprintf(fpFineLog, "starsFormed massFormed ");
        fprintf(fpFineLog, "\n");

        fprintf(fpFineLog, "#Second line:\n#nRung[0] nRung[1] nRung[...] nRung[iCurrMaxRung]\n");
    }
    else {
        /* First, we add a line with relevant statistics */
        fprintf(fpFineLog, "%e %e ", dStep, dTime);
        fprintf(fpFineLog, "%d %e ", starFormed, massFormed);
        fprintf(fpFineLog, "\n");

        /* Second, we add the rung distribution */
        for (int i = 0; i < iCurrMaxRung; i++) fprintf(fpFineLog, "%" PRIu64 " ", nRung[i]);
        fprintf(fpFineLog, "\n");
    }
}

void MSR::EndTimestepIntegration(double dTime, double dDelta) {
    ServiceEndTimestepIntegration::input in;
    in.units = UNITS(parameters, csm->val.h);
#ifdef GRACKLE
    marshal_string(in.achCoolingTable, parameters.get_achCoolingTables());
#endif
    in.dTime = dTime;
    in.dDelta = dDelta;
    in.dConstGamma = parameters.get_dConstGamma();
    in.dTuFac = dTuFacPrimNeutral;
    in.dSFThresholdOD = calc.dSFThresholdOD;
    in.eEOS = get_eEOS();
    in.dBHRadiativeEff = parameters.get_dBHRadiativeEff();
#ifdef STELLAR_EVOLUTION
    in.bChemEnrich = parameters.get_bChemEnrich();
#endif
    double dsec;

    ComputeSmoothing(dTime, dDelta);

    print("Computing primitive variables... ");
    TimerStart(TIMER_ENDINT);
    mdl->RunService<ServiceEndTimestepIntegration>(in);

    TimerStop(TIMER_ENDINT);
    dsec = TimerGet(TIMER_ENDINT);
    print("took {:.5f} seconds\n", dsec);
}

/*
 * For gas, updates predicted velocities to beginning of timestep.
 */
void MSR::KickKDKOpen(double dTime, double dDelta, uint8_t uRungLo, uint8_t uRungHi) {
    TimerStart(TIMER_KICKO);

    auto dDeltaKick = csm->val.bComove ? csmComoveKickFac(csm, dTime, dDelta) : dDelta;
    ServiceKick::input in{.dTime = dTime,
                          .dDelta = dDeltaKick,
                          .dDeltaVPred = 0.0,
                          .dDeltaU = dDelta,
                          .dDeltaUPred = 0.0,
                          .uRungLo = uRungLo,
                          .uRungHi = uRungHi,
                          .bDoGas = DoGas()};
    mdl->RunService<ServiceKick>(in);
    TimerStop(TIMER_KICKO);
}

/*
 * For gas, updates predicted velocities to end of timestep.
 */
void MSR::KickKDKClose(double dTime, double dDelta, uint8_t uRungLo, uint8_t uRungHi) {
    TimerStart(TIMER_KICKC);
    auto dDeltaKick = csm->val.bComove ? csmComoveKickFac(csm, dTime, dDelta) : dDelta;
    ServiceKick::input in{.dTime = dTime,
                          .dDelta = dDeltaKick,
                          .dDeltaVPred = dDeltaKick,
                          .dDeltaU = dDelta,
                          .dDeltaUPred = dDelta,
                          .uRungLo = uRungLo,
                          .uRungHi = uRungHi,
                          .bDoGas = DoGas()};
    mdl->RunService<ServiceKick>(in);
    TimerStop(TIMER_KICKC);
}

int cmpTime(const void* v1, const void* v2) {
    double* d1 = (double*)v1;
    double* d2 = (double*)v2;

    if (*d1 < *d2) return (-1);
    else if (*d1 == *d2) return (0);
    else return (1);
}

void MSR::InitCosmology(CSM csm) {
    ServiceInitLightcone::input in;
    if (parameters.has_h()) csm->val.h = parameters.get_h();
    mdl->RunService<ServiceInitCosmology>(csm->val);
    if (parameters.get_bLightCone()) {
        auto sqdegLCP = parameters.get_sqdegLCP();
        if (sqdegLCP <= 0 || sqdegLCP >= 4 * M_1_PI * 180.0 * 180.0) {
            in.alphaLCP = -1; // indicates we want an all sky lightcone, a bit weird but this is the flag.
        }
        else {
            in.alphaLCP = sqrt(sqdegLCP * M_1_PI) * (M_PI / 180.0);
        }
        in.bBowtie = parameters.get_bBowtie();
        in.bLightConeParticles = parameters.get_bLightConeParticles();
        in.dBoxSize = parameters.get_dBoxSize();
        in.dRedshiftLCP = parameters.get_dRedshiftLCP();
        in.hLCP = parameters.get_hLCP();
        in.mOctants = 0;
        if (parameters.get_bLightConeParticles()) {
            if (parameters.has_lstLightConeOctants()) {
                auto octants = parameters.get_lstLightConeOctants();
                for (auto o : octants) in.mOctants |= (1 << o);
            }
        }
        mdl->RunService<ServiceInitLightcone>(in);
    }
}

void MSR::ZeroNewRung(uint8_t uRungLo, uint8_t uRungHi, int uRung) {
    if (parameters.get_bMemNewRung()) {
        ServiceZeroNewRung::input in(uRung, uRungLo, uRungHi);
        mdl->RunService<ServiceZeroNewRung>(in);
    }
}

/*
 * bGreater = 1 => activate all particles at this rung and greater.
 */
void MSR::ActiveRung(int iRung, int bGreater) {
    ServiceActiveRung::input in(iRung, bGreater);
    mdl->RunService<ServiceActiveRung>(in);

    if (iRung == 0 && bGreater) nActive = N;
    else {
        int i;

        nActive = 0;
        for (i = iRung; i <= (bGreater ? parameters.get_iMaxRung() : iRung); i++) nActive += nRung[i];
    }
}

auto MSR::ActiveOrder() -> uint64_t {
    mdl->RunService<OldDD::ServiceActiveOrder>({}, nActive);
    return nActive;
}

int MSR::CountRungs(uint64_t* nRungs) {
    ServiceCountRungs::output out;
    int i, iMaxRung = 0;
    mdl->RunService<ServiceCountRungs>({}, out);
    for (i = 0; i <= MAX_RUNG; ++i) {
        nRung[i] = out[i];
        if (nRung[i]) iMaxRung = i;
        if (nRungs) nRungs[i] = nRung[i];
    }
    iCurrMaxRung = iMaxRung;

    const uint64_t nDT = (N * parameters.get_dFracDualTree());
    const uint64_t nDD = (N * parameters.get_dFracNoDomainDecomp());
    uint64_t nActive = 0;
    iRungDD = 0;
    iRungDT = 0;
    for (i = iCurrMaxRung; i >= 0; --i) {
        nActive += nRung[i];
        if (nActive > nDT && !iRungDT) iRungDT = i;
        if (nActive > nDD && !iRungDD) iRungDD = i;
    }

    return iMaxRung;
}

void MSR::AccelStep(uint8_t uRungLo, uint8_t uRungHi, double dTime, double dDelta) {
    auto a = csmTime2Exp(csm, dTime);
    auto dVelFac = csm->val.bComove ? 1.0 / (a * a) : 1.0;
    ServiceAccelStep::input in{.dEta = Eta(),
                               .dVelFac = dVelFac,
                               .dAccFac = 1.0 / (a * a * a),
                               .dDelta = dDelta,
                               .dhMinOverSoft = parameters.get_dhMinOverSoft(),
                               .iMaxRung = static_cast<int>(parameters.get_iMaxRung()),
                               .uRungLo = uRungLo,
                               .uRungHi = uRungHi,
                               .bDoGravity = DoGravity(),
                               .bEpsAcc = parameters.get_bEpsAccStep()};
    mdl->RunService<ServiceAccelStep>(in);
}

auto MSR::GetMaxRung() -> uint8_t {
    ServiceMaxRung::output out;
    mdl->RunService<ServiceMaxRung>({}, out);
    return out.uMaxRung;
}

void MSR::SetGlobalRung(uint8_t uRung) {
    ServiceSetGlobalRung::input in{.uMaxRung = uRung};
    mdl->RunService<ServiceSetGlobalRung>(in);
}

void MSR::DensityStep(uint8_t uRungLo, uint8_t uRungHi, double dTime, double dDelta) {

    print_detail("Calculating Rung Densities...\n");
    int bSymmetric = 0;
    Smooth(dTime, dDelta, SMX_DENSITY, bSymmetric, parameters.get_nSmooth());
    auto expand = csmTime2Exp(csm, dTime);

    ServiceDensityStep::input in{.dDelta = dDelta,
                                 .dEta = Eta(),
                                 .dRhoFac = 1.0 / (expand * expand * expand),
                                 .iMaxRung = static_cast<int>(parameters.get_iMaxRung()),
                                 .uRungLo = uRungLo,
                                 .uRungHi = uRungHi};
    mdl->RunService<ServiceDensityStep>(in);
}

/*
 ** Returns the Very Active rung based on the number of very active particles desired,
 ** or the fixed rung that was specified in the parameters.
 */
void MSR::UpdateRung(uint8_t uRung) {
    ServiceUpdateRung::input in;
    ServiceUpdateRung::output nRungCount;
    int iTempRung, iOutMaxRung;

    /* If we are called, it is a mistake -- this happens in analysis mode */
    if (!parameters.get_bMemNewRung()) return;

    in.uRungLo = uRung;
    in.uRungHi = MaxRung();
    in.uMinRung = uRung;
    in.uMaxRung = MaxRung();

    mdl->RunService<ServiceUpdateRung>(in, nRungCount);

    iTempRung = MaxRung() - 1;
    while (nRungCount[iTempRung] == 0 && iTempRung > 0) --iTempRung;
    iOutMaxRung = iTempRung;

    const auto nTruncateRung = parameters.get_nTruncateRung();
    while (nRungCount[iOutMaxRung] <= nTruncateRung && iOutMaxRung > uRung) {
        // print_detail("n_CurrMaxRung = {nCurrMaxRung}  (iCurrMaxRung = {iCurrMaxRung}):  Promoting particles to iCurrMaxrung = {iCurrMaxRung-1}\n");
        print_detail("n_CurrMaxRung = {} (iCurrMaxRung = {}):  Promoting particles to iCurrMaxrung = {}\n", nRungCount[iOutMaxRung], iOutMaxRung,
                     iOutMaxRung - 1);

        in.uMaxRung = iOutMaxRung; /* Note this is the forbidden rung so no -1 here */
        mdl->RunService<ServiceUpdateRung>(in, nRungCount);

        iTempRung = MaxRung() - 1;
        while (nRungCount[iTempRung] == 0 && iTempRung > 0) --iTempRung;
        iOutMaxRung = iTempRung;
    }

    /*
    ** Now copy the rung distribution to the msr structure!
    */
    for (iTempRung = 0; iTempRung < MaxRung(); ++iTempRung) nRung[iTempRung] = nRungCount[iTempRung];

    iCurrMaxRung = iOutMaxRung;

    if (parameters.get_bVRungStat()) {
        print("Rung distribution:\n\n");
        for (iTempRung = 0; iTempRung <= iCurrMaxRung; ++iTempRung) {
            if (nRungCount[iTempRung] == 0) continue;
            print("   rung:{} {}\n", iTempRung, nRungCount[iTempRung]);
        }
        print("\n");
    }
}

/*
 ** Open the healpix output file, and also the particles files if requested.
 */
void MSR::LightConeOpen(int iStep) {
    if (parameters.get_bLightCone()) {
        ServiceLightConeOpen::input lc;
        if (parameters.get_bLightConeParticles()) {
            auto filename = [&]() {
                if (parameters.has_lstLightConeOctants()) {
                    auto octants = parameters.get_lstLightConeOctants();
                    return BuildNameWithOI(iStep, ".lcp");
                }
                else return BuildNameWithI(iStep, ".lcp");
            }();
            marshal_string(lc.achOutFileLcp, filename);
        }
        else lc.achOutFileLcp[0] = 0;
        lc.nSideHealpix = parameters.get_nSideHealpix();
        mdl->RunService<ServiceLightConeOpen>(lc);
    }
}

/*
 ** Close the files for this step.
 */
void MSR::LightConeClose(int iStep) {
    if (parameters.get_bLightCone()) {
        ServiceLightConeClose::input lc;
        // BuildName actually creates directories, so only do so if we are making the maps
        if (parameters.get_nSideHealpix()) {
            auto filename = BuildNameWithI(iStep, ".hpb");
            marshal_string(lc.achOutFileHpb, filename);
        }
        else lc.achOutFileHpb[0] = 0;
        mdl->RunService<ServiceLightConeClose>(lc);
    }
}

/*
 ** Correct velocities from a^2 x_dot to a x_dot (physical peculiar velocities) using the
 ** position dependent scale factor within the light cone. This could be expensive.
 */
void MSR::LightConeVel() {
    double dsec;
    ServiceLightConeVel::input in;

    in.dBoxSize = parameters.get_dBoxSize();
    TimerStart(TIMER_NONE);
    mdl->RunService<ServiceLightConeVel>(in);
    TimerStop(TIMER_NONE);
    dsec = TimerGet(TIMER_NONE);
    print("Converted lightcone velocities to physical, Wallclock: {:.5f} secs.\n", dsec);
}

/* True if we should omit the opening kick */
int MSR::CheckForOutput(int iStep, int nSteps, double dTime, int* pbDoCheckpoint, int* pbDoOutput) {
    int iStop, iCheck;
    long lSec = time(0) - lPrior;

    /*
    ** Check for user interrupt.
    */
    iStop = CheckForStop(STOPFILE);
    iCheck = CheckForStop(CHECKFILE);

    /*
    ** Check to see if the runtime has been exceeded.
    */
    const auto iWallRunTime = parameters.get_iWallRunTime();
    if (!iStop && iWallRunTime > 0) {
        if (iWallRunTime * 60 - (time(0) - lStart) < ((int)(lSec * 1.5))) {
            print("RunTime limit exceeded. Writing checkpoint and exiting.\n"
                  "    iWallRunTime(sec): {}   Time running: {}   Last step: {}\n",
                  iWallRunTime * 60, time(0) - lStart, lSec);
            iStop = 1;
        }
    }

    /* Check to see if there should be an output */
    const auto iSignalSeconds = parameters.get_iSignalSeconds();
    if (!iStop && timeGlobalSignalTime > 0) { /* USR1 received */
        if ((time(0) + (lSec * 1.5)) > timeGlobalSignalTime + iSignalSeconds) {
            print("RunTime limit exceeded.  Writing checkpoint and exiting.\n"
                  "    iSignalSeconds: {}   Time running: {}   Last step: {}\n",
                  iSignalSeconds, time(0) - lStart, lSec);
            iStop = 1;
        }
    }

    /*
    ** Output if 1) we've hit an output time
    **           2) We are stopping
    **           3) we're at an output interval
    */
    if (iCheck || (CheckInterval() > 0 && (bGlobalOutput || iStop || (iStep % CheckInterval() == 0)))) {
        bGlobalOutput = 0;
        *pbDoCheckpoint = 1 | (iStop << 1);
    }

    if ((OutInterval() > 0 && (bGlobalOutput || iStep == iStep_list.front() || iStop || iStep == nSteps || (iStep % OutInterval() == 0)))) {
        bGlobalOutput = 0;
        *pbDoOutput = 1 | (iStop << 1);
    }

    return *pbDoOutput || *pbDoCheckpoint;
}

void MSR::StellarFeedback(double dTime) {
    ServiceSNFeedback::input in{.nSmooth = static_cast<int>(parameters.get_nSmooth()),
                                .bPeriodic = parameters.get_bPeriodic(),
                                .dTime = dTime,
                                .dCCSNFBDelay = parameters.get_dCCSNFBDelay(),
                                .dSNIaFBDelay = parameters.get_dSNIaFBDelay(),
                                .dSNFBDu = calc.dSNFBDu,
                                .dCCSNFBSpecEnergy = calc.dCCSNFBSpecEnergy,
                                .dSNIaFBSpecEnergy = calc.dSNIaFBSpecEnergy,
                                .dConstGamma = parameters.get_dConstGamma()};
    mdl->RunService<ServiceSNFeedback>(in);
}

void MSR::ChemEnrich(double dTime) {
    ServiceChemEnrich::input in{.nSmooth = static_cast<int>(parameters.get_nSmooth()), .bPeriodic = parameters.get_bPeriodic(), .dTime = dTime};
    mdl->RunService<ServiceChemEnrich>(in);
}

int MSR::NewTopStepKDK(double& dTime,                                           /* MODIFIED: Current simulation time */
                       double dDelta, double dTheta, int nSteps, int bDualTree, /* Should be zero at rung 0! */
                       uint8_t uRung,                                           /* Rung level */
                       double* pdStep,                                          /* Current step */
                       uint8_t* puRungMax, int* pbDoCheckpoint, int* pbDoOutput, int* pbNeedKickOpen) {
    const auto bEwald = parameters.get_bEwald();
    const auto bGravStep = parameters.get_bGravStep();
    const auto nPartRhoLoc = parameters.get_nPartRhoLoc();
    const auto iTimeStepCrit = parameters.get_iTimeStepCrit();
    double dDeltaRung, dTimeFixed;
    uint32_t uRoot2 = 0;
    int bKickOpen = 1;
    /*
    ** The iStep variable serves only to give a number to the lightcone and group output files.
    ** We define this to be the output number of the final radius of the lightcone surface.
    */
    int iStep = (int)(*pdStep) + 1;

    /* IA: If the next rung that is the first one with actives particles, we
     *   are sure that all particles are synchronized, thus we can output some statistics with
     *   finer time resolution whilst being accurate.
     *
     *   We assume that there is no 'sandwiched' rung, i.e., the differences in dt are smooth
     */
    if ((nRung[0] != 0 && uRung == 0) || ((nRung[uRung] == 0) && (nRung[uRung + 1] > 0))) OutputFineStatistics(*pdStep, dTime);

    if (uRung == iRungDT + 1) {
        if (parameters.get_bDualTree() && uRung < *puRungMax) {
            /* HACK: FIXME: Don't use the dual tree before z=2; the overlap region is too large */
            /* better would be to construct the tree matching remote processor shape as well as local */
            double a = csmTime2Exp(csm, dTime);
            if (a < (1.0 / 3.0)) bDualTree = 0;
            else {
                bDualTree = 1;
                ServiceDumpTrees::input dump(iRungDT);
                mdl->RunService<ServiceDumpTrees>(dump);
                print_detail("Half Drift, uRung: {rung}\n", "rung"_a = iRungDT);
                dDeltaRung = dDelta / (uintmax_t(1) << iRungDT); // Main tree step
                Drift(dTime, 0.5 * dDeltaRung, FIXROOT);
                dTimeFixed = dTime + 0.5 * dDeltaRung;
                BuildTreeFixed(bEwald, iRungDT);
            }
        }
        else bDualTree = 0;
    }
    if (uRung < *puRungMax) {
        bDualTree = NewTopStepKDK(dTime, dDelta, dTheta, nSteps, bDualTree, uRung + 1, pdStep, puRungMax, pbDoCheckpoint, pbDoOutput, pbNeedKickOpen);
    }

    dDeltaRung = dDelta / (uintmax_t(1) << *puRungMax);
    ActiveRung(uRung, 1);
    if (DoGas() && MeshlessHydro()) {
        MeshlessFluxes(dTime, dDelta);
    }
    ZeroNewRung(uRung, MAX_RUNG, uRung);

    if (parameters.get_bBlackholes() && parameters.get_bBHPlaceSeed()) {
        PlaceBHSeed(dTime, *puRungMax);
    }
    /* Drift the "ROOT" (active) tree or all particle */
    if (bDualTree) {
        print_detail("Drift very actives, uRung: {rung}\n", "rung"_a = *puRungMax);
        Drift(dTime, dDeltaRung, ROOT);
    }
    else {
        print_detail("Drift, uRung: {rung}\n", "rung"_a = *puRungMax);
        Drift(dTime, dDeltaRung, -1);
    }
    dTime += dDeltaRung;
    *pdStep += 1.0 / (uintmax_t(1) << *puRungMax);
    if (parameters.get_bDoCooling()) {
        if (csm->val.bComove) {
            const float a = csmTime2Exp(csm, dTime);
            const float z = 1. / a - 1.;

            CoolingUpdate(z);
        }
        else {
            CoolingUpdate(0.);
        }
    }
    if (parameters.get_bDoStarFormation()) StarForm(dTime, dDelta, uRung);

    ActiveRung(uRung, 1);
    UpdateSoft(dTime);
    if (bDualTree && uRung > iRungDT) {
        uRoot2 = FIXROOT;
        BuildTreeActive(bEwald, iRungDT);
    }
    else {
        DomainDecomp(uRung);
        uRoot2 = 0;

        if (NewSPH()) {
            SelAll(-1, 1);
        }

        if (parameters.get_bBlackholes()) {
            auto bBHMerger = parameters.get_bBHMerger();
            if (bBHMerger) {
                SelActives();
                BHMerger(dTime);
            }
            if (parameters.get_bBHAccretion() && !bBHMerger) {
                RemoveDeleted();
            }
        }

        BuildTree(bEwald);
    }

    if (!uRung) {
        bKickOpen = !CheckForOutput(iStep, nSteps, dTime, pbDoCheckpoint, pbDoOutput);
    }
    else bKickOpen = 1;
    *pbNeedKickOpen = !bKickOpen;

    /*
    ** At this point the particles are in sync. As soon as we call the next gravity it will
    ** advance the particles as part of the opening kick. For this reason we do the requested
    ** analysis at this poing while everything is properly synchronized. This includes
    ** writing the healpix and lightcone particles, as well as measuring P(k) for example.
    */
    if (uRung == 0) {
        runAnalysis(iStep, dTime); // Run any registered Python analysis tasks

        const auto iPkInterval = parameters.get_iPkInterval();
        if (iPkInterval && iStep % iPkInterval == 0) OutputPk(iStep, dTime);

        /*
        ** We need to write all light cone files (healpix and LCP) at this point before the last
        ** gravity is called since it will advance the particles in the light cone as part of the
        ** opening kick! We also need to open
        */
        LightConeClose(iStep);
        if (bKickOpen) LightConeOpen(iStep + 1);

        /* Compute the grids of linear species at main timesteps, before gravity is called */
        auto nGridLin = parameters.get_nGridLin();
        if (csm->val.classData.bClass && parameters.get_achLinSpecies().length() && nGridLin) {
            GridCreateFFT(nGridLin);
            SetLinearGrid(dTime, dDelta, nGridLin, 1, bKickOpen);
            if (parameters.get_bDoLinPkOutput()) OutputLinearPk(*pdStep, dTime);
            LinearKick(dTime, dDelta, 1, bKickOpen);
            GridDeleteFFT();
        }

        if (parameters.get_bFindGroups()) NewFof(parameters.get_dTau(), parameters.get_nMinMembers());
    }

    // We need to make sure we descend all the way to the bucket with the
    // active tree, or we can get HUGE group cells, and hence too much P-P/P-C
    if (NewSPH()) {
        SelAll(-1, 1);
        SPHOptions SPHoptions = initializeSPHOptions(parameters, csm, dTime, dDelta);
        uint64_t nParticlesOnRung = 0;
        for (int i = MaxRung(); i >= uRung; i--) {
            nParticlesOnRung += nRung[i];
        }
        if (nParticlesOnRung / ((float)N) < SPHoptions.FastGasFraction) {
            // Select Neighbors
            SPHoptions.doGravity = 0;
            SPHoptions.doDensity = 0;
            SPHoptions.nPredictRung = uRung;
            SPHoptions.doSPHForces = 0;
            SPHoptions.doSetDensityFlags = 1;
            *puRungMax = Gravity(uRung, MaxRung(), ROOT, uRoot2, dTime, dDelta, *pdStep, dTheta, 1, bKickOpen, bEwald, bGravStep, nPartRhoLoc, iTimeStepCrit,
                                 SPHoptions);
            // Select Neighbors of Neighbors
            if (SPHoptions.doInterfaceCorrection) {
                SPHoptions.dofBallFactor = 1;
                TreeUpdateFlagBounds(bEwald, ROOT, 0, SPHoptions);
                SPHoptions.doSetDensityFlags = 0;
                SPHoptions.doSetNNflags = 1;
                SPHoptions.useDensityFlags = 1;
                *puRungMax = Gravity(uRung, MaxRung(), ROOT, uRoot2, dTime, dDelta, *pdStep, dTheta, 1, bKickOpen, bEwald, bGravStep, nPartRhoLoc,
                                     iTimeStepCrit, SPHoptions);
                SPHoptions.doSetNNflags = 0;
                SPHoptions.useDensityFlags = 0;
            }
        }
        // Calculate Density
        SPHoptions.doSetDensityFlags = 0;
        SPHoptions.doGravity = 0;
        SPHoptions.doDensity = 1;
        SPHoptions.nPredictRung = uRung;
        SPHoptions.doSPHForces = 0;
        SPHoptions.useDensityFlags = 0;
        if (nParticlesOnRung / ((float)N) < SPHoptions.FastGasFraction) {
            SPHoptions.useDensityFlags = 1;
            SPHoptions.dofBallFactor = 1;
            if (SPHoptions.doInterfaceCorrection) {
                SPHoptions.useNNflags = 1;
            }
            TreeUpdateFlagBounds(bEwald, ROOT, 0, SPHoptions);
            *puRungMax = Gravity(uRung, MaxRung(), ROOT, uRoot2, dTime, dDelta, *pdStep, dTheta, 1, bKickOpen, bEwald, bGravStep, nPartRhoLoc, iTimeStepCrit,
                                 SPHoptions);
            if (SPHoptions.doInterfaceCorrection) {
                SPHoptions.doDensity = 0;
                SPHoptions.doDensityCorrection = 1;
                SPHoptions.useDensityFlags = 1;
                SPHoptions.useNNflags = 0;
                SPHoptions.dofBallFactor = 0;
                TreeUpdateFlagBounds(bEwald, ROOT, 0, SPHoptions);
                *puRungMax = Gravity(uRung, MaxRung(), ROOT, uRoot2, dTime, dDelta, *pdStep, dTheta, 1, bKickOpen, bEwald, bGravStep, nPartRhoLoc,
                                     iTimeStepCrit, SPHoptions);
                UpdateGasValues(uRung, dTime, dDelta, *pdStep, 1, bKickOpen, SPHoptions);
                SPHoptions.doDensityCorrection = 0;
                SPHoptions.useDensityFlags = 0;
            }
        }
        else {
            *puRungMax =
                Gravity(0, MaxRung(), ROOT, uRoot2, dTime, dDelta, *pdStep, dTheta, 1, bKickOpen, bEwald, bGravStep, nPartRhoLoc, iTimeStepCrit, SPHoptions);
            if (SPHoptions.doInterfaceCorrection) {
                SPHoptions.doDensity = 0;
                SPHoptions.doDensityCorrection = 1;
                SPHoptions.dofBallFactor = 0;
                TreeUpdateFlagBounds(bEwald, ROOT, 0, SPHoptions);
                *puRungMax = Gravity(0, MaxRung(), ROOT, uRoot2, dTime, dDelta, *pdStep, dTheta, 1, bKickOpen, bEwald, bGravStep, nPartRhoLoc, iTimeStepCrit,
                                     SPHoptions);
                UpdateGasValues(0, dTime, dDelta, *pdStep, 1, bKickOpen, SPHoptions);
                SPHoptions.doDensityCorrection = 0;
            }
        }
        // Calculate Forces
        SelAll(-1, 1);
        SPHoptions.doGravity = parameters.get_bDoGravity();
        SPHoptions.doDensity = 0;
        SPHoptions.nPredictRung = uRung;
        SPHoptions.doSPHForces = 1;
        SPHoptions.useDensityFlags = 0;
        SPHoptions.dofBallFactor = 0;
        TreeUpdateFlagBounds(bEwald, ROOT, 0, SPHoptions);
        *puRungMax =
            Gravity(uRung, MaxRung(), ROOT, uRoot2, dTime, dDelta, *pdStep, dTheta, 1, bKickOpen, bEwald, bGravStep, nPartRhoLoc, iTimeStepCrit, SPHoptions);
    }
    else if (parameters.get_bDoGravity()) {
        SPHOptions SPHoptions = initializeSPHOptions(parameters, csm, dTime, dDelta);
        SPHoptions.doGravity = parameters.get_bDoGravity();
        SPHoptions.nPredictRung = uRung;
        *puRungMax =
            Gravity(uRung, MaxRung(), ROOT, uRoot2, dTime, dDelta, *pdStep, dTheta, 1, bKickOpen, bEwald, bGravStep, nPartRhoLoc, iTimeStepCrit, SPHoptions);
    }

    if (DoStellarFeedback() /* or STELLAR_EVOLUTION */) {
        ActiveRung(uRung, 0);
        if (DoStellarFeedback()) {
            StellarFeedback(dTime);
        }
#ifdef STELLAR_EVOLUTION
        if (parameters.get_bChemEnrich()) {
            ChemEnrich(dTime);
        }
#endif
    }

    ActiveRung(uRung, 1);
    if (DoGas() && MeshlessHydro()) {
        EndTimestepIntegration(dTime, dDelta);
        MeshlessGradients(dTime, dDelta);
    }

    if (DoGas() && MeshlessHydro()) {
        HydroStep(dTime, dDelta);
        UpdateRung(uRung);
        uint8_t iTempRung;
        for (iTempRung = 0; iTempRung <= iCurrMaxRung; ++iTempRung) {
            if (nRung[iTempRung] == 0) continue;
            *puRungMax = iTempRung;
        }
    }

    if (!uRung && parameters.get_bFindGroups()) {
        GroupStats();
        WriteGroupStats(BuildNameWithI(iStep, ".fofstats"));
    }

    if (uRung && uRung < *puRungMax)
        bDualTree = NewTopStepKDK(dTime, dDelta, dTheta, nSteps, bDualTree, uRung + 1, pdStep, puRungMax, pbDoCheckpoint, pbDoOutput, pbNeedKickOpen);
    if (bDualTree && uRung == iRungDT + 1) {
        print_detail("Half Drift, uRung: {rung}\n", "rung"_a = iRungDT);
        dDeltaRung = dDelta / (uintmax_t(1) << iRungDT);
        Drift(dTimeFixed, 0.5 * dDeltaRung, FIXROOT);
    }

    return bDualTree;
}

void MSR::TopStepKDK(double dStep,             /* Current step */
                     double dTime,             /* Current time */
                     double dDeltaRung,        /* Time step */
                     double dTheta, int iRung, /* Rung level */
                     int iKickRung,            /* Gravity on all rungs from iRung to iKickRung */
                     int iAdjust) {            /* Do an adjust? */
    double dDeltaStep = dDeltaRung * (uintmax_t(1) << iRung);
    const auto bEwald = parameters.get_bEwald();
    const auto bGravStep = parameters.get_bGravStep();
    const auto nPartRhoLoc = parameters.get_nPartRhoLoc();

    if (iAdjust && (iRung < MaxRung() - 1)) {
        print_detail("{:>{}}Adjust, iRung: {}\n", ' ', 2 * iRung + 2, iRung);
        /* JW: Note -- can't trash uRungNew here! Force calcs set values for it! */
        ActiveRung(iRung, 1);
        if (parameters.get_bAccelStep()) AccelStep(iRung, MAX_RUNG, dTime, dDeltaStep);
        if (DoGas() && MeshlessHydro()) {
            HydroStep(dTime, dDeltaStep);
        }
        if (parameters.get_bDensityStep()) {
            DomainDecomp(iRung);
            ActiveRung(iRung, 1);
            BuildTree(0);
            DensityStep(iRung, MAX_RUNG, dTime, dDeltaStep);
        }
        if (parameters.get_bBlackholes()) {
            BHStep(dTime, dDeltaStep);
        }
        UpdateRung(iRung);
    }

    print_detail("{:>{}}msrKickOpen  at iRung: {} 0.5*dDelta: {}\n", ' ', 2 * iRung + 2, iRung, 0.5 * dDeltaRung);
    KickKDKOpen(dTime, 0.5 * dDeltaRung, iRung, iRung);
    /*
     ** Save fine-grained statistics, assuming that there is no 'sandwiched' rung,
     ** i.e., the differences in dt are smooth
     */
    if ((nRung[0] != 0 && iRung == 0) || ((nRung[iRung] == 0) && (nRung[iRung + 1] > 0))) OutputFineStatistics(dStep, dTime);

    if (CurrMaxRung() > iRung) {
        /*
        ** Recurse.
        */
        TopStepKDK(dStep, dTime, 0.5 * dDeltaRung, dTheta, iRung + 1, iRung + 1, 0);
        dTime += 0.5 * dDeltaRung;
        dStep += 1.0 / (2 << iRung);

        ActiveRung(iRung, 0); /* is this call even needed? */

        TopStepKDK(dStep, dTime, 0.5 * dDeltaRung, dTheta, iRung + 1, iKickRung, 1);
    }
    else if (CurrMaxRung() == iRung) {
        if (DoGas() && MeshlessHydro()) {
            ActiveRung(iKickRung, 1);
            MeshlessFluxes(dTime, dDeltaStep);
        }

        ZeroNewRung(iKickRung, MAX_RUNG, iKickRung); /* brute force */
        /* This Drifts everybody */
        print("{:>{}}Drift, iRung: {}\n", ' ', 2 * iRung + 2, iRung);
        Drift(dTime, dDeltaRung, ROOT);
        dTime += dDeltaRung;
        dStep += 1.0 / (uintmax_t(1) << iRung);

        if (parameters.get_bDoCooling()) {
            if (csm->val.bComove) {
                const float a = csmTime2Exp(csm, dTime);
                const float z = 1. / a - 1.;

                CoolingUpdate(z);
            }
            else {
                CoolingUpdate(0.);
            }
        }

        ActiveRung(iKickRung, 1);
        DomainDecomp(iKickRung);
        if (parameters.get_bAddDelete()) RemoveDeleted();
        BuildTree(bEwald);

        if (parameters.get_bBlackholes()) {
#ifndef DEBUG_BH_ONLY
            auto bBHAccretion = parameters.get_bBHAccretion();
            if (bBHAccretion || parameters.get_bBHFeedback()) {
                BHEvolve(dTime, dDeltaStep);
            }
            if (bBHAccretion) {
                BHAccretion(dTime);
            }
#endif
            if (parameters.get_bBHMerger()) {
                SelActives();
                BHMerger(dTime);
            }
        }

        if (parameters.get_bDoStarFormation()) StarForm(dTime, dDeltaStep, iKickRung);

#ifdef OPTIM_REORDER_IN_NODES
        if (parameters.get_bStars() || parameters.get_bBlackholes()) ReorderInNodes();
#endif

        if (!iKickRung && parameters.get_bFindGroups()) {
            NewFof(parameters.get_dTau(), parameters.get_nMinMembers());
        }

        if (DoGravity() || DoGas()) {
            print("{:>{}}Forces, iRung: {} to {}\n", ' ', 2 * iRung + 2, iKickRung, iRung);
        }

        if (DoGravity()) {
            UpdateSoft(dTime);
            SPHOptions SPHoptions = initializeSPHOptions(parameters, csm, dTime, dDeltaStep);
            SPHoptions.doGravity = 1;
            Gravity(iKickRung, MAX_RUNG, ROOT, 0, dTime, dDeltaStep, dStep, dTheta, 0, 0, bEwald, bGravStep, nPartRhoLoc, parameters.get_iTimeStepCrit(),
                    SPHoptions);
        }

        if (DoStellarFeedback() /* or STELLAR_EVOLUTION */) {
            ActiveRung(iKickRung, 0);
            double dsec;
            print("Computing feedback... ");
            TimerStart(TIMER_FEEDBACK);
            if (DoStellarFeedback()) {
                StellarFeedback(dTime);
            }
            TimerStop(TIMER_FEEDBACK);
            dsec = TimerGet(TIMER_FEEDBACK);
            print("took {:.5f} seconds\n", dsec);

#ifdef STELLAR_EVOLUTION
            print("Computing stellar evolution... ");
            TimerStart(TIMER_STEV);
            if (parameters.get_bChemEnrich()) {
                ChemEnrich(dTime);
            }
            TimerStop(TIMER_STEV);
            dsec = TimerGet(TIMER_STEV);
            print("took {:.5f} seconds\n", dsec);
#endif
        }
        ActiveRung(iKickRung, 1);

        if (DoGas() && MeshlessHydro()) {
            EndTimestepIntegration(dTime, dDeltaStep);
            MeshlessGradients(dTime, dDeltaStep);
        }

        /*
         * move time back to 1/2 step so that KickClose can integrate
         * from 1/2 through the timestep to the end.
         */
        dTime -= 0.5 * dDeltaRung;
    }
    else {
        abort();
    }

    print("{:>{}}KickClose, iRung: {}, 0.5*dDelta: {}\n", ' ', 2 * iRung + 2, iRung, 0.5 * dDeltaRung);
    KickKDKClose(dTime, 0.5 * dDeltaRung, iRung, iRung); /* uses dTime-0.5*dDelta */

    dTime += 0.5 * dDeltaRung; /* Important to have correct time at step end for SF! */

    if (!iKickRung && !iRung) {
        if (parameters.get_bFindGroups()) GroupStats();
        if (parameters.get_bAddDelete()) RemoveDeleted();
        BuildTree(bEwald);
        if (parameters.get_bBlackholes() && parameters.get_bBHPlaceSeed()) {
            PlaceBHSeed(dTime, CurrMaxRung());
        }
    }
}

void MSR::RemoveDeleted() {
    ServiceRemoveDeleted::output out;
    mdl->RunService<ServiceRemoveDeleted>({}, out);
    N = blitz::sum(out.counts);
    nGas = out.counts[FIO_SPECIES_SPH];
    nDark = out.counts[FIO_SPECIES_DARK];
    nStar = out.counts[FIO_SPECIES_STAR];
    nBH = out.counts[FIO_SPECIES_BH];
    nMaxOrder = out.nMaxOrder;
}

void MSR::CountSpecies() {
    ServiceCountSpecies::output out;
    mdl->RunService<ServiceCountSpecies>({}, out);
    assert(out.counts[FIO_SPECIES_SPH] == nGas);
    assert(out.counts[FIO_SPECIES_DARK] == nDark);
    assert(out.counts[FIO_SPECIES_STAR] == nStar);
    assert(out.counts[FIO_SPECIES_BH] == nBH);
    nMaxOrder = out.nMaxOrder;
}

/* Gas routines */
void MSR::InitSph(double dTime, double dDelta, bool bRestart) {
    if (MeshlessHydro()) {
        if (!bRestart) {
            ActiveRung(0, 1);
            InitBall();
            EndTimestepIntegration(dTime, 0.0);
            MeshlessGradients(dTime, 0.0);
            MeshlessFluxes(dTime, 0.0);
            // We do this twice because we need to have uNewRung for the time
            // limiter of Durier & Dalla Vecchia
            HydroStep(dTime, dDelta);
            HydroStep(dTime, dDelta);
        }
    }
    UpdateRung(0);
}

void MSR::CoolSetup(double dTime) {
}

void MSR::Cooling(double dTime, double dStep, int bUpdateState, int bUpdateTable, int bIterateDt) {
}

void MSR::ChemCompInit() {
    ServiceChemCompInit::input in{
        .H = parameters.get_dInitialH(),
        .He = parameters.get_dInitialHe(),
        .C = parameters.get_dInitialC(),
        .N = parameters.get_dInitialN(),
        .O = parameters.get_dInitialO(),
        .Ne = parameters.get_dInitialNe(),
        .Mg = parameters.get_dInitialMg(),
        .Si = parameters.get_dInitialSi(),
        .Fe = parameters.get_dInitialFe(),
        .Metallicity = parameters.get_dInitialMetallicity(),
    };
    mdl->RunService<ServiceChemCompInit>(in);
}

/* END Gas routines */

void MSR::WriteGroupStats(const std::string_view& fname) {
    WriteGroupStats(fname, parallel_write_count());
}

void MSR::WriteGroupStats(const std::string_view& fname, int nParallel) {
    double dsec;

    if (parameters.get_bVStep()) print("Writing group statistics to {}\n", fname);
    TimerStart(TIMER_IO);

    ServiceSaveGroupStats::input hdr;
    marshal_string(hdr.io.filename, fname);
    hdr.io.nSimultaneous = nParallel;
    hdr.io.nSegment = hdr.io.iThread = 0; // setup later
    hdr.io.iReaderWriter = 0;
    if (hdr.io.nSimultaneous > 1) strcat(hdr.io.filename, ".{i}");
    mdl->RunService<ServiceSaveGroupStats>(hdr);

    TimerStop(TIMER_IO);
    dsec = TimerGet(TIMER_IO);
    if (parameters.get_bVStep()) print("Written statistics, Wallclock: {:.5f} secs\n", dsec);
}

void MSR::Hop(double dTime, double dDelta) {
    int i;
    double sec, dsec, ssec;

    ssec = MSR::Time();

    ServiceHopLink::input il{.nSmooth = 80, .bPeriodic = parameters.get_bPeriodic()};
    ServiceHopJoin::input ij{.dHopTau = parameters.get_dHopTau(), .nSmooth = 80, .bPeriodic = parameters.get_bPeriodic()};

    if (parameters.get_bVStep()) {
        if (ij.dHopTau < 0.0) print("Running Grasshopper with adaptive linking length ({:g} times softening)\n", -ij.dHopTau);
        else print("Running Grasshopper with fixed linking length {:g}\n", ij.dHopTau);
    }

    sec = MSR::Time();
    ServiceSmooth::input in{.nSmooth = 80, .bPeriodic = parameters.get_bPeriodic(), .bSymmetric = false, .iSmoothType = SMX_DENSITY_M3};
    uint64_t mSmoothed;
    mdl->RunService<ServiceSmooth>(in, mSmoothed);
    dsec = MSR::Time() - sec;
    if (parameters.get_bVStep()) print("Density calculation complete in {:.5f} secs, finding chains...\n", dsec);

    sec = MSR::Time();
    ServiceHopLink::output nGroups;
    mdl->RunService<ServiceHopLink>(il, nGroups);
    dsec = MSR::Time() - sec;
    if (parameters.get_bVStep()) print("Chain search complete in {:.5f} secs, building minimal tree...\n", dsec);

    /* Build a new tree with only marked particles */
    sec = MSR::Time();
    BuildTreeMarked();
    dsec = MSR::Time() - sec;
    if (parameters.get_bVStep()) print("Tree build complete in {:.5f} secs, merging {} chains...\n", dsec, nGroups);

    sec = MSR::Time();
    i = 0;
    ServiceHopJoin::output j;
    do {
        ++i;
        assert(i < 100);
        mdl->RunService<ServiceHopJoin>(ij, j);
        if (parameters.get_bVStep()) print("... {} iteration{}, {} chains remain\n", i, i == 1 ? "" : "s", j.nGroups);
    } while (!j.bDone);
    nGroups = j.nGroups;
    dsec = MSR::Time() - sec;
    if (parameters.get_bVStep()) print("Chain merge complete in {:.5f} secs, {} groups\n", dsec, nGroups);
    ServiceHopFinishUp::input inFinish{
        .fPeriod = parameters.get_dPeriod(), .nMinGroupSize = static_cast<int>(parameters.get_nMinMembers()), .bPeriodic = parameters.get_bPeriodic()};
    mdl->RunService<ServiceHopFinishUp>(inFinish, nGroups);
    if (parameters.get_bVStep()) print("Removed groups with fewer than {} particles, {} remain\n", inFinish.nMinGroupSize, nGroups);
#if 0
    if (parameters.get_bVStep()) print("Unbinding\n");

    ServiceHopUnbind::input inUnbind{.dTime = dTime,
                                     .bPeriodic = parameters.get_bPeriodic(),
                                     .fPeriod = parameters.get_dPeriod(),
                                     .nMinGroupSize = static_cast<int>(parameters.get_nMinMembers())};
    ServiceHopUnbind::output outUnbind;
    ServiceHopGravity::input inGravity{.dTime = dTime, .dTheta = dThetaMin};

    int iIteration = 0;
    do {
        sec = MSR::Time();
        ServiceHopTreeBuild::input inTreeBuild{.nBucket = parameters.get_nBucket(), .nGroup = parameters.get_nGroup()};
        mdl->RunService<ServiceHopTreeBuild>(inTreeBuild);
        dsec = MSR::Time() - sec;
        if (parameters.get_bVStep()) print("... group trees built, Wallclock: {:.5f} secs\n", dsec);

        sec = MSR::Time();
        mdl->RunService<ServiceHopGravity>(inGravity);
        dsec = MSR::Time() - sec;
        if (parameters.get_bVStep()) print("... gravity complete, Wallclock: {:.5f} secs\n", dsec);

        sec = MSR::Time();
        mdl->RunService<ServiceHopUnbind>(inUnbind, outUnbind);
        nGroups = outUnbind.nGroups;
        dsec = MSR::Time() - sec;
        if (parameters.get_bVStep())
            print("... unbinding complete in {:.5f} secs, {} particles evaporated, {} groups remain\n", dsec, outUnbind.nEvaporated, nGroups);
    } while (++iIteration < 100 && outUnbind.nEvaporated);
#endif
    /*
    ** This should be done as a separate msr function.
    */
    ServiceGroupStats::input inGroupStats;
    inGroupStats.bPeriodic = parameters.get_bPeriodic();
    inGroupStats.dPeriod = parameters.get_dPeriod();
    inGroupStats.rEnvironment[0] = parameters.get_dEnvironment0();
    inGroupStats.rEnvironment[1] = parameters.get_dEnvironment1();
    inGroupStats.iGlobalStart = 1; /* global id 0 means ungrouped particle on all cpus */
    auto dBoxSize = parameters.get_dBoxSize();
    if (parameters.has_dBoxSize() && dBoxSize > 0.0) {
        inGroupStats.rEnvironment[0] /= dBoxSize;
        inGroupStats.rEnvironment[1] /= dBoxSize;
    }
    mdl->RunService<ServiceSaveGroupStats>(inGroupStats);

    dsec = MSR::Time() - ssec;
    if (parameters.get_bVStep()) print("Grasshopper complete, Wallclock: {:.5f} secs\n", dsec);
}

void MSR::NewFof(double dTau, int nMinMembers) {
    double dsec;

    TimerStart(TIMER_FOF);

    ServiceNewFof::input in{.dTau2 = dTau * dTau,
                            .nMinMembers = nMinMembers,
                            .bPeriodic = parameters.get_bPeriodic(),
                            .nReplicas = static_cast<int>(parameters.get_bPeriodic() ? parameters.get_nReplicas() : 0),
                            .nBucket = static_cast<int>(parameters.get_nBucket())};

    if (parameters.get_bVStep()) {
        print("Running FoF with linking length {}\n", dTau);
    }

    TimerStart(TIMER_NONE);
    mdl->RunService<ServiceNewFof>(in);

    TimerStop(TIMER_NONE);
    dsec = TimerGet(TIMER_NONE);
    if (parameters.get_bVStep()) print("Initial FoF calculation complete in {:.5f} secs\n", dsec);

    TimerStart(TIMER_NONE);
    int i = 0;
    ServiceFofPhases::output out{.bMadeProgress = true};
    do {
        ++i;
        assert(i < 100);
        mdl->RunService<ServiceFofPhases>({}, out);
        if (parameters.get_bVStep()) print("... {} iteration{}\n", i, i == 1 ? "" : "s");
    } while (out.bMadeProgress);

    TimerStop(TIMER_NONE);
    dsec = TimerGet(TIMER_NONE);
    if (parameters.get_bVStep()) print("Global merge complete in {:.5f} secs\n", dsec);

    uint64_t nGroups;
    ServiceFofFinishUp::input inFinish{.nMinGroupSize = nMinMembers};
    mdl->RunService<ServiceFofFinishUp>(inFinish, nGroups);

    if (parameters.get_bVStep()) print("Removed groups with fewer than {} particles, {} remain\n", inFinish.nMinGroupSize, nGroups);
    TimerStop(TIMER_FOF);
    dsec = TimerGet(TIMER_FOF);
    if (parameters.get_bVStep()) print("FoF complete, Wallclock: {:.5f} secs\n", dsec);
}

void MSR::GroupStats() {
    if (parameters.get_bVStep()) print("Generating Group statistics\n");
    TimerStart(TIMER_FOF);
    ServiceGroupStats::input inGroupStats;
    inGroupStats.bPeriodic = parameters.get_bPeriodic();
    inGroupStats.dPeriod = parameters.get_dPeriod();
    inGroupStats.rEnvironment[0] = parameters.get_dEnvironment0();
    inGroupStats.rEnvironment[1] = parameters.get_dEnvironment1();
    inGroupStats.iGlobalStart = 1; /* global id 0 means ungrouped particle on all cpus */
    auto dBoxSize = parameters.get_dBoxSize();
    if (parameters.has_dBoxSize() && dBoxSize > 0.0) {
        inGroupStats.rEnvironment[0] /= dBoxSize;
        inGroupStats.rEnvironment[1] /= dBoxSize;
    }
    mdl->RunService<ServiceGroupStats>(inGroupStats);
    TimerStop(TIMER_FOF);
    auto dsec = TimerGet(TIMER_FOF);
    if (parameters.get_bVStep()) print("Group statistics complete, Wallclock: {:.5f} secs\n", dsec);
}

#ifdef MDL_FFTW
double MSR::GenerateIC(int nGrid, int iSeed, double z, double L, CSM csm) {
    fioSpeciesList nSpecies;
    double sec, dsec;
    double mean, rms;
    uint64_t nTotal;
    int j;

    if (csm) {
        this->csm->val = csm->val;
        if (this->csm->val.classData.bClass) csmClassGslInitialize(this->csm);
    }
    else csm = this->csm;

    // We only support periodic initial conditions
    parameters.set_bPeriodic(true);
    parameters.set_bComove(true);

    parameters.set_iSeed(iSeed);
    parameters.set_dRedFrom(z);
    parameters.set_dBoxSize(L);
    parameters.set_nGrid(nGrid);

    auto inptr = std::make_unique<ServiceGenerateIC::input>(); // This has largish tables
    auto& in = *inptr;
    ServiceGenerateIC::output out;
    in.dBoxSize = L;
    in.iSeed = iSeed;
    in.bFixed = parameters.get_bFixedAmpIC();
    in.fPhase = parameters.get_dFixedAmpPhasePI() * M_PI;
    in.nGrid = nGrid;
    in.iLPT = parameters.get_iLPT();
    in.bICgas = parameters.get_bICgas();
    in.nBucket = parameters.get_nBucket();
    in.dInitialT = parameters.get_dInitialT();
    in.dInitialH = parameters.get_dInitialH();
    in.dInitialHe = parameters.get_dInitialHe();
    in.dInitialC = parameters.get_dInitialC();
    in.dInitialN = parameters.get_dInitialN();
    in.dInitialO = parameters.get_dInitialO();
    in.dInitialNe = parameters.get_dInitialNe();
    in.dInitialMg = parameters.get_dInitialMg();
    in.dInitialSi = parameters.get_dInitialSi();
    in.dInitialFe = parameters.get_dInitialFe();
    in.dInitialMetallicity = parameters.get_dInitialMetallicity();

    nTotal = in.nGrid; /* Careful: 32 bit integer cubed => 64 bit integer */
    nTotal *= in.nGrid;
    nTotal *= in.nGrid;
    in.dBoxMass = csm->val.dOmega0 / nTotal;
    if (in.bICgas) nTotal *= 2;

    for (j = 0; j <= FIO_SPECIES_LAST; j++) nSpecies[j] = 0;
    if (in.bICgas) {
        nSpecies[FIO_SPECIES_ALL] = nTotal;
        nSpecies[FIO_SPECIES_SPH] = nTotal / 2;
        nSpecies[FIO_SPECIES_DARK] = nTotal / 2;
    }
    else {
        nSpecies[FIO_SPECIES_ALL] = nSpecies[FIO_SPECIES_DARK] = nTotal;
    }
    InitializePStore(nSpecies, getMemoryModel(), parameters.get_nMemEphemeral()); // We now need a bit of cosmology to set the maximum lightcone depth here.
    InitCosmology(csm);

    in.dBaryonFraction = csm->val.dOmegab / csm->val.dOmega0;
    SetDerivedParameters();
    in.dTuFac = dTuFacPrimNeutral;

    assert(z >= 0.0);
    in.dExpansion = 1.0 / (1.0 + z);

    N = nSpecies[FIO_SPECIES_ALL];
    nGas = nSpecies[FIO_SPECIES_SPH];
    nDark = nSpecies[FIO_SPECIES_DARK];
    nStar = nSpecies[FIO_SPECIES_STAR];
    nBH = nSpecies[FIO_SPECIES_BH];
    nMaxOrder = N - 1; // iOrder goes from 0 to N-1

    if (parameters.get_bVStart()) print("Generating IC...\nN:{} nDark:{} nGas:{} nStar:{}\n", N, nDark, nGas, nStar);

    /* Read the transfer function */
    in.nTf = 0;
    if (parameters.has_achTfFile()) {
        auto achTfFile = parameters.get_achTfFile();
        FILE* fp = fopen(achTfFile.data(), "r");
        char buffer[PATH_MAX];

        if (parameters.get_bVStart()) print("Reading transfer function from {}\n", achTfFile.data());
        if (fp == NULL) {
            perror(achTfFile.data());
            Exit(1);
        }
        while (fgets(buffer, sizeof(buffer), fp)) {
            assert(in.nTf < ServiceGenerateIC::max_tf);
            if (sscanf(buffer, " %lg %lg\n", &in.k[in.nTf], &in.tf[in.nTf]) == 2) {
                in.k[in.nTf] = log(in.k[in.nTf]);
                in.tf[in.nTf] = log(in.tf[in.nTf]);
                ++in.nTf;
            }
        }
        fclose(fp);
        if (parameters.get_bVStart()) print("Transfer function : {} lines kmin {} kmax {}\n", in.nTf, exp(in.k[0]), exp(in.k[in.nTf - 1]));
    }

    sec = MSR::Time();

    /* Figure out the minimum number of particles */
    ServiceGridCreateFFT::input gridIn{.nGrid = nGrid};
    ServiceGridCreateFFT::output gridOut;
    mdl->RunService<ServiceGridCreateFFT>(gridIn, gridOut);

    print_detail("IC Generation @ a={a} with seed {seed}\n", "a"_a = in.dExpansion, "seed"_a = iSeed);
    in.nPerNode = gridOut.nMaxLocal;
    mdl->RunService<ServiceGenerateIC>(in, out);
    mean = 2 * out.noiseMean / N;
    rms = sqrt(2 * out.noiseCSQ / N);

    print_detail("Transferring particles between/within nodes\n");
    pstMoveIC(pst, &in, sizeof(in), NULL, 0);

    mdl->RunService<ServiceGridDeleteFFT>();

    SetClasses();
    dsec = MSR::Time() - sec;
    print_detail("IC Generation Complete @ a={a:.6}, Wallclock: {seconds:.5} secs\n", "a"_a = in.dExpansion, "seconds"_a = dsec);
    print_detail("Mean of noise sample is {mean:g}, RMS {rms:g}\n", "mean"_a = mean, "rms"_a = rms);
    ReorderDecorator();

    return getTime(in.dExpansion);
}
#endif

double MSR::Read(std::string_view achInFile) {
    double dTime, dExpansion;
    FIO fio;
    double dsec;
    fioSpeciesList nSpecies;
    inReadFileFilename achFilename;
    uint64_t mMemoryModel = 0;

    mMemoryModel = getMemoryModel();

    TimerStart(TIMER_NONE);

    auto nBytes = PST_MAX_FILES * (sizeof(fioSpeciesList) + PATH_MAX);
    std::unique_ptr<char[]> buffer{new char[sizeof(inReadFile) + nBytes]};
    auto read = new (buffer.get()) inReadFile;

    /* Add Data Subpath for local and non-local names. */
    marshal_string(achFilename, achInFile);
    print("Reading from {}\n", achFilename);
    fio = fioOpen(achFilename, csm->val.dOmega0, csm->val.dOmegab);
    if (fio == NULL) {
        print_error("ERROR: unable to open input file\n");
        perror(achFilename);
        Exit(1);
    }
    nBytes = fioDump(fio, nBytes, read + 1);

    // If we have the 'Redshift' field, we take that.
    // If not, we assume that the 'Time' field contains the expansion factor
    if (!fioGetAttr(fio, HDF5_HEADER_G, "Redshift", FIO_TYPE_DOUBLE, &dExpansion)) {
        if (!fioGetAttr(fio, HDF5_HEADER_G, "Time", FIO_TYPE_DOUBLE, &dExpansion)) dExpansion = 0.0;
    }
    else {
        dExpansion = 1.0 / (dExpansion + 1.0);
    }
    if (!fioGetAttr(fio, HDF5_HEADER_G, "dEcosmo", FIO_TYPE_DOUBLE, &dEcosmo)) dEcosmo = 0.0;
    if (!fioGetAttr(fio, HDF5_HEADER_G, "dTimeOld", FIO_TYPE_DOUBLE, &dTimeOld)) dTimeOld = 0.0;
    if (!fioGetAttr(fio, HDF5_HEADER_G, "dUOld", FIO_TYPE_DOUBLE, &dUOld)) dUOld = 0.0;

    if (csm->val.bComove) {
        if (!parameters.has_dOmega0()) fioGetAttr(fio, HDF5_COSMO_G, "Omega_m", FIO_TYPE_DOUBLE, &csm->val.dOmega0);
        if (!parameters.has_dLambda()) fioGetAttr(fio, HDF5_COSMO_G, "Omega_Lambda", FIO_TYPE_DOUBLE, &csm->val.dLambda);
        if (!parameters.has_dBoxSize()) {
            double dBoxSize;
            fioGetAttr(fio, HDF5_HEADER_G, "BoxSize", FIO_TYPE_DOUBLE, &dBoxSize);
            parameters.set_dBoxSize(dBoxSize);
        }
        if (!parameters.has_h()) fioGetAttr(fio, HDF5_COSMO_G, "HubbleParam", FIO_TYPE_DOUBLE, &csm->val.h);
    }

    N = fioGetN(fio, FIO_SPECIES_ALL);
    nGas = fioGetN(fio, FIO_SPECIES_SPH);
    nDark = fioGetN(fio, FIO_SPECIES_DARK);
    nStar = fioGetN(fio, FIO_SPECIES_STAR);
    nBH = fioGetN(fio, FIO_SPECIES_BH);
    nMaxOrder = N - 1; // iOrder goes from 0 to N-1

    read->nProcessors = parallel_read_count();

    auto& nFiles = fio->fileList.nFiles;
    print(" Reading {} particles from {} file{} using {} processor{}\n", N, nFiles, (nFiles == 1 ? "" : "s"), read->nProcessors,
          (read->nProcessors == 1 ? "" : "s"));

    dTime = getTime(dExpansion);
    if (parameters.get_bInFileLC()) read->dvFac = 1.0;
    else read->dvFac = getVfactor(dExpansion);

    if (nGas && !DoGas()) parameters.set_hydro_model(HYDRO_MODEL::SPH);
    if (NewSPH()) mMemoryModel |= (PKD_MODEL_NEW_SPH | PKD_MODEL_ACCELERATION | PKD_MODEL_VELOCITY | PKD_MODEL_DENSITY | PKD_MODEL_BALL | PKD_MODEL_NODE_BOB);
    if (NewSPH() && parameters.get_bShearStrengthModel()) mMemoryModel |= PKD_MODEL_STRENGTH;
    if (parameters.get_bTrackPeakPressure()) mMemoryModel |= PKD_MODEL_PEAKPRESSURE;
    if (nStar) mMemoryModel |= PKD_MODEL_STAR;

    for (auto s = FIO_SPECIES(0); s <= FIO_SPECIES_LAST; s = FIO_SPECIES(s + 1)) {
        nSpecies[s] = fioGetN(fio, s);
        read->nSpeciesStart[s] = 0;
        read->nSpeciesEnd[s] = nSpecies[s];
    }
    InitializePStore(nSpecies, mMemoryModel, parameters.get_nMemEphemeral());

    read->dOmega0 = csm->val.dOmega0;
    read->dOmegab = csm->val.dOmegab;

    SetDerivedParameters();
    read->dTuFac = dTuFac;

    if (read->nProcessors > 1) {
        fioClose(fio);
        // NB: For the parallel read, the nSpeciesEnd is the index of the last
        // element to be read, whereas for the serial read (OneNodeRead) it holds
        // the number of elements to read (always per-species)
        for (auto s = FIO_SPECIES(0); s <= FIO_SPECIES_LAST; s = FIO_SPECIES(s + 1)) {
            if (nSpecies[s] > 0) {
                read->nSpeciesEnd[s] -= 1;
            }
        }
        pstReadFile(pst, read, sizeof(struct inReadFile) + nBytes, NULL, 0);
    }
    else {
        OneNodeRead(read, fio);
        fioClose(fio);
    }

    TimerStop(TIMER_NONE);
    dsec = TimerGet(TIMER_NONE);
    SetClasses();
    print("Input file has been successfully read, Wallclock: {:.5f} secs.\n", dsec);
    ReorderDecorator();

    /*
    ** If this is a non-periodic box, then we must precalculate the bounds.
    ** We throw away the result, but PKD will keep track for later.
    */
    auto period = parameters.get_dPeriod();
    if (!parameters.get_bPeriodic() || blitz::any(period >= FLOAT_MAXVAL)) {
        CalcBound();
    }

    InitCosmology(csm);

    if (NewSPH()) {
        const auto bEwald = parameters.get_bEwald();
        /*
        ** Initialize kernel target with nSmooth
        */
        SetSPHoptions();
        InitializeEOS();

        if (parameters.has_dSoft()) SetSoft(Soft());
        /*
        ** Initialize fBall
        */
        TimerStart(TIMER_NONE);
        int fBallInitnSmooth = int(double(N) / double(nGas) * 2 * parameters.get_nSmooth());
        print("Initializing fBall (Smooth with {} neighbors) ...\n", fBallInitnSmooth);
        Reorder();
        ActiveRung(0, 1); /* Activate all particles */
        DomainDecomp(-1);
        BuildTree(bEwald);
        Smooth(dTime, 0.0f, SMX_BALL, 0, fBallInitnSmooth);
        Reorder();
        TimerStop(TIMER_NONE);
        dsec = TimerGet(TIMER_NONE);
        print("Initializing fBall complete, Wallclock: {:.5f} secs.\n", dsec);

        /*
        ** Convert U
        */
        TimerStart(TIMER_NONE);
        print("Converting u ...\n");
        ActiveRung(0, 1); /* Activate all particles */
        DomainDecomp(-1);
        BuildTree(bEwald);

        const auto iStartStep = parameters.get_iStartStep();
        auto dTheta = set_dynamic(iStartStep, dTime);

        // Calculate Density
        SPHOptions SPHoptions = initializeSPHOptions(parameters, csm, dTime, /*dDelta*/ 0.0);
        SPHoptions.doDensity = 1;
        SPHoptions.doUConversion = 1;
        const auto bGravStep = parameters.get_bGravStep();
        const auto nPartRhoLoc = parameters.get_nPartRhoLoc();
        const auto iTimeStepCrit = parameters.get_iTimeStepCrit();
        Gravity(0, MAX_RUNG, ROOT, 0, dTime, 0.0f, iStartStep, dTheta, 0, 1, bEwald, bGravStep, nPartRhoLoc, iTimeStepCrit, SPHoptions);
        MemStatus();
        if (SPHoptions.doInterfaceCorrection) {
            SPHoptions.doDensity = 0;
            SPHoptions.doDensityCorrection = 1;
            SPHoptions.dofBallFactor = 0;
            TreeUpdateFlagBounds(bEwald, ROOT, 0, SPHoptions);
            Gravity(0, MAX_RUNG, ROOT, 0, dTime, 0.0f, iStartStep, dTheta, 0, 1, bEwald, bGravStep, nPartRhoLoc, iTimeStepCrit, SPHoptions);
            UpdateGasValues(0, dTime, 0.0f, iStartStep, 0, 1, SPHoptions);
        }
        TimerStop(TIMER_NONE);
        dsec = TimerGet(TIMER_NONE);
        print("Converting u complete, Wallclock: {:.5f} secs.\n", dsec);
        if (parameters.get_bWriteIC() || (NoSteps())) {
            Write(BuildIoName(iStartStep).c_str(), dTime, 0);
        }
        if (NoSteps()) std::_Exit(EXIT_FAILURE); // For now. this will be replaced sood.
    }

    return dTime;
}

// This sets the local pkd->bnd.
void MSR::CalcBound(Bound& bnd) {
    mdl->RunService<ServiceCalcBound>({}, bnd);
}
void MSR::CalcBound() {
    Bound bnd;
    CalcBound(bnd);
}

void MSR::OutputGrid(std::string_view filename, bool k, int iGrid, int nParaWrite) {
    struct inOutput out;
    double dsec, sec = MSR::Time();
    out.eOutputType = k ? OUT_KGRID : OUT_RGRID;
    out.iGrid = iGrid;
    out.iPartner = -1;
    out.nPartner = -1;
    out.iProcessor = 0;
    out.nProcessor = nParaWrite > mdl->Procs() ? mdl->Procs() : nParaWrite;
    marshal_string(out.achOutFile, filename);
    print("Writing grid to {}\n", filename);
    pstOutput(pst, &out, sizeof(out), NULL, 0);
    dsec = MSR::Time() - sec;
    print_detail("Grid has been successfully written, Wallclock: {seconds:.5} secs.\n\n", "seconds"_a = dsec);
}

#ifdef MDL_FFTW
void MSR::OutputPk(int iStep, double dTime) {
    double a, z, vfact, kfact;
    std::string filename;
    int i;

    auto nGridPk = parameters.get_nGridPk();
    if (nGridPk == 0) return;

    auto nBinsPk = parameters.get_nBinsPk();

    if (!csm->val.bComove) a = 1.0;
    else a = csmTime2Exp(csm, dTime);

    /* If the Box Size (in mpc/h) was specified, then we can scale the output power spectrum measurement */
    if (parameters.has_dBoxSize() && parameters.get_dBoxSize() > 0.0) kfact = parameters.get_dBoxSize();
    else kfact = 1.0;
    vfact = kfact * kfact * kfact;
    kfact = 1.0 / kfact;

    std::vector<std::int64_t> nFoldPk;
    if (parameters.has_nFoldPk()) nFoldPk = parameters.get_nFoldPk();
    else nFoldPk.push_back(1);

    for (auto fold : nFoldPk) {
        auto [nPk, fK, fPk, fPkAll] = MeasurePk(int(parameters.get_iPkOrder()), parameters.get_bPkInterlace(), nGridPk, a, nBinsPk, fold);

        filename = BuildName(iStep, ".pk");
        if (nFoldPk.size() > 1) filename = fmt::format("{name}.{fold}", "name"_a = filename, "fold"_a = fold);
        std::ofstream fs(filename);
        if (fs.fail()) {
            std::cerr << "Could not create P(k) file: " << filename << std::endl;
            perror(filename.c_str());
            Exit(errno);
        }
        fmt::print(fs, "# k P(k) N(k) P(k)+{linear}\n", "linear"_a = parameters.get_achPkSpecies());
        fmt::print(fs, "# a={a:.8f}  z={z:.8f}\n", "a"_a = a, "z"_a = 1 / a - 1.0);
        for (i = 0; i < nBinsPk; ++i) {
            if (fPk[i] > 0.0)
                fmt::print(fs, "{k:.8e} {pk:.8e} {nk} {all:.8e}\n", "k"_a = kfact * fK[i] * 2.0 * M_PI * fold, "pk"_a = vfact * fPk[i], "nk"_a = nPk[i],
                           "all"_a = vfact * fPkAll[i]);
        }
        fs.close();
    }

    /* Output the k-grid if requested */
    z = 1 / a - 1;
    auto iDeltakInterval = parameters.get_iDeltakInterval();
    auto dDeltakRedshift = parameters.get_dDeltakRedshift();
    if (iDeltakInterval && (iStep % iDeltakInterval == 0) && z < dDeltakRedshift) {
        auto filename = BuildNameWithI(iStep, ".deltak");
        OutputGrid(filename.c_str(), true, 0, parallel_write_count());
    }
}

void MSR::OutputLinearPk(int iStep, double dTime) {
    std::string filename;
    double a, vfact, kfact;
    int i;

    if (parameters.get_nGridLin() == 0) return;
    if (!csm->val.bComove) return;
    if (!parameters.has_dBoxSize()) return;

    a = csmTime2Exp(csm, dTime);

    auto [nPk, fK, fPk] = MeasureLinearPk(parameters.get_nGridLin(), a, parameters.get_dBoxSize());

    if (!csm->val.bComove) a = 1.0;
    else a = csmTime2Exp(csm, dTime);

    if (parameters.get_dBoxSize() > 0.0) kfact = parameters.get_dBoxSize();
    else kfact = 1.0;
    vfact = kfact * kfact * kfact;
    kfact = 1.0 / kfact;

    filename = BuildName(iStep, ".lin_pk");
    std::ofstream fs(filename);
    if (fs.fail()) {
        std::cerr << "Could not create P_lin(k) file: " << filename << std::endl;
        perror(filename.c_str());
        Exit(errno);
    }
    for (i = 0; i < fK.size(); ++i) {
        if (fPk[i] > 0.0) fmt::print(fs, "{k} {pk} {nk}\n", "k"_a = kfact * fK[i] * 2.0 * M_PI, "pk"_a = vfact * fPk[i], "nk"_a = nPk[i]);
    }
    fs.close();
}

#endif

/*
 **  This routine will output all requested files and fields
 */

void MSR::Output(int iStep, double dTime, double dDelta, int bCheckpoint) {
    int bSymmetric;

    // IA: If we allow for adding/deleting particles, we need to recount them to have the
    //  correct number of particles per specie
    if (parameters.get_bAddDelete()) CountSpecies();

    print("Writing output for step {}\n", iStep);

    if (DoGas() && NoSteps()) { /* Diagnostic Gas */
        Reorder();
        OutArray(BuildName(iStep, ".c"), OUT_C_ARRAY);
        OutArray(BuildName(iStep, ".hsph"), OUT_HSPH_ARRAY);
    }

    if (DoDensity()) {
        ActiveRung(0, 1); /* Activate all particles */
        DomainDecomp(-1);
        BuildTree(0);
        bSymmetric = 0; /* should be set in param file! */
        Smooth(dTime, dDelta, SMX_DENSITY, bSymmetric, parameters.get_nSmooth());
    }
    if (parameters.get_bFindGroups()) {
        Reorder();
        // OutArray(BuildNameWithI(iStep,".fof"),OUT_GROUP_ARRAY);
        WriteGroupStats(BuildNameWithI(iStep, ".fofstats"));
    }
    if (parameters.get_bFindHopGroups()) {
        ActiveRung(0, 1); /* Activate all particles */
        DomainDecomp(-1);
        BuildTree(0);
        Hop(dTime, dDelta);
        Reorder();
        // OutArray(BuildName(iStep,".hopgrp"),OUT_GROUP_ARRAY);
        WriteGroupStats(BuildNameWithI(iStep, ".hopstats"));
    }

    if (parameters.get_bDoAccOutput()) {
        Reorder();
        OutVector(BuildName(iStep, ".acc"), OUT_ACCEL_VECTOR);
    }
    if (parameters.get_bDoPotOutput()) {
        Reorder();
        OutArray(BuildName(iStep, ".pot"), OUT_POT_ARRAY);
    }

    if (DoDensity()) {
        Reorder();
        OutArray(BuildName(iStep, ".den"), OUT_DENSITY_ARRAY);
    }
    if (parameters.get_bDoRungOutput()) {
        Reorder();
        OutArray(BuildName(iStep, ".rung"), OUT_RUNG_ARRAY);
    }
    if (parameters.get_bDoRungDestOutput()) {
        Reorder();
        OutArray(BuildName(iStep, ".rd"), OUT_RUNGDEST_ARRAY);
    }
    if (parameters.get_bDoSoftOutput()) {
        Reorder();
        OutArray(BuildName(iStep, ".soft"), OUT_SOFT_ARRAY);
    }
    Write(BuildIoName(iStep).c_str(), dTime, bCheckpoint);
}

uint64_t MSR::CountSelected() {
    uint64_t N;
    mdl->RunService<ServiceCountSelected>({}, N);
    return N;
}
uint64_t MSR::SelSpecies(uint64_t mSpecies, int setIfTrue, int clearIfFalse) {
    uint64_t N;
    ServiceSelSpecies::input in(mSpecies, setIfTrue, clearIfFalse);
    mdl->RunService<ServiceSelSpecies>(in, N);
    return N;
}
uint64_t MSR::SelAll(int setIfTrue, int clearIfFalse) {
    return SelSpecies(1 << FIO_SPECIES_ALL, setIfTrue, clearIfFalse);
}
uint64_t MSR::SelGas(int setIfTrue, int clearIfFalse) {
    return SelSpecies(1 << FIO_SPECIES_SPH, setIfTrue, clearIfFalse);
}
uint64_t MSR::SelStar(int setIfTrue, int clearIfFalse) {
    return SelSpecies(1 << FIO_SPECIES_STAR, setIfTrue, clearIfFalse);
}
uint64_t MSR::SelDark(int setIfTrue, int clearIfFalse) {
    return SelSpecies(1 << FIO_SPECIES_DARK, setIfTrue, clearIfFalse);
}
uint64_t MSR::SelDeleted(int setIfTrue, int clearIfFalse) {
    // The "UNKNWON" species marks deleted particles
    return SelSpecies(1 << FIO_SPECIES_UNKNOWN, setIfTrue, clearIfFalse);
}
uint64_t MSR::SelActives(int setIfTrue, int clearIfFalse) {
    uint64_t N;
    ServiceSelActives::input in(setIfTrue, clearIfFalse);
    mdl->RunService<ServiceSelActives>(in, N);
    return N;
}
uint64_t MSR::SelBlackholes(int setIfTrue, int clearIfFalse) {
    uint64_t N;
    ServiceSelBlackholes::input in(setIfTrue, clearIfFalse);
    mdl->RunService<ServiceSelBlackholes>(in, N);
    return N;
}
uint64_t MSR::SelGroup(int iGroup, int setIfTrue, int clearIfFalse) {
    uint64_t N;
    ServiceSelGroup::input in(iGroup, setIfTrue, clearIfFalse);
    mdl->RunService<ServiceSelGroup>(in, N);
    return N;
}
uint64_t MSR::SelById(uint64_t idStart, uint64_t idEnd, int setIfTrue, int clearIfFalse) {
    uint64_t N;
    ServiceSelById::input in(idStart, idEnd, setIfTrue, clearIfFalse);
    mdl->RunService<ServiceSelById>(in, N);
    return N;
}
uint64_t MSR::SelMass(double dMinMass, double dMaxMass, int setIfTrue, int clearIfFalse) {
    uint64_t N;
    ServiceSelMass::input in(dMinMass, dMaxMass, setIfTrue, clearIfFalse);
    mdl->RunService<ServiceSelMass>(in, N);
    return N;
}
uint64_t MSR::SelPhaseDensity(double dMinPhaseDensity, double dMaxPhaseDensity, int setIfTrue, int clearIfFalse) {
    uint64_t N;
    ServiceSelPhaseDensity::input in(dMinPhaseDensity, dMaxPhaseDensity, setIfTrue, clearIfFalse);
    mdl->RunService<ServiceSelPhaseDensity>(in, N);
    return N;
}
uint64_t MSR::SelBox(blitz::TinyVector<double, 3> center, blitz::TinyVector<double, 3> size, int setIfTrue, int clearIfFalse) {
    uint64_t N;
    ServiceSelBox::input in(center, size, setIfTrue, clearIfFalse);
    mdl->RunService<ServiceSelBox>(in, N);
    return N;
}
uint64_t MSR::SelSphere(blitz::TinyVector<double, 3> r, double dRadius, int setIfTrue, int clearIfFalse) {
    uint64_t N;
    ServiceSelSphere::input in(r, dRadius, setIfTrue, clearIfFalse);
    mdl->RunService<ServiceSelSphere>(in, N);
    return N;
}
uint64_t MSR::SelCylinder(blitz::TinyVector<double, 3> dP1, blitz::TinyVector<double, 3> dP2, double dRadius, int setIfTrue, int clearIfFalse) {
    uint64_t N;
    ServiceSelCylinder::input in(dP1, dP2, dRadius, setIfTrue, clearIfFalse);
    mdl->RunService<ServiceSelCylinder>(in, N);
    return N;
}

double MSR::TotalMass() {
    ServiceTotalMass::output out;
    mdl->RunService<ServiceTotalMass>({}, out);
    return out;
}

void MSR::CalcDistance(blitz::TinyVector<double, 3> dCenter, double dRadius) {
    ServiceCalcDistance::input in{.dCenter = dCenter, .dRadius = dRadius, .bPeriodic = parameters.get_bPeriodic()};
    mdl->RunService<ServiceCalcDistance>(in);
}

void MSR::CalcCOM(blitz::TinyVector<double, 3> dCenter, double dRadius, blitz::TinyVector<double, 3>& com, blitz::TinyVector<double, 3>& vcm,
                  blitz::TinyVector<double, 3>& L, double& M) {
    ServiceCalcCOM::input in{.dCenter = dCenter, .dRadius = dRadius, .bPeriodic = parameters.get_bPeriodic()};
    ServiceCalcCOM::output out;
    mdl->RunService<ServiceCalcCOM>(in, out);
    M = out.M;
    if (out.M > 0.0) {
        com = out.com / out.M;
        vcm = out.vcm / out.M;
        auto T = blitz::cross(com, vcm);
        L = (out.L - out.M * T) / out.M;
    }
}

void MSR::SetSPHoptions() {
    struct inSetSPHoptions in;
    in.SPHoptions = initializeSPHOptions(parameters, csm, 0.0, 0.0);
    pstSetSPHoptions(pst, &in, sizeof(in), NULL, 0);
}

void MSR::ResetCOM() {
    blitz::TinyVector<double, 3> dCenter(0.0), com, vcm, L;
    double M;
    CalcCOM(dCenter, -1.0, com, vcm, L, M);
    print("Before reseting: x_com = {:.5e}, y_com = {:.5e}, z_com = {:.5e}, vx_com = {:.5e}, vy_com = {:.5e}, vz_com = {:.5e}\n", com[0], com[1], com[2],
          vcm[0], vcm[1], vcm[2]);

    ServiceResetCOM::input in{.r_com = com, .v_com = vcm};
    mdl->RunService<ServiceResetCOM>(in);

    CalcCOM(dCenter, -1.0, com, vcm, L, M);
    print("After reseting: x_com = {:.5e}, y_com = {:.5e}, z_com = {:.5e}, vx_com = {:.5e}, vy_com = {:.5e}, vz_com = {:.5e}\n", com[0], com[1], com[2], vcm[0],
          vcm[1], vcm[2]);
}

void MSR::InitializeEOS() {
    double sec, dsec;
    sec = MSR::Time();
    print("Initialize EOS ...\n");
    pstInitializeEOS(pst, NULL, 0, NULL, 0);
    dsec = MSR::Time() - sec;
    print("EOS initialized, Wallclock: {:.5f} secs\n\n", dsec);
}

void MSR::CalculateKickParameters(struct pkdKickParameters* kick, uint8_t uRungLo, double dTime, double dDelta, double dStep, int bKickClose, int bKickOpen,
                                  SPHOptions SPHoptions) {
    uint8_t uRungLoTemp;
    int i;
    double dt;
    if (SPHoptions.doDensity || SPHoptions.doDensityCorrection) {
        uRungLoTemp = uRungLo;
        uRungLo = SPHoptions.nPredictRung;
    }

    /*
    ** Now calculate the timestepping factors for kick close and open if the
    ** gravity should kick the particles. If the code uses bKickClose and
    ** bKickOpen it no longer needs to store accelerations per particle.
    */
    kick->bKickClose = bKickClose;
    kick->bKickOpen = bKickOpen;
    if (SPHoptions.doGravity || SPHoptions.doSPHForces || SPHoptions.doDensity || SPHoptions.doDensityCorrection) {
        for (i = 0, dt = 0.5 * dDelta; i <= parameters.get_iMaxRung(); ++i, dt *= 0.5) {
            kick->dtClose[i] = 0.0;
            kick->dtOpen[i] = 0.0;
            if (i >= uRungLo) {
                if (csm->val.bComove) {
                    if (bKickClose) {
                        kick->dtClose[i] = csmComoveKickFac(csm, dTime - dt, dt);
                    }
                    if (bKickOpen) {
                        kick->dtOpen[i] = csmComoveKickFac(csm, dTime, dt);
                    }
                }
                else {
                    if (bKickClose) kick->dtClose[i] = dt;
                    if (bKickOpen) kick->dtOpen[i] = dt;
                }
            }
        }
    }

    /*
    ** Create the deltas for the on-the-fly prediction of velocity and the
    ** thermodynamical variable.
    */
    if (SPHoptions.doSPHForces || SPHoptions.doDensity || SPHoptions.doDensityCorrection) {
        double substepWeAreAt = dStep - floor(dStep); // use fmod instead
        double stepStartTime = dTime - substepWeAreAt * dDelta;
        for (i = 0; i <= parameters.get_iMaxRung(); ++i) {
            if (i < uRungLo) {
                /*
                ** For particles with a step larger than the current rung, the temporal position of
                ** the velocity in relation to the current time is nontrivial, so we calculate it here
                */
                double substepSize = 1.0 / (uintmax_t(1) << i); // 1.0 / (1 << i);
                double substepsDoneAtThisSize = floor(substepWeAreAt / substepSize);
                double TPredDrift = stepStartTime + (substepsDoneAtThisSize + 0.5) * substepSize * dDelta;
                double dtPredDrift = dTime - TPredDrift;
                /* Now that we know how much we have to drift, we can calculate the corresponding
                ** drift factor
                */
                if (csm->val.bComove) {
                    /*
                    ** This gives the correct result, even if dtPredDrift is negative
                    ** but we still may want to use
                    ** -csmComoveKickFac(csm,TPredDrift + dtPredDrift,-dtPredDrift);
                    ** if dtPredDrift is negative, just to be sure
                    */
                    kick->dtPredDrift[i] = csmComoveKickFac(csm, TPredDrift, dtPredDrift);
                }
                else {
                    kick->dtPredDrift[i] = dtPredDrift;
                }
            }
            else {
                /*
                ** In this case, all particles are synchronized, which means that
                ** velocity and the thermodynamical variable are either a half step behind
                ** or ahead, so all information is contained in dtOpen and dtClose and the
                ** bMarked flag.
                */
                kick->dtPredDrift[i] = 0.0;
            }
        }
    }

    /*
    ** Create the deltas for the on-the-fly prediction in case of ISPH
    */
    if (NewSPH() && (SPHoptions.doSPHForces || SPHoptions.doDensity || SPHoptions.doDensityCorrection) && SPHoptions.useIsentropic) {
        double substepWeAreAt = dStep - floor(dStep); // use fmod instead
        double stepStartTime = dTime - substepWeAreAt * dDelta;
        for (i = 0; i <= parameters.get_iMaxRung(); ++i) {
            double substepSize = 1.0 / (uintmax_t(1) << i); // 1.0 / (1 << i);
            double substepsDoneAtThisSize = floor(substepWeAreAt / substepSize);
            double TSubStepStart, TSubStepKicked;
            /* The start of the step is different if the time step is larger than the current */
            if (i < uRungLo) {
                TSubStepStart = stepStartTime + (substepsDoneAtThisSize)*substepSize * dDelta;
                TSubStepKicked = stepStartTime + (substepsDoneAtThisSize + 0.5) * substepSize * dDelta;
            }
            else {
                TSubStepStart = stepStartTime + (substepsDoneAtThisSize - 1.0) * substepSize * dDelta;
                TSubStepKicked = stepStartTime + (substepsDoneAtThisSize - 0.5) * substepSize * dDelta;
            }
            /* At the beginning we have a special case
            ** If we are not doing the closing kick, we are also in a special case.
            ** This only ever happens in simulate at the beginning and after a write
            ** where we do not have to undo a kick
            */
            if ((dTime == 0.0) || (!bKickClose)) {
                TSubStepStart = stepStartTime;
                TSubStepKicked = stepStartTime;
            }
            double dtPredISPHUndoOpen = TSubStepStart - TSubStepKicked;
            double dtPredISPHOpen = (dTime - TSubStepStart) / 2.0;
            double dtPredISPHClose = (dTime - TSubStepStart) / 2.0;
            if (csm->val.bComove) {
                kick->dtPredISPHUndoOpen[i] = csmComoveKickFac(csm, TSubStepKicked, dtPredISPHUndoOpen);
                kick->dtPredISPHOpen[i] = csmComoveKickFac(csm, TSubStepStart, dtPredISPHOpen);
                kick->dtPredISPHClose[i] = csmComoveKickFac(csm, TSubStepStart + dtPredISPHOpen, dtPredISPHClose);
            }
            else {
                kick->dtPredISPHUndoOpen[i] = dtPredISPHUndoOpen;
                kick->dtPredISPHOpen[i] = dtPredISPHOpen;
                kick->dtPredISPHClose[i] = dtPredISPHClose;
            }
        }
    }

    if (SPHoptions.doDensity || SPHoptions.doDensityCorrection) {
        uRungLo = uRungLoTemp;
    }
}

void MSR::UpdateGasValues(uint8_t uRungLo, double dTime, double dDelta, double dStep, int bKickClose, int bKickOpen, SPHOptions SPHoptions) {
    ServiceUpdateGasValues::input in;
    in.SPHoptions = SPHoptions;
    double sec, dsec;
    sec = MSR::Time();
    print("Update Gas Values ...\n");

    CalculateKickParameters(&in.kick, uRungLo, dTime, dDelta, dStep, bKickClose, bKickOpen, SPHoptions);

    mdl->RunService<ServiceUpdateGasValues>(in);
    dsec = MSR::Time() - sec;
    print("Gas Values updated, Wallclock: {:.5f} secs\n\n", dsec);
}

void MSR::TreeUpdateFlagBounds(int bNeedEwald, uint32_t uRoot, uint32_t utRoot, SPHOptions SPHoptions) {
    ServiceTreeUpdateFlagBounds::input in;
    const double ddHonHLimit = parameters.get_ddHonHLimit();
    LCL* plcl = pst->lcl();
    PKD pkd = plcl->pkd;
    double sec, dsec;

    auto nTopTree = pkd->NodeSize() * (2 * nThreads - 1);
    auto nMsgSize = sizeof(ServiceDistribTopTree::input) + nTopTree;

    auto buffer = std::make_unique<std::byte[]>(nMsgSize);
    auto inDistribTop = std::span(buffer.get(), nMsgSize);
    auto pDistribTop = new (buffer.get()) ServiceDistribTopTree::input;
    auto topTree = std::span(reinterpret_cast<std::byte*>(pDistribTop + 1), nTopTree);
    auto pkdn = pkd->tree[reinterpret_cast<KDN*>(pDistribTop + 1)];

    print("Update local trees...\n\n");

    pDistribTop->uRoot = uRoot;
    pDistribTop->allocateMemory = 0;

    in.base.nBucket = parameters.get_nBucket();
    in.base.nGroup = parameters.get_nGroup();
    in.base.uRoot = uRoot;
    in.base.utRoot = utRoot;
    in.base.ddHonHLimit = ddHonHLimit;
    in.SPHoptions = SPHoptions;
    sec = MSR::Time();
    nTopTree = mdl->RunService<ServiceTreeUpdateFlagBounds>(in, topTree);
    pDistribTop->nTop = nTopTree / pkd->NodeSize();
    assert(pDistribTop->nTop == (2 * nThreads - 1));
    mdl->RunService<ServiceDistribTopTree>(inDistribTop);
    dsec = MSR::Time() - sec;
    print("Tree updated, Wallclock: {:.5f} secs\n\n", dsec);
}

auto MSR::CountDistance(double dRadius2Outer, double dRadius2Inner) -> uint64_t {
    ServiceCountDistance::input in{.dRadius2Inner = dRadius2Inner, .dRadius2Outer = dRadius2Outer};
    ServiceCountDistance::output out;
    auto nOut = mdl->RunService<ServiceCountDistance>(in, out);
    assert(nOut == sizeof(out));
    return out;
}

double MSR::countSphere(double r, void* vctx) {
    auto* ctx = reinterpret_cast<SPHERECTX*>(vctx);
    ctx->nSelected = ctx->msr->CountDistance(r * r);
    return 1.0 * ctx->nSelected - 1.0 * ctx->nTarget;
}

void MSR::profileRootFind(double* dBins, int lo, int hi, int nAccuracy, SPHERECTX* ctx) {
    int nIter;
    int iBin = (lo + hi) / 2;
    if (lo == iBin) return;

    ctx->nTarget = ((ctx->nTotal - ctx->nInner) * ctx->dFrac * iBin + ctx->nInner);
    dBins[iBin] = illinois(countSphere, ctx, dBins[lo], dBins[hi], 0.0, 1.0 * nAccuracy, &nIter);
    profileRootFind(dBins, lo, iBin, nAccuracy, ctx);
    profileRootFind(dBins, iBin, hi, nAccuracy, ctx);
}

double MSR::countShell(double rInner, void* vctx) {
    auto* ctx = reinterpret_cast<SHELLCTX*>(vctx);
    double rOuter;
    local_t nSelected;

    if (rInner == ctx->rMiddle) nSelected = 0;
    else {
        rOuter = pow(10, 2.0 * log10(ctx->rMiddle) - log10(rInner));
        nSelected = ctx->msr->CountDistance(rOuter * rOuter, rInner * rInner);
    }
    return 1.0 * nSelected - 1.0 * ctx->nTarget;
}

/*
** Calculate a profile.
** Bins are of equal size (same number of particles) between dMinRadius and dLogRadius.
** From dLogRadius to dMaxRadius, the binning is done logarithmicly.
** Setting dLogRadius to dMinRadius results in purely logarithmic binning, while
** setting dLogRadius to dMaxRadius results in purely equal sized binning.
*/
void MSR::Profile(const PROFILEBIN** ppBins, int* pnBins, blitz::TinyVector<double, 3> r, double dMinRadius, double dLogRadius, double dMaxRadius, int nPerBin,
                  int nBins, int nAccuracy) {
    SPHERECTX ctxSphere;
    SHELLCTX ctxShell;
    PROFILEBIN* pBins;
    double sec, dsec;
    blitz::TinyVector<double, 3> com, vcm, L;
    double M;
    size_t inSize;
    int i, j;
    int nBinsInner;
    total_t N, n;
    LCL* plcl = pst->lcl();

    assert(dMinRadius <= dLogRadius);
    assert(dLogRadius <= dMaxRadius);
    assert(dLogRadius == dMinRadius || nPerBin > 0);
    assert(dLogRadius == dMaxRadius || nBins > 0);

    if (dLogRadius == dMaxRadius) nBins = 0;

    CalcDistance(r, dMaxRadius);
    CalcCOM(r, dMaxRadius, com, vcm, L, M);

    if (dLogRadius > dMinRadius) {
        /*
        ** The inner radius is calculated such that the logarithmic mid-point
        ** falls on dMinRadius.  This is done so that the profile is plotted
        ** all the way to the inner radius.  The correct radius must be between
        ** dMinRadius and the logrithmic difference between dMinRadius and
        ** dMaxRadius below dMinRadius.
        */
        ctxShell.rMiddle = dMinRadius;
        ctxShell.nTarget = nPerBin;
        ctxShell.msr = this;
        dMinRadius = illinois(countShell, &ctxShell, pow(10, 2.0 * log10(dMinRadius) - log10(dMaxRadius)), dMinRadius, 0.0, 0.0, NULL);
        N = CountDistance(dLogRadius * dLogRadius, dMinRadius * dMinRadius);
        nBinsInner = (N + nPerBin / 2) / nPerBin;
    }
    else {
        double dOuter;

        nBinsInner = 0;

        /*
        ** Calculate the logarithmic mid-point and verify that there are enough particles
        ** in the first bin.  If not, invoke the root finder.
        */
        ctxShell.rMiddle = dMinRadius;
        ctxShell.nTarget = nPerBin;
        ctxShell.msr = this;
        dMinRadius = pow(10, (2.0 * (nBins + 1) * log10(dMinRadius) - log10(dMaxRadius)) / (2 * nBins));
        dOuter = pow(10, 2.0 * log10(ctxShell.rMiddle) - log10(dMinRadius));
        N = CountDistance(dOuter * dOuter, dMinRadius * dMinRadius);
        if (N < nPerBin - nAccuracy) {
            dMinRadius = illinois(countShell, &ctxShell, pow(10, 2.0 * log10(dMinRadius) - log10(dMaxRadius)), dMinRadius, 0.0, 0.0, NULL);
        }
        dLogRadius = dMinRadius;
    }

    inSize = sizeof(ServiceProfile::input) - sizeof(ServiceProfile::input::dRadii[0]) *
                                                 (sizeof(ServiceProfile::input::dRadii) / sizeof(ServiceProfile::input::dRadii[0]) - nBins - nBinsInner - 1);
    auto buffer = std::make_unique<std::byte[]>(inSize);
    auto in = new (buffer.get()) ServiceProfile::input;

    in->dRadii[0] = dMinRadius;

    /*
    ** Inner, fixed size bins
    */
    if (nBinsInner) {
        sec = Time();
        print_detail("Root finding for {bins} bins\n", "bins"_a = nBinsInner);
        ctxSphere.nTotal = CountDistance(dLogRadius * dLogRadius);
        ctxSphere.nInner = CountDistance(dMinRadius * dMinRadius);
        ctxSphere.msr = this;
        ctxSphere.dFrac = 1.0 / nBinsInner;
        in->dRadii[nBinsInner] = dLogRadius;
        profileRootFind(in->dRadii.data(), 0, nBinsInner, nAccuracy, &ctxSphere);
        dsec = Time() - sec;
        print_detail("Root finding complete, Wallclock: {seconds:.5} secs\n", "seconds"_a = dsec);
    }

    /*
    ** Now logarithmic binning for the outer region.  We still obey nPerBin
    ** as the minimum number of particles to include in each bin.
    */
    if (nBins) {
        double dLogMin;
        double dLogMax = log10(dMaxRadius);
        double dRadius;

        ctxSphere.nTotal = SelSphere(r, dMaxRadius, 1, 1);
        ctxSphere.msr = this;

        N = CountDistance(0.0, dLogRadius * dLogRadius);
        for (i = 1; i < nBins; i++) {
            int nBinsRem = nBins - i + 1;

            dLogMin = log10(in->dRadii[nBinsInner + i - 1]);
            dRadius = pow(10, (dLogMax - dLogMin) / nBinsRem + dLogMin);
            n = CountDistance(dRadius * dRadius);
            if (n - N < nPerBin - nAccuracy) {
                ctxSphere.nTarget = N + nPerBin;
                dRadius = illinois(countSphere, &ctxSphere, 0.0, dMaxRadius, 0.0, 1.0 * nAccuracy, NULL);
                n = ctxSphere.nSelected;
            }
            in->dRadii[nBinsInner + i] = dRadius;
            N = n;
        }
    }

    nBins = nBins + nBinsInner;

    in->dRadii[nBins] = dMaxRadius;

    sec = Time();
    print_detail("Profiling\n");
    in->dCenter = r;
    in->com = com;
    in->vcm = vcm;
    in->L = L;
    in->nBins = nBins + 1;
    in->uRungLo = 0;
    in->uRungHi = MaxRung() - 1;
    mdl->RunService<ServiceProfile>(std::span(buffer.get(), inSize));

    /*
    ** Finalize bin values
    */
    pBins = plcl->pkd->profileBins;
    for (i = 0; i < nBins + 1; i++) {
        if (pBins[i].dMassInBin > 0.0) {
            pBins[i].vel_radial /= pBins[i].dMassInBin;
            pBins[i].vel_radial_sigma /= pBins[i].dMassInBin;
            pBins[i].vel_tang_sigma = sqrt(pBins[i].vel_tang_sigma / pBins[i].dMassInBin);
            if (pBins[i].vel_radial_sigma > pBins[i].vel_radial * pBins[i].vel_radial)
                pBins[i].vel_radial_sigma = sqrt(pBins[i].vel_radial_sigma - pBins[i].vel_radial * pBins[i].vel_radial);
            else pBins[i].vel_radial_sigma = 0.0;
            for (j = 0; j < 3; j++) {
                pBins[i].L[j] /= pBins[i].dMassInBin;
            }
        }
    }

    dsec = Time() - sec;
    print_detail("Profiling complete, Wallclock: {:.5} secs\n\n", dsec);

    if (ppBins) *ppBins = plcl->pkd->profileBins;
    if (pnBins) *pnBins = nBins + 1;
}

#ifdef MDL_FFTW
void MSR::GridCreateFFT(int nGrid) {
    ServiceGridCreateFFT::input in{.nGrid = nGrid};
    mdl->RunService<ServiceGridCreateFFT>(in);
}

void MSR::GridDeleteFFT() {
    mdl->RunService<ServiceGridDeleteFFT>();
}

/* Important: call msrGridCreateFFT() before, and msrGridDeleteFFT() after */
std::tuple<std::vector<uint64_t>, std::vector<float>, std::vector<float>, std::vector<float>> // nPk, fK, fPk, fPkAll
MSR::MeasurePk(int iAssignment, int bInterlace, int nGrid, double a, int nBins, int fold) {
    std::vector<uint64_t> nPk;
    std::vector<float> fK, fPk, fPkAll;
    double dsec;

    GridCreateFFT(nGrid);

    if (nGrid / 2 < nBins) nBins = nGrid / 2;
    assert(nBins <= ServiceGridBinK::max_k_bins);

    TimerStart(TIMER_NONE);
    print("Measuring P(k) with grid size {grid} ({bins} bins)...\n", "grid"_a = nGrid, "bins"_a = nBins);

    AssignMass(iAssignment, 0, 0.0, fold);
    DensityContrast(0);
    if (bInterlace) {
        AssignMass(iAssignment, 1, 0.5, fold);
        DensityContrast(1);
        Interlace(0, 1); // We no longer need grid 1
    }
    WindowCorrection(iAssignment, 0);

    std::tie(nPk, fK, fPk) = GridBinK(nBins, 0);
    if (csm->val.classData.bClass && parameters.get_nGridLin() > 0 && parameters.get_achPkSpecies().length() > 0) {
        AddLinearSignal(0, parameters.get_iSeed(), parameters.get_dBoxSize(), a, parameters.get_bFixedAmpIC(), parameters.get_dFixedAmpPhasePI() * M_PI);
        std::tie(nPk, fK, fPkAll) = GridBinK(nBins, 0);
    }
    else {
        fPkAll.resize(nBins);
    }

    GridDeleteFFT();

    TimerStop(TIMER_NONE);
    dsec = TimerGet(TIMER_NONE);
    print("P(k) Calculated, Wallclock: {:.5} secs\n\n", dsec);
    return {nPk, fK, fPk, fPkAll};
}

void MSR::AssignMass(int iAssignment, int iGrid, float fDelta, int fold) {
    static const char* schemes[] = {"Nearest Grid Point (NGP)", "Cloud in Cell (CIC)", "Triangular Shaped Cloud (TSC)", "Piecewise Cubic Spline (PCS)"};
    assert(iAssignment >= 0 && iAssignment <= 3);
    printf("Assigning mass using %s (order %d)\n", schemes[iAssignment], iAssignment);
    ServiceAssignMass::input mass{.iAssignment = iAssignment, .iGrid = iGrid, .fDelta = fDelta, .fold = fold};
    auto sec = MSR::Time();
    mdl->RunService<ServiceAssignMass>(mass);
    printf("Mass assignment complete, Wallclock: %f secs\n", MSR::Time() - sec);
}

void MSR::DensityContrast(int iGrid, bool k) {
    ServiceDensityContrast::input in{.dTotalMass = TotalMass(), .iGrid = iGrid, .k = k};
    mdl->RunService<ServiceDensityContrast>(in);
}

void MSR::Interlace(int iGridTarget, int iGridSource) {
    ServiceInterlace::input in{.iGridTarget = iGridTarget, .iGridSource = iGridSource};
    mdl->RunService<ServiceInterlace>(in);
}

void MSR::WindowCorrection(int iAssignment, int iGrid) {
    assert(iAssignment >= 0 && iAssignment <= 3);
    ServiceWindowCorrection::input in{.iAssignment = iAssignment, .iGrid = iGrid};
    mdl->RunService<ServiceWindowCorrection>(in);
}

void MSR::AddLinearSignal(int iGrid, int iSeed, double Lbox, double a, bool bFixed, float fPhase) {
    ServiceAddLinearSignal::input in{.iGrid = iGrid, .iSeed = iSeed, .bFixed = bFixed, .fPhase = fPhase, .Lbox = Lbox, .a = a};
    mdl->RunService<ServiceAddLinearSignal>(in);
}

std::tuple<std::vector<uint64_t>, std::vector<float>, std::vector<float>> // nPk, fK, fPk
MSR::GridBinK(int nBins, int iGrid) {
    std::vector<uint64_t> nPk(nBins);
    std::vector<float> fK(nBins), fPk(nBins);
    auto out = std::make_unique<ServiceGridBinK::output>();

    assert(sizeof(out->nPower) / sizeof(out->nPower[0]) > nBins);

    ServiceGridBinK::input in{.nBins = nBins, .iGrid = iGrid};
    mdl->RunService<ServiceGridBinK>(in, *out);
    for (int i = 0; i < nBins; i++) {
        if (out->nPower[i] == 0) fK[i] = fPk[i] = 0;
        else {
            nPk[i] = out->nPower[i];
            fK[i] = exp(out->fK[i] / out->nPower[i]);
            fPk[i] = out->fPower[i] / out->nPower[i];
        }
    }

    /* At this point, dPk[] needs to be corrected by the box size */
    return {nPk, fK, fPk};
}

std::tuple<std::vector<uint64_t>, std::vector<float>, std::vector<float>> // nPk, fK, fPk
MSR::MeasureLinearPk(int nGrid, double dA, double dBoxSize) {
    int i;
    double dsec;

    TimerStart(TIMER_NONE);

    struct ServiceMeasureLinearPk::input in{.dA = dA,
                                            .dBoxSize = dBoxSize,
                                            .fPhase = static_cast<float>(parameters.get_dFixedAmpPhasePI() * M_PI),
                                            .iSeed = static_cast<int>(parameters.get_iSeed()),
                                            .nGrid = nGrid,
                                            .nBins = nGrid / 2,
                                            .bFixed = parameters.get_bFixedAmpIC()};
    auto out = std::make_unique<ServiceMeasureLinearPk::output>();
    print("Measuring P_lin(k) with grid size {} ({} bins)...\n", in.nGrid, in.nBins);
    mdl->RunService<ServiceMeasureLinearPk>(in, *out);
    std::vector<uint64_t> nPk(nBins);
    std::vector<float> fK(nBins), fPk(nBins);
    for (i = 0; i < in.nBins; i++) {
        if (out->nPower[i] == 0) fK[i] = fPk[i] = 0;
        else {
            if (nPk.size()) nPk[i] = out->nPower[i];
            fK[i] = out->fK[i] / out->nPower[i];
            fPk[i] = out->fPower[i] / out->nPower[i];
        }
    }
    /* At this point, dPk[] needs to be corrected by the box size */

    TimerStop(TIMER_NONE);
    dsec = TimerGet(TIMER_NONE);
    print("P_lin(k) Calculated, Wallclock: {:.5} secs\n\n", dsec);
    return {nPk, fK, fPk};
}

void MSR::SetLinearGrid(double dTime, double dDelta, int nGrid, int bKickClose, int bKickOpen) {
    print("Setting force grids of linear species with nGridLin = {}\n", nGrid);
    double dsec;
    TimerStart(TIMER_NONE);

    ServiceSetLinearGrid::input in;
    in.nGrid = nGrid;

    int do_DeltaRho_lin_avg = 1;
    in.a0 = in.a1 = in.a = csmTime2Exp(csm, dTime);
    if (do_DeltaRho_lin_avg) {
        if (bKickClose) in.a0 = csmTime2Exp(csm, dTime - 0.5 * dDelta);
        if (bKickOpen) in.a1 = csmTime2Exp(csm, dTime + 0.5 * dDelta);
    }

    in.dBSize = parameters.get_dBoxSize();
    /* Parameters for the grid realization */
    in.iSeed = parameters.get_iSeed();
    in.bFixed = parameters.get_bFixedAmpIC();
    in.fPhase = parameters.get_dFixedAmpPhasePI() * M_PI;
    mdl->RunService<ServiceSetLinearGrid>(in);

    TimerStop(TIMER_NONE);
    dsec = TimerGet(TIMER_NONE);
    print("Force grids of linear species set, Wallclock: {:.5} secs\n\n", dsec);
}

void MSR::BispectrumSelect(int iGridTarget, int iGridSource, double kmin, double kmax) {
    ServiceBispectrumSelect::input in{.iGridTarget = iGridTarget, .iGridSource = iGridSource, .kmin = kmin, .kmax = kmax};
    mdl->RunService<ServiceBispectrumSelect>(in);
}

auto MSR::BispectrumCalculate(int iGrid1, int iGrid2, int iGrid3) -> double {
    ServiceBispectrumCalculate::input in{.iGrid1 = iGrid1, .iGrid2 = iGrid2, .iGrid3 = iGrid3};
    ServiceBispectrumCalculate::output result;
    mdl->RunService<ServiceBispectrumCalculate>(in, result);
    return result;
}

/* First call SetLinearGrid() to setup the grid */
void MSR::LinearKick(double dTime, double dDelta, int bKickClose, int bKickOpen) {
    double dt = 0.5 * dDelta;
    double dsec;

    print("Applying Linear Kick...\n");
    TimerStart(TIMER_NONE);
    ServiceLinearKick::input in{.dtOpen = 0.0, .dtClose = 0.0};
    if (csm->val.bComove) {
        if (bKickClose) in.dtClose = csmComoveKickFac(csm, dTime - dt, dt);
        if (bKickOpen) in.dtOpen = csmComoveKickFac(csm, dTime, dt);
    }
    else {
        if (bKickClose) in.dtClose = dt;
        if (bKickOpen) in.dtOpen = dt;
    }
    mdl->RunService<ServiceLinearKick>(in);
    TimerStop(TIMER_NONE);
    dsec = TimerGet(TIMER_NONE);
    print("Linear Kick Applied, Wallclock: {:.5} secs\n\n", dsec);
}
#endif

void MSR::OutputOrbits(int iStep, double dTime) {
    int i;

    if (parameters.has_lstOrbits()) {
        double dExp, dvFac;

        if (csm->val.bComove) {
            dExp = csmTime2Exp(csm, dTime);
            dvFac = 1.0 / (dExp * dExp);
        }
        else {
            dExp = dTime;
            dvFac = 1.0;
        }

        // Sort the particle IDs to get them in order so we can binary search later
        auto particle_ids = parameters.get_lstOrbits();
        std::ranges::sort(particle_ids);
        // There should be exactly one of each ID, so remove duplicates and resize output vector
        auto last = std::ranges::unique(particle_ids);
        particle_ids.erase(last.begin(), last.end());
        std::vector<ServiceGetOrbits::output> particles(particle_ids.size());
        // Retrieve the particle data and sort by ID
        auto n = mdl->RunService<ServiceGetOrbits>(particle_ids, particles);
        particles.resize(n / sizeof(ServiceGetOrbits::output));
        std::ranges::sort(particles, {}, &ServiceGetOrbits::output::id);

        auto filename = BuildName(iStep, ".orb");
        std::ofstream fs(filename);
        if (fs.fail()) {
            std::cerr << "Could not create orbit file: " << filename << std::endl;
            perror(filename.c_str());
            Exit(errno);
        }
        fmt::print(fs, "{n} {a}\n", "n"_a = particle_ids.size(), "a"_a = dExp);
        for (i = 0; i < particle_ids.size(); ++i) {
            fmt::print(fs, "{id} {mass:.8e} {x:.16e} {y:.16e} {z:.16e} {vx:.8e} {vy:.8e} {vz:.8e} {phi:.8e}\n", "id"_a = particles[i].id,
                       "mass"_a = particles[i].mass, "phi"_a = particles[i].phi, "x"_a = particles[i].r[0], "y"_a = particles[i].r[1],
                       "z"_a = particles[i].r[2], "vx"_a = particles[i].v[0] * dvFac, "vy"_a = particles[i].v[1] * dvFac, "vz"_a = particles[i].v[2] * dvFac);
        }
        fs.close();
    }
}

void MSR::RsHaloLoadIds(const std::string_view filename_template, bool bAppend) {
#ifdef HAVE_ROCKSTAR
    std::vector<uint64_t> counts;

    TimerStart(TIMER_NONE);
    print("Scanning Rockstar halo binary files...\n");
    ServiceRsHaloCount::input hdr;
    marshal_string(hdr.filename, filename_template);
    hdr.nSimultaneous = hdr.nTotalActive = parallel_read_count();
    hdr.iReaderWriter = 0;
    hdr.nElementSize = 1;
    auto out = new ServiceFileSizes::output[ServiceFileSizes::max_files];
    auto out_span = std::span(out, ServiceFileSizes::max_files);
    auto n = mdl->RunService<ServiceRsHaloCount>(hdr, out_span);
    // auto n = mdl->RunService<ServiceRsHaloCount>(hdr, std::span(out, ServiceFileSizes::max_files));
    n /= sizeof(ServiceFileSizes::output);
    counts.resize(n);
    for (auto i = 0; i < n; ++i) {
        assert(out[i].iFileIndex < n);
        counts[out[i].iFileIndex] = out[i].nFileBytes;
    }
    delete[] out;

    if (counts.size() == 0) {
        perror(filename_template.data());
        abort();
    }
    TimerStop(TIMER_NONE);
    auto dsec = TimerGet(TIMER_NONE);

    print("... identified {} halos in {} files, Wallclock: {:.5} secs\n", std::accumulate(counts.begin(), counts.end(), uint64_t(0)), counts.size(), dsec);
    RsLoadIds(ServiceRsHaloLoadIds::name, counts, filename_template, bAppend);
    TimerStop(TIMER_IO);
#endif
}

void MSR::RsLoadIds(std::string_view service_name, std::vector<uint64_t>& counts, const std::string_view filename_template, bool bAppend) {
    using mdl::ServiceBuffer;
    ServiceBuffer msg{ServiceBuffer::Field<ServiceInput::input>(), ServiceBuffer::Field<ServiceRsLoadIds::io_elements>(counts.size()),
                      ServiceBuffer::Field<ServiceRsLoadIds::input>()};
    auto hdr = static_cast<ServiceInput::input*>(msg.data(0));
    auto elements = static_cast<ServiceRsLoadIds::io_elements*>(msg.data(1));
    auto in = static_cast<ServiceRsLoadIds::input*>(msg.data(2));
    in->bAppend = bAppend;
    hdr->nFiles = counts.size();
    hdr->nElements = std::accumulate(counts.begin(), counts.end(), uint64_t(0));
    std::copy(counts.begin(), counts.end(), elements);
    hdr->iBeg = 0;
    hdr->iEnd = hdr->nElements;
    marshal_string(hdr->io.filename, filename_template);
    hdr->io.nSimultaneous = parallel_read_count();
    hdr->io.nSegment = hdr->io.iThread = 0; // setup later
    hdr->io.iReaderWriter = 0;
    print("Loading {} particle IDs from {} files\n", hdr->nElements, hdr->nFiles);
    TimerStart(TIMER_NONE);
    mdl->RunService(service_name, msg);
    TimerStop(TIMER_NONE);
    auto dsec = TimerGet(TIMER_NONE);
    print("... finished reading particles IDs, Wallclock: {:.5} secs\n", dsec);
}

void MSR::RsLoadIds(const std::string_view filename_template, bool bAppend) {
    std::vector<uint64_t> counts;
    TimerStart(TIMER_NONE);
    print("Scanning Particle ID binary files...\n");
    stat_files(counts, filename_template, sizeof(uint64_t));
    TimerStop(TIMER_NONE);
    auto dsec = TimerGet(TIMER_NONE);
    print("... identified {} IDs in {} files, Wallclock: {:.5} secs\n", std::accumulate(counts.begin(), counts.end(), uint64_t(0)), counts.size(), dsec);
    RsLoadIds(ServiceRsLoadIds::name, counts, filename_template, bAppend);
}

void MSR::RsSaveIds(const std::string_view filename_template) {
    ServiceRsSaveIds::input hdr;
    marshal_string(hdr.io.filename, filename_template);
    hdr.io.nSimultaneous = parallel_write_count();
    hdr.io.nSegment = hdr.io.iThread = 0; // setup later
    hdr.io.iReaderWriter = 0;
    print("Saving particle IDS\n");
    TimerStart(TIMER_NONE);
    mdl->RunService<ServiceRsSaveIds>(hdr);
    TimerStop(TIMER_NONE);
    auto dsec = TimerGet(TIMER_NONE);
    print("... finished writing particles IDs, Wallclock: {:.5} secs\n", dsec);
}

void MSR::RsReorderIds() {
    using mdl::ServiceBuffer;
    ServiceBuffer msg{ServiceBuffer::Field<ServiceRsExtract::input>(mdl->Threads() + 1)};
    auto pOrds = static_cast<ServiceRsExtract::input*>(msg.data(0));
    pOrds[mdl->Threads()] = N;
    print("Reordering particle IDS\n");
    TimerStart(TIMER_NONE);
    mdl->RunService<ServiceGetOrdSplits>({}, std::span(reinterpret_cast<std::byte*>(pOrds), msg.size(0)));
    mdl->RunService<ServiceRsReorderIds>(msg);
    TimerStop(TIMER_NONE);
    auto dsec = TimerGet(TIMER_NONE);
    print("... finished reordering particles IDs, Wallclock: {:.5} secs\n", dsec);
}

void MSR::RsExtract(std::string_view filename_template) {
    print("Extracting matching particles\n");
    using mdl::ServiceBuffer;
    ServiceBuffer msg{ServiceBuffer::Field<ServiceRsExtract::header>(), ServiceBuffer::Field<ServiceRsExtract::input>(mdl->Threads() + 1)};
    auto hdr = static_cast<ServiceRsExtract::header*>(msg.data(0));
    marshal_string(hdr->io.filename, filename_template);
    hdr->io.nSimultaneous = parallel_write_count();
    hdr->io.nSegment = hdr->io.iThread = 0; // setup later
    hdr->io.iReaderWriter = 0;

    auto pOrds = static_cast<ServiceRsExtract::input*>(msg.data(1));
    pOrds[mdl->Threads()] = N;
    TimerStart(TIMER_NONE);
    mdl->RunService<ServiceGetOrdSplits>({}, std::span(reinterpret_cast<std::byte*>(pOrds), msg.size(1)));
    mdl->RunService<ServiceRsExtract>(msg);
    TimerStop(TIMER_NONE);
    auto dsec = TimerGet(TIMER_NONE);
    print("... finished extracting particles, Wallclock: {:.5} secs\n", dsec);
}

void MSR::SetEOSParam() {
    const double dHydFrac = parameters.get_dInitialH();
    const double dnHToRho = MHYDR / dHydFrac / units.dGmPerCcUnit;
    calc.dEOSFloorDen = parameters.get_dEOSFloornH() * dnHToRho;
    calc.dEOSFlooru = parameters.get_dEOSFloorTemp() * dTuFacPrimNeutral;
    if (csm->val.bComove) calc.dEOSFloorMinBaryonOD = parameters.get_dEOSFloorMinOD() * csm->val.dOmegab;
    else calc.dEOSFloorMinBaryonOD = 0.;
    calc.dEOSPolyFloorExponent = parameters.get_dEOSPolyFloorIndex() - 1.;
    calc.dEOSPolyFloorDen = parameters.get_dEOSPolyFloornH() * dnHToRho;
    calc.dEOSPolyFlooru = parameters.get_dEOSPolyFloorTemp() * dTuFacPrimNeutral;
    if (csm->val.bComove) calc.dEOSPolyFloorMinBaryonOD = parameters.get_dEOSPolyFloorMinOD() * csm->val.dOmegab;
    else calc.dEOSPolyFloorMinBaryonOD = 0.;
}

auto MSR::ValidateEOSParam() -> int {
    if (!parameters.has_dOmegab() && parameters.has_dEOSFloorMinOD()) {
        fprintf(stderr, "ERROR: dEOSFloorMinOD is specified but dOmegab is not set\n");
        return 0;
    }
    if (!parameters.has_dOmegab() && parameters.has_dEOSPolyFloorMinOD()) {
        fprintf(stderr, "ERROR: dEOSPolyFloorMinOD is specified but dOmegab is not set\n");
        return 0;
    }
    return 1;
}

void MSR::BHAccretion(double dTime) {
    ServiceBHAccretion::input in{.dScaleFactor = csmTime2Exp(csm, dTime)};
    mdl->RunService<ServiceBHAccretion>(in);
}

void MSR::BHGasPin(double dTime, double dDelta) {
    ServiceBHGasPin::input in{
        .nSmooth = static_cast<int>(parameters.get_nSmooth()),
        .bPeriodic = parameters.get_bPeriodic(),
        .a = csmTime2Exp(csm, dTime),
    };
    mdl->RunService<ServiceBHGasPin>(in);
}

void MSR::BHReposition() {
    mdl->RunService<ServiceBHReposition>();
}

void MSR::BHEvolve(double dTime, double dDelta) {
    ServiceBHEvolve::input in{
        .nSmooth = static_cast<int>(parameters.get_nSmooth()),
        .bPeriodic = parameters.get_bPeriodic(),
        .a = csmTime2Exp(csm, dTime),
        .dBHAccretionAlpha = parameters.get_dBHAccretionAlpha(),
        .dBHAccretionCvisc = parameters.get_dBHAccretionCvisc(),
        .dBHAccretionEddFac = calc.dBHAccretionEddFac,
        .dConstGamma = parameters.get_dConstGamma(),
        .dBHFBEff = calc.dBHFBEff,
        .dBHFBEcrit = calc.dBHFBEcrit,
        .bBHAccretion = parameters.get_bBHAccretion(),
        .bBHFeedback = parameters.get_bBHFeedback(),
    };
    mdl->RunService<ServiceBHEvolve>(in);
}

void MSR::SetBlackholeParam() {
    calc.dBHAccretionEddFac = parameters.get_dBHAccretionEddFac() * (1e3 / MSOLG / units.dMsolUnit) / pow(100. / KPCCM / units.dKpcUnit, 3) / units.dSecUnit /
                              parameters.get_dBHRadiativeEff();
    // We precompute the factor such that we only need to multiply
    // AccretionRate by this amount to get E_feed
    calc.dBHFBEff = parameters.get_dBHFBEff() * parameters.get_dBHRadiativeEff() * pow(LIGHTSPEED * 1e-5 / units.dKmPerSecUnit, 2);

    // This, in principle, will not be a parameter
    double n_heat = 1.0;

    // We convert from Delta T to energy per mass.
    // This needs to be multiplied by the mass of the gas particle
    calc.dBHFBEcrit = parameters.get_dBHFBDT() * dTuFacPrimIonised * n_heat;
}

int MSR::ValidateBlackholeParam() {
    if (parameters.get_bBHAccretion()) {
        if (parameters.get_dBHAccretionAlpha() <= 0) {
            fprintf(stderr, "ERROR: dBHAccretionAlpha should be positive."
                            "If you want to avoid boosting the Bondi accretion rate,"
                            "just set dBHAccretionAlpha=1.0\n");
            return 0;
        }
    }
    if (parameters.get_bBHPlaceSeed() && !parameters.get_bFindGroups()) {
        parameters.set_bFindGroups(true);
        fprintf(stderr, "WARNING: Blackhole seeding requires Friends-of-friends group finding. "
                        "Setting bFindGroups to true\n");
    }
    if (!parameters.get_bMemPotential()) {
        parameters.set_bMemPotential(true);
        fprintf(stderr, "WARNING: The blackhole module requires bMemPotential. "
                        "Setting bMemPotential to true\n");
    }
    if (!parameters.get_bMemMass()) {
        parameters.set_bMemMass(true);
        fprintf(stderr, "WARNING: The blackhole module requires bMemMass. "
                        "Setting bMemMass to true\n");
    }
    if ((parameters.get_bBHMerger() || parameters.get_bBHAccretion()) && !parameters.get_bAddDelete()) {
        parameters.set_bAddDelete(true);
        fprintf(stderr, "WARNING: Blackhole mergers and gas accretion require bAddDelete. "
                        "Setting bAddDelete to true\n");
    }
    return 1;
}

void MSR::BlackholeInit(uint8_t uRungMax) {
    struct ServiceBHInit::input in{.uRungMax = uRungMax};
    mdl->RunService<ServiceBHInit>(in);
}

void MSR::PlaceBHSeed(double dTime, uint8_t uRungMax) {
    assert(parameters.get_bFindGroups());
    printf("Planting seeds...\n");
    TimerStart(TIMER_BHS);
    ServicePlaceBHSeed::input in{
        .dTime = dTime,
        .dScaleFactor = csmTime2Exp(csm, dTime),
        .dDenMin = parameters.get_bDoStarFormation() ? calc.dSFThresholdDen * pow(csmTime2Exp(csm, dTime), 3) : 0.0,
        .dBHMhaloMin = parameters.get_dBHMhaloMin(),
        .dTau = parameters.get_dTau(),
        .dBHSeedMass = parameters.get_dBHSeedMass(),
        .uRungMax = uRungMax,
    };
    ServicePlaceBHSeed::output nBH;
    mdl->RunService<ServicePlaceBHSeed>(in, nBH);
    TimerStop(TIMER_BHS);

    nGas -= nBH;
    nBH += nBH;
    printf("Planted %lu BH seeds \n", nBH);
#ifdef OPTIM_REORDER_IN_NODES
    ReorderInNodes();
#endif
}

void MSR::BHMerger(double dTime) {
    TimerStart(TIMER_BHS);

#ifdef DEBUG_BH_ONLY
    // This was not performed in the main loop because there is no gas
    ReorderInNodes();
#endif

    ServiceBHMerger::input inBHMerger{
        .nSmooth = static_cast<int>(parameters.get_nSmooth()),
        .bPeriodic = parameters.get_bPeriodic(),
    };
    ServiceBHMerger::output outBHMerger;
    mdl->RunService<ServiceBHMerger>(inBHMerger, outBHMerger);
    TimerStop(TIMER_BHS);
}

void MSR::BHStep(double dTime, double dDelta) {
#ifndef DEBUG_BH_ONLY
    ServiceBHStep::input in{
        .nSmooth = static_cast<int>(parameters.get_nSmooth()),
        .bPeriodic = parameters.get_bPeriodic(),
    };
    mdl->RunService<ServiceBHStep>(in);
#endif
}

void MSR::ComputeSmoothing(double dTime, double dDelta) {
    uint64_t nSmoothed;
    int it = 0, maxit = 100;

    printf("Computing density... \n");
    TimerStart(TIMER_DENSITY);
    ServiceHydroDensity::input in{
        .dNeighborsStd = parameters.get_dNeighborsStd(),
        .dhMinOverSoft = parameters.get_dhMinOverSoft(),
        .nSmooth = static_cast<int>(parameters.get_nSmooth()),
        .bPeriodic = parameters.get_bPeriodic(),
        .bIterativeSmoothingLength = parameters.get_bIterativeSmoothingLength(),
    };
    if (in.bIterativeSmoothingLength) {
#ifdef OPTIM_AVOID_IS_ACTIVE
        SelActives();
#else
        SelAll(); // We set all particles as "not converged"
#endif
        do {
            mdl->RunService<ServiceHydroDensity>(in, nSmoothed);
        } while (nSmoothed > 0 && ++it <= maxit);
        if (nSmoothed > 0) {
            // If we can not converge, at least make sure that we update the
            // density and other variables with the latest fBall
            ServiceHydroDensityFinal::input inFinal{.nSmooth = static_cast<int>(parameters.get_nSmooth()), .bPeriodic = parameters.get_bPeriodic()};
            mdl->RunService<ServiceHydroDensityFinal>(inFinal, nSmoothed);
            printf("Smoothing length did not converge for %" PRIu64 " particles\n", nSmoothed);
        }

        TimerStop(TIMER_DENSITY);
        double dsec = TimerGet(TIMER_DENSITY);
        printf("Computing h took %d iterations and %.5f seconds \n", it, dsec);
    }
    else {
        mdl->RunService<ServiceHydroDensity>(in, nSmoothed);
    }
}

void MSR::MeshlessFluxes(double dTime, double dDelta) {
    double dsec;
    printf("Computing fluxes... ");
    TimerStart(TIMER_FLUXES);
#ifdef OPTIM_AVOID_IS_ACTIVE
    SelActives();
#endif
    ServiceHydroFluxes::input in{
        .nSmooth = static_cast<int>(parameters.get_nSmooth()),
        .bPeriodic = parameters.get_bPeriodic(),
        .dTime = dTime,
        .dDelta = dDelta,
        .a = Comove() ? csmTime2Exp(csm, dTime) : 1.0,
        .H = Comove() ? csmTime2Hub(csm, dTime) : 0.0,
        .dConstGamma = parameters.get_dConstGamma(),
        .eEOS = get_eEOS(),
    };
    uint64_t nSmoothed;
    mdl->RunService<ServiceHydroFluxes>(in, nSmoothed);

    TimerStop(TIMER_FLUXES);
    dsec = TimerGet(TIMER_FLUXES);
    printf("took %.5f seconds\n", dsec);
}

void MSR::MeshlessGradients(double dTime, double dDelta) {
    double dsec;
    printf("Computing gradients... ");

    TimerStart(TIMER_GRADIENTS);
#ifdef OPTIM_AVOID_IS_ACTIVE
    SelActives();
#endif
    ServiceHydroGradients::input in{
        .nSmooth = static_cast<int>(parameters.get_nSmooth()),
        .bPeriodic = parameters.get_bPeriodic(),
        .bConditionedBarthJespersen = parameters.get_bConditionedBarthJespersen(),
        .bStricterSlopeLimiter = parameters.get_bStricterSlopeLimiter(),
    };
    uint64_t nSmoothed;
    mdl->RunService<ServiceHydroGradients>(in, nSmoothed);

    TimerStop(TIMER_GRADIENTS);
    dsec = TimerGet(TIMER_GRADIENTS);
    printf("took %.5f seconds\n", dsec);
}

void MSR::HydroStep(double dTime, double dDelta) {
    double dsec;

    printf("Computing hydro time step... ");

    TimerStart(TIMER_TIMESTEP);
#ifdef OPTIM_AVOID_IS_ACTIVE
    SelActives();
#endif
    ServiceHydroStep::input in{
        .nSmooth = static_cast<int>(parameters.get_nSmooth()),
        .bPeriodic = parameters.get_bPeriodic(),
        .bDoGravity = parameters.get_bDoGravity(),
        .dEtaCourant = parameters.get_dEtaCourant(),
        .a = Comove() ? csmTime2Exp(csm, dTime) : 1.0,
        .dCFLacc = parameters.get_dCFLacc(),
        .dDelta = dDelta,
    };
    uint64_t nSmoothed;
    mdl->RunService<ServiceHydroStep>(in, nSmoothed);

    if (parameters.get_bGlobalDt()) {
        SetGlobalRung(GetMaxRung());
    }

    TimerStop(TIMER_TIMESTEP);
    dsec = TimerGet(TIMER_TIMESTEP);
    printf("took %.5f seconds\n", dsec);

    if (parameters.get_bWakeUpParticles()) {
        ServiceWakeParticles::input in{dTime};
        mdl->RunService<ServiceWakeParticles>(in);
    }
}

void MSR::SetFeedbackParam() {
    calc.dSNFBDu = parameters.get_dSNFBDT() * dTuFacPrimIonised;

    auto IMF = ChooseIMF(parameters.get_achIMFType().data(), parameters.get_dIMFMinMass(), parameters.get_dIMFMaxMass());
    const double dCCSNNumPerMass = IMF->UnweightedIntegration(parameters.get_dCCSNMinMass(), parameters.get_dCCSNMaxMass());
    calc.dCCSNFBSpecEnergy = (parameters.get_dCCSNEnergy() / MSOLG) * dCCSNNumPerMass / units.dErgPerGmUnit;

    calc.dSNIaFBSpecEnergy = (parameters.get_dSNIaEnergy() / MSOLG) * parameters.get_dSNIaNumPerMass() / units.dErgPerGmUnit;

    calc.dCCSNFBDelay = parameters.get_dCCSNFBDelay() * SECONDSPERYEAR / units.dSecUnit;
    calc.dSNIaFBDelay = parameters.get_dSNIaFBDelay() * SECONDSPERYEAR / units.dSecUnit;

    const double dnHToRho = MHYDR / parameters.get_dInitialH() / units.dGmPerCcUnit;
    calc.dSNFBEffnH0 = parameters.get_dSNFBEffnH0() * dnHToRho;
}

#ifdef GRACKLE
void MSR::GrackleInit(int bComove, double dScaleFactor) {
    ServiceGrackleInit::input in;
    in.bComove = bComove;
    in.dScaleFactor = dScaleFactor;
    marshal_string(in.achCoolingTable, parameters.get_achCoolingTables());
    in.units = units;
    mdl->RunService<ServiceGrackleInit>(in);
}
#endif

#ifdef STELLAR_EVOLUTION
void MSR::SetStellarEvolutionParam() {
    const double dYrToTime = SECONDSPERYEAR / units.dSecUnit;
    auto achSNIaDTDType = parameters.get_achSNIaDTDType();

    if (achSNIaDTDType == "exponential") {
        calc.dSNIaNorm = parameters.get_dSNIaNumPerMass();
        calc.dSNIaScale = parameters.get_dSNIaExpScale() * dYrToTime;
    }
    else if (achSNIaDTDType == "powerlaw") {
        auto dSNIaPLScale = parameters.get_dSNIaPLScale();
        calc.dSNIaNorm = parameters.get_dSNIaNumPerMass() / (pow(parameters.get_dSNIaPLFinalTime() * dYrToTime, dSNIaPLScale + 1.0) -
                                                             pow(parameters.get_dSNIaPLInitTime() * dYrToTime, dSNIaPLScale + 1.0));
        calc.dSNIaScale = dSNIaPLScale;
    }
    else {
        std::cerr << "ERROR: Undefined SNIa DTD type has been given in " << "achSNIaDTDType parameter: " << achSNIaDTDType << std::endl;
        Exit(1);
    }

    calc.dWindSpecificEkin = 0.5 * pow(parameters.get_dStellarWindSpeed() / units.dKmPerSecUnit, 2);
}

void MSR::StellarEvolutionInit(double dTime) {
    char achPath[280];
    ServiceStellarEvolutionInit::input in;
    STEV_RAWDATA *CCSNData, *AGBData, *SNIaData, *LifetimeData;
    auto achStelEvolPath = parameters.get_achStelEvolPath();

    /* Read the tables */
    snprintf(achPath, sizeof(achPath), "%s/CCSN.hdf5", achStelEvolPath.data());
    CCSNData = stev::ReadTable(achPath);
    assert(CCSNData->nZs == STEV_CCSN_N_METALLICITY);
    assert(CCSNData->nMasses == STEV_CCSN_N_MASS);
    assert(CCSNData->nElems == ELEMENT_COUNT);

    snprintf(achPath, sizeof(achPath), "%s/AGB.hdf5", achStelEvolPath.data());
    AGBData = stev::ReadTable(achPath);
    assert(AGBData->nZs == STEV_AGB_N_METALLICITY);
    assert(AGBData->nMasses == STEV_AGB_N_MASS);
    assert(AGBData->nElems == ELEMENT_COUNT);

    snprintf(achPath, sizeof(achPath), "%s/SNIa.hdf5", achStelEvolPath.data());
    SNIaData = stev::ReadSNIaTable(achPath);
    assert(SNIaData->nElems == ELEMENT_COUNT);

    snprintf(achPath, sizeof(achPath), "%s/Lifetimes.hdf5", achStelEvolPath.data());
    LifetimeData = stev::ReadLifetimeTable(achPath);
    assert(LifetimeData->nZs == STEV_LIFETIME_N_METALLICITY);
    assert(LifetimeData->nMasses == STEV_LIFETIME_N_MASS);

    /* NOTE: The lowest value of the initial mass array is set to the corresponding value
       of the Lifetime table that is being used, while its highest value to the maximum
       stellar mass for a CCSN to occur, as set in dCCSNMaxMass. The IMF is still
       normalized in the range [dIMFMinMass,dIMFMaxMass] */
    const double dMinMass = LifetimeData->pfInitialMass[0];
    const double dMaxMass = parameters.get_dCCSNMaxMass();
    double adInitialMass[STEV_INTERP_N_MASS], adIMF[STEV_INTERP_N_MASS];

    auto IMF = ChooseIMF(parameters.get_achIMFType().data(), parameters.get_dIMFMinMass(), parameters.get_dIMFMaxMass());
    IMF->MassWeightedSample(dMinMass, dMaxMass, STEV_INTERP_N_MASS, adInitialMass, adIMF);

    /* Store the data */
    for (auto i = 0; i < STEV_INTERP_N_MASS; ++i) {
        in.StelEvolData.afInitialMass[i] = adInitialMass[i];
        in.StelEvolData.afIMFLogWeight[i] = adIMF[i];
    }
    const double dDeltaLog = (log10(dMaxMass) - log10(dMinMass)) / (STEV_INTERP_N_MASS - 1);
    in.StelEvolData.fDeltaLogMass = dDeltaLog;

    /* Convert CCSN/AGB initial mass arrays to log for interpolation */
    for (auto i = 0; i < CCSNData->nMasses; ++i) CCSNData->pfInitialMass[i] = log10(CCSNData->pfInitialMass[i]);
    for (auto i = 0; i < AGBData->nMasses; ++i) AGBData->pfInitialMass[i] = log10(AGBData->pfInitialMass[i]);

    /* Interpolate yields and ejected masses to the initial mass array */
    stev::InterpToIMFSampling(&in.StelEvolData, CCSNData, AGBData, parameters.get_dCCSNMinMass());

    for (auto i = 0; i < STEV_CCSN_N_METALLICITY; ++i) {
        if (CCSNData->pfMetallicity[i] > 0.0f) in.StelEvolData.afCCSNMetallicity[i] = log10(CCSNData->pfMetallicity[i]);
        else in.StelEvolData.afCCSNMetallicity[i] = STEV_MIN_LOG_METALLICITY;
    }

    for (auto i = 0; i < STEV_AGB_N_METALLICITY; ++i) {
        if (AGBData->pfMetallicity[i] > 0.0f) in.StelEvolData.afAGBMetallicity[i] = log10(AGBData->pfMetallicity[i]);
        else in.StelEvolData.afAGBMetallicity[i] = STEV_MIN_LOG_METALLICITY;
    }

    for (auto i = 0; i < ELEMENT_COUNT; ++i) in.StelEvolData.afSNIaEjectedMass[i] = SNIaData->pfEjectedMass[i];
    in.StelEvolData.fSNIaEjectedMetalMass = *SNIaData->pfMetalYield;

    for (auto i = 0; i < STEV_LIFETIME_N_METALLICITY; ++i) {
        if (LifetimeData->pfMetallicity[i] > 0.0f) in.StelEvolData.afLifetimeMetallicity[i] = log10(LifetimeData->pfMetallicity[i]);
        else in.StelEvolData.afLifetimeMetallicity[i] = STEV_MIN_LOG_METALLICITY;
    }
    for (auto i = 0; i < STEV_LIFETIME_N_MASS; ++i) in.StelEvolData.afLifetimeInitialMass[i] = log10(LifetimeData->pfInitialMass[i]);
    for (auto i = 0; i < STEV_LIFETIME_N_METALLICITY * STEV_LIFETIME_N_MASS; ++i) {
        in.StelEvolData.afLifetime[i] = log10(LifetimeData->pfLifetime[i] * SECONDSPERYEAR / units.dSecUnit);
    }

    strcpy(in.achSNIaDTDType, parameters.get_achSNIaDTDType().data());
    in.dTime = dTime;
    in.dSNIaMaxMass = parameters.get_dSNIaMaxMass();
    in.dCCSNMinMass = parameters.get_dCCSNMinMass();
    in.dCCSNMaxMass = parameters.get_dCCSNMaxMass();

    /* Send data and finish */
    mdl->RunService<ServiceStellarEvolutionInit>(in);

    stev::FreeTable(CCSNData);
    stev::FreeTable(AGBData);
    stev::FreeSNIaTable(SNIaData);
    stev::FreeLifetimeTable(LifetimeData);
}
#endif

void MSR::StarFormInit(double dTime) {
    ServiceStarFormInit::input in{
        .dTime = dTime,
        .dCCSNFBDelay = calc.dCCSNFBDelay,
        .dSNIaFBDelay = calc.dSNIaFBDelay,
        .bCCSNFeedback = parameters.get_bCCSNFeedback(),
        .bSNIaFeedback = parameters.get_bSNIaFeedback(),
        .bFeedback = DoStellarFeedback(),
    };
    ServiceStarFormInit::output out;
    mdl->RunService<ServiceStarFormInit>(in, out);

    printf("%lu star particles are about to explode in this IC\n", out);
}

void MSR::StarForm(double dTime, double dDelta, int iRung) {
    struct outStarForm out;

    TimerStart(TIMER_STARFORM);

    const double a = csmTime2Exp(csm, dTime);

    ServiceStarForm::input in{
        .dHubble = Comove() ? csmTime2Hub(csm, dTime) : 0.0,
        .dDelta = dDelta,
        .dTime = dTime,
        .dScaleFactor = a,
        .dSFindexKS = parameters.get_dSFindexKS(),
        .dSFnormalizationKS = calc.dSFnormalizationKS,
        .dConstGamma = parameters.get_dConstGamma(),
        .dSFGasFraction = parameters.get_dSFGasFraction(),
        .dSFThresholdDen = calc.dSFThresholdDen * a * a * a, // Send in comoving units
        .dSFThresholdOD = calc.dSFThresholdOD,
        .dSFThresholdu = calc.dSFThresholdu,
        .dSFEfficiency = parameters.get_dSFEfficiency(),
        .dSNFBEfficiency = parameters.get_dSNFBEfficiency(),
        .dSNFBMaxEff = parameters.get_dSNFBMaxEff(),
        .dSNFBEffIndex = parameters.get_dSNFBEffIndex(),
        .dSNFBEffnH0 = calc.dSNFBEffnH0,
#ifdef STELLAR_EVOLUTION
        .dSNIaMaxMass = parameters.get_dSNIaMaxMass(),
        .dCCSNMinMass = parameters.get_dCCSNMinMass(),
        .dCCSNMaxMass = parameters.get_dCCSNMaxMass(),
#endif
#ifdef HAVE_METALLICITY
        .bSFThresholdDenSchaye2004 = parameters.get_bSFThresholdDenSchaye2004(),
#endif
        .bCCSNFeedback = parameters.get_bCCSNFeedback(),
        .bSNIaFeedback = parameters.get_bSNIaFeedback(),
        .eEOS = get_eEOS(),
    };
    if (parameters.get_bVDetails()) printf("Star Form (rung: %d) ... ", iRung);

    ActiveRung(iRung, 1);
    mdl->RunService<ServiceStarForm>(in, out);
    if (parameters.get_bVDetails()) printf("%d Stars formed with mass %g, %d gas deleted\n", out.nFormed, out.dMassFormed, out.nDeleted);
    massFormed += out.dMassFormed;
    starFormed += out.nFormed;

    nGas -= out.nFormed;
    nStar += out.nFormed;

    TimerStop(TIMER_STARFORM);
    auto dsec = TimerGet(TIMER_STARFORM);
    printf("Star Formation Calculated, Wallclock: %f secs\n\n", dsec);
}

void MSR::SetStarFormationParam() {
    calc.dSFThresholdu = parameters.get_dSFThresholdT() * dTuFacPrimNeutral;

    if (csm->val.bComove) {
        // As usual, in code units the critical density is 1.0
        calc.dSFThresholdOD = parameters.get_dSFMinOverDensity() * csm->val.dOmegab;
    }
    else {
        calc.dSFThresholdOD = 0.0;
    }

    const double dnHToRho = MHYDR / parameters.get_dInitialH() / units.dGmPerCcUnit;
    calc.dSFThresholdDen = parameters.get_dSFThresholdDen() * dnHToRho; // Code physical density now
    const double Msolpcm2 = 1. / units.dMsolUnit * pow(units.dKpcUnit * 1e3, 2);
    calc.dSFnormalizationKS = parameters.get_dSFnormalizationKS() / units.dMsolUnit * units.dSecUnit / SECONDSPERYEAR * pow(units.dKpcUnit, 2) *
                              pow(Msolpcm2, -parameters.get_dSFindexKS());
}

int MSR::ValidateStarFormationParam() {
    if (!parameters.get_bDoEOS()) {
        print_warning("WARNING: Star formation is active but no eEOS is selected!\n");
    }
    return 1;
}

void MSR::SetCoolingParam() {
    calc.dCoolingMinu = parameters.get_dCoolingMinTemp() * dTuFacPrimNeutral;
}

/**
 * Initialises properties stored in the cooling_function_data struct
 */
void MSR::CoolingInit(float redshift) {
    printf("Initializing cooling \n");

    /* Allocate the needed structs */
    cooling = (struct cooling_function_data*)malloc(sizeof(struct cooling_function_data));

    /* Read model parameters */

    strcpy(cooling->cooling_table_path, parameters.get_achCoolingTables().data());

    /* Despite the names, the values of H_reion_heat_cgs and He_reion_heat_cgs
     * that are read in are actually in units of electron volts per proton mass.
     * We later convert to units just below */

    cooling->H_reion_z = parameters.get_fH_reion_z();
    if (redshift <= cooling->H_reion_z) {
        cooling->H_reion_done = 1;
    }
    else {
        cooling->H_reion_done = 0;
    }
    cooling->H_reion_heat_cgs = parameters.get_fH_reion_eV_p_H();
    cooling->He_reion_z_centre = parameters.get_fHe_reion_z_centre();
    cooling->He_reion_z_sigma = parameters.get_fHe_reion_z_sigma();
    cooling->He_reion_heat_cgs = parameters.get_fHe_reion_eV_p_H();

    /* Optional parameters to correct the abundances */
    cooling->Ca_over_Si_ratio_in_solar = calc.fCa_over_Si_in_Solar;
    cooling->S_over_Si_ratio_in_solar = calc.fS_over_Si_in_Solar;

    /* Convert H_reion_heat_cgs and He_reion_heat_cgs to cgs
     * (units used internally by the cooling routines). This is done by
     * multiplying by 'eV/m_H' in internal units, then converting to cgs units.
     * Note that the dimensions of these quantities are energy/mass = velocity^2
     */
#define EV_CGS 1.602e-12

    cooling->H_reion_heat_cgs *= (EV_CGS) / (MHYDR);

    cooling->He_reion_heat_cgs *= (EV_CGS) / (MHYDR);

    /* Read in the list of redshifts */
    get_cooling_redshifts(cooling);

    /* Read in cooling table header */
    char fname[eagle_table_path_name_length + 13];
    snprintf(fname, sizeof(fname), "%s/z_0.000.hdf5", cooling->cooling_table_path);
    read_cooling_header(fname, cooling);

    /* Allocate space for cooling tables */
    allocate_cooling_tables(cooling);

    /* Compute conversion factors */
    // This is not ideal, and does not follow PKDGRAV3 philosophy... someday
    // it should be reworked
    cooling->internal_energy_to_cgs = units.dErgPerGmUnit;
    cooling->internal_energy_from_cgs = 1. / cooling->internal_energy_to_cgs;
    cooling->number_density_to_cgs = pow(units.dKpcUnit * KPCCM, -3.);
    cooling->units = units;
    cooling->dConstGamma = parameters.get_dConstGamma();
    cooling->dCoolingMinu = calc.dCoolingMinu;

    /* Store some constants in CGS units */
    const double proton_mass_cgs = MHYDR;
    cooling->inv_proton_mass_cgs = 1. / proton_mass_cgs;
    cooling->T_CMB_0 = parameters.get_fT_CMB_0();

    const double compton_coefficient_cgs = 1.0178085e-37;

    /* And now the Compton rate */
    cooling->compton_rate_cgs = compton_coefficient_cgs * cooling->T_CMB_0 * cooling->T_CMB_0 * cooling->T_CMB_0 * cooling->T_CMB_0;

    /* Set the redshift indices to invalid values */
    cooling->z_index = -10;

    /* set previous_z_index and to last value of redshift table*/
    cooling->previous_z_index = eagle_cooling_N_redshifts - 2;

    // We prepare the data to be scattered among the processes
    ServiceCoolingInit::input in;
    in.in_cooling_data = *(cooling);
    for (int i = 0; i < eagle_cooling_N_redshifts; i++) in.Redshifts[i] = cooling->Redshifts[i];
    for (int i = 0; i < eagle_cooling_N_density; i++) in.nH[i] = cooling->nH[i];
    for (int i = 0; i < eagle_cooling_N_He_frac; i++) in.HeFrac[i] = cooling->HeFrac[i];
    for (int i = 0; i < eagle_cooling_N_temperature; i++) in.Temp[i] = cooling->Temp[i];
    for (int i = 0; i < eagle_cooling_N_temperature; i++) in.Therm[i] = cooling->Therm[i];
    for (int i = 0; i < eagle_cooling_N_abundances; i++) in.SolarAbundances[i] = cooling->SolarAbundances[i];
    for (int i = 0; i < eagle_cooling_N_abundances; i++) in.SolarAbundances_inv[i] = cooling->SolarAbundances_inv[i];

    mdl->RunService<ServiceCoolingInit>(in);
}

/**
 * @brief Common operations performed on the cooling function at a
 * given time-step or redshift. Predominantly used to read cooling tables
 * above and below the current redshift, if not already read in.
 *
 * Also calls the additional H reionization energy injection if need be.
 *
 * @param cosmo The current cosmological model.
 * @param cooling The #cooling_function_data used in the run.
 * @param s The space data, including a pointer to array of particles
 *
 *
 */
void MSR::CoolingUpdate(float redshift) {
    printf("Updating cooling.. z=%f\n", redshift);

    /* Current redshift */
    // const float redshift = cosmo->z;

    /* What is the current table index along the redshift axis? */
    int z_index = -1;
    float dz = 0.f;
    cooling::get_redshift_index(redshift, &z_index, &dz, cooling);
    cooling->dz = dz;

    // We update dz in each process
    mdl->RunService<ServiceCoolingUpdateZ>(dz);

    /* Do we already have the correct tables loaded? */
    if (cooling->z_index == z_index) return;

    /* Which table should we load ? */
    if (z_index == eagle_cooling_N_redshifts) {
        /* Before reionization. Just load the corresponding table */
        get_redshift_invariant_table(cooling, /* photodis=*/1);
    }
    else {
        /* After reionization. Just load the corresponding table */
        if (z_index > eagle_cooling_N_redshifts) {
            /* Between reionization and first z-dependent table */
            get_redshift_invariant_table(cooling, /* photodis=*/0);
        }
        else {
            /* Normal case: two tables bracketing the current z */
            get_cooling_table(cooling, z_index, z_index + 1);
        }

        /* Extra energy for reionization */
        if (!cooling->H_reion_done) {
            printf("Hydrogen reionization! Injecting %.1f eV per proton\n", parameters.get_fH_reion_eV_p_H());
            mdl->RunService<ServiceCoolingHydReionization>();

            /* Flag that reionization happened */
            cooling->H_reion_done = 1;
        }
    }

    /* Store the currently loaded index */
    cooling->z_index = z_index;
    // We send the newly loaded table
    printf("Sending cooling tables...\n");
    ServiceCoolingUpdate::input in;
    in.z_index = cooling->z_index;
    in.previous_z_index = cooling->previous_z_index;
    in.dz = cooling->dz;

    for (int i = 0; i < eagle_cooling_N_loaded_redshifts * num_elements_metal_heating; i++) in.metal_heating[i] = cooling->table.metal_heating[i];

    for (int i = 0; i < eagle_cooling_N_loaded_redshifts * num_elements_HpHe_heating; i++) in.H_plus_He_heating[i] = cooling->table.H_plus_He_heating[i];

    for (int i = 0; i < eagle_cooling_N_loaded_redshifts * num_elements_HpHe_electron_abundance; i++)
        in.H_plus_He_electron_abundance[i] = cooling->table.H_plus_He_electron_abundance[i];

    for (int i = 0; i < eagle_cooling_N_loaded_redshifts * num_elements_temperature; i++) in.temperature[i] = cooling->table.temperature[i];

    for (int i = 0; i < eagle_cooling_N_loaded_redshifts * num_elements_electron_abundance; i++)
        in.electron_abundance[i] = cooling->table.electron_abundance[i];

    mdl->RunService<ServiceCoolingUpdate>(in);
}
