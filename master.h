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

#ifndef MASTER_HINCLUDED
#define MASTER_HINCLUDED
#include "pkd_config.h"
#include <stdint.h>
#include <signal.h>
#include <time.h>
#include <vector>
#include <string_view>
#include <Python.h>
#include "fmt/format.h" // This will be part of c++20
#include "fmt/ostream.h"
#include <fmt/color.h>

#include "pst.h"
#include "mdl.h"
#include "core/memory.h"
#include "pkd_parameters.h"
#include "pkd_enumerations.h"
#include "cooling/cooling_struct.h"
#include "gravity/gravity.h"
#define MSR_INIT_E 1
#define MSR_STEP_E 0

extern time_t timeGlobalSignalTime;
extern int bGlobalOutput;

// IMPORTANT: If you change the timers here then you need to change
// their names in master.cxx (timer_names)
enum msrTimers {
    TIMER_GRAVITY = 0,
    TIMER_IO,
    TIMER_TREE,
    TIMER_DOMAIN,
    TIMER_KICKO,
    TIMER_KICKC,
    TIMER_DENSITY,
    TIMER_ENDINT,
    TIMER_GRADIENTS,
    TIMER_FLUXES,
    TIMER_TIMESTEP,
    TIMER_DRIFT,
    TIMER_FOF,
    TIMER_FEEDBACK,
    TIMER_STARFORM,
    TIMER_BHS,
#ifdef STELLAR_EVOLUTION
    TIMER_STEV,
#endif
    TIMER_NONE,
    TOTAL_TIMERS
};

struct MSRINSTANCE {
    PyObject_HEAD class MSR* msr;
};

class MSR {
  protected:
    const PST pst;
    mdl::mdlClass* mdl;
    PyObject* parameter_overrides = nullptr;

  public:
    explicit MSR(MDL mdl, PST pst);
    ~MSR();

  public:
    template <typename... Args>
    void print(Args&&... args) const {
        FMT_TRY {
            fmt::print(std::forward<Args>(args)...);
        }
        FMT_CATCH(fmt::format_error & f) {
            fmt::print(stderr, fg(fmt::color::crimson) | fmt::emphasis::bold, "ERROR in formatting string: {}\n", f.what());
            mdl->show_backtrace_with_snippets(5, 2);
        }
    }
    template <typename... Args>
    void print_error(const char* format, Args&&... args) const {
        print(stderr, fg(fmt::color::crimson) | fmt::emphasis::bold, format, std::forward<Args>(args)...);
    }
    template <typename... Args>
    void print_warning(const char* format, Args&&... args) const {
        if (parameters.get_bVWarnings()) {
            print(stderr, fg(fmt::color::yellow) | fmt::emphasis::bold, format, std::forward<Args>(args)...);
        }
    }
    template <typename... Args>
    void print_notice(const char* format, Args&&... args) const {
        print(fg(fmt::color::green) | fmt::emphasis::bold, format, std::forward<Args>(args)...);
    }
    template <typename... Args>
    void print_detail(const char* format, Args&&... args) const {
        if (parameters.get_bVDetails()) {
            print(format, std::forward<Args>(args)...);
        }
    }
    int Python(int argc, char* argv[]);
    int ValidateParameters();
    void SetDerivedParameters(bool bRestart = false);
    void Hostname();
    void MemStatus();
    int GetLock();
    void IgnoreSIGBUS();
    double LoadOrGenerateIC();
    void Simulate(double dTime, int iStartStep, bool bRestart = false);
    void Simulate(double dTime);

  private:
    void persist(PyObject* file, PyObject* obj);
    void persist(PyObject* file, int n);
    void persist(PyObject* file, double d);
    template <typename T>
    T restore(PyObject* file, PyObject* replace = Py_None);

  protected:
    void stat_files(std::vector<uint64_t>& counts, const std::string_view& filename_template, uint64_t element_size);
    void Restore(const std::string& filename, int nSizeParticle, int nSizeDecorator);
    void PrintStat(const STAT& ps, char const* pszPrefix, int p);
#ifdef INSTRUMENT_MPI
    void PrintStat(const WelfordStats<double>& ps, char const* pszPrefix, int p);
#endif

  public:
    int64_t parallel_count(int64_t nParallel, bool bParallel = true);
    int64_t parallel_read_count();
    int64_t parallel_write_count();

  public:
    // I/O and IC Generation
    double GenerateIC(int nGrid, int iSeed, double z, double L, CSM csm = nullptr);
    void Restart(const char* filename, PyObject* kwargs, PyObject* species = Py_None, PyObject* classes = Py_None, PyObject* step = Py_None,
                 PyObject* steps = Py_None, PyObject* time = Py_None, PyObject* delta = Py_None, PyObject* E = Py_None, PyObject* U = Py_None,
                 PyObject* Utime = Py_None);
    std::tuple<double, double, int64_t, int64_t, int64_t, int64_t> // dTime,dDelta,iStep,nSteps,nSizeParticle,nSizeNode
    ReadCheckpoint(const char* filename, PyObject* kwargs, PyObject* species = Py_None, PyObject* classes = Py_None, PyObject* step = Py_None,
                   PyObject* steps = Py_None, PyObject* time = Py_None, PyObject* delta = Py_None, PyObject* E = Py_None, PyObject* U = Py_None,
                   PyObject* Utime = Py_None);
    double Read(std::string_view achInFile);
    void Checkpoint(int iStep, int nSteps, double dTime, double dDelta);
    void Write(std::string_view FileName, double dTime, int bCheckpoint);
    void OutASCII(std::string_view File, int iType, int nDims, int iFileType);
    void OutArray(std::string_view File, int iType, int iFileType);
    void OutArray(std::string_view File, int iType);
    void OutVector(std::string_view File, int iType, int iFileType);
    void OutVector(std::string_view File, int iType);
    void Output(int iStep, double dTime, double dDelta, int bCheckpoint);
    void OutputFineStatistics(double dStep, double dTime);

    void RecvArray(void* vBuffer, PKD_FIELD field, int iUnitSize, double dTime, bool bMarked = false);

    // Particle order, domains, trees
    void Reorder();
    void DomainDecomp(int iRung = 0);
    void BuildTree(int bNeedEwald);
    void BuildTreeFixed(int bNeedEwald, uint8_t uRungDD);
    void BuildTreeActive(int bNeedEwald, uint8_t uRungDD);
    void BuildTreeMarked(int bNeedEwald = 0);
    void ReorderDecorator();
#ifdef OPTIM_REORDER_IN_NODES
    void ReorderInNodes();
#endif

    // Gravity
    uint8_t Gravity(uint8_t uRungLo, uint8_t uRungHi, int iRoot1, int iRoot2, double dTime, double dDelta, double dStep, double dTheta, int bKickClose,
                    int bKickOpen, int bEwald, int bGravStep, int nPartRhoLoc, int iTimeStepCrit);
    uint8_t Gravity(uint8_t uRungLo, uint8_t uRungHi, int iRoot1, int iRoot2, double dTime, double dDelta, double dStep, double dTheta, int bKickClose,
                    int bKickOpen, int bEwald, int bGravStep, int nPartRhoLoc, int iTimeStepCrit, SPHOptions SPHoptions);

    // Analysis
    void Smooth(double dTime, double dDelta, int iSmoothType, int bSymmetric, int nSmooth);
    uint64_t ReSmooth(double dTime, double dDelta, int iSmoothType, int bSymmetric);
    void NewFof(double dTau, int nMinMembers);
    void Hop(double dTime, double dDelta);
    void GroupStats();
    void WriteGroupStats(const std::string_view& fname, int nParallel);
    void WriteGroupStats(const std::string_view& fname);
    EphemeralMemory EphemeralMemoryGrid(int nGrid, int nCount);
    std::tuple<std::vector<uint64_t>, std::vector<float>, std::vector<float>, std::vector<float>> // nPk, fK, fPk, fPkAll
    MeasurePk(int iAssignment, int bInterlace, int nGrid, double a, int nBins, int fold = 1);
    void AssignMass(int iAssignment = 3, int iGrid = 0, float fDelta = 0.0f, int fold = 1);
    void DensityContrast(int nGrid, bool k = true);
    void WindowCorrection(int iAssignment, int iGrid);
    void Interlace(int iGridTarget, int iGridSource);
    void AddLinearSignal(int iGrid, int iSeed, double Lbox, double a, bool bFixed = false, float fPhase = 0);
    std::tuple<std::vector<uint64_t>, std::vector<float>, std::vector<float>> // nPk, fK, fPk
    GridBinK(int nBins, int iGrid);
    void BispectrumSelect(int iGridTarget, int iGridSource, double kmin, double kmax);
    auto BispectrumCalculate(int iGrid1, int iGrid2, int iGrid3) -> double;
    void GridCreateFFT(int nGrid);
    void GridDeleteFFT();

  private:
    typedef struct {
        double dFrac;     /* Fraction of particles in each bin */
        uint64_t nTotal;  /* Total number of particles in the range */
        uint64_t nInner;  /* Number inside minimum radius */
        uint64_t nTarget; /* Target number of particles */
        uint64_t nSelected;
        MSR* msr;
    } SPHERECTX;
    static double countSphere(double r, void* vctx);
    static void profileRootFind(double* dBins, int lo, int hi, int nAccuracy, SPHERECTX* ctx);

    typedef struct {
        double rMiddle;
        total_t nTarget; /* Target number of particles */
        MSR* msr;
    } SHELLCTX;
    static double countShell(double rInner, void* vctx);

  public:
    struct msr_analysis_callback {
        PyObject* callback;
        EphemeralMemory memory;
        explicit msr_analysis_callback(PyObject* callback, uint64_t per_particle, uint64_t per_process)
            : callback(callback), memory(per_particle, per_process) {
            Py_INCREF(callback);
        }
        ~msr_analysis_callback() { Py_DECREF(callback); }
    };
    void addAnalysis(PyObject* callback, uint64_t per_particle, uint64_t per_process);
    void runAnalysis(int iStep, double dTime);

  protected:
    std::list<msr_analysis_callback> analysis_callbacks;

  protected:
    UNITS units;

  public: // should be private
    pkd_parameters parameters;
    double set_dynamic(int iStep, double dTime) {
        parameters.set_dynamic("step", iStep);
        parameters.set_dynamic("time", dTime);
        if (csm->val.bComove) parameters.set_dynamic("a", csmTime2Exp(csm, dTime));
        auto dTheta = parameters.get_dTheta();
        parameters.set_dynamic("theta", dTheta);
        return dTheta;
    }

  private:
    PyObject* pDill = nullptr;             // The "dill" module
    PyObject* pDill_load = nullptr;        // The "dill.load" function
    PyObject* pDill_dump = nullptr;        // The "dill.dump" function
    PyObject* pDill_load_module = nullptr; // The "dill.load_module" function
    PyObject* pDill_dump_module = nullptr; // The "dill.dump_module" function
    PyObject* particle_dtype = nullptr;    // The numpy dtype for particles, built from the FieldBuilder

  public:
    struct CALC calc;

  public:
    blitz::TinyVector<double, 3> fCenter;

  public:
    /*
    ** Parameters.
    */
    CSM csm;
    /*
    ** Other stuff...
    */
    int nThreads;
    uint64_t N;
    uint64_t nDark;
    uint64_t nGas;
    uint64_t nStar;
    uint64_t nBH;
    uint64_t nMaxOrder; /* Order number of last particle */
    int nClasses;
    int iCurrMaxRung;

    /*
    ** Timers
    */
    struct msrtimer {
        double sec;
        double acc;
    } ti[TOTAL_TIMERS];

    struct cooling_function_data* cooling;
    struct cooling_tables* cooling_table;
    int starFormed;
    double massFormed;
    /*
     * File for the fine-grained output
     */
    FILE* fpFineLog;
    /*
    ** Tree moments (for Ewald)
    */
    MOMC momTreeRoot[IRUNGMAX + 1];
    double momTreeCom[IRUNGMAX + 1][3];

    /*
    ** Comoving coordinate variables.
    */
    double dEcosmo;
    double dUOld;
    double dTimeOld;
    /*
    ** Delta to use for each set of steps
    */
    std::vector<double> dDelta_list;
    std::vector<uint64_t> iStep_list;
    // int nMaxOuts;
    // int nOuts;
    // double *pdOutTime;

    /*
     * Domain Decomposition Done
     */
    std::vector<uint64_t> nRung;
    int iRungDD, iRungDT;
    int iLastRungRT, iLastRungDD;
    uint64_t nActive;
    int nGroups;
    int nBins;
    int bAntiGrav;

    long lStart; /* starting time of job */
    long lPrior; /* starting time of last step */

    /* Gas */
    double dTuFac;
    double dTuFacPrimIonised;
    double dTuFacPrimNeutral;
    int bUpdateBall;

    /* Values for a restore from checkpoint */
    double dCheckpointTime;
    int iCheckpointStep;
    int nCheckpointThreads;
    char achCheckpointName[PATH_MAX];

  protected:
    static double Time();
    static void Leader();
    static void Trailer();
    std::string MakePath(std::string_view base, std::string_view ext = {});

    double getTime(double dExpansion); // Return simulation time
    double getVfactor(double dTime);
    bool getDeltaSteps(double dTime, int iStartStep);

    auto NoSteps() const {
        auto nSteps = parameters.get_nSteps();
        return nSteps.size() == 1 && nSteps[0] == 0;
    }

    auto OutName() const { return parameters.get_achOutName(); }
    int LogInterval() const { return parameters.get_iLogInterval(); }
    int OutInterval() const { return parameters.get_iOutInterval(); }
    int CheckInterval() const { return parameters.get_iCheckInterval(); }
    double Soft() const { return parameters.get_dSoft(); }
    [[nodiscard]] auto DoDensity() const -> bool { return parameters.get_bDoDensity(); }
    [[nodiscard]] auto DoGas() const -> bool { return parameters.get_hydro_model() != HYDRO_MODEL::NONE; }
    [[nodiscard]] auto NewSPH() const -> bool { return parameters.get_hydro_model() == HYDRO_MODEL::SPH; }
    [[nodiscard]] auto MeshlessFiniteMass() const -> bool {
        auto model = parameters.get_hydro_model();
        return model == HYDRO_MODEL::MFM;
    }
    [[nodiscard]] auto MeshlessFiniteVolume() const -> bool {
        auto model = parameters.get_hydro_model();
        return model == HYDRO_MODEL::MFV;
    }
    [[nodiscard]] auto MeshlessHydro() const -> bool { return MeshlessFiniteMass() || MeshlessFiniteVolume(); }
    [[nodiscard]] auto DoStellarFeedback() const -> bool { return parameters.get_bCCSNFeedback() || parameters.get_bSNIaFeedback(); }
    void StellarFeedback(double dTime);
    void ChemEnrich(double dTime);

    [[nodiscard]] auto DoGravity() const -> bool { return parameters.get_bDoGravity(); }
    double Eta() const { return parameters.get_dEta(); }
    int MaxRung() const { return parameters.get_iMaxRung(); }
    int Comove() const { return csm->val.bComove; }
    uint64_t MaxOrder() const { return nMaxOrder; }
    int CurrMaxRung() const { return iCurrMaxRung; }

    std::string BuildName(std::string_view path, int iStep, std::string_view type = {});
    std::string BuildNameWithI(std::string_view path, int iStep, std::string_view type = {});
    std::string BuildName(int iStep, std::string_view type = {});       // With achOutPath
    std::string BuildNameWithI(int iStep, std::string_view type = {});  // With achOutPath: ensure {i}
    std::string BuildNameWithOI(int iStep, std::string_view type = {}); // With achOutPath: ensure {octant} and {i}
    std::string BuildIoName(int iStep, std::string_view type = {});     // With achIoPath
    std::string BuildCpName(int iStep);                                 // With achCpPath

    int ReadOuts(double dTime);
    void Exit(int status);
    uint64_t getMemoryModel();
    std::tuple<int, int, int> InitializePStore(uint64_t* nSpecies, uint64_t mMemoryModel, uint64_t nEphemeral);
    int CheckForStop(const char* achStopFile);
    int CheckForOutput(int iStep, int nSteps, double dTime, int* pbDoCheckpoint, int* pbDoOutput);
    auto GetClasses() -> std::vector<PARTCLASS>;
    void SetClasses();
    void AllNodeWrite(std::string_view FileName, double dTime, double dvFac, int bDouble);
    auto SetTotal() -> uint64_t;
    auto CalcWriteStart() -> uint64_t;
    void SwitchTheta(double);
    double getTheta(double dTime);
    double SwitchDelta(double dTime, int iStep);
    void InitCosmology(CSM csm);
    void BuildTree(int bNeedEwald, uint32_t uRoot, uint32_t utRoot);
    void ActiveRung(int iRung, int bGreater);
    auto ActiveOrder() -> uint64_t;
    void CalcBound(Bound& bnd);
    void CalcBound();
    void CountSpecies();
    void RemoveDeleted();
    double AdjustTime(double aOld, double aNew);
    void UpdateSoft(double dTime);
    void OutputOrbits(int iStep, double dTime);
    double TotalMass();
    void LightConeOpen(int iStep);
    void LightConeClose(int iStep);
    void LightConeVel();
#ifdef MDL_FFTW
    void SetLinearGrid(double dTime, double dDelta, int nGrid, int bKickClose, int bKickOpen);
    void LinearKick(double dTime, double dDelta, int bKickClose, int bKickOpen);
#endif

    // Timers
    void TimerStart(int iTimer);
    void TimerStop(int iTimer);
    double TimerGet(int iTimer);
    double TimerGetAcc(int iTimer);
    void TimerHeader();
    void TimerRestart();
    void TimerDump(int iStep);
    void CalcEandL(int bFirst, double dTime, double& E, double& T, double& U, double& Eth, blitz::TinyVector<double, 3>& L, blitz::TinyVector<double, 3>& F,
                   double& W);
    void Drift(double dTime, double dDelta, int iRoot);

    auto get_eEOS() -> eEOSparam;
    void ZeroNewRung(uint8_t uRungLo, uint8_t uRungHi, int uRung);
    void KickKDKOpen(double dTime, double dDelta, uint8_t uRungLo, uint8_t uRungHi);
    void KickKDKClose(double dTime, double dDelta, uint8_t uRungLo, uint8_t uRungHi);
    void UpdateRung(uint8_t uRung);
    void AccelStep(uint8_t uRungLo, uint8_t uRungHi, double dTime, double dDelta);
    auto GetMaxRung() -> uint8_t;
    void SetGlobalRung(uint8_t uRung);
    void DensityStep(uint8_t uRungLo, uint8_t uRungHi, double dTime, double dDelta);

    void FastGasPhase1(double dTime, double dDelta, int iSmoothType);
    void FastGasPhase2(double dTime, double dDelta, int iSmoothType);
    void CoolSetup(double dTime);
    void Cooling(double dTime, double dStep, int bUpdateState, int bUpdateTable, int bInterateDt);
    void InitSph(double dTime, double dDelta, bool bRestart);
    auto CountDistance(double dRadius2Outer, double dRadius2Inner = 0.0) -> uint64_t;

    // Meshless hydrodynamics
    void MeshlessGradients(double dTime, double dDelta);
    void MeshlessFluxes(double dTime, double dDelta);
    void HydroStep(double dTime, double dDelta);
    void ComputeSmoothing(double dTime, double dDelta);
    void ChemCompInit();
    void EndTimestepIntegration(double dTime, double dDelta);

    // Cooling
    void SetCoolingParam();
    void CoolingUpdate(float redshift);
    void CoolingInit(float redshift);
#ifdef GRACKLE
    void GrackleInit(int bComove, double dScaleFactor);
#endif
    void SetStarFormationParam();
    int ValidateStarFormationParam();
    void StarFormInit(double dTime);
    void StarForm(double dTime, double dDelta, int iRung);
    void SetFeedbackParam();
#ifdef STELLAR_EVOLUTION
    void SetStellarEvolutionParam();
    void StellarEvolutionInit(double dTime);
#endif
    void SetEOSParam();
    int ValidateEOSParam();
    void SetBlackholeParam();
    int ValidateBlackholeParam();
    void BlackholeInit(uint8_t uRungMax);
    void PlaceBHSeed(double dTime, uint8_t uRungMax);
    void BHMerger(double dTime);
    void BHEvolve(double dTime, double dDelta);
    void BHGasPin(double dTime, double dDelta);
    void BHReposition();
    void BHAccretion(double dTime);
    void BHStep(double dTime, double dDelta);

    void Initialize();
    void writeParameters(const std::string& baseName, int iStep, int nSteps, double dTime, double dDelta);
    void DomainDecompOld(int iRung);

    int CountRungs(uint64_t* nRungs);
    void SetSoft(double);
    void InitBall();

    std::tuple<std::vector<uint64_t>, std::vector<float>, std::vector<float>> // nPk, fK, fPk
    MeasureLinearPk(int nGridLin, double a, double dBoxSize);
    void OutputPk(int iStep, double dTime);
    void OutputLinearPk(int iStep, double dTime);

    int NewTopStepKDK(double& dTime,                                           /* MODIFIED: Current simulation time */
                      double dDelta, double dTheta, int nSteps, int bDualTree, /* Should be zero at rung 0! */
                      uint8_t uRung,                                           /* Rung level */
                      double* pdStep,                                          /* Current step */
                      uint8_t* puRungMax, int* pbDoCheckpoint, int* pbDoOutput, int* pbNeedKickOpen);
    void TopStepKDK(double dStep,             /* Current step */
                    double dTime,             /* Current time */
                    double dDelta,            /* Time step */
                    double dTheta, int iRung, /* Rung level */
                    int iKickRung,            /* Gravity on all rungs from iRung
                                              to iKickRung */
                    int iAdjust);             /* Do an adjust? */

    void CalcDistance(blitz::TinyVector<double, 3> dCenter, double dRadius);
    void CalcCOM(blitz::TinyVector<double, 3> dCenter, double dRadius, blitz::TinyVector<double, 3>& com, blitz::TinyVector<double, 3>& vcm,
                 blitz::TinyVector<double, 3>& L, double& M);
    auto initializeSPHOptions(double dTime, double dDelta) -> SPHOptions;
    void SetSPHoptions(double dTime, double dDelta);
    void ResetCOM();
    void InitializeEOS();
    void CalculateKickParameters(struct pkdKickParameters* kick, uint8_t uRungLo, double dTime, double dDelta, double dStep, int bKickClose, int bKickOpen,
                                 SPHOptions SPHoptions);
    void UpdateGasValues(uint8_t uRungLo, double dTime, double dDelta, double dStep, int bKickClose, int bKickOpen, SPHOptions SPHoptions);
    void TreeUpdateFlagBounds(int bNeedEwald, uint32_t uRoot, uint32_t utRoot, SPHOptions SPHoptions);

  public:
    void Profile(const PROFILEBIN** pBins, int* pnBins, blitz::TinyVector<double, 3> r, double dMinRadius, double dLogRadius, double dMaxRadius, int nPerBin,
                 int nBins, int nAccuracy);
    void OutputGrid(std::string_view filename, bool k = false, int iGrid = 0, int nParaWrite = 0);

    uint64_t CountSelected();
    uint64_t SelSpecies(uint64_t mSpecies, int setIfTrue = true, int clearIfFalse = true);
    uint64_t SelAll(int setIfTrue = true, int clearIfFalse = true);
    uint64_t SelGas(int setIfTrue = true, int clearIfFalse = true);
    uint64_t SelStar(int setIfTrue = true, int clearIfFalse = true);
    uint64_t SelDark(int setIfTrue = true, int clearIfFalse = true);
    uint64_t SelDeleted(int setIfTrue = true, int clearIfFalse = true);
    uint64_t SelBlackholes(int setIfTrue = true, int clearIfFalse = true);
    uint64_t SelActives(int setIfTrue = true, int clearIfFalse = true);
    uint64_t SelGroup(int iGroup, int setIfTrue = true, int clearIfFalse = true);
    uint64_t SelMass(double dMinMass, double dMaxMass, int setIfTrue = true, int ClearIfFalse = true);
    uint64_t SelById(uint64_t idStart, uint64_t idEnd, int setIfTrue = true, int clearIfFalse = true);
    uint64_t SelPhaseDensity(double dMinPhaseDensity, double dMaxPhaseDensity, int setIfTrue = true, int clearIfFalse = true);
    uint64_t SelBox(blitz::TinyVector<double, 3> center, blitz::TinyVector<double, 3> size, int setIfTrue = true, int clearIfFalse = true);
    uint64_t SelSphere(blitz::TinyVector<double, 3> r, double dRadius, int setIfTrue = true, int clearIfFalse = true);
    uint64_t SelCylinder(blitz::TinyVector<double, 3> dP1, blitz::TinyVector<double, 3> dP2, double dRadius, int setIfTrue = true, int clearIfFalse = true);

  public:
    void RsLoadIds(std::string_view service_name, std::vector<uint64_t>& counts, const std::string_view filename, bool bAppend = false);
    void RsLoadIds(const std::string_view filename, bool bAppend = false);
    void RsHaloLoadIds(const std::string_view filename, bool bAppend = false);
    void RsSaveIds(const std::string_view filename);
    void RsReorderIds();
    void RsExtract(std::string_view filename_template);
};
#endif
