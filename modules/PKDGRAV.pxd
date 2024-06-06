# distutils: language = c++
# cimport CMSR
# from CMSR cimport PARTCLASS, FIO_SPECIES, PKD_FIELD
# from CMSR cimport OUT_TYPE as CMSR_OUT_TYPE
from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.string cimport string
from libc.stdint cimport int64_t,uint64_t,uint8_t
import cython
import numpy as np
cimport numpy as cnp
cimport cosmo
from cosmology cimport Cosmology
from ephemeral cimport EphemeralMemory
import ephemeral

cdef extern from "blitz/array.h" namespace "blitz" nogil:
    cdef cppclass BLITZ1 "1":
        pass
    cdef cppclass BLITZ2 "2":
        pass
    cdef cppclass BLITZ3 "3":
        pass

    cdef cppclass TinyVector[T,N]:
        TinyVector() except +
        TinyVector(TinyVector[T,N]&) except +
        TinyVector(T) except +
        TinyVector(T,T) except +
        TinyVector(T,T,T) except +


cdef extern from "io/fio.h":
    cpdef enum FIO_SPECIES:
        FIO_SPECIES_DARK
        FIO_SPECIES_SPH
        FIO_SPECIES_STAR
        FIO_SPECIES_BH
        FIO_SPECIES_ALL

cdef extern from "io/outtype.h":
    cpdef enum OUT_TYPE:
        OUT_DENSITY_ARRAY
        OUT_POT_ARRAY
        OUT_AMAG_ARRAY
        OUT_IMASS_ARRAY
        OUT_RUNG_ARRAY
        OUT_DIVV_ARRAY
        OUT_VELDISP2_ARRAY
        OUT_VELDISP_ARRAY
        OUT_PHASEDENS_ARRAY
        OUT_SOFT_ARRAY
        OUT_POS_VECTOR
        OUT_VEL_VECTOR
        OUT_ACCEL_VECTOR
        OUT_MEANVEL_VECTOR
        OUT_IORDER_ARRAY
        OUT_C_ARRAY
        OUT_HSPH_ARRAY
        OUT_RUNGDEST_ARRAY
        OUT_MARKED_ARRAY
        OUT_CACHEFLUX_ARRAY
        OUT_CACHECOLL_ARRAY
        OUT_AVOIDEDFLUXES_ARRAY
        OUT_COMPUTEDFLUXES_ARRAY
        OUT_HOP_STATS
        OUT_GROUP_ARRAY
        OUT_GLOBALGID_ARRAY
        OUT_BALL_ARRAY
        OUT_PSGROUP_ARRAY
        OUT_PSGROUP_STATS

cdef extern from "core/particle.h":
    # ctypedef struct PARTCLASS:
    cdef cppclass PARTCLASS:
        PARTCLASS()
        PARTCLASS(FIO_SPECIES eSpecies,float fMass,float fSoft,int iMat)
        float       fMass
        float       fSoft
        int         iMat
        FIO_SPECIES eSpecies

    cpdef enum PKD_FIELD:
        FIELD_POSITION     "PKD_FIELD::oPosition"
        FIELD_ACCELERATION "PKD_FIELD::oAcceleration"
        FIELD_VELOCITY     "PKD_FIELD::oVelocity"
        FIELD_POTENTIAL    "PKD_FIELD::oPotential"
        FIELD_GROUP        "PKD_FIELD::oGroup"
        FIELD_MASS         "PKD_FIELD::oMass"
        FIELD_SOFTENING    "PKD_FIELD::oSoft"
        FIELD_DENSITY      "PKD_FIELD::oDensity"
        FIELD_BALL         "PKD_FIELD::oBall"
        FIELD_PARTICLE_ID  "PKD_FIELD::oParticleID"
        FIELD_GLOBAL_GID   "PKD_FIELD::oGlobalGid"

cdef extern from "smooth/smoothfcn.h":
    cpdef enum SMOOTH_TYPE:
        SMOOTH_TYPE_DENSITY "SMX_DENSITY"
        SMOOTH_TYPE_F1      "SMX_DENSITY_F1"
        SMOOTH_TYPE_M3      "SMX_DENSITY_M3"
        SMOOTH_TYPE_GRADIENT_M3 "SMX_GRADIENT_M3"
        SMOOTH_TYPE_HOP_LINK "SMX_HOP_LINK"
        SMOOTH_TYPE_BALL    "SMX_BALL"
        SMOOTH_TYPE_PRINTNN "SMX_PRINTNN"
        SMOOTH_TYPE_HYDRO_DENSITY "SMX_HYDRO_DENSITY"
        SMOOTH_TYPE_HYDRO_DENSITY_FINAL "SMX_HYDRO_DENSITY_FINAL"
        SMOOTH_TYPE_HYDRO_GRADIENT "SMX_HYDRO_GRADIENT"
        SMOOTH_TYPE_HYDRO_FLUX "SMX_HYDRO_FLUX"
        SMOOTH_TYPE_HYDRO_STEP "SMX_HYDRO_STEP"
        SMOOTH_TYPE_HYDRO_FLUX_VEC "SMX_HYDRO_FLUX_VEC"
        SMOOTH_TYPE_SN_FEEDBACK "SMX_SN_FEEDBACK"
        SMOOTH_TYPE_BH_MERGER "SMX_BH_MERGER"
        SMOOTH_TYPE_BH_EVOLVE "SMX_BH_EVOLVE"
        SMOOTH_TYPE_BH_STEP "SMX_BH_STEP"
        SMOOTH_TYPE_BH_GASPIN "SMX_BH_GASPIN"
        SMOOTH_TYPE_CHEM_ENRICHMENT "SMX_CHEM_ENRICHMENT"

include "pkd_parameters.pxi"
include "pkd_enumerations.pxi"

cdef extern from "master.h":
    cdef cppclass MSR:
        pkd_parameters parameters

        uint64_t N
        # MSR() except +
        void testv(vector[PARTCLASS] &v)
        #tuple[double,double,int64_t,int64_t,int64_t,int64_t]
        void *ReadCheckpoint(const char *filename, object kwargs,
                 object species, object classes, object step, object steps,
                 object time, object delta, object E, object U, object Utime)
        void Restart(const char *filename,object kwargs)
        void Restart(const char *filename, object kwargs,
                 object species, object classes, object step, object steps,
                 object time, object delta, object E, object U, object Utime)
        double GenerateIC(int nGrid,int iSeed,double z,double L,cosmo.csmContext * csm)
        void AddLinearSignal(int iGrid, int iSeed, double Lbox, double a, bool bFixed, float fPhase)
        double Read(string achInFile)
        void Write(string pszFileName,double dTime,bool bCheckpoint)
        void DomainDecomp(int iRung)
        void BuildTree(bool bNeedEwald)
        void Reorder()
        void Hostname()
        double LoadOrGenerateIC()
        void Simulate(double dTime,double dDelta,int iStartStep,int nSteps)
        void Simulate(double dTime)
        uint8_t Gravity(uint8_t uRungLo, uint8_t uRungHi,int iRoot1,int iRoot2,
                        double dTime,double dDelta,double dStep,double dTheta,
                        bool bKickClose,bool bKickOpen,bool bEwald,bool bGravStep,
                        int nPartRhoLoc,bool iTimeStepCrit)
        EphemeralMemory EphemeralMemoryGrid(int nGrid,int nCount);
        void addAnalysis(object callback,uint64_t per_particle, uint64_t per_process)
        # std::tuple<std::vector<uint64_t>,std::vector<float>,std::vector<float>,std::vector<float>> // nPk, fK, fPk, fPkAll
        void *MeasurePk(int iAssignment,int bInterlace,int nGrid,double a,int nBins)
        # std::tuple<std::vector<uint64_t>,std::vector<float>,std::vector<float>> // nPk, fK, fPk
        void *GridBinK(int nBins, int iGrid)
        void GridCreateFFT(int nGrid)
        void GridDeleteFFT()
        void AssignMass(int iAssignment,int iGrid,float fDelta)
        void DensityContrast(int nGrid,bool k)
        void WindowCorrection(int iAssignment,int iGrid);
        void Interlace(int iGridTarget,int iGridSource);
        void BispectrumSelect(int iGridTarget,int iGridSource,double kmin,double kmax);
        double BispectrumCalculate(int iGrid1,int iGrid2,int iGrid3);
        void OutputGrid(const char *filename, bool k, int iGrid, bool nParaWrite)
        void NewFof(double dTau,int nMinMembers)
        void GroupStats()
        void Smooth(double dTime,double dDelta,int iSmoothType,int bSymmetric,int nSmooth)
        int ReSmooth(double dTime,double dDelta,int iSmoothType,int bSymmetric)
        void OutASCII(const char *pszFile,int iType,int nDims,int iFileType)
        uint64_t CountSelected()
        void RecvArray(void *vBuffer,PKD_FIELD field,int iUnitSize,double dTime,bool bMarked)
        uint64_t SelBox(TinyVector[double,BLITZ3] center, TinyVector[double,BLITZ3] size,int setIfTrue,int clearIfFalse)
        uint64_t SelSphere(TinyVector[double,BLITZ3] r, double dRadius,int setIfTrue,int clearIfFalse)
        uint64_t SelCylinder(TinyVector[double,BLITZ3] dP1, TinyVector[double,BLITZ3] dP2, double dRadius, int setIfTrue, int clearIfFalse)
        uint64_t SelSpecies(uint64_t mSpecies,int setIfTrue,int clearIfFalse)
        void RsLoadIds(string filename,bool bAppend)
        void RsHaloLoadIds(string filename,bool bAppend)
        void RsSaveIds(string filename)
        void RsReorderIds()
        void RsExtract(string filename_template)

cdef public MSR *msr0 "PKDGRAV_msr0"

cdef extern from *:
    """
    struct MeasurePkStruct {
        std::vector<std::uint64_t> nPk;
        std::vector<float>    fK;
        std::vector<float>    fPk;
        std::vector<float>    fPkAll;
    };

    template<typename T>
    inline MeasurePkStruct UnpackMeasurePk(T t) {
        return {std::move(std::get<0>(t)), 
                std::move(std::get<1>(t)), 
                std::move(std::get<2>(t)), 
                std::move(std::get<3>(t))};
    }

    struct GridBinKStruct {
        std::vector<std::uint64_t> nPk;
        std::vector<float>    fK;
        std::vector<float>    fPk;
    };

    template<typename T>
    inline GridBinKStruct UnpackGridBinK(T t) {
        return {std::move(std::get<0>(t)), 
                std::move(std::get<1>(t)), 
                std::move(std::get<2>(t))};
    }
    struct ReadCheckpointStruct {
        double dTime;
        double dDelta;
        int64_t iStep;
        int64_t nSteps;
        int64_t nSizeParticle;
        int64_t nSizeNode;
    };
    template<typename T>
    inline ReadCheckpointStruct UnpackReadCheckpoint(T t) {
        return {std::move(std::get<0>(t)), 
                std::move(std::get<1>(t)), 
                std::move(std::get<2>(t)), 
                std::move(std::get<3>(t)), 
                std::move(std::get<4>(t)), 
                std::move(std::get<5>(t))};
    }
    """

    ctypedef struct MeasurePkStruct:
        vector[uint64_t]  nPk
        vector[float]     fK
        vector[float]     fPk
        vector[float]     fPkAll

    MeasurePkStruct UnpackMeasurePk(void *)

    ctypedef struct GridBinKStruct:
        vector[uint64_t]  nPk
        vector[float]     fK
        vector[float]     fPk

    GridBinKStruct UnpackGridBinK(void *)

    ctypedef struct ReadCheckpointStruct:
        double dTime
        double dDelta
        int64_t iStep
        int64_t nSteps
        int64_t nSizeParticle
        int64_t nSizeNode

    ReadCheckpointStruct UnpackReadCheckpoint(void *)

cdef inline tuple MeasurePk(int iAssignment,int bInterlace,int nGrid,double a,int nBins):
    cdef MeasurePkStruct result = UnpackMeasurePk(msr0.MeasurePk(iAssignment,bInterlace,nGrid,a,nBins))
    cdef:
        cnp.uint64_t[:]  nPk    = <cnp.uint64_t[:result.nPk.size()]>result.nPk.data()
        cnp.float32_t[:] fK     = <cnp.float32_t[:result.fK.size()]>result.fK.data()
        cnp.float32_t[:] fPk    = <cnp.float32_t[:result.fPk.size()]>result.fPk.data()
        cnp.float32_t[:] fPkAll = <cnp.float32_t[:result.fPkAll.size()]>result.fPkAll.data()
    return np.array(nPk), np.array(fK), np.array(fPk), np.array(fPkAll)

cdef inline tuple GridBinK(int nBins, iGrid):
    cdef GridBinKStruct result = UnpackGridBinK(msr0.GridBinK(nBins,iGrid))
    cdef:
        cnp.uint64_t[:]  nPk    = <cnp.uint64_t[:result.nPk.size()]>result.nPk.data()
        cnp.float32_t[:] fK     = <cnp.float32_t[:result.fK.size()]>result.fK.data()
        cnp.float32_t[:] fPk    = <cnp.float32_t[:result.fPk.size()]>result.fPk.data()
    return np.array(nPk), np.array(fK), np.array(fPk)

cdef inline read_checkpoint(filename: str,kwargs,species=None,classes=None,step=None,steps=None,time=None,delta=None,E=None,U=None,Utime=None):
    cdef ReadCheckpointStruct result = UnpackReadCheckpoint(msr0.ReadCheckpoint(filename.encode('utf-8'),kwargs,species,classes,step,steps,time,delta,E,U,Utime))
    return result.dTime,result.dDelta,result.iStep,result.nSteps,result.nSizeParticle,result.nSizeNode

# cpdef load(str filename)
cpdef save(str filename,double time=*)
cpdef domain_decompose(int rung=*)
# cpdef build_tree(bool ewald=*)
cpdef reorder()
#cpdef simulate()
cpdef measure_pk(int grid,int bins=*,double a=*,bool interlace=*,int order=*,double L=*)
cpdef fof(double tau,int minmembers=*)
cpdef smooth(SMOOTH_TYPE type,n=*,time=*,delta=*,symmetric=*,resmooth=*)
cpdef get_array(PKD_FIELD field,double time=*,bool marked=*)
cpdef write_array(filename,OUT_TYPE field)

cdef inline vector[PARTCLASS] new_partclass_vector():
    return vector[PARTCLASS]()

cdef inline double *a2d(a):
    cdef double[::1] v = a
    return &v[0]

cdef inline float *a2f(a):
    cdef float[::1] v = a
    return &v[0]

cdef inline void *a2f2(a):
    cdef float[:,:] v = a
    return &v[0,0]

cdef inline void *a2d2(a):
    cdef double[:,:] v = a
    return &v[0,0]

cdef inline void *a2u2(a):
    cdef uint64_t[:,:] v = a
    return &v[0,0]

cdef inline int64_t *a2i64(a):
    cdef int64_t[::1] v = a
    return &v[0]

cdef inline uint64_t *a2u64(a):
    cdef uint64_t[::1] v = a
    return &v[0]
