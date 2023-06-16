# distutils: language = c++
# cimport CMSR
# from CMSR cimport PARTCLASS, FIO_SPECIES, PKD_FIELD
# from CMSR cimport OUT_TYPE as CMSR_OUT_TYPE
from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.string cimport string
from libc.stdint cimport int64_t,uint64_t,uint8_t
import cython

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

cdef extern from "param.h":
    cdef struct prmContext:
        pass

include "pkd_parameters.pxi"

cdef extern from "master.h":
    cdef cppclass MSR:
        pkd_parameters parameters
        prmContext *prm

        uint64_t N
        # MSR() except +
        void testv(vector[PARTCLASS] &v)
        void Restart(int n, const char *baseName, int iStep, int nSteps, double dTime, double dDelta,
                    size_t nDark, size_t nGas, size_t nStar, size_t nBH,
                    double dEcosmo,double dUOld, double dTimeOld,
                    vector[PARTCLASS] &aClasses,object arguments,object specified)
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
        void MeasurePk(int iAssignment,int bInterlace,int nGrid,double a,int nBins,uint64_t *nPk,float *fK,float *fPk,float *fPkAll)
        void NewFof(double dTau,int nMinMembers)
        void GroupStats()
        uint64_t CountSelected()
        void RecvArray(void *vBuffer,PKD_FIELD field,int iUnitSize,double dTime,bool bMarked)
        uint64_t SelBox(TinyVector[double,BLITZ3] center, TinyVector[double,BLITZ3] size,int setIfTrue,int clearIfFalse)
        uint64_t SelSphere(TinyVector[double,BLITZ3] r, double dRadius,int setIfTrue,int clearIfFalse)
        uint64_t SelCylinder(TinyVector[double,BLITZ3] dP1, TinyVector[double,BLITZ3] dP2, double dRadius, int setIfTrue, int clearIfFalse)

cdef public MSR *msr0 "PKDGRAV_msr0"

cpdef restart(object arguments,object specified,list species,list classes,int n,str name,
    int step,int steps,double time,double delta,double E,double U,double Utime)
cpdef load(str filename)
cpdef save(str filename,double time=*)
cpdef domain_decompose(int rung=*)
cpdef build_tree(bool ewald=*)
cpdef reorder()
#cpdef simulate()
cpdef measure_pk(int grid,int bins=*,double a=*,bool interlace=*,int order=*,double L=*)
cpdef fof(double tau,int minmembers=*)
cpdef get_array(PKD_FIELD field,double time=*,bool marked=*)

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