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

#include <stdint.h>
#include <signal.h>
#include <time.h>
#include <vector>
#include <Python.h>

#include "param.h"
#include "pst.h"
#include "mdl.h"
#include "parameters.h"

#define MSR_INIT_E		1
#define MSR_STEP_E		0

extern time_t timeGlobalSignalTime;
extern int bGlobalOutput;

struct MSRINSTANCE {
    PyObject_HEAD
    class MSR *msr;
    };

class MSR {
protected:
    const PST pst;
    mdl::mdlClass *mdl;
    bool bVDetails;
public:
    explicit MSR(MDL mdl,PST pst) : pst(pst), mdl(static_cast<mdl::mdlClass *>(mdl)), bVDetails(false) {}
    ~MSR();
public:
    int Python(int argc, char *argv[]);
    int ValidateParameters();
    void Hostname();
    void MemStatus();
    int GetLock();
    double LoadOrGenerateIC();
    void Simulate(double dTime,double dDelta,int iStartStep,int nSteps);
    void Simulate(double dTime);
public:
    // Parameters
    bool      wasParameterSpecified(const char *name) const;
    bool      getParameterBoolean(  const char *name) const;
    double    getParameterDouble(   const char *name) const;
    long long getParameterLongLong( const char *name) const;
    void      setParameter(         const char *name,bool v,     int bSpecified=false);
    void      setParameter(         const char *name,double v,   int bSpecified=false);
    void      setParameter(         const char *name,long long v,int bSpecified=false);

    size_t getLocalGridMemory(int nGrid);

    // I/O and IC Generation
    double GenerateIC();
    void Restart(int n, const char *baseName, int iStep, int nSteps, double dTime, double dDelta);
    double Read(const char *achInFile);
    void Checkpoint(int iStep, int nSteps, double dTime, double dDelta);
    void Write(const char *pszFileName,double dTime,int bCheckpoint);
    void OutArray(const char *pszFile,int iType,int iFileType);
    void OutArray(const char *pszFile,int iType);
    void OutVector(const char *pszFile,int iType,int iFileType);
    void OutVector(const char *pszFile,int iType);
    void Output(int iStep, double dTime, double dDelta, int bCheckpoint);

    void RecvArray(void *vBuffer,int field,int iUnitSize,double dTime,bool bMarked=false);

    // Particle order, domains, trees
    void Reorder();
    void DomainDecomp(int iRung=0);
    void BuildTree(int bNeedEwald);
    void BuildTreeFixed(int bNeedEwald,uint8_t uRungDD);
    void BuildTreeActive(int bNeedEwald,uint8_t uRungDD);
    void BuildTreeMarked(int bNeedEwald=0);

    // Gravity
    uint8_t Gravity(uint8_t uRungLo, uint8_t uRungHi,int iRoot1,int iRoot2,
    	double dTime,double dDelta,double dStep,double dTheta,
    	int bKickClose,int bKickOpen,int bEwald,int bGravStep,int nPartRhoLoc,int iTimeStepCrit,int nGroup,SPHOptions SPHoptions);

    // Analysis
    void Smooth(double dTime,double dDelta,int iSmoothType,int bSymmetric,int nSmooth);
    void ReSmooth(double dTime,double dDelta,int iSmoothType,int bSymmetric);
    void NewFof(double exp);
    void Hop(double dTime,double dDelta);
    void GroupStats();
    void HopWrite(const char *fname);
    void MeasurePk(int iAssignment,int bInterlace,int nGrid,double a,int nBins,uint64_t *nPk,float *fK,float *fPk,float *fPkAll);
    void AssignMass(int iAssignment=4,int iGrid=0,float fDelta=0.0f);
    void DensityContrast(int nGrid,bool k=true);
    void WindowCorrection(int iAssignment,int iGrid);
    void Interlace(int iGridTarget,int iGridSource);
    void AddLinearSignal(int iGrid, int iSeed, double Lbox, double a, bool bFixed=false, float fPhase=0);
    void GridBinK(int nBins, int iGrid,uint64_t *nPk,float *fK,float *fPk);
    void BispectrumSelect(int iGridTarget,int iGridSource,double kmin,double kmax);
    double BispectrumCalculate(int iGrid1,int iGrid2,int iGrid3);
    void GridCreateFFT(int nGrid);
    void GridDeleteFFT();

private:
    typedef struct {
	double dFrac;       /* Fraction of particles in each bin */
	uint64_t nTotal;    /* Total number of particles in the range */
	uint64_t nInner;    /* Number inside minimum radius */
	uint64_t nTarget;   /* Target number of particles */
	uint64_t nSelected;
	MSR *msr;
	} SPHERECTX;
    static double countSphere(double r,void *vctx);
    static void profileRootFind( double *dBins, int lo, int hi, int nAccuracy, SPHERECTX *ctx );

    typedef struct {
	double rMiddle;
	total_t nTarget;   /* Target number of particles */
	MSR *msr;
	} SHELLCTX;
    static double countShell(double rInner,void *vctx);

public:
    struct msr_analysis_callback {
	PyObject *callback;
	struct MSRINSTANCE *msr;
	PyObject *memory;
	explicit msr_analysis_callback(PyObject *callback,MSRINSTANCE *msr,PyObject *memory) {
	    this->callback = callback;
	    this->msr = msr;
	    this->memory = memory;
	    Py_INCREF(callback);
	    Py_INCREF(memory);
	    }
//	~msr_analysis_callback() {
//	    Py_DECREF(callback);
//	    }
	};
    void addAnalysis(PyObject *callback,MSRINSTANCE *msr,PyObject *memory);
    void runAnalysis(int iStep,double dTime);
protected:
    std::list<msr_analysis_callback> analysis_callbacks;

// Parameters from Python interface
private:
    int64_t     getScalarInteger(const char *name, PyObject *v);
    double      getScalarNumber(const char *name, PyObject *v);
    std::string getScalarString(const char *name, PyObject *v);
public: // should be private
    PyObject *arguments=nullptr, *specified=nullptr;
public:
    int64_t                  getScalarInteger(const char *name);
    double                   getScalarNumber(const char *name);
    std::string              getScalarString(const char *name);
    std::vector<int64_t>     getVectorInteger(const char *name);
    std::vector<double>      getVectorNumber(const char *name);
    std::vector<std::string> getVectorString(const char *name);
    bool setParameters(PyObject *kwobj,bool bIgnoreUnknown=false);

    PRM prm;
    LCL lcl;
    blitz::TinyVector<double,3> fCenter;
public:
    /*
    ** Parameters.
    */
    struct parameters param;
    CSM csm;
    /*
    ** Other stuff...
    */
    int nThreads;
    uint64_t N;
    uint64_t nDark;
    uint64_t nGas;
    uint64_t nStar;
    uint64_t nMaxOrder;		/* Order number of last particle */
    int iCurrMaxRung;

    /*
    ** Tree moments (for Ewald)
    */
    MOMC momTreeRoot[IRUNGMAX+1];
    double momTreeCom[IRUNGMAX+1][3];

    /*
    ** Comoving coordinate variables.
    */
    double dEcosmo;
    double dUOld;
    double dTimeOld;
    /*
    ** Redshift output points.
    */
    std::vector<double> dOutTimes;
    int iOut;
    // int nMaxOuts;
    // int nOuts;
    // double *pdOutTime;

    /*
     * Domain Decomposition Done
     */
    std::vector<uint64_t> nRung;
    int iRungDD, iRungDT;
    int iLastRungRT,iLastRungDD;
    uint64_t nActive;
    int nGroups;
    int nBins;
    int bAntiGrav;

    long lStart; /* starting time of job */
    long lPrior; /* starting time of last step */

    /* Gas */
    double dTuFac;

    /* Values for a restore from checkpoint */
    double dCheckpointTime;
    int iCheckpointStep;
    int nCheckpointThreads;
    char achCheckpointName[PST_FILENAME_SIZE];
    int nCheckpointClasses;
    PARTCLASS aCheckpointClasses[PKD_MAX_CLASSES];
protected:
    static double Time();
    static void Leader();
    static void Trailer();
    static void MakePath(const char *dir,const char *base,char *path);

    double getTime(double dExpansion); // Return simulation time
    double getVfactor(double dTime);
    bool getDeltaSteps(double dTime,int iStartStep,double &dDelta,int &nSteps);

    const char *OutName() const { return param.achOutName;}
    //int Steps()           const { return param.nSteps; }
    int LogInterval()     const { return param.iLogInterval; }
    int OutInterval()     const { return param.iOutInterval; }
    int CheckInterval()   const { return param.iCheckInterval; }
    double Soft()         const { return param.dSoft; }
    int DoDensity()       const { return param.bDoDensity; }
    int DoGas()           const { return param.bDoGas; }
    int DoGravity()       const { return param.bDoGravity; }
    double Eta()          const { return param.dEta; }
    int MaxRung()         const { return param.iMaxRung; }
    int Comove()          const { return csm->val.bComove; }
    uint64_t MaxOrder()   const { return nMaxOrder; }
    int CurrMaxRung()     const { return iCurrMaxRung; }

    std::string BuildName(const char *path,int iStep,const char *type="");
    std::string BuildName(int iStep,const char *type=""); // With achOutPath
    std::string BuildIoName(int iStep,const char *type="");
    std::string BuildCpName(int iStep,const char *type="");

    void ReadOuts(double dTime,double dDelta);
    void msrprintf(const char *Format, ... ) const;
    void Exit(int status);
    uint64_t getMemoryModel();
    void InitializePStore(uint64_t *nSpecies,uint64_t mMemoryModel);
    int CheckForStop(const char *achStopFile);
    int CheckForOutput(int iStep,int nSteps,double dTime,int *pbDoCheckpoint,int *pbDoOutput);
    bool OutTime(double dTime);
    void SetClasses();
    void SwapClasses(int id);
    void OneNodeRead(struct inReadFile *in, FIO fio);
    void AllNodeWrite(const char *pszFileName, double dTime, double dvFac, int bDouble);
    uint64_t CalcWriteStart();
    void SwitchTheta(double);
    double getTheta(double dTime);
    double SwitchDelta(double dTime,double dDelta,int iStep,int nSteps);
    void InitCosmology();
    void BuildTree(int bNeedEwald,uint32_t uRoot,uint32_t utRoot);
    void ActiveRung(int iRung, int bGreater);
    void ActiveOrder();
    void CalcBound(Bound &bnd);
    void CalcBound();
    void GetNParts();
    double AdjustTime(double aOld, double aNew);
    void UpdateSoft(double dTime);
    int GetParticles(int nIn, uint64_t *ID, struct outGetParticles *out);
    void OutputOrbits(int iStep,double dTime);
    double TotalMass();
    void LightConeOpen(int iStep);
    void LightConeClose(int iStep);
    void LightConeVel();
#ifdef MDL_FFTW
    void SetLinGrid(double dTime, double dDelta, int nGrid, int bKickClose, int bKickOpen);
    void LinearKick(double dTime, double dDelta, int bKickClose, int bKickOpen);
#endif
    void CalcEandL(int bFirst,double dTime,double *E,double *T,double *U,double *Eth,double *L,double *F,double *W);
    void Drift(double dTime,double dDelta,int iRoot);

    void SmoothSetSMF(SMF *smf, double dTime, double dDelta);
    void ZeroNewRung(uint8_t uRungLo, uint8_t uRungHi, int uRung);
    void KickKDKOpen(double dTime,double dDelta,uint8_t uRungLo,uint8_t uRungHi);
    void KickKDKClose(double dTime,double dDelta,uint8_t uRungLo,uint8_t uRungHi);
    void UpdateRung(uint8_t uRung);
    void AccelStep(uint8_t uRungLo,uint8_t uRungHi,double dTime,double dDelta);
    void SphStep(uint8_t uRungLo,uint8_t uRungHi,double dTime,double dDelta);
    void DensityStep(uint8_t uRungLo,uint8_t uRungHi,double dTime,double dDelta);

    void FastGasPhase1(double dTime,double dDelta,int iSmoothType);
    void FastGasPhase2(double dTime,double dDelta,int iSmoothType);
    void CoolSetup(double dTime);
    void Cooling(double dTime,double dStep,int bUpdateState, int bUpdateTable,int bInterateDt);
    void AddDelParticles();
    void StarForm(double dTime, double dDelta, int iRung);
    void InitSph(double dTime,double dDelta);
    void Sph(double dTime, double dDelta, double dStep);
    uint64_t CountDistance(double dRadius2Inner, double dRadius2Outer);

    void Initialize();
    void writeParameters(const char *baseName,int iStep,int nSteps,double dTime,double dDelta);
    void OutASCII(const char *pszFile,int iType,int nDims,int iFileType);
    void DomainDecompOld(int iRung);

    void SaveParameters();
    int CountRungs(uint64_t *nRungs);
    void SetSoft(double);

    void MeasureLinPk(int nGridLin,double a,double dBoxSize, uint64_t *nPk,float *fK,float *fPk);
    void OutputPk(int iStep,double dTime);
    void OutputLinPk(int iStep, double dTime);

    int NewTopStepKDK(
	double &dTime,	/* MODIFIED: Current simulation time */
	double dDelta,
	double dTheta,
	int nSteps,
	int bDualTree,      /* Should be zero at rung 0! */
	uint8_t uRung,	/* Rung level */
	double *pdStep,	/* Current step */
	uint8_t *puRungMax,int *pbDoCheckpoint,int *pbDoOutput,int *pbNeedKickOpen);
    void TopStepKDK(
		    double dStep,	/* Current step */
		    double dTime,	/* Current time */
		    double dDelta,	/* Time step */
		    double dTheta,
		    int iRung,		/* Rung level */
		    int iKickRung,	/* Gravity on all rungs from iRung
					    to iKickRung */
		    int iAdjust);		/* Do an adjust? */

    void InitRelaxation();
    void Relaxation(double dTime,double deltaT,int iSmoothType,int bSymmetric);
    void CalcDistance(const double *dCenter, double dRadius );
    void CalcCOM(const double *dCenter, double dRadius,
		double *com, double *vcm, double *L, double *M);
    void CalcMtot(double *M, uint64_t *N);
    void SetSPHoptions();
    void TreeUpdateFlagBounds(int bNeedEwald,uint32_t uRoot,uint32_t utRoot,SPHOptions SPHoptions);

public:
    void Profile(
	const PROFILEBIN **pBins, int *pnBins, double *r,
	double dMinRadius, double dLogRadius, double dMaxRadius,
	int nPerBin, int nBins, int nAccuracy );
    void OutputGrid(const char *filename, bool k=false, int iGrid=0, int nParaWrite=0);

    uint64_t CountSelected();
    uint64_t SelSpecies(uint64_t mSpecies,bool setIfTrue=true,bool clearIfFalse=true);
    uint64_t SelAll(bool setIfTrue=true,bool clearIfFalse=true);
    uint64_t SelGas(bool setIfTrue=true,bool clearIfFalse=true);
    uint64_t SelStar(bool setIfTrue=true,bool clearIfFalse=true);
    uint64_t SelDark(bool setIfTrue=true,bool clearIfFalse=true);
    uint64_t SelDeleted(bool setIfTrue=true,bool clearIfFalse=true);
    uint64_t SelBlackholes(bool setIfTrue=true,bool clearIfFalse=true);
    uint64_t SelGroup(int iGroup,bool setIfTrue=true,bool clearIfFalse=true);
    uint64_t SelMass(double dMinMass,double dMaxMass,int setIfTrue,int ClearIfFalse);
    uint64_t SelById(uint64_t idStart,uint64_t idEnd,int setIfTrue,int clearIfFalse);
    uint64_t SelPhaseDensity(double dMinPhaseDensity,double dMaxPhaseDensity,int setIfTrue,int clearIfFalse);
    uint64_t SelBox(double *dCenter, double *dSize,bool setIfTrue=true,bool clearIfFalse=true);
    uint64_t SelSphere(double *r, double dRadius,int setIfTrue,int clearIfFalse);
    uint64_t SelCylinder(double *dP1, double *dP2, double dRadius, int setIfTrue, int clearIfFalse );
    };
#endif
