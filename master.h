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

class MSR {
protected:
    const PST pst;
    const MDL mdl;
public:
    explicit MSR(MDL mdl,PST pst) : pst(pst), mdl(mdl) {}
    ~MSR();
public:
    int Python(int argc, char *argv[]);
    int ValidateParameters();
    void Hostname();
    void MemStatus();
    int GetLock();
    double LoadOrGenerateIC();
    void Simulate(double dTime,int iStartStep);
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

    // I/O and IC Generation
    double GenerateIC();
    void Restart(int n, const char *baseName, int iStep, double dTime);
    double Read(const char *achInFile);
    void Checkpoint(int iStep, double dTime);
    void Write(const char *pszFileName,double dTime,int bCheckpoint);
    void OutArray(const char *,int);
    void OutVector(const char *,int);
    void Output(int iStep, double dTime, int bCheckpoint);

    // Particle order, domains, trees
    void Reorder();
    void DomainDecomp(int iRung=0);
    void BuildTree(int bNeedEwald);
    void BuildTreeFixed(int bNeedEwald,uint8_t uRungDD);
    void BuildTreeActive(int bNeedEwald,uint8_t uRungDD);
    void BuildTreeMarked();

    // Gravity
    uint8_t Gravity(uint8_t uRungLo, uint8_t uRungHi,int iRoot1,int iRoot2,
	double dTime,double dStep,int bKickClose,int bKickOpen,int bEwald,int bGravStep,int nPartRhoLoc,int iTimeStepCrit,
	int nGroup,int *piSec,uint64_t *pnActive);

    // Analysis
    void Smooth(double dTime,int iSmoothType,int bSymmetric,int nSmooth);
    void ReSmooth(double dTime,int iSmoothType,int bSymmetric);
    void NewFof(double exp);
    void Hop(double exp);
    void GroupStats();
    void HopWrite(const char *fname);
    void MeasurePk(int iAssignment,int bInterlace,int nGrid,double a,int nBins,uint64_t *nPk,float *fK,float *fPk,float *fPkAll);

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
    PRM prm;
    LCL lcl;
    double fCenter[6];
    /*
    ** Parameters.
    */
    PyObject *arguments, *specified;
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
    double dThetaMin;

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
    double getTime(double dExpansion, double *dvFac);

    const char *OutName() const { return param.achOutName;}
    int Steps()           const { return param.nSteps; }
    double Delta()        const { return param.dDelta; }
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

    char *BuildName(char *achFile,int iStep,char *defaultPath);
    char *BuildName(char *achFile,int iStep);
    char *BuildIoName(char *achFile,int iStep);
    void ReadOuts(double dTime);
    void msrprintf(const char *Format, ... ) const;
    void Exit(int status);
    uint64_t getMemoryModel();
    void InitializePStore(uint64_t *nSpecies);
    int CheckForStop(const char *achStopFile);
    int CheckForOutput(int iStep,double dTime,int *pbDoCheckpoint,int *pbDoOutput);
    bool OutTime(double dTime);
    void SetClasses();
    void SwapClasses(int id);
    void OneNodeRead(struct inReadFile *in, FIO fio);
    void AllNodeWrite(const char *pszFileName, double dTime, double dvFac, int bDouble);
    uint64_t CalcWriteStart();
    void SwitchTheta(double);
    double SwitchDelta(double dTime,int iStep);
    void InitCosmology();
    void SetParameters();
    void BuildTree(int bNeedEwald,uint32_t uRoot,uint32_t utRoot);
    void ActiveRung(int iRung, int bGreater);
    void ActiveOrder();
    void CalcBound(BND *pbnd);
    void CalcVBound(BND *pbnd);
    void GetNParts();
    void ScaleVel(double dvFac);
    double AdjustTime(double aOld, double aNew);
    void UpdateSoft(double dTime);
    int GetParticles(int nIn, uint64_t *ID, struct outGetParticles *out);
    void OutputOrbits(int iStep,double dTime);
    void GridCreateFFT(int nGrid);
    void GridDeleteFFT();
    double TotalMass();
    void LightConeOpen(int iStep);
    void LightConeClose(int iStep);
    void LightConeVel();
#ifdef MDL_FFTW
    void AssignMass(int iAssignment,int nGrid);
    void SetLinGrid(double dTime, int nGrid, int bKickClose, int bKickOpen);
    void LinearKick(double dTime, int bKickClose, int bKickOpen);
#endif
    void CalcEandL(int bFirst,double dTime,double *E,double *T,double *U,double *Eth,double *L,double *F,double *W);
    void Drift(double dTime,double dDelta,int iRoot);

    void SmoothSetSMF(SMF *smf, double dTime);
    void ZeroNewRung(uint8_t uRungLo, uint8_t uRungHi, int uRung);
    void KickKDKOpen(double dTime,double dDelta,uint8_t uRungLo,uint8_t uRungHi);
    void KickKDKClose(double dTime,double dDelta,uint8_t uRungLo,uint8_t uRungHi);
    void UpdateRung(uint8_t uRung);
    void AccelStep(uint8_t uRungLo,uint8_t uRungHi,double dTime);
    void SphStep(uint8_t uRungLo,uint8_t uRungHi,double dTime);
    void DensityStep(uint8_t uRungLo,uint8_t uRungHi,double dTime);

    void FastGasPhase1(double dTime,int iSmoothType);
    void FastGasPhase2(double dTime,int iSmoothType);
    void CoolSetup(double dTime);
    void Cooling(double dTime,double dStep,int bUpdateState, int bUpdateTable,int bInterateDt);
    void AddDelParticles();
    void StarForm(double dTime, int iRung);
    void InitSph(double dTime);
    void Sph(double dTime, double dStep);
    uint64_t CountDistance(double dRadius2Inner, double dRadius2Outer);

    int Initialize();
    void writeParameters(const char *baseName,int iStep,double dTime);
    void OutASCII(const char *pszFile,int iType,int nDims);
    void DomainDecompOld(int iRung);

    void SaveParameters();
    int CountRungs(uint64_t *nRungs);
    void SetSoft(double);

    void MeasureLinPk(int nGridLin,double a,double dBoxSize, uint64_t *nPk,float *fK,float *fPk);
    void OutputPk(int iStep,double dTime);
    void OutputLinPk(int iStep, double dTime);

    int NewTopStepKDK(
	int bDualTree,      /* Should be zero at rung 0! */
	uint8_t uRung,	/* Rung level */
	double *pdStep,	/* Current step */
	double *pdTime,	/* Current time */
	uint8_t *puRungMax,int *piSec,int *pbDoCheckpoint,int *pbDoOutput,int *pbNeedKickOpen);
    void TopStepKDK(
		    double dStep,	/* Current step */
		    double dTime,	/* Current time */
		    double dDelta,	/* Time step */
		    int iRung,		/* Rung level */
		    int iKickRung,	/* Gravity on all rungs from iRung
					    to iKickRung */
		    int iAdjust,		/* Do an adjust? */
		    double *pdActiveSum,
		    int *piSec);

    void InitRelaxation();
    void Relaxation(double dTime,double deltaT,int iSmoothType,int bSymmetric);
    void CalcDistance(const double *dCenter, double dRadius );
    void CalcCOM(const double *dCenter, double dRadius,
		double *com, double *vcm, double *L, double *M);
    void Profile(
	const PROFILEBIN **pBins, int *pnBins, double *r,
	double dMinRadius, double dLogRadius, double dMaxRadius,
	int nPerBin, int nBins, int nAccuracy );

    void SelAll();
    void SelGas();
    void SelStar();
    void SelDeleted();
    uint64_t SelMass(double dMinMass,double dMaxMass,int setIfTrue,int ClearIfFalse);
    uint64_t SelById(uint64_t idStart,uint64_t idEnd,int setIfTrue,int clearIfFalse);
    uint64_t SelPhaseDensity(double dMinPhaseDensity,double dMaxPhaseDensity,int setIfTrue,int clearIfFalse);
    uint64_t SelBox(double *dCenter, double *dSize,int setIfTrue,int clearIfFalse);
    uint64_t SelSphere(double *r, double dRadius,int setIfTrue,int clearIfFalse);
    uint64_t SelCylinder(double *dP1, double *dP2, double dRadius, int setIfTrue, int clearIfFalse );
    };
#endif
