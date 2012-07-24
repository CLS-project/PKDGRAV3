#ifndef PARAMETERS_HINCLUDED
#define PARAMETERS_HINCLUDED

#include "cosmo.h"

#ifdef PLANETS
typedef struct {
    int iOutcomes;
    double dDensity;
    double dBounceLimit;
    int iBounceOption;
    double dEpsN;
    double dEpsT;
    int bFixCollapse;
    } COLLISION_PARAMS;
#endif

/*
** Don't even think about putting a pointer in here!!
*/
struct parameters {
    /*
    ** Parameters for PKDGRAV.
    */
    int nThreads;
    int bDiag;
    int bNoGrav;
    int bOverwrite;
    int bVWarnings;
    int bVStart;
    int bVStep;
    int bVRungStat;
    int bVDetails;
    int bPeriodic;
    int bRestart;
    int bParaRead;
    int bParaWrite;
    int bStandard;
    int iCompress;
    int bHDF5;
    int bDoublePos;
    int bGravStep;
    int bEpsAccStep;
    int bSqrtPhiStep;
    int bAccelStep; /* true if bEpsAccStep or bSqrtPhiStep */
    int bDensityStep;
    int iTimeStepCrit;
    int nPartRhoLoc;
    int nPartColl;
    int nTruncateRung;
    int bDoDensity;
#ifdef USE_PNG
    int nPNGResolution;
#endif
    int bDoRungOutput;
    int bDoRungDestOutput;
    int bDoGravity;
    int bHermite;
    int bAarsethStep;
    int nBucket;
    int nGroup;
    int n2min;
    int iOutInterval;
    int iCheckInterval;
    int iLogInterval;
    int iOrder;
    int bEwald;
    int iEwOrder;
    int nReplicas;
    int iStartStep;
    int nSteps;
    int nSmooth;
    int iMaxRung;
    int nRungVeryActive;
    int nPartVeryActive;
    int nTreeBitsLo;
    int nTreeBitsHi;
    int iWallRunTime;
    int bPhysicalSoft;
    int bSoftMaxMul;
    int nSoftNbr;
    int bSoftByType;
    int bDoSoftOutput;
    int bDoAccOutput;
    int bDoPotOutput;
    int iCacheSize;
    int iWorkQueueSize;
    int iCUDAQueueSize;
    int bAddDelete;
    /* BEGIN Gas Parameters */
    int bDoGas;
    int bGasAdiabatic;
    int bGasIsothermal;
    int bGasCooling;
    int bInitTFromCooling;
    int bStarForm;
    int bFeedback;
    int iViscosityLimiter;
    int iDiffusion;
    int iRungCoolTableUpdate;
    double dEtaCourant;
    double dEtaUDot;
    double dConstAlpha;
    double dConstBeta;
    double dConstGamma;
    double dMeanMolWeight;
    double dGasConst;
    double dTuFac;
    double dMsolUnit;
    double dKpcUnit;
    double ddHonHLimit;
    double dKBoltzUnit;
    double dGmPerCcUnit;
    double dComovingGmPerCcUnit;
    double dErgPerGmUnit;
    double dSecUnit;
    double dKmPerSecUnit;
    double dhMinOverSoft;
    double dMetalDiffusionCoeff;
    double dThermalDiffusionCoeff;
#ifdef COOLING
    COOLPARAM CoolParam;
#endif
    /* StarForm and Feedback */
    double SFdEfficiency;
    double SFdTMax;
    double SFdPhysDenMin;
    double SFdComovingDenMin;
    double SFdESNPerStarMass;
    double SFdtCoolingShutoff;

    double SFdtFeedbackDelay;
    double SFdMassLossPerStarMass;
    double SFdZMassPerStarMass;
    double SFdInitStarMass;
    double SFdMinGasMass;
    double SFdvFB;
    int SFbdivv;
    
    /* END Gas Parameters */
    double dEta;
    double dExtraStore;
    double dSoft;
    double dSoftMax;
    double dDelta;
    double dEwCut;
    double dEwhCut;
    double dTheta;
    double dTheta2;
    double dThetaMax;
    double daSwitchTheta;
    double dPeriod;
    double dxPeriod;
    double dyPeriod;
    double dzPeriod;
    double dPreFacRhoLoc;
    double dFacExcludePart;
    double dEccFacMax;
    CSM csm;
    double dRedTo;
    double dRedFrom;
    double dCentMass;
    char achDigitMask[256];
    char achInFile[256];
    char achOutName[256];
    char achOutPath[256];
    char achIoPath[256];
    char achDataSubPath[256];
    char achOutTypes[256];
    char achCheckTypes[256];
#ifdef USE_PYTHON
    char achScriptFile[256];
#endif
    double dGrowDeltaM;
    double dGrowStartT;
    double dGrowEndT;
    double dFracNoDomainDecomp;
    double dFracNoDomainRootFind;
    double dFracNoDomainDimChoice;
    /*
    ** Additional parameters for group finding.
    */
    int	bFindGroups;
    int	nMinMembers;
    double dTau;
    double dVTau;
    int bTauAbs;
    int	nBins;
    int	iCenterType;
    double binFactor;
    double fMinRadius;
    int bLogBins;
    int	bTraceRelaxation;
#ifdef PLANETS
    int bHeliocentric;
    int bCollision;
    int iCollLogOption;
    char achCollLog[256];
    COLLISION_PARAMS CP;
#ifdef SYMBA
    int bSymba;
#endif
#endif /* PLANETS */

#ifdef USE_MDL_IO
    int nIO;
#endif

    /* IC Generation */
    int bWriteIC;
#ifdef USE_GRAFIC
    double h;
    double dBoxSize;
    int nGrid;
    int iSeed;
#endif

#ifdef USE_LUSTRE
    int nStripeSize;
    int nStripeCount;
#endif

    int nDomainRungs;
    int iDomainMethod;

    /*
    ** Memory models.  Other parameters can force these to be set.
    */
    int bMemAcceleration;
    int bMemVelocity;
    int bMemPotential;
    int bMemGroups;
    int bMemMass;
    int bMemSoft;
    int bMemHermite;
    int bMemRelaxation;
    int bMemVelSmooth;
    int bMemPsMetric;
    int bMemNewDD ;
    int bMemNodeMoment;
    int bMemNodeAcceleration;
    int bMemNodeVelocity;
    int bMemNodeSphBounds;
    int bMemNodeBnd;
    int bMemNodeVBnd;
    };

#endif
