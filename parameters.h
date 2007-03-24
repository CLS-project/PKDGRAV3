#ifndef PARAMETERS_HINCLUDED
#define PARAMETERS_HINCLUDED

#include "cosmo.h"

/*
** Don't even think about putting a pointer in here!!
*/
struct parameters {
    /*
    ** Parameters for PKDGRAV.
    */
    int nThreads;
    int bDiag;
    int bOverwrite;
    int bVWarnings;
    int bVStart;
    int bVStep;
    int bVRungStat;
    int bVDetails;
    int bPeriodic;
    int bParaRead;
    int bParaWrite;
    int bCannonical;
    int bStandard;
    int bDoublePos;
    int bGravStep;
    int bEpsAccStep;
    int bSqrtPhiStep;
    int bAccelStep; /* true if bEpsAccStep or bSqrtPhiStep */
    int bDensityStep;
    int iTimeStepCrit;
    int nPColl;
    int nTruncateRung;
    int bDoDensity;
    int bDodtOutput;
    int bDoRungOutput;
    int bDoGravity;
    int bAntiGrav;
    int nBucket;
    int iOutInterval;
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
    int nGrowMass;
    int iWallRunTime;
    int bPhysicalSoft;  
    int bSoftMaxMul;
    int bVariableSoft;
    int nSoftNbr;
    int bSoftByType;
    int bDoSoftOutput;
    double dEta;
    double dExtraStore;
    double dSoft;
    double dSoftMax;
    double dDelta;
    double dEwCut;
    double dEwhCut;
    double dTheta;
    double dTheta2;
    double daSwitchTheta;
    double dPeriod;
    double dxPeriod;
    double dyPeriod;
    double dzPeriod;
    CSM csm;
    double dRedTo;
    double dCentMass;
    char achDigitMask[256];
    char achInFile[256];
    char achOutName[256];
    char achDataSubPath[256];
    double dGrowDeltaM;
    double dGrowStartT;
    double dGrowEndT;
    double dFracNoTreeSqueeze;
    double dFracNoDomainDecomp;
    double dFracNoDomainDimChoice;
    int    bRungDD;
    double dRungDDWeight;
    /*
    ** Additional parameters for group finding.
    */
    int	nFindGroups;	
    int	nMinMembers;	
    double dTau;
    int bTauAbs;
    int	nBins;	
    int	bUsePotmin;	
    int	nMinProfile;	
    double fBinsRescale;
    double fContrast;
    double Delta;
    double binFactor;
    int bLogBins; 
#ifdef RELAXATION
    int	bTraceRelaxation;	
#endif
    };

#endif
