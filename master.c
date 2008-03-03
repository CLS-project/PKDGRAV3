#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#define _LARGEFILE_SOURCE
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> /* for unlink() */
#include <stddef.h>
#include <string.h>
#include <inttypes.h>
#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif
#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>

#ifdef HAVE_SYS_PARAM_H
#include <sys/param.h> /* for MAXHOSTNAMELEN, if available */
#endif
#include <rpc/types.h>
#include <rpc/xdr.h>

#include "master.h"
#include "tipsydefs.h"
#include "outtype.h"
#include "smoothfcn.h"
#ifdef USE_MDL_IO
#include "io.h"
#endif
#include "ssio.h"

#ifdef USE_BSC
#include "mpitrace_user_events.h"
#endif

#define LOCKFILE ".lockfile"	/* for safety lock */
#define STOPFILE "STOP"			/* for user interrupt */

#define NEWTIME
#ifdef NEWTIME 
double msrTime() {
    struct timeval tv;
    struct timezone tz;

    tz.tz_minuteswest=0; 
    tz.tz_dsttime=0;
    gettimeofday(&tv,NULL);
    return (tv.tv_sec+(tv.tv_usec*1e-6));
    }
#else
double msrTime() {
    return (1.0*time(0));
    }
#endif

void _msrLeader(void) {
    puts("pkdgrav2.2 Joachim Stadel & Doug Potter Sept 2007");
    puts("USAGE: pkdgrav2 [SETTINGS | FLAGS] [SIM_FILE]");
    puts("SIM_FILE: Configuration file of a particular simulation, which");
    puts("          includes desired settings and relevant input and");
    puts("          output files. Settings specified in this file override");
    puts("          the default settings.");
    puts("SETTINGS");
    puts("or FLAGS: Command line settings or flags for a simulation which");
    puts("          will override any defaults and any settings or flags");
    puts("          specified in the SIM_FILE.");
    }


void _msrTrailer(void)
    {
    puts("(see man page for more information)");
    }


void _msrExit(MSR msr,int status)
    {
    MDL mdl=msr->mdl;

    msrFinish(msr);
    mdlFinish(mdl);
    exit(status);
    }


void
_msrMakePath(const char *dir,const char *base,char *path)
    {
    /*
    ** Prepends "dir" to "base" and returns the result in "path". It is the
    ** caller's responsibility to ensure enough memory has been allocated
    ** for "path".
    */

    if (!path) return;
    path[0] = 0;
    if (dir) {
	strcat(path,dir);
	strcat(path,"/");
	}
    if (!base) return;
    strcat(path,base);
    }


void msrInitialize(MSR *pmsr,MDL mdl,int argc,char **argv)
    {
    MSR msr;
    int i,j,ret;
    int id,nDigits;
    struct inSetAdd inAdd;
    struct inLevelize inLvl;
    struct inGetMap inGM;

    msr = (MSR)malloc(sizeof(struct msrContext));
    assert(msr != NULL);
    msr->mdl = mdl;
    msr->pst = NULL;
    msr->lcl.pkd = NULL;
    *pmsr = msr;
    csmInitialize(&msr->param.csm);
    /*
    ** Now setup for the input parameters.
    **
    ** NOTE: nThreads & bDiag are parsed here, but the actual values are
    ** read from the command line via mdlInitialize(). This means the
    ** values of nThreads & bDiag read by prmAddParam() are ignored!
    */
    prmInitialize(&msr->prm,_msrLeader,_msrTrailer);
    msr->param.nThreads = 1;
    prmAddParam(msr->prm,"nThreads",1,&msr->param.nThreads,sizeof(int),"sz",
		"<nThreads>");
    msr->param.bDiag = 0;
    prmAddParam(msr->prm,"bDiag",0,&msr->param.bDiag,sizeof(int),"d",
		"enable/disable per thread diagnostic output");
    msr->param.bOverwrite = 0;
    prmAddParam(msr->prm,"bOverwrite",0,&msr->param.bOverwrite,sizeof(int),
		"overwrite","enable/disable overwrite safety lock = +overwrite");
    msr->param.bVWarnings = 1;
    prmAddParam(msr->prm,"bVWarnings",0,&msr->param.bVWarnings,sizeof(int),
		"vwarnings","enable/disable warnings = +vwarnings");
    msr->param.bVStart = 1;
    prmAddParam(msr->prm,"bVStart",0,&msr->param.bVStart,sizeof(int),
		"vstart","enable/disable verbose start = +vstart");
    msr->param.bVStep = 1;
    prmAddParam(msr->prm,"bVStep",0,&msr->param.bVStep,sizeof(int),
		"vstep","enable/disable verbose step = +vstep");
    msr->param.bVRungStat = 1;
    prmAddParam(msr->prm,"bVRungStat",0,&msr->param.bVRungStat,sizeof(int),
		"vrungstat","enable/disable rung statistics = +vrungstat");
    msr->param.bVDetails = 0;
    prmAddParam(msr->prm,"bVDetails",0,&msr->param.bVDetails,sizeof(int),
		"vdetails","enable/disable verbose details = +vdetails");
    nDigits = 5;
    prmAddParam(msr->prm,"nDigits",1,&nDigits,sizeof(int),"nd",
		"<number of digits to use in output filenames> = 5");
    msr->param.bPeriodic = 0;
    prmAddParam(msr->prm,"bPeriodic",0,&msr->param.bPeriodic,sizeof(int),"p",
		"periodic/non-periodic = -p");
    msr->param.bParaRead = 1;
    prmAddParam(msr->prm,"bParaRead",0,&msr->param.bParaRead,sizeof(int),"par",
		"enable/disable parallel reading of files = +par");
    msr->param.bParaWrite = 1;
    prmAddParam(msr->prm,"bParaWrite",0,&msr->param.bParaWrite,sizeof(int),"paw",
		"enable/disable parallel writing of files = +paw");
    msr->param.bCannonical = 1;
    prmAddParam(msr->prm,"bCannonical",0,&msr->param.bCannonical,sizeof(int),"can",
		"enable/disable use of cannonical momentum = +can");
    msr->param.bDoDensity = 1;
    prmAddParam(msr->prm,"bDoDensity",0,&msr->param.bDoDensity,sizeof(int),
		"den","enable/disable density outputs = +den");
#ifdef USE_PNG
    msr->param.nPNGResolution = 0;
    prmAddParam(msr->prm,"nPNGResolution",0,&msr->param.nPNGResolution,sizeof(int),
		"png","PNG output resolution (zero disables) = 0");
#endif
    msr->param.bDodtOutput = 0;
    prmAddParam(msr->prm,"bDodtOutput",0,&msr->param.bDodtOutput,sizeof(int),
		"dtout","enable/disable dt outputs = -dtout");
    msr->param.nBucket = 8;
    prmAddParam(msr->prm,"nBucket",1,&msr->param.nBucket,sizeof(int),"b",
		"<max number of particles in a bucket> = 8");
    msr->param.dFracNoTreeSqueeze = 0.01;
    prmAddParam(msr->prm,"dFracNoTreeSqueeze",2,&msr->param.dFracNoTreeSqueeze,
		sizeof(double),"fnts",
		"<Fraction of Active Particles for no Tree Squeeze> = 0.01");
    msr->param.iStartStep = 0;
    prmAddParam(msr->prm,"iStartStep",1,&msr->param.iStartStep,
		sizeof(int),"nstart","<initial step numbering> = 0");
    msr->param.nSteps = 0;
    prmAddParam(msr->prm,"nSteps",1,&msr->param.nSteps,sizeof(int),"n",
		"<number of timesteps> = 0");
    msr->param.iOutInterval = 0;
    prmAddParam(msr->prm,"iOutInterval",1,&msr->param.iOutInterval,sizeof(int),
		"oi","<number of timesteps between snapshots> = 0");
    strcpy(msr->param.achOutTypes,"rvmsp");
    prmAddParam(msr->prm,"achOutTypes",3,msr->param.achOutTypes,256,"ot",
		"<output types for snapshort> = \"rvmso\"");
    msr->param.iCheckInterval = 0;
    prmAddParam(msr->prm,"iCheckInterval",1,&msr->param.iCheckInterval,sizeof(int),
		"oc","<number of timesteps between checkpoints> = 0");
    strcpy(msr->param.achCheckTypes,"RVMSP");
    prmAddParam(msr->prm,"achCheckTypes",3,msr->param.achCheckTypes,256,"ct",
		"<output types for checkpoints> = \"RVMSO\"");
    msr->param.iLogInterval = 10;
    prmAddParam(msr->prm,"iLogInterval",1,&msr->param.iLogInterval,sizeof(int),
		"ol","<number of timesteps between logfile outputs> = 10");
    msr->param.bEwald = 1;
    prmAddParam(msr->prm,"bEwald",0,&msr->param.bEwald,sizeof(int),"ewald",
		"enable/disable Ewald correction = +ewald");
    msr->param.iEwOrder = 4;
    prmAddParam(msr->prm,"iEwOrder",1,&msr->param.iEwOrder,sizeof(int),"ewo",
		"<Ewald multipole expansion order: 1, 2, 3 or 4> = 4");
    msr->param.nReplicas = 0;
    prmAddParam(msr->prm,"nReplicas",1,&msr->param.nReplicas,sizeof(int),"nrep",
		"<nReplicas> = 0 for -p, or 1 for +p");
    msr->param.dSoft = 0.0;
    prmAddParam(msr->prm,"dSoft",2,&msr->param.dSoft,sizeof(double),"e",
		"<gravitational softening length> = 0.0");
    msr->param.dSoftMax = 0.0;
    prmAddParam(msr->prm,"dSoftMax",2,&msr->param.dSoftMax,sizeof(double),"eMax",
		"<maximum comoving gravitational softening length (abs or multiplier)> = 0.0");
    msr->param.bPhysicalSoft = 0;
    prmAddParam(msr->prm,"bPhysicalSoft",0,&msr->param.bPhysicalSoft,sizeof(int),"PhysSoft",
		"<Physical gravitational softening length> -PhysSoft");
    msr->param.bSoftMaxMul = 1;
    prmAddParam(msr->prm,"bSoftMaxMul",0,&msr->param.bSoftMaxMul,sizeof(int),"SMM",
		"<Use maximum comoving gravitational softening length as a multiplier> +SMM");
    msr->param.bVariableSoft = 0;
    prmAddParam(msr->prm,"bVariableSoft",0,&msr->param.bVariableSoft,sizeof(int),"VarSoft",
		"<Variable gravitational softening length> -VarSoft");
    msr->param.nSoftNbr = 32;
    prmAddParam(msr->prm,"nSoftNbr",1,&msr->param.nSoftNbr,sizeof(int),"VarSoft",
		"<Neighbours for Variable gravitational softening length> 32");
    msr->param.bSoftByType = 1;
    prmAddParam(msr->prm,"bSoftByType",0,&msr->param.bSoftByType,sizeof(int),"SBT",
		"<Variable gravitational softening length by Type> +SBT");
    msr->param.bDoSoftOutput = 0;
    prmAddParam(msr->prm,"bDoSoftOutput",0,&msr->param.bDoSoftOutput,sizeof(int),
		"softout","enable/disable soft outputs = -softout");
    msr->param.bDoRungOutput = 0;
    prmAddParam(msr->prm,"bDoRungOutput",0,&msr->param.bDoRungOutput,sizeof(int),
		"rungout","enable/disable rung outputs = -rungout");
    msr->param.dDelta = 0.0;
    prmAddParam(msr->prm,"dDelta",2,&msr->param.dDelta,sizeof(double),"dt",
		"<time step>");
    msr->param.dEta = 0.1;
    prmAddParam(msr->prm,"dEta",2,&msr->param.dEta,sizeof(double),"eta",
		"<time step criterion> = 0.1");
    msr->param.bGravStep = 0;
    prmAddParam(msr->prm,"bGravStep",0,&msr->param.bGravStep,sizeof(int),
		"gs","<Gravity timestepping according to iTimeStep Criterion>");
    msr->param.bEpsAccStep = 0;
    prmAddParam(msr->prm,"bEpsAccStep",0,&msr->param.bEpsAccStep,sizeof(int),
		"ea", "<Sqrt(Epsilon on a) timestepping>");
    msr->param.bSqrtPhiStep = 0;
    prmAddParam(msr->prm,"bSqrtPhiStep",0,&msr->param.bSqrtPhiStep,sizeof(int),
		"sphi", "<Sqrt(Phi) on a timestepping>");
    msr->param.bDensityStep = 0;
    prmAddParam(msr->prm,"bDensityStep",0,&msr->param.bDensityStep,sizeof(int),
		"isrho", "<Sqrt(1/Rho) timestepping>");
    msr->param.iTimeStepCrit = 0;
    prmAddParam(msr->prm,"iTimeStepCrit",1,&msr->param.iTimeStepCrit,sizeof(int),
		"tsc", "<Criteria for dynamical time-stepping>");
    msr->param.nPartRhoLoc = 32;
    prmAddParam(msr->prm,"nPartRhoLoc",1,&msr->param.nPartRhoLoc,sizeof(int),
		"nprholoc", "<Number of particles for local density in dynamical time-stepping>");
    msr->param.dPreFacRhoLoc = 4.0*M_PI/3.0;
    prmAddParam(msr->prm,"dPreFacRhoLoc",2,&msr->param.dPreFacRhoLoc,sizeof(double),
		"nprholoc", "<Pre-factor for local density in dynamical time-stepping>");
    msr->param.nPartColl = 0;
    prmAddParam(msr->prm,"nPartColl",1,&msr->param.nPartColl,sizeof(int),
		"npcoll", "<Number of particles in collisional regime>");
    msr->param.nTruncateRung = 0;
    prmAddParam(msr->prm,"nTruncateRung",1,&msr->param.nTruncateRung,sizeof(int),"nTR",
		"<number of MaxRung particles to delete MaxRung> = 0");
    msr->param.iMaxRung = 16;
    prmAddParam(msr->prm,"iMaxRung",1,&msr->param.iMaxRung,sizeof(int),
		"mrung", "<maximum timestep rung>");
    msr->param.nRungVeryActive = 31;
    prmAddParam(msr->prm,"nRungVeryActive",1,&msr->param.nRungVeryActive,
		sizeof(int), "nvactrung", "<timestep rung to use very active timestepping>");
    msr->param.nPartVeryActive = 0;
    prmAddParam(msr->prm,"nPartVeryActive",1,&msr->param.nPartVeryActive,
		sizeof(int), "nvactpart", "<number of particles to use very active timestepping>");
    msr->param.dEwCut = 2.6;
    prmAddParam(msr->prm,"dEwCut",2,&msr->param.dEwCut,sizeof(double),"ew",
		"<dEwCut> = 2.6");
    msr->param.dEwhCut = 2.8;
    prmAddParam(msr->prm,"dEwhCut",2,&msr->param.dEwhCut,sizeof(double),"ewh",
		"<dEwhCut> = 2.8");
    msr->param.dTheta = 0.8;
    msr->param.dTheta2 = msr->param.dTheta;
    prmAddParam(msr->prm,"dTheta",2,&msr->param.dTheta,sizeof(double),"theta",
		"<Barnes opening criterion> = 0.8");
    prmAddParam(msr->prm,"dTheta2",2,&msr->param.dTheta2,sizeof(double),
		"theta2","<Barnes opening criterion for a >= daSwitchTheta> = 0.8");
    msr->param.daSwitchTheta = 1./3.;
    prmAddParam(msr->prm,"daSwitchTheta",2,&msr->param.daSwitchTheta,sizeof(double),"aSwitchTheta",
		"<a to switch theta at> = 1./3.");
    msr->param.dPeriod = 1.0;
    prmAddParam(msr->prm,"dPeriod",2,&msr->param.dPeriod,sizeof(double),"L",
		"<periodic box length> = 1.0");
    msr->param.dxPeriod = 1.0;
    prmAddParam(msr->prm,"dxPeriod",2,&msr->param.dxPeriod,sizeof(double),"Lx",
		"<periodic box length in x-dimension> = 1.0");
    msr->param.dyPeriod = 1.0;
    prmAddParam(msr->prm,"dyPeriod",2,&msr->param.dyPeriod,sizeof(double),"Ly",
		"<periodic box length in y-dimension> = 1.0");
    msr->param.dzPeriod = 1.0;
    prmAddParam(msr->prm,"dzPeriod",2,&msr->param.dzPeriod,sizeof(double),"Lz",
		"<periodic box length in z-dimension> = 1.0");
    msr->param.achInFile[0] = 0;
    prmAddParam(msr->prm,"achInFile",3,msr->param.achInFile,256,"I",
		"<input file name> (file in TIPSY binary format)");
    strcpy(msr->param.achOutName,"pkdgrav");
    prmAddParam(msr->prm,"achOutName",3,msr->param.achOutName,256,"o",
		"<output name for snapshots and logfile> = \"pkdgrav\"");
    strcpy(msr->param.achOutPath,"");
    prmAddParam(msr->prm,"achOutPath",3,msr->param.achOutPath,256,"op",
		"<output path for snapshots and logfile> = \"\"");
    strcpy(msr->param.achIoPath,"");
    prmAddParam(msr->prm,"achIoPath",3,msr->param.achIoPath,256,"iop",
		"<output path for snapshots and logfile> = \"\"");
    msr->param.csm->bComove = 0;
    prmAddParam(msr->prm,"bComove",0,&msr->param.csm->bComove,sizeof(int),
		"cm", "enable/disable comoving coordinates = -cm");
    msr->param.csm->dHubble0 = 0.0;
    prmAddParam(msr->prm,"dHubble0",2,&msr->param.csm->dHubble0, 
		sizeof(double),"Hub", "<dHubble0> = 0.0");
    msr->param.csm->dOmega0 = 1.0;
    prmAddParam(msr->prm,"dOmega0",2,&msr->param.csm->dOmega0,
		sizeof(double),"Om", "<dOmega0> = 1.0");
    msr->param.csm->dLambda = 0.0;
    prmAddParam(msr->prm,"dLambda",2,&msr->param.csm->dLambda,
		sizeof(double),"Lambda", "<dLambda> = 0.0");
    msr->param.csm->dOmegaDE = 0.0;
    prmAddParam(msr->prm,"dOmegaDE",2,&msr->param.csm->dOmegaDE,
		sizeof(double),"OmDE", "Omega for Dark Energy using w0 and wa parameters: <dOmegaDE> = 0.0");
    msr->param.csm->w0 = -1.0;
    prmAddParam(msr->prm,"w0",2,&msr->param.csm->w0,
		sizeof(double),"w0", "w0 parameter for Dark Energy <w0> = -1.0 (pure Lambda)");
    msr->param.csm->wa = 0.0;
    prmAddParam(msr->prm,"wa",2,&msr->param.csm->wa,
		sizeof(double),"wa", "wa parameter for Dark Energy <wa> = 0.0 (pure Lambda)");
    msr->param.csm->dOmegaRad = 0.0;
    prmAddParam(msr->prm,"dOmegaRad",2,&msr->param.csm->dOmegaRad,
		sizeof(double),"Omrad", "<dOmegaRad> = 0.0");
    msr->param.csm->dOmegab = 0.0;
    prmAddParam(msr->prm,"dOmegab",2,&msr->param.csm->dOmegab,
		sizeof(double),"Omb", "<dOmegab> = 0.0");
    strcpy(msr->param.achDataSubPath,".");
    prmAddParam(msr->prm,"achDataSubPath",3,msr->param.achDataSubPath,256,
		NULL,NULL);
    msr->param.dExtraStore = 0.1;
    prmAddParam(msr->prm,"dExtraStore",2,&msr->param.dExtraStore,
		sizeof(double),NULL,NULL);
    msr->param.nSmooth = 64;
    prmAddParam(msr->prm,"nSmooth",1,&msr->param.nSmooth,sizeof(int),"s",
		"<number of particles to smooth over> = 64");
    msr->param.bStandard = 0;
    prmAddParam(msr->prm,"bStandard",0,&msr->param.bStandard,sizeof(int),"std",
		"output in standard TIPSY binary format = -std");
    msr->param.bHDF5 = 0;
    prmAddParam(msr->prm,"bHDF5",0,&msr->param.bHDF5,sizeof(int),"hdf5",
		"output in HDF5 format = -hdf5");
    msr->param.bDoublePos = 0;
    prmAddParam(msr->prm,"bDoublePos",0,&msr->param.bDoublePos,sizeof(int),"dp",
		"input/output double precision positions (standard format only) = -dp");
    msr->param.dRedTo = 0.0;	
    prmAddParam(msr->prm,"dRedTo",2,&msr->param.dRedTo,sizeof(double),"zto",
		"specifies final redshift for the simulation");
    msr->param.nGrowMass = 0;
    prmAddParam(msr->prm,"nGrowMass",1,&msr->param.nGrowMass,sizeof(int),
		"gmn","<number of particles to increase mass> = 0");
    msr->param.dGrowDeltaM = 0.0;
    prmAddParam(msr->prm,"dGrowDeltaM",2,&msr->param.dGrowDeltaM,
		sizeof(double),"gmdm","<Total growth in mass/particle> = 0.0");
    msr->param.dGrowStartT = 0.0;
    prmAddParam(msr->prm,"dGrowStartT",2,&msr->param.dGrowStartT,
		sizeof(double),"gmst","<Start time for growing mass> = 0.0");
    msr->param.dGrowEndT = 1.0;
    prmAddParam(msr->prm,"dGrowEndT",2,&msr->param.dGrowEndT,
		sizeof(double),"gmet","<End time for growing mass> = 1.0");
    msr->param.dFracNoDomainDecomp = 0.001;
    prmAddParam(msr->prm,"dFracNoDomainDecomp",2,&msr->param.dFracNoDomainDecomp,
		sizeof(double),"fndd",
		"<Fraction of Active Particles for no DD> = 0.001");
    msr->param.dFracNoDomainRootFind = 0.02;
    prmAddParam(msr->prm,"dFracNoDomainRootFind",2,&msr->param.dFracNoDomainRootFind,
		sizeof(double),"fndrf",
		"<Fraction of Active Particles for no DD root finding> = 0.02");
    msr->param.dFracNoDomainDimChoice = 0.1;
    prmAddParam(msr->prm,"dFracNoDomainDimChoice",2,&msr->param.dFracNoDomainDimChoice,
		sizeof(double),"fnddc",
		"<Fraction of Active Particles for no DD dimension choice> = 0.1");
    msr->param.bDoGravity = 1;
    prmAddParam(msr->prm,"bDoGravity",0,&msr->param.bDoGravity,sizeof(int),"g",
		"enable/disable interparticle gravity = +g");
#ifdef HERMITE
    msr->param.bHermite = 0;
    prmAddParam(msr->prm,"bHermite",0,&msr->param.bHermite,
		sizeof(int),"hrm","<Hermite integratot>");
    msr->param.bAarsethStep = 0;
    prmAddParam(msr->prm,"bAarsethStep",0,&msr->param.bAarsethStep,sizeof(int),
				"aas","<Aarseth timestepping>");
#endif
    msr->param.iWallRunTime = 0;
    prmAddParam(msr->prm,"iWallRunTime",1,&msr->param.iWallRunTime,
		sizeof(int),"wall",
		"<Maximum Wallclock time (in minutes) to run> = 0 = infinite");
    msr->param.bAntiGrav = 0;
    prmAddParam(msr->prm,"bAntiGrav",0,&msr->param.bAntiGrav,sizeof(int),"antigrav",
		"reverse gravity making it repulsive = -antigrav");
    msr->param.nFindGroups = 0;
    prmAddParam(msr->prm,"nFindGroups",1,&msr->param.nFindGroups,sizeof(int),
		"nFindGroups","<number of iterative FOFs> = 0");
    msr->param.nMinMembers = 16;
    prmAddParam(msr->prm,"nMinMembers",1,&msr->param.nMinMembers,sizeof(int),
		"nMinMembers","<minimum number of group members> = 16");
    msr->param.dTau = 0.164;      
    prmAddParam(msr->prm,"dTau",2,&msr->param.dTau,sizeof(double),"dTau",
                "<linking lenght for first FOF in units of mean particle separation> = 0.164");
    msr->param.dVTau = 1.0e10;
    prmAddParam(msr->prm,"dVTau",2,&msr->param.dVTau,sizeof(double),"dVTau",
		"<velocity space linking lenght for phase-space FOF, set to HUGE for plain FOF> = 1.0e10");
    msr->param.bTauAbs = 0;
    prmAddParam(msr->prm,"bTauAbs",0,&msr->param.bTauAbs,sizeof(int),"bTauAbs",
		"<if 1 use simulation units for dTau, not mean particle separation> = 0");
    msr->param.nBins = 0;
    prmAddParam(msr->prm,"nBins",1,&msr->param.nBins,sizeof(int),"nBins",
		"<number of bin in profiles, no profiles if 0 or negative> = 0");
    msr->param.bUsePotmin = 0;
    prmAddParam(msr->prm,"bUsePotmin",0,&msr->param.bUsePotmin,sizeof(int),"bUsePotmin",
		"<if 1 use minima of pot. instead of com for fist FOF level profiles> = 0");
    msr->param.nMinProfile = 2000;
    prmAddParam(msr->prm,"nMinProfile",1,&msr->param.nMinProfile,sizeof(int),"nMinProfile",
		"<minimum number of particles in a group to make a profile> = 2000");
    msr->param.fBinsRescale = 1.0;
    prmAddParam(msr->prm,"fBinsRescale",2,&msr->param.fBinsRescale,sizeof(double),"fBinsRescale",
		"<if bigger than one, nMinProfile gets smaller in recursive fofs> = 1.0");
    msr->param.fContrast = 10.0;
    prmAddParam(msr->prm,"fContrast",2,&msr->param.fContrast,sizeof(double),"fContrast",
		"<the density contrast used for recursive fofs> = 10.0");
    msr->param.Delta = 200.0;
    prmAddParam(msr->prm,"Delta",2,&msr->param.Delta,sizeof(double),"Delta",
		"<the density contrast whitin the virial radius over simulation unit density> = 200.0");
    msr->param.binFactor = 3.0;
    prmAddParam(msr->prm,"binFactor",2,&msr->param.binFactor,sizeof(double),"binFactor",
		"<ratio of largest spherical bin to fof determined group radius> = 3.0");
    msr->param.bLogBins = 0;
    prmAddParam(msr->prm,"bLogBins",0,&msr->param.bLogBins,
		sizeof(int),"bLogBins","use logaritmic bins instead of linear = -bLogBins");
#ifdef RELAXATION
    msr->param.bTraceRelaxation = 0;
    prmAddParam(msr->prm,"bTraceRelaxation",0,&msr->param.bTraceRelaxation,sizeof(int),
		"rtrace","<enable/disable relaxation tracing> = -rtrace");
#endif /* RELAXATION */

#ifdef PLANETS 
    msr->param.bCollision = 0;
	prmAddParam(msr->prm,"bCollision",0,&msr->param.bCollision,
		    sizeof(int),"hc","<Collisions>");
    msr->param.bHeliocentric = 0;
	prmAddParam(msr->prm,"bHeliocentric",0,&msr->param.bHeliocentric,
		    sizeof(int),"hc","use/don't use Heliocentric coordinates = -hc");  
    msr->param.dCentMass = 1.0;
    prmAddParam(msr->prm,"dCentMass",2,&msr->param.dCentMass,sizeof(double),
		"fgm","specifies the central mass for Keplerian orbits");

#ifdef SYMBA
    msr->param.bSymba = 1;
    prmAddParam(msr->prm,"bSymba",0,&msr->param.bSymba,sizeof(int),
		"sym","use Symba integrator");
#endif 

/* collision stuff */
    msr->param.iCollLogOption = 0;
    prmAddParam(msr->prm,"iCollLogOption",1,&msr->param.iCollLogOption,
		sizeof(int),"clog","<Collision log option> = 0");	
    msr->param.CP.iOutcomes = BOUNCE;
    prmAddParam(msr->prm,"iOutcomes",1,&msr->param.CP.iOutcomes,
		sizeof(int),"outcomes","<Allowed collision outcomes> = 0");    
    msr->param.CP.dBounceLimit = 1.0;
    prmAddParam(msr->prm,"dBounceLimit",2,&msr->param.CP.dBounceLimit,
		sizeof(double),"blim","<Bounce limit> = 1.0");
    msr->param.CP.iBounceOption = ConstEps;
    prmAddParam(msr->prm,"iBounceOption",1,&msr->param.CP.iBounceOption,
		sizeof(int),"bopt","<Bounce option> = 0");
    msr->param.CP.dEpsN = 1.0;
    prmAddParam(msr->prm,"dEpsN",2,&msr->param.CP.dEpsN,
		sizeof(double),"epsn","<Coefficient of restitution> = 1");
    msr->param.CP.dEpsT = 1.0;
    prmAddParam(msr->prm,"dEpsT",2,&msr->param.CP.dEpsT,
		sizeof(double),"epst","<Coefficient of surface friction> = 1");
    msr->param.CP.dDensity = 0.0;
    prmAddParam(msr->prm,"dDensity",2,&msr->param.CP.dDensity,
		sizeof(double),"density","<Merged particle density> = 0");
    msr->param.CP.bFixCollapse = 0;
    prmAddParam(msr->prm,"bFixCollapse",0,&msr->param.CP.bFixCollapse,
		sizeof(int),"overlap","enable/disable overlap fix = -overlap");
#endif /* PLANETS */

    /* IC Generation */
#ifdef USE_GRAFIC
    msr->param.h = 0.0;
    prmAddParam(msr->prm,"h",2,&msr->param.h,
		sizeof(double),"h","<hubble parameter h> = 0");
    msr->param.dBoxSize = 0.0;
    prmAddParam(msr->prm,"dBoxSize",2,&msr->param.dBoxSize,
		sizeof(double),"mpc","<Simulation Box size in Mpc> = 0");
    msr->param.nGrid = 0;
    prmAddParam(msr->prm,"nGrid",1,&msr->param.nGrid,
		sizeof(int),"grid","<Grid size for IC 0=disabled> = 0");	
    msr->param.iSeed = 0;
    prmAddParam(msr->prm,"iSeed",1,&msr->param.iSeed,
		sizeof(int),"seed","<Random seed for IC> = 0");	
    msr->param.bWriteIC = 0;
    prmAddParam(msr->prm,"bWriteIC",1,&msr->param.bWriteIC,
		sizeof(int),"wic","<Write IC after generating> = 0");	
#endif


/*
** Set the box center to (0,0,0) for now!
*/
    for (j=0;j<3;++j) msr->fCenter[j] = 0.0;
    /*
    ** Define any "LOCAL" parameters (LCL)
    */
    msr->lcl.pszDataPath = getenv("PTOOLS_DATA_PATH");
    /*
    ** Process command line arguments. 
    ** (actual parameters read from the param file are set here)
    */ 
    ret = prmArgProc(msr->prm,argc,argv);
    if (!ret) {
	_msrExit(msr,1);
	}

    if (nDigits < 1 || nDigits > 9) {
	(void) fprintf(stderr,"Unreasonable number of filename digits.\n");
	_msrExit(msr,1);
	}

    (void) sprintf(msr->param.achDigitMask,"%%s.%%0%ii",nDigits);
    /*
    ** Make sure that we have some setting for nReplicas if bPeriodic is set.
    */
    if (msr->param.bPeriodic && !prmSpecified(msr->prm,"nReplicas")) {
	msr->param.nReplicas = 1;
	}
    /*
    ** Warn that we have a setting for nReplicas if bPeriodic NOT set.
    */
    if (!msr->param.bPeriodic && msr->param.nReplicas != 0) {
	printf("WARNING: nReplicas set to non-zero value for non-periodic!\n");
	}

    if ( msr->param.nGrid ) {
	if (msr->param.achInFile[0]) {
	    puts("ERROR: do not specify an input file when generating IC");
	    _msrExit(msr,1);
	}

	if ( msr->param.iSeed == 0 ) {
	    puts("ERROR: Random seed for IC not specified");
	    _msrExit(msr,1);
	}
	if ( msr->param.dBoxSize <= 0 ) {
	    puts("ERROR: Box size for IC not specified");
	    _msrExit(msr,1);
	}
	if ( msr->param.h <= 0 ) {
	    puts("ERROR: Hubble parameter (h) was not specified for IC generation");
	    _msrExit(msr,1);
	}
    }
    else if (!msr->param.achInFile[0]) {
	puts("ERROR: no input file specified");
	_msrExit(msr,1);
	}

    if (msr->param.dTheta <= 0) {
	if (msr->param.dTheta == 0 && msr->param.bVWarnings)
	    fprintf(stderr,"WARNING: Zero opening angle may cause numerical problems\n");
	else if (msr->param.dTheta < 0) {
	    fprintf(stderr,"ERROR: Opening angle must be non-negative\n");
	    _msrExit(msr,1);
	    }
	}

    if ( msr->param.dFracNoDomainDecomp > msr->param.dFracNoDomainRootFind
	 || msr->param.dFracNoDomainRootFind > msr->param.dFracNoDomainDimChoice
	 || msr->param.dFracNoDomainDecomp<0.0 || msr->param.dFracNoDomainDimChoice > 1.0 ) {
	puts("ERROR: check that 0 <= dFracNoDomainDecomp <= dFracNoDomainRootFind <= dFracNoDomainDimChoice <= 1");
	_msrExit(msr,1);
	}
   
    msr->nThreads = mdlThreads(mdl);
	
    /*
    ** Always set bCannonical = 1 if bComove == 0
    */
    if (!msr->param.csm->bComove) {
	if (!msr->param.bCannonical)
	    printf("WARNING: bCannonical reset to 1 for non-comoving (bComove == 0)\n");
	msr->param.bCannonical = 1;
	} 
    /* 
     * Softening 
     */
	
    if (msr->param.bPhysicalSoft || msr->param.bVariableSoft) {
	if (msr->param.bPhysicalSoft && !msrComove(msr)) {
	    printf("WARNING: bPhysicalSoft reset to 0 for non-comoving (bComove == 0)\n");
	    msr->param.bPhysicalSoft = 0;
	    }
#ifndef CHANGESOFT
	fprintf(stderr,"ERROR: You must compile with -DCHANGESOFT to use changing softening options\n");
	_msrExit(msr,1);
#endif
	if (msr->param.bVariableSoft && !prmSpecified(msr->prm,"bDoSoftOutput")) msr->param.bDoSoftOutput=1;
  
	if (msr->param.bPhysicalSoft && msr->param.bVariableSoft) {
	    fprintf(stderr,"ERROR: You may only choose one of Physical or Variable softening\n");
	    _msrExit(msr,1);
	    }
	}
    /*
    ** Determine the period of the box that we are using.
    ** Set the new d[xyz]Period parameters which are now used instead
    ** of a single dPeriod, but we still want to have compatibility
    ** with the old method of setting dPeriod.
    */
    if (prmSpecified(msr->prm,"dPeriod") && 
	!prmSpecified(msr->prm,"dxPeriod")) {
	msr->param.dxPeriod = msr->param.dPeriod;
	}
    if (prmSpecified(msr->prm,"dPeriod") && 
	!prmSpecified(msr->prm,"dyPeriod")) {
	msr->param.dyPeriod = msr->param.dPeriod;
	}
    if (prmSpecified(msr->prm,"dPeriod") && 
	!prmSpecified(msr->prm,"dzPeriod")) {
	msr->param.dzPeriod = msr->param.dPeriod;
	}
    /*
    ** Periodic boundary conditions can be disabled along any of the
    ** x,y,z axes by specifying a period of zero for the given axis.
    ** Internally, the period is set to infinity (Cf. pkdBucketWalk()
    ** and pkdDrift(); also the INTERSECT() macro in smooth.h).
    */
    if (msr->param.dPeriod  == 0) msr->param.dPeriod  = FLOAT_MAXVAL;
    if (msr->param.dxPeriod == 0) msr->param.dxPeriod = FLOAT_MAXVAL;
    if (msr->param.dyPeriod == 0) msr->param.dyPeriod = FLOAT_MAXVAL;
    if (msr->param.dzPeriod == 0) msr->param.dzPeriod = FLOAT_MAXVAL;
    /*
    ** Determine opening type.
    */

    msr->dCrit = msr->param.dTheta;
    if (!prmSpecified(msr->prm,"dTheta2")) 
	msr->param.dTheta2 = msr->param.dTheta;
    /*
    ** Initialize comove variables.
    */
    msr->nMaxOuts = 100;
    msr->pdOutTime = malloc(msr->nMaxOuts*sizeof(double));
    assert(msr->pdOutTime != NULL);
    msr->nOuts = 0;

    /*
    ** Check timestepping.
    */

    if (msr->param.iMaxRung < 1) {
	msr->param.iMaxRung = 1;
	if (msr->param.bVWarnings)
	    (void) fprintf(stderr,"WARNING: iMaxRung set to 1\n");
	}

    if (msr->param.bGravStep && !msr->param.bDoGravity) {
	puts("ERROR: need gravity to use gravity stepping...");
	_msrExit(msr,1);
	}
    if (msr->param.bEpsAccStep || msr->param.bSqrtPhiStep) {
	msr->param.bAccelStep = 1;
	}
    else {
	msr->param.bAccelStep = 0;
    }
    if (msr->param.bGravStep) {
	msr->param.bEpsAccStep = 0;   /* we must do this because the meaning of Eta is different */
	}
   
#ifdef PLANETS
	switch (msr->param.iCollLogOption) {
	case COLL_LOG_NONE:
		break;
	case COLL_LOG_VERBOSE:
		(void) strcpy(msr->param.achCollLog,COLL_LOG_TXT);
		break;
	case COLL_LOG_TERSE:
		(void) strcpy(msr->param.achCollLog,COLL_LOG_BIN);
		break;
	default:
		puts("ERROR: Invalid collision log option");
		_msrExit(msr,1);
		}
	if (msr->param.iCollLogOption && msr->param.bVStart)
		printf("Collision log: \"%s\"\n",msr->param.achCollLog);

	COLLISION_PARAMS *CP = &msr->param.CP;
		if (!(CP->iOutcomes & (MERGE | BOUNCE | FRAG))) {
			puts("ERROR: must specify one of MERGE/BOUNCE/FRAG");
			_msrExit(msr,1);
			}

		msr->dEcoll = 0.0;
#ifdef SYMBA
		msr->param.bHeliocentric = 0;
#endif 	
#endif /* PLANETS */


    pstInitialize(&msr->pst,msr->mdl,&msr->lcl);
  
    pstAddServices(msr->pst,msr->mdl);
    /*
    ** Create the processor subset tree.
    */
    for (id=1;id<msr->nThreads;++id) {
	if (msr->param.bVDetails) printf("Adding %d to the pst\n",id);
	inAdd.id = id;
	pstSetAdd(msr->pst,&inAdd,sizeof(inAdd),NULL,NULL);
	}
    if (msr->param.bVDetails) printf("\n");
    /*
    ** Levelize the PST.
    */
    inLvl.iLvl = 0;
    pstLevelize(msr->pst,&inLvl,sizeof(inLvl),NULL,NULL);
    /*
    ** Create the processor mapping array for the one-node output
    ** routines.
    */
    msr->pMap = malloc(msr->nThreads*sizeof(int));
    assert(msr->pMap != NULL);
    inGM.nStart = 0;
    pstGetMap(msr->pst,&inGM,sizeof(inGM),msr->pMap,NULL);
    msr->iCurrMaxRung = 0;
    /*
    ** Mark the Domain Decompositon as not done
    */
    msr->iLastRungRT = -1;
    msr->nRung = malloc((msr->param.iMaxRung+1)*sizeof(uint64_t));
    assert(msr->nRung != NULL);
    for (i=0;i<=msr->param.iMaxRung;++i) msr->nRung[i] = 0;

    msr->iRungVeryActive = msr->param.iMaxRung; /* No very active particles */
    msr->bSavePending = 0;                      /* There is no pending save */
    }


void msrLogParams(MSR msr,FILE *fp)
    {
    double z;
    int i;

#ifdef __DATE__
#ifdef __TIME__
    fprintf(fp,"# Compiled: %s %s\n",__DATE__,__TIME__);
#endif
#endif
    fprintf(fp,"# Preprocessor macros:");
#ifdef CHANGESOFT
    fprintf(fp," CHANGESOFT");
#endif
#ifdef DEBUG
    fprintf(fp," DEBUG");
#endif
#ifdef RELAXATION
    fprintf(fp,"RELAXATION");
#endif
#ifdef _REENTRANT
    fprintf(fp," _REENTRANT");
#endif
#if defined(MAXHOSTNAMELEN) && defined(HAVE_GETHOSTNAME)
	{
	char hostname[MAXHOSTNAMELEN];
	fprintf(fp,"\n# Master host: ");
	if (gethostname(hostname,MAXHOSTNAMELEN))
	    fprintf(fp,"unknown");
	else
	    fprintf(fp,"%s",hostname);
	}
#endif
	fprintf(fp,"\n# N: %"PRIu64,msr->N);
	fprintf(fp," nThreads: %d",msr->param.nThreads);
	fprintf(fp," bDiag: %d",msr->param.bDiag);
	fprintf(fp," Verbosity flags: (%d,%d,%d,%d,%d)",msr->param.bVWarnings,
		msr->param.bVStart,msr->param.bVStep,msr->param.bVRungStat,
		msr->param.bVDetails);
	fprintf(fp,"\n# bPeriodic: %d",msr->param.bPeriodic);
	fprintf(fp," bComove: %d",msr->param.csm->bComove);
	fprintf(fp,"\n# bParaRead: %d",msr->param.bParaRead);
	fprintf(fp," bParaWrite: %d",msr->param.bParaWrite);
	fprintf(fp," bCannonical: %d",msr->param.bCannonical);
	fprintf(fp," bStandard: %d",msr->param.bStandard);
	fprintf(fp," bHDF5: %d",msr->param.bHDF5);
	fprintf(fp," nBucket: %d",msr->param.nBucket);
	fprintf(fp,"\n# dFracNoTreeSqueeze: %g",msr->param.dFracNoTreeSqueeze);
	fprintf(fp," iOutInterval: %d",msr->param.iOutInterval);
	fprintf(fp," iCheckInterval: %d",msr->param.iCheckInterval);
	fprintf(fp," iLogInterval: %d",msr->param.iLogInterval);
	fprintf(fp," iEwOrder: %d",msr->param.iEwOrder);
	fprintf(fp," nReplicas: %d",msr->param.nReplicas);
	fprintf(fp,"\n# dEwCut: %f",msr->param.dEwCut);
	fprintf(fp," dEwhCut: %f",msr->param.dEwhCut);
	fprintf(fp,"\n# iStartStep: %d",msr->param.iStartStep);
	fprintf(fp," nSteps: %d",msr->param.nSteps);
	fprintf(fp," nSmooth: %d",msr->param.nSmooth);
	fprintf(fp," dExtraStore: %f",msr->param.dExtraStore);
	if (prmSpecified(msr->prm,"dSoft"))
	    fprintf(fp," dSoft: %g",msr->param.dSoft);
	else
	    fprintf(fp," dSoft: input");
	fprintf(fp,"\n# bPhysicalSoft: %d",msr->param.bPhysicalSoft);
	fprintf(fp," bVariableSoft: %d",msr->param.bVariableSoft);
	fprintf(fp," nSoftNbr: %d",msr->param.nSoftNbr);
	fprintf(fp," bSoftByType: %d",msr->param.bSoftByType);
	fprintf(fp," bSoftMaxMul: %d",msr->param.bSoftMaxMul);
	fprintf(fp," dSoftMax: %g",msr->param.dSoftMax);
	fprintf(fp," bDoSoftOutput: %d",msr->param.bDoSoftOutput);
	fprintf(fp,"\n# dDelta: %g",msr->param.dDelta);
	fprintf(fp," dEta: %g",msr->param.dEta);
	fprintf(fp," iMaxRung: %d",msr->param.iMaxRung);
	fprintf(fp," nRungVeryActive: %d",msr->param.nRungVeryActive);
	fprintf(fp," bDoRungOutput: %d",msr->param.bDoRungOutput);
	fprintf(fp,"\n# bGravStep: %d",msr->param.bGravStep);
	fprintf(fp," bEpsAccStep: %d",msr->param.bEpsAccStep);
	fprintf(fp," bSqrtPhiStep: %d",msr->param.bSqrtPhiStep);
	fprintf(fp," bDensityStep: %d",msr->param.bDensityStep);
	fprintf(fp," nTruncateRung: %d",msr->param.nTruncateRung);
	fprintf(fp,"\n# iTimeStepCrit: %d",msr->param.iTimeStepCrit);
	fprintf(fp," nPartRhoLoc: %d", msr->param.nPartRhoLoc);
	fprintf(fp," dPreFacRhoLoc: %g", msr->param.dPreFacRhoLoc);
	fprintf(fp," nPartColl: %d", msr->param.nPartColl);
	fprintf(fp,"\n# bDoGravity: %d",msr->param.bDoGravity);
#ifdef HERMITE
	fprintf(fp," bHermite: %d",msr->param.bHermite);
	fprintf(fp," bAarsethStep: %d",msr->param.bAarsethStep);
#endif
	fprintf(fp,"\n# dFracNoDomainDecomp: %g",msr->param.dFracNoDomainDecomp);
	fprintf(fp," dFracNoDomainRootFind: %g",msr->param.dFracNoDomainRootFind);
	fprintf(fp," dFracNoDomainDimChoice: %g",msr->param.dFracNoDomainDimChoice);
	fprintf(fp,"\n# nTruncateRung: %d",msr->param.nTruncateRung);
	fprintf(fp,"\n# nGrowMass: %d",msr->param.nGrowMass);
	fprintf(fp," dGrowDeltaM: %g",msr->param.dGrowDeltaM);
	fprintf(fp," dGrowStartT: %g",msr->param.dGrowStartT);
	fprintf(fp," dGrowEndT: %g",msr->param.dGrowEndT);
	fprintf(fp,"\n# Group Find: nFindGroups: %d",msr->param.nFindGroups);
	fprintf(fp," dTau: %g",msr->param.dTau);
	fprintf(fp," dVTau: %g",msr->param.dVTau);
	fprintf(fp," bTauAbs: %d",msr->param.bTauAbs);
	fprintf(fp," nMinMembers: %d",msr->param.nMinMembers);	
	fprintf(fp," nBins: %d",msr->param.nBins);
	fprintf(fp," bUsePotmin: %d",msr->param.bUsePotmin);
	fprintf(fp," nMinProfile: %d",msr->param.nMinProfile);
	fprintf(fp," fBinsRescale: %g",msr->param.fBinsRescale);
	fprintf(fp," fContrast: %g",msr->param.fContrast);
	fprintf(fp," Delta: %g",msr->param.Delta);
	fprintf(fp," binFactor: %g",msr->param.binFactor);
	fprintf(fp," bLogBins: %d",msr->param.bLogBins);
#ifdef RELAXATION
	fprintf(fp,"\n# Relaxation estimate: bTraceRelaxation: %d",msr->param.bTraceRelaxation);	
#endif
#ifdef PLANETS
#ifdef SYMBA
	fprintf(fp," bSymba: %d",msr->param.bSymba);
#endif	
	fprintf(fp," bCollision: %d",msr->param.bCollision);
	fprintf(fp," bHeliocentric: %d",msr->param.bHeliocentric);
	fprintf(fp," dCentMass: %g",msr->param.dCentMass);
	fprintf(fp,"\n# Collisions...");
	fprintf(fp," iCollLogOption: %d",msr->param.iCollLogOption);
	fprintf(fp,"\n# iOutcomes: %d",msr->param.CP.iOutcomes);   
	fprintf(fp," dBounceLimit: %g",msr->param.CP.dBounceLimit);
	fprintf(fp," dEpsN: %g",msr->param.CP.dEpsN);
	fprintf(fp," dEpsT: %g",msr->param.CP.dEpsT);
	fprintf(fp," dDensity: %g",msr->param.CP.dDensity);
	fprintf(fp,"\n# bFixCollapse: %d",msr->param.CP.bFixCollapse);
#endif /* PLANETS */

	fprintf(fp," dTheta: %f",msr->param.dTheta);
	fprintf(fp,"\n# dPeriod: %g",msr->param.dPeriod);
	fprintf(fp," dxPeriod: %g",
		msr->param.dxPeriod >= FLOAT_MAXVAL ? 0 : msr->param.dxPeriod);
	fprintf(fp," dyPeriod: %g",
		msr->param.dyPeriod >= FLOAT_MAXVAL ? 0 : msr->param.dyPeriod);
	fprintf(fp," dzPeriod: %g",
		msr->param.dzPeriod >= FLOAT_MAXVAL ? 0 : msr->param.dzPeriod);
	fprintf(fp,"\n# dHubble0: %g",msr->param.csm->dHubble0);
	fprintf(fp," dOmega0: %g",msr->param.csm->dOmega0);
	fprintf(fp," dLambda: %g",msr->param.csm->dLambda);
	fprintf(fp," dOmegaDE: %g",msr->param.csm->dOmegaDE);
	fprintf(fp," w0: %g",msr->param.csm->w0);
	fprintf(fp," wa: %g",msr->param.csm->wa);
	fprintf(fp," dOmegaRad: %g",msr->param.csm->dOmegaRad);
	fprintf(fp," dOmegab: %g",msr->param.csm->dOmegab);
	fprintf(fp,"\n# achInFile: %s",msr->param.achInFile);
	fprintf(fp,"\n# achOutName: %s",msr->param.achOutName); 
	fprintf(fp,"\n# achOutPath: %s",msr->param.achOutPath); 
	fprintf(fp,"\n# achIoPath: %s",msr->param.achIoPath); 
	fprintf(fp,"\n# achDataSubPath: %s",msr->param.achDataSubPath);
	if (msr->param.csm->bComove) {
	    fprintf(fp,"\n# RedOut:");
	    if (msr->nOuts == 0) fprintf(fp," none");
	    for (i=0;i<msr->nOuts;i++) {
		if (i%5 == 0) fprintf(fp,"\n#   ");
		z = 1.0/csmTime2Exp(msr->param.csm, msr->pdOutTime[i]) - 1.0;
		fprintf(fp," %f",z);
		}
	    fprintf(fp,"\n");
	    }
	else {
	    fprintf(fp,"\n# TimeOut:");
	    if (msr->nOuts == 0) fprintf(fp," none");
	    for (i=0;i<msr->nOuts;i++) {
		if (i%5 == 0) fprintf(fp,"\n#   ");
		fprintf(fp," %f",msr->pdOutTime[i]);
		}
	    fprintf(fp,"\n");
	    }
    }

int
msrGetLock(MSR msr)
    {
    /*
    ** Attempts to lock run directory to prevent overwriting. If an old lock
    ** is detected with the same achOutName, an abort is signaled. Otherwise
    ** a new lock is created. The bOverwrite parameter flag can be used to
    ** suppress lock checking.
    */

    FILE *fp = NULL;
    char achTmp[256],achFile[256];

    _msrMakePath(msr->param.achDataSubPath,LOCKFILE,achTmp);
    _msrMakePath(msr->lcl.pszDataPath,achTmp,achFile);
    if (!msr->param.bOverwrite && (fp = fopen(achFile,"r"))) {
	(void) fscanf(fp,"%s",achTmp);
	(void) fclose(fp);
	if (!strcmp(msr->param.achOutName,achTmp)) {
	    (void) printf("ABORT: %s detected.\nPlease ensure data is safe to "
			  "overwrite. Delete lockfile and try again.\n",achFile);
	    return 0;
	    }
	}
    if (!(fp = fopen(achFile,"w"))) {
	if (msr->param.bOverwrite && msr->param.bVWarnings) {
	    (void) printf("WARNING: Unable to create %s...ignored.\n",achFile);
	    return 1;
	    }
	else {
	    (void) printf("Unable to create %s\n",achFile);
	    return 0;
	    }
	}
    (void) fprintf(fp,msr->param.achOutName);
    (void) fclose(fp);
    return 1;
    }

int
msrCheckForStop(MSR msr)
    {
    /*
    ** Checks for existence of STOPFILE in run directory. If found, the file
    ** is removed and the return status is set to 1, otherwise 0.
    */

    static char achFile[256];
    static int first_call = 1;

    FILE *fp = NULL;

    if (first_call) {
	char achTmp[256];
	_msrMakePath(msr->param.achDataSubPath,STOPFILE,achTmp);
	_msrMakePath(msr->lcl.pszDataPath,achTmp,achFile);
	first_call = 0;
	}
    if ((fp = fopen(achFile,"r"))) {
	(void) printf("User interrupt detected.\n");
	(void) fclose(fp);
	(void) unlink(achFile);
	return 1;
	}
    return 0;
    }

void msrFinish(MSR msr)
    {
#ifndef USE_MDL_IO
    int id;
#endif

    /*
    ** It is possible that a previous save is still pending on the I/O processor.  Wait for it.
    */
#ifdef USE_MDL_IO
    if ( msr->bSavePending ) {
	mdlSetComm(msr->mdl,1);
	mdlGetReply(msr->mdl,0,NULL,NULL);
	mdlSetComm(msr->mdl,0);
    }
#endif

#ifdef USE_MDL_IO
    mdlStop(msr->mdl);
#else
    for (id=1;id<msr->nThreads;++id) {
	if (msr->param.bVDetails) printf("Stopping thread %d\n",id);		
	mdlReqService(msr->mdl,id,SRV_STOP,NULL,0);
	mdlGetReply(msr->mdl,id,NULL,NULL);
	}
#endif
    pstFinish(msr->pst);
    csmFinish(msr->param.csm);
    /*
    ** finish with parameter stuff, deallocate and exit.
    */
    prmFinish(msr->prm);
    free(msr->pMap);
    free(msr->nRung);
    free(msr->pdOutTime);
    free(msr);
    }

#ifdef USE_HDF5
void msrOneNodeReadHDF5(MSR msr, struct inReadTipsy *in)
    {
    int i,id;
    int *nParts;		/* number of particles for each processor */
    uint64_t nStart;
    PST pst0;
    LCL *plcl;
    char achInFile[PST_FILENAME_SIZE];
    char achOutName[PST_FILENAME_SIZE];
    int nid;
    int inswap;
    struct inSetParticleTypes intype;
    hid_t fileID;
    IOHDF5 io;

    nParts = malloc(msr->nThreads*sizeof(*nParts));
    for (id=0;id<msr->nThreads;++id) {
	nParts[id] = -1;
	}

    pstOneNodeReadInit(msr->pst, in, sizeof(*in), nParts, &nid);
    assert((size_t)nid == msr->nThreads*sizeof(*nParts));
    for (id=0;id<msr->nThreads;++id) {
	assert(nParts[id] > 0);
	}

    pst0 = msr->pst;
    while(pst0->nLeaves > 1)
	pst0 = pst0->pstLower;
    plcl = pst0->plcl;
    /*
    ** Add the local Data Path to the provided filename.
    */
    _msrMakePath(plcl->pszDataPath,in->achInFile,achInFile);
    _msrMakePath(plcl->pszDataPath,in->achOutName,achOutName);

    fileID=H5Fopen(achInFile, H5F_ACC_RDONLY, H5P_DEFAULT);
    assert(fileID >= 0);
    io = ioHDF5Initialize( fileID, 32768, IOHDF5_SINGLE );
    assert( io != NULL );

    nStart = nParts[0];
    assert(msr->pMap[0] == 0);
    for (i=1;i<msr->nThreads;++i) {
	id = msr->pMap[i];
	/* 
	 * Read particles into the local storage.
	 */
	assert(plcl->pkd->nStore >= nParts[id]);
	pkdReadHDF5(plcl->pkd, io, in->dvFac, nStart, nParts[id]);
	nStart += nParts[id];
	/* 
	 * Now shove them over to the remote processor.
	 */
	inswap = 0;
	mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
	pkdSwapAll(plcl->pkd, id);
	mdlGetReply(pst0->mdl,id,NULL,NULL);
    	}
    assert(nStart == msr->N);
    /* 
     * Now read our own particles.
     */
    pkdReadHDF5(plcl->pkd, io, in->dvFac, 0, nParts[0]);
    pstSetParticleTypes(msr->pst,&intype,sizeof(intype),NULL,NULL);

    ioHDF5Finish(io);
    H5Fclose(fileID);
    }
#endif

void msrOneNodeReadTipsy(MSR msr, struct inReadTipsy *in)
    {
    int i,id;
    int *nParts;				/* number of particles for each processor */
    uint64_t nStart;
    PST pst0;
    LCL *plcl;
    char achInFile[PST_FILENAME_SIZE];
    char achOutName[PST_FILENAME_SIZE];
    int nid;
    int inswap;
    struct inSetParticleTypes intype;

    nParts = malloc(msr->nThreads*sizeof(*nParts));
    assert(nParts!=NULL);
    for (id=0;id<msr->nThreads;++id) {
	nParts[id] = -1;
	}

    pstOneNodeReadInit(msr->pst, in, sizeof(*in), nParts, &nid);
    assert((size_t)nid == msr->nThreads*sizeof(*nParts));
    for (id=0;id<msr->nThreads;++id) {
	assert(nParts[id] > 0);
	}

    pst0 = msr->pst;
    while(pst0->nLeaves > 1)
	pst0 = pst0->pstLower;
    plcl = pst0->plcl;
    /*
    ** Add the local Data Path to the provided filename.
    */
    _msrMakePath(plcl->pszDataPath,in->achInFile,achInFile);
    _msrMakePath(plcl->pszDataPath,in->achOutName,achOutName);

    nStart = nParts[0];
    assert(msr->pMap[0] == 0);
    for (i=1;i<msr->nThreads;++i) {
	id = msr->pMap[i];
	/* 
	 * Read particles into the local storage.
	 */
	assert(plcl->pkd->nStore >= nParts[id]);
	pkdReadTipsy(plcl->pkd,achInFile,achOutName,nStart,nParts[id],
		     in->bStandard,in->dvFac,in->bDoublePos);
	nStart += nParts[id];
	/* 
	 * Now shove them over to the remote processor.
	 */
	inswap = 0;
	mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
	pkdSwapAll(plcl->pkd, id);
	mdlGetReply(pst0->mdl,id,NULL,NULL);
    	}
    assert(nStart == msr->N);
    /* 
     * Now read our own particles.
     */
    pkdReadTipsy(plcl->pkd,achInFile,achOutName,0,nParts[0],in->bStandard,in->dvFac,
		 in->bDoublePos);
    pstSetParticleTypes(msr->pst,&intype,sizeof(intype),NULL,NULL);
    }

int xdrHeader(XDR *pxdrs,struct dump *ph)
    {
    if (!xdr_double(pxdrs,&ph->time)) return 0;
    if (!xdr_u_int(pxdrs,&ph->nbodies)) return 0;
    if (!xdr_u_int(pxdrs,&ph->ndim)) return 0;
    if (!xdr_u_int(pxdrs,&ph->nsph)) return 0;
    if (!xdr_u_int(pxdrs,&ph->ndark)) return 0;
    if (!xdr_u_int(pxdrs,&ph->nstar)) return 0;
    if (!xdr_u_int(pxdrs,&ph->pad)) return 0;
    return 1;
    }

double getTime(MSR msr, double dExpansion, double *dvFac) {
    double dTime,aTo,tTo,z;
    if (msr->param.csm->bComove) {
	if(msr->param.csm->dHubble0 == 0.0) {
	    printf("No hubble constant specified\n");
	    _msrExit(msr,1);
	    }
	dTime = csmExp2Time(msr->param.csm,dExpansion);
	z = 1.0/dExpansion - 1.0;
	if (msr->param.bVStart)
	    printf("Input file, Time:%g Redshift:%g Expansion factor:%g iStartStep:%d\n",
		   dTime,z,dExpansion,msr->param.iStartStep);
	if (prmSpecified(msr->prm,"dRedTo")) {
	  if (msr->param.dRedTo <= -1.0) {
	    printf("Badly specified final redshift (zTo <= -1.0), check -zto parameter.\n");
	    _msrExit(msr,1);
	  }
	    if (!prmArgSpecified(msr->prm,"nSteps") &&
		prmArgSpecified(msr->prm,"dDelta")) {
		aTo = 1.0/(msr->param.dRedTo + 1.0);
		tTo = csmExp2Time(msr->param.csm,aTo);
		if (msr->param.bVStart)
		    printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
			   tTo,1.0/aTo-1.0,aTo);
		if (tTo < dTime) {
		    printf("Badly specified final redshift, check -zto parameter.\n");
		    _msrExit(msr,1);
		    }
		msr->param.nSteps = (int)ceil((tTo-dTime)/msr->param.dDelta);
		msr->param.dDelta =
		  (tTo-dTime)/(msr->param.nSteps -
			       msr->param.iStartStep);
		}
	    else if (!prmArgSpecified(msr->prm,"dDelta") &&
		     prmArgSpecified(msr->prm,"nSteps")) {
		aTo = 1.0/(msr->param.dRedTo + 1.0);
		tTo = csmExp2Time(msr->param.csm,aTo);
		if (msr->param.bVStart)
		    printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
			   tTo,1.0/aTo-1.0,aTo);
		if (tTo < dTime) {
		    printf("Badly specified final redshift, check -zto parameter.\n");	
		    _msrExit(msr,1);
		    }
		if(msr->param.nSteps != 0)
		    msr->param.dDelta =
			(tTo-dTime)/(msr->param.nSteps -
				     msr->param.iStartStep);
				
		else
		    msr->param.dDelta = 0.0;
		}
	    else if (!prmSpecified(msr->prm,"nSteps") &&
		     prmFileSpecified(msr->prm,"dDelta")) {
		aTo = 1.0/(msr->param.dRedTo + 1.0);
		tTo = csmExp2Time(msr->param.csm,aTo);
		if (msr->param.bVStart)
		    printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
			   tTo,1.0/aTo-1.0,aTo);
		if (tTo < dTime) {
		    printf("Badly specified final redshift, check -zto parameter.\n");
		    _msrExit(msr,1);
		    }
		msr->param.nSteps = (int)ceil((tTo-dTime)/msr->param.dDelta);
		msr->param.dDelta =
		  (tTo-dTime)/(msr->param.nSteps -
			       msr->param.iStartStep);
		}
	    else if (!prmSpecified(msr->prm,"dDelta") &&
		     prmFileSpecified(msr->prm,"nSteps")) {
		aTo = 1.0/(msr->param.dRedTo + 1.0);
		tTo = csmExp2Time(msr->param.csm,aTo);
		if (msr->param.bVStart)
		    printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
			   tTo,1.0/aTo-1.0,aTo);
		if (tTo < dTime) {
		    printf("Badly specified final redshift, check -zto parameter.\n");	
		    _msrExit(msr,1);
		    }
		if(msr->param.nSteps != 0)
		    msr->param.dDelta =	(tTo-dTime)/(msr->param.nSteps
						     - msr->param.iStartStep);
		else
		    msr->param.dDelta = 0.0;
		}
	    }
	else {
	    tTo = dTime + msr->param.nSteps*msr->param.dDelta;
	    aTo = csmTime2Exp(msr->param.csm,tTo);
	    if (msr->param.bVStart)
		printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
		       tTo,1.0/aTo-1.0,aTo);
	    }
	if (msr->param.bCannonical) {
	    *dvFac = dExpansion*dExpansion;
	    }
	else {
	    *dvFac = 1.0;
	    }
	}
    else {
	dTime = dExpansion;
	if (msr->param.bVStart) printf("Input file, Time:%g iStartStep:%d\n",dTime,msr->param.iStartStep);
	tTo = dTime + (msr->param.nSteps - msr->param.iStartStep)*msr->param.dDelta;
	if (msr->param.bVStart) {
	    printf("Simulation to Time:%g\n",tTo);
	    }
	*dvFac = 1.0;
	}

    return dTime;
}

#ifdef USE_HDF5
#ifdef USE_MDL_IO
static double _msrIORead(MSR msr, const char *achFilename, int iStep )
{
    LCL *plcl = msr->pst->plcl;
    struct inStartLoad load;
    struct inIOLoad in;
    struct inPlanLoad inPlan;
    struct outPlanLoad outPlan;
    struct inSetParticleTypes intype;
    double dTime, dvFac;
    int nPlan;
    int outSize;
    char *p1, *p2;
    int n, i;
    total_t N;

    /* Load the files into the I/O processors */
    strcpy( load.achInName, achFilename );
    p1 = strstr( load.achInName, "&N" );
    if ( p1 ) {
	p2 = strstr( achFilename, "&N" );
	assert( p2 != NULL );
	n = p1 - load.achInName;
	strcpy( p1, msrOutName(msr) );
	strcat( p1, p2+2 );
    }
    p1 = strstr( load.achInName, "&S" );
    if ( p1 ) {
	p2 = strstr( achFilename, "&S" );
	assert( p2 != NULL );
	n = p1 - load.achInName;
	sprintf( p1, "%05d", iStep );
	strcat( p1+5, p2+2 );
    }

    /* First we ask the I/O processors to see how many files we have, and how
       many particles there are in each file. The number of files can be
       different from the number of I/O processors (if, for example, the run
       was started with a different number of processors) */
    strcpy(inPlan.achInName,load.achInName);
    mdlSetComm(msr->mdl,1);
    mdlReqService(msr->mdl,0,IO_PLAN_LOAD,&inPlan,sizeof(inPlan));
    mdlGetReply(msr->mdl,0,&outPlan,&nPlan);
    mdlSetComm(msr->mdl,0);

    /* Count the total number of particles */
    load.nFiles = 0;
    N = 0;
    for( i=0; i<MDL_MAX_IO_PROCS; i++ ) {
	if ( outPlan.nCount[i] ) {
	    assert( load.nFiles == i );
	    load.nFiles++;
	}
	load.nCount[i] = outPlan.nCount[i];
	N += outPlan.nCount[i];
    }

    msr->nDark = N;
    msr->nGas  = 0;
    msr->nStar = 0;
    msr->N = msr->nDark+msr->nGas+msr->nStar;
    msr->nMaxOrder = msr->N;
    msr->nMaxOrderGas = msr->nGas;
    msr->nMaxOrderDark = msr->nGas + msr->nDark;

    dTime = getTime(msr,outPlan.dExpansion,&dvFac);
    if (msr->param.bVStart)
	printf("Reading file...\nN:%"PRIu64" nDark:%"PRIu64" nGas:%"PRIu64" nStar:%"PRIu64"\n",msr->N,
		   msr->nDark,msr->nGas,msr->nStar);

    /* Fire up the load process.  The I/O processors will enter start an
     mdlSend() eventually (after loading the data). */
    mdlSetComm(msr->mdl,1);
    mdlReqService(msr->mdl,0,IO_START_LOAD,&load,sizeof(load));
    mdlSetComm(msr->mdl,0);

    in.nDark = msr->nDark;
    in.nGas = msr->nGas;
    in.nStar = msr->nStar;
    in.dvFac = dvFac;
    in.fExtraStore = msr->param.dExtraStore;
    /*
    ** Provide the period.
    */
    in.fPeriod[0] = msr->param.dxPeriod;
    in.fPeriod[1] = msr->param.dyPeriod;
    in.fPeriod[2] = msr->param.dzPeriod;
    in.nBucket = msr->param.nBucket;

    pstIOLoad( msr->pst, &in, sizeof(in), NULL, NULL );

    mdlSetComm(msr->mdl,1);
    mdlGetReply(msr->mdl,0,NULL,NULL);
    mdlSetComm(msr->mdl,0);

    pstSetParticleTypes(msr->pst,&intype,sizeof(intype),NULL,NULL);

    if (msr->param.bVDetails) puts("Input file has been successfully read.");
    /*
    ** Now read in the output points, passing the initial time.
    ** We do this only if nSteps is not equal to zero.
    */
    if (msrSteps(msr) > 0) msrReadOuts(msr,dTime);
    /*
    ** Set up the output counter.
    */
    for (msr->iOut=0;msr->iOut<msr->nOuts;++msr->iOut) {
	if (dTime < msr->pdOutTime[msr->iOut]) break;
	}


    return dTime;
}
#endif

static double _msrReadHDF5(MSR msr, const char *achFilename)
{
    hid_t fileID;
    IOHDF5 io;

    //struct dump h;
    double dExpansion;
    struct inReadTipsy in;
    char achInFile[PST_FILENAME_SIZE];
    LCL *plcl = msr->pst->plcl;
    double dTime;
    struct inSetParticleTypes intype;
    //uint64_t tlong;

    strcpy(in.achInFile,achFilename);

    /* Add local Data Path. */
    _msrMakePath(plcl->pszDataPath,achFilename,achInFile);

    fileID=H5Fopen(achInFile, H5F_ACC_RDONLY, H5P_DEFAULT);
    if ( fileID < 0 ) {
	printf("Could not open InFile:%s\n",achInFile);
	_msrExit(msr,1);
    }
    io = ioHDF5Initialize( fileID, 32768, IOHDF5_SINGLE );

    assert(ioHDF5ReadAttribute( io, "dTime", H5T_NATIVE_DOUBLE, &dExpansion ));

    msr->nDark = ioHDF5DarkCount(io);
    msr->nGas = ioHDF5GasCount(io);
    msr->nStar = ioHDF5StarCount(io);
    msr->N = msr->nDark+msr->nGas+msr->nStar;

    msr->nMaxOrder = msr->N;
    msr->nMaxOrderGas = msr->nGas;
    msr->nMaxOrderDark = msr->nGas + msr->nDark;
    assert(msr->N == msr->nDark+msr->nGas+msr->nStar);

    ioHDF5Finish(io);
    H5Fclose(fileID);

    dTime = getTime(msr,dExpansion,&in.dvFac);
    if (msr->param.bVStart)
	printf("Reading file...\nN:%"PRIu64" nDark:%"PRIu64" nGas:%"PRIu64" nStar:%"PRIu64"\n",msr->N,
		   msr->nDark,msr->nGas,msr->nStar);

    in.nFileStart = 0;
    in.nFileEnd = msr->N - 1;
    in.nBucket = msr->param.nBucket;
    in.nDark = msr->nDark;
    in.nGas = msr->nGas;
    in.nStar = msr->nStar;
    in.bStandard = msr->param.bStandard;
    in.bDoublePos = msr->param.bDoublePos;
    strcpy(in.achOutName,msr->param.achOutName);    
    /*
    ** Since pstReadTipsy causes the allocation of the local particle
    ** store, we need to tell it the percentage of extra storage it
    ** should allocate for load balancing differences in the number of
    ** particles.
    */
    in.fExtraStore = msr->param.dExtraStore;
    /*
    ** Provide the period.
    */
    in.fPeriod[0] = msr->param.dxPeriod;
    in.fPeriod[1] = msr->param.dyPeriod;
    in.fPeriod[2] = msr->param.dzPeriod;

    if(msr->param.bParaRead)
	pstReadHDF5(msr->pst,&in,sizeof(in),NULL,NULL);
    else
	msrOneNodeReadHDF5(msr, &in);
    pstSetParticleTypes(msr->pst, &intype, sizeof(intype), NULL, NULL);
    if (msr->param.bVDetails) puts("Input file has been successfully read.");
    /*
    ** Now read in the output points, passing the initial time.
    ** We do this only if nSteps is not equal to zero.
    */
    if (msrSteps(msr) > 0) msrReadOuts(msr,dTime);
    /*
    ** Set up the output counter.
    */
    for (msr->iOut=0;msr->iOut<msr->nOuts;++msr->iOut) {
	if (dTime < msr->pdOutTime[msr->iOut]) break;
	}
    return(dTime);
}
#endif

#ifdef USE_GRAFIC
double msrGenerateIC(MSR msr)
{
    struct inGenerateIC in;
    struct outGenerateIC out;
    struct inSetParticleTypes intype;
    int nOut;
    double sec,dsec;
    double dvFac;

    in.h = msr->param.h;
    in.dBoxSize = msr->param.dBoxSize;
    in.iSeed = msr->param.iSeed;
    in.nGrid = msr->param.nGrid;
    in.omegac= msr->param.csm->dOmega0 - msr->param.csm->dOmegab;
    in.omegab= msr->param.csm->dOmegab;
    in.omegav= msr->param.csm->dLambda;
    in.bCannonical = msr->param.bCannonical && msr->param.csm->bComove;
    in.fExtraStore = msr->param.dExtraStore;
    in.fPeriod[0] = msr->param.dxPeriod;
    in.fPeriod[1] = msr->param.dyPeriod;
    in.fPeriod[2] = msr->param.dzPeriod;
    in.nBucket = msr->param.nBucket;

    msr->nDark = in.nGrid * in.nGrid * in.nGrid;
    msr->nGas  = 0;
    msr->nStar = 0;
    msr->N = msr->nDark+msr->nGas+msr->nStar;
    msr->nMaxOrder = msr->N;
    msr->nMaxOrderGas = msr->nGas;
    msr->nMaxOrderDark = msr->nGas + msr->nDark;

    if (msr->param.bVStart)
	printf("Generating IC...\nN:%"PRIu64" nDark:%"PRIu64
	       " nGas:%"PRIu64" nStar:%"PRIu64"\n",
	       msr->N, msr->nDark,msr->nGas,msr->nStar);

    sec = msrTime();
    pstGenerateIC(msr->pst,&in,sizeof(in),&out,&nOut);
    dsec = msrTime() - sec;
    if (msr->param.bVDetails) {
	dsec = msrTime() - sec;
	printf("IC Generation Complete, Wallclock: %f secs\n\n",dsec);
	}

    pstSetParticleTypes(msr->pst, &intype, sizeof(intype), NULL, NULL);

    return getTime(msr,out.dExpansion,&dvFac);
}
#endif

static double _msrReadTipsy(MSR msr, const char *achFilename)
    {
    FILE *fp;
    struct dump h;
    double dExpansion;
    struct inReadTipsy in;
    char achInFile[PST_FILENAME_SIZE];
    LCL *plcl = msr->pst->plcl;
    double dTime,aTo,tTo,z;
    struct inSetParticleTypes intype;
    uint64_t tlong;

    strcpy(in.achInFile,achFilename);


    /* Add local Data Path. */
    _msrMakePath(plcl->pszDataPath,achFilename,achInFile);

	
    fp = fopen(achInFile,"r");
    
    if (!fp) {
	printf("Could not open InFile:%s\n",achInFile);
	_msrExit(msr,1);
    }

    /*
    ** Assume tipsy format for now, and dark matter only.
    */
    if (msr->param.bStandard) {
	XDR xdrs;

	xdrstdio_create(&xdrs,fp,XDR_DECODE);
	xdrHeader(&xdrs,&h);
	xdr_destroy(&xdrs);
	}
    else {
	fread(&h,sizeof(struct dump),1,fp);
	}
    fclose(fp);
    msr->N = h.nbodies;
    msr->nDark = h.ndark;
    msr->nGas = h.nsph;
    msr->nStar = h.nstar;
    dExpansion = h.time;

    if (h.pad) {
	/*
	** This means we are dealing with *large* tipsy files
	** The pad field provides 8 more high bits to each of 
	** nbodies, nsph, ndark and nstar.
	*/
	tlong = h.pad & 0x000000ff;
	tlong = tlong << 32;
	msr->N += tlong;
	tlong = h.pad & 0x0000ff00;
	tlong = tlong << 24;
	msr->nGas += tlong;
	tlong = h.pad & 0x00ff0000;
	tlong = tlong << 16;
	msr->nDark += tlong;
	tlong = h.pad & 0xff000000;
	tlong = tlong << 8;
	msr->nStar += tlong;
	}

    msr->nMaxOrder = msr->N;
    msr->nMaxOrderGas = msr->nGas;
    msr->nMaxOrderDark = msr->nGas + msr->nDark;
    assert(msr->N == msr->nDark+msr->nGas+msr->nStar);
    if (msr->param.csm->bComove) {
	if(msr->param.csm->dHubble0 == 0.0) {
	    printf("No hubble constant specified\n");
	    _msrExit(msr,1);
	    }
	dTime = csmExp2Time(msr->param.csm,dExpansion);
	z = 1.0/dExpansion - 1.0;
	if (msr->param.bVStart)
	    printf("Input file, Time:%g Redshift:%g Expansion factor:%g iStartStep:%d\n",
		   dTime,z,dExpansion,msr->param.iStartStep);
	if (prmSpecified(msr->prm,"dRedTo")) {
	  if (msr->param.dRedTo <= -1.0) {
	    printf("Badly specified final redshift (zTo <= -1.0), check -zto parameter.\n");
	    _msrExit(msr,1);
	  }
	    if (!prmArgSpecified(msr->prm,"nSteps") &&
		prmArgSpecified(msr->prm,"dDelta")) {
		aTo = 1.0/(msr->param.dRedTo + 1.0);
		tTo = csmExp2Time(msr->param.csm,aTo);
		if (msr->param.bVStart)
		    printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
			   tTo,1.0/aTo-1.0,aTo);
		if (tTo < dTime) {
		    printf("Badly specified final redshift, check -zto parameter.\n");
		    _msrExit(msr,1);
		    }
		msr->param.nSteps = (int)ceil((tTo-dTime)/msr->param.dDelta);
		msr->param.dDelta =
		  (tTo-dTime)/(msr->param.nSteps -
			       msr->param.iStartStep);
		}
	    else if (!prmArgSpecified(msr->prm,"dDelta") &&
		     prmArgSpecified(msr->prm,"nSteps")) {
		aTo = 1.0/(msr->param.dRedTo + 1.0);
		tTo = csmExp2Time(msr->param.csm,aTo);
		if (msr->param.bVStart)
		    printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
			   tTo,1.0/aTo-1.0,aTo);
		if (tTo < dTime) {
		    printf("Badly specified final redshift, check -zto parameter.\n");	
		    _msrExit(msr,1);
		    }
		if(msr->param.nSteps != 0)
		    msr->param.dDelta =
			(tTo-dTime)/(msr->param.nSteps -
				     msr->param.iStartStep);
				
		else
		    msr->param.dDelta = 0.0;
		}
	    else if (!prmSpecified(msr->prm,"nSteps") &&
		     prmFileSpecified(msr->prm,"dDelta")) {
		aTo = 1.0/(msr->param.dRedTo + 1.0);
		tTo = csmExp2Time(msr->param.csm,aTo);
		if (msr->param.bVStart)
		    printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
			   tTo,1.0/aTo-1.0,aTo);
		if (tTo < dTime) {
		    printf("Badly specified final redshift, check -zto parameter.\n");
		    _msrExit(msr,1);
		    }
		msr->param.nSteps = (int)ceil((tTo-dTime)/msr->param.dDelta);
		msr->param.dDelta =
		  (tTo-dTime)/(msr->param.nSteps -
			       msr->param.iStartStep);
		}
	    else if (!prmSpecified(msr->prm,"dDelta") &&
		     prmFileSpecified(msr->prm,"nSteps")) {
		aTo = 1.0/(msr->param.dRedTo + 1.0);
		tTo = csmExp2Time(msr->param.csm,aTo);
		if (msr->param.bVStart)
		    printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
			   tTo,1.0/aTo-1.0,aTo);
		if (tTo < dTime) {
		    printf("Badly specified final redshift, check -zto parameter.\n");	
		    _msrExit(msr,1);
		    }
		if(msr->param.nSteps != 0)
		    msr->param.dDelta =	(tTo-dTime)/(msr->param.nSteps
						     - msr->param.iStartStep);
		else
		    msr->param.dDelta = 0.0;
		}
	    }
	else {
	    tTo = dTime + msr->param.nSteps*msr->param.dDelta;
	    aTo = csmTime2Exp(msr->param.csm,tTo);
	    if (msr->param.bVStart)
		printf("Simulation to Time:%g Redshift:%g Expansion factor:%g\n",
		       tTo,1.0/aTo-1.0,aTo);
	    }
	if (msr->param.bCannonical) {
	    in.dvFac = dExpansion*dExpansion;
	    }
	else {
	    in.dvFac = 1.0;
	    }
	}
    else {
	dTime = dExpansion;
	if (msr->param.bVStart) printf("Input file, Time:%g iStartStep:%d\n",dTime,msr->param.iStartStep);
	tTo = dTime + (msr->param.nSteps - msr->param.iStartStep)*msr->param.dDelta;
	if (msr->param.bVStart) {
	    printf("Simulation to Time:%g\n",tTo);
	    printf("Reading file...\nN:%"PRIu64" nDark:%"PRIu64" nGas:%"PRIu64" nStar:%"PRIu64" Time:%g\n",
		   msr->N,msr->nDark,msr->nGas,msr->nStar,dTime);
	    }
	in.dvFac = 1.0;
	}
    in.nFileStart = 0;
    in.nFileEnd = msr->N - 1;
    in.nBucket = msr->param.nBucket;
    in.nDark = msr->nDark;
    in.nGas = msr->nGas;
    in.nStar = msr->nStar;
    in.bStandard = msr->param.bStandard;
    in.bDoublePos = msr->param.bDoublePos;
    strcpy(in.achOutName,msr->param.achOutName);    
    /*
    ** Since pstReadTipsy causes the allocation of the local particle
    ** store, we need to tell it the percentage of extra storage it
    ** should allocate for load balancing differences in the number of
    ** particles.
    */
    in.fExtraStore = msr->param.dExtraStore;
    /*
    ** Provide the period.
    */
    in.fPeriod[0] = msr->param.dxPeriod;
    in.fPeriod[1] = msr->param.dyPeriod;
    in.fPeriod[2] = msr->param.dzPeriod;

    if(msr->param.bParaRead)
	pstReadTipsy(msr->pst,&in,sizeof(in),NULL,NULL);
    else
	msrOneNodeReadTipsy(msr, &in);
    pstSetParticleTypes(msr->pst, &intype, sizeof(intype), NULL, NULL);
    if (msr->param.bVDetails) puts("Input file has been successfully read.");
    /*
    ** Now read in the output points, passing the initial time.
    ** We do this only if nSteps is not equal to zero.
    */
    if (msrSteps(msr) > 0) msrReadOuts(msr,dTime);
    /*
    ** Set up the output counter.
    */
    for (msr->iOut=0;msr->iOut<msr->nOuts;++msr->iOut) {
	if (dTime < msr->pdOutTime[msr->iOut]) break;
	}
    return(dTime);
    }


#ifdef USE_MDL_IO
void msrIOWrite(MSR msr, const char *achOutName, double dTime, int bCheckpoint)
{
    double dExp, dvFac;

    struct inStartIO inStart;
    struct inStartSave save;

    if (msr->param.csm->bComove) {
	dExp = csmTime2Exp(msr->param.csm,dTime);
	if (msr->param.bCannonical) {
	    dvFac = 1.0/(dExp*dExp);
	    }
	else {
	    dvFac = 1.0;
	    }
	}
    else {
	dExp = dTime;
	dvFac = 1.0;
	}

    /* Ask the I/O processors to start a save operation */
    save.dTime       = dExp;
    save.N           = msr->N;                                         /* Total */
    save.bCheckpoint = bCheckpoint;
    save.dEcosmo     = msr->dEcosmo;
    save.dTimeOld    = msr->dTimeOld;
    save.dUOld       = msr->dUOld;

    strcpy(save.achOutName,achOutName); 
    mdlSetComm(msr->mdl,1);
    if ( msr->bSavePending )
	mdlGetReply(msr->mdl,0,NULL,NULL);
    mdlReqService(msr->mdl,0,IO_START_SAVE,&save,sizeof(save));
    msr->bSavePending = 1;
    mdlSetComm(msr->mdl,0);

    /* Execute a save operation on the worker processors */


    inStart.dTime = dExp;
    inStart.dvFac = dvFac;
    inStart.bDoublePos = msr->param.bDoublePos;
    inStart.N = msr->N;                                                /* Total */
    strcpy(inStart.achOutName,achOutName);

    inStart.dEcosmo  = msr->dEcosmo;
    inStart.dTimeOld = msr->dTimeOld;
    inStart.dUOld    = msr->dUOld;

    //_msrMakePath(plcl->pszDataPath,in->achOutFile,achOutFile);


    /*char achOutFile[PST_FILENAME_SIZE];*/

    pstStartIO( msr->pst, &inStart, sizeof(inStart), NULL, NULL );

#if 0
    /* Get the reply from the I/O processor */
    mdlSetComm(msr->mdl,1);
    mdlGetReply(msr->mdl,0,NULL,NULL);
    mdlSetComm(msr->mdl,0);
#endif
}
#endif

#ifdef USE_HDF5
/*
** This function saves all of the input parameters, as well as single-variable
** state information.
*/
void msrSaveParameters(MSR msr, IOHDF5 io)
{
    ioHDF5WriteAttribute( io, "nThreads", H5T_NATIVE_INT, &msr->param.nThreads );
    ioHDF5WriteAttribute( io, "bDiag", H5T_NATIVE_INT, &msr->param.bDiag );
    ioHDF5WriteAttribute( io, "bOverwrite", H5T_NATIVE_INT, &msr->param.bOverwrite );
    ioHDF5WriteAttribute( io, "bVWarnings", H5T_NATIVE_INT, &msr->param.bVWarnings );
    ioHDF5WriteAttribute( io, "bVStart", H5T_NATIVE_INT, &msr->param.bVStart );
    ioHDF5WriteAttribute( io, "bVStep", H5T_NATIVE_INT, &msr->param.bVStep );
    ioHDF5WriteAttribute( io, "bVRungStat", H5T_NATIVE_INT, &msr->param.bVRungStat );
    ioHDF5WriteAttribute( io, "bVDetails", H5T_NATIVE_INT, &msr->param.bVDetails );
    ioHDF5WriteAttribute( io, "bPeriodic", H5T_NATIVE_INT, &msr->param.bPeriodic );
    ioHDF5WriteAttribute( io, "bParaRead", H5T_NATIVE_INT, &msr->param.bParaRead );
    ioHDF5WriteAttribute( io, "bParaWrite", H5T_NATIVE_INT, &msr->param.bParaWrite );
    ioHDF5WriteAttribute( io, "bCannonical", H5T_NATIVE_INT, &msr->param.bCannonical );
    ioHDF5WriteAttribute( io, "bStandard", H5T_NATIVE_INT, &msr->param.bStandard );
    ioHDF5WriteAttribute( io, "bDoublePos", H5T_NATIVE_INT, &msr->param.bDoublePos );
    ioHDF5WriteAttribute( io, "bGravStep", H5T_NATIVE_INT, &msr->param.bGravStep );
    ioHDF5WriteAttribute( io, "bEpsAccStep", H5T_NATIVE_INT, &msr->param.bEpsAccStep );
    ioHDF5WriteAttribute( io, "bSqrtPhiStep", H5T_NATIVE_INT, &msr->param.bSqrtPhiStep );
    ioHDF5WriteAttribute( io, "bAccelStep", H5T_NATIVE_INT, &msr->param.bAccelStep );
    ioHDF5WriteAttribute( io, "bDensityStep", H5T_NATIVE_INT, &msr->param.bDensityStep );
    ioHDF5WriteAttribute( io, "iTimeStepCrit", H5T_NATIVE_INT, &msr->param.iTimeStepCrit );
    //ioHDF5WriteAttribute( io, "nPColl", H5T_NATIVE_INT, &msr->param.nPColl );
    ioHDF5WriteAttribute( io, "nTruncateRung", H5T_NATIVE_INT, &msr->param.nTruncateRung );
    ioHDF5WriteAttribute( io, "bDoDensity", H5T_NATIVE_INT, &msr->param.bDoDensity );
#ifdef USE_PNG
    ioHDF5WriteAttribute( io, "nPNGResolution", H5T_NATIVE_INT, &msr->param.nPNGResolution );
#endif
    ioHDF5WriteAttribute( io, "bDodtOutput", H5T_NATIVE_INT, &msr->param.bDodtOutput );
    ioHDF5WriteAttribute( io, "bDoRungOutput", H5T_NATIVE_INT, &msr->param.bDoRungOutput );
    ioHDF5WriteAttribute( io, "bDoGravity", H5T_NATIVE_INT, &msr->param.bDoGravity );
    ioHDF5WriteAttribute( io, "bAntiGrav", H5T_NATIVE_INT, &msr->param.bAntiGrav );
    ioHDF5WriteAttribute( io, "nBucket", H5T_NATIVE_INT, &msr->param.nBucket );
    ioHDF5WriteAttribute( io, "dFracNoTreeSqueeze", H5T_NATIVE_DOUBLE, &msr->param.dFracNoTreeSqueeze );
    ioHDF5WriteAttribute( io, "iOutInterval", H5T_NATIVE_INT, &msr->param.iOutInterval );
    ioHDF5WriteAttribute( io, "iCheckInterval", H5T_NATIVE_INT, &msr->param.iCheckInterval );
    ioHDF5WriteAttribute( io, "iLogInterval", H5T_NATIVE_INT, &msr->param.iLogInterval );
    ioHDF5WriteAttribute( io, "iOrder", H5T_NATIVE_INT, &msr->param.iOrder );
    ioHDF5WriteAttribute( io, "bEwald", H5T_NATIVE_INT, &msr->param.bEwald );
    ioHDF5WriteAttribute( io, "iEwOrder", H5T_NATIVE_INT, &msr->param.iEwOrder );
    ioHDF5WriteAttribute( io, "nReplicas", H5T_NATIVE_INT, &msr->param.nReplicas );
    ioHDF5WriteAttribute( io, "iStartStep", H5T_NATIVE_INT, &msr->param.iStartStep );
    ioHDF5WriteAttribute( io, "nSteps", H5T_NATIVE_INT, &msr->param.nSteps );
    ioHDF5WriteAttribute( io, "nSmooth", H5T_NATIVE_INT, &msr->param.nSmooth );
    ioHDF5WriteAttribute( io, "iMaxRung", H5T_NATIVE_INT, &msr->param.iMaxRung );
    ioHDF5WriteAttribute( io, "nRungVeryActive", H5T_NATIVE_INT, &msr->param.nRungVeryActive );
    ioHDF5WriteAttribute( io, "nPartVeryActive", H5T_NATIVE_INT, &msr->param.nPartVeryActive );
    ioHDF5WriteAttribute( io, "nGrowMass", H5T_NATIVE_INT, &msr->param.nGrowMass );
    ioHDF5WriteAttribute( io, "iWallRunTime", H5T_NATIVE_INT, &msr->param.iWallRunTime );
    ioHDF5WriteAttribute( io, "bPhysicalSoft", H5T_NATIVE_INT, &msr->param.bPhysicalSoft );  
    ioHDF5WriteAttribute( io, "bSoftMaxMul", H5T_NATIVE_INT, &msr->param.bSoftMaxMul );
    ioHDF5WriteAttribute( io, "bVariableSoft", H5T_NATIVE_INT, &msr->param.bVariableSoft );
    ioHDF5WriteAttribute( io, "nSoftNbr", H5T_NATIVE_INT, &msr->param.nSoftNbr );
    ioHDF5WriteAttribute( io, "bSoftByType", H5T_NATIVE_INT, &msr->param.bSoftByType );
    ioHDF5WriteAttribute( io, "bDoSoftOutput", H5T_NATIVE_INT, &msr->param.bDoSoftOutput );

    ioHDF5WriteAttribute( io, "dEta", H5T_NATIVE_DOUBLE, &msr->param.dEta );
    ioHDF5WriteAttribute( io, "dExtraStore", H5T_NATIVE_DOUBLE, &msr->param.dExtraStore );
    ioHDF5WriteAttribute( io, "dSoft", H5T_NATIVE_DOUBLE, &msr->param.dSoft );
    ioHDF5WriteAttribute( io, "dSoftMax", H5T_NATIVE_DOUBLE, &msr->param.dSoftMax );
    ioHDF5WriteAttribute( io, "dDelta", H5T_NATIVE_DOUBLE, &msr->param.dDelta );
    ioHDF5WriteAttribute( io, "dEwCut", H5T_NATIVE_DOUBLE, &msr->param.dEwCut );
    ioHDF5WriteAttribute( io, "dEwhCut", H5T_NATIVE_DOUBLE, &msr->param.dEwhCut );
    ioHDF5WriteAttribute( io, "dTheta", H5T_NATIVE_DOUBLE, &msr->param.dTheta );
    ioHDF5WriteAttribute( io, "dTheta2", H5T_NATIVE_DOUBLE, &msr->param.dTheta2 );
    ioHDF5WriteAttribute( io, "daSwitchTheta", H5T_NATIVE_DOUBLE, &msr->param.daSwitchTheta );
    ioHDF5WriteAttribute( io, "dPeriod", H5T_NATIVE_DOUBLE, &msr->param.dPeriod );
    ioHDF5WriteAttribute( io, "dxPeriod", H5T_NATIVE_DOUBLE, &msr->param.dxPeriod );
    ioHDF5WriteAttribute( io, "dyPeriod", H5T_NATIVE_DOUBLE, &msr->param.dyPeriod );
    ioHDF5WriteAttribute( io, "dzPeriod", H5T_NATIVE_DOUBLE, &msr->param.dzPeriod );
    ioHDF5WriteAttribute( io, "bComove", H5T_NATIVE_INT, &msr->param.csm->bComove );
    ioHDF5WriteAttribute( io, "dHubble0", H5T_NATIVE_DOUBLE, &msr->param.csm->dHubble0 );
    ioHDF5WriteAttribute( io, "dOmega0", H5T_NATIVE_DOUBLE, &msr->param.csm->dOmega0 );
    ioHDF5WriteAttribute( io, "dLambda", H5T_NATIVE_DOUBLE, &msr->param.csm->dLambda );
    ioHDF5WriteAttribute( io, "dOmegaDE", H5T_NATIVE_DOUBLE, &msr->param.csm->dOmegaDE );
    ioHDF5WriteAttribute( io, "w0", H5T_NATIVE_DOUBLE, &msr->param.csm->w0 );
    ioHDF5WriteAttribute( io, "wa", H5T_NATIVE_DOUBLE, &msr->param.csm->wa );
    ioHDF5WriteAttribute( io, "dOmegaRad", H5T_NATIVE_DOUBLE, &msr->param.csm->dOmegaRad );
    ioHDF5WriteAttribute( io, "dOmegab", H5T_NATIVE_DOUBLE, &msr->param.csm->dOmegab );
    ioHDF5WriteAttribute( io, "dRedTo", H5T_NATIVE_DOUBLE, &msr->param.dRedTo );
    ioHDF5WriteAttribute( io, "dCentMass", H5T_NATIVE_DOUBLE, &msr->param.dCentMass );
    ioHDF5WriteAttribute( io, "dGrowDeltaM", H5T_NATIVE_DOUBLE, &msr->param.dGrowDeltaM );
    ioHDF5WriteAttribute( io, "dGrowStartT", H5T_NATIVE_DOUBLE, &msr->param.dGrowStartT );
    ioHDF5WriteAttribute( io, "dGrowEndT", H5T_NATIVE_DOUBLE, &msr->param.dGrowEndT );
    ioHDF5WriteAttribute( io, "dFracNoDomainDecomp", H5T_NATIVE_DOUBLE, &msr->param.dFracNoDomainDecomp );
    ioHDF5WriteAttribute( io, "dFracNoDomainRootFind", H5T_NATIVE_DOUBLE, &msr->param.dFracNoDomainRootFind );
    ioHDF5WriteAttribute( io, "dFracNoDomainDimChoice", H5T_NATIVE_DOUBLE, &msr->param.dFracNoDomainDimChoice );

    /* Restart information */
    ioHDF5WriteAttribute( io, "dEcosmo", H5T_NATIVE_DOUBLE, &msr->dEcosmo );
    ioHDF5WriteAttribute( io, "dTimeOld", H5T_NATIVE_DOUBLE, &msr->dTimeOld );
    ioHDF5WriteAttribute( io, "dUOld", H5T_NATIVE_DOUBLE, &msr->dUOld );
}
#endif

/*
** This function makes some DANGEROUS assumptions!!!
** Main problem is that it calls pkd level routines, bypassing the
** pst level. It uses plcl pointer which is not desirable.
*/
void msrOneNodeWriteTipsy(MSR msr, struct inWriteTipsy *in, int bCheckpoint)
    {
    int i,id;
    uint64_t nStart;
    PST pst0;
    LCL *plcl;
    char achOutFile[PST_FILENAME_SIZE];
    int inswap;
#ifdef USE_HDF5
    hid_t fileID;
    IOHDF5 io;
    IOHDF5V ioDen, ioPot;
#endif

    pst0 = msr->pst;
    while(pst0->nLeaves > 1)
	pst0 = pst0->pstLower;
    plcl = pst0->plcl;
    /*
    ** Add the local Data Path to the provided filename.
    */
    _msrMakePath(plcl->pszDataPath,in->achOutFile,achOutFile);

    /* 
     * First write our own particles.
     */
#ifdef USE_HDF5
    if ( in->bStandard == 2 ) {
	printf( "Writing HDF5 format file to %s\n", achOutFile );
	fileID=H5Fcreate(achOutFile, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if ( fileID < 0 ) {
	    fprintf(stderr,"Unable to create %s\n", achOutFile);
	    H5assert(fileID);
	}
	io = ioHDF5Initialize( fileID, 32768, bCheckpoint );
	ioDen  = ioHDFF5NewVector( io, "density",  IOHDF5_SINGLE );
	ioPot  = ioHDFF5NewVector( io, "potential",IOHDF5_SINGLE );
	ioHDF5WriteAttribute( io, "dTime", H5T_NATIVE_DOUBLE, &in->dTime );
	msrSaveParameters(msr,io);
	pkdWriteHDF5(plcl->pkd, io, ioDen, ioPot, in->dvFac );
    }
    else
#endif
    {
	pkdWriteTipsy(plcl->pkd,achOutFile,plcl->nWriteStart,in->bStandard,
		  in->dvFac,in->bDoublePos); 
    }
    nStart = plcl->pkd->nLocal;
    assert(msr->pMap[0] == 0);
    for (i=1;i<msr->nThreads;++i) {
	id = msr->pMap[i];
	/* 
	 * Swap particles with the remote processor.
	 */
	inswap = 0;
	mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
	pkdSwapAll(plcl->pkd, id);
	mdlGetReply(pst0->mdl,id,NULL,NULL);
	/* 
	 * Write the swapped particles.
	 */
#ifdef USE_HDF5
	if ( in->bStandard == 2 ) {
	    pkdWriteHDF5(plcl->pkd, io, ioDen, ioPot, in->dvFac);
	}
	else
#endif
	{
	    pkdWriteTipsy(plcl->pkd,achOutFile,nStart, in->bStandard,
			  in->dvFac, in->bDoublePos); 
	}
	nStart += plcl->pkd->nLocal;
	/* 
	 * Swap them back again.
	 */
	inswap = 0;
	mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
	pkdSwapAll(plcl->pkd, id);
	mdlGetReply(pst0->mdl,id,NULL,NULL);
    	}

#ifdef USE_HDF5
	if ( in->bStandard == 2 ) {
	    ioHDF5Finish(io);
	    H5assert(H5Fflush(fileID,H5F_SCOPE_GLOBAL));
	    H5assert(H5Fclose(fileID));
	}
#endif

    assert(nStart == msr->N);
    }


void msrCalcWriteStart(MSR msr) 
    {
    struct outSetTotal out;
    struct inSetWriteStart in;

    pstSetTotal(msr->pst,NULL,0,&out,NULL);
    assert(out.nTotal == msr->N);
    in.nWriteStart = 0;
    pstSetWriteStart(msr->pst,&in,sizeof(in),NULL,NULL);
    }


void _msrWriteTipsy(MSR msr,char *pszFileName,double dTime,int bCheckpoint)
    {
    FILE *fp;
    struct dump h;
    struct inWriteTipsy in;
    char achOutFile[PST_FILENAME_SIZE];
    LCL *plcl = msr->pst->plcl;
    /*
    ** Calculate where to start writing.
    ** This sets plcl->nWriteStart.
    */
    msrCalcWriteStart(msr);
    /*
    ** Add Data Subpath for local and non-local names.
    */
    _msrMakePath(msr->param.achDataSubPath,pszFileName,in.achOutFile);
    /*
    ** Add local Data Path.
    */
    _msrMakePath(plcl->pszDataPath,in.achOutFile,achOutFile);

    in.bStandard = msr->param.bStandard;
    in.bDoublePos = msr->param.bDoublePos;
    /*
    ** Assume tipsy format for now.
    */
    h.nbodies = msr->N;
    h.ndark = msr->nDark;
    h.nsph = msr->nGas;
    h.nstar = msr->nStar;
    if (msr->param.csm->bComove) {
	in.dTime = csmTime2Exp(msr->param.csm,dTime);
	if (msr->param.bCannonical) {
	    in.dvFac = 1.0/(in.dTime*in.dTime);
	    }
	else {
	    in.dvFac = 1.0;
	    }
	}
    else {
	in.dTime = dTime;
	in.dvFac = 1.0;
	}
    h.ndim = 3;
    h.time = in.dTime;
    if (msr->param.bVDetails) {
	if (msr->param.csm->bComove) {
	    printf("Writing file...\nTime:%g Redshift:%g\n",
		   dTime,(1.0/h.time - 1.0));
	    }
	else {
	    printf("Writing file...\nTime:%g\n",dTime);
	    }
	}

#ifdef USE_HDF5
    if ( msr->param.bHDF5 ) {
	in.bStandard = 2;
	msrOneNodeWriteTipsy(msr, &in,bCheckpoint);
    }
    else
#endif
    {
	fp = fopen(achOutFile,"w");
	if (!fp) {
	    printf("Could not open OutFile:%s\n",achOutFile);
	    _msrExit(msr,1);
	}
	if (in.bStandard) {
	    XDR xdrs;
	    
	    xdrstdio_create(&xdrs,fp,XDR_ENCODE);
	    xdrHeader(&xdrs,&h);
	    xdr_destroy(&xdrs);
	}
	else {
	    fwrite(&h,sizeof(struct dump),1,fp);
	}
	fclose(fp);
	
	if(msr->param.bParaWrite)
	    pstWriteTipsy(msr->pst,&in,sizeof(in),NULL,NULL);
	else
	    msrOneNodeWriteTipsy(msr, &in, bCheckpoint);
    }

    if (msr->param.bVDetails)
	puts("Output file has been successfully written.");
    }


void msrSetSoft(MSR msr,double dSoft)
    {
    struct inSetSoft in;
  
    if (msr->param.bVDetails) printf("Set Softening...\n");
    in.dSoft = dSoft;
    pstSetSoft(msr->pst,&in,sizeof(in),NULL,NULL);
    }


void msrDomainDecomp(MSR msr,int iRung,int bGreater,int bSplitVA) {
    struct inDomainDecomp in;
    struct outCalcBound outcb;
    uint64_t nActive;
    const uint64_t nDD = msr->N*msr->param.dFracNoDomainDecomp;
    const uint64_t nRT = msr->N*msr->param.dFracNoDomainRootFind;
    const uint64_t nSD = msr->N*msr->param.dFracNoDomainDimChoice;
    double sec,dsec;
    int iRungDD,iRungRT,iRungSD;
    int i,j;
    int bRestoreActive = 0;

    in.bDoRootFind = 1;
    in.bDoSplitDimFind = 1;
    in.nBndWrap[0] = 0;
    in.nBndWrap[1] = 0;
    in.nBndWrap[2] = 0;

    /*
    ** If we are dealing with a nice periodic volume in all
    ** three dimensions then we can set the initial bounds
    ** instead of calculating them.
    */
    if (msr->param.bPeriodic && 
	msr->param.dxPeriod < FLOAT_MAXVAL &&
	msr->param.dyPeriod < FLOAT_MAXVAL && 
	msr->param.dzPeriod < FLOAT_MAXVAL) {
	for (j=0;j<3;++j) {
	    in.bnd.fCenter[j] = msr->fCenter[j];
	    }
	in.bnd.fMax[0] = 0.5*msr->param.dxPeriod;
	in.bnd.fMax[1] = 0.5*msr->param.dyPeriod;
	in.bnd.fMax[2] = 0.5*msr->param.dzPeriod;
#ifdef USE_NEWDD
	assert(in.bnd.fMax[0] == in.bnd.fMax[1]);
	assert(in.bnd.fMax[1] == in.bnd.fMax[2]);
#endif	
	}
    else {
	pstCalcBound(msr->pst,NULL,0,&outcb,NULL);
	in.bnd = outcb.bnd;
	for (j=1;j<3;++j) {
	    in.bnd.fMax[0] = (in.bnd.fMax[j] > in.bnd.fMax[0])?in.bnd.fMax[j]:in.bnd.fMax[0];
	    }
	/*
	** Make the bounding box 10% larger in this case to accomodate particles drifting out.
	*/
	for (j=0;j<3;++j) {
	    in.bnd.fMax[j] = 1.10*in.bnd.fMax[0];
	    }
	}

    /*
    ** All of this could be calculated once for the case that the number
    ** of particles don't change. Or calculated every time the number of 
    ** particles does change.
    */
    assert(bGreater != 0);
    nActive = 0;
    iRungDD = 0;
    iRungRT = 0;
    iRungSD = 0;
    for (i=msr->iCurrMaxRung;i>=0;--i) {
	nActive += msr->nRung[i];
	if (nActive > nDD && !iRungDD) iRungDD = i;
	if (nActive > nRT && !iRungRT) iRungRT = i;
	if (nActive > nSD && !iRungSD) iRungSD = i;	
    }
    assert(iRungDD >= iRungRT);
    assert(iRungRT >= iRungSD);

    if (msr->iLastRungRT < 0) {
	/*
	** We need to do a full domain decompotition with iRungRT particles being active.
	** However, since I am not sure what the exact state of the domains can be at this point 
	** I had better do a split dim find as well.
	*/
	msr->iLastRungRT = iRungRT;  
	msrActiveRung(msr,iRungRT,bGreater);
	bRestoreActive = 1;
	in.bDoRootFind = 1;
	in.bDoSplitDimFind = 1;
    }
    else if (iRung > iRungDD && !bSplitVA) {
	if (msr->iLastRungRT < iRungRT) {
	    msr->iLastRungRT = iRungRT;
	    msrActiveRung(msr,iRungRT,bGreater); 
	    bRestoreActive = 1;
	    in.bDoRootFind = 1;
	    in.bDoSplitDimFind = 0;
	}
	else {
	    if (msr->param.bVRungStat) {
		printf("Skipping Domain Decomposition (nActive = %"PRIu64"/%"PRIu64", iRung:%d iRungDD:%d iLastRungRT:%d)\n",
		       msr->nActive,msr->N,iRung,iRungDD,msr->iLastRungRT);
	    }
	    return;  /* do absolutely nothing! */
	}
    }
    else if (iRung > iRungRT) {
	if (msr->iLastRungRT < iRungRT) {
	    msr->iLastRungRT = iRungRT;
	    msrActiveRung(msr,iRungRT,bGreater); 
	    bRestoreActive = 1;
	    in.bDoRootFind = 1;
	    in.bDoSplitDimFind = 0;
	}
	else {
	    if (msr->param.bVRungStat) {
		printf("Skipping Root Finder (nActive = %"PRIu64"/%"PRIu64", iRung:%d iRungRT:%d iLastRungRT:%d)\n",
		       msr->nActive,msr->N,iRung,iRungRT,msr->iLastRungRT);
	    }
	    in.bDoRootFind = 0;
	    in.bDoSplitDimFind = 0;
	}
    }
    else if (iRung > iRungSD) {
	if (msr->param.bVRungStat) {
	    printf("Skipping Domain Dim Choice (nActive = %"PRIu64"/%"PRIu64", iRung:%d iRungSD:%d iLastRungRT:%d)\n",
		   msr->nActive,msr->N,iRung,iRungSD,msr->iLastRungRT);
	}
	msr->iLastRungRT = iRung;
	in.bDoRootFind = 1;
	in.bDoSplitDimFind = 0;
    }
    else {
	msr->iLastRungRT = iRung;
	in.bDoRootFind = 1;
	in.bDoSplitDimFind = 1;
    }

    in.nActive = msr->nActive;
    in.nTotal = msr->N;

#ifdef USE_BSC
    MPItrace_event(10000, 2 );
#endif
    in.bSplitVA = bSplitVA;
    if (msr->param.bVDetails) {
	printf("Domain Decomposition: nActive (Rung %d) %"PRIu64" SplitVA:%d\n",
	       iRungDD,msr->nActive,bSplitVA);
	printf("Domain Decomposition... \n");
	sec = msrTime();
	}

    pstDomainDecomp(msr->pst,&in,sizeof(in),NULL,NULL);

#ifdef USE_BSC
    MPItrace_event(10000, 0 );
#endif
    if (msr->param.bVDetails) {
	dsec = msrTime() - sec;
	printf("Domain Decomposition complete, Wallclock: %f secs\n\n",dsec);
	}

    if (bRestoreActive) {
	/* Restore Active data */
	msrActiveRung(msr,iRung,bGreater);
	}
    }

/*
** This the meat of the tree build, but will be called by differently named
** functions in order to implement special features without recoding...
*/
void _BuildTree(MSR msr,double dMass,double dTimeStamp,int bExcludeVeryActive,int bNeedEwald) {
    struct inBuildTree in;
    struct ioCalcRoot root;
    KDN *pkdn;
    int iDum,nCell;

    if (msr->param.bVDetails) printf("Building local trees...\n\n");

    in.nBucket = msr->param.nBucket;
    in.diCrit2 = 1/(msr->dCrit*msr->dCrit);
    nCell = 1<<(1+(int)ceil(log((double)msr->nThreads)/log(2.0)));
    pkdn = malloc(nCell*sizeof(KDN));
    assert(pkdn != NULL);
    in.iCell = ROOT;
    in.nCell = nCell;
    in.bTreeSqueeze = (msr->nActive > (uint64_t)floor(((double)msr->N)*msr->param.dFracNoTreeSqueeze));
    in.bExcludeVeryActive = bExcludeVeryActive;
    in.dTimeStamp = dTimeStamp;
    if (msr->param.bVDetails) {
	double sec,dsec;
	sec = msrTime();
	pstBuildTree(msr->pst,&in,sizeof(in),pkdn,&iDum);
	printf("Done pstBuildTree\n");
	dsec = msrTime() - sec;
	printf("Tree built, Wallclock: %f secs\n\n",dsec);
	}
    else {
	pstBuildTree(msr->pst,&in,sizeof(in),pkdn,&iDum);
	}
    pstDistribCells(msr->pst,pkdn,nCell*(int)sizeof(KDN),NULL,NULL);
    free(pkdn);
    if (!bExcludeVeryActive && bNeedEwald) {
	/*
	** For simplicity we will skip calculating the Root for all particles 
	** with exclude very active since there are missing particles which 
	** could add to the mass and because it probably is not important to 
	** update the root so frequently.
	*/
	pstCalcRoot(msr->pst,NULL,0,&root,&iDum);
	pstDistribRoot(msr->pst,&root,sizeof(struct ioCalcRoot),NULL,NULL);
	}
    }

void msrBuildTree(MSR msr,double dMass,double dTime,int bNeedEwald) {
    const int bExcludeVeryActive = 0;
    _BuildTree(msr,dMass,dTime,bExcludeVeryActive,bNeedEwald);
    }

void msrBuildTreeExcludeVeryActive(MSR msr,double dMass,double dTime) {
    const int bNeedEwald = 0;
    const int bExcludeVeryActive = 1;
    _BuildTree(msr,dMass,dTime,bExcludeVeryActive,bNeedEwald);
    }


#ifdef GASOLINE
void msrCalcBoundBall(MSR msr,double fBallFactor)
    {
    struct inCalcBoundBall in;
    BND *pbnd;
    int iDum,nCell;

    in.fBallFactor = fBallFactor;
    nCell = 1<<(1+(int)ceil(log((double)msr->nThreads)/log(2.0)));
    pbnd = malloc(nCell*sizeof(BND));
    assert(pbnd != NULL);
    in.iCell = ROOT;
    in.nCell = nCell;
    pstCalcBoundBall(msr->pst,&in,sizeof(in),pbnd,&iDum);
    pstDistribBoundBall(msr->pst,pbnd,nCell*sizeof(BND),NULL,NULL);
    free(pbnd);
    }
#endif

void msrReorder(MSR msr) {
    struct inDomainOrder in;

    in.iMaxOrder = msrMaxOrder(msr)-1;
    if (msr->param.bVDetails) {
	double sec,dsec;
	printf("Ordering...\n");
	sec = msrTime();
	pstDomainOrder(msr->pst,&in,sizeof(in),NULL,NULL);
	pstLocalOrder(msr->pst,NULL,0,NULL,NULL);
	dsec = msrTime() - sec;
	printf("Order established, Wallclock: %f secs\n\n",dsec);
	}
    else {
	pstDomainOrder(msr->pst,&in,sizeof(in),NULL,NULL);
	pstLocalOrder(msr->pst,NULL,0,NULL,NULL);
	}
    /*
    ** Mark domain decomp as not done.
    */
    msr->iLastRungRT = -1;
    }


void msrOutArray(MSR msr,char *pszFile,int iType)
    {
    struct inOutArray in;
    char achOutFile[PST_FILENAME_SIZE];
    LCL *plcl;
    PST pst0;
    FILE *fp;
    int id,i;
    int inswap;

    pst0 = msr->pst;
    while(pst0->nLeaves > 1)
	pst0 = pst0->pstLower;
    plcl = pst0->plcl;
    if (pszFile) {
	/*
	** Add Data Subpath for local and non-local names.
	*/
	_msrMakePath(msr->param.achDataSubPath,pszFile,in.achOutFile);
	/*
	** Add local Data Path.
	*/
	_msrMakePath(plcl->pszDataPath,in.achOutFile,achOutFile);
	fp = fopen(achOutFile,"w");

	if (!fp) {
	    printf("Could not open Array Output File:%s\n",achOutFile);
	    _msrExit(msr,1);
	    }
	}
    else {
	printf("No Array Output File specified\n");
	_msrExit(msr,1);
	return;
	}
    /*
    ** Write the Header information and close the file again.
    */
    fprintf(fp,"%"PRIu64"\n",msr->N);
    fclose(fp);
    /* 
     * First write our own particles.
     */
    assert(msr->pMap[0] == 0);
    pkdOutArray(plcl->pkd,achOutFile,iType); 
    for (i=1;i<msr->nThreads;++i) {
	id = msr->pMap[i];
	/* 
	 * Swap particles with the remote processor.
	 */
	inswap = 0;
	mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
	pkdSwapAll(plcl->pkd, id);
	mdlGetReply(pst0->mdl,id,NULL,NULL);
	/* 
	 * Write the swapped particles.
	 */
	pkdOutArray(plcl->pkd,achOutFile,iType); 
	/* 
	 * Swap them back again.
	 */
	inswap = 0;
	mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
	pkdSwapAll(plcl->pkd, id);
	mdlGetReply(pst0->mdl,id,NULL,NULL);
	}
    }


void msrOutVector(MSR msr,char *pszFile,int iType)
    {
    struct inOutVector in;
    char achOutFile[PST_FILENAME_SIZE];
    LCL *plcl;
    PST pst0;
    FILE *fp;
    int id,i;
    int inswap;
    int iDim;

    pst0 = msr->pst;
    while(pst0->nLeaves > 1)
	pst0 = pst0->pstLower;
    plcl = pst0->plcl;
    if (pszFile) {
	/*
	** Add Data Subpath for local and non-local names.
	*/
	_msrMakePath(msr->param.achDataSubPath,pszFile,in.achOutFile);
	/*
	** Add local Data Path.
	*/
	_msrMakePath(plcl->pszDataPath,in.achOutFile,achOutFile);
	fp = fopen(achOutFile,"w");

	if (!fp) {
	    printf("Could not open Vector Output File:%s\n",achOutFile);
	    _msrExit(msr,1);
	    }
	}
    else {
	printf("No Vector Output File specified\n");
	_msrExit(msr,1);
	return;
	}
    /*
    ** Write the Header information and close the file again.
    */
    fprintf(fp,"%"PRIu64"\n",msr->N);
    fclose(fp);

    /* 
     * First write our own particles.
     */
    assert(msr->pMap[0] == 0);
    for (iDim=0;iDim<3;++iDim) {
	pkdOutVector(plcl->pkd,achOutFile,iDim,iType); 
	for (i=1;i<msr->nThreads;++i) {
	    id = msr->pMap[i];
	    /* 
	     * Swap particles with the remote processor.
	     */
	    inswap = 0;
	    mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
	    pkdSwapAll(plcl->pkd, id);
	    mdlGetReply(pst0->mdl,id,NULL,NULL);
	    /* 
	     * Write the swapped particles.
	     */
	    pkdOutVector(plcl->pkd,achOutFile,iDim,iType); 
	    /* 
	     * Swap them back again.
	     */
	    inswap = 0;
	    mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
	    pkdSwapAll(plcl->pkd,id);
	    mdlGetReply(pst0->mdl,id,NULL,NULL);
	    }
	}
    }


void msrSmooth(MSR msr,double dTime,int iSmoothType,int bGasOnly,
	       int bSymmetric, int eParticleTypes)
    {
    struct inSmooth in;

    /*
    ** Make sure that the type of tree is a density binary tree!
    */
    in.nSmooth = msr->param.nSmooth;
    in.bGasOnly = bGasOnly;
    in.bPeriodic = msr->param.bPeriodic;
    in.bSymmetric = bSymmetric;
    in.iSmoothType = iSmoothType;
    in.eParticleTypes = eParticleTypes;
#ifdef SYMBA
    in.smf.dSunMass = msr->dSunMass;
#endif
    if (msrComove(msr)) {
	in.smf.H = csmTime2Hub(msr->param.csm,dTime);
	in.smf.a = csmTime2Exp(msr->param.csm,dTime);
	}
    else {
	in.smf.H = 0.0;
	in.smf.a = 1.0;
	}
    if (msr->param.bVStep) {
	double sec,dsec;
	printf("Smoothing...\n");
	sec = msrTime();
	pstSmooth(msr->pst,&in,sizeof(in),NULL,NULL);
	dsec = msrTime() - sec;
	printf("Smooth Calculated, Wallclock: %f secs\n\n",dsec);
	}
    else {
	pstSmooth(msr->pst,&in,sizeof(in),NULL,NULL);
	}
    }


void msrReSmooth(MSR msr,double dTime,int iSmoothType,int bGasOnly,
		 int bSymmetric,int eParticleTypes)
    {
    struct inReSmooth in;

    /*
    ** Make sure that the type of tree is a density binary tree!
    */
    in.nSmooth = msr->param.nSmooth;
    in.bGasOnly = bGasOnly;
    in.bPeriodic = msr->param.bPeriodic;
    in.bSymmetric = bSymmetric;
    in.iSmoothType = iSmoothType;
    in.eParticleTypes = eParticleTypes;
    if (msrComove(msr)) {
	in.smf.H = csmTime2Hub(msr->param.csm,dTime);
	in.smf.a = csmTime2Exp(msr->param.csm,dTime);
	}
    else {
	in.smf.H = 0.0;
	in.smf.a = 1.0;
	}
    if (msr->param.bVStep) {
	double sec,dsec;
	printf("ReSmoothing...\n");
	sec = msrTime();
	pstReSmooth(msr->pst,&in,sizeof(in),NULL,NULL);
	dsec = msrTime() - sec;
	printf("ReSmooth Calculated, Wallclock: %f secs\n\n",dsec);
	}
    else {
	pstReSmooth(msr->pst,&in,sizeof(in),NULL,NULL);
	}
    }

void msrUpdateSoft(MSR msr,double dTime) {
#ifdef CHANGESOFT
    if (!(msr->param.bPhysicalSoft || msr->param.bVariableSoft)) return;
    if (msr->param.bPhysicalSoft) {
	struct inPhysicalSoft in;

	in.dFac = 1./csmTime2Exp(msr->param.csm,dTime);
	in.bSoftMaxMul = msr->param.bSoftMaxMul;
	in.dSoftMax = msr->param.dSoftMax;

	if (msr->param.bSoftMaxMul && in.dFac > in.dSoftMax) in.dFac = in.dSoftMax;

	pstPhysicalSoft(msr->pst,&in,sizeof(in),NULL,NULL);
	}
    else {
	struct inPostVariableSoft inPost;
	int bSymmetric = 0;
	int bGasOnly;

	pstPreVariableSoft(msr->pst,NULL,0,NULL,NULL);

	if (msr->param.bSoftByType) {
	    if (msr->nDark) {
		msrBuildTree(msr,-1.0,dTime,0);
		bGasOnly = 0;
		assert(0); /* can't do this yet! */
		msrSmooth(msr,dTime,SMX_NULL,bGasOnly,bSymmetric,TYPE_DARK);
		}
	    if (msr->nGas) {
		msrBuildTree(msr,-1.0,dTime,0);
		bGasOnly = 1;
		msrSmooth(msr,dTime,SMX_NULL,bGasOnly,bSymmetric,TYPE_GAS);
		}
	    if (msr->nStar) {
		msrBuildTree(msr,-1.0,dTime,0);
		bGasOnly = 0;
		assert(0); /* can't do this yet! */
		msrSmooth(msr,dTime,SMX_NULL,bGasOnly,bSymmetric,TYPE_STAR);
		}
	    }
	else {
	    msrBuildTree(msr,-1.0,dTime,0);
	    bGasOnly = 0;
	    msrSmooth(msr,dTime,SMX_NULL,bGasOnly,bSymmetric,TYPE_ALL);
	    }

	inPost.dSoftMax = msr->param.dSoftMax;
	inPost.bSoftMaxMul = msr->param.bSoftMaxMul;
	pstPostVariableSoft(msr->pst,&inPost,sizeof(inPost),NULL,NULL);
	}
#endif
    }

#define PRINTGRID(FRM,VAR) {\
    printf("      % 8d % 8d % 8d % 8d % 8d % 8d % 8d % 8d % 8d % 8d\n",0,1,2,3,4,5,6,7,8,9);\
    for (i=0;i<msr->nThreads/10;++i) {\
	printf("%4d: "FRM" "FRM" "FRM" "FRM" "FRM" "FRM" "FRM" "FRM" "FRM" "FRM"\n",i*10,\
	       out[i*10+0].VAR,out[i*10+1].VAR,out[i*10+2].VAR,out[i*10+3].VAR,out[i*10+4].VAR,\
	       out[i*10+5].VAR,out[i*10+6].VAR,out[i*10+7].VAR,out[i*10+8].VAR,out[i*10+9].VAR);\
	}\
    switch (msr->nThreads%10) {\
    case 0: break;\
    case 1: printf("%4d: "FRM"\n",i*10,\
		   out[i*10+0].VAR); break;\
    case 2: printf("%4d: "FRM" "FRM"\n",i*10,\
		   out[i*10+0].VAR,out[i*10+1].VAR); break;\
    case 3: printf("%4d: "FRM" "FRM" "FRM"\n",i*10,\
		   out[i*10+0].VAR,out[i*10+1].VAR,out[i*10+2].VAR); break;\
    case 4: printf("%4d: "FRM" "FRM" "FRM" "FRM"\n",i*10,\
		   out[i*10+0].VAR,out[i*10+1].VAR,out[i*10+2].VAR,out[i*10+3].VAR); break;\
    case 5: printf("%4d: "FRM" "FRM" "FRM" "FRM" "FRM"\n",i*10,\
		   out[i*10+0].VAR,out[i*10+1].VAR,out[i*10+2].VAR,out[i*10+3].VAR,out[i*10+4].VAR); break;\
    case 6: printf("%4d: "FRM" "FRM" "FRM" "FRM" "FRM" "FRM"\n",i*10,\
		   out[i*10+0].VAR,out[i*10+1].VAR,out[i*10+2].VAR,out[i*10+3].VAR,out[i*10+4].VAR,\
		   out[i*10+5].VAR); break;\
    case 7: printf("%4d: "FRM" "FRM" "FRM" "FRM" "FRM" "FRM" "FRM"\n",i*10,\
		   out[i*10+0].VAR,out[i*10+1].VAR,out[i*10+2].VAR,out[i*10+3].VAR,out[i*10+4].VAR,\
		   out[i*10+5].VAR,out[i*10+6].VAR); break;\
    case 8: printf("%4d: "FRM" "FRM" "FRM" "FRM" "FRM" "FRM" "FRM" "FRM"\n",i*10,\
		   out[i*10+0].VAR,out[i*10+1].VAR,out[i*10+2].VAR,out[i*10+3].VAR,out[i*10+4].VAR,\
		   out[i*10+5].VAR,out[i*10+6].VAR,out[i*10+7].VAR); break;\
    case 9: printf("%4d: "FRM" "FRM" "FRM" "FRM" "FRM" "FRM" "FRM" "FRM" "FRM"\n",i*10,\
		   out[i*10+0].VAR,out[i*10+1].VAR,out[i*10+2].VAR,out[i*10+3].VAR,out[i*10+4].VAR,\
		   out[i*10+5].VAR,out[i*10+6].VAR,out[i*10+7].VAR,out[i*10+8].VAR); break;\
    }\
}

void msrGravity(MSR msr,double dTime,double dStep,int bEwald,
		int *piSec,uint64_t *pnActive) {
    struct inGravity in;
    struct outGravity *out;
    int i,id,iDum;
    double sec,dsec,dTotFlop;

    if (msr->param.bVStep) printf("Calculating Gravity, Step:%f\n",dStep);
    in.dTime = dTime;
    in.nReps = msr->param.nReplicas;
    in.bPeriodic = msr->param.bPeriodic;
    in.bEwald = bEwald;
    in.dEwCut = msr->param.dEwCut;
    in.dEwhCut = msr->param.dEwhCut;

    out = malloc(msr->nThreads*sizeof(struct outGravity));
    assert(out != NULL);
    
    sec = msrTime();
    pstGravity(msr->pst,&in,sizeof(in),out,&iDum);
    dsec = msrTime() - sec;

    *piSec = dsec;
    *pnActive = 0;
    for (id=0;id<msr->nThreads;++id) {
	*pnActive += out[id].nActive;
	}

    if (msr->param.bVStep) {
	/*
	** Output some info...
	*/
	dTotFlop = 0.0;
	for (id=0;id<msr->nThreads;++id) {
	    dTotFlop += out[id].dFlop;
	    if (out[id].nActive > 0) {
		out[id].dPartSum /= out[id].nActive;
		out[id].dCellSum /= out[id].nActive;
	    }
	}

	if(dsec > 0.0) {
	    double dGFlops = dTotFlop/dsec*1e-9;
	    printf("Gravity Calculated, Wallclock: %f secs, GFlops:%.1f, Flop:%.3g\n",
		   dsec,dGFlops,dTotFlop);
	}
	else {
	    printf("Gravity Calculated, Wallclock: %f secs, GFlops:unknown, Flop:%.3g\n",
		   dsec,dTotFlop);
	}
	/*
	** Now comes the really verbose output for each processor.
	*/
	printf("Walk Timings:\n");
	PRINTGRID("% 8.2f",dWalkTime);
#ifdef INSTRUMENT
	printf("Compute Percentage:\n");
	PRINTGRID("% 8.2f",dComputing);
	printf("Cache Wait Percentage:\n");
	PRINTGRID("% 8.2f",dWaiting);
	printf("Load Imbalance Percentage:\n");
	PRINTGRID("% 8.2f",dSynchronizing);
#endif
	printf("Number of Active:\n");
	PRINTGRID("% 8d",nActive);
	printf("Average Number of P-P per Active Particle:\n");
	PRINTGRID("% 8.1f",dPartSum);
	printf("Average Number of P-C per Active Particle:\n");
	PRINTGRID("% 8.1f",dCellSum);
    }
    free(out);
}


void msrCalcEandL(MSR msr,int bFirst,double dTime,double *E,double *T,
		  double *U,double *Eth,double L[])
    {
    struct outCalcEandL out;
    double a;
    int k;

    pstCalcEandL(msr->pst,NULL,0,&out,NULL);
    *T = out.T;
    *U = out.U;
    *Eth = out.Eth;
    for (k=0;k<3;k++) L[k] = out.L[k];
    /*
    ** Do the comoving coordinates stuff.
    ** Currently L is not adjusted for this. Should it be?
    */
    a = csmTime2Exp(msr->param.csm,dTime);
    if (!msr->param.bCannonical) *T *= pow(a,4.0);
    /*
     * Estimate integral (\dot a*U*dt) over the interval.
     * Note that this is equal to integral (W*da) and the latter
     * is more accurate when a is changing rapidly.
     */
    if (msr->param.csm->bComove && !bFirst) {
	msr->dEcosmo += 0.5*(a - csmTime2Exp(msr->param.csm, msr->dTimeOld))
	    *((*U) + msr->dUOld);
	}
    else {
	msr->dEcosmo = 0.0;
	}
    msr->dTimeOld = dTime;
    msr->dUOld = *U;
    *U *= a;
    *E = (*T) + (*U) - msr->dEcosmo + a*a*(*Eth);
    }


void msrDrift(MSR msr,double dTime,double dDelta)
    {
    struct inDrift in;
    int j;

    if (msr->param.bCannonical) {
	in.dDelta = csmComoveDriftFac(msr->param.csm,dTime,dDelta);
	}
    else {
	in.dDelta = dDelta;
	}
    for (j=0;j<3;++j) {
	in.fCenter[j] = msr->fCenter[j];
	}
    in.bPeriodic = msr->param.bPeriodic;
    in.fCentMass = msr->param.dCentMass;
    in.dTime = dTime;
    pstDrift(msr->pst,&in,sizeof(in),NULL,NULL);
    }

void msrDriftInactive(MSR msr,double dTime,double dDelta)
    {
    struct inDrift in;
    int j;

    if (msr->param.bCannonical) {
	in.dDelta = csmComoveDriftFac(msr->param.csm,dTime,dDelta);
	}
    else {
	in.dDelta = dDelta;
	}
    for (j=0;j<3;++j) {
	in.fCenter[j] = msr->fCenter[j];
	}
    in.dTime = dTime;
    in.bPeriodic = msr->param.bPeriodic;
    in.fCentMass = msr->param.dCentMass;
    pstDriftInactive(msr->pst,&in,sizeof(in),NULL,NULL);
    }

/*
 * For gasoline, updates predicted velocities to beginning of timestep.
 */
void msrKickKDKOpen(MSR msr,double dTime,double dDelta)
    {
    double H,a;
    struct inKick in;
    struct outKick out;
	
    if (msr->param.bCannonical) {
	in.dvFacOne = 1.0;		/* no hubble drag, man! */
	in.dvFacTwo = csmComoveKickFac(msr->param.csm,dTime,dDelta);
	}
    else {
	/*
	** Careful! For non-cannonical we want H and a at the 
	** HALF-STEP! This is a bit messy but has to be special
	** cased in some way.
	*/
	dTime += dDelta/2.0;
	a = csmTime2Exp(msr->param.csm,dTime);
	H = csmTime2Hub(msr->param.csm,dTime);
	in.dvFacOne = (1.0 - H*dDelta)/(1.0 + H*dDelta);
	in.dvFacTwo = dDelta/pow(a,3.0)/(1.0 + H*dDelta);
	}
    if (msr->param.bAntiGrav) {
	in.dvFacTwo = -in.dvFacTwo;
	in.dvPredFacTwo = -in.dvPredFacTwo;
	}
    pstKick(msr->pst,&in,sizeof(in),&out,NULL);
    if (msr->param.bVDetails)
	printf("KickOpen: Avg Wallclock %f, Max Wallclock %f\n",
	       out.SumTime/out.nSum,out.MaxTime);
    }

/*
 * For gasoline, updates predicted velocities to end of timestep.
 */
void msrKickKDKClose(MSR msr,double dTime,double dDelta)
    {
    double H,a;
    struct inKick in;
    struct outKick out;
	
    if (msr->param.bCannonical) {
	in.dvFacOne = 1.0; /* no hubble drag, man! */
	in.dvFacTwo = csmComoveKickFac(msr->param.csm,dTime,dDelta);
	}
    else {
	/*
	** Careful! For non-cannonical we want H and a at the 
	** HALF-STEP! This is a bit messy but has to be special
	** cased in some way.
	*/
	dTime += dDelta/2.0;
	a = csmTime2Exp(msr->param.csm,dTime);
	H = csmTime2Hub(msr->param.csm,dTime);
	in.dvFacOne = (1.0 - H*dDelta)/(1.0 + H*dDelta);
	in.dvFacTwo = dDelta/pow(a,3.0)/(1.0 + H*dDelta);
	}
    if (msr->param.bAntiGrav) {
	in.dvFacTwo = -in.dvFacTwo;
	in.dvPredFacTwo = -in.dvPredFacTwo;
	}
    pstKick(msr->pst,&in,sizeof(in),&out,NULL);
    if (msr->param.bVDetails)
	printf("KickClose: Avg Wallclock %f, Max Wallclock %f\n",
	       out.SumTime/out.nSum,out.MaxTime);
    }


int msrOutTime(MSR msr,double dTime)
    {	
    if (msr->iOut < msr->nOuts) {
	if (dTime >= msr->pdOutTime[msr->iOut]) {
	    ++msr->iOut;
	    return(1);
	    }
	else return(0);
	}
    else return(0);
    }


int cmpTime(const void *v1,const void *v2) 
    {
    double *d1 = (double *)v1;
    double *d2 = (double *)v2;

    if (*d1 < *d2) return(-1);
    else if (*d1 == *d2) return(0);
    else return(1);
    }

void msrReadOuts(MSR msr,double dTime)
    {
    char achFile[PST_FILENAME_SIZE];
    char ach[PST_FILENAME_SIZE];
    LCL *plcl = &msr->lcl;
    FILE *fp;
    int i,ret;
    double z,a,n;
    char achIn[80];
	
    /*
    ** Add Data Subpath for local and non-local names.
    */
    achFile[0] = 0;
    sprintf(achFile,"%s/%s.red",msr->param.achDataSubPath,
	    msr->param.achOutName);
    /*
    ** Add local Data Path.
    */
    if (plcl->pszDataPath) {
	strcpy(ach,achFile);
	sprintf(achFile,"%s/%s",plcl->pszDataPath,ach);
	}
    fp = fopen(achFile,"r");
    if (!fp) {
	msr->nOuts = 0;
	return;
	}
    i = 0;
    while (1) {
	if (!fgets(achIn,80,fp)) goto NoMoreOuts;
	switch (achIn[0]) {
	case 'z':
	    ret = sscanf(&achIn[1],"%lf",&z);
	    if (ret != 1) goto NoMoreOuts;
	    a = 1.0/(z+1.0);
	    msr->pdOutTime[i] = csmExp2Time(msr->param.csm,a);
	    break;
	case 'a':
	    ret = sscanf(&achIn[1],"%lf",&a);
	    if (ret != 1) goto NoMoreOuts;
	    msr->pdOutTime[i] = csmExp2Time(msr->param.csm,a);
	    break;
	case 't':
	    ret = sscanf(&achIn[1],"%lf",&msr->pdOutTime[i]);
	    if (ret != 1) goto NoMoreOuts;
	    break;
	case 'n':
	    ret = sscanf(&achIn[1],"%lf",&n);
	    if (ret != 1) goto NoMoreOuts;
	    msr->pdOutTime[i] = dTime + (n-0.5)*msrDelta(msr);
	    break;
	default:
	    ret = sscanf(achIn,"%lf",&z);
	    if (ret != 1) goto NoMoreOuts;
	    a = 1.0/(z+1.0);
	    msr->pdOutTime[i] = csmExp2Time(msr->param.csm,a);
	    }
	++i;
	if(i > msr->nMaxOuts) {
	    msr->nMaxOuts *= 2;
	    msr->pdOutTime = realloc(msr->pdOutTime,
				     msr->nMaxOuts*sizeof(double));
	    assert(msr->pdOutTime != NULL);
	    }
	}
    NoMoreOuts:
    msr->nOuts = i;
    /*
    ** Now sort the array of output times into ascending order.
    */
    qsort(msr->pdOutTime,msr->nOuts,sizeof(double),cmpTime);
    fclose(fp);
    }


int msrSteps(MSR msr)
    {
    return(msr->param.nSteps);
    }

char *msrOutName(MSR msr)
    {
    return(msr->param.achOutName);
    }

char *_BuildName(MSR msr,char *achFile,int iStep,char *defaultPath)
{
    char achOutPath[256], *p;
    int n;

    if ( defaultPath[0] ) {
	strcpy( achOutPath, defaultPath );
	p = strstr( achOutPath, "&N" );
	if ( p ) {
	    n = p - achOutPath;
	    strcpy( p, msrOutName(msr) );
	    strcat( p+2, defaultPath + n + 2 );
	}
	else {
	    n = strlen(achOutPath);
	    if ( !n || achOutPath[n-1]!='/' )
		achOutPath[n++] = '/';
	    strcpy(achOutPath+n,msrOutName(msr));
	}
    }
    else {
	strcpy(achOutPath,msrOutName(msr));
    }

    p = strstr( achOutPath, "&S" );
    if ( p ) {
	n = p - achOutPath;
	strncpy( achFile, achOutPath, n );
	achFile += n;
	sprintf( achFile, "%05d", iStep );
	strcat( achFile, p+2 );
    }
    else {
	sprintf(achFile,msr->param.achDigitMask,msrOutName(msr),iStep);
    }
    return achFile;
}

char *msrBuildName(MSR msr,char *achFile,int iStep)
{
    return _BuildName(msr,achFile,iStep, msr->param.achOutPath);
}

char *msrBuildIoName(MSR msr,char *achFile,int iStep)
{
    if ( msr->param.achIoPath[0] )
	return _BuildName(msr,achFile,iStep, msr->param.achIoPath);
    else
	return msrBuildName(msr,achFile,iStep);
}



double msrDelta(MSR msr)
    {
    return(msr->param.dDelta);
    }


int msrLogInterval(MSR msr)
    {
    return(msr->param.iLogInterval);
    }


int msrOutInterval(MSR msr)
    {
    return(msr->param.iOutInterval);
    }


const char *msrOutTypes(MSR msr)
    {
    return(msr->param.achOutTypes);
    }


int msrCheckInterval(MSR msr)
    {
    return(msr->param.iCheckInterval);
    }


const char *msrCheckTypes(MSR msr)
    {
    return(msr->param.achCheckTypes);
    }


int msrComove(MSR msr)
    {
    return(msr->param.csm->bComove);
    }


double msrSoft(MSR msr)
    {
    return(msr->param.dSoft);
    }


void msrSwitchTheta(MSR msr,double dTime)
    {
    double a;

    a = csmTime2Exp(msr->param.csm,dTime);
    if (a >= msr->param.daSwitchTheta) msr->dCrit = msr->param.dTheta2; 
    }


void
msrInitStep(MSR msr) {
    struct inSetRung insr;
    struct inInitStep in;

    /*
    ** Here we can pass down all parameters of the simulation
    ** before any timestepping takes place. This should happen 
    ** just after the file has been read and the PKD structure
    ** initialized for each processor.
    */
    in.param = msr->param;
    in.csm = *msr->param.csm;
    pstInitStep(msr->pst, &in, sizeof(in), NULL, NULL);
    
    /*
    ** Initialize particles to lowest rung. (what for?)
    */
    insr.iRung = msr->param.iMaxRung - 1;
    pstSetRung(msr->pst, &insr, sizeof(insr), NULL, NULL);
    msr->iCurrMaxRung = insr.iRung;
    }


void
msrSetRung(MSR msr, int iRung)
    {
    struct inSetRung in;

    in.iRung = iRung;
    pstSetRung(msr->pst, &in, sizeof(in), NULL, NULL);
    msr->iCurrMaxRung = in.iRung;
    }


int msrMaxRung(MSR msr)
    {
    return msr->param.iMaxRung;
    }


int msrCurrMaxRung(MSR msr)
    {
    return msr->iCurrMaxRung;
    }


double msrEta(MSR msr)
    {
    return msr->param.dEta;
    }

/*
 * bGreater = 1 => activate all particles at this rung and greater.
 */
void msrActiveRung(MSR msr, int iRung, int bGreater) {
    struct inActiveRung in;

    in.iRung = iRung;
    in.bGreater = bGreater;
    pstActiveRung(msr->pst, &in, sizeof(in), NULL, NULL);

    if ( iRung==0 && bGreater )
	msr->nActive = msr->N;
    else {
	int i;

	assert( msr->nRung != NULL );

	msr->nActive = 0;
	for( i=iRung; i<= (bGreater?msr->param.iMaxRung:iRung); i++ )
	    msr->nActive += msr->nRung[i];
    }
}

void msrActiveOrder(MSR msr)
    {
    pstActiveOrder(msr->pst,NULL,0,&(msr->nActive),NULL);
    }

void msrSetRungVeryActive(MSR msr, int iRung) 
    {
    struct inSetRung in;

    msr->iRungVeryActive = iRung;

    in.iRung = iRung;
    pstSetRungVeryActive(msr->pst,&in,sizeof(in),NULL,NULL);
    }

int msrCurrRung(MSR msr, int iRung)
    {
    struct inCurrRung in;
    struct outCurrRung out;

    in.iRung = iRung;
    pstCurrRung(msr->pst, &in, sizeof(in), &out, NULL);
    return out.iCurrent;
    }

void
msrGravStep(MSR msr,double dTime)
    {
    struct inGravStep in;
    double expand;

    in.dEta = msrEta(msr);
    expand = csmTime2Exp(msr->param.csm,dTime);
    in.dRhoFac = 1.0/(expand*expand*expand);
    pstGravStep(msr->pst,&in,sizeof(in),NULL,NULL);
    }

void
msrAccelStep(MSR msr,double dTime)
    {
    struct inAccelStep in;
    double a;

    in.dEta = msrEta(msr);
    a = csmTime2Exp(msr->param.csm,dTime);
    if (msr->param.bCannonical) {
	in.dVelFac = 1.0/(a*a);
	}
    else {
	in.dVelFac = 1.0;
	}
    in.dAccFac = 1.0/(a*a*a);
    in.bDoGravity = msrDoGravity(msr);
    in.bEpsAcc = msr->param.bEpsAccStep;
    in.bSqrtPhi = msr->param.bSqrtPhiStep;
    pstAccelStep(msr->pst,&in,sizeof(in),NULL,NULL);
    }

void
msrDensityStep(MSR msr,double dTime,int eParticleTypes)
    {
    struct inDensityStep in;
    double expand;
    int bGasOnly,bSymmetric;

    if (msr->param.bVDetails) printf("Calculating Rung Densities...\n");
    bGasOnly = 0;
    bSymmetric = 0;
    msrSmooth(msr,dTime,SMX_DENSITY,bGasOnly,bSymmetric,eParticleTypes);
    in.dEta = msrEta(msr);
    expand = csmTime2Exp(msr->param.csm,dTime);
    in.dRhoFac = 1.0/(expand*expand*expand);
    pstDensityStep(msr->pst,&in,sizeof(in),NULL,NULL);
    }

void
msrInitDt(MSR msr)
    {
    struct inInitDt in;
    
    in.dDelta = msrDelta(msr);
    pstInitDt(msr->pst,&in,sizeof(in),NULL,NULL);
    }

/*
** Returns the Very Active rung based on the number of very active particles desired,
** or the fixed rung that was specified in the parameters.
*/
int msrDtToRung(MSR msr, int iRung, double dDelta, int bAll) {
    struct inDtToRung in;
    struct outDtToRung out;
    int iTempRung,iOutMaxRung,iRungVeryActive;
    uint64_t sum;
    char c;

    in.iRung = iRung;
    in.dDelta = dDelta;
    in.iMaxRung = msrMaxRung(msr);
    in.bAll = bAll;

    pstDtToRung(msr->pst, &in, sizeof(in), &out, NULL);

    iTempRung =msrMaxRung(msr)-1;
    while (out.nRungCount[iTempRung] == 0 && iTempRung > 0) --iTempRung;
    iOutMaxRung = iTempRung;

    while (out.nRungCount[iOutMaxRung] <= msr->param.nTruncateRung && iOutMaxRung > iRung) {
	if (msr->param.bVDetails) printf("n_CurrMaxRung = %"PRIu64"  (iCurrMaxRung = %d):  Promoting particles to iCurrMaxrung = %d\n",
					 out.nRungCount[iOutMaxRung],iOutMaxRung,iOutMaxRung-1);

	in.iMaxRung = iOutMaxRung; /* Note this is the forbidden rung so no -1 here */
	pstDtToRung(msr->pst, &in, sizeof(in), &out, NULL);

	iTempRung =msrMaxRung(msr)-1;
	while (out.nRungCount[iTempRung] == 0 && iTempRung > 0) --iTempRung;
	iOutMaxRung = iTempRung;
	}
    /*
    ** Now copy the rung distribution to the msr structure!
    */
    for (iTempRung=0;iTempRung < msrMaxRung(msr);++iTempRung) msr->nRung[iTempRung] = out.nRungCount[iTempRung];
    /*
    ** Now we want to make a suggestion for the current very active rung based on the number of 
    ** particles in the deepest rungs.
    */
    if (msr->param.nPartVeryActive > 0) {
	iRungVeryActive = msrMaxRung(msr);
	sum = 0;
	while (sum < msr->param.nPartVeryActive && iRungVeryActive > 0) {
	    sum += out.nRungCount[--iRungVeryActive];
	    }
	}
    else {
	iRungVeryActive = msr->param.nRungVeryActive;
	}
	
    msr->iCurrMaxRung = iOutMaxRung;

    if (msr->param.bVRungStat) {
    	printf("Rung distribution:\n");
	printf("\n");
	for (iTempRung=0;iTempRung <= msr->iCurrMaxRung;++iTempRung) {
	    if (out.nRungCount[iTempRung] == 0) continue;
	    if (iTempRung > iRungVeryActive) c = 'v';
	    else c = ' ';
	    printf(" %c rung:%d %"PRIu64"\n",c,iTempRung,out.nRungCount[iTempRung]);
	    }
	printf("\n");
	}
    
    /*
     * Set VeryActive particles
     */
    msrSetRungVeryActive(msr, iRungVeryActive);
    return(iRungVeryActive);
    }


void msrTopStepKDK(MSR msr,
		   double dStep,	/* Current step */
		   double dTime,	/* Current time */
		   double dDelta,	/* Time step */
		   int iRung,		/* Rung level */
		   int iKickRung,	/* Gravity on all rungs from iRung
					   to iKickRung */
		   int iRungVeryActive,  /* current setting for iRungVeryActive */
		   /*
		   ** Note that iRungVeryActive is one less than the first rung with VA particles!
		   */
		   int iAdjust,		/* Do an adjust? */
		   double *pdActiveSum,
		   int *piSec)
    {
    double dMass = -1.0;
    uint64_t nActive;
    int bSplitVA;

    if(iAdjust && (iRung < msrMaxRung(msr)-1)) {
	if (msr->param.bVDetails) {
	    printf("%*cAdjust, iRung: %d\n",2*iRung+2,' ',iRung);
	    }
	msrActiveRung(msr, iRung, 1);
	msrInitDt(msr);
	if (msr->param.bGravStep) {
	    msrGravStep(msr,dTime);
	    }
	if (msr->param.bAccelStep) {
	    msrAccelStep(msr,dTime);
	    }
	if (msr->param.bDensityStep) {
	    bSplitVA = 0;
	    msrDomainDecomp(msr,iRung,1,bSplitVA);
	    msrActiveRung(msr,iRung,1);
	    msrBuildTree(msr,dMass,dTime,0);
	    msrDensityStep(msr,dTime,TYPE_ALL);
	    }
	iRungVeryActive = msrDtToRung(msr,iRung,dDelta,1);
	}
    if (msr->param.bVDetails) {
	printf("%*cmsrKickOpen  at iRung: %d 0.5*dDelta: %g\n",
	       2*iRung+2,' ',iRung,0.5*dDelta);
	}
    msrActiveRung(msr,iRung,0); /* EXACT */
    msrKickKDKOpen(msr,dTime,0.5*dDelta);
    if ((msrCurrMaxRung(msr) > iRung) && (iRungVeryActive > iRung)) {
	/*
	** Recurse.
	*/
	msrTopStepKDK(msr,dStep,dTime,0.5*dDelta,iRung+1,iRung+1,iRungVeryActive,0,
		      pdActiveSum,piSec);
	dTime += 0.5*dDelta;
	dStep += 1.0/(2 << iRung);
	msrActiveRung(msr,iRung,0); /* EXACT */
	msrTopStepKDK(msr,dStep,dTime,0.5*dDelta,iRung+1,iKickRung,iRungVeryActive,1,
		      pdActiveSum,piSec);
	}
    else if(msrCurrMaxRung(msr) == iRung) {
	/* This Drifts everybody */
	if (msr->param.bVDetails) {
	    printf("%*cDrift, iRung: %d\n",2*iRung+2,' ',iRung);
	    }
	msrDrift(msr,dTime,dDelta);
	dTime += dDelta;
	dStep += 1.0/(1 << iRung);

	msrActiveRung(msr,iKickRung,1);
	bSplitVA = 0;
	msrDomainDecomp(msr,iKickRung,1,bSplitVA);

	if(msrDoGravity(msr)) {
	    msrActiveRung(msr,iKickRung,1);
	    msrUpdateSoft(msr,dTime);
	    if (msr->param.bVDetails) {
		printf("%*cGravity, iRung: %d to %d\n",2*iRung+2,' ',iKickRung,iRung);
		}
	    msrBuildTree(msr,dMass,dTime,msr->param.bEwald);
	    msrGravity(msr,dTime,dStep,msr->param.bEwald,piSec,&nActive);
	    *pdActiveSum += (double)nActive/msr->N;
	    }

#ifdef PLANETS
	/* Sun's direct and indirect gravity */
	if(msr->param.bHeliocentric) {
	  msrGravSun(msr);
	    }
#endif
	/*
	 * move time back to 1/2 step so that KickClose can integrate
	 * from 1/2 through the timestep to the end.
	 */
	dTime -= 0.5*dDelta;
	}
    else {
	double dDeltaTmp;
	int i;
	
	/*
	 * We have more rungs to go, but we've hit the very active limit.
	 */
	/*
	 * Activate VeryActives
	 */
	msrActiveRung(msr,msr->iRungVeryActive+1,1);

	/*
	 * Drift the non-VeryActive particles forward 1/2 timestep
	 */
	if (msr->param.bVDetails) {
	    printf("%*cInActiveDrift at iRung: %d, 0.5*dDelta: %g\n",
		   2*iRung+2,' ',iRung,0.5*dDelta);
	    }
	msrDriftInactive(msr, dTime, 0.5*dDelta);
	/*
	 * Build a tree out of them for use by the VeryActives
	 */
	if(msrDoGravity(msr)) {
	    msrUpdateSoft(msr,dTime + 0.5*dDelta);
	    /*
	    ** Domain decomposition for parallel exclude very active is going to be 
	    ** placed here shortly.
	    */
	    bSplitVA = 1;
	    msrDomainDecomp(msr,iRung,1,bSplitVA);

	    if (msr->param.bVDetails) {
		printf("%*cBuilding exclude very active tree: iRung: %d\n",
		       2*iRung+2,' ',iRung);
		}
	    msrBuildTreeExcludeVeryActive(msr,dMass,dTime + 0.5*dDelta);
	    }
	/*
	 * Perform timestepping on individual processors.
	 */
	msrStepVeryActiveKDK(msr, dStep, dTime, dDelta, iRung);
	dTime += dDelta;
	dStep += 1.0/(1 << iRung);
	/*
	 * Move Inactives to the end of the step.
	 */
	if (msr->param.bVDetails) {
	    printf("%*cInActiveDrift at iRung: %d, 0.5*dDelta: %g\n",
		   2*iRung+2,' ',iRung,0.5*dDelta);
	    }
	msrActiveRung(msr,msr->iRungVeryActive+1,1);
	/* 
	** The inactives are half time step behind the actives. 
	** Move them a half time step ahead to synchronize everything again.
	*/
	msrDriftInactive(msr, dTime - 0.5*dDelta, 0.5*dDelta);

	/*
	 * Regular Tree gravity
	 */
	msrActiveRung(msr,iKickRung,1);
	bSplitVA = 0;
	msrDomainDecomp(msr,iKickRung,1,bSplitVA);

	if(msrDoGravity(msr)) {
	    msrActiveRung(msr,iKickRung,1);
	    msrUpdateSoft(msr,dTime);
	    if (msr->param.bVDetails) {
		printf("%*cGravity, iRung: %d to %d\n",
		       2*iRung+2,' ',iKickRung,msrCurrMaxRung(msr));
		}
	    msrBuildTree(msr,dMass,dTime,msr->param.bEwald);
	    msrGravity(msr,dTime,dStep,msr->param.bEwald,piSec,&nActive);
	    *pdActiveSum += (double)nActive/msr->N;
	    }

#ifdef PLANETS
	/* Sun's direct and indirect gravity */
	if(msr->param.bHeliocentric) {
	  msrGravSun(msr);
	    }
#endif

	dDeltaTmp = dDelta;
	for(i = msrCurrMaxRung(msr); i > iRung; i--)
	    dDeltaTmp *= 0.5;
	
	for(i = msrCurrMaxRung(msr); i > iRung; i--) { /* close off all
							 the VeryActive Kicks
						      */
	    if (msr->param.bVDetails) {
		printf("%*cVeryActive msrKickClose at iRung: %d, 0.5*dDelta: %g\n",
		       2*iRung+2,' ',i, 0.5*dDeltaTmp);
		}
	    msrActiveRung(msr,i,0); /* EXACT */
	    msrKickKDKClose(msr,dTime - 0.5*dDeltaTmp,0.5*dDeltaTmp);
	    dDeltaTmp *= 2.0;
	    }
	/*
	 * move time back to 1/2 step so that KickClose can integrate
	 * from 1/2 through the timestep to the end.
	 */
	dTime -= 0.5*dDelta;
	}
    
    if (msr->param.bVDetails) {
	printf("%*cKickClose, iRung: %d, 0.5*dDelta: %g\n",
	       2*iRung+2,' ',iRung, 0.5*dDelta);
	}
    msrActiveRung(msr,iRung,0); /* EXACT */
    msrKickKDKClose(msr,dTime,0.5*dDelta);
    }

void
msrStepVeryActiveKDK(MSR msr, double dStep, double dTime, double dDelta,
		     int iRung)
    {
    struct inStepVeryActive in;
    struct outStepVeryActive out;
    
#ifdef PLANETS
    struct inSunIndirect ins;
    struct outSunIndirect outs;

    if(msr->param.bHeliocentric){
      int k;        

      ins.iFlag = 2; /* for inactive particles */ 
      pstSunIndirect(msr->pst,&ins,sizeof(ins),&outs,NULL); 
	for (k=0;k<3;k++){ 
        in.aSunInact[k] = outs.aSun[k];
	in.adSunInact[k] = outs.adSun[k];
	}
      }
#endif

    in.dStep = dStep;
    in.dTime = dTime;
    in.dDelta = dDelta;
    in.iRung = iRung;
    in.diCrit2 = 1/(msr->dCrit*msr->dCrit);   /* could set a stricter opening criterion here */
    in.nMaxRung = msrCurrMaxRung(msr);
#ifdef PLANETS 
    in.dSunMass = msr->dSunMass;
#endif
    /*
     * Start Particle Cache on all nodes (could be done as part of
     * tree build)
     */
    pstROParticleCache(msr->pst, NULL, 0, NULL, NULL);
    
    pstStepVeryActiveKDK(msr->pst, &in, sizeof(in), &out, NULL);
    /*
     * Finish Particle Cache on all nodes
     */
    pstParticleCacheFinish(msr->pst, NULL, 0, NULL, NULL);
    msr->iCurrMaxRung = out.nMaxRung;
    }

#ifdef HERMITE
void msrTopStepHermite(MSR msr,
				   double dStep,	/* Current step */
				   double dTime,	/* Current time */
				   double dDelta,	/* Time step */
				   int iRung,		/* Rung level */
				   int iKickRung,	/* Gravity on all rungs from iRung
									   to iKickRung */
		                   int iRungVeryActive,  /* current setting for iRungVeryActive */
				   int iAdjust,		/* Do an adjust? */
				   double *pdActiveSum,
				   int *piSec)
{
    double dMass = -1.0;
    uint64_t nActive;
    int bSplitVA;
    const int bNeedEwald = 0;

    if(iAdjust && (iRung < msrMaxRung(msr)-1)) {
        if (msr->param.bVDetails) {
	  printf("%*cAdjust, iRung: %d\n",2*iRung+2,' ',iRung);
           }
	msrActiveRung(msr, iRung, 1);
	msrInitDt(msr);
	if (msr->param.bGravStep) {
	    msrGravStep(msr,dTime);
	    }
	if (msr->param.bAarsethStep) {
	    msrAarsethStep(msr);
	    }
	if (msr->param.bAccelStep) {
	    msrAccelStep(msr,dTime);
	    }
	if (msr->param.bDensityStep) {
	    bSplitVA = 0;
	    msrDomainDecomp(msr,iRung,1,bSplitVA);
	    msrActiveRung(msr,iRung,1);
	    msrBuildTree(msr,dMass,dTime,bNeedEwald);
	    msrDensityStep(msr,dTime,TYPE_ALL);
	    }
	iRungVeryActive = msrDtToRung(msr,iRung,dDelta,1);
    }
	
    if ((msrCurrMaxRung(msr) > iRung) && (iRungVeryActive > iRung)) 
      {
		/*
		 ** Recurse.
		 */
		msrTopStepHermite(msr,dStep,dTime,0.5*dDelta,iRung+1,iRung+1,iRungVeryActive,0,
					  pdActiveSum,piSec);
		dTime += 0.5*dDelta;
		dStep += 1.0/(2 << iRung);
		/*msrActiveRung(msr,iRung,0);*/
		msrTopStepHermite(msr,dStep,dTime,0.5*dDelta,iRung+1,iKickRung,iRungVeryActive,1,
					  pdActiveSum,piSec);
		}
    else if(msrCurrMaxRung(msr) == iRung) {
   
                dTime += dDelta;
		dStep += 1.0/(1 << iRung);

         /* 
           This Predicts everybody 
           dt = dTime - dTime0(pkd)     
           x = x0 + v0*dt + 0.5*a0*dt*dt+ ad0*dt*dt*dt/6.0
           v = v0 + a0*dt + 0.5*ad0*dt*dt
          */
		if (msr->param.bVDetails){
		  printf("%*cPredict, iRung: %d\n",2*iRung+2,' ',iRung);                          
		}
                msrActiveRung(msr,0,1);/* activate everybody*/
                msrPredictor(msr,dTime);     
            
		msrActiveRung(msr,iKickRung,1);
		bSplitVA = 0;
		msrDomainDecomp(msr,iKickRung,1,bSplitVA);
	
		if(msrDoGravity(msr)) {
		  msrActiveRung(msr,iKickRung,1);
		  msrUpdateSoft(msr,dTime);
		  if (msr->param.bVDetails) {
		    printf("%*cGravity, iRung: %d to %d\n",2*iRung+2,' ',iKickRung,iRung);
		  }
		  msrBuildTree(msr,dMass,dTime,bNeedEwald);
		  msrGravity(msr,dTime,dStep,0,0,piSec,&nActive);
		  *pdActiveSum += (double)nActive/msr->N;
		}
#ifdef PLANETS
		/* Sun's direct and indirect gravity */
		if(msr->param.bHeliocentric) {
		  msrGravSun(msr);
		}
#endif                               
		 msrActiveRung(msr,iKickRung,1); /*just in case*/
		/* 
		  Corrector step 
		    (see Kokubo and Makino 2004 for the optimal coefficients) 
		     add = -6.D0*(a0-a)/(dt*dt)-2.0*(2.0*ad0+ad)/dt
                    addd = 12.0*(a0-a)/(dt*dt*dt)+6.0*(ad0+ad)/(dt*dt) 
		     x = x + add*dt**4/24.0+addd*dt**5*alpha/120.0 
		     v = v + add*dt**3/6.0+addd*dt**4/24.0 
	        */
		 msrCorrector(msr,dTime);    
 
		 if (msr->param.bVDetails){
		  printf("%*cCorrect, iRung: %d\n",2*iRung+2,' ',iRung);                          
		}
#ifdef PLANETS
                if(msr->param.bHeliocentric){                      		
		/*
		Here is a step for correcting the Sun's direct gravity.
		a is recalculated using the correctors' position and 
		velocity.  
		*/          
		int nite = 0;
                int nitemax = 3; 
                /* number of interation for correcting the Sun's gravity 
                   P(EC)'^(nitemax) scheme (see Kokubo et al 1998)*/           
                do{    
                nite += 1;
                msrSunCorrector(msr,dTime);  		
                }while(nite < nitemax);
                }
		if(msr->param.bCollision){   
		  msrDoCollision(msr,dTime,dDelta);  
		  }
#endif 
		/* 
		Copy the present values of activated particles as the initial values 
		x_0 = x, v_0 = v, a_0 = a, dTime0 = dTime
		*/ 
                msrCopy0(msr,dTime); 
    }

else {        
   	double dDeltaTmp;
	int i;

	/*printf("iRungVeryActive: %d CurrMaxrung: %d  iRung: %d, 0.5*dDelta: %g n/",iRungVeryActive, msrCurrMaxRung(msr),iRung,0.5*dDelta);*/

	/*
	 * We have more rungs to go, but we've hit the very active limit.
	 */
	/*
	 * Activate VeryActives
	 */
	msrActiveRung(msr,msr->iRungVeryActive+1,1);
	/*
	 * Predict the non-VeryActive particles forward 1/2 timestep
	 */
	if (msr->param.bVDetails)
	  	if (msr->param.bVDetails) {
	    printf("%*cInActivePredict at iRung: %d, 0.5*dDelta: %g\n",
		   2*iRung+2,' ',iRung,0.5*dDelta);
	    }
	    
	msrPredictorInactive(msr, dTime + 0.5*dDelta);
	/*
	 * Build a tree out of them for use by the VeryActives
	 */
	if(msrDoGravity(msr)) {
	    msrUpdateSoft(msr,dTime + 0.5*dDelta);
	    /*
	    ** Domain decomposition for parallel exclude very active is going to be 
	    ** placed here shortly.
	    */
	    bSplitVA = 1;
	    msrDomainDecomp(msr,iRung,1,bSplitVA);

	    if (msr->param.bVDetails) {
		printf("%*cBuilding exclude very active tree: iRung: %d\n",
		       2*iRung+2,' ',iRung);
		}
	    msrBuildTreeExcludeVeryActive(msr,dMass,dTime + 0.5*dDelta);
	    }
	/*
	 * Perform timestepping on individual processors.
	 */
	if (msr->param.bVDetails)
	    printf("VeryActive at iRung: %d\n", iRung);
	msrStepVeryActiveHermite(msr, dStep, dTime, dDelta, iRung);
	dTime += dDelta;
	dStep += 1.0/(1 << iRung);
	/*
	 * Move Inactives to the end of the step.
	 */
	if (msr->param.bVDetails) {
	    printf("%*cInActivePredictor at iRung: %d, 0.5*dDelta: %g\n",
		   2*iRung+2,' ',iRung,0.5*dDelta);
	    }
	msrActiveRung(msr,msr->iRungVeryActive+1,1);
	/* 
	** The inactives are half time step behind the actives. 
	** Move them a half time step ahead to synchronize everything again.
	*/	
	msrPredictorInactive(msr, dTime);

	/*
	 * Regular Tree gravity
	 */
	msrActiveRung(msr,iKickRung,1);
	bSplitVA = 0;
	msrDomainDecomp(msr,iKickRung,1,bSplitVA);

	if(msrDoGravity(msr)) {
	    msrActiveRung(msr,iKickRung,1);
	    msrUpdateSoft(msr,dTime);
	    if (msr->param.bVDetails) {
		printf("%*cGravity, iRung: %d to %d\n",
		       2*iRung+2,' ',iKickRung,msrCurrMaxRung(msr));
		}
	    msrBuildTree(msr,dMass,dTime,bNeedEwald);
	    msrGravity(msr,dTime,dStep,0,0,piSec,&nActive);
	    *pdActiveSum += (double)nActive/msr->N;
	    }
#ifdef PLANETS
	/* Sun's direct and indirect gravity */
		if(msr->param.bHeliocentric) {
		  msrGravSun(msr);
		}
#endif 
	if (msr->param.bVDetails) {
		printf("%*cVeryActive msrCorrector at iRung: %d, 0.5*dDelta: %g\n",
		       2*iRung+2,' ',i, 0.5*dDeltaTmp);
		}	
	        msrCorrector(msr,dTime);   

#ifdef PLANETS         
                if(msr->param.bHeliocentric){          		
		int nite = 0;
                int nitemax = 3;                                                     
                do{    
                nite += 1;
                msrSunCorrector(msr,dTime);  		
                }while(nite < nitemax);
                }
		if(msr->param.bCollision){   
		  msrDoCollision(msr,dTime,dDelta);  
		}
#endif
                msrCopy0(msr,dTime);
          }
}

void
msrStepVeryActiveHermite(MSR msr, double dStep, double dTime, double dDelta,
		     int iRung)
    {
    struct inStepVeryActiveH in;
    struct outStepVeryActiveH out;    
 
#ifdef PLANETS
    if(msr->param.bHeliocentric){
      int k;        
      struct inSunIndirect ins;
      struct outSunIndirect outs;

      ins.iFlag = 2; /* for inactive particles */ 
      pstSunIndirect(msr->pst,&ins,sizeof(ins),&outs,NULL); 
	for (k=0;k<3;k++){ 
        in.aSunInact[k] = outs.aSun[k];
	in.adSunInact[k] = outs.adSun[k];
	}
      }
#endif

    in.dStep = dStep;
    in.dTime = dTime;
    in.dDelta = dDelta;
    in.iRung = iRung;
    in.diCrit2 = 1/(msr->dCrit*msr->dCrit);   /* could set a stricter opening criterion here */
    in.nMaxRung = msrCurrMaxRung(msr);
#ifdef PLANETS 
    in.dSunMass = msr->dSunMass;
#endif
    /*
     * Start Particle Cache on all nodes (could be done as part of
     * tree build)
     */
    pstROParticleCache(msr->pst, NULL, 0, NULL, NULL);
    
    pstStepVeryActiveHermite(msr->pst, &in, sizeof(in), &out, NULL);

#ifdef PLANETS 
/* maybe we should use collision flag to determine if we call 
   the following function.  */
    if(msr->param.bCollision){
    struct outGetVariableVeryActive outGet;
    pstGetVariableVeryActive(msr->pst, NULL, 0, &outGet, NULL);
    msr->dEcoll += outGet.dDeltaEcoll;
    outGet.dDeltaEcoll = 0.0;/*just in case */
    }
#endif
    /*
     * Finish Particle Cache on all nodes
     */
    pstParticleCacheFinish(msr->pst, NULL, 0, NULL, NULL);
    msr->iCurrMaxRung = out.nMaxRung;
    }

void msrCopy0(MSR msr,double dTime)
{
	struct inCopy0 in;

	in.dTime = dTime;
	pstCopy0(msr->pst,&in,sizeof(in),NULL,NULL);
	}

void msrPredictor(MSR msr,double dTime)
{
	struct inPredictor in;
	
	in.dTime = dTime;
	pstPredictor(msr->pst,&in,sizeof(in),NULL,NULL);
	}

void msrCorrector(MSR msr,double dTime)
{
	struct inCorrector in;

	in.dTime = dTime;
	pstCorrector(msr->pst,&in,sizeof(in),NULL,NULL);
	}

#ifdef PLANETS 
void msrSunCorrector(MSR msr,double dTime)
{
	struct inSunCorrector in;
	
	in.dTime = dTime;
	in.dSunMass = msr->dSunMass;
	pstSunCorrector(msr->pst,&in,sizeof(in),NULL,NULL);
	}
#endif

void msrPredictorInactive(MSR msr,double dTime)
{
	struct inPredictorInactive in;
	
	in.dTime = dTime;
	pstPredictorInactive(msr->pst,&in,sizeof(in),NULL,NULL);
	}

void
msrAarsethStep(MSR msr)
{
    struct inAarsethStep in;

    in.dEta = msrEta(msr);
    pstAarsethStep(msr->pst,&in,sizeof(in),NULL,NULL);
    }

void
msrFirstDt(MSR msr)
{  
    pstFirstDt(msr->pst,NULL,0,NULL,NULL);
    }

#endif /* Hermite*/

uint64_t msrMaxOrder(MSR msr)
{
    return msr->nMaxOrder;
    }

void
msrAddDelParticles(MSR msr)
    {
    struct outColNParts *pColNParts;
    uint64_t *pNewOrder;
    struct inSetNParts in;
    struct inSetParticleTypes intype;
#ifdef PLANETS
    struct inHandSunMass inh;
#endif
    int iOut;
    int i;
    
    if (msr->param.bVDetails) printf("Changing Particle number\n");
    pColNParts = malloc(msr->nThreads*sizeof(*pColNParts));
    assert(pColNParts!=NULL);
    pstColNParts(msr->pst, NULL, 0, pColNParts, &iOut);
    /*
     * Assign starting numbers for new particles in each processor.
     */
    pNewOrder = malloc(msr->nThreads*sizeof(*pNewOrder));
    assert(pNewOrder!=NULL);
    for(i=0;i<msr->nThreads;i++) {
	/*
	 * Detect any changes in particle number, and force a tree
	 * build.
	 */
     if (pColNParts[i].nNew != 0 || pColNParts[i].nDeltaGas != 0 ||
	    pColNParts[i].nDeltaDark != 0 || pColNParts[i].nDeltaStar != 0) {
	  /*printf("Particle assignments have changed!\n");
	    printf("need to rebuild tree, code in msrAddDelParticles()\n");
	    printf("needs to be updated. Bailing out for now...\n");
	    exit(-1); */	    
	pNewOrder[i] = msr->nMaxOrder;
	msr->nMaxOrder += pColNParts[i].nNew;
	msr->nGas += pColNParts[i].nDeltaGas;
	msr->nDark += pColNParts[i].nDeltaDark;
	msr->nStar += pColNParts[i].nDeltaStar;      
     }
     }
    msr->N = msr->nGas + msr->nDark + msr->nStar;

    msr->nMaxOrderDark = msr->nMaxOrder;

    pstNewOrder(msr->pst,pNewOrder,(int)sizeof(*pNewOrder)*msr->nThreads,NULL,NULL);

    if (msr->param.bVDetails)
	printf("New numbers of particles: %"PRIu64" gas %"PRIu64" dark %"PRIu64" star\n",
	       msr->nGas, msr->nDark, msr->nStar);

    in.nGas = msr->nGas;
    in.nDark = msr->nDark;
    in.nStar = msr->nStar;
    in.nMaxOrderGas = msr->nMaxOrderGas;
    in.nMaxOrderDark = msr->nMaxOrderDark;    
    pstSetNParts(msr->pst,&in,sizeof(in),NULL,NULL);
    pstSetParticleTypes(msr->pst,&intype,sizeof(intype),NULL,NULL);

#ifdef PLANETS
    inh.dSunMass = msr->dSunMass;
    pstHandSunMass(msr->pst,&inh,sizeof(inh),NULL,NULL);
#endif
    free(pNewOrder);
    free(pColNParts);
    }

int msrDoDensity(MSR msr)
    {
    return(msr->param.bDoDensity);
    }

#ifdef USE_PNG
int msrPNGResolution(MSR msr)
    {
    return(msr->param.nPNGResolution);
    }
#endif

int msrDoGravity(MSR msr)
    {
    return(msr->param.bDoGravity);
    }

void msrInitTimeSteps(MSR msr,double dTime,double dDelta) 
    {
    double dMass = -1.0;
    const int bNeedEwald = 0;

    msrActiveRung(msr,0,1); /* Activate all particles */
    msrInitDt(msr);
    if (msr->param.bGravStep) {
	msrGravStep(msr,dTime);
	}
    if (msr->param.bAccelStep) {
	msrAccelStep(msr,dTime);
	}
    if (msr->param.bDensityStep) {
	msrDomainDecomp(msr,0,1,0);
	msrActiveRung(msr,0,1); /* Activate all particles */
	msrBuildTree(msr,dMass,dTime,bNeedEwald);
	msrDensityStep(msr,dTime,TYPE_ALL);
	}
    msrDtToRung(msr,0,dDelta,1);
    }


void msrFof(MSR msr,int nFOFsDone,int iSmoothType,int bSymmetric,int eParticleTypes, double exp)
    {
    struct inFof in;
    in.nFOFsDone = nFOFsDone;
    in.nSmooth = msr->param.nSmooth;
    in.bPeriodic = msr->param.bPeriodic;
    in.bSymmetric = bSymmetric;
    in.iSmoothType = iSmoothType;
    in.eParticleTypes = eParticleTypes;
    in.smf.a = exp;
    in.smf.dTau2 = pow(msr->param.dTau,2.0);
    in.smf.dVTau2 = pow(msr->param.dVTau,2.0);
    if(msr->param.bTauAbs == 0) {
      in.smf.dTau2 *= pow(msr->param.csm->dOmega0,-0.6666);
    } else {
      in.smf.dTau2 /= exp*exp;
      in.smf.dVTau2 *= exp*exp;
    }
    in.smf.bTauAbs = msr->param.bTauAbs;
    in.smf.nMinMembers = msr->param.nMinMembers;
    in.smf.fContrast = msr->param.fContrast;
    in.smf.Delta = msr->param.Delta;
    if (msr->param.bVStep) {
	double sec,dsec;
	printf("Doing FOF with dTau2=%e , dVTau=%e ...\n", in.smf.dTau2, in.smf.dVTau2);
	sec = msrTime();
	pstFof(msr->pst,&in,sizeof(in),NULL,NULL);
	dsec = msrTime() - sec;
	printf("FOF Calculated, Wallclock: %f secs\n\n",dsec);
	}
    else {
	pstFof(msr->pst,&in,sizeof(in),NULL,NULL);
	}
    }

void msrGroupMerge(MSR msr, double exp)
    {
    struct inGroupMerge in;
    int nGroups;
    in.bPeriodic = msr->param.bPeriodic;
    in.smf.nMinMembers = msr->param.nMinMembers;
    in.smf.Delta = msr->param.Delta;
    in.smf.a = exp;
    if (msr->param.bVStep) {
	double sec,dsec;
	printf("Doing GroupMerge...\n");
	sec = msrTime();
	pstGroupMerge(msr->pst,&in,sizeof(in),&nGroups,NULL);
	dsec = msrTime() - sec;
	printf("GroupMerge done, Wallclock: %f secs\n",dsec);
	}
    else {
	pstGroupMerge(msr->pst,&in,sizeof(in),&nGroups,NULL);
	}
    msr->nGroups = nGroups;
    printf("MASTER: TOTAL groups: %i \n" ,nGroups);
    }

void msrGroupProfiles(MSR msr,int nFOFsDone,int iSmoothType,int bSymmetric,int eParticleTypes, double exp)
    {
    int nBins;
    struct inGroupProfiles in;
    in.bPeriodic = msr->param.bPeriodic;
    in.nTotalGroups = msr->nGroups;
    in.nSmooth = msr->param.nSmooth;
    in.bSymmetric = bSymmetric;
    in.iSmoothType = iSmoothType;
    in.eParticleTypes = eParticleTypes;
    in.nFOFsDone = nFOFsDone;
    in.smf.nMinMembers = msr->param.nMinMembers;
    in.smf.nBins = msr->param.nBins;
    in.smf.nMinProfile = msr->param.nMinProfile/pow(msr->param.fBinsRescale, nFOFsDone);	
    in.bLogBins = msr->param.bLogBins;
    in.smf.bUsePotmin = msr->param.bUsePotmin;
    in.smf.Delta = msr->param.Delta;
    in.smf.binFactor = msr->param.binFactor;
    in.smf.a = exp;
    if (msr->param.bVStep) {
	double sec,dsec;
	printf("Doing GroupProfiles...\n");
	sec = msrTime();
	pstGroupProfiles(msr->pst,&in,sizeof(in),&nBins,NULL);
	dsec = msrTime() - sec;
	printf("GroupProfiles done, Wallclock: %f secs\n",dsec);
	}
    else {
	pstGroupProfiles(msr->pst,&in,sizeof(in),&nBins,NULL);
	}
    msr->nBins = nBins;
    printf("MASTER: TOTAL bins: %i TOTAL groups: %i \n" ,nBins,msr->nGroups);
    }

void msrOutGroups(MSR msr,char *pszFile,int iOutType, double dTime){

    struct inOutArray in;
    char achOutFile[PST_FILENAME_SIZE];
    LCL *plcl;
    PST pst0;
    FILE *fp;
    double dvFac,time;
   
    pst0 = msr->pst;
    while(pst0->nLeaves > 1)
	pst0 = pst0->pstLower;
    plcl = pst0->plcl;
    if (pszFile) {
	/*
	** Add Data Subpath for local and non-local names.
	*/
	_msrMakePath(msr->param.achDataSubPath,pszFile,in.achOutFile);
	/*
	** Add local Data Path.
	*/
	_msrMakePath(plcl->pszDataPath,in.achOutFile,achOutFile);

	fp = fopen(achOutFile,"w");
	if (!fp) {
	    printf("Could not open Group Output File:%s\n",achOutFile);
	    _msrExit(msr,1);
	    }
	}
    else {
	printf("No Group Output File specified\n");
	_msrExit(msr,1);
	return;
	}
    if (msrComove(msr)) {
	time = csmTime2Exp(msr->param.csm,dTime);
	if (msr->param.bCannonical) {
	    dvFac = 1.0/(time*time);
	    }
	else {
	    dvFac = 1.0;
	    }
	}
    else {
	time = dTime;
	dvFac = 1.0;
	}

    if(iOutType == OUT_GROUP_TIPSY_NAT || iOutType == OUT_GROUP_TIPSY_STD){
	/*
	** Write tipsy header.
	*/
	struct dump h;
	h.nbodies = msr->nGroups;
	h.ndark = 0;
	h.nsph = 0;
	h.nstar = msr->nGroups;
	h.time = time;
	h.ndim = 3;

	if (msrComove(msr)) {
	    printf("Writing file...\nTime:%g Redshift:%g\n",
		   dTime,(1.0/h.time - 1.0));
	    }
	else {
	    printf("Writing file...\nTime:%g\n",dTime);
	    }
		
	if (iOutType == OUT_GROUP_TIPSY_STD) {
	    XDR xdrs;
	    xdrstdio_create(&xdrs,fp,XDR_ENCODE);
	    xdrHeader(&xdrs,&h);
	    xdr_destroy(&xdrs);
	    }
	else {
	    fwrite(&h,sizeof(struct dump),1,fp);
	    }	
	}	
    fclose(fp);	
    /* 
     * Write the groups.
     */
    assert(msr->pMap[0] == 0);
    pkdOutGroup(plcl->pkd,achOutFile,iOutType,0,dvFac);
    }

void msrDeleteGroups(MSR msr){

    LCL *plcl;
    PST pst0;
   
    pst0 = msr->pst;
    while(pst0->nLeaves > 1)
	pst0 = pst0->pstLower;
    plcl = pst0->plcl;

    if(plcl->pkd->groupData)free(plcl->pkd->groupData); 
    if(plcl->pkd->groupBin)free(plcl->pkd->groupBin);
    plcl->pkd->nBins = 0;
    plcl->pkd->nGroups = 0;
    }
#ifdef RELAXATION
void msrInitRelaxation(MSR msr)
    {
    pstInitRelaxation(msr->pst,NULL,0,NULL,NULL);
    }	
void msrRelaxation(MSR msr,double dTime,double deltaT,int iSmoothType,int bSymmetric,int eParticleTypes)
    {
    struct inSmooth in;
    in.nSmooth = msr->param.nSmooth;
    in.bPeriodic = msr->param.bPeriodic;
    in.bSymmetric = bSymmetric;
    in.iSmoothType = iSmoothType;
    in.eParticleTypes = eParticleTypes;
    in.dfBall2OverSoft2 = (msr->param.bLowerSoundSpeed ? 0 :
			   4.0*msr->param.dhMinOverSoft*msr->param.dhMinOverSoft);
    if (msrComove(msr)) {
	in.smf.H = csmTime2Hub(msr->param.csm,dTime);
	in.smf.a = csmTime2Exp(msr->param.csm,dTime);
	}
    else {
	in.smf.H = 0.0;
	in.smf.a = 1.0;
	}
    in.smf.dDeltaT = deltaT;
    if (msr->param.bVStep) {
	double sec,dsec;
	printf("Smoothing for relaxation...dDeltaT = %f \n",deltaT);
	sec = msrTime();
	pstSmooth(msr->pst,&in,sizeof(in),NULL,NULL);
	dsec = msrTime() - sec;
	printf("Relaxation Calculated, Wallclock: %f secs\n\n",dsec);
	}
    else {
	pstSmooth(msr->pst,&in,sizeof(in),NULL,NULL);
	}
    }
#endif /* RELAXATION */

#ifdef PLANETS 
void
msrOneNodeReadSS(MSR msr,struct inReadSS *in)
{
    int i,id;
    int *nParts;
    int nStart;
    PST pst0;
    LCL *plcl;
    char achInFile[PST_FILENAME_SIZE];
    int nid;
    int inswap;

    nParts = malloc(msr->nThreads*sizeof(*nParts));
    for (id=0;id<msr->nThreads;++id) {
		nParts[id] = -1;
		}

    pstOneNodeReadInit(msr->pst,in,sizeof(*in),nParts,&nid);
    assert(nid == msr->nThreads*sizeof(*nParts));
    for (id=0;id<msr->nThreads;++id) {
		assert(nParts[id] > 0);
		}

    pst0 = msr->pst;
    while(pst0->nLeaves > 1)
		pst0 = pst0->pstLower;
    plcl = pst0->plcl;
    /*
     ** Add the local Data Path to the provided filename.
     */
	_msrMakePath(plcl->pszDataPath,in->achInFile,achInFile);

    nStart = nParts[0];
	assert(msr->pMap[0] == 0);
    for (i=1;i<msr->nThreads;++i) {
		id = msr->pMap[i];
		/* 
		 * Read particles into the local storage.
		 */
		assert(plcl->pkd->nStore >= nParts[id]);
		pkdReadSS(plcl->pkd,achInFile,nStart,nParts[id]);
		nStart += nParts[id];
		/* 
		 * Now shove them over to the remote processor.
		 */
		inswap = 0;
		mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
		pkdSwapAll(plcl->pkd, id);
		mdlGetReply(pst0->mdl,id,NULL,NULL);
    	}
    assert(nStart == msr->N);
    /* 
     * Now read our own particles.
     */
    pkdReadSS(plcl->pkd,achInFile,0,nParts[0]);
    }

double
msrReadSS(MSR msr)
{
	SSIO ssio;
	SSHEAD head;
	struct inReadSS in;
	struct inHandSunMass inh;

	struct inSetParticleTypes intype;
	char achInFile[PST_FILENAME_SIZE];
	LCL *plcl = msr->pst->plcl;
	double dTime;

	if (msr->param.achInFile[0]) {
		/*
		 ** Add Data Subpath for local and non-local names.
		 */
		_msrMakePath(msr->param.achDataSubPath,msr->param.achInFile,in.achInFile);
		/*
		 ** Add local Data Path.
		 */
		_msrMakePath(plcl->pszDataPath,in.achInFile,achInFile);

		if (ssioOpen(achInFile,&ssio,SSIO_READ)) {
			printf("Could not open InFile:%s\n",achInFile);
			_msrExit(msr,1);
			}
		}
	else {
		printf("No input file specified\n");
		_msrExit(msr,1);
		}

	/* Read header */

	if (ssioHead(&ssio,&head)) {
		printf("Could not read header of InFile:%s\n",achInFile);
		_msrExit(msr,1);
		}
	if (ssioClose(&ssio)) {
		printf("Could not close InFile:%s\n",achInFile);
		_msrExit(msr,1);
		}

	msr->N = msr->nDark = head.n_data;
	msr->nGas = msr->nStar = 0;
	msr->nMaxOrder = msr->N;
	msr->nMaxOrderGas = msr->nGas; /* NOW ALWAYS ZERO : was always -1 */
	msr->nMaxOrderDark = msr->nDark;
        msr->nPlanets = head.n_planets;        
	msr->dEcoll = head.dEcoll;        
        msr->dSunMass = head.dSunMass;        


	dTime = head.time;
	if (msr->param.bVStart) {
		double tTo;
		printf("Input file...N=%i,Time=%g\n",msr->N,dTime);
		tTo = dTime + msr->param.nSteps*msr->param.dDelta;
		printf("Simulation to Time:%g\n",tTo);
		}

	in.nFileStart = 0;
	in.nFileEnd = msr->N - 1;
	in.nBucket = msr->param.nBucket;
	in.nDark = msr->nDark;
	in.nGas = msr->nGas;	/* always zero */
	in.nStar = msr->nStar;	/* always zero */
	in.iOrder = msr->param.iOrder;

	/*
	 ** Since pstReadSS causes the allocation of the local particle
	 ** store, we need to tell it the percentage of extra storage it
	 ** should allocate for load balancing differences in the number of
	 ** particles.
	 */
	in.fExtraStore = msr->param.dExtraStore;

	in.fPeriod[0] = msr->param.dxPeriod;
	in.fPeriod[1] = msr->param.dyPeriod;
	in.fPeriod[2] = msr->param.dzPeriod;

	if (msr->param.bParaRead)
	    pstReadSS(msr->pst,&in,sizeof(in),NULL,NULL);
	else
	    msrOneNodeReadSS(msr,&in);
	pstSetParticleTypes(msr->pst,&intype,sizeof(intype),NULL,NULL);
	if (msr->param.bVDetails) puts("Input file successfully read.");

	inh.dSunMass = msr->dSunMass;       
	pstHandSunMass(msr->pst,&inh,sizeof(inh),NULL,NULL);
	

	/*
	 ** Now read in the output points, passing the initial time.
	 ** We do this only if nSteps is not equal to zero.
	 */
	if (msrSteps(msr) > 0) msrReadOuts(msr,dTime);
	/*
	 ** Set up the output counter.
	 */
	for (msr->iOut=0;msr->iOut<msr->nOuts;++msr->iOut) {
		if (dTime < msr->pdOutTime[msr->iOut]) break;
		}
	return(dTime);
	}

void
msrOneNodeWriteSS(MSR msr,struct inWriteSS *in)
{
    int i,id;
    int nStart;
    PST pst0;
    LCL *plcl;
    char achOutFile[PST_FILENAME_SIZE];
    int inswap;

    pst0 = msr->pst;
    while(pst0->nLeaves > 1)
		pst0 = pst0->pstLower;
    plcl = pst0->plcl;
    /*
     ** Add the local Data Path to the provided filename.
     */
	_msrMakePath(plcl->pszDataPath,in->achOutFile,achOutFile);

    /* 
     * First write our own particles.
     */
    pkdWriteSS(plcl->pkd,achOutFile,plcl->nWriteStart);
    nStart = plcl->pkd->nLocal;
	assert(msr->pMap[0] == 0);
    for (i=1;i<msr->nThreads;++i) {
		id = msr->pMap[i];
		/* 
		 * Swap particles with the remote processor.
		 */
		inswap = 0;
		mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
		pkdSwapAll(plcl->pkd,id);
		mdlGetReply(pst0->mdl,id,NULL,NULL);
		/* 
		 * Write the swapped particles.
		 */
		pkdWriteSS(plcl->pkd,achOutFile,nStart);
		nStart += plcl->pkd->nLocal;
		/* 
		 * Swap them back again.
		 */
		inswap = 0;
		mdlReqService(pst0->mdl,id,PST_SWAPALL,&inswap,sizeof(inswap));
		pkdSwapAll(plcl->pkd, id);
		mdlGetReply(pst0->mdl,id,NULL,NULL);
    	}
    assert(nStart == msr->N);
    }

void
msrWriteSS(MSR msr,char *pszFileName,double dTime)
{
	SSIO ssio;
	SSHEAD head;
	struct inWriteSS in;
	char achOutFile[PST_FILENAME_SIZE];
	LCL *plcl = msr->pst->plcl;

	/*
	 ** Calculate where each processor should start writing.
	 ** This sets plcl->nWriteStart.
	 */
	msrCalcWriteStart(msr);
	/*
	 ** Add Data Subpath for local and non-local names.
	 */
	_msrMakePath(msr->param.achDataSubPath,pszFileName,in.achOutFile);
	/*
	 ** Add local Data Path.
	 */
	_msrMakePath(plcl->pszDataPath,in.achOutFile,achOutFile);
	
	if (ssioOpen(achOutFile,&ssio,SSIO_WRITE)) {
		printf("Could not open OutFile:%s\n",achOutFile);
		_msrExit(msr,1);
		}

	/* Write header */

	head.time = dTime;
	head.n_data = msr->N;
	head.n_planets = msr->nPlanets;
	head.dEcoll = msr->dEcoll;
	head.dSunMass = msr->dSunMass;

	if (ssioHead(&ssio,&head)) {
		printf("Could not write header of OutFile:%s\n",achOutFile);
		_msrExit(msr,1);
		}
	if (ssioClose(&ssio)) {
		printf("Could not close OutFile:%s\n",achOutFile);
		_msrExit(msr,1);
		}

	if(msr->param.bParaWrite)
	    pstWriteSS(msr->pst,&in,sizeof(in),NULL,NULL);
	else
	    msrOneNodeWriteSS(msr,&in);

	if (msr->param.bVDetails) puts("Output file successfully written.");
	}

void msrGravSun(MSR msr)
{
	struct inGravSun in;

	struct inSunIndirect ins;
	struct outSunIndirect outs;

	int j;

	/* Calculate Sun's indirect gravity */
	ins.iFlag = 0;	/* for all particles */
	pstSunIndirect(msr->pst,&ins,sizeof(ins),&outs,NULL); 
	for (j=0;j<3;++j){ 
                in.aSun[j] = outs.aSun[j];
		in.adSun[j] = outs.adSun[j];
		}
	/* printf("asun = %e %e %e adsun = %e %e %e \n",in.aSun[0],in.aSun[1],in.aSun[2],in.adSun[0],in.adSun[1],in.adSun[2]); */

	in.dSunMass = msr->dSunMass;

	pstGravSun(msr->pst,&in,sizeof(in),NULL,NULL);

}

static char *
_msrParticleLabel(MSR msr,int iColor)
{
	switch (iColor) {
	case SUN:
		return "SUN";
	case JUPITER:
		return "JUPITER";
	case SATURN:
		return "SATURN";
	case URANUS:
		return "URANUS";
	case NEPTUNE:
		return "NEPTUNE";
	case PLANETESIMAL:
		return "PLANETESIMAL";
	default:
		return "UNKNOWN";
		}
	}

void
msrDoCollision(MSR msr,double dTime,double dDelta)
{

	struct outNextCollision next;
	struct inGetColliderInfo inGet;
	struct outGetColliderInfo outGet;
	struct inDoCollision inDo;
	struct outDoCollision outDo;
	struct outCheckHelioDist outCh;

	COLLIDER *c1 = &inDo.Collider1,*c2 = &inDo.Collider2,*c;
	double sec;
	unsigned int nCol=0,nMis=0,nMrg=0,nBnc=0,nFrg=0;
		 
        inDo.bPeriodic = msr->param.bPeriodic;

	/* we first check heliocentric distance */
	pstCheckHelioDist(msr->pst,NULL,0,&outCh,NULL);
	msr->dEcoll += outCh.dT;
	msr->dSunMass += outCh.dSM;

	do{	
	    pstNextCollision(msr->pst,NULL,0,&next,NULL);
	    /*printf("%i,%i\n",next.iOrder1,next.iOrder2);*/
	    
	    /* process the collision */
	    if (COLLISION(next.dt)) {
		
		assert(next.iOrder1 >= 0);
		assert(next.iOrder2 >= 0);
		
		inDo.dt = next.dt;
		inGet.iOrder = next.iOrder1;
		pstGetColliderInfo(msr->pst,&inGet,sizeof(inGet),&outGet,NULL);
		*c1 = outGet.Collider; /* struct copy */
		
		/*printf("%i,%i\n",c1->id.iOrder,inGet.iOrder);*/
		
		assert(c1->id.iOrder == inGet.iOrder);
		
		inGet.iOrder = next.iOrder2;
		pstGetColliderInfo(msr->pst,&inGet,sizeof(inGet),&outGet,NULL);
		*c2 = outGet.Collider;
		/*printf("%i,%i\n",c2->id.iOrder,inGet.iOrder);*/
		assert(c2->id.iOrder == inGet.iOrder);
		inDo.CP = msr->param.CP; /* copy collisional parmas */
		
		pstDoCollision(msr->pst,&inDo,sizeof(inDo),&outDo,NULL);
		
		msr->dEcoll += outDo.dT; /* account for kinetic energy loss + (potential)*/
		
		++nCol;
		switch (outDo.iOutcome) {
		    case MISS:
			++nMis;
			--nCol;
			break;
		    case MERGE:
			++nMrg;
			break;
		    case BOUNCE:
			++nBnc;
			break;
		    case FRAG:
			++nFrg;
			break;
		    default:
			assert(0); /* unknown outcome */
		}	
		
#ifdef IGNORE_FOR_NOW/*DEBUG*/
		if (outDo.iOutcome & FRAG) {
		    /* see Zoe's version */
		}
#endif
		switch (msr->param.iCollLogOption) { /* log collision if requested */
		    case COLL_LOG_NONE:
			break;
		    case COLL_LOG_VERBOSE:
		    {
			FILE *fp;
			int i;
			
			fp = fopen(msr->param.achCollLog,"a");
			assert(fp != NULL);
			
			/* for (i=0;i<3;i++) {
			   c1->r[i] += c1->v[i]*next.dt;
			   c2->r[i] += c2->v[i]*next.dt;
			   } */
			
			fprintf(fp,"%s-%s COLLISION:T=%e\n",
				_msrParticleLabel(msr,c1->iColor),
				_msrParticleLabel(msr,c2->iColor),dTime);
			
			fprintf(fp,"***1:p=%i,o=%i,i=%i,oi=%i,M=%e,R=%e,dt=%e,rung=%i,"
				"r=(%e,%e,%e),v=(%e,%e,%e),w=(%e,%e,%e)\n",
				c1->id.iPid,c1->id.iOrder,c1->id.iIndex,c1->id.iOrgIdx,
				c1->fMass,c1->fRadius,c1->dt,c1->iRung,
				c1->r[0],c1->r[1],c1->r[2],
				c1->v[0],c1->v[1],c1->v[2],
				c1->w[0],c1->w[1],c1->w[2]);
			
			fprintf(fp,"***2:p=%i,o=%i,i=%i,oi=%i,M=%e,R=%e,dt=%e,rung=%i,"
				"r=(%e,%e,%e),v=(%e,%e,%e),w=(%e,%e,%e)\n",
				c2->id.iPid,c2->id.iOrder,c2->id.iIndex,c2->id.iOrgIdx,
				c2->fMass,c2->fRadius,c2->dt,c2->iRung,
				c2->r[0],c2->r[1],c2->r[2],
				c2->v[0],c2->v[1],c2->v[2],
				c2->w[0],c2->w[1],c2->w[2]);
			fprintf(fp,"***OUTCOME=%s dT=%e\n",
				outDo.iOutcome == MISS ? "MISS" :
				outDo.iOutcome == MERGE ? "MERGE" :
				outDo.iOutcome == BOUNCE ? "BOUNCE" :
				outDo.iOutcome == FRAG ? "FRAG" : "UNKNOWN",outDo.dT);
			for (i=0;i<(outDo.nOut < MAX_NUM_FRAG ? outDo.nOut : MAX_NUM_FRAG);i++) {
			    c = &outDo.Out[i];
			    
			    fprintf(fp,"***out%i:p=%i,o=%i,i=%i,oi=%i,M=%e,R=%e,rung=%i,"
				    "r=(%e,%e,%e),v=(%e,%e,%e),w=(%e,%e,%e)\n",i,
				    c->id.iPid,c->id.iOrder,c->id.iIndex,c->id.iOrgIdx,
				    c->fMass,c->fRadius,c->iRung,
				    c->r[0],c->r[1],c->r[2],
				    c->v[0],c->v[1],c->v[2],
				    c->w[0],c->w[1],c->w[2]);
			}
			fclose(fp);
			break;
		    }
		    case COLL_LOG_TERSE:
		    {
			/*
			** FORMAT: For each event, time (double), collider 1 iOrgIdx
			** (int), collider 2 iOrgIdx (int), number of post-collision
			** particles (int), iOrgIdx for each of these (n * int).
			*/
			
			FILE *fp;
			XDR xdrs;		       
			int i;
			
			if (outDo.iOutcome != MERGE && outDo.iOutcome != FRAG)
			    break; /* only care when particle indices change */
			fp = fopen(msr->param.achCollLog,"a");
			assert(fp != NULL);
			xdrstdio_create(&xdrs,fp,XDR_ENCODE);
			
			(void) xdr_double(&xdrs,&dTime);
			/* MERGE =1, BOUNCE =2*/
			(void) xdr_int(&xdrs,&outDo.iOutcome); 
			
			/* info for c1*/
			(void) xdr_int(&xdrs,&c1->iColor);
			(void) xdr_int(&xdrs,&c1->id.iOrgIdx);			
			(void) xdr_double(&xdrs,&c1->fMass);
			(void) xdr_double(&xdrs,&c1->fRadius);
			for (i=0;i<N_DIM;i++)
			    (void)xdr_double(&xdrs,&c1->r[i]);
			for (i=0;i<N_DIM;i++)
			    (void)xdr_double(&xdrs,&c1->v[i]);
			for (i=0;i<N_DIM;i++)
			    (void)xdr_double(&xdrs,&c1->w[i]);
			
			/* info for c2*/
			(void) xdr_int(&xdrs,&c2->iColor);
			(void) xdr_int(&xdrs,&c2->id.iOrgIdx);			
			(void) xdr_double(&xdrs,&c2->fMass);
			(void) xdr_double(&xdrs,&c2->fRadius);
			for (i=0;i<N_DIM;i++)
			    (void)xdr_double(&xdrs,&c2->r[i]);
			for (i=0;i<N_DIM;i++)
			    (void)xdr_double(&xdrs,&c2->v[i]);
			for (i=0;i<N_DIM;i++)
			    (void)xdr_double(&xdrs,&c2->w[i]);
			
			xdr_destroy(&xdrs);
			(void) fclose(fp);
			break;
		    }
		    default:
			assert(0); /* invalid collision log option */
		} /* logging */
	    } /* if collision */
	} while (COLLISION(next.dt));	

	msrAddDelParticles(msr); /* clean up any deletions */
	
	if (msr->param.bVStep) {
	    double dsec = msrTime() - sec;
	    printf("%i collision%s: %i miss%s, %i merger%s, %i bounce%s, %i frag%s\n",
		   nCol,nCol==1?"":"s",nMis,nMis==1?"":"es",nMrg,nMrg==1?"":"s",
		   nBnc,nBnc==1?"":"s",nFrg,nFrg==1?"":"s");
	    printf("Collision search completed, time = %g sec\n",dsec);
	}
}

#ifdef SYMBA
void msrTopStepSymba(MSR msr,
		     double dStep,	/* Current step */
		     double dTime,	/* Current time */
		     double dDelta,	/* Time step */
		     int iRung,		/* Rung level */
		     int iKickRung,	/* Gravity on all rungs from iRung
					   to iKickRung */
		     int iRungVeryActive,  /* current setting for iRungVeryActive */
		     int iAdjust,		/* Do an adjust? */
		     double *pdActiveSum,
		     int *piSec)
{
    double dMass = -1.0;
    uint64_t nActive;
    int bSplitVA;
    double dDeltaTmp;
    int i;
    const int bNeedEwald = 0;
  
    msrActiveRung(msr,0,1); /* activate all particles */ 		 
    /* F_0 for particles at iRung = 0 and 1 */		
    msrKickKDKOpen(msr,dTime,0.5*dDelta);   
    /*	
    ** drift all particles to the end of time step
    */
    msrKeplerDrift(msr, dDelta);
    /*
     * check min.distance (drmin2) during drift 
     */ 	
    msrSmooth(msr,dTime,SMX_SYMBA,0,0,TYPE_ALL);
    /*
     * Determine p->iRung from p->drmin 
     ** If drmin2 < 3Hill, (but drmin > 3Hill), this interacting pair 
     ** is sent to iRung = 1. Positions and velocities before the drift 
     ** are retrived (x = xb, v = vb) for particles with p_iRung >=1
     */ 
    msrDrminToRung(msr,iRung);
    /*
     * Activate VeryActives
     */
    msrActiveRung(msr,1,1); /* msr->iRungVeryActive+1 = 1 */
  
    /*assert(msr->nActive == 0); temporaryly */
    if(msr->nActive){ /* if any particles in close encounters */
    if(msrDoGravity(msr)) {	   
      /*
      ** Domain decomposition for parallel exclude very active is going to be 
      ** placed here shortly.
      */
      bSplitVA = 1;
      msrDomainDecomp(msr,iRung,1,bSplitVA);	    
    }
    /*
     * Perform timestepping of particles in close encounters on 
     * individual processors.
     */
    msrStepVeryActiveSymba(msr,dStep,dTime, dDelta, iRung);    
    }

    dTime += dDelta;
    dStep += 1.0;
    
    /*
     * Regular Tree gravity
     */
    msrActiveRung(msr,iKickRung,1); /* iKickRung = 0 */ 
    bSplitVA = 0;
    msrDomainDecomp(msr,iKickRung,1,bSplitVA);
    
    if(msrDoGravity(msr)) {
      msrActiveRung(msr,iKickRung,1);	   
      if (msr->param.bVDetails) {
	printf("%*cGravity, iRung: %d to %d\n",
	       2*iRung+2,' ',iKickRung,msrCurrMaxRung(msr));
      }
      msrBuildTree(msr,dMass,dTime,bNeedEwald);    
      msrGravity(msr,dTime,dStep,0,0,piSec,&nActive);
      *pdActiveSum += (double)nActive/msr->N;
    }
        
    msrActiveRung(msr,iRung,1);			
    msrKickKDKClose(msr,dTime-0.5*dDelta,0.5*dDelta);
    
    /* linear shifts due to Sun's velocity. 
       The next half-step shifts are included */
    msrDriftSun(msr,dTime-0.5*dDelta,dDelta); 
    
    if(msr->param.bCollision){   
      msrDoCollision(msr,dTime,dDelta);     
    }
    
}

void
msrStepVeryActiveSymba(MSR msr, double dStep, double dTime, double dDelta,
		     int iRung){
    struct inStepVeryActiveS in;
    struct outStepVeryActiveS out;    
    
    in.dStep = dStep;
    in.dTime = dTime;
    in.dDelta = dDelta;
    in.iRung = iRung;
    /* could set a stricter opening criterion here */
    in.diCrit2 = 1/(msr->dCrit*msr->dCrit);  
    in.nMaxRung = msrCurrMaxRung(msr);
    in.dSunMass = msr->dSunMass;
       
    pstStepVeryActiveSymba(msr->pst, &in, sizeof(in), &out, NULL);

    if(msr->param.bCollision){
    struct outGetVariableVeryActive outGet;
    pstGetVariableVeryActive(msr->pst, NULL, 0, &outGet, NULL);
    msr->dEcoll += outGet.dDeltaEcoll;
    outGet.dDeltaEcoll = 0.0;/*just in case */
    }
    msr->iCurrMaxRung = out.nMaxRung;
}

void msrDrminToRung(MSR msr,int iRung){
  struct inDrminToRung in;
  struct outDrminToRung out;
  int iTempRung,iOutMaxRung;
  char c;
 
  in.iRung = iRung;
  in.iMaxRung = msrMaxRung(msr);

  pstDrminToRung(msr->pst, &in, sizeof(in), &out, NULL);
  iTempRung =msrMaxRung(msr)-1;
  while (out.nRungCount[iTempRung] == 0 && iTempRung > 0) --iTempRung;
  iOutMaxRung = iTempRung;

  /*
    ** Now copy the rung distribution to the msr structure!
    */
  for (iTempRung=0;iTempRung < msrMaxRung(msr);++iTempRung){
    msr->nRung[iTempRung] = out.nRungCount[iTempRung];
  }
  msr->nRung[msrMaxRung(msr)] = 0; /* just for sure */
  msr->iCurrMaxRung = iOutMaxRung;
  
  if (msr->param.bVRungStat) {
    	printf("Rung distribution:\n");
	for (iTempRung=0;iTempRung <= msr->iCurrMaxRung;++iTempRung) {
	    if (out.nRungCount[iTempRung] == 0) continue;
	    if (iTempRung > 0 ) c = 'v'; /* iRungVeryActive = 0*/
	    else c = ' ';
	    printf(" %c rung:%d %d\n",c,iTempRung,out.nRungCount[iTempRung]);
	    }
	printf("\n");
	}
  
  /*
   * Set VeryActive particles
   * Remember, the first very active particle is at iRungVeryActive + 1 
   */
  msrSetRungVeryActive(msr, 0); /* iRungVeryActive = 0*/
}

void msrDriftSun(MSR msr,double dTime,double dDelta){
	struct inDriftSun in;
	struct outMomSun outm;

	int j;
	/* Calculate Sun's momentum */
	pstMomSun(msr->pst,NULL,0,&outm,NULL); 
       
	for (j=0;j<3;++j){ 	    
                in.vSun[j] = outm.momSun[j]/msr->dSunMass;	
	}
	in.dDelta = dDelta;
	pstDriftSun(msr->pst,&in,sizeof(in),NULL,NULL);

} 

void msrKeplerDrift(MSR msr,double dDelta){
    struct inKeplerDrift in;
    
    in.dDelta = dDelta;
    in.dSunMass = msr->dSunMass;
    pstKeplerDrift(msr->pst,&in,sizeof(in),NULL,NULL);
	
} 


#endif /* SYMBA */
#endif /* PLANETS*/

void msrWrite(MSR msr,char *pszFileName,double dTime,int bCheckpoint)
{
#ifdef PLANETS 
    msrWriteSS(msr,pszFileName,dTime);
#else
#ifdef USE_MDL_IO
    /* If we are using I/O processors, then we do it totally differently */
    if ( mdlIO(msr->mdl) ) {
	msrIOWrite(msr,pszFileName,dTime,bCheckpoint);
    }
    else
#endif
    {
	_msrWriteTipsy(msr,pszFileName,dTime,bCheckpoint);
    }
#endif
}


double msrRead(MSR msr, int iStep)
{
    double dTime;
#ifdef PLANETS    
    dTime = msrReadSS(msr); /* must use "Solar System" (SS) I/O format... */
#else
    char achFilename[PST_FILENAME_SIZE];

    if ( ! msr->param.achInFile[0] ) {
	printf("No input file specified\n");
	_msrExit(msr,1);
	return -1.0;
    }

    /* Add Data Subpath for local and non-local names. */
    _msrMakePath(msr->param.achDataSubPath,msr->param.achInFile,achFilename);

#ifdef USE_MDL_IO
    /* If we have MDL I/O, and the MDL pattern is part of the file */
    if ( mdlIO(msr->mdl) && strstr(achFilename,"&I")!=NULL ) {
	dTime = _msrIORead(msr,achFilename, iStep);
    }
    else {
#endif


#ifdef USE_HDF5
    /* We can automatically detect if a given file is in HDF5 format */
    if ( H5Fis_hdf5(achFilename) ) {
	dTime = _msrReadHDF5(msr,achFilename);
    }
    else
#endif
    {
	dTime = _msrReadTipsy(msr,achFilename);
    }
#endif

#ifdef USE_MDL_IO
    /* If we are using I/O processors, then preallocate space to save */
    if ( mdlIO(msr->mdl) ) {
	struct inIOSetup setup;
	setup.N = msr->N;
	mdlSetComm(msr->mdl,1);
	mdlReqService(msr->mdl,0,IO_SETUP,&setup,sizeof(setup));
	mdlGetReply(msr->mdl,0,NULL,NULL);
	mdlSetComm(msr->mdl,0);
    }
#endif

#ifdef USE_MDL_IO
    }
#endif

    return dTime;
}
