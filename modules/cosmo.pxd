import cython

cdef extern from "cosmo.h":

    ctypedef struct classDataStruct:
        int bClass

    ctypedef struct csmVariables:
        int bComove
        double dHubble0
        double dOmega0
        double dLambda
        double dOmegaRad
        double dOmegab
        double dOmegaDE
        double w0
        double wa
        double dSigma8
        double dNormalization
        double dSpectral
        double dRunning
        double dPivot
        double h
        classDataStruct classData

    ctypedef struct csmContext:
        csmVariables val

    void csmInitialize(csmContext **pcsm)
    void csmFinish(csmContext *csm)
    void csmClassRead(csmContext *csm, const char *achFilename, double dBoxSize, double h,
                  int nLinear, const char **aLinear, int nPower, const char **aPower)
    void csmClassGslInitialize(csmContext *csm)

    double csmRadMatEquivalence(csmContext *csm)
    double csmExp2Hub(csmContext *csm, double dExp)
    double csmTime2Hub(csmContext *csm, double dTime)
    double csmExp2Time(csmContext *csm, double dExp)
    double csmTime2Exp(csmContext *csm, double dTime)

    double csmComoveDriftInt(csmContext *csm, double dIExp);
    double csmComoveKickInt(csmContext *csm, double dIExp);
    double csmComoveDriftFac(csmContext *csm, double dTime, double dDelta);
    double csmComoveKickFac(csmContext *csm, double dTime, double dDelta);
    double csmComoveLookbackTime2Exp(csmContext *csm, double dComoveTime);
    void csmComoveGrowth(csmContext *csm, double a, double *D1LPT, double *D2LPT, double *f1LPT, double *f2LPT);

    double csmRhoBar_m    (csmContext *csm, double a)
    double csmRhoBar_lin  (csmContext *csm, double a)
    double csmRhoBar_pk   (csmContext *csm, double a)
    double csmDelta_m     (csmContext *csm, double a,                double k)
    double csmTheta_m     (csmContext *csm, double a,                double k)
    double csmDelta_lin   (csmContext *csm, double a,                double k)
    double csmDelta_pk    (csmContext *csm, double a,                double k)
    double csmDeltaRho_lin(csmContext *csm, double a, double a_next, double k)
    double csmDeltaRho_pk (csmContext *csm, double a,                double k)
    double csmZeta        (csmContext *csm,                          double k)
