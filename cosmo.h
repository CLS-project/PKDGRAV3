#ifndef COSMO_HINCLUDED
#define COSMO_HINCLUDED

typedef struct csmContext {
    int bComove;
    double dHubble0;
    double dOmega0;
    double dLambda;
    double dOmegaRad;
    double dOmegab;
    double dOmegaDE;
    double w0;
    double wa;
    double dSigma8;
    double dSpectral;
    } * CSM;

void csmInitialize(CSM *pcsm);
void csmFinish(CSM csm);
double csmExp2Hub(CSM csm, double dExp);
double csmTime2Hub(CSM csm,double dTime);
double csmExp2Time(CSM csm,double dExp);
double csmTime2Exp(CSM csm,double dTime);
double csmComoveDriftInt(CSM csm, double dIExp);
double csmComoveKickInt(CSM csm, double dIExp);
double csmComoveDriftFac(CSM csm,double dTime,double dDelta);
double csmComoveKickFac(CSM csm,double dTime,double dDelta);
double csmComoveLookbackTime2Exp(CSM csm,double dComoveTime);

/*
 ** returns the speed of light in simulation units, given
 ** the simulation length unit in h^-1 Mpc.
 */
static inline double dLightSpeedSim(double dMpcUnit)
{
	/*
	 ** Find the speed of light in simulation units.
	 **
	 ** c[Mpc/Gyr] = c[cm/s] * Julian Year[s] / pc[cm] * 1000 
	 ** c_sim = c[Mpc/Gyr] * (x Gyrs/ 1 sim time) * ( 1 sim length/Boxsize (Mpc))
	 ** x = 1/sqrt(4.498*h*h*2.776e-4)
	 */
	return(8676.85/dMpcUnit);
	}


#endif
