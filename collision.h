#ifdef PLANETS
#ifndef COLLISION_HINCLUDED
#define COLLISION_HINCLUDED

#include "pkd.h" /*for PARTICLE struct */
#include "ssio.h"  /*for SSDATA struct */
#include "parameters.h" /* for COLLISION_PARAMS struct */

#define NORMAL 0	/* drift types */
#define KEPLER 1

int REJECT(PARTICLE *p);
#define REJECT(p) ((p)->dtCol < 0)

#define DBL_MAX 1.7976931348623157E+308
int COLLISION(double t);
#define COLLISION(t) ((t) < DBL_MAX)

unsigned int BIT(unsigned int n);
#define BIT(n) (1 << (n))

#define COLL_LOG_NONE		0
#define COLL_LOG_VERBOSE	1
#define COLL_LOG_TERSE		2

#define MISS	0
#define MERGE	BIT(0)
#define BOUNCE	BIT(1)
#define FRAG	BIT(2)

#define MAX_NUM_FRAG 4

enum {ConstEps,Frosty200,Frosty120,Compacted,Glancing}; /* bounce options */

enum {EscVel,MaxTrv}; /* slide options */

typedef struct {
    int iPid;
    int iOrder;
    int iIndex;
    int iOrgIdx;
    } PARTICLE_ID;

typedef struct {
    PARTICLE_ID id;
    FLOAT fMass;
    FLOAT fRadius;
    FLOAT r[3];
    FLOAT v[3];
    FLOAT w[3];
#ifdef HERMITE
    FLOAT a[3];
    FLOAT ad[3];
#endif
#ifdef SYMBA
    FLOAT drmin;
    int   iOrder_VA[5]; /* iOrder's of particles within 3 hill radius*/
    int   i_VA[5];    /* pointers of particles */
    int   n_VA;       /* number of particles */
    double  hill_VA[5]; /* mutual hill radius calculated in grav.c */
#endif
    int iColor;
    FLOAT dt;
    int iRung;
    } COLLIDER;

void PutColliderInfo(const COLLIDER *c,int iOrder2,PARTICLE *p,double dt);
void pkdNextCollision(PKD pkd, double *dtCol, int *iOrder1, int *iOrder2);
void pkdGetColliderInfo(PKD pkd, int iOrder, COLLIDER *c);
void pkdDoCollision(PKD pkd, double dt, const COLLIDER *c1, const COLLIDER *c2,
		    int bPeriodic, const COLLISION_PARAMS *CP, int *piOutcome,
		    double *dT,	COLLIDER *cOut, int *pnOut);
void pkdDoCollisionVeryActive(PKD pkd,double dTime);
static char * _pkdParticleLabel(PKD pkd,int iColor);
void pkdGetVariableVeryActive(PKD pkd, double *dDeltaEcoll);
void pkdCheckHelioDist(PKD pkd,double *dT,double *dSM);

#endif

#endif /* PLANETS */

