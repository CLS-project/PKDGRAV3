#ifdef STELLAR_EVOLUTION


#include "master.h"


/* NOTAS:
   - El segundo argumento de la funcion prmAddParam define el nombre con el cual se busca un
   parametro en el archivo de parametros. El sexto argumento define como se puede especificar
   el parametro desde la linea de comandos.
   - Al momento de distribuir las tablas, todas las cantidades y parametros relevantes (tiempos
   y masas) han sido convertidos a unidades internas del codigo.
*/


/* TODOs y DUDAS:
   - Armar archivos hdf5
   - Implementar la inicializacion cuando se empieza de un archivo restart o de una snapshot
   - Considerar cambiar a float parametros que se usan bastante
   - Convertir parametros de masa a log?
   - Por que no estoy utilizando las funciones de interpolacion de interpolate.h?
   - Adicionar free(pkd->StelEvolData) en algun lado
*/


#define STEV_CCSN_N_METALLICITY               5
#define STEV_CCSN_N_MASS                     11
#define STEV_AGB_N_METALLICITY                3
#define STEV_AGB_N_MASS                      23
#define STEV_LIFETIMES_N_METALLICITY          6
#define STEV_LIFETIMES_N_MASS                30
#define STEV_N_ELEM                           9

#define STEV_INTERP_N_MASS                  200


/*
 * ---------------------
 * STRUCTURE DEFINITIONS
 * --------------------- 
 */

typedef struct inStellarEvolution {
   /* Pointer to the function that gives the number of SNIa in [dTime, dTime + dDelta] */
   float (*fcnNumSNIa)(PKD pkd, STARFIELDS *pStar, float fTime);

   /* Initial mass array for CCSN and AGB tables */
   float afMasses[STEV_INTERP_N_MASS];

   /* Initial Mass Function values array */
   float afIMF[STEV_INTERP_N_MASS];

   /* Core Collapse SNe arrays */
   float afCCSN_Zs[STEV_CCSN_N_METALLICITY];
   float afCCSN_Yields[STEV_N_ELEM * STEV_CCSN_N_METALLICITY * STEV_INTERP_N_MASS];
   float afCCSN_MetalYield[STEV_CCSN_N_METALLICITY * STEV_INTERP_N_MASS];
   float afCCSN_EjectedMass[STEV_CCSN_N_METALLICITY * STEV_INTERP_N_MASS];

   /* AGB arrays */
   float afAGB_Zs[STEV_AGB_N_METALLICITY];
   float afAGB_Yields[STEV_N_ELEM * STEV_AGB_N_METALLICITY * STEV_INTERP_N_MASS];
   float afAGB_MetalYield[STEV_AGB_N_METALLICITY * STEV_INTERP_N_MASS];
   float afAGB_EjectedMass[STEV_AGB_N_METALLICITY * STEV_INTERP_N_MASS];

   /* Type Ia SNe */
   float afSNIa_EjectedMass[STEV_N_ELEM];
   float fSNIa_EjectedMetalMass;

   /* Stellar lifetimes arrays */
   float afLifetimes_Zs[STEV_LIFETIMES_N_METALLICITY];
   float afLifetimes_Masses[STEV_LIFETIMES_N_MASS];
   float afLifetimes[STEV_LIFETIMES_N_METALLICITY * STEV_LIFETIMES_N_MASS];
} STEV_DATA;


typedef struct StellarEvolutionRawData {
   int nZs;
   int nSpecs;
   int nMasses;
   float *pfZs;
   float *pfMasses;
   float *pfYields;
   float *pfMetalYield;
   float *pfEjectedMass;
   float *pfLifetimes;
} STEV_RAWDATA;



/*
 * --------------
 * MAIN FUNCTIONS
 * --------------
 */

void msrStellarEvolutionInit(MSR);
int pstStellarEvolutionInit(PST, void *, int, void *, int);
int pkdStellarEvolutionInit(PKD, STEV_DATA *);

void smChemEnrich(PARTICLE *p, float fBall, int nSmooth, NN *nnList, SMF *smf);
void initChemEnrich(void *vpkd, void *vp);
void combChemEnrich(void *vpkd, void *vp1, void *vp2);



/*
 * ----------------
 * HELPER FUNCTIONS
 * ----------------
 */

STEV_RAWDATA *stevReadTable(char *, char **, int);
STEV_RAWDATA *stevReadSNIaTable(char *);
STEV_RAWDATA *stevReadLifetimesTable(char *);
void stevFreeTable(STEV_RAWDATA *);
void stevFreeSNIaTable(STEV_RAWDATA *);
void stevFreeLifetimesTable(STEV_RAWDATA *);

void stevChabrierIMF(float *, int, double, double, float *);
float stevExponentialNumSNIa(PKD, STARFIELDS *, float);
float stevPowerlawNumSNIa(PKD, STARFIELDS *, float);


static inline void stevGetIndex1D(const float *restrict pfTable, const int nSize,
				  const float fVal, int *piIdx, float *restrict pfDelta) {
   const float epsilon = 1e-4f;

   if (fVal < pfTable[0] + epsilon) {
      /* We are below the first element */
      *piIdx = 0;
      *pfDelta = 0.0f;
   }
   else if (fVal > pfTable[nSize - 1] - epsilon) {
      /* We are beyond the last element */
      *piIdx = nSize - 2;
      *pfDelta = 1.0f;
   }
   else {
      /* Normal case */
      int i;
      for (i = 1; (i < nSize - 1) && (fVal > pfTable[i]); i++)
	 ;

      *piIdx = --i;
      *pfDelta = (fVal - pfTable[i]) / (pfTable[i + 1] - pfTable[i]);
   }
}

/* Por ahora voy a asumir que todo esta en log, excepto fMass */
static inline float stevLifetimeFunction(PKD pkd, STARFIELDS *pStar, float fMass) {
   int iIdxMass;
   float fDeltaMass;
   stevGetIndex1D(pkd->StelEvolData->afLifetimes_Masses, STEV_LIFETIMES_N_MASS,
		  log10(fMass), &iIdxMass, &fDeltaMass);

   int oLowerZ = pStar->Lifetimes.iIdxZ * STEV_LIFETIMES_N_MASS;
   float fDeltaZ = pStar->Lifetimes.fDeltaZ;

   float *afTimesLowerZ = pkd->StelEvolData->afLifetimes + oLowerZ;
   float *afTimesUpperZ = afTimesLowerZ + STEV_LIFETIMES_N_MASS;

   float fLogTime = afTimesLowerZ[iIdxMass] * (1.0f - fDeltaZ) * (1.0f - fDeltaMass);
   fLogTime += afTimesLowerZ[iIdxMass + 1] * (1.0f - fDeltaZ) * fDeltaMass;
   fLogTime += afTimesUpperZ[iIdxMass] * fDeltaZ * (1.0f - fDeltaMass);
   fLogTime += afTimesUpperZ[iIdxMass + 1] * fDeltaZ * fDeltaMass;

   return pow(10.0, fLogTime);
}

/* Por ahora voy a asumir que todo esta en log, excepto fTime */
static inline float stevInverseLifetimeFunction(PKD pkd, STARFIELDS *pStar, float fTime) {
   int oLowerZ = pStar->Lifetimes.iIdxZ * STEV_LIFETIMES_N_MASS;
   float fDeltaZ = pStar->Lifetimes.fDeltaZ;

   float *afTimesLowerZ = pkd->StelEvolData->afLifetimes + oLowerZ;
   float *afTimesUpperZ = afTimesLowerZ + STEV_LIFETIMES_N_MASS;

   int iIdxTimeLowerZ;
   float fDeltaTimeLowerZ;
   stevGetIndex1D(afTimesLowerZ, STEV_LIFETIMES_N_MASS, log10(fTime),
		  &iIdxTimeLowerZ, &fDeltaTimeLowerZ);

   int iIdxTimeUpperZ;
   float fDeltaTimeUpperZ;
   stevGetIndex1D(afTimesUpperZ, STEV_LIFETIMES_N_MASS, log10(fTime),
		  &iIdxTimeUpperZ, &fDeltaTimeUpperZ);

   float *afMasses = pkd->StelEvolData->afLifetimes_Masses;

   float fLogMass = afMasses[iIdxTimeLowerZ] * (1.0f - fDeltaZ) * (1.0f - fDeltaTimeLowerZ);
   fLogMass += afMasses[iIdxTimeLowerZ + 1] * (1.0f - fDeltaZ) * fDeltaTimeLowerZ;
   fLogMass += afMasses[iIdxTimeUpperZ] * fDeltaZ * (1.0f - fDeltaTimeUpperZ);
   fLogMass += afMasses[iIdxTimeUpperZ + 1] * fDeltaZ * fDeltaTimeUpperZ;

   return pow(10.0, fLogMass);
}


#endif	/* STELLAR_EVOLUTION */
