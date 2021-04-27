#ifdef STELLAR_EVOLUTION

#include "master.h"


/* NOTAS:
   - El segundo argumento de la funcion prmAddParam define el nombre con el cual se busca un
   parametro en el archivo de parametros. El sexto argumento define como se puede especificar
   el parametro desde la linea de comandos.
   - Estoy convirtiendo los parametros a unidades internas en msrInitialize.
   - Los parametros son accesibles a traves de pkd->param
   - Include de stellarevolution.h en archivos .c en pst y pkd
 */


/* TODOs y DUDAS:
   - Armar archivos hdf5. Excluyo S y Ca?
   - Entender como se guardan y leen archivos de restart.
   - Imponer conservacion de energia cinetica y momento lineal al distribuir la masa?
   - Preguntar a Isaac sobre si la funcion stevFreeTable esta liberando la memoria
   correctamente
   - Es correcto normalizar la IMF entre 0.1 y 100.0 Msol y despues utilizar tablas desde
   120 Msun?
   - Volver a verificar que las tablas se estan leyendo correctamente
   - Transponer buffers?
 */


#define STEV_CCSN_N_METALLICITY               5
#define STEV_CCSN_N_MASS                     11
#define STEV_AGB_N_METALLICITY                3
#define STEV_AGB_N_MASS                      23
#define STEV_LIFETIMES_N_METALLICITY          6
#define STEV_LIFETIMES_N_MASS                30
#define STEV_N_ELEM     chemistry_element_count   /* Defined in pkd.h */

#define STEV_INTERP_N_MASS                  200


/*
 * ---------------------
 * STRUCTURE DEFINITIONS
 * --------------------- 
 */

typedef struct inStellarEvolution {
   /* Pointer to the function that gives the number of SNIa in [dTime, dTime + dDelta] */
   float (*fcnNumSNIa)(double dTime, double dDelta, double norm, double scale);

   /* Initial mass array for CCSN and AGB tables */
   float afMasses[STEV_INTERP_N_MASS];

   /* Initial Mass Function values array */
   float afIMF[STEV_INTERP_N_MASS];

   /* Core Collapse SNe arrays */
   float afCCSN_Zs[STEV_CCSN_N_METALLICITY];
   float afCCSN_Yields[STEV_CCSN_N_METALLICITY * STEV_N_ELEM * STEV_INTERP_N_MASS];
   float afCCSN_MetalYield[STEV_CCSN_N_METALLICITY * STEV_INTERP_N_MASS];
   float afCCSN_EjectedMass[STEV_CCSN_N_METALLICITY * STEV_INTERP_N_MASS];

   /* AGB arrays */
   float afAGB_Zs[STEV_AGB_N_METALLICITY];
   float afAGB_Yields[STEV_AGB_N_METALLICITY * STEV_N_ELEM * STEV_INTERP_N_MASS];
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
int pkdStellarEvolutionInit(PKD, /* ADD NECESSARY ARGUMENTS */);

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

void stevChabrierIMF(float *, float *, int, double, double);
float stevExponentialNumSNIa(double dTime, double dDelta, double norm, double scale);
float stevPowerlawNumSNIa(double dTime, double dDelta, double norm, double scale);


static inline stevGetIndex1D(const float *restrict pfTable, const int nSize, const float fVal,
			     int *piIdx, float *restrict pfDelta) {
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


#endif	/* STELLAR_EVOLUTION */
