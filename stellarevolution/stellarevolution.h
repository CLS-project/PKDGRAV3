#ifdef STELLAR_EVOLUTION

#include "master.h"


/* NOTAS:
   - El segundo argumento de la funcion prmAddParam define el nombre con el cual se busca un
   parametro en el archivo de parametros. El sexto argumento define como se puede especificar
   el parametro desde la linea de comandos.
   - Estoy convirtiendo los parametros a unidades internas en msrInitialize.
   - Los parametros son accesibles a traves de pkd->param
 */


/* TODOs:
   - Armar archivos hdf5. Excluyo S y Ca?
   - Entender como se guardan y leen archivos de restart.
   - Imponer conservacion de energia cinetica y momento lineal al distribuir la masa?
   - Preguntar a Isaac sobre si la funcion stevFreeTable esta liberando la memoria
   correctamente
 */


#define STEV_CCSN_N_METALLICITY               5
#define STEV_CCSN_N_MASS                     11
#define STEV_AGB_N_METALLICITY                3
#define STEV_AGB_N_MASS                      23
#define STEV_LIFETIMES_N_METALLICITY          6
#define STEV_LIFETIMES_N_MASS                30
#define STEV_N_ELEM     chemistry_element_count   /* AM: Defined in pkd.h */

#define STEV_INTERP_N_MASS                  200


/*
 * ---------------------
 * STRUCTURE DEFINITIONS
 * --------------------- 
 */

typedef struct StellarEvolutionBuffer {
   // Initial mass array for CCSN and AGB tables
   float Masses[STEV_INTERP_N_MASS];

   // Initial Mass Function values array
   float IMF[STEV_INTERP_N_MASS];

   // Core Collapse SNe arrays
   float CCSN_Zs[STEV_CCSN_N_METALLICITY];
   float CCSN_Yields[STEV_N_ELEM * STEV_INTERP_N_MASS * STEV_CCSN_N_METALLICITY];
   float CCSN_MetalYield[STEV_INTERP_N_MASS * STEV_CCSN_N_METALLICITY];
   float CCSN_EjectedMass[STEV_INTERP_N_MASS * STEV_CCSN_N_METALLICITY];

   // AGB arrays
   float AGB_Zs[STEV_AGB_N_METALLICITY];
   float AGB_Yields[STEV_N_ELEM * STEV_INTERP_N_MASS * STEV_AGB_N_METALLICITY];
   float AGB_MetalYield[STEV_INTERP_N_MASS * STEV_AGB_N_METALLICITY];
   float AGB_EjectedMass[STEV_INTERP_N_MASS * STEV_AGB_N_METALLICITY];

   // Type Ia SNe
   float SNIa_EjectedMass[STEV_N_ELEM];
   float SNIa_EjectedMetalMass;

   // Stellar lifetimes arrays
   float Lifetimes_Masses[STEV_LIFETIMES_N_MASS];
   float Lifetimes_Zs[STEV_LIFETIMES_N_METALLICITY];
   float Lifetimes[STEV_LIFETIMES_N_MASS * STEV_LIFETIMES_N_METALLICITY];
} STEV_BUFFER;


typedef struct StellarEvolutionData {
   float CCSN_MinMass;
   float SNIa_MaxMass;

   float SNIa_Norm;
   float SNIa_Scale;
   float (*N_SNIa)(double dTime, double dDelta, float norm, float scale);

   float *Masses;
   float *IMF;

   float *CCSN_Zs;
   float *CCSN_Yields;
   float *CCSN_MetalYield;
   float *CCSN_EjectedMass;

   float *AGB_Zs;
   float *AGB_Yields;
   float *AGB_MetalYield;
   float *AGB_EjectedMass;

   float *SNIa_EjectedMass;
   float *SNIa_EjectedMetalMass;

   float *Lifetimes_Masses;
   float *Lifetimes_Zs;
   float *Lifetimes;
} STEV_DATA;

typedef struct StellarEvolutionRawData {
   int nZs;
   int nSpecs;
   int nMasses;
   float *Zs;
   float *Masses;
   float *Yields;
   float *MetalYield;
   float *EjectedMass;
   float *Lifetimes;
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

void stevReadTable(STEV_RAWDATA *, char *, char **, int);
void stevFreeTable(STEV_RAWDATA *);
void stevReadSNIaTable(STEV_RAWDATA *, char *);
void stevFreeSNIaTable(STEV_RAWDATA *);
void stevReadLifetimesTable(STEV_RAWDATA *, char *);
void stevFreeLifetimesTable(STEV_RAWDATA *);


#endif	/* STELLAR_EVOLUTION */
