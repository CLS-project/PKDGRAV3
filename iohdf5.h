#ifndef IOHDF5_H
#define IOHDF5_H
#include <hdf5.h>

#ifdef __cplusplus
extern "C" {
#endif

#define H5assert(rc) assert( rc >= 0 )

    typedef struct {
	FLOAT v[3];
	} ioV3;

    typedef uint64_t PINDEX;

    typedef struct {
	PINDEX iStart;            /* Start iOrder number */
	PINDEX iNext;             /* Next iOrder number (starts at iStart) */
	PINDEX *iOrder;           /* Buffered iOrder numbers (or NULL) */
	hid_t  setOrder_id;
	} IOORDER;

    typedef struct {
	uint8_t iClass;
	PINDEX iOrderStart;       /* Start index of this class */
	FLOAT  fMass;             /* Particle mass */
	FLOAT  fSoft;             /* Softening */
	PINDEX nCount;            /* Number of particles with this class */
	} classEntry;

    typedef struct {
	uint_fast32_t nClasses;   /* Number of different classes */
	classEntry Class[256];
	hid_t   setClass_id;
	uint8_t *piClass;         /* Class index (or NULL) */
	FLOAT   *fMass;
	FLOAT   *fSoft;
	} IOCLASS;

    typedef struct {
	char    szGroupName[32];

	PINDEX  nTotal;           /* Total number of particles in the file */
	PINDEX  iOffset;          /* Particle offset into the file */
	uint_fast32_t iIndex;     /* Index of next particle in memory */
	uint_fast32_t nBuffered;  /* Number of buffered particles */
	hid_t   group_id;         /* Group /dark, /gas, /star, etc. */
	hid_t   setR_id;          /* set of Positions */
	hid_t   setV_id;          /* set of Velocities */
	IOORDER Order;
	IOCLASS Class;

	ioV3    *R;               /* Positions (always present) */
	ioV3    *V;               /* Velocities (always present) */

	} IOBASE;


    struct ioHDF5;
    typedef struct ioHDF5v {
	struct ioHDF5v *next;
	struct ioHDF5 *io;
	hid_t   diskFloat;        /* FLOAT on disk: float or double */
	hid_t   set_id;           /* vector */

	PINDEX  iOffset;          /* Particle offset into the file */
	PINDEX  nTotal;
	uint_fast32_t iIndex;     /* Next particle in memory to read */
	uint_fast32_t nBuffered;  /* Number of buffered particles */

	char    name[32];

	float   *s;
	double  *d;

	} *IOHDF5V;

    typedef struct ioHDF5 {
	hid_t   fileID;           /* HDF5 file handle */
	uint_fast32_t iChunkSize;
	hid_t   parametersID;

	hid_t   memFloat;         /* FLOAT in memory: float or double */
	hid_t   diskFloat_R;      /* FLOAT on disk: float or double */
	hid_t   diskFloat_V;      /* FLOAT on disk: float or double */

	IOBASE  darkBase;
	IOBASE  gasBase;
	IOBASE  starBase;

	unsigned char bRead;
	unsigned char bWrite;

	IOHDF5V vectorList;
	} *IOHDF5;

#define IOHDF5_SINGLE 0
#define IOHDF5_DOUBLE 1

    IOHDF5 ioHDF5Initialize( hid_t fileID, hid_t iChunkSize, int bDouble );

    void ioHDF5Finish( IOHDF5 io );

    PINDEX ioHDF5DarkCount( IOHDF5 io );
    PINDEX ioHDF5GasCount( IOHDF5 io );
    PINDEX ioHDF5StarCount( IOHDF5 io );

    void ioHDF5AddDark( IOHDF5 io, PINDEX iOrder,
			const FLOAT *r, const FLOAT *v,
			FLOAT fMass, FLOAT fSoft, float fPot );

    void ioHDF5AddGas(  IOHDF5 io, PINDEX iOrder,
			const FLOAT *r, const FLOAT *v,
			FLOAT fMass, FLOAT fSoft, float fPot,
			FLOAT fTemp, FLOAT fMetals);

    void ioHDF5AddStar( IOHDF5 io, PINDEX iOrder,
			const FLOAT *r, const FLOAT *v,
			FLOAT fMass, FLOAT fSoft, float fPot,
			FLOAT fMetals, FLOAT fTForm);

    IOHDF5V ioHDFF5OpenVector( IOHDF5 io, const char *name, int bDouble );
    IOHDF5V ioHDFF5NewVector( IOHDF5 io, const char *name, int bDouble );
    void ioHDF5AddVector( IOHDF5V iov, PINDEX iOrder, FLOAT v );
    FLOAT ioHDF5GetVector( IOHDF5V iov );
    void ioHDF5SeekVector( IOHDF5V iov, PINDEX Offset );


    void ioHDF5WriteAttribute( IOHDF5 io, const char *name,
			       hid_t dataType, void *data );

    int ioHDF5ReadAttribute( IOHDF5 io, const char *name,
			     hid_t dataType, void *data );

    int  ioHDF5GetDark( IOHDF5 io, PINDEX *iOrder,
			FLOAT *r, FLOAT *v,
			FLOAT *fMass, FLOAT *fSoft, float *fPot );

    int  ioHDF5GetGas(  IOHDF5 io, PINDEX *iOrder,
			FLOAT *r, FLOAT *v,
			FLOAT *fMass, FLOAT *fSoft, float *fPot,
			FLOAT *fTemp, FLOAT *fMetals);
    int  ioHDF5GetStar( IOHDF5 io, PINDEX *iOrder,
			FLOAT *r, FLOAT *v,
			FLOAT *fMass, FLOAT *fSoft, float *fPot,
			FLOAT *fMetals, FLOAT *fTForm);

    void ioHDF5SeekDark( IOHDF5 io, PINDEX Offset );
    void ioHDF5SeekGas(  IOHDF5 io, PINDEX Offset );
    void ioHDF5SeekStar( IOHDF5 io, PINDEX Offset );
#ifdef __cplusplus
    }
#endif

#endif
