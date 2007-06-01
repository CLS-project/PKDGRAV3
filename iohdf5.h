#ifndef IOHDF5_H
#define IOHDF5_H
#include <hdf5.h>

#define H5assert(rc) assert( rc >= 0 )

typedef struct {
    FLOAT v[3];
} ioV3;

typedef uint32_t PINDEX;

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
} classEntry;

typedef struct {
    int nClasses;             /* Number of different classes */
    classEntry class[256];
    hid_t   setClass_id;
    uint8_t *iClass;          /* Class index (or NULL) */
    FLOAT   *fMass;
    FLOAT   *fSoft;
} IOCLASS;

typedef struct {
    char    szGroupName[32];

    PINDEX  iOffset;          /* Particle offset into the file */
    int     iIndex;           /* Index of next particle in memory */
    int     nTotal;           /* Total number */
    int     nBuffered;        /* Number of buffered particles */
    hid_t   group_id;         /* Group /dark, /gas, /star, etc. */
    hid_t   setR_id;          /* set of Positions */
    hid_t   setV_id;          /* set of Velocities */
    IOORDER Order;
    IOCLASS Class;

    ioV3    *R;               /* Positions (always present) */
    ioV3    *V;               /* Velocities (always present) */

} IOBASE;

typedef struct ioHDF5 {
    hid_t   fileID;           /* HDF5 file handle */
    int     iChunkSize;
    hid_t   parametersID;

    hid_t   memFloat;         /* FLOAT in memory: float or double */
    hid_t   diskFloat;        /* FLOAT on disk: float or double */

    IOBASE  darkBase;
    IOBASE  gasBase;
    IOBASE  starBase;

} *IOHDF5;

IOHDF5 ioHDF5Initialize( hid_t fileID, hid_t iChunkSize, int bSingle );

void ioHDF5Finish( IOHDF5 io );

void ioHDF5AddDark( IOHDF5 io, PINDEX iOrder,
		    const FLOAT *r, const FLOAT *v,
		    FLOAT fMass, FLOAT fSoft, FLOAT fPot );

void ioHDF5AddGas(  IOHDF5 io, PINDEX iOrder,
		    const FLOAT *r, const FLOAT *v,
		    FLOAT fMass, FLOAT fSoft, FLOAT fPot,
		    FLOAT fTemp, FLOAT fMetals);

void ioHDF5AddStar( IOHDF5 io, PINDEX iOrder,
		    const FLOAT *r, const FLOAT *v,
		    FLOAT fMass, FLOAT fSoft, FLOAT fPot,
		    FLOAT fMetals, FLOAT fTForm);

void ioHDF5WriteAttribute( IOHDF5 io, const char *name,
			   hid_t dataType, void *data );

int  ioHDF5GetDark( IOHDF5 io, PINDEX *iOrder,
		    FLOAT *r, FLOAT *v,
		    FLOAT *fMass, FLOAT *fSoft, FLOAT *fPot );

int  ioHDF5GetGas(  IOHDF5 io, PINDEX *iOrder,
		    FLOAT *r, FLOAT *v,
		    FLOAT *fMass, FLOAT *fSoft, FLOAT *fPot,
		    FLOAT *fTemp, FLOAT *fMetals);
int  ioHDF5GetStar( IOHDF5 io, PINDEX *iOrder,
		    FLOAT *r, FLOAT *v,
		    FLOAT *fMass, FLOAT *fSoft, FLOAT *fPot,
		    FLOAT *fMetals, FLOAT *fTForm);

void ioHDF5SeekDark( IOHDF5 io, PINDEX Offset );
void ioHDF5SeekGas(  IOHDF5 io, PINDEX Offset );
void ioHDF5SeekStar( IOHDF5 io, PINDEX Offset );

#endif
