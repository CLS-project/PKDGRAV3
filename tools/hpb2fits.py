import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
from sys import argv,exit

hpb_type = np.dtype([('grouped', '=i4'),('ungrouped', '=i4'), ('potential', '=f4')])

bPotential = True

NSIDE = 64
N = 12 * NSIDE**2

hdr = fits.Header()
hdr['EXTEND'] = 'T'
primary_hdu = fits.PrimaryHDU(header=hdr)
primary_hdu.writeto('test.fits')


# Use BinTableHDU as a template
hdr = fits.BinTableHDU(Table(names=['SIGNAL'],dtype=['=f4' if bPotential else '=i4']),name='BINTABLE').header
hdr['ORDERING'] = ("RING","Pixel ordering scheme, either RING or NESTED")
hdr['INDXSCHM'] = ("IMPLICIT","Pixel indexing scheme (IMPLICIT or EXPLICIT)")
hdr['NSIDE']    = (NSIDE,"Resolution parameter for HEALPIX")
hdr['COORDSYS'] = ("C","Pixelisation coordinate system")
hdr['PIXTYPE']  = ("HEALPIX", "HEALPIX Pixelisation")
hdr['NAXIS']    = 2
hdr['NAXIS2']   = N
hdr['NAXIS1']   = 1
hdr['BITPIX']   = -32 if bPotential else 32

hdu = fits.StreamingHDU('test.fits',hdr)
for file in argv[1:]:
    with open(file,'rb') as hpb:
        data = np.fromfile(hpb,dtype=hpb_type)
        data = pd.DataFrame(data,columns=data.dtype.names)
        if bPotential:
            hdu.write(data['potential'].values)
        else:
            data = (data.grouped+data.ungrouped).to_frame()
            data.columns=['signal']
            hdu.write(data['signal'].values)
hdu.close()
