# Download an example FITS file, create a 2D cutout, and save it to a
# new FITS file, including the updated cutout WCS.
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.utils.data import download_file
from astropy.wcs import WCS

# Download the image
filename = 'exampleFile.fits'
position = (477,632)
size = 250

# Load the image and the WCS
hdu = fits.open(filename)[0]
wcs = WCS(hdu.header)

# Make the cutout, including the WCS
cutout = Cutout2D(hdu.data, position=position, size=size, wcs=wcs)

 # Put the cutout image in the FITS HDU
hdu.data = cutout.data

# Update the FITS header with the cutout WCS
hdu.header.update(cutout.wcs.to_header())

# Write the cutout to a new FITS file
cutout_filename = 'example_cutout.fits'
hdu.writeto(cutout_filename, overwrite=True)
