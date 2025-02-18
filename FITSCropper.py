from astropy.io import fits
import numpy as np
from astropy.nddata import Cutout2D
from astropy import units as u
from astropy.wcs import WCS
import os

fits_files = [
    #enter the name of your FITS files as a list, seperated by commas
]
base_path = "/enter/your/base/path"

for i, file in enumerate(fits_files):
    file_path = base_path + file

    with fits.open(file_path) as hdul:
        print(f"Processing {file}...")
        hdu = hdul[0]  # Open the first HDU
        img = hdu.data  # Extract image data
        wcs = WCS(hdu.header)  # Extract WCS
        print(f"{file} wcs: ")
        # Define cutout parameters (ensure position is in pixels)
        position = (#x_center, y_center)  # Pixel position (must be integers)
        size = (#height, width)  # Size in pixels (height, width)

        # Create cutout
        cutout = Cutout2D(img, position=position, size=size, wcs=wcs)

        # Update FITS HDU
        hdu.data = cutout.data
        hdu.header.update(cutout.wcs.to_header())

        # Save cutout to new FITS file
        cutout_filename = base_path + 'your desired file name' + file

        hdu.writeto(cutout_filename, overwrite=True)

    print(f"Saved cutout to {cutout_filename} in {os.path.abspath(cutout_filename)}")
