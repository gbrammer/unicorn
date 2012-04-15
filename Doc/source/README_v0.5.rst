==========================
3D-HST v0.5 - Apr 12, 2012
==========================

Introduction
============

This release (v0.5) is intended to provide interested users with an idea of the data quality and format of the 3D-HST WFC3/G141 grism products.  Full catalogs of the spectra, photometry, and products derived from these (redshifts, stellar-population properties, etc.) are still under development and are expected to be released by Fall 2012.  This README provides information on the specific products released here.  For more general information on the survey, see the summary paper: 

    `"3D-HST: A Wide-field Grism Spectroscopic Survey with the Hubble Space Telescope"`, Brammer et al., 2012, ApJS-submitted (`arXiv:xxx <http://www.arXiv.org | arXiv>`_)

Please cite the above reference in any publications that make use of 3D-HST data products.


Data products
=============

Region files
------------

Region files suitable for plotting with DS9 that define the layout of the  3D-HST primary WFC3/G141 and parallel ACS/G800L pointings are provided `here xxx <http://3dhst.research.yale.edu/Data.html>`_.

Field mosaic images
-------------------

The v0.5 release includes FITS mosaics of the F140W and G141 images covering the 3D-HST survey fields.  Each individual pointing is registered to an astrometric reference image by matching objects detected in the pointing and the reference with the IRAF task ``xyxymatch`` and shifts and a rotation necessary to align the frames are determined with the task ``geomap``.  A 1st-order polynomial background is subtracted from the F140W image after aggressively masking sources.  The G141 background is subtracted following the method described in detail by Brammer et al. (2012).  Briefly, we first divide the G141 images by the F140W *imaging* flat-field and then subtract a scaled average background image, chosen from a set that reasonably spans the variation in 2D structure of the background.

The final mosaics are combined with the MultiDrizzle software, with ``PIXFRAC=0.8`` and 0.06 arcsec output pixels.  The second FITS extension of the mosaics is the "inverse variance" weight map computed by MultiDrizzle.

+--------------------------------+-----------------------------------+
|  **Field**                     |  **WCS reference**                |
+--------------------------------+-----------------------------------+
| AEGIS                          |  EGS / ACS-F814W                  |
+--------------------------------+-----------------------------------+
| COSMOS                         |  CANDELS F160W                    |
+--------------------------------+-----------------------------------+
| GOODS-N                        |  GOODS ACS-F850LP                 |
+--------------------------------+-----------------------------------+
| GOODS-S                        |  CANDELS F160W, HUDF09 F160W      |
+--------------------------------+-----------------------------------+
| UDS                            |  CANDELS F160W                    |
+--------------------------------+-----------------------------------+


Field Browser
-------------

The main products of the v0.5 release are interfaces for browsing the full mosaics and for browsing the spectra in individual pointings as described below.  The interface uses the Google maps API to display the large FITS images and includes separate tabs to show the broad multi-wavelength properties of objects in the F140W images that have G141 spectra. 

In the maps, one can

- Pan and zoom with the navigation buttons at left.  Simply clicking and dragging the main map also pans, and double-clicking zooms in one level of the map.

- Toggle the overlay of the pointing footprints with the green square at right.  Turning the overlays off can help browser performance navigating the maps.

- Center the field on particular RA/Dec coordinates (must be in DD:MM:SS format).  The map recenters after values have been entered in the tabs at the top center of the display and the center position is briefly flashed by a yellow box.  The coordinate values themselves are updated after manually panning the map.

- Search online VizieR catalogs around the center position of the map (default 1" search radius) by clicking on the small magnifying glass icon next to the RA/Dec coordinate tab.

- Change the displayed image to different observed wavelengths (Xray, Radio, mid-IR, etc), depending on the availability in each field, by clicking on the labeled tabs at upper right.

 
Pointing Browser
----------------

Similar to the interface of the Field Browser, the Pointing Browser allows navigating a single 3D-HST pointing and finding detailed information about objects detected in that pointing, including their extracted 1D and 2D spectra and a redshift fit to the combined fit of the spectrum and ancillary matched photometry, when available.  The result of the fit is shown in the figure at lower right, and the redshift fitting is described in more detail by Brammer et al. (2012).  These figures are intended primarily to demonstrate how the spectral shapes are typically consistent with the broad-band SEDs and how the frequent emission lines of H-alpha and OIII can be robustly identified with the redshift constraints provided by the photometry.  Note that for the current browse pages, the spectral extraction has been limited to objects brighter than **H(140)<24**.

Most of the navigation features of the Field Browser are replicated for the individual pointings.  Rather than multi-wavelength images, however, tabs are provided for the F140W direct image, the G141 grism image, and the model image produced by the ``aXe`` software (see Brammer et al. 2012 for details on how the model is produced).  Emission lines show up as compact features in the G141 image.  However, some compact features are zeroth-order grism spectra, which can be identified in the model image.  

- Clicking on the circular object markers will load that object's information in the bottom panels.  This includes the ID and coordinates from the SExtractor catalog extracted for that pointing (the full catalog is linked at the top of the page).  Also shown is the F140W thumbnail, the (optimally) extracted 1D spectrum, and the 2D spectrum.  Clicking on first and last of these downloads individual FITS files of the thumbnails and 2D spectra, and the 1D spectrum image is linked to an ASCII version of the spectrum.  Tarfiles of all of the F140W thumbnails and 1D/2D spectra are linked at the top of the page.

- The right-most figure for each object demonstrates the redshift fit to the combined spectrum and photometry.  The left-most panel of that figure shows the full broad-band SED, while the middle panel shows the SED centered around the G141 spectrum.  The right-most panel shows the redshift probability distribution (yellow = photometry only, blue = spectrum only, purple = combined spectrum + photometry).  When a previously-measured spectroscopic redshift is available for a given object, it is indicated by a vertical green line.

- Clicking on the small "=" icon at upper left changes the layout of the display to show the products for a single object (default) or for 25 objects simultaneously, ordered by magnitude.  In both cases, clicking on the "ID" field for a given object recenters the map on that object.  In the multi-object view, one can page to the next set of brighter or fainter objects by clicking on the "+" and "-" buttons in the header of the ID column.


Contact
=======
For data questions, please contact Gabriel Brammer (gbrammer@eso.org), Ivelina Momcheva (ivelina.momcheva@yale.edu).  Additional general contact information can be found at http://3dhst.research.yale.edu/Team.html.


