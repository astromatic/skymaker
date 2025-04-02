# SkyMaker

a program that simulates astronomical images.

Check out the [official web page].

[official web page]: https://astromatic.net/software/skymaker

The general SYNTAX is similar to that of SExtractor:

> sky [<list_file>] [-c <Configuration_file>] [-<keyword> <value>] ...
> to dump a default configuration file: sky -d 
> to dump a default extended configuration file: sky -dd 

- A list file is an ASCII file containing a list of objects that can be added to
  the simulated image. An example is provided in the config/sample.list file.
  Note that only stars (code = 100) galaxies (code = 200), and image rasters
  (code = 300) are recognized
  in this version.
- Keyword parameters given in the command line override those from the
  configuration file.
- If the list-file is given as unique argument, SkyMaker searches for a
  default configuration file called "sky.conf".
- SkyMaker creates 2 files in output: the image itself, and a catalog containing
  the objects it contains (with name toto.list if IMAGE_NAME was set to
  toto.fits).
- Currently, the following TYPEs can be used with the IMAGE_TYPE keyword:
  PUPIL_REAL, PUPIL_IMAGINARY, PUPIL_MODULUS, PUPIL_PHASE, PUPIL_MTF,
  PSF_MTF, PSF_FULLRES, PSF_FINALRES, SKY_NONOISE and SKY. The two latter
  keywords should be used for creating actual instrument images.
- A FITS header (any FITS image, or even an ASCII dump) can be provided through
  the IMAGE_HEADER keyword: simply replace "INTERNAL" by the file name. SkyMaker
  will then make a copy of this header for the simulated image, enabling the
  latter to be easily processed through your usual reduction tools.
- Thanks to Pascal Fouque, parameters describing common optical
  aberrations (including defocus, spheric, astigmatism and coma) have been
  included in the description of the pupil phase-plane. Their normalisation
  follow the ESO convention (equivalent angular diameter of a circle, in the
  focal plane, which encloses 80% of the PSF flux; this is generally slightly
  more than the FWHM). However the user is invited to check this normalisation,
  and report any unexpected result.
- If a SEED_* parameter is set to 0, the corresponding random generator is
  initialized to a "random" (function of time) value.
- Beware of large AUREOLE_RADIUS values: during the calculation of the image,
  a temporary border of <AUREOLE_RADIUS> pixels in thickness is added all
  around the image, and can significantly affect the computation time.

