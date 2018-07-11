cleantools Backend Routines
===========================

Introduction
------------
In order to make the deconvolution routine as widely applicable as possible, the routines have been developed in a modular 
fashion.  The main GUI application for handling the cleaning, CLEANVIEW, does not handle the beam modeling or deconvolution on 
its own.  Rather, it calls the command line utilities, build_beam5, build_beam_grid2, and alfa_clean7, to do the computational 
"heavy lifting".  The purpose of this document is to detail the inner workings of these routines and explain how they can be 
integrated into other GUI applications or analysis routines.  This document also explains some of the associated routines used by 
build_beam5 and alfa_clean7.



build_beam5
-----------
Overview:

build_beam5 is the utility used to model the effective beam pattern found at a particular location in an ALALFA data cube.  This 
is the utility called by BEAMVIEW and CLEANVIEW that does the actual work.  

Methods:

build_beam5 works by first finding all of the drifts that contribute to that requested grid point and then simulating observations 
of a point source located at that point.  The simulated drifts use the az=180° drifts from the A1963 project (Irwin et al. 2008) that 
have been re-sampled to N' x N' @ 0.05'/px resolution and are carried out on a sky that is  N' x N' @ 0.05'/px.  The 0.05'/px 
resolution was chosen to minimize the error associated with mapping the declination of each drift onto the grid and to insure 
that 1 second in RA involved an integer number of pixels.  Values of N are 21, 25, 31, 35, 41, 45, and 51.  The wide variety map 
sizes available is motivated by limitations of the CLEAN method, namely that a NxN image needs a 2Nx2N beam to fully clean it.  
Map sizes larger than 51’ are not created because the 51’ map contains virtually all of the effective beam.

The details of the simulation of the effective beam at grid point (i,j) are:
  1)	The grid's grid_makeup structures are parsed to find all drifts that contribute to the point (i,j).
	2)	The simulated sky with central point source is created and the beams are allowed to drift across.  A record is created for every second in RA for each beam and each drift.  These drifts are fairly straight forward and are not aware of the two polarizations or missing data (via bad boxes set in the grid's pos structures).  They are, however, aware of the RA values where the individual drifts start and stop.
	3)	The simulated drifts are then mapped to a 1'/px grid similar to how grip_prep.pro makes grids.  For the weight function we read the grid.wf_fwhm tag and assume that a Gaussian weighting function is being used.  We also give beam 0 a higher (factor of 1.27) weight that comes from my reading of grid_prep.pro.
	4)	The resulting map is the beam associated with the map at that point.  This map is then normalized and returned.

The beams that build_beam5 produces agree well with isolated bright continuum sources in several grids, both in terms of FWHM 
sizes and the overall shape of the beam.

Usage:

To use build_beam5, type:

```
build_beam5, grid, [120, 98], beam, MapSize=MapSize, AntPat=AntPat, /CExts, /Silent, Output=Output
```

Inputs:

There are three main inputs to this utility, the first is a grid structure while the second is a two-elements integer array that 
specifies the (i,j) coordinates of the point in the grid where the beam needs to be modeled.  The final input, MapSize, is used 
to set the size of the beam map that is modeled.  The default value is 25 so as to be backward compatible with older versions of
build_beam.  AntPat is an optional variable that is used to pass a collection of the seven ALFA feed response patterns to 
build_beam5.  This is useful if build_beam5 is being used as part of an automated method to build beams that all have the same 
MapSize.  

There are also two optional keyword, SILENT and CEXTS.  SILENT can be set to turn off verbose output to the terminal.  CEXTS 
causes build_beam5 to use the external shared library clean_tools.so for all of the computations.  Since this library is written 
in C, it is significantly faster than the IDL version.  This library can optionally be complied with OpenMP to further increase 
its speed on multiple cpu/core systems.  Note:  To use this feature, the library functions need to be linked into IDL.  See 
below for instructions on how to do this.

Outputs:

There are two outputs from this utility.  The first is the beam pattern itself, which is a NxN element double array that has a 
1 arc minute pixel scale.  The other output is a ~10-element string array that holds the verbose output.  This array is created 
even if the SILENT keyword is set and holds information about the modeling process.  The format of the log is documented in the 
CLN structure documentation.

Known Issues:

  * build_beam5 currently assumes that the grid has a 1’ pixel scale.  It does not read the GRID.DELTADEC tag to find out what the 
    real pixel scale is.  Adapting the current routine for an arbitrary pixel scale would be difficult because of the 0.05’ pixel 
    scale of the ALFA beam maps imposed by the current version of a1963_beams.  
  * build_beam5 also assumes that the weighting function used to generate the grid is a Gaussian but it does use the GRID.WF_FWHM
    to get the correct FWHM.  However, my readings of grid_prep indicates that a Gaussian is currently the only option for a 
    weighting function.
  * build_beam5 also assumes that the data quality is uniform, both across and channels and polarization, and ignores flagging 
    information stored in the GRID structure.  
  * The routine can only handle drifts where the azimuth arm is at either ~0° or ~180°.  This is fine for drifts currently being 
    produced that are decidedly north or south of zenith at Arecibo.  For future drifts near zenith, other azimuth arm angles may 
    be used which will cause a problem.  The A1963 maps do not include beam maps at azimuth angles centered around ~90° or ~270°.

Limitations:

  * build_beam5 depends on seven files being located in the directory specified by the AGCDIR variable.  These files are 
    a1963_N_180.sav, where N is 21, 25, 31, 35, 41, 45, and 51 and contain the az=180° beams interpolated to a resolution of 
    0.05'/px from the data collected by the A1963  project (Irwin et al. 2008).  These files are generated by the routine 
    A1963_BEAMS which is documented below.
  * build_beam5 does not currently have an option to return higher resolution beam maps that may be of interest in studying the 
    effective beam pattern.  This limitation mainly exists because of the amount of computer time needed to generate these maps.
  * In order for the CEXTS keyword to be used, IDL needs to know how to interface with the clean_tools.so library.  To do this, 
    run:
    ```
    linkimage, 'build_beam_worker', '/path/to/clean_tools.so', 1, 'build_beam_worker', max_args=6, min_args=6
    linkimage, 'lib_version', '/path/to/clean_tools.so', 1, 'lib_version', max_args=1, min_args=0
    ```
    Information about the clean_tools.so library is detailed in Appendix C.


build_beam_grid2
----------------
Overview:

Since the effective beam pattern at every grid point in a data cube is different, a new beam needs to be computed for each point 
in the region that needs to be cleaned.  build_beam_grid2 serves as a way to automate build_beam5 into filling a region with 
beam patterns, compute the flux correction factor for each map point, and determine the best-match FWHM and PA for the Gaussian 
restore beam.

Methods:

The routine works as follow:
  1) Holder arrays are initialized for both beams and FluxCorr based on the values of llx, lly, urx, ury, and MapSize.
	2) The modeling keywords (DecOnly, Smart) are checked for.  If no keywords are found, full beam modeling is used.
	3) The modeling loop begins in one of three modes:
	     * DecOnly – Beams are calculated for all dec. values for the RA value specified by (urx+llx)/2.  The computed beams are 
         then copied to all other RA values.  This is the fastest method.
	     * Smart – Similar to DecOnly, but the weight map is also checked for regions where the weight changes by more than 10% 
         relative to where the DecOnly beams are computed.  This method is usually two to three times slower than DecOnly, but 
         still completes in a reasonable (<3 minutes) time.
	     * Full – All beams are computed.  This is the slowest but most exact method of modeling the beams needed for cleaning.
	4) After the beams are computed, the FWHM and PA is determined by fitting a 2-D Gaussian to the inner 11’x11’ of the central 
     beam located at ((urx+llx)/2, (ury+lly)/2).  After the FWHM is know, the FluxCorr correction factor is computed for every
     beam.  This factor is based off the area ratio of the dirty and clean beams.

Usage:

To use build_beam_grid2, type:

```
build_beam_grid2, grid, llx, lly, urx, ury, beams, MapSize=35 FluxCorr=FluxCorr, FWHM=RestoreFWHM, PA=RestorePA, /Silent, Output=output, /Smart, /CExts
```

Inputs:

There are a variety of inputs needed by build_beam_grid2.  These are:
  * grid - A grid structure.
  * llx - Lower left-hand x coordinate of the region where beams need to be built.
  * lly - Lower left-hand y coordinate of the region where beams need to be built.
  * urx - Upper right-hand x coordinate of the region where beams need to be built.
  * ury - Upper right-hand y coordinate of the region where beams need to be built.
  * MapSize - The spatial extent of each of the beams that are built.  Valid values are 21, 25, 31, 35, 41, 45, and 51 arc minutes
   on a side.  The default value is 25 ‘.
  * DecOnly  - Optional keyword to change how many beams are computed.  If both this and Smart are set, DecOnly is used.  If neither 
    is set, the full set of beams is computed.  See above for details.
  * Smart - Optional keyword to change how many beams are computed.  If neither this nor DecOnly is set, the full set of beams is 
    computed.  See above for details.
  * CExts - CEXTS causes the build_beam5 routine used by build_beam_grid2 to use the external shared library clean_tools.so for all 
    of the computations.  Since this library is written in C, it is significantly faster than the IDL version.  This library can 
    optionally be complied with OpenMP to further increase its speed on multiple cpu/core systems.  Note:  To use this, the library 
    functions need to be linked into IDL.  See below for how to do this.
  * Silent - Keyword to turn off information statements.  By not setting this keyword, the statements are printed to the terminal.

Outputs:

There are also a variety of outputs returned by build_beam_grid2.  These are:
  * beams - A 4-D array that contains the effective beam pattern at each point in the map.  The dimensions are map RA, map dec., 
    beam RA, and beam dec.  The sizes of first two dimensions are set by llx, lly, urx, and ury.  The sizes of the last two 
    dimensions are set by MapSize.
  * FluxCorr - A 2-D array that contains the correction factor needed to go from mJy/"dirty" beam to mJy/"clean" beam for each 
    point.  This correction is applied to the fluxes via 1/FluxCorr.
  * FWHM - The recommended FWHM for the Gaussian restore beam.  This is based off fitting a 2-D Gaussian to the central beam.  
  * PA - The recommended position angle, in degrees east of north, for the Gaussian restore beam.  This is based off fitting a 
    2-D Gaussian to the central beam.  
  * Output - All logging/information statements are saved to a string array.  This array is returned using this variable.  This 
    information is useful for interfacing with other programs.

Known Issues:

  * There are currently no know issues with build_beam_grid2.

Limitations:

  * Although setting the Smart keyword does something to try to estimate where the shape of the beam changes, there is not a clear 
   link between changes in the beam shape and changes in the weight map.  Or more specifically, not all ≥10% changes in the weight 
   map have different beam shapes.


alfa_clean7
-----------
Overview:

alfa_clean7 is the routine responsible for the deconvolution of the data cube and is the routine called by CLEANVIEW when the 
"Clean Region" button is pushed.  This routine also has an associated helper routine, combine_chans, that is documented in below.

Methods:

The procedure works as follows:
	1) The data cube is read into the program as are all of the CLEAN control options (flux limit, gain, etc.).
	2) The cleaning loop selects the maximum point in each map for cleaning.  If more than one point fits this selection criterion 
     then the one closest to the lower left-hand corner is used. The cleaning process removes a portion of the total flux from that 
     map that is set by the gain and proportional to the value of the peak.  The removed flux is added to the CLEANed map as a 
     delta-function.  The loop continues until:
	     * the RMS limit is hit (Sigma*RMS > map maximum point),
	     * the Flux limit is hit (Flux > map maximum point), or
	      * the iteration limit is hit.
     If more than one flux limit is defined, the larger limit is used.
	3) The CLEANed data the convolved with a Gaussian beam with a FWHM of 4.6' by 3.1' or whatever is specified by the FWHM flag.  
     This re-convolved map is then added back to the residual flux to create the output map.  

Usage:

To use alfa_clean7, type:

```
alfa_clean7, data_cube, beams, cleaned_map, Continuum=cont_map, ,C_Cleand=cont_map_cleaned, MapRMS=3.12, NIter=2000L, Flux=5.1, FWHM=3.0, Sigma=5.0, Gain=0.1, FWHM=[3.8,4.3], PA=-87.5, Range=[564,679] ,/AllowSum, Exit_Status=exit_status, /CExts, /Silent, Output=output
```

Inputs:

There are a variety of inputs that alfa_clean7 needs in order to clean the correct channels from the data cube and to the correct 
depth.  These are:
  * data_cube - Data cube.  This is a 3-D double array with dimensions [v, ra, dec].
  * beams - N' by N' at 1'/px resolution beams from build_beam5.pro for each pixel in the region to be cleaned.  This is a 4-D array 
    with dimensions [map RA, map dec., beam RA, beam dec.].
  * Continuum - Map of continuum sources for cleaning.  This is a 2-D double array with dimensions [ra, dec].
  * MapRMS - RMS of the grid that the map is from.  In CLEANVIEW, this is computed using the robust_sigma routine from the IDL 
    Astronomy Users' Library.
  * Niter - Maximum number of iterations in the CLEAN loop.  The default is 2,000.  This should be a long integer to avoid overflowing
    the integer limit.
  * Flux - The flux level, in mJy/beam, to clean down to.
  * Sigma - The sigma level with which to clean down to.  This is multiplied by MapRMS to convert to a flux.
  * Range - A two elements integer or long integer array that specifies the channel range to clean.  Note:  This range is inclusive.
  * Gain - The gain used for the CLEAN loop.  The default value is 0.05.
  * AllowSum - Keyword set to combine channels to improve the per-channel signal-to-noise ratio.  For details of how the channel 
    combining process works, see below.
  * FWHM - The FWHM of the Gaussian restore beam.  The default value is 4.6' by 3.1'.
  * PA - The position angle of the Gaussian restore beam.  The default value is 0 degrees.
  * CExts - CEXTS causes alfa_clean7 to use the external shared library clean_tools.so for all of the computations.  Since this 
    library is written in C, it is significantly faster than the IDL version.  This library can optionally be complied with OpenMP 
    to further increase its speed on multiple cpu/core systems.  Note:  To use this, the library functions need to be linked into 
    IDL.  See below for how to do this.
  * Silent - Keyword to turn off information statements.  By not setting this keyword, the statements are printed to the terminal.

Ouputs:

There are also a variety of outputs that alfa_clean7 generates that store not only the deconvolved data but also information about 
the deconvolution process.  These output variables are:
  * cleaned_map - The CLEANed version of data_cube.  The dimensions are the same as those of data_cube.
  * c_cleaned - The CLEANed version of the continuum map.  The dimensions are the same as those of Continuum.
  * Exit_Status - Numerical exist status for each channel’s cleaning loop.  The codes are:
      * -1:  Exit caused by the iteration limit.
      * +1:  Exit caused by the RMS limit.
  * Output - All logging/information statements are saved to a string array.  This array is returned using this variable.  This 
    information is useful for interfacing with other programs.

Known Issues:

  * There are currently no known issues with alfa_clean7.

Limitations:  

  * When the AllowSum option is used, a limitation is encountered with the channel range contains a prime number of channels.  For 
    a more detailed description of the problem, see below.
  * In order for the CEXTS keyword to be used, IDL needs to know how to interface with the clean_tools.so library.  To do this, run:
    ```
    linkimage, 'alfa_clean_worker', '/path/to/clean_tools.so', 1, 'alfa_clean_worker', max_args=13, min_args=13
    linkimage, 'lib_version', '/path/to/clean_tools.so', 1, 'lib_version', max_args=1, min_args=0
    ```
    Information about the clean_tools.so library is detailed below.


A1963_BEAMS
-----------
Overview:

A1963_BEAMS combines and interpolates the beam maps from the A1963 project (Irwin et al. 2008) to generate beam maps for azimuths 
between 104° and 255°.  If the IDL function “reverse” is used, these range of azimuths can be extended to include 284° to 75°.  

Methods:

The procedure works are follows:
  1) The input azimuth is compared with the azimuths mapped by the A1963 project.  The two azimuths from the project that bracket 
     the input azimuth are then used to create the new maps.  The FITS files corresponding to the two azimuths are loaded into 
     memory.  
  2) The beams are re-sampled to the correct spatial resolution using a bicubic interpolation routine that is included with the 
     a1963_beams.pro file.  The new resolution is 0.05'/px.
  3) Once the interpolation is complete, the beams are averaged together using a linear combination of the two azimuths.
  4) This procedure is then repeated for the remaining six ALFA beams.  The results are saved to a 3-D double array where the first 
     dimension indexes across the ALFA beam numbers.

Usage:

To run a1963_beams, type:

```
a1963_beams, final_beams, Width=25, /Recenter, Azimuth=185, File=’out.sav’, /Verbose
```

Inputs:

There are five inputs to this routine.  Azimuth is used to specify the azimuth where the maps should be computed.  Width sets the 
spatial extent of the beams in arc minutes.  Recenter is a keyword that independently recenters each beam each beam after the 
interpolation step.  If this keyword is set, the maximum value in the map is shifted to be located at the center of the array.  If 
this keyword is not set, the original centers are used.  File can be used as either a keyword or a variable.  When it is used as 
a keyword, the beams are saved to an IDL sav file named “a1963_<width in arc minutes>_<azimuth in degree>.sav”.  When it is used 
as a variable, the value is used as the file name.  The Verbose keyword turns on information statements and creates plots of the 
individual interpolated beams as the routine runs.

Once the beams have been re-sampled, the user is prompted for whether or not they wish to inspect and clean the maps.  If yes, then 
each beam map is displayed with a logarithmic stretch.  Cleaning is accomplished by using the mouse to select polygonal regions that 
are likely to be artifacts and setting the region inside the polygon to zero.  A left click adds a new vertex to the current polygon, 
while a right click adds a new vertex and closes the current polygon.  Additional polygons can be added by using another left click.  
Once all of the artifacts have been masked, a middle click applies the mask and moves to the next beam.

Outputs:

There is only one output to this routine.  It is a 3-D double array with dimension [7, (Width/0.05), (Width/0.05)].  The first 
dimension indexes the ALFA beam number.


Known Issues:
  * There are currently no known issues with this utility.

Limitations:

  * The routines depends on the FITS files being located in a sub-directory of the current directory called A1963 and having the 
    same names as when they were downloaded.  About 90MB of FITS files from http://www.astro.queensu.ca/~irwin/pub/ngc2903/ are 
    need to cover the full azimuth range.  
  * This routine is not currently distributed with the other cleaning tools.

COMBINE_CHANS
-------------
Overview:

combine_chans is the routine responsible for combining channels to improve the per-channel signal-to-noise ratio.  This is the 
command used when the AllowSum keyword is set in alfa_clean7.  

Methods:

This routine works as follows:
  1) The number of channels in the channel range is found and then factored.  These factors are then arranged in increasing order.
  2) The routine steps through each of the factors combining that number of channels.
  3) The stepping continues until 75% or more of the channels in current range have a peak flux above the flux limit.
  4) Once the deconvolution is complete, the pseudo-channels are split into the correct number of real channels.  

Usage:

combine_chans should only be called from inside alfa_clean7.  In light of this, no further documentation on its usage is given here.

Known Issues:  

  * There are currently no known issues with this routine.

Limitations:
  * Since combine_chans combines channels in integer multiples, channel ranges with a prime number of channels pose a problem. 
    For example, if your channel range has 50 channels, this factors into:  1, 2, 5, 10, 25, and 50.  If the channel range is 53, 
    however, this factors into: 1 and 53.  Thus, combining the channels may result in the entire range being put into only one 
    pseudo-channel.  Care should be taken to not to set the number of channels to a prime or a number that has only a few factors.


The clean_tools.so Library
--------------------------
Overview:

The clean_tools.so library contains versions of the computational sections of build_beam5 and alfa_clean7 that have been 
reimplemented in C to improve their execution speed and take advantage of multiple cpu/core systems via parallelization with OpenMP.  
In many cases, the speed up of switching to C alone is up to a factor of 10.  Parallelization further increases this margin, 
particularly in the case of the deconvolution routine.

Methods:

There are three functions currently included in the library:  build_beam_worker, alfa_clean_worker, and lib_version.  The purposes 
of the first two are self explanatory.  The final routine, lib_version, is used simply to report the version of the library that 
is being used.  As of this writing, the library version is 20080730.

Usage:

With the exception of lib_version, the routines found in the clean_tools.so library are not intended to be used outside of 
build_beam5 and alfa_clean7.  

Known Issues:  

  * In order for OpenMP support to be enabled, clean_tools.so needs to be complied with a complier that has an OpenMP library that 
    can be dlopen()ed and supports the -openmp complier flag.  GCC versions ≥4.3.0 meet these requirements as does the Intel C++ 
    compiler version 10.1.  
  * When complied with OpenMP support and static scheduling, not all of the calls that build_beam5 makes to build_beam_worker 
    produce valid results.  The reason for this is not understood and there  is nothing in common between failed calls, i.e., a 
    failed setup can be run again and finish successfully.  In the current version of build_beam_worker avoids this problem by 
    setting the scheduling to dynamic and using a fixed chunk size of 200.  

Limitations:  

  * The beam modeling routine in clean_tools.so, build_beam_worker, is not as efficient as it could be.  On a dual core machine 
    with two threads, the peak processor usage is only ~150% with a mean ~130%.  This is likely due to issues that arise from 
    the relatively small computational load.
  * The deconvolution routine in clean_tools.so, alfa_clean_worker, does not handle the continuum map at the same time as the 
    spectral data.  This leads to a slight loss in efficiency.  This loss, however, is not so great as to warrant not using the 
    clean_tool.so library.
