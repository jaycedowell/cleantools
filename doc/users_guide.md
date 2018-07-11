A GRIDVIEW2/CLEANFLUX Session
=============================
In order to use the deconvolution routine on a data cube, the cube will need to be viewed in GRIDVIEW2.  GRIDVIEW2 is 
an expanded version of the original GRIDVIEW that includes all of the functionality of the original, plus the 
ability to deconvolve sources.  To start a GRIDVIEW2 session, run:

  ```@wasinit2
  @alfinit
  @cleaninit2
  restore, 'gridbf_1236+15a.sav'
  gridview2, grid, Cat='ex3d_1236+15a.cat'
  ```

The major differences between GRIDVIEW and GRIDVIEW2 are in the layout of the buttons near the bottom middle of the window.  
The buttons have been rearranged to make space for a new "CLEAN REGION" button that is responsible for running the 
deconvolution.  Only the "CLEAN REGION" button can be used to deal with the deconvolution;  "MEASURE FLUX" still calls GALFLUX
on the raw data.  

Selecting a region for deconvolution is analogous to measuring the flux.  All you need to do is click the "CLEAN REGION" 
button and draw a box around the source you want to clean.  You should draw a larger box than you normally would for measuring 
the flux so that the regions around the source can be cleaned and remove any side lobe artifacts from nearby sources.  The 
selected region  is further padded by 12’ on each side so that sidelobes down to the 5% level can be accounted for at the edge.  
Although this extra padding is included, it is “hidden” from view and is not displayed in the CLEANVIEW window.  In the event 
that the CLEANVIEW cannot properly pad the selected cleaning region due to a grid edge, an error will appear and the routine 
will immediately exit.  To clean these edge source, you may need to use the adjacent grid.

For further information about usage of GRIDVIEW/GRIDVIEW2 see the "Instruction on how to generate, baseline, and 'flatfield' a grid"
(http://caborojo.astro.cornell.edu/alfalfalog/idldocs/instruct_gridproc.txt), Section 4.1 for details.  Known issues and current 
limitations of GRIDVIEW2 are listed in Appendix A.

Once that region box has been draw, CLEANVIEW starts.  CLEANVIEW is the utility that handle the setup and control of the 
deconvolution.  The main parts of the routine are:
  1) Menus
    * The menus in CLEANVIEW are similar to those in GRIDVIEW2, i.e., 'Scaling' controls the scaling of the channel map 
      and 'Catalogs controls the overlay of the AGC, HVC, 3D, Checklist, and Continuum catalogs.  For convenience, the catalogs 
      selected in GRIDVIEW2 are automatically selected in CLEANVIEW.  The 'Settings' menu, however, behaves in a different manner
      than in GRIDVIEW2.  This menu is allows you to select between displaying the spectral and continuum maps and allows for the 
      setting of spatial and spectral smoothing parameters for the map being displayed.  There is also an additional 'Color' 
      menu where the display color map can be changed if desired.  The default color map is a blue-white linear stretch.
  2) Image Info
    * This panel display information about the location of the cursor on the current map, as well as the intensity 
      at that pixel.  Information about the channel number, velocity, and cleaning status are also displayed at the bottom of 
      the display window.
  3) CLEAN Control
    * This panel controls the cleaning process, i.e., how deep the each channel could be cleaned, what channel 
      range to clean, and if channels can be combined into wider pseudo-channels to improve the signal-to-noise.  The Grid RMS 
      value cannot be changed and is set at the start of the routine.  The other values can be changed.  Clicking on the "Combine
      to Improve S/N?" button changes the value from "No" to "Yes" and vice versa.  
    * The "Combine to Improve S/N?" button needs additional explanation.  On weaker sources the per channel signal-to-noise ratio is 
      too low for the cleaning to be effective, i.e., few channels have sufficient flux for the cleaning method to iterate.  In order
      to overcome this, this option can be enabled to combine adjacent channels together.  This is done as follows:
	      1)	The number of channels in the channel range is found and then factored.  These factors are then arranged in increasing order.
	      2)	The routine steps through each of the factors combining that number of channels.
	      3)	The stepping continues until 75% or more of the channels in current range have a peak flux above the flux limit.
      The deconvolution routine then cleans this new set of pseudo-channels in the same was as it does real channels.  Once the 
      deconvolution has completed, the pseudo-channels are split into the appropriate number of real channels.  Each of the split 
      channels associated with a pseudo-channel are identical.
  4) CLEAN Log
    * After the deconvolution has finished, this panel contains all of the information about how each channel was cleaned.  The top two labels (“Flux” and “Iteration”) give a summary of how many channels had an exit status of each type.  “Flux” means that the channel was cleaned to the flux limit set by max([Flux Limit, Grid RMS*Flux Sigma Limit]), while “Iteration” means that the channel was cleaned until the iteration limit was reached.  The text box below these labels stores detailed information about how each channel was clean.  If the “Combine to Improve S/N?” button is set to “Yes”, it also contains information about how many channel were combined to form the higher S/N pseudo-channels.  The details of reading this log file can be found in the reference guide to the CLN structure generated by CLEANFLUX.  A subset of this information can also be found in Section 1.2.6
  5) Buttons
    * The buttons located at the bottom of the window control various aspects of CLEANVIEW.  The buttons are:
      1) "View Beam" computes and launches the BEAMVIEW utility that displays the effective beam pattern as a function of position 
         in the data and allows you to save a JPEG of the central beam.  BEAMVIEW is a relatively straight forward program that 
         is documented in Appendix C.  
      2) "Clean Region" starts the deconvolution of the currently displayed region.  When this button is clicked, the values 
         stored in "CLEAN Control" panel are read in and the command line utility `alfa_clean7` is called to carry out the 
         deconvolution.  Details about `alfa_clean7` can be found in the guide to the general-purpose deconvolution utilities.  
      3) "Velocity Field" launches the VELOCITY FIELDVIEW window that allows you to get a quick look at the velocity field around 
         the object that you have drawn a box around.  The operation of this button in analogous to the "MEASURE FLUX" buttons in 
         both GRIDVIEW2 and CLEANFLUX.  This routine is not fully functional yet and should be used with caution.  What 
         documentation that exists on this routine can be found in Appendix D.  
      4) "Measure Flux" launches CLEANFLUX for measuring the flux of an object on the deconvolved map.  Its operation is analogous 
         to that of GAFLUX version 3 and is documented below.


Usage
-----
The window opens with a solid blue background and displays "Initializing..." on the map plot.  During this message, the sky variance 
of the data cube is being computed through an adaptive clipping method from the GSFC's IDL Astronomy User's Library and the beams 
needed to clean the selected region are generated.  Depending on the computer and the size of the selected region, this will take 
about 60 seconds.  Once this has completed, CLEANVIEW can be configured and the deconvolution process started.  The main controls 
are under the “CLEAN Control” panel (see Section 1.2.3).  The parameters flux limit, flux sigma limit, and iteration limit control
how deeply each channel is cleaned.  The default values for these parameters should be sufficient for cleaning most sources.  For 
particularly high-flux sources, i.e., ones with obvious sidelobes, the iteration limit may need to be increased to ensure that all 
channels are cleaned sufficiently.  Channel range and Combine to Improve S/N control which channels are cleaned.  Channel range 
needs to be a two-element comma-separated entry.  Only channels in this range are deconvolved; the remaining channels are left 
untouched. 

Deconvolution, like measuring fluxes, is an iterative process.  The routines should first be run with the default values here in 
order to gauge which values need to be tweaked in subsequent runs.  Once the initial run has been completed, the log is displayed 
in the "CLEAN log" panel.  To number to first look at is the total number of Iteration limited channels.  If there are a large 
number of iteration limited channels it may be best to increase the iteration limit and try again.  Another way to determine 
the quality of the cleaning is to used the “Velocity Field” button to generate a rough velocity field map.  This brings up the 
VELOCITY FIELDVIEW window that has options to define the rough spectral extent of the source.  However, this analysis method 
works best for only the brightest sources or sources that have relatively large widths.  To get a better idea of how the channels 
were cleaned, it can be informative to read through the log output.  Although this output is wordy and rather lengthy, it does 
provide the most complete view of the deconvolution process.  The structure of the log is fully documented in the section that 
describes the CLN structure.  For convenience, the most relevant portion of the documentation is reproduced above.

Once a flux box has been drawn in CLEANVIEW, CLEANFLUX is started and used to measure the source position, velocity, flux, 
etc.  CLEANFLUX is basically GALFLUX version 3 that has been modified to account for the fact that (1) the deconvolved data 
already has its polarizations combined and (2) that the clean beam has a different integrated area for each source.  See the
"Instruction on how to generate, baseline, and 'flatfield' a grid", Section 4.2 for details on the usage of CLEANFLUX/GALFLUX.  
If the region needs to be cleaned again, all you need to do with tweak the relevant parameters and push "Clean Region" again.  
This will automatically delete the previous results and restart the process. 

Information and the beam building process, known issues and limitations of the current version of CLEANVIEW can be found in 
Appendix B.

Brief Description of the CLEAN log format
-----------------------------------------
The most important sections of this log for determining the quality of the cleaning are in the format of:

```
...
>  Channel  378   
>  exiting on flux limit   
>   used 137 iterations  
>   peak at loop exit 16.411306 mJy   
>   residual RMS is 2.5879821 mJy
...
```

The first line gives the channel (or pseudo-channel if CLEANVIEW is configured to with the "Combine to improve S/N?" option 
is set to yes) number.  The second contains what condition caused the cleaning loop to exit on that channel (flux or iteration 
limit).  The next line gives the number of iterations used in cleaning.  For channels with no peaks about the flux density 
thresholds set by FLUX_LIMIT and SIGMA_list, this value is zero.  The fourth line gives the flux density of the last peak 
cleaned.  These two lines are particularly useful in determining to what depth each channel should be cleaned.  The final 
line gives the RMS of the remaining (un-cleaned) component of the channel.  For a well cleaned channel, this value should 
be comparable to the value stored in ROBUST_RMS.

Known Issues and Limitations of GRIDVIEW2
-----------------------------------------
Know Issues:

  * There are currently no known issues with this routine.

Limitations:

  * There are currently no limitations, relative to GRIDVIEW, with this routine.


The Beam Building Process, Known Issues, and Limitations of CLEANVIEW
---------------------------------------------------------------------
The Beam Building Process:

Although the effective beam pattern is position-dependent, CLEANVIEW uses a variety of approximations to reduce the 
computational load associated with generating a beam pattern for each pixel that needs to be cleaned.  These approximations are 
based on the fact that beam changes shape faster in dec. than it does in RA.  This allows CLEANVIEW to take a "smart" approach 
to modeling the beams.  First, beams for the entire dec. range are created for the RA value that lies at the center of the region 
to be cleaned.  Next, the weight map for the regions is examined to find pixels where the weight changes relative to the central 
strip of beams.  For each pixel where the changes is ≥10%, a new beam is computed.  Finally, the remaining area of the grid is 
filled in by copying the computed beams to areas that do not have a beam.  The copying is carried out from east to west.  

Know Issues:

  * There is an issue with the coordinate system the cursor sees when switching between CLEANVIEW and BEAMVIEW/VELOCITY FIELDVIEW.  
    This causes the coordinate system in the channel map window to be reset and the Image Info panel to not update.  To fix this, 
    simply left click once on the channel map window.  This will reset the coordinates and start the Image Info panel updating.
  * Altering the scaling/colorbar options before deconvolution causes CLEANVIEW to crash.  This problem happens because the arrays 
    that hold the cleaned map are initialized to zero at the start of CLEANVIEW.  Typing “retall” at the IDL command prompt should 
    clear this up.  
  
Limitations:

  * Not all of the catalogs listed in GRIDVIEW2 are currently available in CLEANVIEW.  The catalogs that are missing are 
    Multi-AGC and NVSS.
  * There is currently no option that allows you to jump to a specific channel.  The only way to move between channels is 
    by moving the channel slider.  
  * When using logarithmic or histogram equalization scaling, the colorbar is forced to automatic scaling.  Linear scaling is 
    currently the only option that supports constant (-5 to 10 mJy/beam) scaling.
  * The grid RMS computed and displayed by CLEANVIEW is not a rigourous RMS.  Rather it is the RMS associated with ~1,000,000 
    pixels (~1/16th of the total pixels) from the grid.d structure.  This method is used because it is much faster and less
    memory intensive than using the full grid.d array.  The two values are also in good agreement with each other, with a 
    difference less than 0.01 mJy/beam.
  * The file format saved using "File->Save Data" is currently only documented in the file data_structures.md file.
  * When the "Combine to Improve S/N?" option is selected, there may be an issue with using a channel range that specifies a 
    prime number of channels.  This is due to the way that channels are combined in only integer multiples.  Care should be 
    taken not to specify a channel range that has a prime number of channels.  Channel ranges that are multiples of two are 
    best.  For an alternative description of this particular limitation, see the backend_guide.md file.


BEAMVIEW
--------
Overiew:

BEAMVIEW is a utility for displaying the effective beam patterns computed by build_beam5 for the deconvolution.  The five main 
parts to this routine are:  the "Effective Beam Pattern" plot window, the "Unique Beam Regions" plot window, the "Detail" panel, 
the “Properties” panel, and the “Beam Modeling Log” panel.  The Detail panel gives information about the beam strength located 
under the current cursor position on the Effective Beam Pattern plot window while the Properties panel gives details about the 
beam as a whole.  This includes the location (RA/dec and pixel coordinates on the grid displayed in GRIDVIEW2), the fitted FWHP 
of the beam, and the position angle of the beam in degrees East of North.  The Unique Beam Regions plot window displays the area 
to be cleaned and plots each unique beam as a different color.  Clicking on this window updates the Effective Beam Pattern window 
with the beam used at that point.  The Beam Modeling Log panel shows the logs generated when the beams were modeled.  This 
information includes how many drifts were used to create the model, if the clean_tools.so library was used for the computations, 
the position error incurred by modeling the drifts on a discrete grid, and the total time required to model the beams.  

The BEAMVIEW menu bar as a variety of options for displaying and saving the beam pattern.  Under the "File" menu there is an 
option to export the central beam as a JPEG.  Under the "Color" menu, different color maps can be chosen along with the option 
to invert the colors.  The default is a rainbow+white color map.  The "Scaling" menu has options to change the stretch displayed
from linear to logarithmic and vice versa.

Know Issues:

  * It should be noted that after BEAMVIEW has been closed, it may be necessary to click the map window in CLEANVIEW to 
    re-initialize its coordinate system.  Otherwise, the Image Info panel will not update properly.


VELOCITY FIELDVIEW
------------------
Overview:

VELOCITY FIELDVIEW is a utility for creating a rough velocity field of an object.  It is invoked using the "Velocity Field"  
button in CLEANVIEW and drawing a box around the region to examine.  The main parts of this routine are:
  1) Menus:  VELOCITY FIELDVIEW has a variety of menu options for controlling the creation of the velocity field.  
       * "Color" is similar to the color menu in BEAMVIEW and allows for the color velocity field map to be changed.
       * "Contours" controls how many and of what type of contours are displayed.  By default it selects 10 contour 
         levels and overplots white lines at the boundaries of the contours.  The overplotting is changed by using the 
	 "Overplot Contours" option.
       * "Clipping" tells the routine which column density values should be displayed in the field map.  The default is 
         10^19 atoms/cm^2.

  2) Spectrum:  The spectrum consists of the upper three panels in VELOCITY FIELDVIEW.  The uppermost panel deals with adjusting 
     the scale displayed in the spectrum window.  This, combined with the right-click zoom, is useful for centering up on a 
     particular feature and adjusting the display range.  To apply the new range, simplify click "Rescale".
     
     The next lower panel displays the the spectrum of the center pixel of the source box.  This spectrum can be interacted with 
     using the left, middle, and right mouse buttons.  The left button is used in a manner similar to GALFLUX/CLEANFLUX to select 
     the region of interest in the spectrum.  The middle button resets the spectrum to its original size and clears all channel 
     markers.  The right button allows for zooming in.

    Below the spectrum is a panel that displays information about the cursors current location.  

  3) Buttons:  There are two buttons located at the bottom of the window, "Reset" and "Compute".  
       * "Reset" resets the spectrum window to display all channels at a flux density range of -10 to 30 mJy/beam.  It also 
         erases the velocity field map window. 
       * "Compute" is used to generate the velocity field map after the upper and lower channel ranges have been selected in 
         the spectrum window using the left mouse button.  

Usage:

Specifying the channel range to use in VELOCITY FIELDVIEW is similar to that in GALFLUX/CLEANFLUX;  you use the left mouse button
to define the range, and the right mouse button to zoom in.  Once the range has been picked, you need to click on Compute to 
generate and display the map.  If you need to redo your channel range selection you can either left click a third time to clear 
the marks (this maintains the current zoom), middle click to clear the marks and reset the view window, or click.

Once a map has been computed, the color, number and type of contours, and clipping levels can to altered to change the display.  
These options are accessible from the menu bar.  It is also possible to save a JPEG version of the velocity field map.

Know Issues:

  * When VELOCITY FIELDVIEW is first started, the coordinate system for the spectrum is not initialized, i.e., the cursor info 
    will not update.  To fix this,  simply click once on the spectrum window.  This should bring up the dashed cross-hairs and 
    begin updating the cursor info.  This is sometimes necessary after computing the field map as well.
  * If "Compute" is clicked before the spectral range is defined, the routine will crash.  To recover from this, type 'retall' at 
    the IDL command prompt.
  * Sometimes when changing clipping levels, the previous velocity field map will not be erased and the new map will be simply 
    overplotted.
  * It may be necessary to click the map window in CLEANVIEW to re-initialize its coordinate system after VELOCITY FIELDVIEW has 
    been exited.  Otherwise, CLEANVIEW's Image Info panel will not update properly.
