;+
; NAME:
;      GRIDView
; PURPOSE:
;       Visualization tool for astronomical 3D data cubes.
;
; EXPLANATION:
;
;       GRIDView displays output from the grid_prep
;       procedure.  It shows the 2-D view of the RA/Dec plane 
;       and allows the user to move through the cube in cz, changing
;       polarization, etc.  The user can also switch to a continuum
;       view of data in the grid structure  A clickable color bar can scale the
;       intensity/contrast of the image.  A weights map shows the
;       spectral/continuum weights.  The GRIDinfo and controls are
;       found in the lower left of the screen.  They include sliders
;       to move the spectral cube and choose the level of boxcar
;       smoothing.  The user can use the keyboard, NEXT/PREVIOUS
;       buttons, or arrow keys to "move" through the grid.  The
;       default is for single channel display.  Selected the
;       integration button will sum plus/minus the number of channels
;       specified.  Plus/minus 5 channels is recommended as a given
;       default.  In spectral mode, the AGC can be overlaid on the
;       map, and clicking on the map will give a DSS 2 blue image of
;       specified size in the lower right corner.  In continuum mode,
;       the NVSS will be overlaid, color-coded by intensity, and 
;       clicking on the map will give an NVSS image of specified size
;       in the lower right corner.  The user can customize the color
;       scheme under the settings menu, and scalling (linear, log,
;       hist_eq) is chosen under the scaling menu (linear is default).
;       Log scaling is determined by calculating an offset and
;       subtracting before taking the log of the pixels.  Passing a
;       keyword with a filename for a 3D extracted catalog will give
;       the user the option of overlaying the catalog contents.
;
;
; CALLING SEQUENCE:
;       gridview2, grid
;
; INPUTS:
;       grid - grid structure as output by grid_prep
;
;
; OPTIONAL INPUT:
;
;
; OPTIONAL INPUT KEYWORD:
;
;    cat3D=filenamestring - allows user to overlay a 3D extractor catalog
;          
;
; OUTPUTS:
;       none
;
;
; RESTRICTIONS:
;
;        Must have output from grid_prep
;
; EXAMPLE:
;
;       Pretty straight forward:  gridview2, grid, cat3D='my_awesome_galaxies.cat'
;
; PROCEDURES USED:
;         DSSQUERY, GSFC IDL USER'S LIBRARY
;         HOR.  See pjp AO library
;         VER.  See pjp AO library
;         TVCIRCLE
;         ADSTRING
;         TVBOXBK
;         NEDQUERY
;         RATICKNAME
;         DECTICKNAME
;         SIGFIGS
;         DECIMALS
;         GETDSS
;         
;
;
;
; MODIFICATION HISTORY:
;       WRITTEN, October 2, 2005
;       October  26, 2005, BK - weights window added, GUI
;                                     reorganized
;       November  3, 2005, BK - spectra display option
;       November  8, 2005. BK - smoothing and scaling options
;                                     added.
;       November 26, 2005. BK - Colorbar autoscaling
;       November 27, 2005. BK - weights map doesn't scale
;                                  - Simple spectrum display
;                                  - cube pixel values shown
;                                  - allowed constant value scaling
;                                    option
;       November 28, 2005. BK - Zoom function added to spectrum
;                                    display.  Aperture option
;                                    removed.
;       November 29, 2005. BK - Zoom function added to main
;                                    window.  Spectral weights added
;                                    to spectrum display
;       November 30, 2005. BK - AGC info display fixed.
;       December  1, 2005. BK - New version distributed.
;       December  2, 2005. BK - Added NED Search Query
;       December 13, 2005. BK - custom catalog overlay
;                                  - cz defaults to 21 cm
;                                    measurements and widths,
;                                    otherwise then goes to optical.
;       December 21, 2005. BK - Detcode added to AGC display
;        January 07, 2006. BK - fixed AGC (no cz) option -
;                                     thanks to MPH!
;        January 08, 2006. BK - added Sloan to cutout optical
;                                    image option.  Access through
;                                    spawn and wget. 
;        January 15, 2006. BK - Added AGC-multi mode - display
;                                    all AGC galaxies with different
;                                    colors based on current cz
;                                    status.
;        January 20, 2006. BK - Fixed coordinate display problem
;                                    with new procedure.  Clicking on
;                                    optical image display opens a
;                                    larger view via getdss.pro
;
;        February 14,2006. BK - Added xs,ys,xregion,yregion to
;                                     state structure to prevent
;                                     coordinate reset when opening
;                                     auxillary windows.
;
;                                  - Zoom featured fixed - cube
;                                    coordinates now correct
;
;        February 28,2006. BK - Added button for measuring flux
;           March 08,2006. BK - Added 'g' keyboard for getting
;                                    DSS images
;           March 10,2006. BK - Added check boxes for currently
;                                    selected menu items
;           March 11,2006. BK - Added checklist option for
;                                    measuring fluxes, so one can
;                                    remember what was measured
;           March 12,2006. BK - Fixed zoombox range error for
;                                    spectrum plotting window
;
;           June  01,2006. BK - Simple error checking if user enters
;                               bad 3D cat file name
;
;
;
;           July  16,2006  BK/AS - added Strong continuum sources (>
;                                  250 mJy/beam from NVSS)
;
;           July  25,2006  BK   - Rearranged GUI to fit on smaller
;                                 screens. Functionality remains the
;                                 same
;
;           August 8,2006  BK   - Spectral weights now displayed at
;                                 constant scale ALWAYS.  Update
;                                 released.
;
;           Sept. 23,2006  BK   - modified weights scaling (larger
;                                 range)
;
;           Dec.  12,2006  BK   - Added postscript output
;
;           Jan.  06,2007  BK   - Can view pixel spectra with both pols
;
;----------------------------------------------------------

pro gridview2_initgrid

common gridstate, grid

end


;-------------------------------------------
; Initialize common block
pro gridview2_initcommon, cat3D=cat3D

common gridview2_state, state
common gridstate

;grid=congrid(grid,420,420,735)

;gridswitch, reform(grid.d[*,0,*,*]), gridout

;gridout=congrid(gridout,420,420,735)

ramin=grid.ramin
ramax=grid.ramin+grid.(6)*(grid.deltara/3600.0)
decmin=grid.decmin
decmax=grid.decmin+grid.(7)*(grid.deltadec/60.0)

;racenter=ramin+(ramax-ramin)/2.0
;deccenter=decmin+(decmax-decmin)/2.0

;print, racenter
;print, deccenter

;print, double(decmax-decmin)*60.

;print, 'Obtaining Optical image.  Please wait...'
;queryDSS, [racenter*15.0,deccenter], image, header, imsize=double(decmax-decmin)*60., survey='2b'

;restore, 'imagetest.sav'

;image=randomn(1,1000,1000)



;optimagesize=size(image)
;image=congrid(image, optimagesize[1], optimagesize[1])

;raimagearr=(dindgen(optimagesize[1])/(optimagesize[1]-1))*(ramax-ramin)+ramin
;decimagearr=(dindgen(optimagesize[1])/(optimagesize[1]-1))*(decmax-decmin)+decmin

;raimagearr=(dindgen(optimagesize[1])/(optimagesize[1]-1))*(2.0)+12.0
;decimagearr=(dindgen(optimagesize[1])/(optimagesize[1]-1))*(2.0)+9.0


;raimagearr=reverse(raimagearr)

;raarr=findgen(420)/419.0*(ramax-ramin)+ramin
;raarr=reverse(raarr)
;decarr=findgen(420)/419.0*(decmax-decmin)+decmin
;czarr=findgen(735)/734.0*(grid.velarr[0]-grid.velarr[n_elements(grid.velarr)-1])+grid.velarr[n_elements(grid.velarr)-1]
;czarr=reverse(czarr)

imagesize=string(15.0*float(ramax-ramin), format='(f5.1)')+' X '+$
        string(float(decmax-decmin), format='(f5.1)')+' degrees'
velrange=strcompress(float(grid.velarr[n_elements(grid.velarr)-1]), /remove_all)+ ' to '+$
 strcompress(float(grid.velarr[0]), /remove_all)+' km/s'


;restore, '/home/dorado4/bkent/a2010/virgo3D/table3D.sav'
;restore, '/home/caborojo4/bkent/a2010/virgosouth/table2D.sav'

;Open catalogs
common agcshare, agcdir

;If passed as a keyword, open and rearrange the input 3D catalog

if (n_elements(cat3D) ne 0) then begin
  ;print, '3D catalog!'
  restore, agcdir+'catalog_template.sav' 
  table3D=read_ascii(cat3D, template=catalog_template)
  ;convert3Dcat, catin, table3D
  extracats=1
endif else begin
  ;opens catalog_template and dummy file
  restore, agcdir+'catalog_template.sav'
  extracats=0
endelse

;restore, '/home/dorado4/bkent/a2010/virgo3D/table3D.sav'
;extracats=1



restore, agcdir+'agcnorth.sav'
restore, agcdir+'nvsscat_alfalfa.sav'
restore, agcdir+'dhbb.sav'

;Shrink down the NVSS - only save what is needed for this grid
 invss=where(nvsscat.ra_2000_ lt ramax*15.0 AND $
             nvsscat.ra_2000_ gt ramin*15.0 AND $
             nvsscat.dec_2000_ lt decmax AND $
             nvsscat.dec_2000_ gt decmin)



agc=temporary(agcnorth)
rec=n_elements(agc.agcnumber) 
rahr=dblarr(rec) 
decdeg=dblarr(rec) 

;Reverse velocity correction of de Heij et al.
;Solar apex information - all from J2000 coords

            vapex=19.5  ;km/s
            lapex=56.157337   ;Galactic Long
            bapex=22.765049   ;Galactic Lat

            glactc, double(dhbb._raj2000/15.0),double(dhbb._dej2000), 2000, l,b,1
            onedeg=3.14159265359/180.0
            cos_p = sin(b*onedeg)*sin(bapex*onedeg) + cos(b*onedeg)*cos(bapex*onedeg)*cos(l*onedeg-lapex*onedeg)
            vcorrect=vapex*cos_p
            dhbb_vhelio=dhbb.lrv+vcorrect   ;To go from LSR to HELIO  ADD!!!!

rahr=double(agc.rah)+double(agc.ram)/double(60.0)+double(agc.ras10)/(10.0D)/3600.0D
decdeg=abs(double(agc.decd))+double(agc.decm)/double(60.0)+double(agc.decs)/3600.0D
signs=where(agc.sign eq '-')
if (signs[0] ne -1) then decdeg[signs]=-decdeg[signs]


state={baseID:0L, $                           ; ID of top base
       fmenu:0L, $                            ; ID of File menu
       helpmenu:0L, $                         ; ID of Help menu
       settingsmenu:0L, $                     ; ID of Settings menu
       scalingmenu:0L, $                      ; ID of image scaling menu
       plotwindowone:0L, $                    ; ID of main window where image is shown
       plotwindowtwo:0L, $                    ; ID of window for weights display
       plotwindowcolorbar:0L, $               ; ID of window for colorbar
       plotwindowimage:0L, $                  ; ID of DSS/NVSS image display
       button_linear:0L, $
       button_log:0L, $
       button_histeq:0L, $
       button_autoscale:0L, $
       button_constantscale:0L, $
       ;grid:grid, $                           ; Structure for the grid (output of grid_prep)
       ramin:ramin, $
       ramax:ramax, $
       decmin:decmin, $
       decmax:decmax, $
       xcubemin:0, $
       xcubemax:n_elements(grid.d[0,0,*,0])-1, $
       ycubemin:0, $
       ycubemax:n_elements(grid.d[0,0,0,*])-1, $
       current_xcube:0L, $
       current_ycube:0L, $
       xval:0L, $                             ; ID for RA label
       yval:0L, $                             ; ID for Dec label
       zval:0L, $                             ; ID for cz label
       xpix:0L, $                             ; ID for x pixel label
       ypix:0L, $                             ; ID for y pixel label
       intensitylevel:0L, $                   ; ID for pixel intensity level
       mousestatus:0L, $                      ; mousestatus for colorbar control
       imagesize:imagesize, $                 ; String for GRIDinfo display for image size
       velrange:velrange, $                   ; String for GRIDinfo display for velocity range
       projection:'Projection: '+grid.map_projection, $
       epoch:'Epoch: '+strcompress(grid.epoch, /remove_all), $
       dataset:'GRID: '+grid.name, $
       velslider:0L, $
       xs:dblarr(2), $             ;Added BK Feb 2005 to prevent coord conversion from having probs
       xregion:dblarr(2), $        ;Added BK Feb 2005 to prevent coord conversion from having probs
       ys:dblarr(2), $             ;Added BK Feb 2005 to prevent coord conversion from having probs
       yregion:dblarr(2), $        ;Added BK Feb 2005 to prevent coord conversion from having probs
       smoothslider:0L, $
       gotochannel:0L, $
       colortable:1, $
       polbutton:lonarr(3), $
       currentpol:2, $
       polAcolor:3, $
       polBcolor:8, $
       avgpolcolor:1, $
       polAselect:0L, $
       polBselect:0L, $
       avgpolselect:0L, $
       velintegselect:0L, $
       velintegselectstatus:0, $
       velintegtext:0L, $
       px:fltarr(2), $
       py:fltarr(2), $
       sx:0.0, $
       sy:0.0, $
       multiagcstatus:0, $
       agcoverlaystatus:0, $ 
       agcoverlaystatus_nocz:0,$
       overlaystatus3D:0,$
       overlaystatus2D:0,$
       dhbboverlaystatus:0, $
       nvssoverlaystatus:0, $
       agc:temporary(agc), $
       rahr:temporary(rahr), $
       decdeg:temporary(decdeg), $
       currentagc:lonarr(3,5000), $
       currentagc_nocz:lonarr(3,5000), $
       currentdhbb:lonarr(3, 1000), $
       agcinfo:'', $
       nvsscat:nvsscat[invss], $
       dhbb:temporary(dhbb), $
       dhbb_vhelio:dhbb_vhelio, $
       exportstart:0L, $
       exportnumframes:0L, $
       exportnumsteps:0L, $
       exportdirectory:0L, $
       apertureFWHM:0L, $
       plotwindowspectrum:0L, $
       plotwindowspectrumweights:0L, $
       spectrumwindowwidth:0L, $
       getspectrumstatus:0, $
       fluxmeasurestatus:0, $
       zoomstatus:1, $    ;Zoom is on by default
       spectrumaxisbutton:lonarr(2), $
       spectrumaxisstatus:0, $    ;0 for velocity, and 1 for channels
       spectrum:dblarr(n_elements(grid.velarr)), $
       spectrumweights:dblarr(n_elements(grid.velarr)), $
       maxweight:median(grid.w), $
       spectrum_xmin:0, $
       spectrum_xmax:0, $
       spectrum_ymin:0, $
       spectrum_ymax:0, $
       spectrumon:0, $
       stretchlower:0, $
       stretchupper:200, $
       minval:0.0D, $
       maxval:0.0D, $
       colorbar_y:0, $
       opticalsize:0L, $
       imageoptions:[' DSS ', 'Sloan'], $
       currentimage:0, $
       currentimage_rahr:0.0D, $
       currentimage_decdeg:0.0D, $
       dssimage:dblarr(210,210), $
       sloanimage:dblarr(3,210,210), $
       mapchoice:lonarr(2), $
       mapcurrent:'spectral', $ 
       imagelabel:0L, $
       copyright:0L, $
       keystatus:0, $
       currentscaling:'linearscaling', $
       colorbarscaling:'colorbar_autoscale', $
       extracats:extracats, $
       table3D:table3D, $
       colorbarmin:-5.0, $     ;Default minimum
       colorbarmax:10.0, $     ;Default maximum
       colorbarmin_text:0L, $
       colorbarmax_text:0L, $
         wid:0L, $          ; The window index number.
         drawID:0L, $    ; The draw widget identifier.
         pixID:-1, $         ; The pixmap identifier (undetermined now).
         xsize:0, $      ; The X size of the graphics window.
         ysize:0, $      ; The Y size of the graphics window.
         zoomsx:-1, $            ; The X value of the static corner of the box.
         zoomsy:-1, $            ; The Y value of the static corner of the box.
         boxColor:0L, $ ; The rubberband box color.}
       checkliststatus:0L, $
       cntsrcstatus:0L, $
       overlayboxbutton:0L, $
       overlayboxstatus:0, $
       overlaybox:dblarr(10,5), $
       ps_centerpixel_x:0L, $
       ps_centerpixel_y:0L, $
       ps_pixelrange:0L, $
       ps_centerchannel:0L, $
       ps_channelrange:0L, $
       ps_checkbox:0, $
       drift_draw:0L, $
       posfilelist:strarr(10000L), $
       driftfiles:strarr(n_elements(grid.grid_makeup[0].driftname)), $
       filelist:strarr(n_elements(grid.grid_makeup[0].driftname)), $
       makeupindex:0L, $
       listnumber:-1}


    ;   table2D:table2D, $
    ;   table3D:table3D}

;Remove unneeded variables
delvarx, agcnorth, agc, nvsscat, table3D, cattable, table2D, catalog_template


end






;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;STARTUP;;;;;;;;;;;;;
pro gridview2_startup, cat3D=cat3D

if (not (xregistered('gridview2', /noshow))) then begin

;print, 'Gridview startup'

common gridview2_state
common gridstate

gridview2_initcommon, cat3D=cat3D

set_plot, 'x'

loadct, 1, /silent
stretch, state.stretchlower,state.stretchupper

;Reset plots
hor
ver
!p.multi=0
device, decomposed=0

;Create widgets

;Top level base
tlb=widget_base(/column, title='GRIDview2', $
                tlb_frame_attr=1, xsize=1065, ysize=790, mbar=top_menu, $
                uvalue='tlb')

;Menu options
state.fmenu=widget_button(top_menu, value=' File ')
 
  buttonjpegoutput=widget_button(state.fmenu, value=' Export JPEG ', event_pro='gridview2_jpeg_output')
 ;; buttonpsoutput=widget_button(state.fmenu, value=' Export Single Postcript ', event_pro='gridview2_ps_output_gui')
  buttonexit=widget_button(state.fmenu, value=' Exit ', /separator, $
                           event_pro='gridview2_exit')

state.settingsmenu=widget_button(top_Menu, value=' Settings ')
  buttonpolsettings=widget_button(state.settingsmenu, value=' Pol Colors ', uvalue='polcolors', $
                           event_pro='gridview2_polsettings')
  buttonscalesettings=widget_button(state.settingsmenu, value=' Colorbar scale ', uvalue='scalesettings', $
                           event_pro='gridview2_scalesettings')
  state.overlayboxbutton=widget_button(state.settingsmenu, value= ' Add box marker ', uvalue='overlaybox', /Checked_Menu)
  clearboxmarkers=widget_buttoN(state.settingsmenu, value= ' Clear box markers ', uvalue='clearoverlaybox')

state.scalingmenu=widget_button(top_menu, value=' Scaling ')
  state.button_linear=widget_button(state.scalingmenu, value=' Linear ', uvalue='linearscaling',/Checked_Menu) 
  state.button_log=widget_button(state.scalingmenu, value=' Logarithmic ', uvalue='logscaling',/Checked_Menu) 
  state.button_histeq=widget_button(state.scalingmenu, value=' Histogram EQ ', uvalue='histeqscaling',/Checked_Menu) 

widget_control, state.button_linear, set_button=1

colorbarmenu=widget_button(top_menu, value=' Colorbar ')
  state.button_autoscale=widget_button(colorbarmenu, value=' Autoscale ', uvalue='colorbar_autoscale',/Checked_Menu)
  state.button_constantscale=widget_button(colorbarmenu, value=' Constant scale ', uvalue='colorbar_constantscale',/Checked_Menu)

widget_control, state.button_autoscale, set_button=1

catalogmenu=widget_button(top_menu, value=' Catalogs ')
  buttonmultiagc=widget_button(catalogmenu, value= ' Multi-AGC mode', uvalue='multiagc',/Checked_Menu)
  buttonagccatalog=widget_button(catalogmenu, value=' AGC Catalog (known cz)', uvalue='agcoverlay',/Checked_Menu)
  buttonagccatalog_nocz=widget_button(catalogmenu, value=' AGC Catalog (no known cz)', uvalue='agcoverlay_nocz',/Checked_Menu) 
  buttondhvvcatalog=widget_button(catalogmenu, value= ' HVC catalog (de Heij 2002)', uvalue='dhbboverlay',/Checked_Menu)
  buttonnvsscatalog=widget_button(catalogmenu, value=' NVSS Catalog ', uvalue='nvssoverlay',/Checked_Menu) 
  buttoncatalog3D=widget_button(catalogmenu, value=' 3D Catalog ', uvalue='catalog3D',/Checked_Menu) 
  buttoncatalog2D=widget_button(catalogmenu, value=' 2D Catalog (not active)', uvalue='catalog2D',/Checked_Menu)

;Special case of uvalue - an array that contains the currently marked values
  buttonchecklist=widget_button(catalogmenu, value=' Checklist (GALflux) ', uvalue='checklist', /Checked_Menu)

  buttoncntsrc=widget_button(catalogmenu, value=' Continuum (>250mJy/bm) ', uvalue='cntsrcs', /Checked_Menu)

state.helpmenu=widget_button(top_menu, value=' Help ')
  buttonquickstart=widget_button(state.helpmenu, value= ' Quickstart ', event_pro='gridview2_quickstart')
  buttonhelp=widget_button(state.helpmenu, value= ' About ', event_pro='gridview2_help')
  

;RIGHT BASE - upper base on the screen

rightbase=widget_base(tlb, xsize=1125, ysize=500, /row)
     state.plotwindowone=widget_draw(rightbase, xsize=495, ysize=495, frame=1, $
               uvalue='plotwindowone', event_pro='gridview2_display_event', /motion_events, /button_events, retain=2,$
               keyboard_events=1)
     state.plotwindowcolorbar=widget_draw(rightbase, xsize=50, ysize=495, frame=1, uvalue='plotwindowcolorbar', $
                                       event_pro='gridview2_display_event', /motion_events, /button_events)
     state.plotwindowtwo=widget_draw(rightbase, xsize=495, ysize=495, frame=1, $
               uvalue='plotwindowone', event_pro='gridview2_display_event', retain=2)



;LEFT BASE - Information display / controls

leftbase=widget_base(tlb, xsize=1000, ysize=375, /row)
 
infocontrolbase=widget_base(leftbase, xsize=400, ysize=370, /column)

     lowerbase=widget_base(infocontrolbase,  xsize=350, ysize=110,frame=1, /row)
     infobase=widget_base(lowerbase,  xsize=180, ysize=120, /column)
        ;label=widget_label(infobase, value='GRIDinfo', frame=1)
        label=widget_label(infobase, value=state.dataset, frame=1)
        label=widget_label(infobase, value=state.imagesize)
        label=widget_label(infobase, value=state.velrange)
        label=widget_label(infobase, value=state.epoch)
        ;label=widget_label(infobase, value=state.projection)
        
  mapchoicebase=widget_base(infobase, /row, /align_left, /exclusive)
          state.mapchoice[0]=widget_button(mapchoicebase, value='Spectral    ', uvalue='spectral')
          state.mapchoice[1]=widget_button(mapchoicebase, value='Continuum', uvalue='continuum')
           widget_control, state.mapchoice[0], set_button=1

     coordbase=widget_base(lowerbase,  xsize=210, ysize=120, /column)
        xbase=widget_base(coordbase, /row)
           label=widget_label(xbase, value='RA/DEC: ')
           state.xval=widget_label(xbase, value='12 23 24.5')
;        ybase=widget_base(coordbase, /row)
;           label=widget_label(ybase, value='DEC value: ')
           state.yval=widget_label(xbase, value='+12 23 24.5')
        zbase=widget_base(coordbase, /row)
           label=widget_label(zbase, value='cz value:  ')
           state.zval=widget_label(zbase, value=strcompress(grid.velarr[0]))
        xpixbase=widget_base(coordbase, /row)
           label=widget_label(xpixbase, value='(X,Y) pix:  ')
           state.xpix=widget_label(xpixbase, value='999')
;        ypixbase=widget_base(coordbase, /row)
;           label=widget_label(ypixbase, value='Y pixel:  ')
           state.ypix=widget_label(xpixbase, value='999')
        brightnessbase=widget_base(coordbase, /row)
           label=widget_label(brightnessbase, value='Intensity: ')
           state.intensitylevel=widget_label(brightnessbase, value='999.999999')
       

     sliderbase=widget_base(infocontrolbase, xsize=390, ysize=50, /row)
         buttonprevious=widget_button(sliderbase, xsize=40, ysize=15, $
                                 value='PREV', uvalue='previouschannel')
        buttonnext=widget_button(sliderbase, xsize=40, ysize=15, $
                                 value='NEXT', uvalue='nextchannel')

        state.velslider=widget_slider(sliderbase, min=0, max=n_elements(grid.velarr)-1, $
                                      title='Channel', value=0, xsize=230, uvalue='velslider', /drag)
        gotobase=widget_base(sliderbase, xsize=180, ysize=35, /row)
;        label=widget_label(gotobase, value='     Goto channel:')
        state.gotochannel=widget_text(gotobase, xsize=6, value='0',$ 
                                      /editable, uvalue='gotochannel')
     buttoncontrolbase=widget_base(infocontrolbase, xsize=390, ysize=40, /row, /frame)
        
        

        label=widget_label(buttoncontrolbase, value='Image smooth')
        state.smoothslider=widget_slider(buttoncontrolbase, min=0, max=10, $
                                        value=3, xsize=120, uvalue='smoothslider', /drag)

         integrationbase=widget_base(buttoncontrolbase, ysize=30, /row, /align_left, /frame)
        integrationbuttonbase=widget_base(integrationbase, ysize=25, /row, /align_left, /nonexclusive)
          state.velintegselect=widget_button(integrationbuttonbase, value='Boxcar (+/-)', $
                                             uvalue='velintegselect')
          state.velintegtext=widget_text(integrationbase, value='5', xsize=4, /editable, uvalue='velintegtext')
          label=widget_label(integrationbase, value='chan')



        choicebase=widget_base(infocontrolbase, /row, /align_left)
        

        polchoice=widget_droplist(choicebase, value=['Pol A', 'Pol B', 'Avg Pol'], $
                         event_pro='gridview2_event', uvalue='polchoice')

        widget_control, polchoice, SET_DROPLIST_SELECT=2

          ;state.polbutton[0]=widget_button(polbase, value='Pol A   ', uvalue='pola')
          ;state.polbutton[1]=widget_button(polbase, value='Pol B',    uvalue='polb')
          ;state.polbutton[2]=widget_button(polbase, value='Avg Pol', uvalue='avgpol')
          ;widget_control, state.polbutton[2], set_button=1
          state.currentpol=2   ;2 for averaged pols
          state.colortable=1

        ;choicebase2=widget_base(choicebase, /column, /align_left)
       

        agcinfobase=widget_base(choicebase, ysize=20, /row, /align_left, /frame)
        state.agcinfo=widget_label(agcinfobase, value='        AGC Information Displayed Here          ')
        
        ;specbuttonbase=widget_base(choicebase2, ysize=30, /row, /align_left, /frame)

       
        

        ;agcbase=widget_base(choicebase2, ysize=35, /row, /align_left)
        ;agcbuttonbase=widget_base(agcbase, ysize=25, /row, /align_left, /nonexclusive, /frame)
        ;state.agcoverlay=widget_button(agcbuttonbase, value='Overlay AGC galaxies', $
        ;                                     uvalue='agcoverlay')
        ;nvssbuttonbase=widget_base(agcbase, ysize=25, /row, /align_left, /nonexclusive, /frame)
        ;state.nvssoverlay=widget_buttON(nvssbuttonbase, value='Overlay NVSS sources', $
        ;                                     uvalue='nvssoverlay')
        
;SPECTRUM DISPLAY

spectrumbase=widget_base(leftbase, xsize=423, ysize=360, /column)

;label=widget_label(spectrumbase, value='Spectrum Display', frame=1, /align_center)
state.plotwindowspectrum=widget_draw(spectrumbase, xsize=420, ysize=210, frame=1, $
             uvalue='plotwindowspectrum', retain=2, /button_events, event_pro='gridview2_display_event')
;state.plotwindowspectrumweights=widget_draw(spectrumbase, xsize=420, ysize=100, frame=1, $
;             uvalue='plotwindowspectrumweights', retain=2, /button_events, event_pro='gridview2_display_event')

spectrumbuttonbase=widget_base(spectrumbase, xsize=420, ysize=25, /row)
 

  ;label=widget_label(spectrumbuttonbase, value='     Aperture FWHM: ')
  ;state.apertureFWHM=widget_text(spectrumbuttonbase, xsize=8, value='4.0', /editable, uvalue='apertureFWHM')
  ;label=widget_label(spectrumbuttonbase, value=' arcminutes')

;spectrumchannelbase=widget_base(spectrumbase, xsize=420, ysize=30, /row)
;label=widget_label(spectrumbuttonbase, value='Spectral window width: ')
;state.spectrumwindowwidth=widget_text(spectrumbuttonbase, xsize=8, value='250', /editable, uvalue='spectrumwindowwidth')
;label=widget_label(spectrumbuttonbase, value=' channels')

;spectrumaxisbase=widget_base(spectrumbase, xsize=420,
;ysize=30,/row,/align_left)

;buttonbase_upper = widget_base(spectrumbuttonbase, /Row, /Align_Left)
 getspectrum=widget_button(spectrumbuttonbase, xsize=100, ysize=25, value='PIXEL SPECTRUM', uvalue='getspectrum')
 cleanregion=widget_button(spectrumbuttonbase, xsize=100, ysize=25, value='CLEAN REGION', uvalue='cleanregion')
 fluxmeasure=widget_button(spectrumbuttonbase, xsize=100, ysize=25, value='MEASURE FLUX', uvalue='fluxmeasure')

 spectrumaxisbuttonbase=widget_base(spectrumbuttonbase, /row, /align_left, /exclusive)
  state.spectrumaxisbutton[0]=widget_button(spectrumaxisbuttonbase, value='Vel.', uvalue='axisvelocity')
  state.spectrumaxisbutton[1]=widget_button(spectrumaxisbuttonbase, value='Chan.', uvalue='axischannel')
  widget_control, state.spectrumaxisbutton[0], set_button=1

bothpolsbase=widget_base(spectrumbase, xsize=100, ysize=25, /row)
    button=widget_button(bothpolsbase, xsize=100, ysize=25, value='View POLS', uvalue='bothpols')
;buttonbase_lower = widget_base(spectrumbuttonbase, /Row, /Align_Left)
; polsbutton=widget_button(buttonbase_lower, xsize=120, ysize=25, value='View POLS', uvalue='bothpols')
;  spectrumaxisbuttonbase=widget_base(buttonbase_lower, /row, /align_left, /exclusive)
;   state.spectrumaxisbutton[0]=widget_button(spectrumaxisbuttonbase, value='Velocity', uvalue='axisvelocity')
;   state.spectrumaxisbutton[1]=widget_button(spectrumaxisbuttonbase, value='Channel', uvalue='axischannel')
;   widget_control, state.spectrumaxisbutton[0], set_button=1

;bothpolsbase=widget_base(spectrumbase, xsize=120, ysize=25, /row)
;    button=widget_button(buttonbase_lower, xsize=120, ysize=25, value='View POLS', uvalue='bothpols')


;DSS/NVSS Image Display

imagebase=widget_base(leftbase, xsize=230, ysize=360, /column)
 
  state.plotwindowimage=widget_draw(imagebase, xsize=210, ysize=210, frame=1, uvalue='plotwindowimage', /button_events)

  sizetextbase=widget_base(imagebase, xsize=300, ysize=30, /row)

  state.imagelabel=widget_droplist(sizetextbase, value=state.imageoptions, event_pro='gridview2_event', uvalue='imageoption')
  state.opticalsize=widget_text(sizetextbase, xsize=4, value='10', /editable, uvalue='opticalsize')
  label=widget_label(sizetextbase, value='arcmin')


;   info.ncalibcalmaskbeamselect=widget_droplist(ncalibcalmaskheader, value=beamselections, $
;           event_pro='bpdgui_calmaskselect')

;   state.imagelabel=widget_label(imagebase, value=' DSS Image Display', frame=1, /align_center)
  ;label=widget_label(imagebase, value='Right-click/press g for image')
  
  ;label=widget_label(sizetextbase, value='Image size:  ')
  
  ;spacer=widget_base(imagebase, xsize=200, ysize=20)
  ;state.copyright=widget_draw(imagebase, xsize=150, ysize=20, uvalue='copyright', $
;                         /align_right, /button_events, event_pro='gridview2_display_event')
   
;Realization
widget_control, tlb, /realize

state.baseID=tlb

;Xmanager startup
xmanager, 'gridview2', state.baseID, /no_block

 ;widget_control, state.copyright, get_value=index
 ;wset, index
 ;xyouts, 10, 6, 'LOVEDATA, Inc.  Ithaca, NY', /device

endif 

end 


;------------------------------------------------------------------
;DISPLAY EVENT HANDLER

pro gridview2_display_event, event

common gridview2_state
common gridstate

widget_control, event.id, get_uvalue=uvalue

case uvalue of

    'plotwindowone': begin
        
        



;--------------------------------------------------------------------------
;--------------------------------------------------------------------------

       ;First section of case choice is reserved for zoom box
       ;Begin preliminary case statement for mouse control and motion 

if (event.type le 2 AND event.press ne 4 AND state.getspectrumstatus eq 0) then begin

;print, 'EVENT TYPE ', event.type
;print, 'EVENT PRESS ' , event.press

eventTypes = ['DOWN', 'UP', 'MOTION']
thisEvent = eventTypes[event.type]
state.xsize=495   ;Hardwired
state.ysize=495   ;Hardwired

widget_control, state.gotochannel, get_value=channelstring
               channel=long(channelstring[0])

CASE thisEvent OF

   'DOWN': BEGIN
      if (event.press eq 1) then begin
       ; Turn motion events on for the draw widget.

       state.mousestatus=1

      Widget_Control, state.plotwindowone, Draw_Motion_Events=1
      widget_control, state.plotwindowone, get_value=index
      state.wid=index
      state.drawID=state.plotwindowone

         ; Create a pixmap. Store its ID. Copy window contents into it.

      Window, 19, /Pixmap, XSize=state.xsize, YSize=state.ysize
      state.pixID = !D.Window
      Device, Copy=[0, 0, state.xsize, state.ysize, 0, 0, state.wid]

         ; Get and store the static corner of the box.

      state.zoomsx = event.x
      state.zoomsy = event.y

      xpos=round(event.x-state.px[0]-1) 
      ypos=round(event.y-state.py[0]-1)

     

      endif
  END


   'UP': BEGIN
      
      if (state.mousestatus eq 1) then begin

       state.mousestatus=0

         ; Erase the last box drawn. Destroy the pixmap.

      widget_control, state.plotwindowone, get_value=index
      state.wid=index
      state.drawID=state.plotwindowone

      WSet, state.wid
      Device, Copy=[0, 0, state.xsize, state.ysize, 0, 0, state.pixID]
      
      ; Order the box coordinates.

      sx = Min([state.zoomsx, event.x], Max=dx)
      sy = Min([state.zoomsy, event.y], Max=dy)

        ;New min and max for zooming
   
      ;print, sx,sy,dx,dy

     WDelete, state.pixID

     ;Determine coordinates FOR SQUARE BOX to zoom in on

           ;xpos=round(event.x-state.px[0]-1) 
           ;ypos=round(event.y-state.py[0]-1)
           
          
           ;xcube=-round((float(state.xcubemax-state.xcubemin+1)/state.sx)*xpos)+(state.xcubemax)
           ;ycube=round((float(state.ycubemax-state.ycubemin+1)/state.sy)*ypos)+(state.ycubemin)

     width=abs(dx-sx)

     xpos_min=sx
     xpos_max=dx
     ypos_min=sy

     

    
          
    
;If status codes are correct, go ahead and zoom in
if (state.zoomstatus eq 1 AND state.fluxmeasurestatus eq 0) then begin

    
     if (sy lt dy) then ypos_max=sy+width
     if (sy gt dy) then ypos_max=sy-width
     if (sy eq dy) then ypos_max=dy
     ;Convert from device window coords to
     ;   enlarged cube coords, and then to original cube cooridnates

     ;print, 'Device ',xpos_min, xpos_max, ypos_min, ypos_max
     ;print, 'S and D ', sx,dx,sy,dy

     xpos_min=round(xpos_min-state.px[0]-1) 
     xpos_max=round(xpos_max-state.px[0]-1) 
     ypos_min=round(ypos_min-state.py[0]-1)
     ypos_max=round(ypos_max-state.py[0]-1)




     if ((sx eq dx AND sy eq dy) OR $
         (xpos_min lt 0  OR ypos_min lt 0 OR $
          xpos_max gt state.sx-1 OR ypos_max gt state.sy-1) $
        ) then begin
         ;print, 'Zooming out'
         state.xcubemax=n_elements(grid.d[0,0,*,0])-1
         state.xcubemin=0
         state.ycubemin=0
         state.ycubemax=n_elements(grid.d[0,0,0,*])-1
         
         state.ramin=grid.ramin
         state.ramax=grid.ramin+grid.(6)*(grid.deltara/3600.0)
         state.decmin=grid.decmin
         state.decmax=grid.decmin+grid.(7)*(grid.deltadec/60.0) 

     endif else begin
     ;print, 'Cube LARGE ',xpos_min, xpos_max, ypos_min, ypos_max

     newxmax=-round((float(state.xcubemax-state.xcubemin+1)/state.sx)*xpos_min)+(state.xcubemax)
     newxmin=-round((float(state.xcubemax-state.xcubemin+1)/state.sx)*xpos_max)+(state.xcubemax)
     newymin= round((float(state.ycubemax-state.ycubemin+1)/state.sy)*ypos_min)+(state.ycubemin)
     newymax= round((float(state.ycubemax-state.ycubemin+1)/state.sy)*ypos_max)+(state.ycubemin)

     state.xcubemax=newxmax
     state.xcubemin=newxmin
     state.ycubemin=newymin
     state.ycubemax=newymax


     ;Protect against going 'outside' the cube
     if (state.xcubemin lt 0) then state.xcubemin=0
     if (state.ycubemin lt 0) then state.ycubemin=0
     if (state.xcubemax gt n_elements(grid.d[0,0,*,0])-1) then state.xcubemax=n_elements(grid.d[0,0,*,0])-1
     if (state.ycubemax gt n_elements(grid.d[0,0,0,*])-1) then state.ycubemax=n_elements(grid.d[0,0,0,*])-1

     ;print, 'Cube Small ',state.xcubemin,state.xcubemax,state.ycubemin,state.ycubemax

     index_min=where(grid.grid_makeup.i eq state.xcubemin AND $
                     grid.grid_makeup.j eq state.ycubemin)
     index_max=where(grid.grid_makeup.i eq state.xcubemax AND $
                     grid.grid_makeup.j eq state.ycubemax)

     state.ramin=grid.grid_makeup[index_min].ra-(grid.deltara/3600.0)/2.0
     state.ramax=grid.grid_makeup[index_max].ra+(grid.deltara/3600.0)/2.0
     state.decmin=grid.grid_makeup[index_min].dec-(grid.deltadec/3600.0)/2.0
     state.decmax=grid.grid_makeup[index_max].dec+(grid.deltadec/3600.0)/2.0

     ;print, state.ramin, state.ramax

     if (state.ramin gt state.ramax) then state.ramin=state.ramin-24.0

     ;print, state.ramin, state.ramax

     ;state.ramin=grid.ramin+state.xcubemin*(grid.deltara/3600.0)
     ;state.ramax=grid.ramin+state.xcubemax*(grid.deltara/3600.0)
     ;state.decmin=grid.decmin+state.ycubemin*(grid.deltadec/60.0)
     ;state.decmax=grid.decmin+state.ycubemax*(grid.deltadec/60.0)
     
     endelse

;     print, xpos_min,xpos_max, ypos_min, ypos_max


  ; Turn draw motion events off. Clear any events queued for widget.

    Widget_Control, state.drawID, Clear_Events=1

     gridview2_display, channel
     gridview2_display_weights, channel


endif


;If fluxmeasure is turned on, don't zoom - open the flux measuring
;                                          program and set the correct
;                                          windowingn parameters
if (state.zoomstatus eq 0 AND state.fluxmeasurestatus eq 1) then begin

     ypos_max=dy

     xpos_min=round(xpos_min-state.px[0]-1) 
     xpos_max=round(xpos_max-state.px[0]-1) 
     ypos_min=round(ypos_min-state.py[0]-1)
     ypos_max=round(ypos_max-state.py[0]-1)



     xmax=-round(float(state.xcubemax-state.xcubemin+1)/state.sx*xpos_min)+(state.xcubemax)
     xmin=-round(float(state.xcubemax-state.xcubemin+1)/state.sx*xpos_max)+(state.xcubemax)
     ymin=round(float(state.ycubemax-state.ycubemin+1)/state.sy*ypos_min)+(state.ycubemin)
     ymax=round(float(state.ycubemax-state.ycubemin+1)/state.sy*ypos_max)+(state.ycubemin)

     ;Protect against going 'outside' the cube
     if (xmin lt 0) then xmin=0
     if (ymin lt 0) then ymin=0
     if (xmax gt n_elements(grid.d[0,0,*,0])-1) then xmax=n_elements(grid.d[0,0,*,0])-1
     if (ymax gt n_elements(grid.d[0,0,0,*])-1) then ymax=n_elements(grid.d[0,0,0,*])-1

     ;status=dialog_message('COMING SOON!  Measure flux: ('+strcompress(xmin)+','+strcompress(ymin)+') to ('+$
     ;                                                      strcompress(xmax)+','+strcompress(ymax)+')')
     
     ;Call the flux measureing routine GALFLUX

     state.fluxmeasurestatus=0
     state.zoomstatus=1

     Widget_Control, state.drawID, Clear_Events=1

     galflux,llx=xmin, lly=ymin, urx=xmax, ury=ymax

 ;    !p.clip=[90,60,573,370,0,1000]

 ;    !x.region=[ 0.00000 ,     1.00000]
 ;    !x.s=[0.15000500,   0.00078690126]
 ;    !x.margin=[ 10.0000 ,     3.00000]
 ;    !x.window=[ 0.150005 ,    0.955005]
 ;    !x.crange=[0, n_elements(grid.velarr)-1]
;     !x.range=[0, n_elements(grid.velarr)-1]


;     !y.region=[ 0.00000 ,     1.00000]
;     !y.s=[0.40833835,     0.025833334]   ;set for vertical scaling between -10 and 20
;     !y.margin=[ 4.0,     2.00000]
;     !y.window=[ 0.150005 ,    0.925005]

 endif

if (state.zoomstatus eq 0 AND state.fluxmeasurestatus eq 2) then begin
     ypos_max=dy

     xpos_min=round(xpos_min-state.px[0]-1)
     xpos_max=round(xpos_max-state.px[0]-1)
     ypos_min=round(ypos_min-state.py[0]-1)
     ypos_max=round(ypos_max-state.py[0]-1)

     xmax=-round(float(state.xcubemax-state.xcubemin+1)/state.sx*xpos_min)+(state.xcubemax)
     xmin=-round(float(state.xcubemax-state.xcubemin+1)/state.sx*xpos_max)+(state.xcubemax)
     ymin=round(float(state.ycubemax-state.ycubemin+1)/state.sy*ypos_min)+(state.ycubemin)
     ymax=round(float(state.ycubemax-state.ycubemin+1)/state.sy*ypos_max)+(state.ycubemin)

     ;Protect against going 'outside' the cube
     if (xmin lt 0) then xmin=0
     if (ymin lt 0) then ymin=0
     if (xmax gt n_elements(grid.d[0,0,*,0])-1) then xmax=n_elements(grid.d[0,0,*,0])-1
     if (ymax gt n_elements(grid.d[0,0,0,*])-1) then ymax=n_elements(grid.d[0,0,0,*])-1

     ;Call the flux measureing routine GALFLUX
     state.fluxmeasurestatus=0
     state.zoomstatus=1

     Widget_Control, state.drawID, Clear_Events=1

     ;galflux,llx=xmin, lly=ymin, urx=xmax, ury=ymax
     cleanview, llx=xmin, lly=ymin, urx=xmax, ury=ymax

endif
	

;If zoom status is on and the user wished to create a box overlay,
;then save the box they have just outlined.

if (state.zoomstatus eq 0 AND state.overlayboxstatus eq 1) then begin

;Stored as (status, sx,sy,dx,dy) in the intarr state.boxoverlay

    boxindex=where(state.overlaybox[*,0] eq 0)

    result1=convert_coord(sx, sy, /device, /double, /to_data)
    result2=convert_coord(dx,dy, /device, /double, /to_data)

    if (boxindex[0] ne -1) then state.overlaybox[boxindex[0],*]=[1,result1[0], result1[1], result2[0], result2[1]]

    ;Reset control parameters
    state.overlayboxstatus=0
    state.zoomstatus=1
    widget_control, state.overlayboxbutton, set_Button=0

   

    Widget_Control, state.drawID, Clear_Events=1

     gridview2_display, channel
     gridview2_display_weights, channel


endif



  
    endif





      END

   'MOTION': BEGIN
       if (event.press eq 0 AND state.mousestatus eq 1) then begin

         ; Here is where the actual box is drawn and erased.
         ; First, erase the last box.

      widget_control, state.plotwindowone, get_value=index
      state.wid=index
      state.drawID=state.plotwindowone

      WSet, state.wid
      Device, Copy=[0, 0, state.xsize, state.ysize, 0, 0, state.pixID]

         ; Get the coodinates of the new box and draw it.

      sx = state.zoomsx
      sy = state.zoomsy
      dx = event.x
      dy = event.y
      loadct, 1, /silent
      state.boxcolor=!D.N_Colors-1

     
      width=abs(dx-sx)
     

      ;xcenter=sx+(dx-sx)/2.0
      ;ycenter=sy+(dy-sy)/2.0

      ;tvboxbk, width, xcenter,ycenter, $
      ;               linestyle=0, /device, color='0000FF'XL, thick=2.0    ;RED

     ;PlotS, [sx, sx, dx, dx, sx], [sy, dy, dy, sy, sy], /Device, $
     ;    Color=state.boxColor, thick=1.5


      device, decomposed=1


      ;plot a flux measurement box of arbitrary size
      if (state.fluxmeasurestatus eq 1 AND state.zoomstatus eq 0) then begin
           color='0000FF'XL  ;RED
      
           PlotS, [sx, sx, dx, dx, sx], [sy, dy, dy, sy, sy], /Device, $
          Color=color, thick=1.5

      endif

      ;plot a CLEAN region box of arbitrary size
      if (state.fluxmeasurestatus eq 2 AND state.zoomstatus eq 0) then begin
	   color='FF00FF'XL  ;RED

	   PlotS, [sx, sx, dx, dx, sx], [sy, dy, dy, sy, sy], /Device, $
          Color=color, thick=1.5

      endif


      ;plot a zoom box (GREEN) of square proportions
      if (state.fluxmeasurestatus eq 0 AND state.zoomstatus eq 1) then begin
           color='00FF00'XL  ;GREEN


         if (sy lt dy) then begin
          PlotS, [sx, sx, dx, dx, sx], [sy, sy+width, sy+width, sy, sy], /Device, $
            color=color, thick=1.5 
            endif else begin 
          PlotS, [sx, sx, dx, dx, sx], [sy, sy-width, sy-width, sy, sy], /Device, $
            color=color, thick=1.5    ;RED
        endelse

      endif

      if (state.fluxmeasurestatus eq 1 AND state.zoomstatus eq 0) then begin
          xyouts, sx+10, sy-10, 'FLUX BOX', /device, color=color
      endif
      
      if (state.fluxmeasurestatus eq 2 AND state.zoomstatus eq 0) then begin
	  xyouts, sx+10, sy-10, 'CLEAN Region', /device, color=color
      endif


      ;Plot a fiducial box of arbitrary size that will stay on screen
      if (state.overlayboxstatus eq 1) then begin
          plotS, [sx, sx, dx, dx, sx], [sy, dy, dy, sy, sy], /Device, $
          Color='FFFFFF'XL, thick=1.5 ;WHITE

      endif



      device, decomposed=0

      loadct, state.colortable, /silent

      endif

  END

ELSE:

ENDCASE

endif



        
;--------------------------------------------------------------------------
;--------------------------------------------------------------------------

       widget_control, state.plotwindowone, get_value=index
       wset, index

       xdevice=event.x
       ydevice=event.y

       ramin=state.ramin
       ramax=state.ramax
       decmin=state.decmin
       decmax=state.decmax

       widget_control, state.gotochannel, get_value=channelstring
         channel=long(channelstring[0]) 

       hor, ramax, ramin
       ver, decmin, decmax

       result=convert_coord(xdevice, ydevice, /device, /double, /to_data)

       xdata=result[0]
       ydata=result[1]

       if (xdata lt ramin OR xdata gt ramax OR ydata lt decmin OR ydata gt decmax) then begin
           widget_control, state.xval, set_value='----'
           widget_control, state.yval, set_value='----'
           widget_control, state.xpix, set_value='----'
           widget_control, state.ypix, set_value='----'
           widget_control, state.intensitylevel, set_value='----'

        endif else begin
           radecconvert, xdata, ydata, rastring, decstring
           widget_control, state.xval, set_value=string(rastring[0], format='(a11)')+' ,'
           widget_control, state.yval, set_value=' '+string(decstring[0], format='(a11)')

                                ;Update window for optical image
                                ;display if there is a mouse click

           xpos=round(event.x-state.px[0]-1) 
           ypos=round(event.y-state.py[0]-1)
           
          
           xcube=-round((float(state.xcubemax-state.xcubemin+1)/state.sx)*xpos)+(state.xcubemax)
           ycube=round((float(state.ycubemax-state.ycubemin+1)/state.sy)*ypos)+(state.ycubemin)

           if (xcube lt 0) then xcube=0
           if (ycube lt 0) then ycube=0
           if (xcube gt n_elements(grid.d[0,0,*,0])-1) then xcube=n_elements(grid.d[0,0,*,0])-1
           if (ycube gt n_elements(grid.d[0,0,0,*])-1) then ycube=n_elements(grid.d[0,0,0,*])-1

            widget_control, state.xpix, set_value=strcompress(xcube, /remove_all)
            widget_control, state.ypix, set_value=strcompress(ycube, /remove_all)

         if (state.mapcurrent eq 'spectral') then begin

            if (state.currentpol eq 0 OR state.currentpol eq 1) then $
                 intensityvalue=grid.d[channel,state.currentpol,xcube,ycube]
            if (state.currentpol eq 2) then $
                 intensityvalue=(grid.d[channel,0,xcube,ycube]+grid.d[channel,1,xcube,ycube])/2.0

            widget_control, state.intensitylevel, set_value=strcompress(intensityvalue, /remove_all)

         endif

          if (state.mapcurrent eq 'continuum') then begin

             if (state.currentpol eq 0 OR state.currentpol eq 1) then $
                 intensityvalue=grid.cont[state.currentpol,xcube,ycube]
             if (state.currentpol eq 2) then $
                 intensityvalue=(grid.cont[0,xcube,ycube]+grid.cont[1,xcube,ycube])/2.0

             widget_control, state.intensitylevel, set_value=strcompress(intensityvalue, /remove_all)

         endif
        
         ;If AGC catalog is turned on, check mouse coordinates
         
         agcindex=where(state.currentagc[1,*] lt event.x+5 AND $
                        state.currentagc[1,*] gt event.x-5 AND $
                        state.currentagc[2,*] lt event.y+5 AND $
                        state.currentagc[2,*] gt event.y-5)

         if (agcindex[0] ne -1) then begin
               agcnum=state.currentagc[0,agcindex]
               agcnum=agcnum[0]
               catindex=where(state.agc.agcnumber eq agcnum)
               stringout='A'+strcompress(state.agc.agcnumber[catindex], /remove_all)+',Type='+ $
                         strcompress(state.agc.DESCRIPTION[catindex], /remove_all)+',Vopt='+ $
                         strcompress(state.agc.Vopt[catindex], /remove_all)+',V21='+ $
                         strcompress(state.agc.V21[catindex], /remove_all)+',detcode='+ $
                         strcompress(state.agc.detcode[catindex], /remove_all)
               widget_control, state.agcinfo, set_value=stringout[0]
           endif else begin 
               widget_control, state.agcinfo, set_value='--------------' 
           endelse 
             

           ;IF N key is hit, query NED and bring up text box

          if (event.type eq 5 AND (StrTrim(event.ch,2) eq 'n' OR StrTrim(event.ch,2) eq 'N') $
                 AND state.keystatus eq 0) then begin

              gridview2_nedgui, xdata,ydata

              state.keystatus=1
               

          endif

          ;IF P key is hit, open the postscript dialog box

           if (event.type eq 5 AND (StrTrim(event.ch,2) eq 'p' OR StrTrim(event.ch,2) eq 'P') $
                 AND state.keystatus eq 0) then begin

               xcube=-round((float(state.xcubemax-state.xcubemin+1)/state.sx)*xpos)+(state.xcubemax)
               ycube=round((float(state.ycubemax-state.ycubemin+1)/state.sy)*ypos)+(state.ycubemin)

                

                  if (xcube lt 0) then xcube=0
                  if (ycube lt 0) then ycube=0
                  if (xcube gt n_elements(grid.d[0,0,*,0])-1) then xcube=n_elements(grid.d[0,0,*,0])-1
                  if (ycube gt n_elements(grid.d[0,0,0,*])-1) then ycube=n_elements(grid.d[0,0,0,*])-1

              gridview2_ps_output_gui, event, xcube, ycube

              state.keystatus=0
               

          endif

          ;If H key is hit, open the drift viewing window

           if (event.type eq 5 AND (StrTrim(event.ch,2) eq 'h' OR StrTrim(event.ch,2) eq 'H') $
                 AND state.keystatus eq 0) then begin

               xcube=-round((float(state.xcubemax-state.xcubemin+1)/state.sx)*xpos)+(state.xcubemax)
               ycube=round((float(state.ycubemax-state.ycubemin+1)/state.sy)*ypos)+(state.ycubemin)

                  if (xcube lt 0) then xcube=0
                  if (ycube lt 0) then ycube=0
                  if (xcube gt n_elements(grid.d[0,0,*,0])-1) then xcube=n_elements(grid.d[0,0,*,0])-1
                  if (ycube gt n_elements(grid.d[0,0,0,*])-1) then ycube=n_elements(grid.d[0,0,0,*])-1

                

              gridview2_fetchdrifts, xcube, ycube

              state.keystatus=1
               

          endif



          ;If G key is hit, query DSS 2 Blue and Sloan

          if (event.type eq 5 AND (StrTrim(event.ch,2) eq 'g' OR StrTrim(event.ch,2) eq 'G') $
                 AND state.keystatus eq 0) then begin

               case state.mapcurrent of 

               'spectral': begin
               
               if (state.getspectrumstatus eq 0) then begin
               widget_control, state.baseID, hourglass=1

               
        
               widget_control, state.opticalsize, get_value=imagesizestring
               imagesize=float(imagesizestring[0])
               
               xsizepixel=state.sx
               ysizepixel=state.sy

               xsizearcmin=15.0*(ramax-ramin)*60.0
               ysizearcmin=(decmax-decmin)*60.0

               pixscalex=xsizepixel/xsizearcmin
               pixscaley=ysizepixel/ysizearcmin

               device, decomposed=1
               tvboxbk, imagesize*pixscalex, event.x,event.y, $
                     linestyle=0, /device, color='0000FF'XL ;RED
              
               state.currentimage_rahr=xdata
               state.currentimage_decdeg=ydata               

               print, xdata
               if (xdata lt 0.0) then xdata=24.0+xdata
               print, xdata
                         
               queryDSS, [xdata*15.0,ydata], image, header, imsize=imagesize, survey='2b'
               state.dssimage=congrid(image, 210,210)

               osfamily = strupcase(!version.os_family)
               if (osfamily eq 'UNIX') then begin

                 url='http://casjobs.sdss.org/ImgCutoutDR6/getjpeg.aspx?ra='+$
                    strcompress(xdata*15.0, /remove_all)+$
                    '&dec='+strcompress(ydata, /remove_all)+$
                    '&scale='+strcompress(imagesize/6.67,/remove_all)+$
                    '&opt=GI&width=400&height=400'

                 filename='~/12junksdss.jpg'

                 spawn, 'wget -q -O '+ filename + " '" + url + "'"

                 read_jpeg, filename, image, true=1

                 state.sloanimage=congrid(image, 3,210,210)

                 spawn, '/bin/rm -r ~/12junksdss.jpg'
               endif


               gridview2_cutout_display

               widget_control, state.baseID, hourglass=0
        

       
               endif

           end

           'continuum': begin

               
               widget_control, state.baseID, hourglass=1

               widget_control, state.opticalsize, get_value=imagesizestring
               imagesize=float(imagesizestring[0])

               xsizepixel=state.sx
               ysizepixel=state.sy

               xsizearcmin=15.0*(ramax-ramin)*60.0
               ysizearcmin=(decmax-decmin)*60.0

               pixscalex=xsizepixel/xsizearcmin
               pixscaley=ysizepixel/ysizearcmin

               device, decomposed=1
               tvboxbk, imagesize*pixscalex, event.x,event.y, $
                     linestyle=0, /device, color='0000FF'XL ;RED

               osfamily = strupcase(!version.os_family)
               if (osfamily eq 'UNIX') then begin
 
                url='http://skys.gsfc.nasa.gov/cgi-bin/pskcall?VCOORD='+strcompress(xdata*15.0, /remove_all)+$
                ','+strcompress(ydata, /remove_all)+'&SURVEY=VLA+NVSS+(1.4+Ghz)&SFACTR='+$
                   strcompress(imagesize/60.0, /remove_all)+$
                   '&PIXELX=400&PIXELY=400&COLTAB=Rainbow&RETURN=GIF'

                spawn, 'wget -q -O ~/temp1357_IDL.gif ' + "'" + url + "'"
                read_gif, '~/temp1357_IDL.gif', image, r,g,b

                spawn, '/bin/rm -r ~/temp1357_IDL.gif'
            
                widget_control, state.plotwindowimage, get_value=index
                wset, index
                device, decomposed=0
                loadct, 13, /silent
                tvscl, congrid(image, 210,210)
               endif else begin
                xyouts, 100,100, 'Available only on Linux', alignment=0.5, /device
               endelse

               widget_control, state.baseID, hourglass=0

           end 

           else:

       endcase


          endif

          


          

       endelse

       ;Update window for optical/nvss image display

       
   
       ;if (state.getspectrumstatus eq 1 and state.mapcurrent eq 'spectral') then begin

       ;    if (xdata gt ramin AND xdata lt ramax AND ydata gt decmin AND ydata lt decmax) then begin

       ;        xsizepixel=state.sx
       ;        ysizepixel=state.sy

       ;        xsizearcmin=15.0*(ramax-ramin)*60.0
       ;        ysizearcmin=(decmax-decmin)*60.0

       ;        pixscalex=xsizepixel/xsizearcmin
       ;        pixscaley=ysizepixel/ysizearcmin

       ;        widget_control, state.apertureFWHM, get_value=aperturestring

        ;       radiuspixel=double(aperturestring[0])*pixscalex/2.0

         ;      gridview2_display, channel

        ;       red='0000FF'XL

     ;          device, decomposed=1
     ;          tvboxbk, 31, event.x,event.y, /device, color=red
     ;          tvcircle, radiuspixel, event.x, event.y, /device, color=red
     ;          device, decomposed=0

     ;      endif else begin
     ;          gridview2_display, channel
     ;          gridview2_display_weights, channel
     ;      endelse


     ; endif


      


   end


   'plotwindowcolorbar': begin

       case event.type of
           0: begin    ;Mouse button press

               state.mousestatus = 1 ;Mouse button is down
               state.colorbar_y=event.y

             end

           1: begin    ;Mouse button release

                state.mousestatus = 0   ;Mouse button is up
                state.colorbar_y=0   ;Reset value
             end

           2: begin           ; Button press
               
               if (state.mousestatus eq 1) then begin

                   y_diff_percent=abs(event.y-state.colorbar_y)/float(state.sy)

               if (event.y lt state.colorbar_y) then $
                 state.stretchupper=state.stretchupper-(y_diff_percent*255)
               if (event.y gt state.colorbar_y) then $
                 state.stretchupper=state.stretchupper+(y_diff_percent*255)

               widget_control, state.gotochannel, get_value=channelstring
               channel=long(channelstring[0])

               gridview2_refresh, channel
               
                 state.colorbar_y=event.y

                 endif

             end

        else:

         endcase


    end

'copyright': begin
 
               if (event.press gt 0) then gridview2_utility
 
             end


'plotwindowspectrum': begin

;For this window deal only with up, down, and mouse movement for box control

if (event.type le 2) then begin

eventTypes = ['DOWN', 'UP', 'MOTION']
thisEvent = eventTypes[event.type]
state.xsize=420
state.ysize=210

widget_control, state.gotochannel, get_value=channelstring
               channel=long(channelstring[0])

CASE thisEvent OF

   'DOWN': BEGIN
       ; Turn motion events on for the draw widget.

      Widget_Control, state.plotwindowspectrum, Draw_Motion_Events=1
      widget_control, state.plotwindowspectrum, get_value=index
      state.wid=index
      state.drawID=state.plotwindowspectrum

         ; Create a pixmap. Store its ID. Copy window contents into it.

      Window, 19, /Pixmap, XSize=state.xsize, YSize=state.ysize
      state.pixID = !D.Window
      Device, Copy=[0, 0, state.xsize, state.ysize, 0, 0, state.wid]

         ; Get and store the static corner of the box.

      state.zoomsx = event.x
      state.zoomsy = event.y

      ;print, event.y

      gridview2_plotspectrum


     END


   'UP': BEGIN
       
         ; Erase the last box drawn. Destroy the pixmap.

      widget_control, state.plotwindowspectrum, get_value=index
      state.wid=index
      state.drawID=state.plotwindowspectrum

      WSet, state.wid
      Device, Copy=[0, 0, state.xsize, state.ysize, 0, 0, state.pixID]
     
         ; Order the box coordinates.

      sx = Min([state.zoomsx, event.x], Max=dx)
      sy = Min([state.zoomsy, event.y], Max=dy)

        ;New min and max for zooming
   
      ;print, sx,sy,dx,dy
   
;NOTE:  Use this scaling because there are TWO plots in the window
sy=1.44056*sy-90.7
dy=1.44056*dy-90.7
     
spectrum=state.spectrum
spectrumweights=state.spectrumweights

;Plot velocity on x-axis
if (state.spectrumaxisstatus eq 0) then begin

;Redundant, but just for consistency
widget_control, state.plotwindowspectrum, get_value=index
wset, index

hor, grid.velarr[state.spectrum_xmax], grid.velarr[state.spectrum_xmin]
ver, state.spectrum_ymin, state.spectrum_ymax

      resultcoords_min=convert_coord(sx, sy, /device, /double, /to_data)
      resultcoords_max=convert_coord(dx, dy, /device, /double, /to_data)

      startvel=round(resultcoords_min[0])
      stopvel=round(resultcoords_max[0])
      minflux=round(resultcoords_min[1])
      maxflux=round(resultcoords_max[1])

;print, startvel, stopvel


chanmin_array=where(grid.velarr ge stopvel)
startchannel=chanmin_array[n_elements(chanmin_array)-1]
chanmax_array=where(grid.velarr le startvel)
stopchannel=chanmax_array[0]

;Added BK March 13, 2006
if (stopchannel gt n_elements(grid.velarr)-1) then stopchannel=n_elements(grid.velarr)-1
if (startchannel lt 0) then starthchannel=0

;Situation where user clicks to zoom out
if (sx eq dx OR sy eq dy) then begin
startchannel=0
stopchannel=n_elements(grid.velarr)-1
minflux=min(spectrum)
maxflux=max(spectrum)
endif

if (stopchannel eq -1) then stopchannel=n_elements(grid.velarr)-1
if (startchannel eq -1) then startchannel=0

;Repeat for spectral weights

;widget_control, state.plotwindowspectrumweights, get_value=index
;wset, index

device, decomposed=1

hor, grid.velarr[stopchannel], grid.velarr[startchannel]
ver, 0.0,1.1

color='0000FF'XL   ;RED
plot,grid.velarr[startchannel:stopchannel],spectrumweights[startchannel:stopchannel]/state.maxweight,$
           xstyle=1, ystyle=1, xtitle='Velocity [km/s]', position=[0.07,0.1,0.98,0.3], $
           charsize=0.7, ytitle='Weight', /nodata

oplot,grid.velarr[startchannel:stopchannel], spectrumweights[startchannel:stopchannel]/state.maxweight, color=color

flag, grid.velarr[channel], color='00FF00'XL  ; GREEN current channel (velocity)


hor, grid.velarr[stopchannel], grid.velarr[startchannel]
ver, minflux,maxflux


color='0000FF'XL   ;RED
plot,grid.velarr[startchannel:stopchannel],spectrum[startchannel:stopchannel],$
      xstyle=1, ystyle=1, position=[0.07,0.3, 0.98,0.98], /noerase, $
      charsize=0.7, ytitle='Flux Density [mJy/beam]', /nodata, xtickn=[' ',' ',' ' ,' ' ,' ' ,' ',' ']

oplot,grid.velarr[startchannel:stopchannel], spectrum[startchannel:stopchannel], color=color


flag, grid.velarr[channel], color='00FF00'XL  ; GREEN current channel (velocity)



device, decomposed=0

endif else begin   ;Plot channel on x-axis

widget_control, state.plotwindowspectrum, get_value=index
wset, index

hor, state.spectrum_xmin, state.spectrum_xmax
ver, state.spectrum_ymin, state.spectrum_ymax

      resultcoords_min=convert_coord(sx, sy, /device, /double, /to_data)
      resultcoords_max=convert_coord(dx, dy, /device, /double, /to_data)

      startchannel=round(resultcoords_min[0])
      stopchannel=round(resultcoords_max[0])
      minflux=round(resultcoords_min[1])
      maxflux=round(resultcoords_max[1])

;Repeat for spectral weights


if (startchannel lt 0) then startchannel=0
if (stopchannel eq -1) then stopchannel=n_elements(grid.velarr)-1

hor, startchannel, stopchannel
ver, 0.0,1.1

device, decomposed=1

color='0000FF'XL   ;RED
plot,spectrumweights/state.maxweight,xstyle=1, ystyle=1, position=[0.07,0.1,0.98,0.3], $
    xtitle='Channel', charsize=0.7, ytitle='Weight', /nodata
oplot, spectrumweights/state.maxweight, color=color
flag, channel, color='00FF00'XL  ; GREEN centroid for redshift



hor, startchannel, stopchannel
ver, minflux,maxflux


color='0000FF'XL   ;RED
plot,spectrum,xstyle=1, ystyle=1,position=[0.07,0.3, 0.98,0.98], $
     charsize=0.7, ytitle='Flux Density [mJy/beam]', /nodata, /noerase, xtickn=[' ',' ',' ',' ',' ',' ',' ']
oplot, spectrum, color=color

flag, channel, color='00FF00'XL  ; GREEN centroid for redshift



device, decomposed=0

endelse

;Store away for any changes
state.spectrum=spectrum
state.spectrumweights=spectrumweights
state.spectrum_xmin=startchannel
state.spectrum_xmax=stopchannel
state.spectrum_ymin=minflux
state.spectrum_ymax=maxflux

WDelete, state.pixID

  ; Turn draw motion events off. Clear any events queued for widget.

    Widget_Control, state.drawID, Draw_Motion_Events=0, Clear_Events=1

   gridview2_display, channel
   gridview2_display_weights, channel

    END

   'MOTION': BEGIN

         ; Here is where the actual box is drawn and erased.
         ; First, erase the last box.

      widget_control, state.plotwindowspectrum, get_value=index
      state.wid=index
      state.drawID=state.plotwindowspectrum

      WSet, state.wid
      Device, Copy=[0, 0, state.xsize, state.ysize, 0, 0, state.pixID]

         ; Get the coodinates of the new box and draw it.

      sx = state.zoomsx
      sy = state.zoomsy
      dx = event.x
      dy = event.y
      loadct, 1, /silent
      state.boxcolor=!D.N_Colors-1
      PlotS, [sx, sx, dx, dx, sx], [sy, dy, dy, sy, sy], /Device, $
         Color=state.boxColor
      loadct, state.colortable, /silent

  END

  ELSE:

ENDCASE

endif

end

else:

endcase


case event.type of

1: begin
;BK Note 1
   if (state.getspectrumstatus eq 1) then $
         gridview2_fetch_spectrum, $
           round(event.x-state.px[0]-1), round(event.y-state.py[0]-1)
           
   end

   else:

endcase



end





;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;EVENT HANDLER FOR GENERAL BUTTON EVENTS
pro gridview2_event, event

common gridview2_state
common gridstate


widget_control, event.id, get_uvalue=uvalue



case uvalue of

'velslider': begin

    widget_control, event.id, get_value=channel
    
    

    widget_control, state.zval, set_value=strcompress(grid.velarr[channel])
    widget_control, state.gotochannel, set_value=strcompress(channel, /remove_all)

    ;gridview2_update_spectrum_box

    gridview2_refresh, channel   ; 1 means mouse drag
    gridview2_refresh_weights, channel
    ;if (event.drag eq 0) then gridview2_display, channel ; 0 means end of drag
    

end

'gotochannel': begin

widget_control, event.id, get_value=channelstring

channel=long(channelstring)
channel=channel[0]

if (channel gt n_elements(grid.velarr)-1) then $
                channel=n_elements(grid.velarr)-1
if (channel lt 0) then channel=0

widget_control, state.zval, set_value=strcompress(grid.velarr[channel])

widget_control, state.velslider, set_value=channel
widget_control, state.gotochannel, set_value=strcompress(channel, /remove_all)

gridview2_display, channel
gridview2_display_weights, channel

end


'nextchannel': begin

widget_control, state.velslider, get_value=channel

if (channel eq n_elements(grid.velarr)-1) then channel=n_elements(grid.velarr)-2



widget_control, state.zval, set_value=strcompress(grid.velarr[channel+1])
widget_control, state.gotochannel, set_value=strcompress(channel+1, /remove_all)
widget_control, state.velslider, set_value=channel+1



gridview2_refresh, channel+1
gridview2_refresh_weights, channel+1

end


'previouschannel': begin

widget_control, state.velslider, get_value=channel

if (channel eq 0) then channel=1

widget_control, state.zval, set_value=strcompress(grid.velarr[channel-1])
widget_control, state.gotochannel, set_value=strcompress(channel-1, /remove_all)
widget_control, state.velslider, set_value=channel-1



gridview2_refresh, channel-1
gridview2_refresh_weights, channel-1

end

'polchoice':begin

    state.currentpol=event.index
    case event.index of
        0: state.colortable=state.polAcolor
        1: state.colortable=state.polBcolor
        2: state.colortable=state.avgpolcolor
        else:
    endcase
    widget_control, state.velslider, get_value=channel
    gridview2_display, channel
    gridview2_display_weights, channel
    

end


'pola': begin

    state.currentpol=0
    state.colortable=state.polAcolor
    widget_control, state.velslider, get_value=channel
    gridview2_display, channel
    gridview2_display_weights, channel

end

'polb': begin

    state.currentpol=1
    state.colortable=state.polBcolor
    widget_control, state.velslider, get_value=channel
    gridview2_display, channel
    gridview2_display_weights, channel

end

'avgpol': begin

    state.currentpol=2
    state.colortable=state.avgpolcolor
    widget_control, state.velslider, get_value=channel
    gridview2_display, channel
    gridview2_display_weights, channel

end


'velintegselect': begin

    if (state.velintegselectstatus eq 0) then begin 
          state.velintegselectstatus=1
    endif else begin
          state.velintegselectstatus=0
    endelse

   
end

'multiagc': begin

if (state.multiagcstatus eq 0) then begin
        state.multiagcstatus=1
        Widget_Control, event.id, Set_Button=1
    endif else begin
        state.multiagcstatus=0
        Widget_Control, event.id, Set_Button=0
endelse

widget_control, state.gotochannel, get_value=channelstring

   channel=long(channelstring)
   channel=channel[0]

   gridview2_display, channel
   gridview2_display_weights, channel


end




'agcoverlay': begin

   if (state.agcoverlaystatus eq 0) then begin
        state.agcoverlaystatus=1
        Widget_Control, event.id, Set_Button=1
    endif else begin
        state.agcoverlaystatus=0
        Widget_Control, event.id, Set_Button=0
    endelse

   widget_control, state.gotochannel, get_value=channelstring

   channel=long(channelstring)
   channel=channel[0]

   gridview2_display, channel
   gridview2_display_weights, channel

end

'agcoverlay_nocz': begin

   if (state.agcoverlaystatus_nocz eq 0) then begin
        state.agcoverlaystatus_nocz=1
        Widget_Control, event.id, Set_Button=1
    endif else begin
        state.agcoverlaystatus_nocz=0
        Widget_Control, event.id, Set_Button=0
    endelse

   widget_control, state.gotochannel, get_value=channelstring

   channel=long(channelstring)
   channel=channel[0]

   gridview2_display, channel
   gridview2_display_weights, channel

end


'dhbboverlay': begin

if (state.dhbboverlaystatus eq 0) then begin
        state.dhbboverlaystatus=1
        Widget_Control, event.id, Set_Button=1
    endif else begin
        state.dhbboverlaystatus=0
        Widget_Control, event.id, Set_Button=0
    endelse




end



'catalog3D': begin

if (state.overlaystatus3D eq 0) then begin

    if (state.extracats eq 1) then begin
        state.overlaystatus3D=1
        Widget_Control, event.id, Set_Button=1
    endif else begin
         ;cat3D=dialog_pickfile(/read) 
        
         ;restore, '/home/dorado3/galaxy/idl_alfa/catalog_template.sav' 
         ;table3D=read_ascii(cat3D, template=catalog_template)
         ;convert3Dcat, catin, table3D
         ;state.table3D=table3D
         ;state.extracats=1
         ;state.overlaystatus3D=1

        status=dialog_message('You need to pass a file name as a keyword upon startup.')


    endelse



    endif else begin
        state.overlaystatus3D=0
        Widget_Control, event.id, Set_Button=0
    endelse

   widget_control, state.gotochannel, get_value=channelstring

   channel=long(channelstring)
   channel=channel[0]

   gridview2_display, channel
   gridview2_display_weights, channel


end

'checklist': begin

;check=findfile('checklist.txt', count=count)

if (state.checkliststatus eq 1) then begin
    state.checkliststatus=0
    Widget_Control, event.id, Set_Button=0
endif else begin
    state.checkliststatus=1
    Widget_Control, event.id, Set_Button=1
endelse
    

end

'cntsrcs': begin

if (state.cntsrcstatus eq 1) then begin
    state.cntsrcstatus=0
    Widget_Control, event.id, Set_Button=0
endif else begin
    state.cntsrcstatus=1
    Widget_Control, event.id, Set_Button=1
endelse
    
end



'overlaybox': begin

if (state.overlayboxstatus eq 1) then begin
   state.overlayboxstatus=0
   widget_control, event.id, set_button=0
   state.zoomstatus=1
endif else begin
    state.overlayboxstatus=1
    widget_control, event.id, set_button=1
    state.zoomstatus=0
endelse

end


'clearoverlaybox': begin

widget_control, state.gotochannel, get_value=channelstring

   channel=long(channelstring)
   channel=channel[0]

state.overlaybox[*,*]=0

   gridview2_display, channel
   gridview2_display_weights, channel

end


'nvssoverlay': begin

    if (state.nvssoverlaystatus eq 0) then begin
        state.nvssoverlaystatus=1
        Widget_Control, event.id, Set_Button=1
    endif else begin
        state.nvssoverlaystatus=0
        Widget_Control, event.id, Set_Button=0
    endelse

   widget_control, state.gotochannel, get_value=channelstring

   channel=long(channelstring)
   channel=channel[0]

   gridview2_display, channel
   gridview2_display_weights, channel


end


'getspectrum': begin

  if (state.mapcurrent eq 'spectral') then begin

   if (state.getspectrumstatus eq 0) then begin
        state.getspectrumstatus=1
    endif else begin
        state.getspectrumstatus=0
   endelse

   endif

end

'fluxmeasure': begin

   if (state.mapcurrent eq 'spectral') then begin

      if (state.fluxmeasurestatus eq 0) then begin
        state.fluxmeasurestatus=1
        state.zoomstatus=0      ;turn zoom off
      endif else begin
        state.fluxmeasurestatus=0
        state.zoomstatus=1      ;turn zoom on
      endelse 

   endif



end


'cleanregion': begin
	if (state.mapcurrent eq 'spectral') then begin
		if (state.fluxmeasurestatus eq 0) then begin
			state.fluxmeasurestatus=2
			state.zoomstatus=0
		endif else begin
			state.fluxmeasurestatus=0
			state.zoomstatus=1
		endelse
	endif
end

'axisvelocity': begin

state.spectrumaxisstatus=0

if (state.spectrumon eq 1) then gridview2_plotspectrum

end


'axischannel': begin

state.spectrumaxisstatus=1

if (state.spectrumon eq 1) then gridview2_plotspectrum

end


'spectral': begin
   state.mapcurrent='spectral'
   widget_control, state.gotochannel, get_value=channelstring

   channel=long(channelstring)
   channel=channel[0]

   gridview2_display, channel
   gridview2_display_weights, channel


    widget_control, state.opticalsize, set_value='10'
    widget_control, state.imagelabel, set_value=[' DSS ', 'Sloan']
    state.currentimage=0

end

'continuum': begin

   state.mapcurrent='continuum'
   widget_control, state.gotochannel, get_value=channelstring

   channel=long(channelstring)
   channel=channel[0]
   
   gridview2_display, channel
   gridview2_display_weights, channel

    widget_control, state.opticalsize, set_value='20' 
    widget_control, state.imagelabel, set_value=' NVSS '

end

'linearscaling': begin

widget_control, state.button_linear, set_button=1
widget_control, state.button_log,    set_button=0
widget_control, state.button_histeq, set_button=0

state.currentscaling='linearscaling'



   widget_control, state.gotochannel, get_value=channelstring

   channel=long(channelstring)
   channel=channel[0]
   
   gridview2_display, channel
   gridview2_display_weights, channel

end

'logscaling': begin

widget_control, state.button_linear, set_button=0
widget_control, state.button_log,    set_button=1
widget_control, state.button_histeq, set_button=0

state.currentscaling='logscaling'

widget_control, state.gotochannel, get_value=channelstring

   channel=long(channelstring)
   channel=channel[0]
   
   gridview2_display, channel
   gridview2_display_weights, channel

end


'histeqscaling': begin

widget_control, state.button_linear, set_button=0
widget_control, state.button_log,    set_button=0
widget_control, state.button_histeq, set_button=1

state.currentscaling='histeqscaling'

widget_control, state.gotochannel, get_value=channelstring

   channel=long(channelstring)
   channel=channel[0]
   
   gridview2_display, channel
   gridview2_display_weights, channel


end

'smoothslider': begin

widget_control, state.gotochannel, get_value=channelstring

   channel=long(channelstring)
   channel=channel[0]

gridview2_refresh, channel
gridview2_refresh_weights, channel


end


'colorbar_autoscale': begin

state.colorbarscaling='colorbar_autoscale'
Widget_Control, event.id, Set_Button=1

widget_control, state.button_constantscale, set_button=0

widget_control, state.gotochannel, get_value=channelstring

   channel=long(channelstring)
   channel=channel[0]
   
   gridview2_display, channel
   gridview2_display_weights, channel

end


'colorbar_constantscale': begin

state.colorbarscaling='colorbar_constantscale'
Widget_Control, event.id, Set_Button=1

widget_control, state.button_autoscale, set_button=0


widget_control, state.gotochannel, get_value=channelstring

   channel=long(channelstring)
   channel=channel[0]
   
   gridview2_display, channel
   gridview2_display_weights, channel

end


'imageoption': begin

state.currentimage=event.index

gridview2_cutout_display

end

'plotwindowimage': begin


widget_control, state.opticalsize, get_value=imagesizestring
imagesize=float(imagesizestring[0])


ra=state.currentimage_rahr
if (ra lt 0.0) then ra=ra+24

getdss, ra, state.currentimage_decdeg, imagesize=imagesize


end

'bothpols': begin

spectrum=grid.d[*,*,state.current_xcube,state.current_ycube]
spectrumweights=grid.w[*,*,state.current_xcube,state.current_ycube]

plot_both_pols, spectrum, spectrumweights, grid.velarr, state.spectrum_xmin,$
                state.spectrum_xmax, state.maxweight, state.baseID


end



else:

endcase

end


;-------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;WINDOW display procedure - resets entire window

pro gridview2_display, channel

common gridview2_state
common gridstate

widget_control, state.plotwindowone, get_value=index
wset, index


device, decomposed=0
loadct, state.colortable, /silent
stretch, state.stretchlower,state.stretchupper

ramin=state.ramin
ramax=state.ramax
decmin=state.decmin
decmax=state.decmax

;print, ramin, ramax, decmin, decmax

hor, ramin, ramax
ver, decmin, decmax

posarray=[0.15,0.15,0.95,0.95]
xstyle=1
ystyle=1

color='00FFFF'XL
charsize=1.0
xtitle='!3RA [HMS] J2000'
ytitle='!3Dec [DMS] J2000'
ticklen=-0.01
if (state.currentpol eq 0) then poltitle='Pol A'
if (state.currentpol eq 1) then poltitle='Pol B'
if (state.currentpol eq 2) then poltitle='Avg Pol'

device, decomposed=1 
plot, [0,0], /nodata, xstyle=xstyle, ystyle=ystyle, position=posarray, $
             color=color, charsize=charsize, xtitle=xtitle, ytitle=ytitle, $
             xtick_get=xvals, ytick_get=yvals, ticklen=ticklen

device, decomposed=0

nxticklabels=n_elements(xvals)
nyticklabels=n_elements(yvals)

xvals=xvals
yvals=yvals

xspacing=((xvals[n_elements(xvals)-1]-xvals[0])*60.0)/(nxticklabels-1)
yspacing=((yvals[n_elements(yvals)-1]-yvals[0])*60.0)/(nyticklabels-1)

;ticlabels, xvals[0]*15.0, nxticklabels, xspacing, xticlabs,/ra,delta=1
;ticlabels, yvals[0], nyticklabels, yspacing, yticlabs, delta=1


xticlabs=ratickname(xvals*15.0)
yticlabs=dectickname(yvals)


state.PX = !X.WINDOW * !D.X_VSIZE 
state.PY = !Y.WINDOW * !D.Y_VSIZE 
state.SX = state.PX[1] - state.PX[0] + 1 
state.SY = state.PY[1] - state.PY[0] + 1

erase

hor, ramax, ramin
ver, decmin, decmax
map=dblarr(state.sx, state.sy)


case state.mapcurrent of

'spectral': begin

if (state.velintegselectstatus eq 1) then begin
    widget_control, state.velintegtext, get_value=step
    step=long(step[0])
    lowerchan=channel-step
    upperchan=channel+step
    if (lowerchan lt 0) then lowerchan=0
    if (upperchan gt n_elements(grid.velarr)-1) then upperchan=n_elements(grid.velarr)-1
    if (state.currentpol eq 0 OR state.currentpol eq 1) then begin

       map=total(grid.d[lowerchan:upperchan,state.currentpol,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax],1)/(float(2*step+1))
               
    endif else begin
       
        map=reform(total(grid.d[lowerchan:upperchan,*,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax],1)/(float(2*step+1)))
        map=total(map,1)/2.
 
        
      

    endelse

      infostring='Channels '+strcompress(lowerchan, /remove_all)+' to '+strcompress(upperchan, /remove_all)+ $
             '  Smoothed over -> '+strcompress(grid.velarr[lowerchan], /remove_all)+$
             ' to '+strcompress(grid.velarr[upperchan], /remove_all)+' km/s   '+poltitle
       
 
endif else begin  

        if (state.currentpol eq 0 OR state.currentpol eq 1) then begin
          map=reform(grid.d[channel,state.currentpol,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax])
        endif else begin
          map=reform(grid.d[channel,*,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax])
          map=total(map,1)/2.
          
        endelse

infostring='Channel '+strcompress(channel)+ $
             '    Velocity= '+strcompress(grid.velarr[channel])+' km/s   '+poltitle

endelse


end

'continuum': begin

if (state.currentpol eq 0 OR state.currentpol eq 1) then begin
          map=reform(grid.cont[state.currentpol,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax])
        endif else begin
          map=reform(grid.cont[*,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax])
          map=total(map,1)/2.
          
        endelse

infostring='Continuum    '+poltitle

end


else:

endcase


gridview2_displaymap, reverse(congrid(reform(map), state.sx, state.sy))


device, decomposed=1
plot, [0,0], /nodata, xstyle=1, ystyle=1, /noerase, position=posarray, $ 
             xtickn=reverse(xticlabs), ytickn=yticlabs, $
             color=color, charsize=charsize, xtitle=xtitle, ytitle=ytitle, ticklen=ticklen
             



mask=fltarr(600,25)
tvscl, mask, 10,5

xyouts, 10, 10, infostring, /device, color=color, charsize=1.1

device, decomposed=0

;Get region and "s" info for plotting parameters - IMPORTANT!!
state.xs=!x.s
state.ys=!y.s
state.xregion=!x.region
state.yregion=!y.region


gridview2_loadColorBar

wset, index

gridview2_agcdisplay
gridview2_nvss

end



;----------------------------------------------------------------------------
;Resets weights window

pro gridview2_display_weights, channel

common gridview2_state
common gridstate

widget_control, state.plotwindowtwo, get_value=index
wset, index

device, decomposed=0
loadct, state.colortable, /silent
stretch, state.stretchlower,state.stretchupper

ramin=state.ramin
ramax=state.ramax
decmin=state.decmin
decmax=state.decmax

hor, ramin, ramax
ver, decmin, decmax

posarray=[0.15,0.15,0.95,0.95]
xstyle=1
ystyle=1


color='00FFFF'XL
charsize=1.0
xtitle='RA [HMS] J2000'
ytitle='Dec [DMS] J2000'
ticklen=-0.01
if (state.currentpol eq 0) then poltitle='Pol A'
if (state.currentpol eq 1) then poltitle='Pol B'
if (state.currentpol eq 2) then poltitle='Avg Pol'

device, decomposed=1 
plot, [0,0], /nodata, xstyle=xstyle, ystyle=ystyle, position=posarray, $
             color=color, charsize=charsize, xtitle=xtitle, ytitle=ytitle, $
             xtick_get=xvals, ytick_get=yvals, ticklen=ticklen



device, decomposed=0

nxticklabels=n_elements(xvals)
nyticklabels=n_elements(yvals)



xvals=xvals
yvals=yvals

xspacing=((xvals[n_elements(xvals)-1]-xvals[0])*60.0)/(nxticklabels-1)
yspacing=((yvals[n_elements(yvals)-1]-yvals[0])*60.0)/(nyticklabels-1)

;ticlabels, xvals[0]*15.0, nxticklabels, xspacing, xticlabs,/ra,delta=1
;ticlabels, yvals[0], nyticklabels, yspacing, yticlabs, delta=1

xticlabs=ratickname(xvals*15.0)
yticlabs=dectickname(yvals)


state.PX = !X.WINDOW * !D.X_VSIZE 
state.PY = !Y.WINDOW * !D.Y_VSIZE 
state.SX = state.PX[1] - state.PX[0] + 1 
state.SY = state.PY[1] - state.PY[0] + 1

erase

hor, ramax, ramin
ver, decmin, decmax


case state.mapcurrent of

'spectral':  begin

if (state.velintegselectstatus eq 1) then begin
    widget_control, state.velintegtext, get_value=step
    step=long(step[0])
    lowerchan=channel-step
    upperchan=channel+step
    if (lowerchan lt 0) then lowerchan=0
    if (upperchan gt n_elements(grid.velarr)-1) then upperchan=n_elements(grid.velarr)-1
    if (state.currentpol eq 0 OR state.currentpol eq 1) then begin

       map=total(grid.w[lowerchan:upperchan,state.currentpol,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax],1)/(float(2*step))
               
    endif else begin
       
        map=reform(total(grid.w[lowerchan:upperchan,*,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax],1)/(float(2*step)))
        map=total(map,1)/2.
       

    endelse

      
      
       
 
endif else begin  

        if (state.currentpol eq 0 OR state.currentpol eq 1) then begin
          map=reform(grid.w[channel,state.currentpol,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax])
        endif else begin
          map=reform(grid.w[channel,*,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax])
          map=total(map,1)/2.
         
        endelse


    endelse

   infostring='Spectral weights display'

end

'continuum':  begin

if (state.currentpol eq 0 OR state.currentpol eq 1) then begin
          map=reform(grid.cw[state.currentpol,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax])
        endif else begin
          map=reform(grid.cw[*,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax])
          map=total(map,1)/2.
          
        endelse

   infostring='Continuum weights display'


end 

else:

endcase

;Next line removed to keep linear scaling for weights
;gridview2_displaymap, reverse(congrid(reform(map), state.sx, state.sy))

widget_control, state.smoothslider, get_value=smoothval
smoothval=smoothval[0]
;tvscl, smooth(reverse(congrid(reform(map), state.sx, state.sy)),$
;         smoothval, /edge_truncate), state.px[0], state.py[0]

tv, bytscl(smooth(reverse(congrid(reform(map), state.sx, state.sy)),smoothval, /edge_truncate), min=0, max=30), state.px[0], state.py[0]


device, decomposed=1
plot, [0,0], /nodata, xstyle=1, ystyle=1, /noerase, position=posarray, $ 
             xtickn=reverse(xticlabs), ytickn=yticlabs, $
             color=color, charsize=charsize, xtitle=xtitle, ytitle=ytitle, ticklen=ticklen
            

 xyouts, 225,15, infostring, /device, color=color, charsize=1.3

end



;---------------------------------------------------------------------------
;Refresh only the images - faster and less flickering

pro gridview2_refresh, channel

common gridview2_state
common gridstate

widget_control, state.plotwindowone, get_value=index
wset, index

ramin=state.ramin
ramax=state.ramax
decmin=state.decmin
decmax=state.decmax

hor, ramax, ramin
ver, decmin, decmax

;Reset region and "s" info for plotting parameters to state values- IMPORTANT!!
!x.s=state.xs
!y.s=state.ys
!x.region=state.xregion
!y.region=state.yregion

device, decomposed=0
loadct, state.colortable, /silent
stretch, state.stretchlower,state.stretchupper

if (state.currentpol eq 0) then poltitle='Pol A'
if (state.currentpol eq 1) then poltitle='Pol B'
if (state.currentpol eq 2) then poltitle='Avg Pol'


case state.mapcurrent of

'spectral': begin

if (state.velintegselectstatus eq 1) then begin
    widget_control, state.velintegtext, get_value=step
    step=long(step[0])
    lowerchan=channel-step
    upperchan=channel+step
    if (lowerchan lt 0) then lowerchan=0
    if (upperchan gt n_elements(grid.velarr)-1) then upperchan=n_elements(grid.velarr)-1
    if (state.currentpol eq 0 OR state.currentpol eq 1) then begin

       map=total(grid.d[lowerchan:upperchan,state.currentpol,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax],1)/(float(2*step+1))
               
    endif else begin
       
        map=reform(total(grid.d[lowerchan:upperchan,*,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax],1)/(float(2*step+1)))
        map=total(map,1)/2.
       

    endelse

      
        infostring='Channels '+strcompress(lowerchan, /remove_all)+' to '+strcompress(upperchan, /remove_all)+ $
             '  Smoothed over -> '+strcompress(grid.velarr[lowerchan], /remove_all)+$
             ' to '+strcompress(grid.velarr[upperchan], /remove_all)+' km/s   '+poltitle
       
 
endif else begin  

        if (state.currentpol eq 0 OR state.currentpol eq 1) then begin
          map=reform(grid.d[channel,state.currentpol,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax])
        endif else begin
          map=reform(grid.d[channel,*,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax])
          map=total(map,1)/2.
        endelse

infostring='!3Channel '+strcompress(channel)+ $
             '    Velocity= '+strcompress(grid.velarr[channel])+' km/s   '+poltitle

endelse

end

'continuum': begin

if (state.currentpol eq 0 OR state.currentpol eq 1) then begin
          map=reform(grid.cont[state.currentpol,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax])
        endif else begin
          map=reform(grid.cont[*,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax])
          map=total(map,1)/2.
          
        endelse

infostring='Continuum    '+poltitle

end

else:

endcase



gridview2_displaymap, reverse(congrid(reform(map), state.sx, state.sy))


color='00FFFF'XL    ;YELLOW

device, decomposed=1
plots, findgen(long(state.sx))+state.px[0], $
       fltarr(long(state.sx))+state.px[0], color=color, /device
plots, fltarr(long(state.sx))+state.px[0],  $
       findgen(long(state.sx))+state.px[0], color=color, /device

plots, findgen(long(state.sx))+state.px[0], $
       fltarr(long(state.sx))+state.px[0]+state.sx-1, color=color, /device
plots, fltarr(long(state.sx))+state.px[0]+state.sx-1, $
       findgen(long(state.sx))+state.px[0], color=color, /device

mask=fltarr(600,25)
tvscl, mask, 10,5

xyouts, 10, 10, infostring, /device, color=color, charsize=1.1

device, decomposed=0

gridview2_loadColorBar

wset, index

gridview2_agcdisplay
gridview2_nvss

end


;---------------------------------------------------------------------------
;Refresh only the images - faster and less flickering

pro gridview2_refresh_weights, channel

common gridview2_state
common gridstate

widget_control, state.plotwindowtwo, get_value=index
wset, index

ramin=state.ramin
ramax=state.ramax
decmin=state.decmin
decmax=state.decmax

hor, ramax, ramin
ver, decmin, decmax

;Reset region and "s" info for plotting parameters to state values- IMPORTANT!!
!x.s=state.xs
!y.s=state.ys
!x.region=state.xregion
!y.region=state.yregion


device, decomposed=0
loadct, state.colortable, /silent
stretch, state.stretchlower,state.stretchupper


if (state.currentpol eq 0) then poltitle='Pol A'
if (state.currentpol eq 1) then poltitle='Pol B'
if (state.currentpol eq 2) then poltitle='Avg Pol'


case state.mapcurrent of

'spectral': begin

if (state.velintegselectstatus eq 1) then begin
    widget_control, state.velintegtext, get_value=step
    step=long(step[0])
    lowerchan=channel-step
    upperchan=channel+step
    if (lowerchan lt 0) then lowerchan=0
    if (upperchan gt n_elements(grid.velarr)-1) then upperchan=n_elements(grid.velarr)-1
    if (state.currentpol eq 0 OR state.currentpol eq 1) then begin

       map=total(grid.w[lowerchan:upperchan,state.currentpol,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax],1)/(float(2*step+1))
               
    endif else begin
       
        map=reform(total(grid.w[lowerchan:upperchan,*,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax],1)/(float(2*step+1)))
        map=total(map,1)/2.

    endelse

endif else begin  

        if (state.currentpol eq 0 OR state.currentpol eq 1) then begin
          map=reform(grid.w[channel,state.currentpol,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax])
        endif else begin
          map=reform(grid.w[channel,*,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax])
          map=total(map,1)/2.
        
        endelse

    endelse

end

'continuum':  begin

if (state.currentpol eq 0 OR state.currentpol eq 1) then begin
          map=reform(grid.cw[state.currentpol,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax])
        endif else begin
          map=reform(grid.cw[*,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax])
          map=total(map,1)/2.
          
        endelse

end

else:
 
endcase


;gridview2_displaymap, reverse(congrid(reform(map), state.sx, state.sy))
widget_control, state.smoothslider, get_value=smoothval
smoothval=smoothval[0]
;tvscl, smooth(reverse(congrid(reform(map), state.sx, state.sy)),smoothval, /edge_truncate), state.px[0], state.py[0]
tv, bytscl(smooth(reverse(congrid(reform(map), state.sx, state.sy)),smoothval, /edge_truncate), min=0, max=30), state.px[0], state.py[0]


color='00FFFF'XL  ;YELLOW

device, decomposed=1
plots, findgen(long(state.sx))+state.px[0], fltarr(long(state.sx))+state.px[0], color=color, /device
plots, fltarr(long(state.sx))+state.px[0], findgen(long(state.sx))+state.px[0], color=color, /device

plots, findgen(long(state.sx))+state.px[0], fltarr(long(state.sx))+state.px[0]+state.sx-1, color=color, /device
plots, fltarr(long(state.sx))+state.px[0]+state.sx-1, findgen(long(state.sx))+state.px[0], color=color, /device

device, decomposed=0

end

;-----------------------------------------------------------------------
;Procedure that is passed map, and map (spectral, continuum or
;weights) is shown, with appropriate scaling applied

pro gridview2_displaymap, map

common gridview2_state
common gridstate

;Fetch the smoothing state
widget_control, state.smoothslider, get_value=smoothval
smoothval=smoothval[0]


case state.currentscaling of

   'linearscaling': begin

     if (state.colorbarscaling eq 'colorbar_autoscale') then begin
        tvscl, smooth(map,smoothval, /edge_truncate), state.px[0], state.py[0]
     endif else begin
        tv, bytscl(smooth(map,smoothval,/edge_truncate),min=state.colorbarmin,max=state.colorbarmax),$
             state.px[0], state.py[0]   
      endelse

     state.minval=min(map)
     state.maxval=max(map)

    end

    'logscaling': begin
       offset=min(map)-(max(map)-min(map))*0.01
       tvscl, smooth(alog10(map-offset), smoothval, /edge_truncate), state.px[0], state.py[0] 
       state.minval=alog10(min(map)-offset)
       state.maxval=alog10(max(map)-offset)

    end

    'histeqscaling':begin

        tvscl, smooth(hist_equal(map, minv=min(map), maxv=max(map)), smoothval, /edge_truncate), state.px[0], state.py[0] 
        state.minval=min(map)
        state.maxval=max(map)

    end

else:

endcase

end




;-------------------------------------------------------------
;Load colorbar

pro gridview2_loadColorBar

common gridview2_state
common gridstate

widget_control, state.gotochannel, get_value=channelstring
channel=channelstring[0]

widget_control, state.plotwindowcolorbar, get_value=index
wset, index

erase

device, decomposed=0
loadct, state.colortable, /silent
stretch, state.stretchlower,state.stretchupper

yellow='00FFFF'XL


device, decomposed=1

if (state.maxval eq 0.0) then state.maxval=10.0

if (state.colorbarscaling eq 'colorbar_autoscale') then begin
Colorbar, range=[state.minval,state.maxval], vertical=1, position=[0.5,0.15,0.9,0.95], ticklen=-0.05, color=yellow
endif else begin
Colorbar, range=[state.colorbarmin,state.colorbarmax], vertical=1, $
        position=[0.5,0.15,0.9,0.95], ticklen=-0.05, color=yellow
endelse

xyouts, 12,50, 'mJy/', color=yellow, /device
xyouts, 15,40, 'beam', color=yellow, /device

device, decomposed=0




end


;-----------------------------------------------------
;Overlay AGC galaxies if user so chooses
pro gridview2_agcdisplay

common gridview2_state
common gridstate


ramin=state.ramin
ramax=state.ramax
decmin=state.decmin
decmax=state.decmax


widget_control, state.gotochannel, get_value=channelstring

channel=long(channelstring)
channel=channel[0] 


    lowerchan=channel-14
    upperchan=channel+14


if (lowerchan lt 0) then lowerchan=0
if (upperchan gt n_elements(grid.velarr)-1) then upperchan=n_elements(grid.velarr)-1


if (state.agcoverlaystatus eq 1 AND state.mapcurrent eq 'spectral') then begin

iagc=where(state.rahr lt ramax AND $
           state.rahr gt ramin AND $
           state.decdeg lt decmax AND $
           state.decdeg gt decmin AND $
           ((state.agc.v21 lt grid.velarr[channel]+state.agc.width/2.0 AND $
             state.agc.v21 gt grid.velarr[channel]-state.agc.width/2.0 AND $
             state.agc.v21 ne 0) OR $
            (state.agc.vopt lt grid.velarr[lowerchan] AND $
             state.agc.vopt gt grid.velarr[upperchan] AND $
             state.agc.vopt ne 0)))




color='00FFFF'XL    ;YELLOW squares

hor, ramax, ramin
ver, decmax, decmin



if (iagc[0] ne -1) then begin
device, decomposed=1
resultcoords=convert_coord(state.rahr[iagc], state.decdeg[iagc], /data, /double, /to_device)

plots, state.rahr[iagc], state.decdeg[iagc], psym=6,symsize=2, color=color, /data
xyouts,resultcoords[0,*]+10, resultcoords[1,*]-10, $
             strcompress(state.agc.agcnumber[iagc], /remove_all), color=color, /device

for i=0,n_elements(resultcoords[0,*])-1 do state.currentagc[*,i]=[state.agc.agcnumber[iagc[i]],resultcoords[0,i], resultcoords[1,i]]

device, decomposed=0
endif else begin
   state.currentagc=lonarr(3,5000)
endelse

endif


;Plot AGC galaxies with no known redshifts in the area

if (state.agcoverlaystatus_nocz eq 1 AND state.mapcurrent eq 'spectral') then begin

iagc=where(state.rahr lt ramax AND $
           state.rahr gt ramin AND $
           state.decdeg lt decmax AND $
           state.decdeg gt decmin AND $
           state.agc.v21 eq 0 AND $
           state.agc.vopt eq 0 )



color='0000FF'XL    ;RED DIAMONDS

hor, ramax, ramin
ver, decmax, decmin

if (iagc[0] ne -1) then begin
device, decomposed=1
resultcoords=convert_coord(state.rahr[iagc], state.decdeg[iagc], /data, /double, /to_device)

plots, state.rahr[iagc], state.decdeg[iagc], psym=4,symsize=3, color=color, /data
xyouts,resultcoords[0,*]+10, resultcoords[1,*]-10, $
             strcompress(state.agc.agcnumber[iagc], /remove_all), color=color, /device



device, decomposed=0
endif

endif




;-------------------------------------------------------------------

if (state.overlaystatus2D eq 1) then begin

;;;;;FOR BRIAN KENT ONLY - VIRGO SOUTH STUFF - NOT GENERALIZED!!

;Load from 2D extractions



i2D=where(state.table2D.ra lt ramax AND $
           state.table2D.ra gt ramin AND $
           state.table2D.dec lt decmax AND $
           state.table2D.dec gt decmin AND $
           state.table2D.velocity lt grid.velarr[lowerchan] AND $
           state.table2D.velocity gt grid.velarr[upperchan])



color='00FF00'XL    ;GREEN Triangle

hor, ramax, ramin
ver, decmax, decmin

if (i2D[0] ne -1) then begin
device, decomposed=1
resultcoords=convert_coord(state.table2D[i2D].ra, state.table2D[i2D].dec, /data, /double, /to_device)

plots, state.table2D[i2D].ra, state.table2D[i2D].dec, psym=5,symsize=2, color=color, /data
xyouts,resultcoords[0,*]-25, resultcoords[1,*]-10, $
             strcompress(state.table2D[i2D].sourcenumber, /remove_all), color=color, /device

device, decomposed=0
endif


endif




if (state.overlaystatus3D eq 1 AND state.mapcurrent eq 'spectral' AND state.extracats eq 1) then begin

;Load from 3D extractions

;modified to accept a larger width
i3D=where(state.table3D.ra lt ramax AND $
           state.table3D.ra gt ramin AND $
           state.table3D.dec lt decmax AND $
           state.table3D.dec gt decmin AND $
           state.table3D.cz lt grid.velarr[lowerchan] AND $
           state.table3D.cz gt grid.velarr[upperchan])



color='0000FF'XL    ;RED Plus sign

hor, ramax, ramin
ver, decmax, decmin

if (i3D[0] ne -1) then begin
device, decomposed=1
resultcoords=convert_coord(state.table3D.ra[i3D], state.table3D.dec[i3D], /data, /double, /to_device)

plots, state.table3D.ra[i3D], state.table3D.dec[i3D], psym=1,symsize=2, color=color, thick=2.0, /data
xyouts,resultcoords[0,*]-25, resultcoords[1,*]+10, $
             strcompress(state.table3D.id[i3D], /remove_all), color=color, /device

device, decomposed=0
endif

endif

;-------------------------------------
;New multi-color agc plotting

if (state.multiagcstatus eq 1 AND state.mapcurrent eq 'spectral') then begin
symsize=1.5


;Plot All AGC galaxies in the range in red

iagc=where(state.rahr lt ramax AND $
           state.rahr gt ramin AND $
           state.decdeg lt decmax AND $
           state.decdeg gt decmin)

color='0000FF'XL    ;RED squares

hor, ramax, ramin
ver, decmax, decmin

if (iagc[0] ne -1) then begin
device, decomposed=1
resultcoords=convert_coord(state.rahr[iagc], state.decdeg[iagc], /data, /double, /to_device)

plots, state.rahr[iagc], state.decdeg[iagc], psym=6,symsize=symsize, color=color, /data
xyouts,resultcoords[0,*]+10, resultcoords[1,*]-10, $
             strcompress(state.agc.agcnumber[iagc], /remove_all), color=color, /device

for i=0,n_elements(resultcoords[0,*])-1 do state.currentagc[*,i]=[state.agc.agcnumber[iagc[i]],resultcoords[0,i], resultcoords[1,i]]

device, decomposed=0
endif else begin
   state.currentagc=lonarr(3,5000)
endelse


;Plot those that have known redshift or have been searched in HI,
;yielding positive OR negative results in Cyan

iagc=where(state.rahr lt ramax AND $
           state.rahr gt ramin AND $
           state.decdeg lt decmax AND $
           state.decdeg gt decmin AND $
          (state.agc.v21 ne 0 OR state.agc.vopt ne 0))

color='FFFF00'XL    ;CYAN squares

hor, ramax, ramin
ver, decmax, decmin

if (iagc[0] ne -1) then begin
device, decomposed=1
resultcoords=convert_coord(state.rahr[iagc], state.decdeg[iagc], /data, /double, /to_device)

plots, state.rahr[iagc], state.decdeg[iagc], psym=6,symsize=symsize, color=color, /data
xyouts,resultcoords[0,*]+10, resultcoords[1,*]-10, $
             strcompress(state.agc.agcnumber[iagc], /remove_all), color=color, /device
endif

;Plot those that are in the redshift range in yellow

iagc=where(state.rahr lt ramax AND $
           state.rahr gt ramin AND $
           state.decdeg lt decmax AND $
           state.decdeg gt decmin AND $
           ((state.agc.v21 lt grid.velarr[channel]+state.agc.width/2.0 AND $
             state.agc.v21 gt grid.velarr[channel]-state.agc.width/2.0 AND $
             state.agc.v21 ne 0) OR $
            (state.agc.vopt lt grid.velarr[lowerchan] AND $
             state.agc.vopt gt grid.velarr[upperchan] AND $
             state.agc.vopt ne 0)))


color='00FFFF'XL    ;YELLOW squares

hor, ramax, ramin
ver, decmax, decmin

if (iagc[0] ne -1) then begin
device, decomposed=1
resultcoords=convert_coord(state.rahr[iagc], state.decdeg[iagc], /data, /double, /to_device)

plots, state.rahr[iagc], state.decdeg[iagc], psym=6,symsize=symsize, color=color, /data
xyouts,resultcoords[0,*]+10, resultcoords[1,*]-10, $
             strcompress(state.agc.agcnumber[iagc], /remove_all), color=color, /device

endif



endif

;-----------------------------------------
;Plot de Heij et al HVC catalog

if (state.dhbboverlaystatus eq 1 AND state.mapcurrent eq 'spectral') then begin

idhbb=where(state.dhbb._raj2000/15.0 lt ramax AND $
            state.dhbb._raj2000/15.0 gt ramin AND $
            state.dhbb._dej2000 lt decmax AND $
            state.dhbb._dej2000 gt decmin AND $
            (state.dhbb_vhelio lt grid.velarr[channel]+state.dhbb.fwhm/2.0 AND $
             state.dhbb_vhelio gt grid.velarr[channel]-state.dhbb.fwhm/2.0))




color='0000FF'XL    ;Red Asterisks

hor, ramax, ramin
ver, decmax, decmin

if (idhbb[0] ne -1) then begin
device, decomposed=1
resultcoords=convert_coord(state.dhbb._raj2000[idhbb]/15.0, state.dhbb._dej2000[idhbb], /data, /double, /to_device)

plots, state.dhbb._raj2000[idhbb]/15.0, state.dhbb._dej2000[idhbb], psym=6,symsize=2, color=color, /data
xyouts,resultcoords[0,*]+10, resultcoords[1,*]-10, $
             'HVC '+strcompress(state.dhbb.seq[idhbb], /remove_all), color=color, /device

;for i=0,n_elements(resultcoords[0,*])-1 do state.currentagc[*,i]=[state.agc.agcnumber[iagc[i]],resultcoords[0,i], resultcoords[1,i]]

device, decomposed=0
endif else begin
   state.currentdhbb=lonarr(3,1000)
endelse

endif


if (state.checkliststatus eq 1 AND state.mapcurrent eq 'spectral') then begin

 check=findfile('checklist.txt', count=count)

    if (count eq 1) then begin
       data=read_ascii('checklist.txt', delimiter=',')
       ;help, data
       rahr=data.field1[0,*]
       decdeg=data.field1[1,*]
       vel=data.field1[2,*]
       width=data.field1[3,*]
       
   
     icheck=where(rahr lt ramax AND $
           rahr gt ramin AND $
           decdeg lt decmax AND $
           decdeg gt decmin AND $
           (vel lt grid.velarr[channel]+1.2*width/2.0 AND $
             vel gt grid.velarr[channel]-1.2*width/2.0))

          if (icheck[0] ne -1) then begin

              device, decomposed=1
              plots, rahr[icheck], decdeg[icheck], psym=4, symsize=2, thick=2, color='0000FF'XL, /data
              device, decomposed=0


          endif                 ;icheck if statement
           
           
       

     endif                     ;count if statement




 endif      ;checkliststatus if statement


;over plot strong continuum sources

if (state.cntsrcstatus eq 1 AND (state.mapcurrent eq 'spectral' OR state.mapcurrent eq 'continuum')) then begin


      
      ramin=state.ramin
      ramax=state.ramax
      decmin=state.decmin
      decmax=state.decmax



;Peak ranges in Jy/beam
peaklower=[.25]    
peakupper=[1.0e6]
peakstring=['(5-50)', '(50-100)', '(100-200)', '(200-300)', '(>300)']



for i=0, n_elements(peakupper)-1 do begin

   invss=where(state.nvsscat.ra_2000_ lt ramax*15.0 AND $
           state.nvsscat.ra_2000_ gt ramin*15.0 AND $
           state.nvsscat.dec_2000_ lt decmax AND $
           state.nvsscat.dec_2000_ gt decmin AND $
           state.nvsscat.peak_int ge peaklower[i] AND $
           state.nvsscat.peak_int lt peakupper[i]) 

   hor, ramax, ramin
   ver, decmax, decmin

   ;overplot white squares

  

   if (invss[0] ne -1) then begin
      device, decomposed=1
      plots,state.nvsscat[invss].ra_2000_/15.0, state.nvsscat[invss].dec_2000_ , $
         psym=6,symsize=2.0,thick=2.0, color='FFFFFF'XL, /data
      device, decomposed=0

   endif  


 endfor  


 endif      ;cntsrcstatus if statement




;Plot overlay boxes if any are set in ra/dec (hr/deg) decimal coordinates

boxindex=where(state.overlaybox[*,0] eq 1)

if (boxindex[0] ne -1) then begin

    for i=0, n_elements(boxindex)-1 do begin
        sx=state.overlaybox[i,1]
        sy=state.overlaybox[i,2]
        dx=state.overlaybox[i,3]
        dy=state.overlaybox[i,4]
        plotS, [sx, sx, dx, dx, sx], [sy, dy, dy, sy, sy], /data, $
          Color='FFFFFF'XL, thick=1.5 ;WHITE 
    endfor


endif


end

;---------------------------------------------------------------------------------
;----Overlay sources from the NVSS Catalog (Condon et al. 1998)
pro gridview2_nvss

common gridview2_state
common gridstate

;Will plot in spectral or continuum mode

if (state.nvssoverlaystatus eq 1) then begin

ramin=state.ramin
ramax=state.ramax
decmin=state.decmin
decmax=state.decmax



;Peak ranges in mJy/beam
peaklower=[0.005, 0.05, 0.1, 0.2, 0.3]    
peakupper=[0.05,   0.1, 0.2, 0.3, 100]
peakstring=['(5-50)', '(50-100)', '(100-200)', '(200-300)', '(>300)']

;CYAN, GREEN, MAGENTA, PINK, YELLOW, RED
colors=['FFFF00'XL, '00FF00'XL,'FF00FF'XL,'7F7FFF'XL, '0000FF'XL]

symsize=[0.9,1.2,2.0,2.0,2.2]
psym=[1,1,5,5,6]
xlevs=[50,130,210,290,370]
ylevs=480

for i=0, n_elements(peakupper)-1 do begin




   invss=where(state.nvsscat.ra_2000_ lt ramax*15.0 AND $
           state.nvsscat.ra_2000_ gt ramin*15.0 AND $
           state.nvsscat.dec_2000_ lt decmax AND $
           state.nvsscat.dec_2000_ gt decmin AND $
           state.nvsscat.peak_int ge peaklower[i] AND $
           state.nvsscat.peak_int lt peakupper[i])

if (ramin lt 0.0) then begin

    invss=where((state.nvsscat.ra_2000_ lt ramax*15.0 AND $
           state.nvsscat.ra_2000_ gt ramin*15.0 AND $
           state.nvsscat.dec_2000_ lt decmax AND $
           state.nvsscat.dec_2000_ gt decmin AND $
           state.nvsscat.peak_int ge peaklower[i] AND $
           state.nvsscat.peak_int lt peakupper[i]) OR $
               (state.nvsscat.ra_2000_ gt (ramin+24.0D)*15.0 AND $
                state.nvsscat.dec_2000_ lt decmax AND $
                state.nvsscat.dec_2000_ gt decmin AND $
                state.nvsscat.peak_int ge peaklower[i] AND $
                state.nvsscat.peak_int lt peakupper[i]))    

endif

 
   hor, ramax, ramin
   ver, decmax, decmin

   if (invss[0] ne -1) then begin
      device, decomposed=1
      plots,state.nvsscat[invss].ra_2000_/15.0, state.nvsscat[invss].dec_2000_ , $
         psym=psym[i],symsize=symsize[i], color=colors[i], /data
      device, decomposed=0

      if (ramin lt 0.0) then begin
          device, decomposed=1
      plots,(state.nvsscat[invss].ra_2000_/15.0)-24.0, state.nvsscat[invss].dec_2000_ , $
         psym=psym[i],symsize=symsize[i], color=colors[i], /data
      device, decomposed=0
      print, (state.nvsscat[invss].ra_2000_/15.0)
      endif

  endif





;Display Legend
device, decomposed=1
plots, xlevs[i], ylevs+3, psym=psym[i], symsize=1.0, color=colors[i], /device
xyouts,xlevs[i]+8, ylevs, peakstring[i], color=colors[i], charsize=1.0, /device
device, decomposed=0

endfor


xyouts, 425, ylevs, 'mJy/beam', charsize=1.0, color='FFFFFF'XL, /device



endif

end



;-----------------------------------------------------
;Exit event handler
pro gridview2_exit, event

common gridview2_state
common gridstate
hor
ver
!p.multi=0

widget_control, state.baseID, /destroy
loadct, 0, /silent
stretch, 0,100

delvarx, state
delvarx, grid

print, 'Exiting Gridview2...'

end


;------------------------------------------------------------
;Help
pro gridview2_help, event

   common gridview2_state
common gridstate

h=['GRIDview2        ', $
   'Written, December 2007', $
   'Based off gridview', $
   ' ', $
   'Last update, Wednesday, January 30, 2008']


if (not (xregistered('gridview2_help', /noshow))) then begin

helptitle = strcompress('Gridview HELP')

    help_base =  widget_base(group_leader = state.baseID, $
                             /column, /base_align_right, title = helptitle, $
                             uvalue = 'help_base')

    help_text = widget_text(help_base, /scroll, value = h, xsize = 85, ysize = 15)
    
    help_done = widget_button(help_base, value = ' Done ', uvalue = 'help_done')

    widget_control, help_base, /realize
    xmanager, 'gridview2_help', help_base, /no_block
    
endif


end

;----------------------------------------------------------------------

pro gridview2_help_event, event

widget_control, event.id, get_uvalue = uvalue

case uvalue of
    'help_done': widget_control, event.top, /destroy
    else:
endcase

end

;----------------------------------------------------------------------

pro gridview2_utility
  common gridview2_state
common gridstate

if (not (xregistered('gridview2_utility', /noshow))) then begin

utiltitle = strcompress('...')

    util_base =  widget_base(group_leader = state.baseID, $
                             /column, /base_align_right, title = utiltitle, $
                             uvalue = 'util_base')

    util_draw = widget_draw(util_base, xsize=800, ysize=600)
    
    util_done = widget_button(util_base, value = ' Do you LOVE your data? ', uvalue = 'util_done')

    widget_control, util_base, /realize
    xmanager, 'gridview2_utility', util_base, /no_block
    
    widget_control, util_draw, get_value=index
    wset, index
    url='http://www.astro.cornell.edu/~bkent/images/datautil'
    spawn, 'wget -q -O ~/datautil3754 ' + "'" + url + "'"
    read_jpeg, '~/datautil3754', testjunk
    tvscl, testjunk, true=1
    spawn, '/bin/rm -r ~/datautil3754'
endif 



end

;----------------------
pro gridview2_utility_event, event

widget_control, event.id, get_uvalue = uvalue

case uvalue of
    'util_done': widget_control, event.top, /destroy
    else:
endcase

end



;-----------------------------------------------------------
;cutout Image display procedure
pro gridview2_cutout_display

common gridview2_state
common gridstate

widget_control, state.plotwindowimage, get_value=index
wset, index
device, decomposed=0
loadct, 1, /silent     


case state.currentimage of 

0: begin
              ;USE FOR DSS2 BLUE
             
               stretch, state.stretchupper, state.stretchlower
               tvscl, state.dssimage
               stretch, state.stretchlower, state.stretchupper

   end


1: begin

               ;USE TO GET SDSS IMAGES

               tv,state.sloanimage, true=1 
               osfamily = strupcase(!version.os_family)
               if (osfamily ne 'UNIX') then begin
                   xyouts, 100,100, 'Available only on Linux',alignment=0.5, /device


               endif    


               
   end

else:

endcase






end



;;;;;;EVENT HANDLERS FOR SETTINGS

;------------------------------------------------------------
;Settings for Polarization colors
pro gridview2_polsettings, event

   common gridview2_state
common gridstate

if (not (xregistered('gridview2_polsettings', /noshow))) then begin

colortablenames=[ 'B-W LINEAR', 'BLUE/WHITE', 'RED TEMPERATURE', 'RAINBOW', 'GREEN/WHITE LINEAR', 'STERN SPECIAL']
colortableselections=['0','1','3','13','8','14']


settingstitle = strcompress('Gridview Pol Settings')

    settings_base =  widget_base(group_leader = state.baseID, $
                             /column, /base_align_right, title = settingstitle, $
                             uvalue = 'settings_base')
     PolAbase=widget_base(settings_base, /row)
     label=widget_label(PolAbase, value=' Polarization A ')
     state.polAselect=widget_droplist(PolAbase, value=colortablenames, uvalue=colortableselections)

     PolBbase=widget_base(settings_base, /row)
     label=widget_label(PolBbase, value=' Polarization B ')
     state.polBselect=widget_droplist(PolBbase, value=colortablenames, uvalue=colortableselections)

     avgpolbase=widget_base(settings_base, /row)
     label=widget_label(avgpolbase, value=' Average Pol ')
     state.avgpolselect=widget_droplist(avgpolbase, value=colortablenames, uvalue=colortableselections)


    settings_defaults=widget_button(settings_base, value = ' Defaults ', uvalue = 'polsettings_defaults', event_pro='gridview2_polsettings_done')
    settings_done = widget_button(settings_base, value = ' Done ', uvalue = 'polsettings_done', event_pro='gridview2_polsettings_done')

    widget_control, settings_base, /realize
    xmanager, 'gridview2_polsettings', settings_base, /no_block
    
endif


end

;----------------------------------------------------------------------
;Dummy procedure for dropdown lists
pro gridview2_polsettings_event, event

common gridview2_state
common gridstate

end

;-------------------------------------------------------------------
pro gridview2_polsettings_done, event

common gridview2_state
common gridstate

widget_control, event.id, get_uvalue = uvalue

case uvalue of
    'polsettings_done': begin
     
        polA=widget_info(state.polAselect, /DropList_Select)
        polB=widget_info(state.polBselect, /DropList_Select)
        avgpol=widget_info(state.avgpolselect, /DropList_Select)

        widget_control, state.polAselect, get_uvalue=colortableselections

        colortableselections=long(colortableselections)

        state.polAcolor=colortableselections[polA]
        state.polBcolor=colortableselections[polB]
        state.avgpolcolor=colortableselections[avgpol]

        if (state.currentpol eq 0) then state.colortable=state.polAcolor
        if (state.currentpol eq 1) then state.colortable=state.polBcolor
        if (state.currentpol eq 2) then state.colortable=state.avgpolcolor

     end

     'polsettings_defaults': begin

        state.polAcolor=3
        state.polBcolor=8
        state.avgpolcolor=1

        if (state.currentpol eq 0) then state.colortable=state.polAcolor
        if (state.currentpol eq 1) then state.colortable=state.polBcolor
        if (state.currentpol eq 2) then state.colortable=state.avgpolcolor

     end


    else:
 endcase

        widget_control, state.velslider, get_value=channel
        gridview2_display, channel
        gridview2_display_weights, channel

     widget_control, event.top, /destroy


end


;-------------------------------------------------------------------------
;Settings for colorbar scale
pro gridview2_scalesettings, event

   common gridview2_state
common gridstate
if (not (xregistered('gridview2_scalesettings', /noshow))) then begin




settingstitle = strcompress('Gridview Colorbar scale settings')

    settings_base =  widget_base(group_leader = state.baseID, $
                             /column, /base_align_right, title = settingstitle, $
                             uvalue = 'settings_base')
     minbase=widget_base(settings_base, /row)
     label=widget_label(minbase, value=' Colorbar minimum ')
     state.colorbarmin_text=widget_text(minbase, value=strcompress(state.colorbarmin), $
                     uvalue='colorbarminbox', /editable)

     maxbase=widget_base(settings_base, /row)
     label=widget_label(maxbase, value=' Colorbar maximum ')
     state.colorbarmax_text=widget_text(maxbase, value=strcompress(state.colorbarmax), $
                     uvalue='colorbarmaxbox', /editable)

    
    settings_defaults=widget_button(settings_base, value = ' Defaults ', uvalue = 'scalesettings_defaults', event_pro='gridview2_scalesettings_done')
    settings_done = widget_button(settings_base, value = ' Done ', uvalue = 'scalesettings_done', event_pro='gridview2_scalesettings_done')

    widget_control, settings_base, /realize
    xmanager, 'gridview2_scalesettings', settings_base, /no_block
    
endif


end

;----------------------------------------------------------------------
;Event handler for colobar scale settings window
pro gridview2_scalesettings_done, event

common gridview2_state
common gridstate
widget_control, event.id, get_uvalue = uvalue

case uvalue of
    'scalesettings_done': begin
     
        widget_control, state.colorbarmin_text, get_value=minstring
        state.colorbarmin=float(minstring[0])
        widget_control, state.colorbarmax_text, get_value=maxstring
        state.colorbarmax=float(maxstring[0])


     end

     'scalesettings_defaults': begin

       state.colorbarmin=-5.0
       state.colorbarmax=10.0  

     end


    else:
 endcase

        widget_control, state.velslider, get_value=channel
        gridview2_display, channel
        gridview2_display_weights, channel

     widget_control, event.top, /destroy


end


;---------------------------------------------------------------------
;RA dec conversion - takes ra in decimal hours, and dec in decimal degrees
;                    and converts to strings - uses adstring from GSFC
;                                              ASTRO LIB
      
pro radecconvert, rahr, decdeg, rastring, decstring
        
   radeg=(rahr/24.0)*360.0
      
   result=strtrim(adstring([radeg,decdeg], 2, /truncate))
      
;Use strmatch here?
;Get position of plus or minus sign for dec, and cut the string at
;that position
    
   signpos=strpos(result, '+')
   if (signpos eq -1) then signpos=strpos(result, '-')
           
   rastring=strmid(result, 0, signpos-1)
   decstring=strmid(result, signpos, strlen(result)-1)
  
end





;------------------------------------------------------------------
;Procedure that actually writes out the JPEG files
pro export_jpeg


common gridview2_state
common gridstate
widget_control, state.baseID, hourglass=1

widget_control, state.exportstart, get_value=startstring
widget_control, state.exportnumframes, get_value=numframesstring
widget_control, state.exportnumsteps, get_value=numstepsstring
widget_control, state.exportdirectory, get_value=exportdirectory

start=long(startstring[0])
numframes=long(numframesstring[0])
step=long(numstepsstring[0])

res=700

stop=start+(numframes-1)*step

device, decomposed=0

stretch, state.stretchlower, state.stretchupper

;Beginning for each file
for i=start, stop, step do begin

;GOTO Z device for output
thisDevice=!d.name
set_plot, 'Z', /COPY

Device, Set_Resolution=[res,res], Z_Buffer=0
Erase




lowerchan=i-step
upperchan=i+step
if (lowerchan lt 0) then lowerchan=0
if (upperchan gt n_elements(grid.velarr)-1) $
               then upperchan=n_elements(grid.velarr)-1

filename='vel_'+strcompress(long(grid.velarr[upperchan]), /remove_all)+'_'+$
         strcompress(long(grid.velarr[lowerchan]), /remove_all)+'.jpg'

outputfilename=exportdirectory+filename


loadct, state.colortable, /silent


ramin=state.ramin
ramax=state.ramax
decmin=state.decmin
decmax=state.decmax

hor, ramin, ramax
ver, decmin, decmax

posarray=[0.15,0.15,0.95,0.95]
xstyle=1
ystyle=1

charsize=1.0
xtitle='RA [HMS] J2000'
ytitle='Dec [DMS] J2000'
ticklen=-0.01
if (state.currentpol eq 0) then poltitle='Pol A'
if (state.currentpol eq 1) then poltitle='Pol B'
if (state.currentpol eq 2) then poltitle='Avg Pol'


plot, [0,0], /nodata, xstyle=xstyle, ystyle=ystyle, position=posarray, $
             charsize=charsize, xtitle=xtitle, ytitle=ytitle, $
             xtick_get=xvals, ytick_get=yvals, ticklen=ticklen



nxticklabels=n_elements(xvals)
nyticklabels=n_elements(yvals)



xvals=xvals
yvals=yvals

xspacing=((xvals[n_elements(xvals)-1]-xvals[0])*60.0)/(nxticklabels-1)

yspacing=((yvals[n_elements(yvals)-1]-yvals[0])*60.0)/(nyticklabels-1)

;ticlabels, xvals[0]*15.0, nxticklabels, xspacing, xticlabs,/ra,delta=1
;ticlabels, yvals[0], nyticklabels, yspacing, yticlabs, delta=1

xticlabs=ratickname(xvals*15.0)
yticlabs=dectickname(yvals)

erase

hor, ramax, ramin
ver, decmin, decmax

    
 loadct, state.colortable, /silent   
    
        xsize=561
        ysize=561
        px=105
        py=105

        map=reform(total(grid.d[lowerchan:upperchan,*,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax],1)/(float(2*step+1)))
        map=total(map,1)/2.
        map=reverse(congrid(reform(map), xsize,ysize))

;List of cuts made to get the subgrid structure
;d[lowerchan:upperchan,*,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax]
;w[lowerchan:upperchan,*,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax]
;cont[*,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax]
;cw[*,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax]
;;;;;;;;;;velarr[lowerchan:upperchan]
;;;;;;nx=state.xcubemax-state.xcubemin+1
;;;;;;ny=state.ycubemax-state.ycubemin+1
;;;;;;;;czmin=min(grid.velarr[lowerchan:upperchan])
;;;;;;;;;nz=n_elements(grid.velarr[lowerchan:upperchan])
;;;;ramin=(state.xcubemin*grid.deltara)/60.0+grid.ramin
;;;;;;decmin=(state.ycubemin*grid.deltadec)/60.0+grid.decmin

subgrid={name:grid.name, $
      RAmin:(state.xcubemin*grid.deltara)/3600.0+grid.ramin, $
      Decmin:(state.ycubemin*grid.deltadec)/60.0+grid.decmin, $
      Epoch:grid.epoch, $
      DeltaRA:grid.deltara, $
      DeltaDec:grid.deltadec, $
      NX:state.xcubemax-state.xcubemin+1, NY:state.ycubemax-state.ycubemin+1, $
      map_projection:grid.map_projection, $
      czmin:min(grid.velarr[lowerchan:upperchan]), $
      NZ:n_elements(grid.velarr[lowerchan:upperchan]), $
      velarr:grid.velarr[lowerchan:upperchan], $
      wf_type:grid.wf_type, $
      wf_fwhm:grid.wf_fwhm, $
      han:grid.han, $
      medsubtract:grid.medsubtract, $
      baseline:grid.baseline, $
      calib_facs:grid.calib_facs, $
      grms:grid.grms, $
      date:systime(0), $
      who:grid.who, $
      pos:grid.pos, $
      drift_list:grid.drift_list, $
      grid_makeup:grid.grid_makeup, $
      d:grid.d[lowerchan:upperchan,*,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax], $
      w:grid.w[lowerchan:upperchan,*,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax], $
      cont:grid.cont[*,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax], $
      cw:grid.cw[*,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax]} 

save, subgrid, filename='subgrid.sav'
delvarx, subgrid

;Fetch the smoothing state
widget_control, state.smoothslider, get_value=smoothval
smoothval=smoothval[0]

case state.currentscaling of

   'linearscaling': begin
       tvscl, smooth(map,smoothval, /edge_truncate), px, py
   end

    'logscaling': begin
       offset=min(map)-(max(map)-min(map))*0.01
       tvscl, smooth(alog10(map-offset), smoothval, /edge_truncate), px,py
    end

    'histeqscaling':begin
        tvscl, smooth(hist_equal(map, minv=min(map), maxv=max(map)), smoothval, /edge_truncate), px, py
    end

else:

endcase
      
       infostring=' Channels '+strcompress(lowerchan)+' to '+strcompress(upperchan)+ $
             '  Smoothed over -> '+strcompress(grid.velarr[lowerchan])+$
             ' to '+strcompress(grid.velarr[upperchan])+' km/s   '+poltitle
       
plot, [0,0], /nodata, xstyle=1, ystyle=1, /noerase, position=posarray, $ 
             xtickn=reverse(xticlabs), ytickn=yticlabs, $
             charsize=charsize, xtitle=xtitle, ytitle=ytitle, ticklen=ticklen
             
xyouts, 60, 30, infostring, /device, charsize=charsize

snapshot=tvrd()
tvlct, r,g,b, /get
device, Z_Buffer=1
set_plot, thisDevice
        
image = BytArr(3, res, res)
image[0,*,*] = r[snapshot]
image[1,*,*] = g[snapshot]
image[2,*,*] = b[snapshot]

write_jpeg, outputfilename, image, true=1, quality=100


endfor

widget_control, state.baseID, hourglass=0


end




;---------------------------------------------------------------------
;JPEG output - gui for procedure that outputs integrated velocity maps

pro gridview2_jpeg_output, event

common gridview2_state
common gridstate



if (not (xregistered('gridview2_jpeg_output', /noshow))) then begin

jpeg_outputtitle = strcompress('GRIDview2 JPEG Output')

    jpeg_output_base =  widget_base(group_leader = state.baseID, $
                             /row, /base_align_right, title = jpeg_outputtitle, $
                             uvalue = 'jpeg_output_base')
    configbase=widget_base(jpeg_output_base, /column, /align_left)
       startbase=widget_base(configbase, /row)
       widget_control, state.gotochannel, get_value=currentchannel
          state.exportstart=widget_text(startbase, xsize=8, value=currentchannel,$ 
                                      /editable, uvalue='exportstart')
           label=widget_label(startbase, value='Start Channel        ')
       numframesbase=widget_base(configbase, /row)
          state.exportnumframes=widget_text(numframesbase, xsize=8, value='1', $
                                      /editable, uvalue='exportnumframes')
          label=widget_label(numframesbase, value='# of frames')
       numstepsbase=widget_base(configbase, /row)
          state.exportnumsteps=widget_text(numstepsbase, xsize=8, value='5', $
                                      /editable, uvalue='exportnumsteps')
          label=widget_label(numstepsbase, value='# of chans')

          cd, current=pwd  ;Get current working directory into a string

    filesbase=widget_base(jpeg_output_base, /column, /align_left)
          directorybase=widget_base(filesbase, /row)
             state.exportdirectory=widget_text(directorybase, xsize=40, value=pwd+'/', $
                                               /editable, uvalue='exportdirectory')
             label=widget_label(directorybase, value='Output directory')

     buttonbase=widget_base(filesbase, xsize=100,ysize=100, /align_right, /column)

    jpeg_output_export = widget_button(buttonbase, value = ' Export ', $
                uvalue = 'jpeg_output_export', event_pro='gridview2_jpeg_output_event')
    cancel=widget_button(buttonbase, value=' Cancel ', uvalue='jpeg_output_cancel', $
                                       event_pro='gridview2_jpeg_output_event')

    widget_control, jpeg_output_base, /realize
    xmanager, 'gridview2_jpeg_output', jpeg_output_base, /no_block
    
endif


end

;----------------------------------------------------------------------
;Event handler for JPEG output gui
pro gridview2_jpeg_output_event, event

widget_control, event.id, get_uvalue = uvalue

case uvalue of
    'jpeg_output_export': begin

      export_jpeg

      widget_control, event.top, /destroy

  end

    'jpeg_output_cancel': widget_control, event.top, /destroy


    else:
endcase

end

;---------------------------------------------------------------------
;Postscript output
pro gridview2_ps_output_gui, event, xpixel, ypixel

common gridview2_state
common gridstate

;Result optical image option
state.ps_checkbox=0

if (not (xregistered('gridview2_ps_output_gui', /noshow))) then begin

ps_outputtitle = strcompress('GRIDview2 Postscript Output')

    ps_output_base =  widget_base(group_leader = state.baseID, $
                             /row, /base_align_right, title = ps_outputtitle, $
                             uvalue = 'ps_output_base')

     widget_control, state.gotochannel, get_value=currentchannel
     widget_control, state.velintegtext, get_value=channelrange
    configbase=widget_base(ps_output_base, /column, /align_left)
       pixelbase=widget_base(configbase, /row)
          label=widget_label(pixelbase,   value='Center pixel (X,Y)=')
          state.ps_centerpixel_x=widget_text(pixelbase, xsize=4, value=strcompress(xpixel, /remove_all), /editable, uvalue='ps_centerpixel_x')
          state.ps_centerpixel_y=widget_text(pixelbase, xsize=4, value=strcompress(ypixel, /remove_all), /editable, uvalue='ps_centerpixel_y')
          label=widget_label(pixelbase, value=' +/- ')
          state.ps_pixelrange=widget_text(pixelbase, xsize=4, value='15', /editable, uvalue='ps_pixelrange')
          label=widget_label(pixelbase, value=' pixels')

       channelbase=widget_base(configbase, /row)
          label=widget_label(channelbase, value='Center Channel=    ')
          state.ps_centerchannel=widget_text(channelbase, xsize=4, value=currentchannel, /editable, uvalue='ps_centerchannel')
          label=widget_label(channelbase, value='        +/- ')
          state.ps_channelrange=widget_text(channelbase, xsize=4, value=channelrange, /editable, uvalue='ps_channelrange')
          label=widget_label(channelbase, value=' channels')
          cd, current=pwd  ;Get current working directory into a string

       directorybase=widget_base(configbase, /row)
             label=widget_label(directorybase, value='File:')
             state.exportdirectory=widget_text(directorybase, xsize=50, value=pwd+'/gridview2.ps', $
                                               /editable, uvalue='exportdirectory')
       buttonbase=widget_base(configbase, xsize=340,ysize=25, /align_left, /row)

    ps_output_export = widget_button(buttonbase, value = ' Export ', $
                uvalue = 'ps_output_export', event_pro='gridview2_ps_output_event')
    cancel=widget_button(buttonbase, value=' Cancel ', uvalue='ps_output_cancel', $
                                       event_pro='gridview2_ps_output_event')
    label=widget_label(buttonbase, value='    ')
             checkboxbase=widget_base(buttonbase, ysize=25, /row, /align_left, /nonexclusive)
               checkbox=widget_button(checkboxbase, value='Check this for a DSS image', $
                                                    uvalue='ps_checkbox', event_pro='gridview2_ps_output_event')
    

    widget_control, ps_output_base, /realize
    xmanager, 'gridview2_ps_output_gui', ps_output_base, /no_block
    
endif


end


;----------------------------------------------------------------------
;Event handler for JPEG output gui
pro gridview2_ps_output_event, event

common gridview2_state
common gridstate

widget_control, event.id, get_uvalue = uvalue

case uvalue of
    'ps_output_export': begin

      export_postscript

      widget_control, event.top, /destroy

    

  end

    'ps_checkbox': begin

        if (state.ps_checkbox eq 0) then begin
            state.ps_checkbox=1
        endif else begin
            state.ps_checkbox=0
        endelse

    end

    'ps_output_cancel': widget_control, event.top, /destroy

    


    else:
endcase


 state.keystatus=0

end

;-------------------------------------------------------------------
;Export publication quality postscript files
pro export_postscript

common gridview2_state
common gridstate



rah=grid.ramin+(dindgen(n_elements(grid.d[0,0,*,0]))+0.5)*grid.deltara/3600.   ;units of hours
dec=grid.decmin+(dindgen(n_elements(grid.d[0,0,0,*]))+0.5)*grid.deltadec/(60.) ;units of degrees

  widget_control, state.ps_centerpixel_x, get_value=xpixstring
  widget_control, state.ps_centerpixel_y, get_value=ypixstring
  widget_control, state.ps_pixelrange, get_value=pixelrangestring
  widget_control, state.ps_centerchannel, get_value=centerchannelstring
  widget_control, state.ps_channelrange, get_value=channelrangestring
  widget_control, state.exportdirectory, get_value=outputfile


;Parameters used to make the map
  xcenter=long(xpixstring)
    xcenter=xcenter[0]
  ycenter=long(ypixstring)
    ycenter=ycenter[0]
  pixelrange=long(pixelrangestring)
    pixelrange=pixelrange[0]
  xmin=xcenter-pixelrange
  ymin=ycenter-pixelrange
  xmax=xcenter+pixelrange
  ymax=ycenter+pixelrange
  
  centerchannel=long(centerchannelstring)
    centerchannel=centerchannel[0]

  channelrange=long(channelrangestring)
    channelrange=channelrangestring[0]
  chanmin=centerchannel-channelrange
  chanmax=centerchannel+channelrange

;check for bad input
  if (chanmin lt channelrange) then begin
      print, 'first if'
      chanmin=0
      chanmax=2*pixelrange
  endif

  if (chanmax gt n_elements(grid.velarr)-1) then begin
     chanmax=n_elements(velarr)-1
     chammin=chanmax-2*pixelrange
  endif

  if (xmin lt 0) then begin
    xmin=0
    xmax=2*pixelrange
  endif

  if (xmax gt n_elements(grid.d[0,0,*,0])-1) then begin
      xmax=n_elements(grid.d[0,0,*,0])-1
      xmin=xmax-2*pixelrange
  endif

  if (ymin lt 0) then begin
    ymin=0
    ymax=2*pixelrange
  endif

  if (ymax gt n_elements(grid.d[0,0,0,*])-1) then begin 
      ymax=n_elements(grid.d[0,0,0,*])-1
      ymin=ymax-(2*pixelrange)    
  endif

; xmin=xcenter-pixelrange
; ymin=ycenter-pixelrange
;Redefine xcenter and ycenter after possible shift
xcenter=xmin+pixelrange
ycenter=ymin+pixelrange


;Make the map
map=reform((grid.d[*,0,xmin:xmax,ymin:ymax]+grid.d[*,0,xmin:xmax,ymin:ymax])/2.0)

map=total(reform(map[chanmin:chanmax,*,*]),1)/(2*channelrange+1)

;Get the ra and dec for the axes
rahr=rah[xmin:xmax]
decdeg=dec[ymin:ymax]
velocity=grid.velarr[chanmin:chanmax]

decnr=dec[ycenter]  ;center of image
cosdec=cos(decnr*!dpi/180.)
deltaram=0.25*grid.deltara*cosdec   ;needed

ramin=min(rahr)-0.5*(deltaram/60.0)
ramax=max(rahr)+0.5*(deltaram/60.0)
decmin=min(decdeg)-0.5/60.0
decmax=max(decdeg)+0.5/60.0

rahr_opt=(ramax+ramin)/2.0
decdeg_opt=(decmax+decmin)/2.0

;If chosen by the user, then download a DSS image
if (state.ps_checkbox eq 1) then begin

  widget_control, state.baseID, hourglass=1
  print, 'Please wait.  Obtaining optical image from service...' 
  ;querydss, [rah[xcenter]*15.0, dec[ycenter]], optimage, header, imsize=(decmax-decmin)*60.0, survey='2b'
  querydss, [rahr_opt*15.0, decdeg_opt], optimage, header, imsize=(decmax-decmin)*60.0, survey='2b'
  widget_control, state.baseID, hourglass=0

endif


rms=stddev(map)
c_levels=rms*(dindgen(20)+1)
n_levels=rms*[-2.0, -1.0]

;c_levels=[2,3,4,5,6,8,10,12,14,16,18,20,22,24,25,30,40,50]

psopen, outputfile, set_font='Times', /tt_font, xsize=7.0, ysize=7.0, /inches, /isolatin1

charsize=1.0
position=[0.2,0.2,0.95,0.95]
title=strcompress(grid.velarr[chanmax], /remove_all)+' to '+strcompress(grid.velarr[chanmin], /remove_all)+' km/s'

yvalues=dindgen(241)/6.0
ylabels=strarr(241)
xvalues=dindgen(1441)/60.0
xlabels=strarr(1441)

;Create nice TeX style sexigesimal labels for ra
for i=0, n_elements(xvalues)-1 do begin

xlabels[i]='+'+strcompress(floor(xvalues[i]), /remove_all)+'^h'+strcompress(round((xvalues[i]-floor(xvalues[i]))*60.0), /remove_all)+$
         '^{m}00^{s}'

endfor

;Create nice TeX style sexigesimal labels for dec
for i=0, n_elements(yvalues)-1 do begin

ylabels[i]='+'+strcompress(floor(yvalues[i]), /remove_all)+'^o'+strcompress(round((yvalues[i]-floor(yvalues[i]))*60.0), /remove_all)+$
         '^{\prime}00^{\prime\prime}'

endfor

xindex=where(xvalues gt ramin and xvalues lt ramax)
yindex=where(yvalues gt decmin and yvalues lt decmax)


contour, map, rahr,decdeg, xstyle=1, ystyle=1, $
         xrange=[rahr[n_elements(rahr)-1],rahr[0]],yrange=[decdeg[0], decdeg[n_elements(decdeg)-1]], charsize=charsize, $
         xtitle='!6 Right Ascension (J2000)', ytitle='Declination (J2000)', $
         xticks=n_elements(xindex)-1, yticks=n_elements(yindex)-1, $
         xtickv=xvalues[xindex], ytickv=yvalues[yindex], $
         xtickn=textoidl(xlabels[xindex]), ytickn=textoidl(ylabels[yindex]), $
         ticklen=-0.015,  title=title, $
         position=position,levels=c_levels, xminor=6, yminor=4

;overplot optical image if selected

if (state.ps_checkbox eq 1) then begin
  TVSCL, -optimage, !X.WINDOW(0), !Y.WINDOW(0), $
    XSIZE = !X.WINDOW(1) - !X.WINDOW(0), $
    YSIZE = !Y.WINDOW(1) - !Y.WINDOW(0), /NORM 
 endif


;replot contours over optical image
contour, map, rahr,decdeg, xstyle=1, ystyle=1, $
         xrange=[rahr[n_elements(rahr)-1],rahr[0]],yrange=[decdeg[0], decdeg[n_elements(decdeg)-1]], charsize=charsize, $
         xtitle='!6 Right Ascension (J2000)', ytitle='Declination (J2000)', $
         xticks=n_elements(xindex)-1, yticks=n_elements(yindex)-1, $
         xtickv=xvalues[xindex], ytickv=yvalues[yindex], $
         xtickn=textoidl(xlabels[xindex]), ytickn=textoidl(ylabels[yindex]), $
         ticklen=-0.015,  title=title, $
         position=position,levels=c_levels, xminor=6, yminor=4, /noerase


;Negative contours
contour, map, rahr,decdeg, xstyle=1, ystyle=1, $
         xrange=[rahr[n_elements(rahr)-1],rahr[0]],yrange=[decdeg[0], decdeg[n_elements(decdeg)-1]], charsize=charsize, $
         xtitle='!6 Right Ascension (J2000)', ytitle='Declination (J2000)', $
         xticks=n_elements(xindex)-1, yticks=n_elements(yindex)-1, $
         xtickv=xvalues[xindex], ytickv=yvalues[yindex], $
         xtickn=textoidl(xlabels[xindex]), ytickn=textoidl(ylabels[yindex]), $
         ticklen=-0.015,  title=title, $
         position=position,levels=n_levels, /noerase, c_linestyle=2, xminor=6, yminor=4

xyouts, 7000, 1500,'Map rms: '+strcompress(rms, /remove_all)+' mJy/beam', /device
xyouts, 7000, 1000,'Positive Contours: 1 through 20 times rms', /device
xyouts, 7000, 500,'Negative Contours: -2 and -1 times rms', /device


psclose
print, 'Postscript saved to '+outputfile

end



;----------------------------------------------------------------------
pro gridview2_fetch_spectrum, xpos, ypos

common gridview2_state
common gridstate

state.spectrumon=1

if (xpos lt 0) then xpos=0
if (ypos lt 0) then ypos=0

if (xpos gt state.sx-1) then xpos=state.sx-1
if (ypos gt state.sy-1) then ypos=state.sy-1

;xsizepixel=state.sx
;ysizepixel=state.sy

;ramin=state.ramin
;ramax=state.ramax
;decmin=state.decmin
;decmax=state.decmax

;xsizearcmin=15.0*(ramax-ramin)*60.0
;ysizearcmin=(decmax-decmin)*60.0

;pixscalex=xsizepixel/xsizearcmin
;pixscaley=ysizepixel/ysizearcmin

;widget_control, state.apertureFWHM, get_value=FWHMstring
;FWHM=double(FWHMstring[0])

widget_control, state.gotochannel, get_value=channelstring
channel=long(channelstring[0])

;FOLLOWING TWO LINES NO LONGER NEEDED BECAUSE OF ZOOM FUNCTION
;widget_control, state.spectrumwindowwidth, get_value=widthstring
;width=long(widthstring[0])

;FWHM=4.0  ; units of arcmin

;sigma=FWHM*pixscalex/2.3548    ;sigma in pixels
;print, sigma

;nx=21
;ny=21

;X = DINDGEN(nx) # REPLICATE(1.0, ny)
;X=X-10
;Y = REPLICATE(1.0, nx) # DINDGEN(ny)
;Y=Y-10

;U=X^2+Y^2

;A=1/(2*sigma^2)

;B=A*EXP(-U/(2.0*sigma^2))

width=250

startchannel=channel-width/2
stopchannel=channel+width/2

if (startchannel lt 0) then startchannel=0
if (stopchannel gt n_elements(grid.velarr)-1) then stopchannel=n_elements(grid.velarr)-1

spectrum=dblarr(n_elements(grid.velarr))
spectrumweights=dblarr(n_elements(grid.velarr))

widget_control, state.baseID, hourglass=1

;print, xpos, ypos

;Convert from screen coordinates back to current cube coordinates   
;xcube=-round((float(n_elements(grid.d[0,0,*,0]))/state.sx)*xpos)+n_elements(grid.d[0,0,*,0])-1
;ycube=round((float(n_elements(grid.d[0,0,0,*]))/state.sy)*ypos)

xcube=-round((float(state.xcubemax-state.xcubemin+1)/state.sx)*xpos)+(state.xcubemax-1)+1
ycube=round((float(state.ycubemax-state.ycubemin+1)/state.sy)*ypos)+(state.ycubemin)



;print, xcube, ycube

state.current_xcube=xcube
state.current_ycube=ycube

spectrum=(grid.d[*,0,xcube,ycube]+grid.d[*,1,xcube,ycube])/2.0
spectrumweights=(grid.w[*,0,xcube,ycube]+grid.w[*,1,xcube,ycube])/2.0

;for i=startchannel,stopchannel do begin

;    print, i

;   map=reform(grid.d[i,*,state.xcubemin:state.xcubemax, state.ycubemin:state.ycubemax])
;   map=total(map,1)/2.0
   
;   image=reverse(congrid(map, state.sx,state.sy))
  
   ;Take care of clicks near the edge
;   lowerx=long(xpos-15)
;   upperx=long(xpos+15)
;   lowery=long(ypos-15)
;   uppery=long(ypos+15)


;   if (lowerx lt 0) then lowerx=0
;   if (upperx gt state.sx-1) then upperx=long(state.sx-1)
;   if (lowery lt 0) then lowery=0
;   if (uppery gt state.sy-1) then uppery=long(state.sy-1)

;   imagecutout=image[lowerx:upperx, lowery:uppery]   ;mJy/beam
  
;   convol_result=convol(imagecutout, B, total(B))    ;mJy

   ;weightedmean=total(B*imagecutout)
;   fluxdensity=total(convol_result)/total(B)

;   spectrum[i]=fluxdensity/10.0

;endfor

;BRIAN'S ATTEMPT AT CALIRBATION!!   :)
;which isn't really working

;map=reform(grid.cont[*,state.xcubemin:state.xcubemax,state.ycubemin:state.ycubemax ])
;map=total(map,1)/2.
;xpos=212
;ypos=266
;image=reverse(congrid(map, state.sx,state.sy))

;lowerx=long(xpos-15)
;   upperx=long(xpos+15)
;   lowery=long(ypos-15)
;   uppery=long(ypos+15)

;imagecutout=image[lowerx:upperx, lowery:uppery]   ;mJy/beam

;convol_result=convol(imagecutout, B, total(B)) 

;fluxdensity=total(convol_result)/total(B)

;print, fluxdensity

;FOR CUBE 8 ONLY - ONLY A TEST
;calfactor=fluxdensity/(10^.56)

;spectrum=spectrum/calfactor


widget_control, state.baseID, hourglass=0






;Plot velocity on x-axis
if (state.spectrumaxisstatus eq 0) then begin

widget_control, state.plotwindowspectrum, get_value=index
wset, index

hor, grid.velarr[stopchannel], grid.velarr[startchannel]
ver,0.0,1.1

color='0000FF'XL   ;RED
plot,grid.velarr[startchannel:stopchannel],spectrumweights[startchannel:stopchannel]/state.maxweight,xstyle=1, ystyle=1, xtitle='Velocity [km/s]', charsize=0.7, ytitle='Weight', /nodata, position=[0.07,0.1,0.98,0.3]

oplot,grid.velarr[startchannel:stopchannel], spectrumweights[startchannel:stopchannel]/state.maxweight, color=color

flag, grid.velarr[channel], color='00FF00'XL  ; GREEN current channel (velocity)



hor, grid.velarr[stopchannel], grid.velarr[startchannel]
ver, min(spectrum[startchannel:stopchannel]), max(spectrum[startchannel:stopchannel])

device, decomposed=1
color='0000FF'XL   ;RED
plot,grid.velarr[startchannel:stopchannel],spectrum[startchannel:stopchannel],xstyle=1, ystyle=1, charsize=0.7, ytitle='Flux Density [mJy/beam]', /nodata, /noerase, position=[0.07,0.3,0.98,0.98], $
xtickn=[' ',' ',' ',' ',' ',' ',' ']
oplot,grid.velarr[startchannel:stopchannel], spectrum[startchannel:stopchannel], color=color
;x=indgen(n_elements(spectrum))
;xcentroid=total(spectrum*x)/total(spectrum)

flag, grid.velarr[channel], color='00FF00'XL  ; GREEN current channel (velocity)







device, decomposed=0


endif else begin   ;Plot channel on x-axis

widget_control, state.plotwindowspectrum, get_value=index
wset, index

hor, startchannel, stopchannel
ver, 0.0,1.1

color='0000FF'XL   ;RED
plot,spectrumweights/state.maxweight,xstyle=1, ystyle=1, xtitle='Channel', charsize=0.7, ytitle='Weight', /nodata,  position=[0.07,0.1,0.98,0.3]
oplot, spectrumweights/state.maxweight, color=color
;x=indgen(n_elements(spectrum))
;flag, total(spectrum*x)/total(spectrum), color='00FF00'XL  ; GREEN current channel
flag, channel, color='00FF00'XL  ; GREEN centroid for redshift


hor, startchannel, stopchannel
ver, min(spectrum[startchannel:stopchannel]), max(spectrum[startchannel:stopchannel])

device, decomposed=1
color='0000FF'XL   ;RED
plot,spectrum,xstyle=1, ystyle=1, charsize=0.7, ytitle='Flux Density [mJy]', /nodata, $
position=[0.07,0.3,0.98,0.98], /noerase, $
xtickn=[' ',' ',' ',' ',' ',' ',' ']
oplot, spectrum, color=color
;x=indgen(n_elements(spectrum))
;flag, total(spectrum*x)/total(spectrum), color='00FF00'XL  ; GREEN current channel
flag, channel, color='00FF00'XL  ; GREEN centroid for redshift






device, decomposed=0



endelse

state.spectrum=spectrum
state.spectrumweights=spectrumweights
state.spectrum_xmin=startchannel
state.spectrum_xmax=stopchannel
state.spectrum_ymin=min(spectrum[startchannel:stopchannel])
state.spectrum_ymax=max(spectrum[startchannel:stopchannel])

hor
ver

gridview2_display, channel
gridview2_display_weights, channel

state.getspectrumstatus=0

end






;-----------------------------------------------------------------
;Displays spectrum in spectrum window box

pro gridview2_plotspectrum


common gridview2_state
common gridstate


widget_control, state.gotochannel, get_value=channelstring
channel=long(channelstring[0])

spectrum=state.spectrum
spectrumweights=state.spectrumweights

;Plot velocity on x-axis
if (state.spectrumaxisstatus eq 0) then begin


;widget_control, state.plotwindowspectrumweights, get_value=index
;wset, index

widget_control, state.plotwindowspectrum, get_value=index
wset, index

device, decomposed=1

hor, grid.velarr[state.spectrum_xmax], grid.velarr[state.spectrum_xmin]
ver, 0.0,1.1

color='0000FF'XL   ;RED
plot,grid.velarr[state.spectrum_xmin:state.spectrum_xmax],spectrumweights[state.spectrum_xmin:state.spectrum_xmax]/state.maxweight,$
           xstyle=1, ystyle=1, xtitle='Velocity [km/s]', $
           charsize=0.7, ytitle='Weight', /nodata, position=[0.07,0.1,0.98,0.3]
oplot,grid.velarr[state.spectrum_xmin:state.spectrum_xmax], spectrumweights[state.spectrum_xmin:state.spectrum_xmax]/state.maxweight, color=color

flag, grid.velarr[channel], color='00FF00'XL  ; GREEN current channel (velocity)



hor, grid.velarr[state.spectrum_xmax], grid.velarr[state.spectrum_xmin]
ver, state.spectrum_ymin, state.spectrum_ymax



color='0000FF'XL   ;RED
plot,grid.velarr[state.spectrum_xmin:state.spectrum_xmax],spectrum[state.spectrum_xmin:state.spectrum_xmax],$
xstyle=1, $
 ystyle=1, charsize=0.7, ytitle='Flux Density [mJy/beam]', /nodata, /noerase, position=[0.07,0.3, 0.98,0.98], $
xtickn=[' ',' ',' ',' ',' ',' ',' ',' ']
oplot,grid.velarr[state.spectrum_xmin:state.spectrum_xmax],spectrum[state.spectrum_xmin:state.spectrum_xmax], color=color
;x=indgen(n_elements(spectrum))
;xcentroid=total(spectrum*x)/total(spectrum)

flag, grid.velarr[channel], color='00FF00'XL  ; GREEN current channel (velocity)


device, decomposed=0


endif else begin   ;Plot channel on x-axis

widget_control, state.plotwindowspectrum, get_value=index
wset, index

device, decomposed=1

hor, state.spectrum_xmin, state.spectrum_xmax
ver, 0.0,1.1

color='0000FF'XL   ;RED
plot,spectrumweights/state.maxweight,xstyle=1, ystyle=1, $
    xtitle='Channel', charsize=0.7, ytitle='Weight', /nodata, position=[0.07,0.1,0.98,0.3]
oplot, spectrumweights/state.maxweight, color=color
flag, channel, color='00FF00'XL  ; GREEN centroid for redshift





hor, state.spectrum_xmin, state.spectrum_xmax
ver, state.spectrum_ymin, state.spectrum_ymax

color='0000FF'XL   ;RED
plot,spectrum,xstyle=1, ystyle=1, charsize=0.7, ytitle='Flux Density [mJy/beam]', /nodata,  /noerase, $
position=[0.07,0.3, 0.98,0.98], xtickn=[' ',' ',' ',' ',' ',' ',' ',' ']
oplot, spectrum, color=color
;x=indgen(n_elements(spectrum))
;flag, total(spectrum*x)/total(spectrum), color='00FF00'XL  ; GREEN current channel
flag, channel, color='00FF00'XL  ; GREEN centroid for redshift


device, decomposed=0

endelse    

end









;------------------------------------------------------------------
;Update the box in the spectrum window - NOT USED AS OF OCT 18.
pro gridview2_update_spectrum_box

common gridview2_state
common gridstate
widget_control, state.plotwindowspectrum, get_value=index
wset, index

hor, 0,n_elements(grid.velarr)-1
ver, 0,1.0

plot, [0,0], /nodata, xstyle=1, ystyle=1, xtitle='Channel', charsize=0.7

widget_control, state.gotochannel, get_value=channelstring

channel=long(channelstring[0])
widget_control, state.spectrumwindowwidth, get_value=widthstring

width=long(widthstring[0])

device, decomposed=1
color='0000FF'XL   ;RED
tvboxbk, [width, 1.0], channel, 0.5, /data, color=color
device, decomposed=0




hor
ver

end

;Added January 2006 - view both pols

pro plot_both_pols, spectrum, spectrumweights, velarr, chanmin, chanmax, maxweight, group_leader

;help, spectrum

  if (not (xregistered('plot_both_pols', /noshow))) then begin

       chooseinterpbase =  widget_base(group_leader = group_leader, /column, /base_align_right, $
                     title = 'Spectrum - both polarizations')       

        bothpolsdraw=widget_draw(chooseinterpbase, xsize=800, ysize=400, uvalue='plot_both_pols_draw', /button_events)

      bottombase=widget_base(chooseinterpbase, /row, /base_align_left)
       ;  scalebase=widget_base(bottombase, /column, /base_align_left)
;	 xscalebase=widget_base(scalebase, column=4, /align_left, /grid_layout)
;         label=widget_label(xscalebase, value=' Channels min:')
;           xmin=widget_text(xscalebase, /editable, xsize=8)
;         label=widget_label(xscalebase, value='     Channels max:') 
;           xmax=widget_text(xscalebase, /editable, xsize=8)
         
;         yscalebase=widget_base(scalebase, column=4, /align_left, /grid_layout)
;         label=widget_label(yscalebase, value='Intensity min:')
;           ymin=widget_text(yscalebase, /editable, xsize=8)
;         label=widget_label(yscalebase, value='    Intensity max:')
 ;          ymax=widget_text(yscalebase, /editable, xsize=8)
            

         
         buttonbase=widget_base(bottombase, column=5, /base_align_right, /grid_layout)
         ;chooseinterp_zoomin=widget_button(bottombase, value= 'Zoom In', uvalue='chooseinter_zoomin')
         ;chooseinterp_zoomout=widget_button(bottombase, value= 'Zoom Out', uvalue='chooseinterp_zoomout')
         ;chooseinterp_reset=widget_button(bottombase, value= 'Reset', uvalue='chooseinterp_reset')
         ;chooseinterp_default=widget_button(bottombase, value= 'Default', uvalue='chooseinterp_default')
         chooseinterp_done = widget_button(bottombase, value = ' Close Window ', uvalue = 'plot_both_pols_done')

       ;Realization
       widget_control, chooseinterpbase, /realize
       xmanager, 'plot_both_pols', chooseinterpbase, /no_block
   endif

;make the plot of both pols

widget_control, bothpolsdraw, get_value=index

wset, index

hor, velarr[chanmax], velarr[chanmin]
ver, 0.0,1.2

;Plot weights
plot,velarr[chanmin:chanmax],spectrumweights[chanmin:chanmax,0]/maxweight,$
    xstyle=1, ystyle=1, $
    xtitle='cz [km/s]', charsize=1.0, ytitle='Weight', /nodata, position=[0.07,0.1,0.98,0.3]

device, decomposed=1

oplot,velarr[chanmin:chanmax],spectrumweights[chanmin:chanmax,0]/maxweight, color='0000FF'XL  ;A is red
oplot,velarr[chanmin:chanmax],spectrumweights[chanmin:chanmax,1]/maxweight, color='00FF00'XL  ;B is green

;Plot spectrum

hor, velarr[chanmax], velarr[chanmin]
ver, -10, max(spectrum[chanmin:chanmax])*1.2

plot, velarr[chanmin:chanmax], spectrum[chanmin:chanmax,0], $
      xstyle=1, ystyle=1, charsize=1.0, $
      ytitle='Flux Density [mJy/beam]', /nodata, /noerase,  position=[0.07,0.3,0.98,0.98], $
      xtickn=[' ',' ',' ',' ',' ',' ',' ']

oplot, velarr[chanmin:chanmax], spectrum[chanmin:chanmax,0], color='0000FF'XL
oplot, velarr[chanmin:chanmax], spectrum[chanmin:chanmax,1], color='00FF00'XL

oplot, velarr, fltarr(n_elements(velarr)), linestyle=2, color='FFFFFF'XL

xyouts, 50, 6, 'Pol A', color='0000FF'XL, /device, charsize=2.0
xyouts, 150, 6, 'Pol B', color='00FF00'XL, /device, charsize=2.0

device, decomposed=0

end

;------------------------------------------------
;Events for plots of both polarizations
pro plot_both_pols_event, event

widget_control, event.id, get_uvalue = uvalue

case uvalue of
    'plot_both_pols_done': begin
        widget_control, event.top, /destroy
        
    end

else:

endcase

end

;------------------------------------------------------------------
;------------------------------------------------------------------
;HELP QUICKSTART BLOCK

;------------------------------------------------------------
;Help
pro gridview2_quickstart, event

   common gridview2_state
common gridstate
h=['GRIDview2TM Quickstart',$
'A sanctioned product of LOVEDATA, Inc. Ithaca, NY',$
'CONTENTS',$
'',$
'    * Overview',$
'    * Menu Options',$
'    * Main Window',$
'    * Colorbar',$
'    * Weights Window',$
'    * Information Display',$
'    * Controls',$
'    * Spectrum Window',$
'    * Imaging Window',$
'',$
'Overview',$
'GRIDview2 is a set of procedures written in IDL to assist in the ',$
'visualization of 3D data cubes, primarily those of the ALFALFA ',$
'extragalactic survey. The data from meridian transit drifts scans ',$
'are grouped and combined via procedures (posfind and grid_prep) into ',$
'regularly gridded data, with approprimate headers, coordinates, ',$
'and all pertinent information concerning the creation of the cube. ',$
'GRIDView acts as a tool that can read this data and display/average ',$
'cube spectral and continuum planes, extract single pixel spectra, ',$
'show the weights function derived in the creation of the cube, access ',$
'images, and export for use in other media.',$
'',$
'',$
'Menu Options',$
'File 	',$
'',$
'    * Export JPEG - allows user to export JPEGs in ',$
'        single channel or averaged channel maps in ',$
'        the directory of their chosing. Files are named by velocity range.',$
'    * Export Single Postscript - create a publication quality postscript of a map portion.', $
'    * Exit - exit the program!',$
'',$
'Settings 	',$
'',$
'    * Pol Colors - Choose color table to use with polarziation options',$
'    * Colorbar Scale - If a constant, fixed scaling option is used, then ',$
'        the program will use the intensity range specified here.',$
'',$
'Scaling 	',$
'',$
'    * Linear - Fixed minimum and maximum scaling for auto or constant scaling.',$
'    * Log - Scaling where offset is calculated and subtracted before taking a ',$
'               base 10 logarithm.',$
'    * Histogram EQ - Histogram equalization of the image.',$
'',$
'Colorbar 	',$
'',$
'    * Autoscale - minimum and maximum values are determined using the full ',$
'                     range of the image',$
'    * Constant Scale - minimum and maximum values are taken from input ',$
'                          in the settings menu (default is -5 to +10).',$
'',$
'Catalogs 	',$
'',$
'    * Multi-AGC mode - NEW!  For a given channel, this mode will display', $
'                       all AGC galaxies in the range.  RED for those that do not', $
'                       have any redshift information, LIGHT BLUE for those that have', $
'                       been searched or have a known optical or 21cm derived velocity,', $
'                       but are not in the current channel range (determined by the', $
'                       velocity width), and YELLOW for those with known redshifts', $
'                       or that have been search within the currently selected channel', $
'                       range.', $
'    * AGC (known cz) - known cz catalog overlay.  If the galaxy has 21 cm line',$
'                         information, OR has been searched before for HI, ', $
'                         it will be displayed as a yellow box.  Hover your mouse', $
'                         over the box to see the DETCODE at the bottom of the ',$
'                         GRIDview2 window', $
'    * AGC (no known cz) - unknown cz catalog overlay',$ 
'    * NVSS - continuum mode option. Sources above a peak intensity value ',$ 
'                of 5 mJy/beam are display and color coded.',$ 
'    * 2D extractor - not current active.  Ask Brian Kent for more information',$ 
'    * 3D extractor - pass a keyword at startup using:  ',$
'                                 gridview2, grid, cat3D=filename_string ',$
'                     Sources will be displayed as a red cross with the source number', $
'',$
'Help 	',$
'',$
'    * Quickstart - Text version of this guide.',$
'    * About - Information, update information, and modification history.',$
'',$
'Keyboard Shortcuts ',$
'', $
'  g - Get a DSS/Sloan/NVSS image',$ 
'  h - View component drifts of the grid',$
'  p - Postscript file output',$
'  n - NED database query', $
'',$
'',$
'Main Window',$
'The main window display shows a color image of the current spectral channel ',$
'   (in single channel mode) or an averaged channel map (if ',$
'   the velocity "integration") option is selected. ',$
'    If continuum mode is selected, the continuum for the cube will be shown. ',$
'    Any selected catalogs will be overlaid on the map. The user may LEFT click',$ 
'    with the mouse and drag a square box to zoom to a selected region. ',$
'    A single LEFT mouse click will restore the display to the full range of ',$
'    the cube in RA and Dec. RIGHT clicking will display a DSS image ',$
'    (spectral mode) or NVSS image (continuum mode), overlaying a box showing ',$
'    the image size.  The user may hover the mouse over a point and press ',$
'    the "n" key to query NED within 10 arcminutes of that point on the sky.', $
'',$
'Colorbar',$
'The colorbar window shows the current dynamic range for the currently selected',$ 
'colortable (available in the Settings menu). ',$
'Colorbar can autoscale for the currently displayed map, or use a fixed ',$
'constant scale that can be set in the settings menu. The user can LEFT ',$
'click and drag the colorbar to refine the dynamic range.',$
'',$
'Weights Window',$
'The weights window shows the same area currently seen in the main window. ',$
'It shows a linearly scaled image of the weights for the current ',$
'(single or averaged) spectral channel or continuum map.',$
'',$
'Information Display',$
'The left column of the information display show information from the ',$
'current grid structure - The grid name, size, velocity range and ',$
'coordinate epoch. The switch for spectral and continuum mode is also ',$
'located is this display. The right column of the information display ',$
'shows current coordinates for the mouse pointer in the main window - ',$
'RA/DEC, velocity of the current channel, and x/y/intensnity of the',$
'current cube pixel.',$
'',$
'Controls',$
'In spectral mode, the user can access different spectral channels via ',$
'the channel slider, the NEXT and PREVIOUS buttons, or by typing a desired ',$
'channel into the channel box. If the slider is active, the left and right',$
'arrow keys may be used to scroll through the data cube. A second slider ',$
'can be used to boxcar smooth the currently displayed map, with a SCREEN ',$
'pixel width determined from the slider value. A simple boxcar smooth can ',$
'also be used in the spectral direction by clicking the velocity integration',$  
'option. A average will be taken over +/- the number of channels the user ',$ 
'specifies in the box (default is 5 channels). The box below the velocity ',$
'integration is the Arecibo General Catalog Display. If the AGC is turned ',$
'"on" (catalogs menu), then hovering the mouse over an AGC galaxy will ',$
'display the number, description, and optical and 21cm measured cz values.',$
'',$
'Spectrum Window',$
'Clicking the GET SPECTRUM button, followed by a click on the main window. ',$
'A spectrum will be shown (x axis units are chosen by user - velocity ',$
'or channels). A single LEFT click in the spectrum window will zoom ',$
'out to the full spectrum range. A LEFT click and dragging a box will ',$
'zoom into the selected box range.',$
'',$
'Imaging Window',$
'The user can input an image size; left clicking on the main window will ',$
'retreive an image of that size and display it in the imaging window. ',$
'DSS 2 Blue images are retreived in spectral mode, and NVSS image are ',$
'retrieved in continuum mode. ',$
'',$
'Brian Kent, Cornell University.']


if (not (xregistered('gridview2_quickstart', /noshow))) then begin

helptitle = strcompress('Gridview Quickstart')

    help_base =  widget_base(group_leader = state.baseID, $
                             /column, /base_align_right, title = helptitle, $
                             uvalue = 'help_base')

    help_text = widget_text(help_base, /scroll, value = h, xsize = 90, ysize = 50)
    
    help_done = widget_button(help_base, value = ' Done ', uvalue = 'quickstart_done')

    widget_control, help_base, /realize
    xmanager, 'gridview2_quickstart', help_base, /no_block
    
endif


end

;----------------------------------------------------------------------

pro gridview2_quickstart_event, event

 common gridview2_state
common gridstate
widget_control, event.id, get_uvalue = uvalue

case uvalue of
    'quickstart_done': widget_control, event.top, /destroy
    else:
endcase

end




;----------------------------------------------------------------------
pro gridview2_nedgui, rahr, decdeg

 common gridview2_state
common gridstate
widget_control, state.baseID, hourglass=1

nedquery, rahr*15.0, decdeg,10.0, numberinfo=numberinfo, string_array=string_array

h=[numberinfo, '  Name 			         RA(J2000)   DEC(J2000)	TYPE      VEL    Z	   SKY DIST(arcmin)',string_array]


if (not (xregistered('gridview2_nedgui', /noshow))) then begin

helptitle = strcompress('Gridview NED Results')

    help_base =  widget_base(group_leader = state.baseID, $
                             /column, /base_align_right, title = helptitle, $
                             uvalue = 'ned_base')

    help_text = widget_text(help_base, /scroll, value = h, xsize = 110, ysize = n_elements(h)+2)
    
    help_done = widget_button(help_base, value = ' Done ', uvalue = 'nedgui_done')

    widget_control, help_base, /realize
    xmanager, 'gridview2_nedgui', help_base, /no_block
    
endif



widget_control, state.baseID, hourglass=0



end



;----------------------------------------------------------------------

pro gridview2_nedgui_event, event

 common gridview2_state
common gridstate
widget_control, event.id, get_uvalue = uvalue

case uvalue of
    'nedgui_done': begin
        state.keystatus=0
        widget_control, event.top, /destroy
    end


    else:
endcase

end

;----------------------------------------------------------------------
;Rearrange structure for 3D extraction catalogs
;Not used as of Dec. 14

pro convert3Dcat, catin, catout

table=catin

catout=   {id:table.(0)[0], $
           ra:table.(1)[0], $
           dec:table.(2)[0], $
           radec:table.(3)[0], $
           channel:table.(4)[0], $
           rapixel:table.(5)[0], $
           decpixel:table.(6)[0], $
           width:table.(7)[0], $   
           ara:table.(8)[0], $     
           adec:table.(9)[0], $    
           peak_flux:table.(10)[0], $
           int_flux:table.(11)[0], $ 
           sn:table.(12)[0], $
           rms:table.(13)[0], $
           cz:table.(14)[0], $ 
           agc:table.(15)[0], $
           comment:table.(16)[0]}


for j=0,n_elements(table.(0))-1 do begin


tabletemp={id:table.(0)[j], $ 
           ra:table.(1)[j], $
           dec:table.(2)[j], $
           radec:table.(3)[j], $   
           channel:table.(4)[j], $ 
           rapixel:table.(5)[j], $ 
           decpixel:table.(6)[j], $
           width:table.(7)[j], $
           ara:table.(8)[j], $
           adec:table.(9)[j], $
           peak_flux:table.(10)[j], $
           int_flux:table.(11)[j], $
           sn:table.(12)[j], $ 
           rms:table.(13)[j], $  
           cz:table.(14)[j], $
           agc:table.(15)[j], $
           comment:table.(16)[j]}

catout=[catout, tabletemp]

endfor


catout=catout[1:n_elements(catout)-1]




end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro gridview2_fetchdrifts, xcube, ycube

common gridview2_state
common gridstate


 widget_control, state.baseID, hourglass=1

;Get indicies of contributing elements
index=where(grid.grid_makeup.i eq xcube AND grid.grid_makeup.j eq ycube)

;Get contributing file names
driftindex=where(grid.grid_makeup[index].driftname ne '')
state.driftfiles[0:n_elements(driftindex)-1]=grid.grid_makeup[index].driftname[driftindex]

state.makeupindex=index

url='http://caborojo.astro.cornell.edu/alfalfalog/posmaster.list'

spawn, 'wget -q -O posmaster.list ' + "'" + url + "'"

filename='posmaster.list'

openr,lun,filename,/get_lun
        maxfiles=10000L
        templist=strarr(maxfiles)
        ifound=0L
        on_ioerror,done
       
        while (1) do begin
          inpline=''
          readf,lun,inpline
          strar=strsplit(inpline,' ',/extract)
          templist[ifound]=strtrim(strar[0],2)
          ifound=ifound+1
        end
        done:
        close, lun
        free_lun,lun

   ;print, ''
   ;print, templist[0]   ;Informs user of last update
   ;print, ''
        
   posfilelist=templist[1:ifound-1]  ; trim list to number of files found
   
   ;store for later use
   state.posfilelist[0:n_elements(posfilelist)-1]=posfilelist
   
driftfilesindex=where(state.driftfiles ne '')
driftfiles=state.driftfiles[driftfilesindex]


for i=0,n_elements(driftfiles)-1 do begin
       decrange=strmid(driftfiles[i], 7,7)
      
    for j=0,n_elements(posfilelist)-1 do begin

        result=strpos(posfilelist[j], decrange)
        if (result[0] ne -1) then begin
            ;print, driftfiles[i], '   ', posfilelist[j]
            
            restore, posfilelist[j]
                                ;Have to check the individual posfiles
                                ;to make sure the file is associated
                                ;with the date
            filecheck=where(pos.name eq driftfiles[i])

            if (filecheck[0] ne -1) then begin
            directoryindex=strpos(posfilelist[j], 'pos')
            
            state.filelist[i]=strmid(posfilelist[j],0,directoryindex)+driftfiles[i]
            
            ;print, '  ',filelist[i]
            endif

        endif      

    endfor
endfor


;Start GUI


if (not (xregistered('gridview2_fetchdrifts', /noshow))) then begin

drifttitle = strcompress('Gridview Drift Viewer - Selected Pixel History')

    drift_base =  widget_base(group_leader = state.baseID, $
                             /row, /base_align_left, /align_top, title = drifttitle, $
                             uvalue = 'drift_base')

  control_base=widget_base(drift_base, /column, /base_align_left)
  filelistwidget=widget_list(control_base, value=state.driftfiles[0:n_elements(driftindex)-1], xsize=36, ysize=10, uvalue='filelist', event_pro='gridview2_fetchdrifts_event')
    
    state.drift_draw=widget_draw(drift_base, xsize=600, ysize=570, frame=1, $
               uvalue='drift_draw', event_pro='gridview2_fetchdrifts_event', retain=2)

    label=widget_label(control_base, value=' Selected pixel is (X,Y) = ('+strcompress(xcube, /remove_all)+','+strcompress(ycube, /remove_all)+')')
    view_drift = widget_button(control_base, value = ' View  ', uvalue = 'view_drift', event_pro='gridview2_fetchdrifts_event')
    drift_done = widget_button(control_base, value = ' Close ', uvalue = 'drift_done', event_pro='gridview2_fetchdrifts_event')

    label=widget_label(control_base, value='Choose VIEW to open the selected drift')
    label= widget_label(control_base, value='with FLAGBB.PRO')
    label= widget_label(control_base, value=' ')
    label= widget_label(control_base, value='NOTE: ')
    label= widget_label(control_base, value='Stopping FLAGBB mid-program will not')
    label= widget_label(control_base, value=' remove the drift files from memory,')
    label= widget_label(control_base, value=' and may cause font problems and')
    label= widget_label(control_base, value=' differences in GRIDVIEW.  Choose')
    label= widget_label(control_base, value=' a PIXEL SPECTRUM to reset.')

    widget_control, drift_base, /realize
    xmanager, 'gridview2_fetchdrifts', drift_base, /no_block
    
endif

gridview2_fetchdrifts_plot, -1

 widget_control, state.baseID, hourglass=0

end


;-----------------------------------------
pro gridview2_fetchdrifts_event, event

common gridview2_state
common gridstate

widget_control, event.id, get_uvalue = uvalue

case uvalue of
    'drift_done': begin
        state.keystatus=0
        widget_control, event.top, /destroy
        state.driftfiles[*]=''
        state.filelist[*]=''
        state.posfilelist[*]=''
        listnumber=-1

        widget_control, state.gotochannel, get_value=channelstring

        channel=long(channelstring)
        channel=channel[0]
   
        gridview2_display, channel
        gridview2_display_weights, channel

    end

    'filelist': begin
        listnumber=event.index
        state.listnumber=listnumber
        gridview2_fetchdrifts_plot,listnumber 

    end

    'view_drift': begin
        listnumber=state.listnumber
        ;Trim the list for plotting
        driftfilesindex=where(state.driftfiles ne '')
        driftfiles=state.driftfiles[driftfilesindex]
        filelist=state.filelist[driftfilesindex]

        ;full path and filename to open
        filename=filelist[listnumber]

        ;get correct pos file
        decrange=strmid(driftfiles[listnumber], 7,7)
        posfilelist=state.posfilelist
        posfilename=''

         for j=0,n_elements(posfilelist)-1 do begin

             result=strpos(posfilelist[j], decrange)
             if (result[0] ne -1) then begin

                restore, posfilelist[j]
                                ;Have to check the individual posfiles
                                ;to make sure the file is associated
                                ;with the date
                filecheck=where(pos.name eq driftfiles[listnumber])

                if (filecheck[0] ne -1) then begin
                    posfilename=posfilelist[j]
                endif
             endif      

         endfor

        

         if (posfilename ne '' AND listnumber ne -1) then begin
             
             
             print, 'Opening pos file '+posfilename
             restore, posfilename
             print, 'Opening drift '+filename
             restore, filename

             ;start the flagging session
             flagbb, dred, cont_pt, mask, pos, gausav=11, han=3, agc=1

         endif
        
    end

    else:
endcase


end

;-------------------------------------------
pro gridview2_fetchdrifts_plot, listnumber

common gridview2_state
common gridstate

ramin=min(grid.grid_makeup.ra)
ramax=max(grid.grid_makeup.ra)
decmin=min(grid.grid_makeup.dec)
decmax=max(grid.grid_makeup.dec)

widget_control, state.drift_draw, get_value=window

wset, window

device, decomposed=1

plot, dist(5), /nodata, xrange=[ramax+0.5/15.0, ramin-0.5/15.0], yrange=[decmin-0.5, decmax+0.5], $
      xstyle=1, ystyle=1, xtitle='RA [Hrs]', ytitle='Dec [Deg]', $
        title='Map center at RA: '+$
        strcompress((ramin+ramax)/2.0)+$
        ', Dec: '+strcompress((decmin+decmax)/2.0), $
        charsize=1.4

tvboxbk, [ramax-ramin, decmax-decmin], $
       ramin+(ramax-ramin)/2.0, decmin+(decmax-decmin)/2.0, $ 
       /data, color='0000FF'XL

;Trim the list for plotting
driftfilesindex=where(state.driftfiles ne '')
driftfiles=state.driftfiles[driftfilesindex]


for i=0, n_elements(driftfiles)-1 do begin

    posindex=where(grid.pos.name eq driftfiles[i])

    if (posindex[0] ne -1) then begin

        oplot, grid.pos[posindex].rahr, grid.pos[posindex].decdeg, psym=3, color='00FF00'XL

    endif 

endfor

;Plot currently selected drift

if (listnumber ne -1) then begin
posindex=where(grid.pos.name eq driftfiles[listnumber])
oplot, grid.pos[posindex].rahr, grid.pos[posindex].decdeg, psym=3, color='7F7FFF'XL
endif


index=state.makeupindex

plotsym, 0, /fill
plots, grid.grid_makeup[index].ra, grid.grid_makeup[index].dec, psym=8, color='0000FF'XL, /data

charsize=1.0
xyouts, 10, 40,'- Selected drift', /device, color='7F7FFF'XL, charsize=charsize
xyouts, 10, 25,'- Contributing drifts', /device, color='00FF00'XL, charsize=charsize
xyouts, 10, 10,'- Grid boundary', /device, color='0000FF'XL, charsize=charsize

device, decomposed=0


end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Main program block;;;;;;;;;;;;
pro gridview2, data, cat3D=cat3D


 common gridstate

grid=data

    if (n_params() lt 1) then begin
         print, ''
         print, 'USAGE:  gridview2, grid, cat3D=cat3D'
         print, '           grid  = 3D data cube output from grid_prep.pro'
         print, '           cat3D = optional keyword string of 3D extraction catalog'
         print, ''
         return
    endif

if (n_elements(cat3D) ne 0) then begin

    check=findfile(cat3D, count=count)

    if (count eq 0) then begin
        print, ''
        print, 'The 3D extracted catalog specified does not exist.'
        print, ''
        return
    endif

endif


  gridview2_startup, cat3D=cat3D 

  common gridview2_state 
 
  gridview2_display, 0 
  gridview2_display_weights, 0 

end 
