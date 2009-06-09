; version = 18072008
pro cleanview, llx=llx, lly=lly, urx=urx, ury=ury

; Common blocks we need from gridview2
common gridview2_state
common gridstate

; Intialize common blocks
common cleanview_state,		viewstate, viewdisp
common clean,      		cleanstate, dbox2, cbox2, cleand, c_cleand

; New plan for this part.  I have problems with the way the boundary conditions are handled and
; am running into flux troubles in the outer 12' of the cube where it is unclear whether a peak 
; is at the actual peak of the beam or at some other part.  To get around this, I want to add 
; 12' to each side of the selected region, if possible, and clean that.  12' to each side pads 
; the selected region by 6', which is the HP width of the beam at ~5%.  This should ensure 
; that what is cleaned doesn't run into boundary issues and is usuable.  The map buffer size will 
; be stored in the variable map_buffer and in the cleanstate structure.
map_buffer = 6
orig_llx = llx & orig_lly = lly & orig_urx = urx & orig_ury = ury
llx -= map_buffer
urx += map_buffer
lly -= map_buffer
ury += map_buffer
if llx LT 0 OR urx GT (grid.nx-1) then begin
	ok = Error_Message(['Expansion of the region selection box in the', $
			    'RA direction encountered a grid edge. Consider', $
			    'using the adjacent cube for this source.'])
endif
if lly LT 0 OR ury GT (grid.ny-1) then begin
	ok = Error_Message(['Expansion of the region selection box in the', $
			    'Dec direction encountered a grid edge. ', $
			    'Consider using the adjacent cube for this ', $
			    'source.'])
endif

; Define, grid, dbox2, cbox2, cleand, and c_cleand
rbox = grid.d[*, *, llx:urx, lly:ury]
wbox = grid.w[*, *, llx:urx, lly:ury]
dbox2 = total( rbox*wbox, 2, /Double) / total( wbox, 2, /Double)
  bad = where( finite(dbox2) NE 1 )
  if bad[0] NE -1 then dbox2[bad] = 0.0
;+ Define the variables used to deal with the map buffer
lower_x = map_buffer
lower_y = map_buffer
upper_x = (urx-llx)-map_buffer
upper_y = (ury-lly)-map_buffer
cleand = (dbox2 * 0.0)[*, lower_x:upper_x, lower_y:upper_y]

rbox = grid.cont[*, llx:urx, lly:ury]
wbox = grid.cw[*, llx:urx, lly:ury]
cbox2 = total( rbox*wbox, 1, /Double) / total( wbox, 1, /Double)
  bad = where( finite(cbox2) NE 1 )
  if bad[0] NE -1 then cbox2[bad] = 0.0
c_cleand = (cbox2 * 0.0)[lower_x:upper_x, lower_y:upper_y]

;print,size(dbox2),size(cleand)

; Define the viewstate structure (cleanview main window)
viewstate = { baseID: 0L, display_win: 0L, colorbar_win: 0L, slider: 0L, summary: 0L, clog: 0L, $
		DispSpec: 0L, DispCont: 0L, ImgSmthSlider: 0L, BxcSmthValue: 0L, BxcSmthOnOff: 0L, $
		linear: 0L, logari: 0L, histeq: 0L, acscale: 0L, ccscale: 0L, $
		agccz: 0L, agcncz: 0L, cathvc: 0L, cat3d: 0L, catchk: 0L, catcnt: 0L, $
		colorBLU: 0L, colorGRN: 0L, colorRED: 0L, $
		xval: 0L, yval: 0L, vval: 0L, ival: 0L, $
		xpix: 0L, ypix: 0L, vpix: 0L, $
		ramin: (grid.ramin+(llx+map_buffer)*grid.deltara/3600.0), $
		ramax: (grid.ramin+(urx-map_buffer)*grid.deltara/3600.0), $
		decmin: (grid.decmin+(lly+map_buffer)*grid.deltadec/60.0), $
		decmax: (grid.decmin+(ury-map_buffer)*grid.deltadec/60.0), $
		fmin: -5.0, fmax: 10.0, $
		drawfluxbox: 0, cornerx: 0.0, cornery: 0.0, px: [0.0, 0.0], py: [0.0, 0.0], sx: 0.0, sy: 0.0, mousestate: 0, $
		clean_rms: 0L, clean_flux: 0L, clean_sigma: 0L, clean_niter: 0L, clean_chn: 0L, clean_sum: 0L, $
		save_dir: 0L, save_name: 0L, $
		dir: '', frmsta: 0L, frmnum: 0L, frmste: 0L, $
		xs: dblarr(2), ys: dblarr(2), xregion: dblarr(2), yregion: dblarr(2)}

; Define viewdisp structure (controls how the maps are displayed in the cleanview window)
viewdisp = { UseCons: 0, UseAuto: 0, $
		DispType: 1, ImgSmth: 3, BxcSmth: 0, $
		UseLin: 0, UseLog: 0, UseHeq: 0, $
		CBMin: -5.0, CBMax: 10.0, colortab: 1, $
		UseKA: 0, UseUA: 0, UseHV: 0, Use3D: 0, UseCk: 0, UseCt: 0 }
; Setup the variables based off the gridview2_state
viewdisp.UseKA = state.agcoverlaystatus
viewdisp.UseUA = state.agcoverlaystatus_nocz
viewdisp.UseHV = state.dhbboverlaystatus
viewdisp.Use3D = state.overlaystatus3D
viewdisp.UseCk = state.checkliststatus
viewdisp.UseCT = state.cntsrcstatus

; Define cleanstate structur that holds of the info about cleaning
;+ This structure is initialzed after the display window 

; Begin the Gui
main_base = widget_base(Title='CLEANview', /Column, mbar=top_menu)
viewstate.baseID = main_base
; Menu items
filemenu=widget_button(top_menu, value=' File')
buttonjpeg=widget_button(filemenu, value=' Export JPEG ', event_pro='cleanview_jpeg')
buttonsave=widget_button(filemenu, value=' Save Data ', uvalue='dsave')
buttonexit=widget_button(filemenu, value=' Exit ', uvalue='exit', /Separator)
colormenu=widget_button(top_menu, value=' Color ')
viewstate.colorBLU = widget_button(colormenu, value=' Blue/White ',  uvalue='blue', /Checked_Menu)
viewstate.colorGRN = widget_button(colormenu, value=' Green/White ', uvalue='green', /Checked_Menu)
viewstate.colorRED = widget_button(colormenu, value=' Red/White ',   uvalue='red', /Checked_Menu)
settingmenu=widget_button(top_menu, value=' Settings ')
viewstate.DispSpec=widget_button(settingmenu, value=' Show Spectral Map   ', uvalue='disps', /Checked_Menu)
viewstate.DispCont=widget_button(settingmenu, value=' Show Continuum Map  ', uvalue='dispc', /Checked_Menu)
junk =             widget_button(settingmenu, value=' Set Image Smoothing ', uvalue='smth', /Separator)
scalingmenu=widget_button(top_menu, value=' Scaling ')
viewstate.linear=widget_button(scalingmenu, value=' Linear ', uvalue='linscl', /Checked_Menu) 
viewstate.logari=widget_button(scalingmenu, value=' Logarithmic ', uvalue='logscl', /Checked_Menu) 
viewstate.histeq=widget_button(scalingmenu, value=' Histogram EQ ', uvalue='heqscl', /Checked_Menu) 
colorbarmenu=widget_button(top_menu, value=' Colorbar ')
viewstate.acscale=widget_button(colorbarmenu, value=' Autoscale ', uvalue='colorbar_ascale', /Checked_Menu)
viewstate.ccscale=widget_button(colorbarmenu, value=' Constant scale ', uvalue='colorbar_cscale', /Checked_Menu)
catalogmenu=widget_button(top_menu, value=' Catalogs ')
viewstate.agccz=widget_button(catalogmenu, value=' AGC Catalog (known cz)', uvalue='agc_cz', /Checked_Menu)
viewstate.agcncz=widget_button(catalogmenu, value=' AGC Catalog (no known cz)', uvalue='agc_ncz', /Checked_Menu) 
viewstate.cathvc=widget_button(catalogmenu, value= ' HVC catalog (de Heij 2002)', uvalue='cathvc',/Checked_Menu)
viewstate.cat3d=widget_button(catalogmenu, value=' 3D Catalog ', uvalue='cat3d', /Checked_Menu)
viewstate.catchk=widget_button(catalogmenu, value=' Checklist (GALflux) ', uvalue='catchk', /Checked_Menu)
viewstate.catcnt=widget_button(catalogmenu, value=' Continuum (>250mJy/bm) ', uvalue='catcnt', /Checked_Menu)
helpmenu=widget_button(top_menu, value=' Help ')
buttonstart=widget_button(helpmenu, value=' Quick Start ', uvalue='qstart')
buttonabout=widget_button(helpmenu, value=' About ', uvalue='about',/Separator)

; Plot Base
upper_base = widget_base(main_base, /Row)
	viewstate.display_win = widget_draw(upper_base, uvalue='draw', XSize=512, YSize=512, $
		Frame=1, Event_Pro='cleanview_display_event', /Motion_Events, /Button_Events, $
		Retain=2, /Keyboard_Events)
	viewstate.colorbar_win = widget_draw(upper_base, uvalue='color', XSize=100, YSize=512, $
		Frame=1, Event_Pro='cleanview_display_event', /Motion_Events, /Button_Events)
; Controls and info base
middle_base = widget_base(main_base, /Row)
	; Info on data window
	info_base = widget_base(middle_base, /Column, /Frame)
		info_sub1 = widget_base(info_base, /Align_Center)
			label = widget_label(info_sub1, value='Image Info:')
		info_sub2 = widget_base(info_base, /Row, /Align_Left)
			info_sub3 = widget_base(info_sub2, /Column, /Align_Left)
			info_sub4 = widget_base(info_sub2, /Column, /Align_Left)
				label = widget_label(info_sub3, value='Pixel X:  ')
				viewstate.xpix = widget_label(info_sub3, value='---')
				label = widget_label(info_sub4, value='RA:  ')
				viewstate.xval = widget_label(info_sub4, value='00 00 00.0')
				label = widget_label(info_sub3, value='Pixel Y:  ')
				viewstate.ypix = widget_label(info_sub3, value='---')
				label = widget_label(info_sub4, value='DEC: ')
				viewstate.yval = widget_label(info_sub4, value='+00 00 00.0')
				label = widget_label(info_sub3, value='Channel:  ')
				viewstate.vpix = widget_label(info_sub3, value='   0')
				label = widget_label(info_sub4, value='cz:  ')
				viewstate.vval = widget_label(info_sub4, $
					value=string(grid.velarr[0],Format='(F7.1)'))
				label = widget_label(info_sub3, value='Intensity: ')
				viewstate.ival = widget_label(info_sub4, value='999.999999')
		info_sub5 = widget_base(info_base, /Row, /Align_Center)
			viewstate.slider = widget_slider(info_sub5, min=0, max=n_elements(grid.velarr)-1, $
				title='Channel', value=0, uvalue='slider', /drag)
	; Control of CLEAN
	clean_base = widget_base(middle_base, /Column, /Frame)
		clean_sub1 = widget_base(clean_base, /Align_Center)
			label = widget_label(clean_sub1, value='CLEAN Control:')
		clean_sub2 = widget_base(clean_base, /Align_Left, Row=6)
			label = widget_label(clean_sub2, value='Grid RMS:         ')
			viewstate.clean_rms = widget_text(clean_sub2, value=string(0.0, Format='(F5.2)'), $
				XSize=7)
			label = widget_label(clean_sub2, value='Flux Limit:       ')
			viewstate.clean_flux = widget_text(clean_sub2, value='0.0', $
				uvalue='limitf', /Editable, XSize=7)
			label = widget_label(clean_sub2, value='Flux Sigma Limit: ')
			viewstate.clean_sigma = widget_text(clean_sub2, value='1.0', $
				uvalue='limits', /Editable, XSize=7)
			label = widget_label(clean_sub2, value='Iteration Limit:  ')
			viewstate.clean_niter = widget_text(clean_sub2, value='5000', $
				uvalue='limiti', /Editable, XSize=7)
			label = widget_label(clean_sub2, value='Channel Range: ')
			viewstate.clean_chn = widget_text(clean_sub2, value='0, 1023', $
				uvalue='limitc', /Editable, XSize=10)
			label = widget_label(clean_sub2, value='Combine to Improve S/N? ')
			viewstate.clean_sum = widget_button(clean_sub2, value='No ', uvalue='cleansum')
			
	; Output of CLEAN
	log_base = widget_base(middle_base, /Column, /Frame)
		label = widget_label(log_base, value='CLEAN Log:', /Align_Center)
		viewstate.summary = widget_label(log_base, value='Flux: ----     Iteration: ----', $
			/Align_Left)
		viewstate.clog = widget_text(log_base, value='----', /Scroll, YSize=12, XSize=30)
; Button base
bottom_base = widget_base(main_base, /Row)
	beam_button = widget_button(bottom_base, value='View Beam', uvalue='beam')
	clean_button = widget_button(bottom_base, value='Clean Region', uvalue='clean')
	vfeld_button = widget_button(bottom_base, value='Velocity Field', uvalue='vfeld')
	flux_button = widget_button(bottom_base, value='Measure Flux', uvalue='flux')

; Go for it
widget_control, main_base, /realize
xmanager, 'cleanview', main_base, /no_block

; Setup default settings - settings, scaling, and colorbar
viewdisp.DispType = 1
widget_control, viewstate.DispSpec, set_button=1
viewdisp.UseCons = 1
viewdisp.UseAuto = 0
widget_control, viewstate.linear, set_button=1
viewdisp.UseLin = 1
viewdisp.UseLog = 0
viewdisp.UseHeq = 0
widget_control, viewstate.ccscale, set_button=1
viewdisp.colortab = 1
widget_control, viewstate.colorBLU, Set_Button=1

; Setup default settings - catalog overlays
if viewdisp.UseKA EQ 1 then widget_control, viewstate.agccz, Set_Button=1
if viewdisp.UseUA EQ 1 then widget_control, viewstate.agcncz, Set_Button=1
if viewdisp.UseHV EQ 1 then widget_control, viewstate.cathvc, Set_Button=1
if viewdisp.Use3D EQ 1 then widget_control, viewstate.cat3d, Set_Button=1
if viewdisp.UseCk EQ 1 then widget_control, viewstate.catchk, Set_Button=1
if viewdisp.UseCT EQ 1 then widget_control, viewstate.catcnt, Set_Button=1

; Setup default settings - velocity slider
widget_control, state.velslider, get_value=channel
channel = long(channel)
widget_control, viewstate.slider, set_value=string(channel)

;+ Fake initialization of cleanstate here (enough to get us through)
cleanstate = {  llx: (llx+map_buffer), lly: (lly+map_buffer), urx: (urx-map_buffer), ury: (ury-map_buffer), map_buffer: map_buffer, $
	beam: 1, beamlog: 1, $
	rms: 0.0, flux: 0.0, sigma: 1.0, niter: 5000L, crange: [0L,1023L], allowsum: 0, $
	restore_fwhm: 1, restore_pa: 1, area_corr: 1, $
	cleanlog: strarr(5150), cleanexit: intarr(1024), cleansum: intarr(2)}

cleanview_updatemap, channel, /Reset

; Setup default settings - grid rms
widget_control, viewstate.display_win, get_value=index
wset, index
mask = bytarr(125, 50)
tv, mask, 110, 340
plots, [110, 110, 235, 235, 110], [340, 390, 390, 340, 340], Color=FSC_Color('Yellow'), /Device
xyouts, 120, 360,'Initializing....', /Device, Charsize=1.5, Color=FSC_Color('Yellow')
widget_control, /Hourglass
; Start the beam modeling here
t_center = [(urx+llx)/2, (ury+lly)/2]
t_width = [urx-llx+1, ury-lly+1]
t_mopt = [25, 31, 35, 41, 45, 51]
t_maps = 2.0*(max([(urx-llx+1), (ury-lly+1)]) - 2.0*map_buffer)
if min( where(t_maps LE t_mopt) ) EQ -1 then begin
	t_maps = max(t_mopt)
endif else begin
	t_maps = t_mopt[ min( where(t_maps LE t_mopt) ) ]
endelse
build_beam_grid2, grid, llx, lly, urx, ury, beams, MapSize=t_maps, FluxCorr=flux_factor, FWHM=beam_fwhm, PA=beam_pa, /Silent, Output=beamlog, /CExts, /Smart
beam_fwhm = round(beam_fwhm*1000.0)/1000.0
beam_pa = round(beam_pa*1000.0)/1000.0

; Initialize cleanstate here
cleanstate = {  llx: (llx+map_buffer), lly: (lly+map_buffer), urx: (urx-map_buffer), ury: (ury-map_buffer), map_buffer: map_buffer, $
	beam: beams, beamlog: beamlog, $
	rms: 0.0, flux: 0.0, sigma: 1.0, niter: 5000L, crange: [0L,1023L], allowsum: 0, $
	restore_fwhm: beam_fwhm, restore_pa: beam_pa, area_corr: flux_factor, $
	cleanlog: strarr(5150), cleanexit: intarr(1024), cleansum: intarr(2)}
;+ Cheap shot method to get the sky RMS by not using
;+ all the points.
n_grids = n_elements(grid.d)
values = floor( randomu(seed, n_grids/16L)*n_grids )
values = values[ sort(values) ]
values = values[ uniq(values) ]
grid_rms = robust_sigma(grid.d[values])

widget_control, viewstate.clean_rms, set_value=string(grid_rms, Format='(F5.2)')
cleanstate.rms = grid_rms

cleanview_updatemap, channel

end


pro cleanview_display_event, event

common cleanview_state
common gridstate
common clean

widget_control, event.id, get_uvalue=uvalue

case uvalue of
	'draw': begin
		Event_Types = ['DOWN', 'UP', 'MOTION', '?', '?', '?', '?']
		This_Event = Event_Types[event.type]
		case This_Event of
			'DOWN'  : begin
				if event.press eq 1 then begin
					viewstate.mousestate=1
					; Turn motion events on for the draw widget.
					widget_control, viewstate.display_win, Draw_Motion_Events=1
					widget_control, viewstate.display_win, get_value=index
										
					; Create a pixmap. Store its ID. Copy window contents into it.
					Window, 19, /Pixmap, XSize=512, YSize=512
					Device, Copy=[0, 0, 512, 512, 0, 0, index]
					
					; Get and store the static corner of the box.
					viewstate.cornerx = event.x
					viewstate.cornery = event.y
				endif
			end
			'UP'    : begin
				if viewstate.mousestate then begin
					viewstate.mousestate=0
					
					; Erase the last box drawn. Destroy the pixmap.
					widget_control, viewstate.display_win, get_value=index
					wSet, index
					Device, Copy=[0, 0, 512, 512, 0, 0, 19]
					
					; Order the box coordinates.
					sx = Min([viewstate.cornerx, event.x], Max=dx)
					sy = Min([viewstate.cornery, event.y], Max=dy)
					
					WDelete, 19
					
					;Determine coordinates FOR SQUARE BOX to zoom in on
					width=abs(dx-sx)
					xpos_min=sx
					xpos_max=dx
					ypos_min=sy
					
					if viewstate.drawfluxbox NE 0 then begin
						ypos_max=dy

						xpos_min = round( xpos_min - viewstate.px[0]-1 )
						xpos_max = round( xpos_max - viewstate.px[0]-1 )
						ypos_min = round( ypos_min - viewstate.py[0]-1 )
						ypos_max = round( ypos_max - viewstate.py[0]-1 )
						
						xmax = -round(float(cleanstate.urx-cleanstate.llx+1) / viewstate.sx*xpos_min) + $
							(cleanstate.urx)
						xmin = -round(float(cleanstate.urx-cleanstate.llx+1) / viewstate.sx*xpos_max) + $
							(cleanstate.urx)
						ymin =  round(float(cleanstate.ury-cleanstate.lly+1) / viewstate.sy*ypos_min) + $
							(cleanstate.lly)
						ymax =  round(float(cleanstate.ury-cleanstate.lly+1) / viewstate.sy*ypos_max) + $
							(cleanstate.lly)
						
						;Protect against going 'outside' the cube
						if (xmin lt cleanstate.llx) then xmin = cleanstate.llx
						if (ymin lt cleanstate.lly) then ymin = cleanstate.lly
						if (xmax gt cleanstate.urx) then xmax = cleanstate.urx
						if (ymax gt cleanstate.ury) then ymax = cleanstate.ury
						
						if viewstate.drawfluxbox EQ 1 then begin
							loadct, viewdisp.colortab, /Silent
						
							;Call the flux measureing routine GALFLUX
							viewstate.drawfluxbox = 0
						
							widget_control, viewstate.display_win, Clear_Events=1
							cleanflux, llx=xmin, lly=ymin, urx=xmax, ury=ymax
						endif
						if viewstate.drawfluxbox EQ 2 then begin
							loadct, viewdisp.colortab, /Silent
						
							;Call the flux measureing routine cleanvfield
							viewstate.drawfluxbox = 0
						
							widget_control, viewstate.display_win, Clear_Events=1
							xmin = xmin + cleanstate.llx
							xmax = xmax + cleanstate.llx
							ymin = ymin + cleanstate.lly
							ymax = ymax + cleanstate.lly
							fieldview, llx=xmin, lly=ymin, urx=xmax, ury=ymax
						endif
					endif
				endif
				widget_control, viewstate.slider, get_value=channel
				channel = long(channel[0])

				cleanview_updatemap, channel
			end
			'MOTION': begin
				if event.press eq 0 AND viewstate.mousestate eq 1 then begin
					; Here is where the actual box is drawn and erased.
					; First, erase the last box.

					widget_control, viewstate.display_win, get_value=index
					wset, index
					Device, Copy=[0, 0, 512, 512, 0, 0, 19]

					; Get the coodinates of the new box and draw it.
					sx = viewstate.cornerx
					sy = viewstate.cornery
					dx = event.x
					dy = event.y
					loadct, viewdisp.colortab, /silent
					width=abs(dx-sx)

					;plot a flux measurement box of arbitrary size
					if viewstate.drawfluxbox EQ 1 then begin
						color=FSC_Color('Red')
					
						PlotS, [sx, sx, dx, dx, sx], [sy, dy, dy, sy, sy], /Device, $
							Color=color, thick=1.5

						xyouts, sx+10, sy-10, 'FLUX BOX', /device, color=color
					endif
					if viewstate.drawfluxbox EQ 2 then begin
						color=FSC_Color('Magenta')
					
						PlotS, [sx, sx, dx, dx, sx], [sy, dy, dy, sy, sy], /Device, $
							Color=color, thick=1.5

						xyouts, sx+10, sy-10, 'V. FIELD BOX', /device, color=color
					endif
				endif
			end
			else    :
		endcase

		widget_control, viewstate.display_win, get_value=index
		wset, index
	
		hor, viewstate.ramax, viewstate.ramin
		ver, viewstate.decmin, viewstate.decmax
		
		xdevice=event.x
		ydevice=event.y
	
		widget_control, viewstate.slider, get_value=chnstring
		channel=long(chnstring[0])
		
		result=convert_coord(xdevice, ydevice, /Device, /Double, /To_Data)
		xdata=result[0]
		ydata=result[1]
	
		if (xdata lt viewstate.ramin OR xdata gt viewstate.ramax OR $
			ydata lt viewstate.decmin OR ydata gt viewstate.decmax) then begin
			widget_control, viewstate.xval, set_value='----'
			widget_control, viewstate.yval, set_value='----'
			widget_control, viewstate.xpix, set_value='----'
			widget_control, viewstate.ypix, set_value='----'
			widget_control, viewstate.ival, set_value='----'
		endif else begin
			radecconvert, xdata, ydata, rastring, decstring
			xpix = round((xdata - viewstate.ramin)*3600.0 / grid.deltara + cleanstate.llx)
			ypix = round((ydata - viewstate.decmin)*60.0 / grid.deltadec + cleanstate.lly)
			if viewdisp.DispType then begin
				ival = cleand[channel, xpix-cleanstate.llx, ypix-cleanstate.lly]
			endif else begin
				ival = c_cleand[xpix-cleanstate.llx, ypix-cleanstate.lly]
			endelse
			widget_control, viewstate.xval, set_value=string(rastring[0], format='(a11)')
			widget_control, viewstate.yval, set_value=string(decstring[0], format='(a11)')
			widget_control, viewstate.xpix, set_value=string(xpix,Format='(I3)')
			widget_control, viewstate.ypix, set_value=string(ypix,Format='(I3)')
			widget_control, viewstate.ival, set_value=string(ival,Format='(F10.3)')
		endelse
	end
	else:	

endcase


end


pro cleanview_event, event

common gridview2_state
common cleanview_state
common gridstate
common clean

widget_control, viewstate.slider, get_value=channel
channel = long(channel[0])

widget_control, event.id, get_uvalue=uvalue, get_value=value

change_to = 1
case uvalue of
	; Update channel stuff
	'slider': begin
		; Get current chan
		widget_control, event.id, get_value=chnstring
		channel=long(chnstring[0])

		; Update the info base
		widget_control, viewstate.vpix, set_value=string(chnstring,Format='(I4)')
		widget_control, viewstate.vval, set_value=string(grid.velarr[channel],Format='(F7.1)')

		cleanview_updatemap, channel
	end	

	'disps': begin
		if viewdisp.DispType NE 1 then begin
			viewdisp.DispType = 1
			widget_control, viewstate.DispSpec, Set_Button=1
			widget_control, viewstate.DispCont, Set_Button=0
			cleanview_updatemap, channel, /reset
		endif
	end
	'dispc': begin
		if viewdisp.DispType EQ 1 then begin
			viewdisp.DispType = 0
			widget_control, viewstate.DispSpec, Set_Button=0
			widget_control, viewstate.DispCont, Set_Button=1
			cleanview_updatemap, channel, /reset
		endif
	end
	'smth': begin
		if NOT xregistered('cleanview_imgsmth', /NoShow) then begin
			smthtitle = strcompress('Image Smoothing')
			smth_base =  widget_base(group_leader=viewstate.baseID, /column, /base_align_right, $
				title=smthtitle, uvalue = 'smth_base')
			smth_base1 = widget_base(smth_base, /Row)
				label=widget_label(smth_base1, value='Image smooth: ')
				viewstate.ImgSmthSlider=widget_slider(smth_base1, min=0, max=10, $
                                        value=viewdisp.ImgSmth, xsize=120, uvalue='smthslider', /Drag, $
					Event_Pro='cleanview_event')
			smth_base2 = widget_base(smth_base, /Row)
				smth_base3=widget_base(smth_base2, /row, /align_left, /NonExclusive)
				viewstate.BxcSmthOnOff=widget_button(smth_base3, value='Boxcar (+/-)', $
                                             uvalue='vsmthselect', Event_Pro='cleanview_event')
				viewstate.BxcSmthValue=widget_text(smth_base2, value='5', xsize=4, /Editable, $
						uvalue='vsmthtext', Event_Pro='cleanview_event')
				label=widget_label(smth_base2, value='channels  ')
			smth_done = widget_button(smth_base, value = ' Done ', uvalue = 'about_done', $
				Event_Pro='cleanview_event')
			widget_control, smth_base, /realize
			xmanager, 'cleanview_imgsmth', smth_base, /no_block

			if viewdisp.BxcSmth NE 0 then $
				widget_control, viewstate.BxcSmthOnOff, Set_Button=1
		endif
	end
	'smthslider': begin
		widget_control, event.id, get_value=smthlevel
		viewdisp.ImgSmth=long(smthlevel[0])
		cleanview_updatemap, channel
	end
	'vsmthselect': begin
		if viewdisp.BxcSmth EQ 0 then begin
			widget_control, viewstate.BxcSmthValue, get_value=boxc
			viewdisp.BxcSmth = long(boxc[0])
		endif else begin
			viewdisp.BxcSmth = 0
		endelse
		cleanview_updatemap, channel
	end
	'vsmthtext': begin
		widget_control, event.id, get_value=boxc
		viewdisp.BxcSmth = long(boxc[0])
		cleanview_updatemap, channel
	end

	'blue': begin
		if viewdisp.colortab NE 1 then begin
			viewdisp.colortab = 1
			widget_control, viewstate.colorBLU, Set_Button=1
			widget_control, viewstate.colorGRN, Set_Button=0
			widget_control, viewstate.colorRED, Set_Button=0
			cleanview_updatemap, channel, /reset
		endif
	end
	'green': begin
		if viewdisp.colortab NE 8 then begin
			viewdisp.colortab = 8
			widget_control, viewstate.colorBLU, Set_Button=0
			widget_control, viewstate.colorGRN, Set_Button=1
			widget_control, viewstate.colorRED, Set_Button=0
			cleanview_updatemap, channel, /reset
		endif
	end
	'red': begin
		if viewdisp.colortab NE 3 then begin
			viewdisp.colortab = 3
			widget_control, viewstate.colorBLU, Set_Button=0
			widget_control, viewstate.colorGRN, Set_Button=0
			widget_control, viewstate.colorRED, Set_Button=1
			cleanview_updatemap, channel, /reset
		endif
	end

	; Update catalog information
	'agc_cz': begin
		if viewdisp.UseKA EQ 1 then change_to = 0
		viewdisp.UseKA = change_to
		widget_control, event.id, Set_Button=change_to
		cleanview_updatemap, channel
	end
	'agc_ncz': begin
		if viewdisp.UseUA EQ 1 then change_to = 0
		viewdisp.UseUA = change_to
		widget_control, event.id, Set_Button=change_to
		cleanview_updatemap, channel
	end
	'cathvc': begin
		if viewdisp.UseHV EQ 1 then change_to = 0
		viewdisp.UseHV = change_to
		widget_control, event.id, Set_Button=change_to
		cleanview_updatemap, channel
	end
	'cat3d': begin
		if viewdisp.Use3D EQ 1 then change_to = 0
		viewdisp.Use3D = change_to
		widget_control, event.id, Set_Button=change_to
		cleanview_updatemap, channel
	end
	'catchk': begin
		if viewdisp.UseCk EQ 1 then change_to = 0
		viewdisp.UseCk = change_to
		widget_control, event.id, Set_Button=change_to
		cleanview_updatemap, channel
	end
	'catcnt': begin
		if viewdisp.UseCt EQ 1 then change_to = 0
		viewdisp.UseCt = change_to
		widget_control, event.id, Set_Button=change_to
		cleanview_updatemap, channel
	end

	; Update Colorbar and display
	'colorbar_ascale': begin
		if viewdisp.UseAuto EQ 0 then begin
			; Set parameters in viewdisp to avoid calls to widget_info over and over again
			viewdisp.UseCons = 0
			viewdisp.UseAuto = 1
			; Update the buttons
			widget_control, event.id, Set_Button=1
			widget_control, viewstate.ccscale, set_button=0
			; Render the map with the new settings
			cleanview_updatemap, channel, /reset
		endif
	end
	'colorbar_cscale': begin
		if viewdisp.UseCons EQ 0 then begin
			; Set parameters in viewdisp to avoid calls to widget_info over and over again
			viewdisp.UseCons = 1
			viewdisp.UseAuto = 0
			; Switch over to linear scale
			viewdisp.UseLin = 1
			viewdisp.UseLog = 0
			viewdisp.UseHeq = 0
			; Update the buttons
			widget_control, event.id, Set_Button=1
			widget_control, viewstate.acscale, set_button=0
			widget_control, viewstate.linear, set_button=1
			widget_control, viewstate.logari, set_button=0
			widget_control, viewstate.histeq, set_button=0
			; Render the map with the new settings
			cleanview_updatemap, channel, /reset
		endif
	end

	;Update scaling and display
	'linscl': begin
		if viewdisp.UseLin EQ 0 then begin
			viewdisp.UseLin = 1
			viewdisp.UseLog = 0
			viewdisp.UseHeq = 0

			widget_Control, event.id, Set_Button=1
			widget_control, viewstate.logari, set_button=0
			widget_control, viewstate.histeq, set_button=0
			
			cleanview_updatemap, channel, /reset
		endif
	end
	'logscl': begin
		if viewdisp.UseLog EQ 0 then begin
			viewdisp.UseLin = 0
			viewdisp.UseLog = 1
			viewdisp.UseHeq = 0
			; Switch over to an automatic scale for this one
			viewdisp.UseCons = 0
			viewdisp.UseAuto = 1

			widget_Control, event.id, Set_Button=1
			widget_control, viewstate.linear, set_button=0
			widget_control, viewstate.histeq, set_button=0
			; Force to autoscale
			widget_control, viewstate.ccscale, set_button=0
			widget_control, viewstate.acscale, set_button=1

			cleanview_updatemap, channel, /reset
		endif
	end
	'heqscl': begin
		if viewdisp.UseHeq EQ 0 then begin
			viewdisp.UseLin = 0
			viewdisp.UseLog = 0
			viewdisp.UseHeq = 1
			; Switch over to an automatic scale for this one
			viewdisp.UseCons = 0
			viewdisp.UseAuto = 1

			widget_Control, event.id, Set_Button=1
			widget_control, viewstate.logari, set_button=0
			widget_control, viewstate.linear, set_button=0
			; Force to autoscale
			widget_control, viewstate.ccscale, set_button=0
			widget_control, viewstate.acscale, set_button=1
		
			cleanview_updatemap, channel, /reset
		endif
	end	

	; Update the variables that control cleaning
	'limitf': begin
		cleanstate.flux = float(value)
	end
	'limits': begin
		cleanstate.sigma = float(value)
	end
	'limiti': begin
		cleanstate.niter = float(value)
	end
	'limitc': begin
		cleanstate.crange = long( strsplit(value, ',', /Extract) )
	end
	'cleansum': begin
		if strcmp(value, 'No ') then begin
			cleanstate.allowsum = 1
			widget_control, event.id, set_value='Yes'
		endif else begin
			cleanstate.allowsum = 0
			widget_control, event.id, set_value='No '
		endelse
	end

	; Handle the beam thing
	'beam': begin
		widget_control, /hourglass
		Widget_Control, viewstate.display_win, Clear_Events=1
		beamview
	end

	'clean': begin
		; Retrive value in text boxes
		widget_control, viewstate.clean_flux, get_value=value
			cleanstate.flux = value
		widget_control, viewstate.clean_sigma, get_value=value
			cleanstate.sigma = value
		widget_control, viewstate.clean_niter, get_value=value
			cleanstate.niter = value
		widget_control, viewstate.clean_chn, get_value=value
			cleanstate.crange = long( strsplit(value, ',', /Extract) )
		widget_control, viewstate.clean_sum, get_value=value
			if strcmp(value,'No ') then cleanstate.allowsum = 0 else cleanstate.allowsum = 1

		; Clear our the previous log
		widget_control, viewstate.summary, set_value='Flux: ---- Iteration: ---- Residual: ----'
		widget_control, viewstate.clog, set_value='----'

		widget_control, viewstate.display_win, get_value=index
		wset, index
		mask = bytarr(110, 50)
		tv, mask, 110, 340		
		plots, [110, 110, 220, 220, 110], [340, 390, 390, 340, 340], Color=FSC_Color('Yellow'), /Device
		xyouts, 120, 360, 'Cleaning....', /Device, CharSize=1.5, Color=FSC_Color('Yellow')

		; Set and hourglass and run with it
		widget_control, /Hourglass

		; Reset cleanexit
		cleanstate.cleanexit *= 0
		; And clean it...
; 		alfa_clean7, dbox2, cleanstate.beam, cleand_temp, $
; 			Continuum=cbox2, c_cleand=c_cleand_temp, $
; 			MapRMS=cleanstate.rms, $
; 			NIter=cleanstate.niter, Flux=cleanstate.flux,  Sigma=cleanstate.sigma, Gain=0.1, $
; 			Range=cleanstate.crange, AllowSum=cleanstate.allowsum, Exit_Status=exit_status, $
; 			FWHM=cleanstate.restore_fwhm, PA=cleanstate.restore_pa, /Silent, Output=log

		alfa_clean7, dbox2, cleanstate.beam, cleand_temp, $
			Continuum=cbox2, c_cleand=c_cleand_temp, $
			MapRMS=cleanstate.rms, $
			NIter=cleanstate.niter, Flux=cleanstate.flux,  Sigma=cleanstate.sigma, Gain=0.1, $
			Range=cleanstate.crange, AllowSum=cleanstate.allowsum, Exit_Status=exit_status, $
			FWHM=cleanstate.restore_fwhm, PA=cleanstate.restore_pa, /Silent, Output=log, /CExts

		; Apply the F3 correction here
		for i=cleanstate.crange[0],cleanstate.crange[1] do $
			cleand_temp[i,*,*] = cleand_temp[i,*,*] / cleanstate.area_corr
		c_cleand_temp = c_cleand_temp / cleanstate.area_corr

		; Because we have added a buffer between the selected box and cleaned box, we need to 
		; make sure that we make dimensions match in the end.  That is why alfa_clean6 returns
		; to cleand_temp and c_cleand_temp instead of to the actual common variables.
		;+ Extraction variable decleration here
		lower_x = cleanstate.map_buffer
		lower_y = cleanstate.map_buffer
		upper_x = n_elements(cleand_temp[0,*,0])-1-cleanstate.map_buffer
		upper_y = n_elements(cleand_temp[0,0,*])-1-cleanstate.map_buffer
		;+ Extraction here
		cleand = cleand_temp[*, lower_x:upper_x, lower_y:upper_y]
		c_cleand = c_cleand_temp[lower_x:upper_x, lower_y:upper_y]
		;+ Freeing temporary variables
		DelVarX, cleand_temp, c_cleand_temp

		; Save output parameters to the correct places
		cleanstate.cleanlog = log
		cleanstate.cleanexit[cleanstate.crange[0]:cleanstate.crange[1]] = exit_status
		junk = where(exit_status EQ -1, junk1)
		junk = where(exit_status EQ +1, junk2)
		cleanstate.cleansum = [junk1, junk2]

		; Update main window and display it
		;+ Update log window
		widget_control, viewstate.summary, set_value='Flux: '+string(junk2,Format='(I4)')+ $
			'     Iteration: '+string(junk1,Format='(I4)')
		widget_control, viewstate.clog, set_value=log
		;+ Update display
		cleanview_updatemap, channel
	end

	'vfeld': viewstate.drawfluxbox = 2

	'flux': viewstate.drawfluxbox = 1

	'qstart': begin
		h=['CLEANView Quickstart',$
		'CONTENTS',$
		'',$
		'    * Overview',$
		'    * Menu Options',$
		'    * Main Window',$
		'    * Colorbar',$
		'    * Information Display',$
		'    * Controls',$
		'    * Log Window',$
		'    * Other Notes',$
		'',$
		'Overview',$
		'Clean is a set of procedures written in IDL to assist in the ',$
		'deconvolution of ALFA data, primarily those of the ALFALFA ',$
		'extragalactic survey. CleanView acts as a tool that interfaces with ',$
		'GRIDview2, build_beam3, and alfa_clean5 to select data for deconvolution, ',$
		'model the local beam shape, and deconvolve the data.  The deconvolution ',$
		'process can also be fully controlled from this program.',$
		'',$
		'',$
		'Menu Options',$
		'File 	',$
		'',$
		'    * Export JPEG - allows user to export JPEGs in single channel or averaged ',$
		'                       channel maps in ',$
		'    * Save Data - Saves the deconvolved data along with the local beam shape',$
		'                     to an IDL save file that can be specified.',$
		'    * Exit - Exit the program.',$
		'',$
		'Color 	',$
		'',$
		'    * Blue/White - Set the color table to a linear blue-white stretch.',$
		'    * Green/White - Set the color table to a linear green-white stretch.',$
		'    * Red/White - Set the color table to a linear red-white stretch.',$
		'',$
		'Scaling 	',$
		'',$
		'    * Linear - Fixed minimum and maximum scaling for auto or constant scaling.',$
		'    * Logarithmic - Scaling where offset is calculated and subtracted before taking a ',$
		'                       base 10 logarithm.  This forced auto scaling to be turned on.',$
		'    * Histogram EQ - Histogram equalization of the image.  This forces auto scaling to ',$
		'                        turned on',$
		'',$
		'Colorbar 	',$
		'',$
		'    * Autoscale - minimum and maximum values are determined using the full ',$
		'                     range of the image.  This option is available for all options ',$
		'                     in the Scaling menu.',$
		'    * Constant Scale - minimum and maximum values are taken from input ',$
		'                          in the settings menu (default is -5 to +10).  This ',$
		'                          is only available with the Linear scale option.',$
		'',$
		'Catalogs 	',$
		'',$
		'    * AGC (known cz) - known cz catalog overlay.  If the galaxy has 21 cm line',$
		'                         information, OR has been searched before for HI, ', $
		'                         it will be displayed as a yellow box.  Hover your mouse', $
		'                         over the box to see the DETCODE at the bottom of the ',$
		'                         GRIDview2 window', $
		'    * AGC (no known cz) - unknown cz catalog overlay',$ 
		'    * 3D extractor - pass a keyword at start-up using:  ',$
		'                                 gridview2, grid, cat3D=filename_string ',$
		'                     Sources will be displayed as a red cross with the source number', $
		'',$
		'Help 	',$
		'',$
		'    * Quickstart - Text version of this guide.',$
		'    * About - Information, update information, and modification history.',$
		'',$
		'',$
		'Main Window',$
		'The main window display shows a color image of the current spectral channel ',$
		'   if the map has been deconvolved.  If not, the display are will show a solid ',$
		'   color.  Any selected catalogs will be overlaid on the map. ',$
		'',$
		'Colorbar',$
		'The colorbar window shows the current dynamic range for the currently selected',$ 
		'colortable (available in the Settings menu). ',$
		'Colorbar can autoscale for the currently displayed map, or use a fixed ',$
		'constant scale that can be set in the settings menu.',$
		'',$
		'Information Display',$
		'The left column of the information display show information from the ',$
		'current grid structure - The current coordinates for the mouse pointer ',$
		'in the main window - RA/DEC, velocity of the current channel, and ',$
		'x/y/intensity of the current cube pixel.',$
		'',$
		'Controls',$
		'This column controls all aspects of the deconvolution process.  The grid RMS ',$
		'is computed using GSFC IDL routine robust_sigma and cannot be changed.  The ',$
		'flux, flux sigma, and iteration options control how deeply each channel is ',$
		'cleaned.  "Flux" causes the cleaning loop to exit when the peak flux is below ',$
		'this value.  "Flux sigma" exits when the peak flux is below some multiple of the ',$
		'grid RMS.  "Iteration" is a simple iteration limit.  The channel selection option ',$
		'confines the deconvolution to a set range of channels.  For some sources, the ',$
		'signal strength is too low for a through cleaning.  To get around this, the ',$
		'"combine channels" can be set.  This combines channels until at least 75% of ',$
		'channels can be cleaned at least once.  Once the cleaning is done, the combined ',$
		'channels are split into constituent channels again.',$
		'',$
		'Log Window',$
		'Once the deconvolution has finished, details of the process are display in ',$
		'the log window.  These details include the how many channels were combined ',$
		'for the cleaning process, the exit status of each channel (flux limited, ',$
		'iteration limited, or residuals limited), and the total run time.  This ',$
		'information is also saved when the CLEANed data cube is saved.',$
		'',$
		'Other Notes',$
		'* For most sources, the default values are good.  However, for bright sources or ',$
		'     sources near Galactic HI a higher iteration limit of 5,000 to 10,000 should ',$
		'     be used.  This will slow down the deconvolution but yield better results',$
		'* If the channel range is a prime number and the "combine channel" options is set ',$
		'     then the routine will combine all channels into one.  Be sure to check the ',$
		'     log window to see if this has happened.',$
		'',$
		'',$
		'Jayce Dowell, Indiana University.']
			
		if NOT xregistered('cleanview_quickstart', /NoShow) then begin
			helptitle = strcompress('CLEANview Quickstart')
			help_base =  widget_base(group_leader=viewstate.baseID, /column, /base_align_right, $
				title=helptitle, uvalue = 'help_base')
			help_text = widget_text(help_base, /scroll, value = h, xsize = 90, ysize = 50)
			help_done = widget_button(help_base, value = ' Done ', uvalue = 'about_done', $
				event_pro='cleanview_event')
			widget_control, help_base, /realize
			xmanager, 'cleanview_quickstart', help_base, /no_block
		endif
	end
	'about': begin
		h = ['CleanView - Graphical deconvolution utility', $
		     '  Written, November 2007', $
		     ' ', $
		     '  Last update, Tuesday, March 18, 2008']
		
		if NOT xregistered('cleanview_help', /NoShow) then begin
			about_base =  widget_base(group_leader=viewstate.baseID, title='About CleanView', $
				/Column, /Base_Align_Right, uvalue = 'about_base')
			about_text = widget_text(about_base, /Scroll, value = h, xsize = 45, ysize = 10)
			about_done = widget_button(about_base, value = ' Done ', uvalue = 'about_done', $
				event_pro='cleanview_event')
			
			widget_control, about_base, /realize
			xmanager, 'cleanview_help', about_base, /no_block
		endif
	end
	'about_done': begin
		widget_control, event.top, /destroy
	end

	'dsave': begin
		if NOT xregistered('cleanview_save', /NoShow) then begin
			cd, current=dir
			save_base = widget_base(group_leader=viewstate.baseID, title=' CLEANview Save Data ', /Column, $
				uvalue='save_base')
			  top_base = widget_base(save_base, /Column, /Base_Align_Left)
			    save_text = widget_label(top_base, value='Output Directory:    ')
			    viewstate.save_dir = widget_text(top_base, value=dir+'/', uvalue='save_dir', $
				/Editable, XSize=40)
			    save_text = widget_label(top_base, value='Output File Name:    ')
			    viewstate.save_name = widget_text(top_base, value='cleanview.sav', uvalue='save_name', $
				/Editable, XSize=40)
			  bottom_base = widget_base(save_base, /Row, /Base_Align_Right)
			    save_save = widget_button(bottom_base, value = ' Save ',   uvalue = 'save_save')
			    save_done = widget_button(bottom_base, value = ' Cancel ', uvalue = 'save_cancel')

			widget_control, save_base, /realize
			xmanager, 'cleanview_save', save_base, /no_block
		endif
	end 

	'exit': begin
		widget_control, event.top, /destroy

		hor
		ver
		loadct, 1, /silent
		
		delvarx, viewstate, viewdisp
		delvarx, cleanstate, dbox2, cbox2, cleand, c_cleand

		widget_control, state.BaseID, /sensitive
		
		print, 'Exiting CleanView...'
	end
else:
endcase
		
end


pro cleanview_save_event, event

common cleanview_state
common gridstate
common clean

widget_control, event.id, get_uvalue=uvalue, get_value=value

case uvalue of
	'save_save': begin
		widget_control, viewstate.save_dir,  get_value=dir
		widget_control, viewstate.save_name, get_value=filename
		
		widget_control, /Hourglass
		filename = dir+filename
		save, cleanstate, cleand, c_cleand, filename=filename
		
		widget_control, event.top, /destroy
	end
	'save_cancel': begin
		widget_control, event.top, /destroy
	end
	else:	
endcase

end


pro cleanview_updatemap, channel, Reset=Reset

common cleanview_state
common gridstate
common clean

; Update the color table
device, decomposed=0
loadct,viewdisp.colortab,/Silent

; Update axis
widget_control, viewstate.display_win, get_value=index
wset, index

hor, viewstate.ramax, viewstate.ramin
ver, viewstate.decmin, viewstate.decmax 

if Keyword_Set(Reset) then begin
	plot, [0,0], /nodata, xstyle=1, ystyle=1, xtick_get=xvals, ytick_get=yvals, $
		ticklen=ticklen, charsize=1.0, Position=[0.15,0.15,0.95,0.95]
	nxticklabels=n_elements(xvals)
	nyticklabels=n_elements(yvals)
	xspacing=((xvals[n_elements(xvals)-1]-xvals[0])*60.0)/(nxticklabels-1)
	yspacing=((yvals[n_elements(yvals)-1]-yvals[0])*60.0)/(nyticklabels-1)
	xticlabs=ratickname(xvals*15.0)
	yticlabs=dectickname(yvals)
	
	viewstate.px = !X.WINDOW * !D.X_VSIZE
	viewstate.py = !Y.WINDOW * !D.Y_VSIZE
	viewstate.sx = viewstate.px[1] - viewstate.px[0] - 1
	viewstate.sy = viewstate.py[1] - viewstate.py[0] - 1
	
	erase
	
	hor, viewstate.ramax, viewstate.ramin
	ver, viewstate.decmin, viewstate.decmax
endif else begin
	!X.s = viewstate.xs
	!Y.s =viewstate.ys
	!X.region = viewstate.xregion
	!Y.region = viewstate.yregion
endelse

; Check display scaling and colorbar mapping
if viewdisp.DispType then begin
	lowerchan = ((channel-viewdisp.BxcSmth)>0)
	upperchan = ((channel+viewdisp.BxcSmth)<(n_elements(grid.velarr)-1))
	
	to_disp = reform(total(cleand[lowerchan:upperchan,*,*],1))/(2.0*viewdisp.BxcSmth+1)
endif else begin
	to_disp = c_cleand
endelse
scale_rule = viewdisp.UseLin + viewdisp.UseLog*2 + viewdisp.UseHeq*4
case scale_rule of
	1:	begin
			title_cb = 'Flux Density [mJy/beam]'
			format_cb = '(I4)'
			
			delta = max(to_disp) - min(to_disp)
			d_min = min(to_disp) + delta/5.0
			d_max = max(to_disp) - delta/5.0
	end
	2:	begin
			offset = min(to_disp) - (max(to_disp)-min(to_disp))*0.01
			to_disp = alog10(to_disp - offset)
		
			title_cb = 'Log Flux Density [log(mJy/beam)]'
			format_cb = '(F5.2)'
		
			delta = max(to_disp) - min(to_disp)
			d_min = ((min(to_disp) + delta/5.0)>0)
			d_max = (0.1>(max(to_disp) - delta/5.0))
	end
	4:	begin
			to_disp = hist_equal(to_disp, minv=min(to_disp), maxv=max(to_disp))
		
			title_cb = 'Flux Density [mJy/beam]'
			format_cb = '(I3)'
		
			delta = max(to_disp) - min(to_disp)
			d_min = min(to_disp) + delta/5.0
			d_max = max(to_disp) - delta/5.0
	end
	else:
endcase

if viewdisp.UseCons EQ 1 then begin
	d_min = -5.0
	d_max = 10.0
endif

tv, reverse( bytscl(smooth(congrid(to_disp, viewstate.sx, viewstate.sy), viewdisp.ImgSmth, /Edge_Truncate), min=d_min, max=d_max) ), $
	viewstate.px[0]+2, viewstate.py[0]+2

if Keyword_Set(Reset) then begin
	device, decomposed=1
	plot, [0,0], /nodata, /NoErase, xtitle='RA [HMS] J2000', ytitle='Dec [DMS] J2000', $
		xstyle=1, ystyle=1, Position=[0.15,0.15,0.95,0.95], xtickn=xticlabs, ytickn=yticlabs, $
		ticklen=-0.01, charsize=1.0, Color='00FFFF'XL
	device, decomposed=0
endif

if viewdisp.DispType then begin
	if viewdisp.BxcSmth EQ 0 then begin
		infostring = 'Channel '+strcompress(channel, /Remove_All)+ $
			'    Velocity= '+strcompress(grid.velarr[channel])+' km/s   '
		
			case cleanstate.cleanexit[channel] of
				0:	infostring = infostring+'  NOT Cleaned'
				-1:     infostring = infostring+'  Cleaned; iteration limited'
				1:	infostring = infostring+'  Cleaned; flux limited'
			endcase
	endif else begin
		infostring = 'Channels '+strcompress(lowerchan, /remove_all)+' to '+strcompress(upperchan, /remove_all)+ $
             			'  Smoothed over -> '+strcompress(string(grid.velarr[lowerchan],Format='(F7.1)'), /remove_all)+' to '+ $
				strcompress(string(grid.velarr[upperchan],Format='(F7.1)'), /remove_all)+' km/s   '

		case cleanstate.cleanexit[channel] of
			0:	infostring = infostring+'  NOT Cleaned'
			-1:     infostring = infostring+'  Iteration limited'
			1:	infostring = infostring+'  Flux limited'
		endcase
	endelse
endif else begin
	infostring = 'Contiuum'
endelse

device, decomposed=1
mask=fltarr(600,25)
tvscl, mask, 10,5

xyouts, 10, 10, infostring, /device, Color='00FFFF'XL, charsize=1.1
device, decomposed=0

; Save window parameters
viewstate.xs = !X.s
viewstate.ys = !Y.s
viewstate.xregion = !X.region
viewstate.yregion = !Y.region

viewdisp.CBMin = d_min
viewdisp.CBMax = d_max
widget_control, viewstate.colorbar_win, get_value=index
wset, index
erase

device, decomposed=1
colorbar, Range=[d_min, d_max], Position=[0.50, 0.10, 0.90, 0.90], /Ver, $
	Title=title_cb, Format=format_cb, Color='00FFFF'XL
device, decomposed=0


; Update catalog information
cleanview_updatecat
end




















pro cleanview_updatecat

common cleanview_state
common gridstate
common clean
common gridview2_state

device, decomposed=1

widget_control, viewstate.slider, get_value=chnstring
channel=long(chnstring[0])

lc = ((channel-14) > 0 )
uc = ((channel+14) < (n_elements(grid.velarr)-1))

widget_control, viewstate.display_win, get_value=index
wset, index

hor, viewstate.ramax, viewstate.ramin
ver, viewstate.decmin, viewstate.decmax

if viewdisp.UseKA AND viewdisp.DispType then begin
	iagc=where(state.rahr GT viewstate.ramin AND $
		   state.rahr LT viewstate.ramax AND $
		   state.decdeg GT viewstate.decmin AND $
		   state.decdeg LT viewstate.decmax AND $
		   ((state.agc.v21 lt grid.velarr[channel]+state.agc.width/2.0 AND $
		     state.agc.v21 gt grid.velarr[channel]-state.agc.width/2.0 AND $
		     state.agc.v21 ne 0) OR $
		    (state.agc.vopt lt grid.velarr[lc] AND $
		     state.agc.vopt gt grid.velarr[uc] AND $
		     state.agc.vopt ne 0)))
	color = '00FFFF'XL
	;color=FSC_Color('Yellow')

	if (iagc[0] ne -1) then begin
		plots, state.rahr[iagc], state.decdeg[iagc], PSym=6,SymSize=2, Color=color, /Data

		resultcoords=convert_coord(state.rahr[iagc], state.decdeg[iagc], /data, /double, /to_device)
		xyouts,resultcoords[0,*]+10, resultcoords[1,*]-10, strcompress(state.agc.agcnumber[iagc], /remove_all), $
			Color=color, /Device
	endif
endif

;Plot AGC galaxies with no known redshifts in the area
if viewdisp.UseUA AND viewdisp.DispType  then begin
	iagc=where(state.rahr GT viewstate.ramin AND state.rahr LT viewstate.ramax AND $
		state.decdeg GT viewstate.decmin AND state.decdeg LT viewstate.decmax AND $
           	state.agc.v21 eq 0 AND state.agc.vopt eq 0 )
	
	color = '0000FF'XL
	;color=FSC_Color('Red')
; 
	if (iagc[0] ne -1) then begin
		plots, state.rahr[iagc], state.decdeg[iagc], PSym=4, SymSize=3, Color=Color, /Data
		
		resultcoords=convert_coord(state.rahr[iagc], state.decdeg[iagc], /data, /double, /to_device)
		xyouts,resultcoords[0,*]+10, resultcoords[1,*]-10, strcompress(state.agc.agcnumber[iagc], /remove_all), $
			Color=color, /Device
	endif
endif

; Plot HVC data
if viewdisp.UseHV AND viewdisp.DispType then begin
	idhbb=where(state.dhbb._raj2000/15.0 lt viewstate.ramax AND state.dhbb._raj2000/15.0 gt viewstate.ramin AND $
            state.dhbb._dej2000 lt viewstate.decmax AND state.dhbb._dej2000 gt viewstate.decmin AND $
            (state.dhbb_vhelio lt grid.velarr[channel]+state.dhbb.fwhm/2.0 AND $
             state.dhbb_vhelio gt grid.velarr[channel]-state.dhbb.fwhm/2.0))

	color = '0000FF'XL
	;color=FSC_Color('Red')

	if (idhbb[0] ne -1) then begin
		plots, state.dhbb._raj2000[idhbb]/15.0,  state.dhbb._dej2000[idhbb], PSym=6, SymSize=2, $
			Color=color, /Data

		resultcoords=convert_coord(state.dhbb._raj2000[idhbb]/15.0,state.dhbb._dej2000[idhbb], /data, /double, /to_device)
		xyouts,resultcoords[0,*]+10, resultcoords[1,*]-10, 'HVC '+strcompress(state.dhbb.seq[idhbb], /remove_all), $
			color=color, /device
	endif
endif

;Load from 3D extractions
if viewdisp.Use3D AND viewdisp.DispType then begin
	;modified to accept a larger width
	i3D=where(state.table3D.ra lt viewstate.ramax AND state.table3D.ra gt viewstate.ramin AND $
           state.table3D.dec lt viewstate.decmax AND state.table3D.dec gt viewstate.decmin AND $
           state.table3D.cz lt grid.velarr[lc] AND state.table3D.cz gt grid.velarr[uc])

	color = '0000FF'XL
	;color=FSC_Color('Red')
	
	if (i3D[0] ne -1) then begin
		plots, state.table3D.ra[i3D], state.table3D.dec[i3D], PSym=1, SymSize=2, Color=color, Thick=2.0, /Data
		
		resultcoords=convert_coord(state.table3D.ra[i3D],state.table3D.dec[i3D], /data, /double, /to_device)
		xyouts,resultcoords[0,*]-25, resultcoords[1,*]+10, strcompress(state.table3D.id[i3D], /remove_all), color=color, /device
	endif
endif

; Load in checklist.txt
if viewdisp.UseCk AND viewdisp.DispType then begin
	; Find checklist.txt file, if it exists
	check=findfile('checklist.txt', count=count)
    	if (count eq 1) then begin
		; Read in data
		data=read_ascii('checklist.txt', delimiter=',')

		rahr=data.field1[0,*]
		decdeg=data.field1[1,*]
		vel=data.field1[2,*]
		width=data.field1[3,*]
		
		icheck=where(rahr lt viewstate.ramax AND rahr gt viewstate.ramin AND $
			decdeg lt viewstate.decmax AND decdeg gt viewstate.decmin AND $
			( vel lt (grid.velarr[lc]+1.2*width/2.0) AND vel gt (grid.velarr[uc]-1.2*width/2.0) ) )

		color = '0000FF'XL
		;color=FSC_Color('Red')

		if (icheck[0] NE -1) then $
			plots, rahr[icheck], decdeg[icheck], PSym=4, SymSize=2, Thick=2.0, Color=color, /Data
	endif
endif

; Strong continuum sources
if viewdisp.UseCt then begin
	;Peak ranges in Jy/beam
	peaklower=.25
	peakupper=1.0e6
	
	invss=where(state.nvsscat.ra_2000_ lt viewstate.ramax*15.0 AND state.nvsscat.ra_2000_ gt viewstate.ramin*15.0 AND $
		state.nvsscat.dec_2000_ lt viewstate.decmax AND state.nvsscat.dec_2000_ gt viewstate.decmin AND $
		state.nvsscat.peak_int ge peaklower AND state.nvsscat.peak_int lt peakupper ) 
	
	color = 'FFFFFF'XL
	;color=FSC_Color('White')
	
	if (invss[0] NE -1) then $
		plots,state.nvsscat[invss].ra_2000_/15.0, state.nvsscat[invss].dec_2000_ , $
			PSym=6,SymSize=2.0,Thick=2.0, Color=color, /Data
endif

device, decomposed=0

end



pro cleanview_jpeg, event

common cleanview_state
common gridstate

widget_control, viewstate.slider, get_value=chnstring

if NOT xregistered('cleanview_jpeg', /NoShow) then begin
	jpeg_output_base = widget_base(group_leader=viewstate.baseID, /row, /base_align_right, $
		title='CLEANview JPEG Output', uvalue = 'jpeg_output_base')
	configbase=widget_base(jpeg_output_base, /column, /align_left)
	startbase=widget_base(configbase, /row)
	label=widget_label(startbase,     value='Start Channel:    ')
	viewstate.frmsta=widget_text(startbase, xsize=8, value=string(chnstring[0], Format='(I4)'), /Editable, uvalue='start')
	numframesbase=widget_base(configbase, /row)
	label=widget_label(numframesbase, value='Number of Frames: ')
	viewstate.frmnum=widget_text(numframesbase, xsize=8, value='1', /Editable, uvalue='number')
	numstepsbase=widget_base(configbase, /row)
	label=widget_label(numstepsbase,  value='Number of steps:  ')
	viewstate.frmste=widget_text(numstepsbase, xsize=8, value='5', /Editable, uvalue='step')
	
	cd, current=pwd  ;Get current working directory into a string
	filesbase=widget_base(jpeg_output_base, /column, /align_left)
	directorybase=widget_base(filesbase, /Column)
	label=widget_label(directorybase, value='Output Directory:')
	viewstate.dir=widget_text(directorybase, xsize=40, value=pwd+'/', /Editable, uvalue='dir')
	buttonbase=widget_base(filesbase, xsize=100,ysize=100, /align_right, /column)
	jpeg_output_export = widget_button(buttonbase, value = ' Export ', $
		uvalue = 'export', event_pro='cleanview_jpeg_event')
	cancel=widget_button(buttonbase, value=' Cancel ', uvalue='cancel', $
					event_pro='cleanview_jpeg_event')
	
	widget_control, jpeg_output_base, /realize
	xmanager, 'cleanview_jpeg', jpeg_output_base, /no_block
endif

end




pro cleanview_jpeg_event, event

common cleanview_state
common gridstate
common clean

widget_control, event.id, get_uvalue = uvalue

case uvalue of 
	'export': begin
		widget_control, viewstate.slider, get_value=chnstring
		widget_control, viewstate.frmsta, get_value=starts
		widget_control, viewstate.frmnum, get_value=nums
		widget_control, viewstate.frmste, get_value=steps
		widget_control, viewstate.dir,    get_value=dir

		start = long(starts[0])
		num = long(nums[0])
		step = long(steps[0])
		stop=start+(num-1)*step

		;Beginning for each file
		widget_control, /Hourglass
		for channel=start, stop, step do begin
			thisDevice=!d.name
			set_plot, 'Z', /COPY
			
			Device, Set_Resolution=[700, 800], Z_Buffer=0
			Erase
						
			lowerchan=channel-step
			upperchan=channel+step
			if (lowerchan lt 0) then lowerchan=0
			if (upperchan gt n_elements(grid.velarr)-1) $
				then upperchan=n_elements(grid.velarr)-1
			
			filename='clean_vel_'+strcompress(long(grid.velarr[upperchan]), /remove_all)+'_'+$
				strcompress(long(grid.velarr[lowerchan]), /remove_all)+'.jpg'
			
			filename=dir+filename
			
			loadct, viewdisp.colortab, /silent
			to_disp = total( reform(cleand[lowerchan:upperchan, *, *]), 1)
			to_disp = to_disp / float(step+1)
			
			; Check display scaling and colorbar mapping
			if viewdisp.UseLin EQ 1 then begin
				to_disp = to_disp
			
				title_cb = 'Flux Density [mJy/beam]'
				format_cb = '(I4)'
				
				delta = max(to_disp) - min(to_disp)
				d_min = min(to_disp) + delta/5.0
				d_max = max(to_disp) - delta/5.0
			endif
			
			if viewdisp.UseLog EQ 1 then begin
				offset = min(to_disp)-(max(to_disp)-min(to_disp))*0.01
				to_disp = alog10(to_disp - offset)
			
				title_cb = 'Log Flux Density [log(mJy/beam)]'
				format_cb = '(F5.2)'
			
				delta = max(to_disp) - min(to_disp)
				d_min = min(to_disp) + delta/5.0
				d_max = max(to_disp) - delta/5.0
			endif
			
			if viewdisp.UseHeq EQ 1 then begin
				to_disp = hist_equal(to_disp, minv=min(to_disp), maxv=max(to_disp))
			
				title_cb = 'Flux Density [mJy/beam]'
				format_cb = '(I3)'
			
				delta = max(to_disp) - min(to_disp)
				d_min = min(to_disp) + delta/5.0
				d_max = max(to_disp) - delta/5.0
			endif
			
			if viewdisp.UseCons EQ 1 then begin
				d_min = -5.0
				d_max = 10.0
			endif
			
			xrange = [viewstate.ramax, viewstate.ramin]
			yrange = [viewstate.decmin, viewstate.decmax]
			plot, [0,0], /nodata, xrange=xrange, yrange=yrange, xstyle=1, ystyle=1, $
				xtick_get=xvals, ytick_get=yvals, ticklen=ticklen, charsize=1.0, $
				Position=[0.15,0.203,0.95,0.953]
			nxticklabels=n_elements(xvals)
			nyticklabels=n_elements(yvals)
			xspacing=((xvals[n_elements(xvals)-1]-xvals[0])*60.0)/(nxticklabels-1)
			yspacing=((yvals[n_elements(yvals)-1]-yvals[0])*60.0)/(nyticklabels-1)
			xticlabs=ratickname(xvals*15.0)
			yticlabs=dectickname(yvals)
			plot, [0,0], /nodata, xrange=xrange, yrange=yrange, $
				xtitle='RA [HMS] J2000', ytitle='Dec [DMS] J2000', $
				xstyle=1, ystyle=1, Position=[0.15,0.256,0.95,0.956], $
				xtickn=xticlabs, ytickn=yticlabs, ticklen=-0.01, charsize=1.0
			px = !X.WINDOW * !D.X_VSIZE
			py = !Y.WINDOW * !D.Y_VSIZE
			sx = px[1] - px[0] - 1
			sy = py[1] - py[0] - 1
			
			hor, viewstate.ramax, viewstate.ramin
			ver, viewstate.decmin, viewstate.decmax
			
			tv, reverse( bytscl(congrid(reform(to_disp), sx, sy), min=d_min, max=d_max) ), $
				px[0]+1, py[0]+2	

			infostring='Channels '+strcompress(lowerchan)+' to '+strcompress(upperchan)+ $
				'  Smoothed over -> '+string(grid.velarr[lowerchan],Format='(F7.1)')+$
				' to '+string(grid.velarr[upperchan],Format='(F7.1)')+' km/s'
			xyouts, 40, 140, infostring, /Device, charsize=1.0

			infostring='Clean Settings:'
			xyouts, 40, 110, infostring, /Device, CharSize=1.0

			infostring='Local Beam:'
			xyouts, 588, 110, infostring, /Device, CharSize=1.0

			flux_lim = max([cleanstate.flux, cleanstate.rms*cleanstate.sigma], type)
			flux_type = 'Flux'
			if type EQ 1 then flux_type = 'Sigma'
			infostring=' Grid RMS: '+string(cleanstate.rms, Format='(F4.2)')+' mJy    '+ $
				' Flux Limit: '+string(flux_lim, Format='(F5.2)')+' mJy     '+ $
				' Flux Limit Type: '+flux_type
			xyouts, 40, 90, infostring, /Device, charsize=0.9

			comb_chan = 'No'
			if cleanstate.allowsum EQ 1 then comb_chan = 'Yes'
			infostring=' Iteration Limit: '+string(cleanstate.niter, Format='(I5)')+'    '+ $
				'Channel Range: '+string(cleanstate.crange[0], Format='(I4)')+' to '+ $
				string(cleanstate.crange[1], Format='(I4)')+'     '+ $
				'Combine Channels: '+comb_chan
			xyouts, 40, 70, infostring, /Device, charsize=0.9

			beam_cont = cleanstate.beam / max(cleanstate.beam)
			beam_cont = 10.0*alog10(beam_cont)
			

			hor,-12,12
			ver,-12,12
			tvimage, bytscl(beam_cont, min=-24.0, max=0.0), Pos=[0.843,0.03,0.95,0.124], $
				/Minus_One
			plot, [0,0], /NoData, XRange=[-12,12], YRange=[-12,12], /NoErase, $
				Pos=[0.843,0.03,0.95,0.124], CharSize=0.75
			
			infostring = 'Chanel '+strcompress(channel, /Remove_All)+ $
				'    Velocity= '+strcompress(grid.velarr[channel])+' km/s   '
			case cleanstate.cleanexit[channel] of
				 0:	infostring = infostring+'  NOT Cleaned'
				-1:     infostring = infostring+'  Cleaned; iteration limited'
				 1:	infostring = infostring+'  Cleaned; flux limited'
				 2:     infostring = infostring+'  Cleaned; residuals limited'
			endcase

			xyouts, 10, 10, infostring, /device, charsize=1.1

			snapshot=tvrd()
			tvlct, r,g,b, /get
			device, Z_Buffer=1
			set_plot, thisDevice
				
			image24 = BytArr(3, 700, 800)
			image24[0,*,*] = r[snapshot]
			image24[1,*,*] = g[snapshot]
			image24[2,*,*] = b[snapshot]
			
			write_jpeg, filename, image24, true=1, quality=100
		endfor

		widget_control, event.top, /destroy
	end
	'cancel': begin
		widget_control, event.top, /destroy
	end
	else: 
endcase

end
