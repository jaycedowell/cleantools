; version = 04132008
pro fieldview, llx=llx, lly=lly, urx=urx, ury=ury

; Include common blocks from cleanview and gridview2
common cleanview_state
common gridstate
common clean

; Common block initialization for fieldview
common fieldview_state, fieldstate

; Fieldstate structure
fieldstate = { baseID: 0L, specwin: 0L, fieldwin: 0L, barwin: 0L, $
		xpix: 0, ypix: 0, xval: 0, yval: 0, vpix: 0, vval: 0, $
		cpix: 0, cval: 0, fval: 0, $
		llx: llx, lly: lly, urx: urx, ury: ury, $
		xmin: 0L, ymin: 0L, xmax: 0L, ymax: 0L, $
		XRange: fltarr(2), YRange: fltarr(2), CRange: [-1.0, -1.0], $
		px: [0.0, 0.0], py: [0.0, 0.0], sx: 0.0, sy: 0.0, $
		ramin: (grid.ramin+llx*grid.deltara/3600.0), $
		ramax: (grid.ramin+urx*grid.deltara/3600.0), $
		decmin: (grid.decmin+lly*grid.deltadec/60.0), $
		decmax: (grid.decmin+ury*grid.deltadec/60.0), $
		chmin:    0, velmin: max(grid.velarr), $
		chmax: 1023, velmax: min(grid.velarr), $
		flxmin: -10.0, flxmax: 30.0, $
		mousestate: 0, cornerx: 0, cornery: 0, $
		Spec: grid.velarr, $
		zeroth: fltarr(urx-llx+1, ury-lly+1), $
		first:  fltarr(urx-llx+1, ury-lly+1), $
		colorBLU: 0L, colorGRN: 0L, colorRED: 0L, colorRAI: 0L, colortab: 0, $
		colorINV: 0L, colorinvert: 0, $
		showct: 0L, cont05: 0L, cont10: 0L, cont15: 0L, overp: 1, levels: 10, $
		clip518: 0L, clip119: 0L, clip519: 0L, clip120: 0L, cliplevel: 1d19, $
		vmin: -100.0, vmax: 100.0, $
		ClickCount: -1, $
		xs: dblarr(2), ys: dblarr(2), xregion: dblarr(2), yregion: dblarr(2), $
		dir: 0L, nam: 0L }

; Move into the coordinate system of the cleand box
dx = fieldstate.urx - fieldstate.llx
dy = fieldstate.ury - fieldstate.lly
fieldstate.llx = fieldstate.llx - cleanstate.llx
fieldstate.lly = fieldstate.lly - cleanstate.lly
fieldstate.urx = fieldstate.llx + dx
fieldstate.ury = fieldstate.lly + dy

field_main_base = widget_base(group_leader=group_leader, Title='Velocity FIELDview', mbar=top_menu, $
	/Column, /Base_Align_Left)
fieldstate.baseID = field_main_base
filemenu=widget_button(top_menu, value=' File')
  buttonjpeg=widget_button(filemenu, value=' Export JPEG ', event_pro='fieldview_jpeg')
  buttonexit=widget_button(filemenu, value=' Exit ', uvalue='close', /Separator)
colormenu=widget_button(top_menu, value=' Color ')
  fieldstate.colorBLU=widget_button(colormenu, value=' Blue/White ',  uvalue='blue', /Checked_Menu)
  fieldstate.colorGRN=widget_button(colormenu, value=' Green/White ', uvalue='green', /Checked_Menu)
  fieldstate.colorRED=widget_button(colormenu, value=' Red/White ',   uvalue='red', /Checked_Menu)
  fieldstate.colorRAI=widget_button(colormenu, value=' Rainbow ',     uvalue='rainbow', /Checked_Menu)
  fieldstate.colorINV=widget_button(colormenu, value=' Invert Colors', uvalue='invert', /Separator, /Checked_Menu)
contourmenu=widget_button(top_menu, Value=' Contours ')
  fieldstate.showct=widget_button(contourmenu, value=' Overplot Contours ', uvalue='show', /Checked_Menu)
  fieldstate.cont05=widget_button(contourmenu, value='  5 Levels ', uvalue='level05', /Checked_Menu, /Separator)
  fieldstate.cont10=widget_button(contourmenu, value=' 10 Levels ', uvalue='level10', /Checked_Menu)
  fieldstate.cont15=widget_button(contourmenu, value=' 15 Levels ', uvalue='level15', /Checked_Menu)
clipmenu=widget_button(top_menu, value=' Clipping ')
  fieldstate.clip518=widget_button(clipmenu, value=' log(N) > 18.70 ', uvalue='clip518', /Checked_Menu)
  fieldstate.clip119=widget_button(clipmenu, value=' log(N) > 19.00 ', uvalue='clip119', /Checked_Menu)
  fieldstate.clip519=widget_button(clipmenu, value=' log(N) > 19.70 ', uvalue='clip519', /Checked_Menu)
  fieldstate.clip120=widget_button(clipmenu, value=' log(N) > 20.00 ', uvalue='clip120', /Checked_Menu)
helpmenu=widget_button(top_menu, value=' Help ')
  buttonabout=widget_button(helpmenu, value=' About ', uvalue='about')

field_top_base = widget_base(field_main_base, /Row, /Frame, /Align_Center)
	label = widget_label(field_top_base, value='Xmin: ')
	fieldstate.xmin = widget_text(field_top_base, xsize=5, value='0', /Editable)
	label = widget_label(field_top_base, value='  Xmax: ')
	fieldstate.xmax = widget_text(field_top_base, xsize=5, value='1023', /Editable)
	label = widget_label(field_top_base, value='  Ymin: ')
	fieldstate.ymin = widget_text(field_top_base, xsize=5, value='-10', /Editable)
	label = widget_label(field_top_base, value='  Ymax: ')
	fieldstate.ymax = widget_text(field_top_base, xsize=5, value='30', /Editable)
	label = widget_label(field_top_base, value='  ')
	field_rescale = widget_button(field_top_base, value=' Rescale ', uvalue='rescale')

fieldstate.specwin  = widget_draw(field_main_base, XSize=525, YSize=220, Frame=1, uvalue='spect', $
	Event_Pro='fieldview_display_event', /Motion_Events, /Button_Events, /Align_Center)

field_upper_base = widget_base(field_main_base, /Row, /Frame, /Align_Center)
	label = widget_label(field_upper_base, value='Channel: ')
	fieldstate.cpix = widget_label(field_upper_base, value='----')
	label = widget_label(field_upper_base, value='   Velocity: ')
	fieldstate.cval = widget_label(field_upper_base, value='-------')
	label = widget_label(field_upper_base, value=' km/s   Flux Density: ')
	fieldstate.fval = widget_label(field_upper_base, value='-------')
	label = widget_label(field_upper_base, value=' mJy/beam')

field_middle_base = widget_base(field_main_base, /Row)
	field_left_base = widget_base(field_middle_base, /Frame, /Align_Left)
		fieldstate.fieldwin = widget_draw(field_left_base, XSize=425, YSize=425, Frame=1, $
			uvalue='field', /Align_Center)

	field_right_base = widget_base(field_middle_base, /Frame, /Align_Left)
		fieldstate.barwin = widget_draw(field_right_base, XSize=100, YSize=425, Frame=1, $
			uvalue='clbar', /Align_Center)
field_bottom_base = widget_base(field_main_base, /Row)
	field_reset = widget_button(field_bottom_base, value=' Reset ',   uvalue='reset')
	field_comp  = widget_button(field_bottom_base, value=' Compute ', uvalue='comp')

widget_control,field_main_base, /realize
xmanager, 'fieldview', field_main_base, /No_Block

fieldstate.colortab=39
widget_control, fieldstate.colorRAI, Set_Button=1
fieldstate.overp=1
widget_control, fieldstate.showct, Set_Button=1
fieldstate.levels=10
widget_control, fieldstate.cont10, Set_Button=1
fieldstate.cliplevel=1d19
widget_control, fieldstate.clip119, Set_Button=1

fieldview_updatemap, /Reset

fieldstate.crange = [-1.0, -1.0]
fieldview_updatespec, /Reset, /Init
fieldstate.clickcount = -1

!P.clip = [78, 33, 498, 208, 0, 1000]
!X.s = [0.15000001, 0.00078201367]
!X.region = [0.0357143, 0.984286]
!Y.s = [0.31000000, 0.031999999]
!Y.region = [0.00000, 0.999990]

fieldstate.xs = !X.s
fieldstate.ys = !Y.s
fieldstate.xregion = !X.region
fieldstate.yregion = !Y.region

end




pro fieldview_display_event, event

common fieldview_state
common cleanview_state
common gridstate
common clean

; The cross hair cursor from galflux3
result=convert_coord(event.x, event.y, /device, /double, /to_data)
device, decomposed=1
widget_control, fieldstate.specwin, get_value=index 

wset, index
Device, Copy=[0, 0, 575, 220, 0, 0, 20] 

y=[-1000,20000]
x=[0,16000]
oplot, intarr(2)+long(round(result[0])), y,color='7F7F7F'XL, linestyle=2   ; gray dashed line
oplot, x, fltarr(2)+result[1], color='7F7F7F'XL, linestyle=2  ; gray dashed line
device, decomposed=0

widget_control, event.id, get_uvalue=uvalue

Event_Types = ['DOWN', 'UP', 'MOTION']
This_Event = Event_Types[event.type]
case This_Event of
	'DOWN'  : begin
		if event.press eq 4 then begin
			widget_control, fieldstate.specwin, Draw_Motion_Events=1
			widget_control, fieldstate.specwin, get_value=index
								
			; Create a pixmap. Store its ID. Copy window contents into it.
			Window, 19, /Pixmap, XSize=575, YSize=220
			Device, Copy=[0, 0, 575, 220, 0, 0, index]
			
			; Get and store the static corner of the box.
			fieldstate.mousestate = 1
			fieldstate.cornerx = event.x
			fieldstate.cornery = event.y
		endif

		if event.press eq 2 then begin
			fieldstate.crange = [-1.0, -1.0]
			fieldstate.clickcount = 0

			fieldview_updatespec, /Reset, /Init
		endif

		if event.press eq 1 then begin
			result=convert_coord(event.x, event.y, /device, /double, /to_data)
			case fieldstate.clickcount of
				-1:	begin
						fieldstate.crange = [-1.0, -1.0]
						fieldstate.clickcount = 0
						fieldview_updatespec, /Reset
						break
				end
				2:	begin
						fieldstate.crange = [-1.0, -1.0]
						fieldstate.clickcount = 0
						fieldview_updatespec, /Reset
						break
				end
				else: 	begin
					fieldstate.crange[fieldstate.clickcount] = round(result[0])
					fieldview_updatespec
					fieldstate.clickcount = fieldstate.clickcount + 1
				end
			endcase
		endif
	end
	'UP'	: begin
		if fieldstate.mousestate EQ 1 then begin
			fieldstate.mousestate=0
			
			; Erase the last box drawn. Destroy the pixmap.
			widget_control, fieldstate.specwin, get_value=index
			wSet, index
			Device, Copy=[0, 0, 575, 220, 0, 0, 19]
			
			; Order the box coordinates.
			sx = Min([fieldstate.cornerx, event.x])
			dx = Max([fieldstate.cornerx, event.x])
			sy = Min([fieldstate.cornery, event.y])
			dy = Max([fieldstate.cornery, event.y])
			
			WDelete, 19
			
			widget_control, fieldstate.specwin, get_value=index
			wSet, index

			hor, fieldstate.chmin, fieldstate.chmax
			ver, fieldstate.flxmin, fieldstate.flxmax

			result1=convert_coord(sx, sy, /device, /double, /to_data)
			result2=convert_coord(dx, dy, /device, /double, /to_data)
			xmin = result1[0]
			xmax = result2[0]
			ymin = result1[1]
			ymax = result2[1]
				
			; Protect against going 'outside' the cube
			if (xmin lt 0) then xmin = cleanstate.llx
			if (xmax gt n_elements(grid.velarr)-1) then xmax = n_elements(grid.velarr)-1
			if (ymin lt fieldstate.flxmin) then ymin = fieldstate.flxmin
			if (ymax gt fieldstate.flxmax) then ymax = fieldstate.flxmax
				
			widget_control, fieldstate.specwin, Clear_Events=1
			fieldstate.chmin  = xmin
			fieldstate.chmax  = xmax
			fieldstate.flxmin = ymin
			fieldstate.flxmax = ymax

			fieldview_updatespec
		endif
	end
	'MOTION': begin
		if event.press eq 0 AND fieldstate.mousestate eq 1 then begin
			; Here is where the actual box is drawn and erased.
			; First, erase the last box.

			widget_control, fieldstate.specwin, get_value=index
			wset, index
			Device, Copy=[0, 0, 575, 220, 0, 0, 19]

			; Get the coodinates of the new box and draw it.
			sx = fieldstate.cornerx
			sy = fieldstate.cornery
			dx = event.x
			dy = event.y
			width=abs(dx-sx)

			color=FSC_Color('Green')
			
			PlotS, [sx, sx, dx, dx, sx], [sy, dy, dy, sy, sy], /Device, $
				Color=color, thick=1.5

				xyouts, sx+10, sy-10, 'Zoom Box', /device, color=color
		endif
	end
endcase

widget_control, fieldstate.specwin, get_value=index
wset, index

xdevice = event.x
ydevice = event.y

hor, fieldstate.chmin, fieldstate.chmax
ver, fieldstate.flxmin, fieldstate.flxmax

result=convert_coord(xdevice, ydevice, /Device, /Double, /To_Data)
xdata=result[0]
ydata=result[1]

if (xdata lt fieldstate.chmin OR xdata gt fieldstate.chmax OR $
	ydata lt fieldstate.flxmin OR ydata gt fieldstate.flxmax) then begin
	widget_control, fieldstate.cpix, set_value='----'
	widget_control, fieldstate.cval, set_value='------'
	widget_control, fieldstate.fval, set_value='------'
endif else begin
	cpix = round(xdata)
	cval = grid.velarr[cpix]
	flux = fieldstate.spec[cpix]
	widget_control, fieldstate.cpix, set_value=string(cpix,Format='(I4)')
	widget_control, fieldstate.cval, set_value=string(cval,Format='(F7.1)')
	widget_control, fieldstate.fval, set_value=string(flux,Format='(F7.2)')
endelse

end



pro fieldview_event, event

common fieldview_state
common cleanview_state
common gridstate
common clean

widget_control, event.id, get_uvalue=uvalue, get_value=value

case uvalue of
	'blue': begin
		if fieldstate.colortab NE 1 then begin
			fieldstate.colortab = 1
			widget_control, fieldstate.colorBLU, set_button=1
			widget_control, fieldstate.colorGRN, set_button=0
			widget_control, fieldstate.colorRED, set_button=0
			widget_control, fieldstate.colorRAI, set_button=0
			fieldview_updatemap
		endif
	end
	'green': begin
		if fieldstate.colortab NE 8 then begin
			fieldstate.colortab = 8
			widget_control, fieldstate.colorBLU, set_button=0
			widget_control, fieldstate.colorGRN, set_button=1
			widget_control, fieldstate.colorRED, set_button=0
			widget_control, fieldstate.colorRAI, set_button=0
			fieldview_updatemap
		endif	
	end
	'red': begin
		if fieldstate.colortab NE 3 then begin
			fieldstate.colortab = 3
			widget_control, fieldstate.colorBLU, set_button=0
			widget_control, fieldstate.colorGRN, set_button=0
			widget_control, fieldstate.colorRED, set_button=1
			widget_control, fieldstate.colorRAI, set_button=0
			fieldview_updatemap
		endif	
	end
	'rainbow': begin
		if fieldstate.colortab NE 39 then begin
			fieldstate.colortab = 39
			widget_control, fieldstate.colorBLU, set_button=0
			widget_control, fieldstate.colorGRN, set_button=0
			widget_control, fieldstate.colorRED, set_button=0
			widget_control, fieldstate.colorRAI, set_button=1
			fieldview_updatemap
		endif	
	end
	'invert': begin
		if fieldstate.colorinvert NE 1 then begin
			fieldstate.colorinvert = 1
			widget_control,fieldstate.colorINV, set_button=1
		endif else begin
			fieldstate.colorinvert = 0
			widget_control,fieldstate.colorINV, set_button=0
		endelse
		fieldview_updatemap
	end

	'show': begin
		if fieldstate.overp NE 1 then begin
			fieldstate.overp = 1
			widget_control,fieldstate.showct, set_button=1
		endif else begin
			fieldstate.overp = 0
			widget_control,fieldstate.showct, set_button=0
		endelse
		fieldview_updatemap
	end
	'level05': begin
		if fieldstate.levels NE  5 then begin
			fieldstate.levels = 5
			widget_control,fieldstate.cont05, set_button=1
			widget_control,fieldstate.cont10, set_button=0
			widget_control,fieldstate.cont15, set_button=0
			fieldview_updatemap
		endif
	end
	'level10': begin
		if fieldstate.levels NE 10 then begin
			fieldstate.levels = 10
			widget_control,fieldstate.cont05, set_button=0
			widget_control,fieldstate.cont10, set_button=1
			widget_control,fieldstate.cont15, set_button=0
			fieldview_updatemap
		endif
	end
	'level15': begin
		if fieldstate.levels NE 15 then begin
			fieldstate.levels = 15
			widget_control,fieldstate.cont05, set_button=0
			widget_control,fieldstate.cont10, set_button=0
			widget_control,fieldstate.cont15, set_button=1
			fieldview_updatemap
		endif
	end

	'clip518': begin
		if fieldstate.cliplevel NE 5d18 then begin
			fieldstate.cliplevel = 5d18
			widget_control,fieldstate.clip518, set_button=1
			widget_control,fieldstate.clip119, set_button=0
			widget_control,fieldstate.clip519, set_button=0
			widget_control,fieldstate.clip120, set_button=0	
			fieldview_updatemap, /reset
			fieldview_updatemap
		endif
	end
	'clip119': begin
		if fieldstate.cliplevel NE 1d19 then begin
			fieldstate.cliplevel = 1d19
			widget_control,fieldstate.clip518, set_button=0
			widget_control,fieldstate.clip119, set_button=1
			widget_control,fieldstate.clip519, set_button=0
			widget_control,fieldstate.clip120, set_button=0
			fieldview_updatemap, /reset
			fieldview_updatemap
		endif
	end
	'clip519': begin
		if fieldstate.cliplevel NE 5d19 then begin
			fieldstate.cliplevel = 5d19
			widget_control,fieldstate.clip518, set_button=0
			widget_control,fieldstate.clip119, set_button=0
			widget_control,fieldstate.clip519, set_button=1
			widget_control,fieldstate.clip120, set_button=0
			fieldview_updatemap, /reset
			fieldview_updatemap
		endif
	end
	'clip120': begin
		if fieldstate.cliplevel NE 1d20 then begin
			fieldstate.cliplevel = 1d20
			widget_control,fieldstate.clip518, set_button=0
			widget_control,fieldstate.clip119, set_button=0
			widget_control,fieldstate.clip519, set_button=0
			widget_control,fieldstate.clip120, set_button=1
			fieldview_updatemap, /reset
			fieldview_updatemap
		endif
	end

	'rescale': begin
		widget_control,fieldstate.xmin, get_value=xmin
		widget_control,fieldstate.xmax, get_value=xmax
		widget_control,fieldstate.ymin, get_value=ymin
		widget_control,fieldstate.ymax, get_value=ymax

		fieldstate.chmin  = double(xmin)
		fieldstate.chmax  = double(xmax)
		fieldstate.flxmin = double(ymin)
		fieldstate.flxmax = double(ymax)

		fieldview_updatespec
	end

	'reset': begin
		fieldstate.first = fieldstate.first * 0.0
		fieldstate.crange = [-1.0, -1.0]
		fieldstate.clickcount = 0
		
		fieldview_updatemap, /Reset
		fieldview_updatespec, /Reset, /Init
	end

	'comp': begin
		if fieldstate.crange[0] EQ 0 OR fieldstate.crange[1] EQ 0 then begin
			break
		endif
		; Flip over if needed
		fieldstate.crange = fieldstate.crange[ sort(fieldstate.crange) ]
		fieldstate.vmin = grid.velarr[fieldstate.crange[1]]
		fieldstate.vmax = grid.velarr[fieldstate.crange[0]]

		zeroth = total( cleand[fieldstate.crange[0]:fieldstate.crange[1],*,*], 1 )
		  zeroth *= 3.40895e16 * (mean(grid.velarr[fieldstate.crange])/2.99792458d5 + 1.0)^2.0
		  zeroth *= (max(grid.velarr)-min(grid.velarr))/n_elements(grid.velarr)
		temp = cleand*0.0
		for i=0,(n_elements(grid.velarr)-1) do begin
			if i GE fieldstate.crange[0] AND i LE fieldstate.crange[1] then begin
				temp[i,*,*] = cleand[i,*,*] * grid.velarr[i]
			endif
		endfor
		  temp *= 3.40895e16 * (mean(grid.velarr[fieldstate.crange])/2.99792458d5 + 1.0)^2.0
		  temp *= (max(grid.velarr)-min(grid.velarr))/n_elements(grid.velarr)
		first = total(temp,1) / zeroth

		xmin = fieldstate.llx - cleanstate.llx
		ymin = fieldstate.lly - cleanstate.lly
		xmax = fieldstate.urx - cleanstate.llx
		ymax = fieldstate.ury - cleanstate.lly
		fieldstate.zeroth = zeroth[xmin:xmax, ymin:ymax]
		fieldstate.first =   first[xmin:xmax, ymin:ymax]

		fieldview_updatemap
	end

	'about': begin
		h = ['FieldView - Velocity field utility        ', $
		     '  Written, December 2007', $
		     ' ', $
		     '  Last update, Tuesday, March 18, 2008']
		
		if NOT xregistered('fieldview_help', /noshow) then begin
			about_base =  widget_base(group_leader=fieldstate.baseID, title='About FieldView', $
				/Column, /Base_Align_Right, uvalue = 'about_base')
			about_text = widget_text(about_base, /Scroll, value = h, xsize = 45, ysize = 10)
			about_done = widget_button(about_base, value = ' Done ', uvalue = 'about_done', $
				event_pro='fieldview_event')
			
			widget_control, about_base, /realize
			xmanager, 'fieldview_help', about_base, /no_block
		endif
	end
	'about_done': begin
		widget_control, event.top, /destroy
	end

	'close': begin
		widget_control, event.top, /destroy

		hor
		ver
		loadct, 1, /silent
		
		delvarx, fieldstate
		
		widget_control, viewstate.baseID, /sensitive

		print, 'Exiting FieldView...'
	end
	else:
endcase

end


pro fieldview_updatespec, Reset=Reset, Init=Init

common fieldview_state
common cleanview_state
common gridstate
common clean

widget_control, fieldstate.specwin, get_value=index
wset,index

if Keyword_Set(Reset) then begin
	erase

	xmin = fieldstate.llx - cleanstate.llx
	ymin = fieldstate.lly - cleanstate.lly
	xmax = fieldstate.urx - cleanstate.llx
	ymax = fieldstate.ury - cleanstate.lly

	fieldstate.spec = cleand[*, (xmax+xmin)/2, (ymax+ymin)/2]
	if Keyword_Set(Init) then begin
		fieldstate.chmin = 0L
		fieldstate.chmax = 1023L
		fieldstate.flxmin = -10.0
		fieldstate.flxmax = 30.0
	endif

	xrange = [fieldstate.chmin, fieldstate.chmax]
	yrange = [fieldstate.flxmin, fieldstate.flxmax]	
	plot, [0,0], /NoData, Position=[0.15,0.15,0.95,0.95], $
		XRange=XRange, XStyle=1, YRange=YRange, YStyle=1, $
		XTitle='Channel Number', YTitle='Flux Density [mJy/beam]'

	fieldstate.clickcount = 0
	fieldstate.mousestate = 0

	fieldstate.xs = !X.s
	fieldstate.ys = !Y.s
	fieldstate.xregion = !X.region
	fieldstate.yregion = !Y.region
endif

xrange = [fieldstate.chmin, fieldstate.chmax]
yrange = [fieldstate.flxmin, fieldstate.flxmax]	
xplot = lindgen(fieldstate.chmax-fieldstate.chmin+1) + fieldstate.chmin
plot, xplot, fieldstate.spec[fieldstate.chmin:fieldstate.chmax], $
	XRange=XRange, XStyle=1, YRange=YRange, YStyle=1, $
	XTitle='Channel Number', YTitle='Flux Density [mJy/beam]', $
	Position=[0.15,0.15,0.95,0.95]

if fieldstate.crange[0] NE -1.0 then $
	flag,fieldstate.crange[0],color=FSC_Color('Red')
if fieldstate.crange[1] NE -1.0 then $
	flag,fieldstate.crange[1],color=FSC_Color('Red')

Window, 20, /Pixmap, XSize=575, YSize=220
wset, 20
Device, Copy=[0, 0, 575, 220, 0, 0, index]

!X.s = fieldstate.xs
!Y.s =fieldstate.ys
!X.region = fieldstate.xregion
!Y.region = fieldstate.yregion

end



pro fieldview_updatemap, Reset=Reset

common fieldview_state
common cleanview_state
common gridstate
common clean

loadct,fieldstate.colortab,/Silent

widget_control, fieldstate.fieldwin, get_value=index
wset, index

if Keyword_Set(Reset) then begin
	erase

	xrange = [fieldstate.ramax, fieldstate.ramin]
	yrange = [fieldstate.decmin, fieldstate.decmax]
	plot, [0,0], /nodata, xrange=xrange, yrange=yrange, xstyle=1, ystyle=1, $
		xtick_get=xvals, ytick_get=yvals, ticklen=ticklen, charsize=1.0, $
		Position=[0.20,0.20,0.95,0.95]
	nxticklabels=n_elements(xvals)
	nyticklabels=n_elements(yvals)
	xspacing=((xvals[n_elements(xvals)-1]-xvals[0])*60.0)/(nxticklabels-1)
	yspacing=((yvals[n_elements(yvals)-1]-yvals[0])*60.0)/(nyticklabels-1)
	xticlabs=ratickname(xvals*15.0)
	yticlabs=dectickname(yvals)
	plot, [0,0], /nodata, xrange=reverse(xrange), yrange=yrange, $
		xtitle='Right Ascension (J2000)', ytitle='Declination (J2000)', $
		xstyle=1, ystyle=1, Position=[0.20,0.20,0.95,0.95], $
		xtickn=reverse(xticlabs), ytickn=yticlabs, ticklen=-0.01, charsize=1.0
	fieldstate.px = !X.WINDOW * !D.X_VSIZE
	fieldstate.py = !Y.WINDOW * !D.Y_VSIZE
	fieldstate.sx = fieldstate.px[1] - fieldstate.px[0] - 1
	fieldstate.sy = fieldstate.py[1] - fieldstate.py[0] - 1

	widget_control, fieldstate.barwin, get_value=index
	wset, index
	
	erase
	if fieldstate.colorinvert EQ 1 then begin
		colorbar, Range=[fieldstate.vmin, fieldstate.vmax], Title='Velocity [km/s]', $
			/InvertColors, /Ver, Pos=[0.55,0.10,0.90,0.90]
	endif else begin
		colorbar, Range=[fieldstate.vmin, fieldstate.vmax], Title='Velocity [km/s]', $
			/Ver, Pos=[0.55,0.10,0.90,0.90]
	endelse
endif else begin
	c_levels = findgen(fieldstate.levels)*(fieldstate.vmax-fieldstate.vmin)/(fieldstate.levels-1.0) + fieldstate.vmin

	first = fieldstate.first
	bad = where( finite(first) EQ 0 OR fieldstate.zeroth LT fieldstate.cliplevel )
	  if bad[0] NE -1 then first[bad] = -3000.0

	if fieldstate.colorinvert EQ 1 then begin
		tvlct, r, g, b, /Get
		tvlct, Reverse(r), Reverse(g), Reverse(b)
	endif
	contour, reverse(first), levels=c_levels, c_labels=[0], pos=[0.20,0.20,0.95,0.95], $
		XRange=[0,(size(first))[1]-1],YRange=[0,(size(first))[2]-1], $
		XStyle=5,YStyle=5, /Fill, /NoErase
	if fieldstate.colorinvert EQ 1 then begin
		tvlct, r, g, b, /Get
		tvlct, Reverse(r), Reverse(g), Reverse(b)
	endif

	if fieldstate.overp then begin
		contour, reverse(first), levels=c_levels, c_labels=[0], pos=[0.20,0.20,0.95,0.95], $
			XRange=[0,(size(first))[1]-1],YRange=[0,(size(first))[2]-1], $
			XStyle=5,YStyle=5, /Follow, /Overplot
	endif

	widget_control, fieldstate.barwin, get_value=index
	wset, index

	erase
	if fieldstate.colorinvert EQ 1 then begin
		colorbar, Range=[fieldstate.vmin, fieldstate.vmax], Title='Velocity [km/s]', $
			/InvertColors, /Ver, Pos=[0.55,0.10,0.90,0.90]
	endif else begin
		colorbar, Range=[fieldstate.vmin, fieldstate.vmax], Title='Velocity [km/s]', $
			/Ver, Pos=[0.55,0.10,0.90,0.90]
	endelse
endelse

fieldview_updatespec

end



pro fieldview_jpeg, event

common fieldview_state

if NOT xregistered('fieldview_jpeg', /NoShow) then begin
	jpeg_output_base = widget_base(group_leader=fieldstate.baseID, /row, /base_align_right, $
		title='Velocity FIELDview JPEG Output', uvalue = 'jpeg_output_base')
	
	cd, current=pwd  ;Get current working directory into a string
	cra = (fieldstate.ramax + fieldstate.ramin) / 2.0 	; RA  of box center
	 crah = string(floor(cra), Format='(I02)')
	 cram = string(cra*60.0   mod 60, Format='(I02)')
	 cras = string(cra*3600.0 mod 60, Format='(F04.1)')
	cdc = (fieldstate.decmax + fieldstate.decmin) / 2.0	; Dec of box center
	 if cdc LT 0.0 then cdce = '-' else cdce = '+'
	 cdcd = string(floor(cdc), Format='(I02)')
	 cdcm = string(cdc*60.0   mod 60, Format='(I02)')
	 cdcs = string(cdc*3600.0 mod 60, Format='(I02)')
	fname = 'vfield_'+crah+cram+cras+cdce+cdcd+cdcm+cdcs+'.jpg'	; File name suggestion

	  filesbase=widget_base(jpeg_output_base, /column, /align_left)
	    directorybase=widget_base(filesbase, /Column)
	    label=widget_label(directorybase, value='Output Directory:  ')
	    fieldstate.dir=widget_text(directorybase, xsize=40, value=pwd+'/', /Editable, uvalue='dir')
	    label=widget_label(directorybase, value='Output File Name:  ')
	    fieldstate.nam=widget_text(directorybase, xsize=40, value=fname, /Editable, uvalue='nam')
	  buttonbase=widget_base(filesbase, xsize=100,ysize=100, /align_right, /column)
	    jpeg_output_export = widget_button(buttonbase, value = ' Export ', $
		uvalue = 'export', event_pro='fieldview_jpeg_event')
	    cancel=widget_button(buttonbase, value=' Cancel ', uvalue='cancel', $
					event_pro='fieldview_jpeg_event')
	
	widget_control, jpeg_output_base, /realize
	xmanager, 'fieldview_jpeg', jpeg_output_base, /no_block
endif

end


pro fieldview_jpeg_event, event

common fieldview_state
common cleanview_state
common gridstate
common clean

widget_control, event.id, get_uvalue = uvalue

case uvalue of 
	'export': begin
		widget_control, fieldstate.dir, get_value=dir
		widget_control, fieldstate.nam, get_value=filename

		widget_control, /Hourglass

		filename = dir+filename

		thisDevice = !D.Name
   		set_plot, 'Z', /Copy
		device, set_resolution=[700,800], Z_Buffer=0
		erase

		xrange = [fieldstate.ramax, fieldstate.ramin]
		yrange = [fieldstate.decmin, fieldstate.decmax]
		plot, [0,0], /nodata, xrange=xrange, yrange=yrange, xstyle=1, ystyle=1, $
			xtick_get=xvals, ytick_get=yvals, ticklen=ticklen, charsize=1.0, $
			Position=[0.20,0.20,0.95,0.95]
		nxticklabels=n_elements(xvals)
		nyticklabels=n_elements(yvals)
		xspacing=((xvals[n_elements(xvals)-1]-xvals[0])*60.0)/(nxticklabels-1)
		yspacing=((yvals[n_elements(yvals)-1]-yvals[0])*60.0)/(nyticklabels-1)
		xticlabs=ratickname(xvals*15.0)
		yticlabs=dectickname(yvals)
		plot, [0,0], /nodata, xrange=reverse(xrange), yrange=yrange, $
			xtitle='Right Ascension (J2000)', ytitle='Declination (J2000)', $
			xstyle=1, ystyle=1, Position=[0.20,0.20,0.95,0.95], $
			xtickn=reverse(xticlabs), ytickn=yticlabs, ticklen=-0.01, charsize=1.0
; 		px = !X.WINDOW * !D.X_VSIZE
; 		py = !Y.WINDOW * !D.Y_VSIZE
; 		sx = px[1] - px[0] - 1
; 		sy = py[1] - py[0] - 1
		
		c_levels = findgen(fieldstate.levels)*(fieldstate.vmax-fieldstate.vmin)/(fieldstate.levels-1.0) + fieldstate.vmin

		first = fieldstate.first
		bad = where( finite(first) EQ 0 OR fieldstate.zeroth LT fieldstate.cliplevel )
		if bad[0] NE -1 then first[bad] = -3000.0
	
		if fieldstate.colorinvert EQ 1 then begin
			tvlct, r, g, b, /Get
			tvlct, Reverse(r), Reverse(g), Reverse(b)
		endif
		contour, reverse(first), levels=c_levels, c_labels=[0], pos=[0.20,0.20,0.95,0.95], $
			XRange=[0,(size(first))[1]-1],YRange=[0,(size(first))[2]-1], $
			XStyle=5,YStyle=5, /Fill, /NoErase
		if fieldstate.colorinvert EQ 1 then begin
			tvlct, r, g, b, /Get
			tvlct, Reverse(r), Reverse(g), Reverse(b)
		endif
	
		if fieldstate.overp then begin
			contour, reverse(first), levels=c_levels, c_labels=[0], pos=[0.20,0.20,0.95,0.95], $
				XRange=[0,(size(first))[1]-1],YRange=[0,(size(first))[2]-1], $
				XStyle=5,YStyle=5, /Follow, /Overplot
		endif
	
		if fieldstate.colorinvert EQ 1 then begin
			colorbar, Range=[fieldstate.vmin, fieldstate.vmax], Title='Velocity [km/s]', $
				/InvertColors, Pos=[0.10,0.05,0.90,0.10]
		endif else begin
			colorbar, Range=[fieldstate.vmin, fieldstate.vmax], Title='Velocity [km/s]', $
				Pos=[0.10,0.05,0.90,0.10]
		endelse
		
		snapshot = tvrd()
   		tvlct, r, g, b, /Get
  		device, Z_Buffer=1
  		set_plot, thisDevice
		
		image24 = BytArr(3, 700, 800)
		image24[0,*,*] = r[snapshot]
		image24[1,*,*] = g[snapshot]
		image24[2,*,*] = b[snapshot]
		write_jpeg, filename, image24, true=1, quality=100

		widget_control, event.top, /destroy
	end

	'cancel': widget_control, event.top, /destroy
	else: 
endcase

end
