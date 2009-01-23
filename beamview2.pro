pro beamview

; Information we need from cleanview and gridview
common cleanview_state
common gridstate
common clean

;t_center = [cleanstate.cx, cleanstate.cy]
t_center = [(cleanstate.urx+cleanstate.llx)/2, (cleanstate.ury+cleanstate.lly)/2]
t_width = [cleanstate.urx-cleanstate.llx+1, cleanstate.ury-cleanstate.lly+1]+2*cleanstate.map_buffer

beams = cleanstate.beam
log = cleanstate.beamlog
flux_factor = cleanstate.area_corr

fitl = (size(beams))[3]/2-6
fitu = (size(beams))[3]/2+6

areas = total(total(beams,3),3)
out_posi = strarr(n_elements(areas))
out_posw = strarr(n_elements(areas))
out_beamf = strarr(n_elements(areas))
out_beamp = strarr(n_elements(areas))
out_beamc = strarr(n_elements(areas))
for i=0,(n_elements(areas)-1) do begin
	ubeam_x = i mod (size(areas))[1]
	ubeam_y = i / (size(areas))[1]

	c_beam = reform( beams[ubeam_x, ubeam_y, *, *] )
	xc = t_center[0] - t_width[0]/2 + ubeam_x
	yc = t_center[1] - t_width[1]/2 + ubeam_y
	; String to display for the location of the central beam
	out_posi[i] = '('+strtrim(string(xc,Format='(I3)'),2)+','+ $
		strtrim(string(yc,Format='(I3)'),2)+')'
	; String to display for the location of the central beam
	rac = grid.ramin + xc*grid.deltara/3600.0
		rac_h = fix(rac)
		rac_m = fix( (rac - rac_h)*60.0 )
		rac_s = rac*3600.0 - rac_h*3600.0 - rac_m*60.0
		if rac_s GE 60 then begin
			rac_s = rac_s - 60.0
			rac_m = rac_m + 1
		endif
		rac = string(rac_h,Format='(I02)')+':'+string(rac_m,Format='(I02)')+ $
			':'+string(rac_s,Format='(F05.2)')
	decc = grid.decmin + yc*grid.deltadec/60.0
		decc_d = fix(abs(decc))
		decc_m = fix( (abs(decc) - decc_d)*60.0 )
		decc_s = abs(decc)*3600.0 - decc_d*3600.0 - decc_m*60.0
		if decc_s GE 60 then begin
			decc_s = decc_s - 60.0
			decc_m = decc_m + 1
		endif
		if decc LT 0 then decc_d = -decc_d
		decc = string(decc_d,Format='(I+03)')+':'+ $
			string(decc_m,Format='(I02)')+':'+ $
			string(decc_s,Format='(F04.1)')
	out_posw[i] = '('+rac+', '+decc+')'
	; Now figure out beam FWHM and PA
	junk = gauss2dfit( c_beam[fitl:fitu, fitl:fitu], beam_fit, /Tilt)
	;+ Load up the results
	beam_fwhm = 2.0*sqrt(2.0*alog(2.0))*[beam_fit[3], beam_fit[2]]
	beam_pa = beam_fit[6]*!radeg
	out_beamf[i] = string(beam_fwhm[0],Format='(F4.2)')+"' x "+$
		string(beam_fwhm[1],Format='(F4.2)')+"'"
	out_beamp[i] = string(beam_pa,Format='(F+6.1)')+' degrees'	
	beam_corr = flux_factor[ubeam_x, ubeam_y]
	out_beamc[i] = string(1.0/beam_corr, Format='(F8.6)')
endfor

; Information we need only for beamview
common beamview_state, beamstate

; Define beamstate (beamview window)
beamstate = {baseID: 0L, beam_win: 0L, map_win: 0L, pos_gt: 0L, pos_st: 0L, fwhmt: 0L, pat: 0L, corrt: 0L, $
	xval: 0L, yval: 0L, dval: 0L, $
	cbIndex: 0, scaleLIN: 0L, scaleLOG: 0L, UseLog: 1, $
	colorBLU: 0L, colorGRN: 0L, colorRED: 0L, colorRAI: 0L, colorINV: 0L, colortab: 39, colorinvert: 0, $
	dir: 0L, fnm: 0L, $
	beams: beams, fluxfact: flux_factor, pos_g: out_posi, pos_s: out_posw, fwhm: out_beamf, pa: out_beamp, corr: out_beamc }

; Now for the window
beam_main_base = widget_base(group_leader=viewstate.baseID, Title='BEAMview', mbar=top_menu, $
	/Column, /Base_Align_Left)
beamstate.baseID = beam_main_base
filemenu = widget_button(top_menu, value=' File ')
	jpegbutton = widget_button(filemenu, value=' Export JPEG ', uvalue='jpeg', event_pro='beamview_jpeg')
	exitbutton = widget_button(filemenu, value=' Exit ', uvalue='close', /Separator)
colormenu = widget_button(top_menu, value=' Color ')
	beamstate.colorBLU = widget_button(colormenu, value=' Blue/White ',  uvalue='blue', /Checked_Menu)
	beamstate.colorGRN = widget_button(colormenu, value=' Green/White ', uvalue='green', /Checked_Menu)
	beamstate.colorRED = widget_button(colormenu, value=' Red/White ',   uvalue='red', /Checked_Menu)
	beamstate.colorRAI = widget_button(colormenu, value=' Rainbow ',     uvalue='rainbow', /Checked_Menu)
	beamstate.colorINV = widget_button(colormenu, value='Invert Colors', uvalue='invert', /Checked_Menu, /Separator)
scalingmenu = widget_button(top_menu, value=' Scaling ')
	beamstate.scaleLIN = widget_button(scalingmenu, value=' Linear ', uvalue = 'linear', /Checked_Menu)
	beamstate.scaleLOG = widget_button(scalingmenu, value=' Logarithmic ', uvalue = 'loga', /Checked_Menu)
helpmenu = widget_button(top_menu, Value=' Help ')
	aboutbutton = widget_button(helpmenu, value=' About ', uvalue='about')
	
beam_upper_base = widget_base(beam_main_base, /Row)
	beam_upper_sub1 = widget_base(beam_upper_base, /Column)
		label = widget_label(beam_upper_sub1, value='Effective Beam Pattern:')
		beamstate.beam_win = widget_draw(beam_upper_sub1, XSize=400, YSize=400, Frame=1, $
			Event_Pro='beamview_display_event', /Motion_Events, /Button_Events, /Align_Center)

		beam_upper_sub11 = widget_base(beam_upper_sub1, /Column, /Frame)
			label = widget_label(beam_upper_sub11, value='Detail at Cursor:', /Align_Center)
			left_sub1 = widget_base(beam_upper_sub11, /Row, /Align_Center)
			label = widget_label(left_sub1, value='Delta X: ')
			beamstate.xval = widget_label(left_sub1, value='----')
			label = widget_label(left_sub1, value='  ')
			label = widget_label(left_sub1, value='Delta Y: ')
			beamstate.yval = widget_label(left_sub1, value='----')
			label = widget_label(left_sub1, value='  ')
			label = widget_label(left_sub1, value='Strength [dB]: ')
			beamstate.dval = widget_label(left_sub1, value='----  ')

		label = widget_label(beam_upper_sub1, value='Beam Modeling Logs: ', /Align_Center)
		label = widget_text(beam_upper_sub1, value=log, /Scroll, YSize=10, XSize=60)

	beam_upper_sub2 = widget_base(beam_upper_base, /Column)
		label = widget_label(beam_upper_sub2, value='Unique Beam Regions:')
		beamstate.map_win = widget_draw(beam_upper_sub2, XSize=400, YSize=400, Frame=1, $
			Event_Pro='beamview_map_event', /Motion_Events, /Button_Events, /Align_Center)

		beam_upper_sub21 = widget_base(beam_upper_sub2, /Column, /Frame)
			label = widget_label(beam_upper_sub21, value='Beam Properties:', /Align_Center)	
			right_sub0 = widget_base(beam_upper_sub21, /Row, /Align_Left)
				label = widget_label(right_sub0, value=' ')
			right_sub1 = widget_base(beam_upper_sub21, /Row, /Align_Left)
				label = widget_label(right_sub1, value='Center on Grid [px]:')
				beamstate.pos_gt = widget_label(right_sub1, value=out_posi[n_elements(areas)/2]+'  ',  /Align_Left)
			right_sub2 = widget_base(beam_upper_sub21, /Row, /Align_Left)
				label = widget_label(right_sub2, value='Center on Sky [RA, Dec.]:')
				beamstate.pos_st = widget_label(right_sub2, value=out_posw[n_elements(areas)/2], /Align_Left)
			right_sub3 = widget_base(beam_upper_sub21, /Row, /Align_Left)
				label = widget_label(right_sub3, value='FWHM:')
				beamstate.fwhmt = widget_label(right_sub3, value=out_beamf[n_elements(areas)/2], /Align_Left)
			right_sub4 = widget_base(beam_upper_sub21, /Row, /Align_Left)
				label = widget_label(right_sub4, value='PA:')
				beamstate.pat = widget_label(right_sub4, value=out_beamp[n_elements(areas)/2], /Align_Left)
			right_sub5 = widget_base(beam_upper_sub21, /Row, /Align_Left)
				label = widget_label(right_sub5, value='Area Correction:')
				beamstate.corrt = widget_label(right_sub5, value=out_beamc[n_elements(areas)/2], /Align_Left)
			right_sub6 = widget_base(beam_upper_sub21, /Row, /Align_Left)
				label = widget_label(right_sub6, value=' ')
			right_sub7 = widget_base(beam_upper_sub21, /Row, /Align_Left)
				label = widget_label(right_sub7, value='Note: The white box is the field-of-view shown in CLEANview')

widget_control, beam_main_base, /realize
xmanager, 'beamview', beam_main_base, /No_Block

widget_control, beamstate.colorRAI, Set_Button=1
beamstate.colortab = 39
widget_control, beamstate.scaleLOG, Set_Button=1
beamstate.UseLog = 1

beamstate.cbIndex = n_elements(areas)/2
beamview_update, /Init
end



pro beamview_display_event, event

common beamview_state
common clean

widget_control, beamstate.beam_win, get_value=index
wset, index

tempx = beamstate.cbIndex mod (size(beamstate.fluxfact))[1]
tempy = beamstate.cbIndex / (size(beamstate.fluxfact))[1]
beam = beamstate.beams[tempx, tempy, *, *]
beam = reform(beam)

xdevice=event.x
ydevice=event.y
xdata=round( -(xdevice/400.0 - 0.15)/0.8* n_elements(beam[*,0]) + n_elements(beam[*,0])/2.0 )
ydata=round( (ydevice/400.0 - 0.15)/0.8* n_elements(beam[0,*]) - n_elements(beam[0,*])/2.0 )

if (xdata lt -(n_elements(beam[*,0]))/2 OR xdata gt (n_elements(beam[*,0]))/2 OR $
	ydata lt -(n_elements(beam[0,*]))/2 OR ydata gt (n_elements(beam[0,*]))/2) then begin
	widget_control, beamstate.xval, set_value='----'
	widget_control, beamstate.yval, set_value='----'
	widget_control, beamstate.dval, set_value='----'
endif else begin
	old_except = !Except
	!Except = 0
	xpix = round( xdata + n_elements(beam[*,0])/2 )
	ypix = round( ydata + n_elements(beam[0,*])/2 )
	dval = (10.0d * (alog10(beam) -  alog10(max(beam))))[xpix, ypix]
	widget_control, beamstate.xval, set_value=string(xdata,Format='(I3)')
	widget_control, beamstate.yval, set_value=string(ydata,Format='(I3)')
	if dval GE -50.0 then begin
		widget_control, beamstate.dval, set_value=string(dval,Format='(F+6.2)')
	endif else begin
		widget_control, beamstate.dval, set_value='<-50.'
	endelse
	junk=check_math()
	!Except = old_except
endelse
end



pro beamview_map_event, event

common beamview_state
common clean

widget_control, beamstate.map_win, get_value=index
wset, index

Event_Types = ['DOWN', 'UP', 'MOTION', '?', '?', '?', '?']
This_Event = Event_Types[event.type]
case This_Event of
	'DOWN'  : begin
		xdevice=event.x
		ydevice=event.y
		xdata = floor( -(xdevice/400.0 - 0.15)/0.8* n_elements(cleanstate.area_corr[*,0]) + n_elements(cleanstate.area_corr[*,0]))
		ydata = floor( (ydevice/400.0 - 0.15)/0.8* n_elements(cleanstate.area_corr[0,*]) )
		
		if (xdata GE 0 AND xdata LT (size(beamstate.fluxfact))[1] AND ydata GE 0 AND ydata LT (size(beamstate.fluxfact))[2]) then begin
			new_region = xdata+(size(beamstate.fluxfact))[1]*ydata
			; Since the uniques of regions is decided based on the flux correction, only update if the flux correction where 
			; we have clicked is different than the flux correction of where we were previously at.
			beamstate.cbIndex = new_region
				
			widget_control, beamstate.pos_gt, set_value=((beamstate.pos_g)[new_region])[0]
			widget_control, beamstate.pos_st, set_value=((beamstate.pos_s)[new_region])[0]
			widget_control, beamstate.fwhmt, set_value=((beamstate.fwhm)[new_region])[0]
			widget_control, beamstate.pat, set_value=((beamstate.pa)[new_region])[0]
			widget_control, beamstate.corrt, set_value=((beamstate.corr)[new_region])[0]
				
			beamview_update
		endif
	end
	else:
endcase

end



pro beamview_event, event

common cleanview_state
common beamview_state

widget_control, event.id, get_uvalue=uvalue, get_value=value
to_change = 1
case uvalue of
	'linear': begin
		if beamstate.UseLog NE 0 then begin
			beamstate.UseLog = 0
			widget_control, beamstate.scaleLIN, Set_Button=1
			widget_control, beamstate.scaleLOG, Set_Button=0
			beamview_update
		endif
	end
	'loga': begin
		if beamstate.UseLog NE 1 then begin
			beamstate.UseLog = 1
			widget_control, beamstate.scaleLIN, Set_Button=0
			widget_control, beamstate.scaleLOG, Set_Button=1
			beamview_update
		endif
	end

	'blue': begin
		if beamstate.colortab NE 1 then begin
			beamstate.colortab = 1
			widget_control, beamstate.colorBLU, Set_Button=1
			widget_control, beamstate.colorGRN, Set_Button=0
			widget_control, beamstate.colorRED, Set_Button=0
			widget_control, beamstate.colorRAI, Set_Button=0
			beamview_update
		endif
	end
	'green': begin
		if beamstate.colortab NE 8 then begin
			beamstate.colortab = 8
			widget_control, beamstate.colorBLU, Set_Button=0
			widget_control, beamstate.colorGRN, Set_Button=1
			widget_control, beamstate.colorRED, Set_Button=0
			widget_control, beamstate.colorRAI, Set_Button=0
			beamview_update
		endif
	end
	'red': begin
		if beamstate.colortab NE 3 then begin
			beamstate.colortab = 3
			widget_control, beamstate.colorBLU, Set_Button=0
			widget_control, beamstate.colorGRN, Set_Button=0
			widget_control, beamstate.colorRED, Set_Button=1
			widget_control, beamstate.colorRAI, Set_Button=0
			beamview_update
		endif
	end
	'rainbow': begin
		if beamstate.colortab NE 39 then begin
			beamstate.colortab = 39
			widget_control, beamstate.colorBLU, Set_Button=0
			widget_control, beamstate.colorGRN, Set_Button=0
			widget_control, beamstate.colorRED, Set_Button=0
			widget_control, beamstate.colorRAI, Set_Button=1
			beamview_update
		endif
	end
	'invert': begin
		if beamstate.colorinvert EQ 1 then begin
			beamstate.colorinvert = 0
			widget_control, beamstate.colorINV, Set_Button=0
		endif else begin
			beamstate.colorinvert = 1
			widget_control, beamstate.colorINV, Set_Button=1
		endelse
		beamview_update
	end

	'about': begin
		h = ['CleanView - Graphical beam utility        ', $
		     '  Written, November 2007', $
		     ' ', $
		     '  Last update, Tuesday, March 18, 2008']
		
		if NOT xregistered('beamview_help', /noshow) then begin
			about_base =  widget_base(group_leader=beamstate.baseID, title='About BeamView', $
				/Column, /Base_Align_Right, uvalue = 'about_base')
			about_text = widget_text(about_base, /Scroll, value = h, xsize = 45, ysize = 10)
			about_done = widget_button(about_base, value = ' Done ', uvalue = 'about_done', $
				event_pro='beamview_event')
			
			widget_control, about_base, /realize
			xmanager, 'beamview_help', about_base, /no_block
		endif
	end
	'about_done': begin
		widget_control, event.top, /destroy
	end

	'close': begin
		widget_control, event.top, /destroy
		
		hor
		ver
		loadct,1,/Silent
		device, decomposed=0

		delvarx, beamstate

		widget_control, viewstate.baseID, /sensitive
		
		print, 'Exiting BeamView...'
	end
	else:
endcase
end



pro beamview_update, Init=Init

common cleanview_state
common beamview_state
common gridstate
common clean

if Keyword_Set(Init) then begin
	widget_control, beamstate.map_win, get_value=index
	wset,index

	loadct,39,/Silent
	temp = beamstate.fluxfact
	utemp = (temp[sort(temp)])[ uniq( temp[ sort(temp) ] ) ]
	for i=0,(n_elements(temp)-1) do $
		temp[i] = where( temp[i] EQ utemp )

	hor, viewstate.ramax+cleanstate.map_buffer*grid.deltara/3600.0, viewstate.ramin-cleanstate.map_buffer*grid.deltara/3600.0
	ver, viewstate.decmin-cleanstate.map_buffer*grid.deltadec/60.0, viewstate.decmax+cleanstate.map_buffer*grid.deltadec/60.0 

	plot, [0,0], /nodata, xstyle=1, ystyle=1, xtick_get=xvals, ytick_get=yvals, $
		ticklen=ticklen, charsize=1.0, Position=[0.15,0.15,0.95,0.95]
	nxticklabels=n_elements(xvals)
	nyticklabels=n_elements(yvals)
	xspacing=((xvals[n_elements(xvals)-1]-xvals[0])*60.0)/(nxticklabels-1)
	yspacing=((yvals[n_elements(yvals)-1]-yvals[0])*60.0)/(nyticklabels-1)
	xticlabs=ratickname(xvals*15.0)
	yticlabs=dectickname(yvals)
	
	erase

	to_disp = bytscl( congrid(reverse(temp), 320, 320), min = -1, max=n_elements(utemp) )
	tv,reverse(to_disp),60,60

	device, decomposed=1
	plot, [0,0], /nodata, /NoErase, xtitle='RA [HMS] J2000', ytitle='Dec [DMS] J2000', $
		xstyle=1, ystyle=1, Position=[0.15,0.15,0.95,0.95], xtickn=xticlabs, ytickn=yticlabs

	device, decomposed=0

	x1 = viewstate.ramax
	x2 = viewstate.ramin
	y1 = viewstate.decmin
	y2 = viewstate.decmax
	oplot,[x1, x2, x2, x1, x1], [y1, y1, y2, y2, y1], Color='0000FF'XL, Thick=2.0
endif

widget_control, beamstate.beam_win, get_value=index
wset,index

tempx = beamstate.cbIndex mod (size(beamstate.fluxfact))[1]
tempy = beamstate.cbIndex / (size(beamstate.fluxfact))[1]
beam = beamstate.beams[tempx, tempy, *, *]
beam = reform(beam)

beam_cont = beam / max(beam)
beam_cont = 10.0*alog10(beam_cont)

beam_min = -24.0
beam_max =  00.0
beam_disp = beam_cont
if beamstate.UseLog NE 1 then begin
	beam_min = 0.005
	beam_max = 0.90
	beam_disp = beam / max(beam)
endif

loadct,beamstate.colortab,/Silent
offset = 0B
if beamstate.colorinvert EQ 1 then begin
	tvlct, r, g, b, /Get
	tvlct, Reverse(r), Reverse(g), Reverse(b)
endif

tvimage, bytscl(reverse(beam_disp), min=beam_min, max=beam_max), Pos=[0.15,0.15,0.95,0.95], $
	/Minus_One

if beamstate.colorinvert EQ 1 then begin
	tvlct, r, g, b, /Get
	tvlct, Reverse(r), Reverse(g), Reverse(b)
endif

contour, reverse(beam_cont), levels=[-21, -18, -15, -12, -9, -6, -3], $
	c_labels=[1, 1, 1, 1, 1, 1, 1], Position=[0.15,0.15,0.95,0.95], $
	XTitle='X Offset [arc min]', YTitle='Y Offset [arc min]', $
	XRange=[0,(n_elements(beam_cont[*,0])-1)],YRange=[0,(n_elements(beam_cont[0,*])-1)],XStyle=5,YStyle=5, /NoErase

plot,[0,0], /NoData, XRange=[(n_elements(beam_cont[*,0])-1)/2, -(n_elements(beam_cont[*,0])-1)/2],XStyle=1, $
	YRange=[-(n_elements(beam_cont[0,*])-1)/2, (n_elements(beam_cont[0,*])-1)/2],YStyle=1, $
	XTitle='X Offset [arc min]', YTitle='Y Offset [arc min]', $
	Position=[0.15,0.15,0.95,0.95], /NoErase

end



pro beamview_jpeg, event

common cleanview_state
common beamview_state
common gridstate
common clean

if NOT xregistered('beamview_jpeg', /NoShow) then begin
	jpeg_output_base =  widget_base(group_leader=beamstate.baseID, /Column, /Base_Align_Right, $
		Title = 'BEAMview JPEG Output', uvalue = 'jpeg_output_base')
	
	filename = 'beam_'+strcompress(grid.name)+ $
			'_'+string((cleanstate.urx+cleanstate.llx)/2, Format='(I03)')+ $
			'_'+string((cleanstate.ury+cleanstate.lly)/2, Format='(I03)')+'.jpg'
	cd, current=pwd  ;Get current working directory into a string
	filesbase=widget_base(jpeg_output_base, /Column, /Align_Left)
	directorybase=widget_base(filesbase, /Column)
	label=widget_label(directorybase, value='Output Directory: ', /Align_Left)
	beamstate.dir=widget_text(directorybase, xsize=40, value=pwd+'/', /Editable, $
		uvalue='directory')
	label=widget_label(directorybase, value='Output Filename: ', /Align_Left)
	beamstate.fnm = widget_text(directorybase, xsize=40, value=filename, /Editable, $
		uvalue='filename')
	
	buttonbase=widget_base(filesbase, /Align_Right, /Row)
	jpeg_output_export = widget_button(buttonbase, value = ' Export ', $
		uvalue = 'export', event_pro='beamview_jpeg_event')
	cancel=widget_button(buttonbase, value=' Cancel ', uvalue='cancel', $
		event_pro='beamview_jpeg_event')
	
	widget_control, jpeg_output_base, /realize
	xmanager, 'beamview_jpeg', jpeg_output_base, /no_block 
endif
end



pro beamview_jpeg_event, event

common cleanview_state
common beamview_state
common gridstate
common clean

widget_control, event.id, get_uvalue = uvalue

case uvalue of
	'export': begin
		widget_control, beamstate.dir, get_value=dir
		widget_control, beamstate.fnm, get_value=filename
		widget_control, beamstate.beam_win, get_value=index

		widget_control, /Hourglass

		filename = dir+filename

		thisDevice = !D.Name
   		set_plot, 'Z', /Copy
		device, set_resolution=[700,700], Z_Buffer=0
		erase
		
		beam = beamstate.beams[(size(cleanstate.area_corr))[1]/2, (size(cleanstate.area_corr))[2]/2, *, *]
		beam = reform(beam)

		bindex = (size(cleanstate.area_corr))[1]/2 + (size(cleanstate.area_corr))[2]*(size(cleanstate.area_corr))[1]/2
		beam_posg = beamstate.pos_g[bindex]
		beam_poss = beamstate.pos_s[bindex]
		beam_fwhm = beamstate.fwhm[bindex]
		beam_pa   = beamstate.pa[bindex]
		beam_corr = beamstate.corr[bindex]

		beam_cont = beam / max(beam)
		beam_cont = 10.0*alog10(beam_cont)
		
		beam_min = -24.0
		beam_max =  00.0
		beam_disp = beam_cont
		if beamstate.UseLog NE 1 then begin
			beam_min = 0.005
			beam_max = 0.90
			beam_disp = beam / max(beam)
		endif
		
		loadct,beamstate.colortab,/Silent
		offset = 0B
		if beamstate.colorinvert EQ 1 then begin
			tvlct, r, g, b, /Get
			tvlct, Reverse(r), Reverse(g), Reverse(b)
		endif
		
		tvimage, bytscl(reverse(beam_disp), min=beam_min, max=beam_max), Pos=[0.15,0.15,0.95,0.95], $
			/Minus_One
		
		if beamstate.colorinvert EQ 1 then begin
			tvlct, r, g, b, /Get
			tvlct, Reverse(r), Reverse(g), Reverse(b)
		endif
		
		contour, reverse(beam_cont), levels=[-21, -18, -15, -12, -9, -6, -3], $
			c_labels=[1, 1, 1, 1, 1, 1, 1], Position=[0.15,0.15,0.95,0.95], $
			XTitle='X Offset [arc min]', YTitle='Y Offset [arc min]', $
			XRange=[0,(n_elements(beam_cont[*,0])-1)],YRange=[0,(n_elements(beam_cont[0,*])-1)], XStyle=5,YStyle=5, /NoErase
		
		plot,[0,0], /NoData, XRange=[(n_elements(beam_cont[*,0])-1)/2, -(n_elements(beam_cont[*,0])-1)/2],XStyle=1, $
			YRange=[-(n_elements(beam_cont[0,*])-1)/2, (n_elements(beam_cont[0,*])-1)/2],YStyle=1, $
			XTitle='X Offset [arc min]', YTitle='Y Offset [arc min]', $
			Position=[0.15,0.15,0.95,0.95], /NoErase

		xyouts, 0.18, 0.91, 'Location (i,j):     '+beam_posg, /Norm
		xyouts, 0.18, 0.88, 'Location (ra, dec): '+beam_poss, /Norm
		xyouts, 0.18, 0.85, 'FWHM: '+beam_fwhm, /Norm
		xyouts, 0.18, 0.82, 'PA:   '+beam_pa, /Norm
		xyouts, 0.18, 0.79, 'Correction Factor: '+beam_corr, /Norm
		
		snapshot = tvrd()
   		tvlct, r, g, b, /Get
  		device, Z_Buffer=1
  		set_plot, thisDevice
		
		image24 = BytArr(3, 700, 700)
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
