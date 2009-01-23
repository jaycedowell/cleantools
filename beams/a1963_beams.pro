; bicubic_spline - Bicubic interpolation with a spline.  
function bicubic_spline, Z, X, Y, Xp, Yp
; Create output arraus
new = dblarr((size(X))[1], (size(Yp))[1])
new2 = dblarr((size(Xp))[1], (size(Yp))[1])
; Loop 1:  Loop over columns (x,y) -> (x,y')
val = reform(Y[0,*])
unk = reform(Yp[0,*])
for i=0L,(n_elements(Y[*,0])-1) do begin
	column = reform(Z[i,*])
	temp = interpol(column, val, unk, /Spline)
	new[i,*] = temp
endfor
; Loop 2: Loop over rows (x,y')->(x',y')
val = reform(X[*,0])
unk = reform(Xp[*,0])
for i=0L,(n_elements(Yp[0,*])-1) do begin
	row = reform(new[*,i])
	temp = interpol(row, val, unk, /Spline)
	new2[*,i] = temp
endfor
; Return and exit
return,new2
end



pro a1963_beams, final, Width=Width, ReCenter=ReCenter, Azimuth=Azimuth, Verbose=Verbose, File=File
; Check if we are using the Verbose flag.  If so, display stuff
if Keyword_Set(Verbose) then begin
	window,0,Title='Original TV Image',XSize=512,YSize=512
	window,1,Title='Final TV Image',XSize=512,YSize=512
	loadct, 39, /Silent
endif

if n_elements(Azimuth) EQ 0 then Azimuth = 153.0
az_mstr = [104.0, 107.0, 109.0, 116.0, 128.0, 153.0, 200.0, 230.0, 245.0, 250.0, 253.0, 255.0]
if Azimuth LT min(az_mstr) or Azimuth GT max(az_mstr) then begin
	print,'% Azimuth is out of range. Valid range is '+string(min(az_mstr),Format='(I3)')+$
		' -> '+string(max(az_mstr),Format='(I3)')+' degrees'
	goto,TheEnd
endif
if Azimuth EQ min(az_mstr) or Azimuth EQ max(az_mstr) then begin
	az_use = [min(az_mstr), max(az_mstr)]
endif else begin
	temp1 = az_mstr - Azimuth
	temp2 = temp1*0+10000.0
	for i=0L,(n_elements(temp1)-2) do temp2[i] = temp1[i]*temp1[i+1]
	good = min( where( temp2 LE 0 ))
	az_use = az_mstr[[good,good+1]]
endelse
if Keyword_Set(Verbose) then begin
	print,'% Build beams for Azimuth '+string(Azimuth,Format='(F5.1)')+' degrees'
	print,'% Averaging beams ',az_use
endif

; Setup result variables.  These stay the same for each run so
; we really only need to do this once.
Width = fix( Width / 0.05 )		; arc min -> pixels
if Width mod 2 EQ 0 then Width++	; give us a well defined center
new = dblarr(Width,Width)		; new array (Width+.05)'x(Width+0.05)' @ 0.05'/px
new_ra = new
new_dec = new
for j=0L,(Width-1) do begin
	new_ra[j,*] = (j-Width/2)*0.05d0
	new_dec[*,j] = (j-Width/2)*0.05d0
endfor
final1 = fltarr(2, 7, Width, Width)	; final1 array for two azimuths of 7 ALFA beams
final = fltarr(7, Width, Width)		; the real final

; Make sure bicubic_spline is loaded
if total(strcmp( routine_info(/Functions), 'BICUBIC_SPLINE')) EQ 0 then $
	resolve_routine,'a1963_beams',/Compile_Full_file, /Either

; The main loop over the azimuths
for a=0,1 do begin
	az_name = '_AZ'+string(az_use[a],Format='(I3)')

	; The main loop over beams
	for i=0,6 do begin
		; Clear new
		new = new *0.0
		; Load the FITS file and convert to double
		if Keyword_Set(Verbose) then $
			print,'% Loading Azimuth '+string(az_use[a],Format='(I3)')+$
				', Beam '+strtrim(string(i),2)
		current_beam = 'A1963/FINALBM'+strtrim(string(i),2)+az_name+'.FITS'
		beam = readfits(current_beam, hdr, /Silent)
		beam = double(beam) ;double(beam[384:639,*])
	
		; Get the size and location of the central peak through a simple Gaussian fit
		junk = gauss2dfit(beam, GaussOut, /Tilt)
		x = (size(beam))[1]
		y = (size(beam))[2]
		center_x = round(GaussOut[4])
		center_y = round(GaussOut[5])
	
		; Get the RA and Dec scale.  These came from AIPS so the units are
		; deg / px.  
		ra_scale = sxpar(hdr, 'CDELT1') * 60.0d0		; arc min
		dec_scale = sxpar(hdr, 'CDELT2') * 60.0d0		; arc min
	
		; Computer the coordinates of each pixed in 'beam' relative to the central peak
		ra = dblarr(x,y)
		dec = dblarr(x,y)
		for j=0,(x-1) do ra[j,*] = j
			ra = (ra - center_x) * ra_scale
		for j=0,(y-1) do dec[*,j] = j
			dec = (dec - center_y) *dec_scale
	
		; Re-grid the data using bicubic_spline.  
		if Keyword_Set(Verbose) then $
			print,'% Re-griding data'
		new = bicubic_spline(beam, ra, dec, new_ra, new_dec)
		; Normalize the output
		new = new / total(new)
	
		; Display if requested.
		if Keyword_Set(Verbose) then begin
			; Footwork to get the aspect ratio work
			if( x*abs(ra_scale) GE y*dec_scale ) then begin
				cgx = 512
				cgy = floor(512.0*(y*dec_scale)/abs(x*ra_scale))
			endif else begin
				cgx = floor(512.0*abs(x*ra_scale)/(y*dec_scale))
				cgy = 512
			endelse

			; Erase each window, slap up the plot, and print a label
			wset,0
			erase
			tvscl,congrid(beam,cgx,cgy)
			xyouts,0.1,0.9,'AZ: '+string(az_use[a],Format='(I3)')+', Beam: '+string(i,Format='(I1)'),/Norm, $
				CharSize=1.1

			wset,1
			erase
			tvscl,congrid(new,512,512)
			xyouts,0.1,0.9,'AZ: '+string(az_use[a],Format='(I3)')+', Beam: '+string(i,Format='(I1)'),/Norm, $
				CharSize=1.1
		endif
	
		; Load into the output array to the correct location
		final1[a,i,*,*] = new
	endfor
endfor

; Clean up display windows
if Keyword_Set(Verbose) then begin
	wdelete,0
	wdelete,1
endif


; If the center flag is set, center the max of each beam at the center
if Keyword_Set(ReCenter) then begin
	for a=0,1 do begin
		for i=0,6 do begin
			c = reform(final1[a,i,*,*])
			l = where( c EQ max(c) )
			x = l mod Width
			y = l / Width
			
			c = shift(c, (Width/2-x), (Width/2-y))
			final1[a,i,*,*] = c
			l = where( c EQ max(c) )
			x = l mod Width
			y = l / Width
		endfor
	endfor
endif

; Average the two beam together to get the correct azimuth angle. 
z = 1.0							; two equations, three unknowns, fix one of them
y = (az_use[1] - Azimuth) / (Azimuth - az_use[0])	; solve for y/z
if finite(y) NE 1 then begin				; For some edge cases, this is inifnite.
	z = 0.0
	y = 1.0
endif
; Check to see if we can make some rational numbers out of these
temp = y*(findgen(100)+1)
good = min( where( temp EQ fix(temp), count) )
if y NE 1 and count NE 0 then begin
	y = double(temp[good])
	z = double(good+1)
endif
x = (y + z)
; Tell us about it.
if Keyword_Set(Verbose) then begin
	print,'% Beam weighting is '+string(y,Format='(I3)')+' for '+string(az_use[0],Format='(I3)')+$
		' and '+string(z,Format='(I3)')+' for '+string(az_use[1],Format='(I3)')
endif
; Average everything
for i=0,6 do final[i,*,*] = (y*final1[0,i,*,*]+z*final1[1,i,*,*])/x

; Clean 
ans='N'
read,ans,Prompt='Do you wish to inspect the maps and clean them? (Y/[N])'
if strcmp(ans,'Y') or strcmp(ans,'y') then begin
	; Create output window and print a few simple instructions
	window,0,Title='Cleaning Window',XSize=512,YSize=512
	print,''
	print,'Left Click - Add verticies to the masking polygon'
	print,'Right Click - Close the current polygon'
	print,'Middle Click - Apply mask and move to the next beam'
	print,''

	for i=0,6 do begin
		cur = reform(final[i,*,*])
		mask = fix(cur*0.0+1.0)
		out = cur / max(cur)
		loadct,39,/Silent
		tvimage,bytscl(congrid(alog10(out)*10.0, 512, 512), min=-24.0, max=0.0)
		xyouts,0.90,0.90,'Beam '+string(i,Format='(I1)'),/Norm, $
			Color=FSC_Color('White'),Alignment=1.0

		Window, 19, /Pixmap, XSize=512, YSize=512
		Device, Copy=[0, 0, 512, 512, 0, 0, 0]

		WSet,0
		while !Mouse.Button NE 2 do begin
			cursor, x_temp, y_temp, /Down, /Norm

			Device, Copy=[0, 0, 512, 512, 0, 0, 19]

			plots,x_temp,y_temp,PSYm=1,Color=FSC_Color('Cyan'), /Norm
			polygon_x = [x_temp]
			polygon_y = [y_temp]
			polygon_i = 1

			while !Mouse.Button NE 4 and !Mouse.Button NE 2 do begin
				cursor, x_temp, y_temp, /Down, /Norm
	
				Device, Copy=[0, 0, 614, 614, 0, 0, 19]

				polygon_x = [polygon_x, x_temp]
				polygon_y = [polygon_y, y_temp]
				plots,polygon_x,polygon_y,PSym=1,Color=FSC_Color('Cyan'), /Norm
				plots,polygon_x,polygon_y,LineStyle=2,Color=FSC_Color('Cyan'), /Norm
				polygon_i = polygon_i + 1


				polyfill, [polygon_x, polygon_x[0]], [polygon_y, polygon_y[0]], $
					/Line_Fill, Orientation=45, /Norm, Color=FSC_Color('Red')
				polyfill, [polygon_x, polygon_x[0]], [polygon_y, polygon_y[0]], $
					/Line_Fill, Orientation=-45, /Norm, Color=FSC_Color('Red')
			endwhile

			if polygon_i GT 2 then begin
				temp1 = [0, polygon_i-1]
				plots,polygon_x[temp1],polygon_y[temp1],LineStyle=2,Color=FSC_Color('Red'), /Norm

				polyfill, [polygon_x, polygon_x[0]], [polygon_y, polygon_y[0]], $
					/Line_Fill, Orientation=45, /Norm, Color=FSC_Color('Red')
				polyfill, [polygon_x, polygon_x[0]], [polygon_y, polygon_y[0]], $
					/Line_Fill, Orientation=-45, /Norm, Color=FSC_Color('Red')

				Window, 19, /Pixmap, XSize=512, YSize=512
				Device, Copy=[0, 0, 512, 512, 0, 0, 0]
				WSet,0

				polygon_x = round(polygon_x*Width)
				polygon_y = round(polygon_y*Width)
				good = polyfillv(polygon_x, polygon_y, Width, Width)
				print,'Beam '+string(i,Format='(I1)')+': Masking ',string(n_elements(good),Format='(I5)')+' pixels inside a '+$
					string(n_elements(polygon_x),Format='(I3)')+'-sided polygon'
				; Clever trick to build a 1-D index
				mask[good] = 0
			endif
		endwhile
		final[i,*,*] = cur*mask

		WDelete, 19

		!Mouse.Button=0
	endfor

	wdelete,0
endif

; Save the output to a file
ant_pat = final
if n_elements(File) NE 0 then begin 
	filename = File
endif else begin
	filename = 'a1963_'+string(Width*0.05,Format='(I2)')+'_'+$
		string(Azimuth,Format='(I3)')+'.sav'
endelse
save, ant_pat, filename=filename

TheEnd:

end
