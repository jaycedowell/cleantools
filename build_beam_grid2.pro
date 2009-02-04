; build_beam_grid2 - Given a box where we need to deconvolve, divide the box into sub-
; regions and compute the effective beam pattern in each.  Also, return a map with the 
; flux corrections to handle the area miss-match between the dirty and clean beams.
;
; Usage:  build_beam_grid2, grid, llx, lly, urx, ury, beams [, FluxCorr=<flux correction factor>]
;         	[, FWHM=<FWHM to use for restore beam>] [, PA=<PA to use for restore beam>] 
;		[, /Silent] [, Output=<output from beam builds] [, /DecOnly] [, /Smart] 
;		[, /SingleBeam] [, /CExts]
;
; Inputs:
;  grid   -> grid structure
;  llx    -> x-value of the lower left-hand for the box where beams need to be built
;  lly    -> y-value of the lower left-hand for the box where beams need to be built
;  urx    -> x-value of the upper right-hand for the box where beams need to be built
;  ury    -> y-value of the upper right-hand for the box where beams need to be built
;  FWHM   -> The FWHM of the Gaussian restore beam computed from the central beam.
;  PA     -> The position angle (E of N) for the Gaussian restore beam.
;  Silent -> Turn off information statements.  The default is to display them.
;
; Modeling Keywords:
;  DecOnly -> When this keyword is set, beams are modeled for every declination in the map.  
;             By doing this, the method captures most of the changes in the beam shape that 
;             are associated with change in the beam that contributes most of the flux to a 
;	      given pixel.  This method is fairly fast as only (ury-lly+1) beams are needed.
;  Smart   -> This is similar to DecOnly, but the routine try to detect additional changes in
;             the beam shape that occur along in the RA direction.  It does this by looking 
;             at the weight map computing a new beam where the pixel weight changes by >=5%.
;             This is not as fast as DecOnly because it models more beams.  It is also not
;             entirely clear if doing this is really necessary, i.e., if changes in the beam
;             shape due to the addition of loss of beams really cause 5% changes in the weight
;
; Outputs:
;  beams    -> 4-D array (grid x, grid y, beam x, beam y) containing a beam map for each point
;              in the box specified by llx->urx, lly->ury
;  FluxCorr -> 2-D array 
;  Output   -> All logging statements are saved to a string array.  This array is returned with 
;              this varaible and is useful for interfacing with other programs.
; Change Log:
; 07/11/08 - JD - Initial version 
; 07/16/08 - JD - Added Silet and Ouput keywords/variables
;               - Created version 2 to compute at every grid point.
; 07/25/08 - JD - Added documentation on how to use it
;		- Added option for variable beam map size
; 07/29/08 - JD - Modified to accomidate build_beam5
; 07/31/08 - JD - Added `Dec Only` option
; 09/23/08 - JD - Added `Single Beam` option to make the coding of the automatic deconvolution
;                 routine easier.  See alfalfa_deconv2 for details of the new scheme that lead
;                 to this keyword being added.
;               - Updated the internal documentation to explain the various modeling keywords 
;                 and the beam replication part
;
pro build_beam_grid2, grid, llx, lly, urx, ury, beams, MapSize=MapSize, FluxCorr=FluxCorr, FWHM=FWHM, PA=PA, Silent=Silent, Output=Output, CExts=CExts, DecOnly=DecOnly, Smart=Smart, SingleBeam=SingleBeam

common agcshare

; Start timer
t_start = systime(1)

if Not Keyword_Set(Silent) then begin
	print,''
	print,'build_beam_grid2 started:'
	print,'> setup'
endif
Output = 'build_beam_grid2 started:'
Output = [Output, '> setup']

; Look for the MapSize option and handle it accordingly
if n_elements(MapSize) EQ 0 then MapSize=25

; Create results array `beams` that is a 4-D array.  It will store the beam at each 
; location.  The first two dimensions correspond to the RA,Dec (x,y) in the flux box.
; The last two store the i,j points of the MapSize' x Mapsize' beam.
beams = dblarr(urx-llx+1, ury-lly+1, MapSize, MapSize)
; Create the flux correction array that will store the corrections needed to handle
; the area miss-match
FluxCorr = dblarr(urx-llx+1, ury-lly+1)

; Create a variable (internal to this routine) to hold the results of the Gaussian fits.
GaussOut = dblarr(urx-llx+1, ury-lly+1, 7)

; Time after setup of variables
t_setup = systime(1)

; Restore the appropriate beam maps
case MapSize of 
	21:	restore,agcdir+'../deconvolution/beams/a1963_21_180.sav'
	25:	restore,agcdir+'../deconvolution/beams/a1963_25_180.sav'
	31:	restore,agcdir+'../deconvolution/beams/a1963_31_180.sav'
	35:	restore,agcdir+'../deconvolution/beams/a1963_35_180.sav'
	41:	restore,agcdir+'../deconvolution/beams/a1963_41_180.sav'
	45:	restore,agcdir+'../deconvolution/beams/a1963_45_180.sav'
	51:	restore,agcdir+'../deconvolution/beams/a1963_51_180.sav'
	else:	restore,agcdir+'../deconvolution/beams/a1963_25_180.sav'
endcase

fitl = MapSize/2 - 6
fitu = MapSize/2 + 6

if Keyword_Set(DecOnly) then begin
	; Now loop through and build the beams we need
	if Not Keyword_Set(Silent) then begin
		print,'> building beams using DecOnly Mapping'
	endif
	Output = [Output, '> building beams using DecOnly Mapping']

	grid_i = (urx+llx)/2
	box_i = grid_i - llx

	for m=lly,ury do begin
		grid_j = m
		box_j = m - lly
		
		; This next part is meant to be a basic error check.  For some unknown reason (compiler flags, coding problem, etc.)
		; build_beam_worker returns a beam that is not kosher (FluxCorr < 1).  So, check for this and redo beams that are
		; suspect.  
		;while FluxCorr[box_i, box_j] LT 1 do begin
			build_beam5, grid, [grid_i, grid_j], tile_beam, MapSize=MapSize, /Silent, Output=OutputTemp, AntPat=ant_pat, CExts=CExts
			beams[box_i, box_j, *, *] = tile_beam / max(tile_beam)
			
			junk = gauss2dfit( tile_beam[fitl:fitu, fitl:fitu], temp_beam_out, /Tilt)		; Fit Gaussian to central beam
			GaussOut[box_i, box_j, *] = temp_beam_out
	
			FluxCorr[box_i, box_j] = total(tile_beam) / max(tile_beam)	; Area of beam for FluxCorr
		;endwhile

		; Save output to the Output variable
		Output = [Output, OutputTemp]
	endfor
endif else begin
	if Keyword_Set(Smart) then begin
		; Now loop through and build the beams we need
		if Not Keyword_Set(Silent) then begin
			print,'> building beams using Smart Mapping'
		endif
		Output = [Output, '> building beams using Smart Mapping']

		Weights = total(total(grid.w,1,/Double),1,/Double)/2.0/n_elements(grid.nz)

		; Compute where the beam changes by 5% or more
		x_size = urx-llx+1
		y_size = ury-lly+1
		
		refs = 0.0d
		to_comp = fix(Weights*0)
		index = 0
		
		for j=lly,ury do begin
			ref = Weights[(urx+llx)/2-llx, j]
			to_comp[(urx+llx)/2-llx, j-lly] = index+1
			index = index+1
			
			build_beam5, grid, [(urx+llx)/2, j], tile_beam, MapSize=MapSize, /Silent, Output=OutputTemp, AntPat=ant_pat, CExts=CExts
			beams[(urx+llx)/2-llx, j-lly, *, *] = tile_beam / max(tile_beam)
			Gjunk = gauss2dfit( tile_beam[fitl:fitu, fitl:fitu], temp_beam_out, /Tilt)	; Fit Gaussian to central beam
			GaussOut[(urx+llx)/2-llx, j-lly, *] = temp_beam_out
			FluxCorr[(urx+llx)/2-llx, j-lly] = total(tile_beam) / max(tile_beam)		; Area of beam for FluxCorr
			Output = [Output, OutputTemp]

			BFill = tile_beam/max(tile_beam)
			FFill = total(tile_beam) / max(tile_beam)
			
			for i=llx,urx do begin
				pdiff = 100.0 * (Weights[i,j] - ref)/ref 
				if abs(pdiff) GE 10 then begin
					to_comp[i,j] = index+1
					ref = Weights[i,j]
					index = index+1

					build_beam5, grid, [(urx+llx)/2, j], tile_beam, MapSize=MapSize, /Silent, Output=OutputTemp, $
						AntPat=ant_pat, CExts=CExts
					beams[i-llx, j-lly, *, *] = tile_beam / max(tile_beam)
					Gjunk = gauss2dfit( tile_beam[fitl:fitu, fitl:fitu], temp_beam_out, /Tilt)
					GaussOut[i-llx, j-lly, *] = temp_beam_out
					FluxCorr[i-llx, j-lly] = total(tile_beam) / max(tile_beam)
					Output = [Output, OutputTemp]

					BFill = tile_beam/max(tile_beam)
					FFill = total(tile_beam) / max(tile_beam)
				endif else begin
					to_comp[i,j] = max(to_comp)

					beams[i-llx, j-lly, *, *] = BFill
					FluxCorr[i-llx, j-lly] = FFill
				endelse
			endfor
			
		endfor

		; Report time savings
		n_comp_beams = n_elements( uniq(to_comp[sort(to_comp)]) )
		n_chain_beams = 0
		junk = (to_comp[sort(to_comp)])[ uniq(to_comp[sort(to_comp)]) ]
		for i=0,(n_elements(junk)-1) do begin
			junk2 = where( to_comp EQ junk[i], count )
			if count GT n_chain_beams then n_chain_beams = count
		endfor
		if Not Keyword_Set(Silent) then begin
			print,'>  '+strtrim(string(n_comp_beams),2)+' beams to compute'
			print,'>  longest beam chain is '+strtrim(string(n_chain_beams),2)
			print,'>  area savings is '+strtrim(string(100.-100.0*n_comp_beams/x_size/y_size,Format='(F4.1)'),2)+'%'
		endif
		Output = [Output, '>  '+strtrim(string(n_comp_beams),2)+' beams to compute', $
				'>  longest beam chain is '+strtrim(string(n_chain_beams),2), $
				'>  area savings is '+strtrim(string(100.0*n_comp_beams/x_size/y_size,Format='(F4.1)'),2)+'%' ]

	endif else begin
		if n_elements(SingleBeam) NE 0 then begin
			if Not Keyword_Set(Silent) then begin
				print,'> building beams using SingleBeam Mapping'
			endif
			Output = [Output, '> building beams using SingleBeam Mapping']

			; Interpret the keyword.  If it is just set, then we will just do the central pixel.
			; If a two-element array is supplied, then we use that as the central pixel.  This is 
			; needed to ensure the the center of the galaxy is modeled when the edge expansion 
			; encounters a grid edge.
			if n_elements(SingleBeam) EQ 1 then begin
				grid_i = (llx+urx)/2
				grid_j = (lly+ury)/2
			endif else begin
				grid_i = SingleBeam[0]
				grid_j = SingleBeam[1]
			endelse
			; This next part looks a little strange:  Here's the reason.  We want to know the beam 
			; at a given grid point but we are going to make all pixels in the output beam map have the
			; same shape.  So, it really doesn't matter where we map the output beam to.  By using the
			; center we are making it easier to combine this method with the beam fitting methods around
			; line 280.
			box_i = (llx+urx)/2 - llx
			box_j = (lly+ury)/2 - lly

			build_beam5, grid, [grid_i, grid_j], tile_beam, MapSize=MapSize, /Silent, Output=OutputTemp, AntPat=ant_pat, CExts=CExts
			beams[box_i, box_j, *, *] = tile_beam / max(tile_beam)
			
			junk = gauss2dfit( tile_beam[fitl:fitu, fitl:fitu], temp_beam_out, /Tilt)		; Fit Gaussian to central beam
			GaussOut[box_i, box_j, *] = temp_beam_out
	
			FluxCorr[box_i, box_j] = total(tile_beam) / max(tile_beam)	; Area of beam for FluxCorr

			; Save output to the Output variable
			Output = [Output, OutputTemp]
		
		endif else begin
			; Else, do the full modeling process
			if Not Keyword_Set(Silent) then begin
				print,'> building beams using Full Mapping'
			endif
			Output = [Output, '> building beams using Full Mapping']
	
			for l=llx,urx do begin
				grid_i = l
				box_i = l - llx
			
				for m=lly,ury do begin
					grid_j = m
					box_j = m - lly
					
					build_beam5, grid, [grid_i, grid_j], tile_beam, MapSize=MapSize, /Silent, Output=OutputTemp, AntPat=ant_pat, CExts=CExts
					beams[box_i, box_j, *, *] = tile_beam / max(tile_beam)
					
					junk = gauss2dfit( tile_beam[fitl:fitu, fitl:fitu], temp_beam_out, /Tilt)		; Fit Gaussian to central beam
					GaussOut[box_i, box_j, *] = temp_beam_out
			
					FluxCorr[box_i, box_j] = total(tile_beam) / max(tile_beam)	; Area of beam for FluxCorr
			
					; Save output to the Output variable
					Output = [Output, OutputTemp]
				endfor
			endfor
		endelse
	endelse
endelse

if Not Keyword_Set(Silent) then begin
	for i=3,(n_elements(Output)-1) do $
		print, '>  '+Output[i]
endif

; Find the best-fit Gaussian to the beam ensemble
if Not Keyword_Set(Silent) then begin
	print,'> determining clean beam FWHM'
endif
Output = [Output, '> determining clean beam FWHM']
center_gauss = reform(GaussOut[(urx+llx)/2-llx, (ury+lly)/2-lly, *])
center_area = 2.0*!Pi*center_gauss[2]*center_gauss[3]
; Finish the primary flux correction array using these values
FluxCorr = FluxCorr / center_area
; Setup Gaussian return values
FWHM = 2.0*sqrt(2.0*alog(2.0))*[center_gauss[2], center_gauss[3]]
PA = center_gauss[6]*!radeg
if Not Keyword_Set(Silent) then begin
	print, '>  clean beam is '+string(FWHM[0],Format='(F5.3)')+''''+' x '+string(FWHM[1],Format='(F5.3)')+''''
endif
Output = [Output, '>  clean beam is '+string(FWHM[0],Format='(F5.3)')+''''+' x '+string(FWHM[1],Format='(F5.3)')+'''']

; Time after primary beam building/flux correction
t_end = systime(1)

; Fill in the map, if needed.  SMART fills its own map because of how it determines where the 
; beam needs to be recomputed.  FULL (the DEFAULT) fills its own map because it computes the 
; beam for every point in the map.  That leaves DECONLY and SINGLEBEAM that need to be helped
; out.  The two methods used below are very similar but differ in how they assign the beam to the
; rest of the map.  DECONLY copy only along lines of contant RA while SINGLEBEAM copies the 
; central beam to all points.  
if Keyword_Set(DecOnly) then begin
	if Not Keyword_Set(Silent) then begin
		print, '> filling in RA'
	endif
	Output = [Output, '> filling in RA']

	for l=llx,urx do begin
		box_i = l - llx
		
		for m=lly,ury do begin
			box_j = m - lly

			beams[box_i, box_j, *, *] = beams[(urx+llx)/2-llx, box_j, *, *]
			FluxCorr[box_i, box_j] = FluxCorr[(urx+llx)/2-llx, box_j]
		endfor
	endfor
endif else begin
	if n_elements(SingleBeam) NE 0 then begin
		if Not Keyword_Set(Silent) then begin
			print, '> filling in beam map'
		endif
		Output = [Output, '> filling in beam map']
		
		for l=llx,urx do begin
			box_i = l - llx
			for m=lly,ury do begin
				box_j = m - lly
	
				beams[box_i, box_j, *, *] = beams[(urx+llx)/2-llx, (ury+lly)/2-lly, *, *]
				FluxCorr[box_i, box_j] = FluxCorr[(urx+llx)/2-llx, (ury+lly)/2-lly]
			endfor
		endfor
			
	endif
endelse

if Not Keyword_Set(Silent) then begin
	print, '> total time '+strtrim(string(t_end-t_start),2)+' seconds'
endif
Output = [Output, '> total time '+strtrim(string(t_end-t_start),2)+' seconds']

end
