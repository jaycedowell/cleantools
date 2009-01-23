; combine_chans - Combine channels different ways to improve the per channel S/N.  The factor
; that determines how many to add together is CleanLimit.  Combine_chans adds channels together
; until at least 75% of channels have a peak > CleanLimit.
;
; Usage:  combine_chans, data_cube, CleanLimit [,NComb = <# of channels in new psuedo channel>] $
;	[,/Silent] [,Output=<strarr of output>]
;
; 
function combine_chans, data, CleanLimit, NComb=NComb, AcceptFrac=AcceptFrac, Silent=Silent, Output=Output
	; If the acceptable fraction is not defined, set it to 0.75
	if n_elements(AcceptFrac) EQ 0 then $
		AcceptFrac = 0.75

	; Define size of data cube
	vs = (size(data))[1]	; Size of velocity array
	xs = (size(data))[2]	; RA size
	ys = (size(data))[3]	; Dec size

	; Define cube RMS
	data_rms = stddev(data)
	
	; Find most of the ways to split up the channels
	temp = vs / (findgen(vs)+1)	; start at one so that if we have a really bright source
					; we get the full channel resolution
	smth = where( temp EQ fix(temp) ) + 1

	; Do the sum and figure out the new velocities associated with each
	; new channel
	for t=0L,(n_elements(smth)-1) do begin
		data2 = dblarr(vs/smth[t], xs, ys)
		metric = 0
		for i=0,(vs/smth[t]-1) do begin
			for j=(i*smth[t]),((i+1)*smth[t]-1) do data2[i,*,*] += data[j,*,*]
			if max(data2[i,*,*]) GT CleanLimit then metric += 1
		endfor
		if metric GE AcceptFrac*vs/smth[t] then goto, GoodEnough

	endfor
	
	GoodEnough:
	; Define the new cube RMS
	data2_rms = stddev(data2)
	; Store the number of channels we are combining
	NComb = smth[(t<(n_elements(smth)-1))]

	; Print some interesting information
	if Not Keyword_Set(Silent) then begin
		print,'> combine_chan: combining '+$
			string(smth[(t<(n_elements(smth)-1))],Format='(I3)')+' channels'
	endif
	Output = [Output, '> combine_chan: combining '+$
			string(smth[(t<(n_elements(smth)-1))],Format='(I3)')+' channels']

	; Return the data cube
	return,data2
end

; alfa_clean7 - Deconvolve a data cube with a CLEAN-type (Hogbom 1974) method.  The 
; procedure accepts a number of inputs that control the CLEANing process, including map 
; RMS, the total number of iterations, and the FWHM of the restored beam.  The procedure
; also assumes that the beams provided are N' x N' @ 1'/px where N is an odd number, 
; i.e., a beam sequence created by build_beam_grid.pro or build_beam_grid2.pro.
;
; The procedure works as follows:
; (1) The data cube is read into the program as are all of the CLEAN control options.
; (2) The cleaning loop selects the maximum point in each map for cleaning.  If more 
;     than one point fits this selection criterion then the one closest to the lower 
;     left-hand corner is used.  The cleaning process removes a total flux from the 
;     map that is 5% of the maximum value.  The loop continues until
;      (a) the RMS limit is hit (Sigma*RMS > map maximum point)
;      (b) the Flux limit is hit (Flux > map maximum point)
;      (c) the iteration limit is hit
;    If more than one flux flag (Flux,Sigma) is selected, the larger limit is used.
; (3) The CLEANed data is re-convolved with a Gaussian beam with a FWHM of 3.5' or 
;     whatever is specified by the FWHM flag.  This re-convolved map is then added back
;     to the residual flux to create the output map.  
;
; Usage:  alfa_clean7, data_cube, beams, cleaned_map [, Continuum=<continuum_map>] [,C_Cleand=C_Cleand] 
;		[, MapRMS=<map rms from grid.grms>] NIter=<max # of iterations>] $
;         	[,Flux = <Flux limit>] [, Gain = <CLEAN loop gain>] [,FWHM = <FWHM of Gaussian restore beam>] $
;               [, PA = <PA of Guassian restore beam>] [,Sigma = <Sigma limit>] [,Range = <channel range to clean>] $
;               [,/AllowSum] [,Exit_Status = <numererical exit status] [,/Silent] [,Output=<strarr of output>]
;
; Inputs:
;  data_cube -> [v, ra, dec] data cube
;  beam -> 4-D cube contining N' x N' @ 1'/px beam from build_beam_grid2.pro
;  Continuum -> [ra, dec] map of continuum sources for cleaning
;  MapRMS -> RMS of the grid that the map is from
;  NIter -> Maximum number of iterations in the CLEAN loop.  The default is 2,000.
;  Flux -> The flux level, in mJy, to clean down to.
;  Sigma -> Sigma level with which to clean down to.  
;  Gain - > Gain to use in the CLEAN loop.  The default is 0.05.
;  Range -> Two element array of channel number to clean
;  AllowSum -> Allow alfa_clean7 to sum channels together to improve the per channel S/N
;  FWHM -> The FWHM of the Gaussian restore beam.  The default is 3.1' x 4.6'
;  PA -> The position angle (E of N) for the Gaussian restore beam.  The default is 0 deg.
;  Silent -> Turn off information statements.  The default is to display them.
;  Output -> All logging statements are saved to a string array.  This array is returned with 
;            this varaible and is useful for interfacing with other programs.
;
; Outputs:
;  cleaned_map -> The velocity integrated CLEANed map.  This map as the same RA and
;                 dec dimensions as the src.srccube.dbox data
;  c_cleand -> Cleaned version of continuum map (if supplied)
;  Exit_Status -> Numerical exit status for each channel's cleaning loop.  The codes are:
;                  -1 => exit caused by iteration limit
;                  +1 => exit caused by RMS limit
;  Output -> All logging statements are saved to a string array.  This array is returned with 
;            this varaible and is useful for interfacing with other programs.
;
; Information Statements:
;  This procedure produces a variety of information statements.  These include current 
;  task (setup, cleaning, and restore), total run time is seconds, name of the current
;  source, CLEAN loop control limits, the loop exit status, and the flux contained in 
;  the map.  These can be turned off with the Silent flag
;
pro alfa_clean7, dbox, beams, cleand, Continuum=Continuum, C_Cleand=C_cleand, MapRMS=MapRMS, $
	NIter=NIter, Flux=Flux, Sigma=Sigma, Gain=Gain, Range=Range, AllowSum=AllowSum, $
	Exit_Status=Exit_Status, FWHM=FWHM, PA=PA, Silent=Silent, Output=Output, CExts=CExts
t_start = systime(1)

; Determine box size
v_size = (size(dbox))[1]
x_size = (size(dbox))[2]
y_size = (size(dbox))[3]
if Not Keyword_Set(Silent) then begin
	print,''
	print,'alfa_clean7 started:'
	print,'> setup'
	print,'>  loaded '+string(v_size,Format='(I4)')+'x'+string(x_size,Format='(I3)')+'x'+ $
		string(y_size,Format='(I3)')+' data cube'
endif
Output = ['alfa_clean7 started:','> setup']
Output = [Output, '>  loaded '+string(v_size,Format='(I4)')+'x'+string(x_size,Format='(I3)')+ $
	'x'+string(y_size,Format='(I3)')+' data cube']

; Setup variabled that control the cleaning process.  
if n_elements(NIter) EQ 0 then NIter = 2000
; Check to see if the Sigma flag is set.  If not, set sigma to something
; rediculous so that we don't exit on this.
if n_elements(Sigma) EQ 0 then Sigma=0.0
; Same thing for flux. 
if n_elements(Flux) EQ 0 then Flux = 0.0
; If niether are set, default to a 10 Sigma limit (Method 1)
if Flux EQ 0.0 AND Sigma EQ 0.0 then Sigma = 10.0
; If both are set, pick the maximum and set flux to that.  Then, set Sigma to 0
if Flux NE 0.0 AND Sigma NE 0.0 then begin
	Flux = max([Flux, Sigma*MapRMS])
	Sigma = 0.0
endif
; Gain Settings for the CLEAN loop
if n_elements(Gain) EQ 0 then Gain = 0.05

; Select range to operate on
if n_elements(Range) NE 2 then Range=[0L,1023L]

; Combine channels if asked
if KeyWord_Set(AllowSum) then begin
	dbox_old = dbox
	dbox = combine_chans(dbox_old[Range[0]:Range[1], *, *], (Flux+Sigma*MapRMS), NComb=NComb, $
		Silent=Silent, Output=Output)
	cleand = dbox_old*0.0
	C_Range = [0L, (n_elements(dbox[*,0,0])-1L)]
endif else C_Range = Range

v_size = (size(dbox))[1]
x_size = (size(dbox))[2]
y_size = (size(dbox))[3]

if Not Keyword_Set(Silent) then begin
	print,'> clean'
	print,'>  cleaning data cube to a peak flux of '+strtrim(string(Flux+Sigma*MapRMS),2)+ $
		' mJy or '+strtrim(string(NIter),2)+' iterations'
endif
Output = [Output, '> clean', '>  cleaning data cube to a peak flux of '+ $
	strtrim(string(Flux+Sigma*MapRms),2)+' mJy or '+ $
	strtrim(string(NIter),2)+' iterations']

; Setup the holding variables.  Working starts off the same as data but slowly has peaks 
; removed.  cleand starts off as zeros and slowly has peaks added to it.  
cleand = dbox*0.0
working = dbox

; Determine Beam half-size
b_size = (size(beams))[3]/2
; Determine beam normalizations
b_norm = total(total(beams, 3, /Double), 3, /Double)

; Run through the CLEANing loop
if Keyword_Set(CExts)  then begin
	if Not Keyword_Set(Silent) then print,'>  using C library '+lib_version()
	Output = [Output, '>  using C library '+lib_version()]

	exit_status = lonarr( v_size ) - 2
	cleand2 = dblarr(x_size, y_size, v_size)
	working2 = dblarr(x_size, y_size, v_size)
	cleand2 = transpose(cleand, [1,2,0])
	working2 = transpose(working, [1,2,0])

	cleand2=alfa_clean_worker(cleand2, working2, beams, C_Range[0], C_Range[1], Flux+Sigma*MapRms, NIter, Gain, $
			b_norm, (size(beams))[3], x_size, v_size, exit_status)

	cleand = transpose(cleand2, [2,0,1])
	working = transpose(working2, [2,0,1])

	; Strip out excess
	exit_status = fix(exit_status[C_Range[0]:C_Range[1]])

endif else begin
	exit_status = intarr( C_Range[1]-C_Range[0]+1 ) -1

	; Loop over channels
	for c=C_Range[0],C_Range[1] do begin
		for l=0L,(NIter-1) do begin
			; Select only one pixel at a time
			peak = (where(max(working[c,*,*]) EQ working[c,*,*]))[0]
			; Amount of flux at peak
			peak_value = (working[c,*,*])[peak]
; 			if l mod 100 EQ 0 then begin
; 				print,'Channel '+strtrim(string(c),2)+', Iteration '+strtrim(string(l),2)+ $
; 					': peak of '+strtrim(string(peak_value),2)+' at '+strtrim(string(peak),2)
; 			endif
	
			; Check to see if we can leave yet.  This is to test for the
			; flux/sigma limits
			if peak_value LT (Flux+Sigma*MapRms) then begin
				exit_status[c-C_Range[0]] = 1+l
				goto, CleanDone	
			endif
	
			; Convert from 1D -> 2D index
			peak_x = peak mod x_size
			peak_y = peak / x_size
	
			; Setup boundaries of removal to handle the edges.
			work_x_lo = max([(peak_x - b_size),0])
			work_x_hi = min([(peak_x + b_size),x_size-1])
			work_y_lo = max([(peak_y - b_size),0])
			work_y_hi = min([(peak_y + b_size),y_size-1])
			beam_x_lo = b_size - (peak_x - work_x_lo)
			beam_x_hi = b_size - (peak_x - work_x_hi)
			beam_y_lo = b_size - (peak_y - work_y_lo)
			beam_y_hi = b_size - (peak_y - work_y_hi)
		
			; Scale beam and fix boundaries to find what we need to remove
			to_remove = Gain*peak_value *  reform(beams[peak_x, peak_y, beam_x_lo:beam_x_hi, beam_y_lo:beam_y_hi])
		
			; Take that amount out of the working data array.  
			working[c, work_x_lo:work_x_hi, work_y_lo:work_y_hi] = reform(working[c, work_x_lo:work_x_hi, work_y_lo:work_y_hi]) $
										- to_remove
			; And add it into the cleaned map
			;+ Modified 07/16/08 - Though about what happens at the edge.   
			cleand[c,peak_x, peak_y] = cleand[c,peak_x, peak_y] + Gain*peak_value * b_norm[peak_x, peak_y]
		endfor
	
		; Exit point if needed	
		CleanDone:

	endfor
endelse

; Print out the results of the deconvolution all at once (for either method)
for c=C_Range[0],C_Range[1] do begin
	resid_rms = stddev(working[c,*,*])
	last_peak = max(working[c,*,*])
	if exit_status[c-C_Range[0]] NE -1 then exit_str='exiting on flux limit' else $
		exit_str='exiting on iteration limit'

	if Not Keyword_Set(Silent) then begin
		print,'>  Channel '+string(c,Format='(I4)')
		print,'>  '+exit_str
		print,'>   used '+strtrim(string((exit_status[c-C_Range[0]]-1)<NIter),2)+' iterations'
		print,'>   peak at loop exit '+strtrim(string(last_peak),2)+' mJy/beam'
		print,'>   residual RMS is '+strtrim(string(resid_rms),2)+' mJy/beam'
	endif
	Output = [Output, '>  Channel '+string(c,Format='(I4)'), '>  '+exit_str, $
		'>   used '+strtrim(string((exit_status[c-C_Range[0]]-1)<NIter),2)+' iterations', $
		'>   peak at loop exit '+strtrim(string(last_peak),2)+' mJy/beam', $
		'>   residual RMS is '+strtrim(string(resid_rms),2)+' mJy/beam']

endfor
; And convert exit_status codes to the standard type...
exit_status = (exit_status < 1)

; Restore the results to a gaussian beam with FWHM, uh, FWHM
if n_elements(FWHM) NE 2 then begin
	if n_elements(FWHM) EQ 1 then begin
		FWHM = replicate(FWHM,2)
	endif else begin
		FWHM = [3.1, 4.6]
	endelse
endif
if n_elements(PA) EQ 0 then PA = 0.0
PA = PA/!radeg
restore_sigma = double(FWHM) / 2.0 / sqrt( 2.0 * alog(2.0) )
if Not Keyword_Set(Silent) then begin
	print,'> restore'
	print,'>  using Gaussian with'+string(FWHM[0], Format='(F5.2)')+"'x"+$
	       string(FWHM[1], Format='(F5.2)')+"' FWHM"
endif
Output = [Output, '> restore', '>  using Gaussian with'+string(FWHM[0], Format='(F5.2)')+"'x"+string(FWHM[1], Format='(F5.2)')+"' FWHM"]

; Build restore kernel up
restore_beam = dblarr(21,21)
is = dblarr(21,21)
js = dblarr(21,21)
for i=0,20 do begin 
	is[i,*] = i-10
	js[*,i] = i-10
endfor
isp = is*cos(PA) + js*sin(PA)
jsp = -is*sin(PA) + js*cos(PA)

restore_beam = exp(-isp^2.0/(2.0*restore_sigma[0]^2.0) - jsp^2.0/(2.0*restore_sigma[1]^2.0))
restore_beam = restore_beam / total(restore_beam, /Double)

; Convolve with restore beam
for c=C_Range[0],C_Range[1] do begin
	temp = reform(cleand[c,*,*])
	cleand[c,*,*] = convolve(temp, restore_beam)
endfor
if Not Keyword_Set(Silent) then begin
	print,'>  clean map flux '+strtrim(string(total(cleand)/1000.0),2)+' Jy/beam'
	print,'>  residual flux  '+strtrim(string(total(working/1000.0)),2)+' Jy/beam'
endif
Output = [Output, '>  clean map flux '+strtrim(string(total(cleand)/1000.0),2)+' Jy/beam', $
	'>  residual flux  '+strtrim(string(total(working/1000.0)),2)+' Jy/beam']

; And add back in residuals
cleand = cleand + working
if Not Keyword_Set(Silent) then begin
	print,'>  total flux     '+strtrim(string(total(cleand)/1000.0),2)+' Jy/beam'
	print,'>  input flux     '+strtrim(string(total(dbox)/1000.0),2)+' Jy/beam'
endif
Output = [Output, '>  total flux     '+strtrim(string(total(cleand)/1000.0),2)+' Jy/beam', $
	'>  input flux     '+strtrim(string(total(dbox)/1000.0),2)+' Jy/beam']

; If we have summed channels together then resample back to the original velocity resolution
if KeyWord_Set(AllowSum) then begin
	cleand_old = cleand
	cleand = dbox_old
	for i=Range[0],Range[1] do $
		cleand[i, *, *] = cleand_old[(i-Range[0])/NComb, *, *] / double(NComb)
	dbox = dbox_old
endif

; Now for the Continuum map, if present
if n_elements(Continuum) NE 0 then begin
	if Not Keyword_Set(Silent) then begin
		print,'> clean'
		print,'>  cleaning continuum map to a peak flux of '+ $
			strtrim(string(Flux+Sigma*MapRms),2)+ ' mJy or '+ $
			strtrim(string(NIter),2)+' iterations'
	endif
	Output = [Output, '> clean', '>  cleaning continuum map to a peak flux of '+ $
	strtrim(string(Flux+Sigma*MapRms),2)+' mJy or '+ $
	strtrim(string(NIter),2)+' iterations']

	if Keyword_Set(CExts) then begin
		if Not Keyword_Set(Silent) then print,'>  using C library '+lib_version()
		Output = [Output, '>  using C library '+lib_version()]

			c_exit_status = lonarr( 1 ) - 2
			c_cleand2 = dblarr(x_size, y_size, 1)
			c_working2 = dblarr(x_size, y_size, 1)
			c_cleand2[*,*,0] = Continuum*0.0
			c_working2[*,*,0] = Continuum

			c_cleand2=alfa_clean_worker(c_cleand2, c_working2, beams, 0, 0, Flux+Sigma*MapRms, NIter, Gain, $
					b_norm, (size(beams))[3], x_size, 1, c_exit_status)
		
			c_cleand = reform(c_cleand2)
			c_working = reform(c_working2)

	endif else begin	
		; Setup the holding variables.  Working starts off the same as data but slowly has peaks 
		; removed.  cleand starts off as zeros and slowly has peaks added to it.  
		c_cleand = Continuum*0.0
		c_working = Continuum
		c_exit_status = -1
		
		for l=0L,(NIter-1) do begin
			; Select only one pixel at a time
			peak = (where(max(c_working) EQ c_working))[0]
			; Amount of flux at peak
			peak_value = c_working[peak]
	
			; Check to see if we can leave yet.  This is to test for the
			; flux/sigma limits
			if peak_value LE (Flux+Sigma*MapRms) then begin
				c_exit_status[0] = 1+l
				goto, C_CleanDone	
			endif
	
			; Convert from 1D -> 2D index
			peak_x = peak mod x_size
			peak_y = peak / x_size
	
			; Setup boundaries of removal to handle the edges.
			work_x_lo = max([(peak_x - b_size),0])
			work_x_hi = min([(peak_x + b_size),x_size-1])
			work_y_lo = max([(peak_y - b_size),0])
			work_y_hi = min([(peak_y + b_size),y_size-1])
			beam_x_lo = b_size - (peak_x - work_x_lo)
			beam_x_hi = b_size - (peak_x - work_x_hi)
			beam_y_lo = b_size - (peak_y - work_y_lo)
			beam_y_hi = b_size - (peak_y - work_y_hi)
		
			; Scale beam and fix boundaries to find what we need to remove
			to_remove = reform(gain * peak_value * beams[peak_x, peak_y, beam_x_lo:beam_x_hi, beam_y_lo:beam_y_hi])
		
			; Take that amount out of the working data array.  We don't update working 
			; and cleand until we know that removing this won't push total(working) < 0.
			temp = c_working
			temp[work_x_lo:work_x_hi, work_y_lo:work_y_hi] = $
				temp[work_x_lo:work_x_hi, work_y_lo:work_y_hi] - to_remove
			
			; We made it this far so update working
			c_working = temp
			; And add it into the cleaned map
			c_cleand[peak_x, peak_y] = c_cleand[peak_x, peak_y] + $
				total(gain * peak_value * beams[peak_x, peak_y, *, *])
		endfor

		; Exit point if needed	
		C_CleanDone:

	endelse

	; Report of exit status and progress
	c_resid_rms = stddev(c_working)
	c_last_peak = max(c_working)
	if c_exit_status[0] NE -1 then c_exit_str='exiting on flux limit' else $
		c_exit_str='exiting on iteration limit'

	if Not Keyword_Set(Silent) then begin
		print,'>  Continuum Map'
		print,'>  '+c_exit_str
		print,'>   used '+strtrim(string((c_exit_status[0]-1)<NIter),2)+' iterations'
		print,'>   peak at loop exit '+strtrim(string(c_last_peak),2)+' mJy/beam'
		print,'>   residual RMS is '+strtrim(string(c_resid_rms),2)+' mJy/beam'
	endif
	Output = [Output, '>  Continuum Map', '>  '+c_exit_str, $
		'>   used '+strtrim(string((c_exit_status[0]-1)<NIter),2)+' iterations', $
		'>   peak at loop exit '+strtrim(string(c_last_peak),2)+' mJy/beam', $
		'>   residual RMS is '+strtrim(string(c_resid_rms),2)+' mJy/beam']
	; And convert exit_status codes to the standard type...
	exit_status = (exit_status < 1)

	if Not Keyword_Set(Silent) then begin
		print,'> restore'
		print,'>  using Gaussian with'+string(FWHM[0], Format='(F5.2)')+"'x"+string(FWHM[1], Format='(F5.2)')+"' FWHM"
	endif
	Output = [Output, '> restore', '>  using Gaussian with'+string(FWHM[0], Format='(F5.2)')+"'x"+string(FWHM[1], Format='(F5.2)')+"' FWHM"]

	temp = c_cleand
	c_cleand = convolve(temp, restore_beam)

	if Not Keyword_Set(Silent) then begin
		print,'>  clean continuum map flux '+strtrim(string(total(c_cleand)/1000.0),2)+' Jy/beam'
		print,'>  residual flux  '+strtrim(string(total(c_working/1000.0)),2)+' Jy/beam'
	endif
	Output = [Output, '>  clean continuum map flux '+strtrim(string(total(c_cleand)/1000.0),2)+' Jy/beam', $
		'>  residual flux  '+strtrim(string(total(c_working/1000.0)),2)+' Jy/beam']
	
	; And add back in residuals
	c_cleand = c_cleand + c_working
	if Not Keyword_Set(Silent) then begin
		print,'>  continuum total flux     '+strtrim(string(total(c_cleand)/1000.0),2)+' Jy/beam'
		print,'>  continuum input flux     '+strtrim(string(total(Continuum)/1000.0),2)+' Jy/beam'
	endif
	Output = [Output, '>  continuum total flux     '+strtrim(string(total(c_cleand)/1000.0),2)+' Jy/beam', $
		'>  continuum input flux     '+strtrim(string(total(Continuum)/1000.0),2)+' Jy/beam']
endif

t_end = systime(1)
if Not Keyword_Set(Silent) then begin
	print,'> total time '+strtrim(string(t_end-t_start),2)+' seconds'
endif
Output = [Output, '> total time '+strtrim(string(t_end-t_start),2)+' seconds']

end
