pro beam_test, grid, Pass=Pass, Fail=Fail, NoExit=NoExit, Extended=Extended

Pass = 0
Fail = 0

Catch, Error_status

if Error_status NE 0 then begin
	print,'Error '+strtrim(string(Error_status),2)+': '+!Error_State.Msg
	if NOT Keyword_Set(NoExit) then $
		exit, status=abs(Error_status)
endif

if Keyword_Set(Extended) then BeamSizes = [21, 25, 31, 35, 41, 45, 51] else BeamSizes = [25] 

for b=0,(n_elements(BeamSizes)-1) do begin
	build_beam5, grid, [72,72], beam1, MapSize=BeamSizes[b], /Silent
	build_beam5, grid, [72,72], beam2, MapSize=BeamSizes[b], /CExts, /Silent
	
	beamr_total = 4.1021862397d-2
	beam1 = beam1 / total(beam1, /Double)
	beam2 = beam2 / total(beam2, /Double)
	beamd = beam1 - beam2
	
	if BeamSizes[b] EQ 25 then begin
		if abs(beamr_total - max(beam1)) GT 1.0d-9 then begin
			print,'Error: log(Differnce between expected and IDL method peaks) > -9'
			Fail = Fail + 1
	
			if NOT Keyword_Set(NoExit) then $
				exit, status=254
		endif else begin
			Pass = Pass + 1
		endelse
	endif
	
	if abs(max(beam1)-max(beam2)) GT 1.0d-9 then begin
		print,'Error: log(Differnce between IDL and C method peaks) > -9'
		Fail = Fail + 1

		if NOT Keyword_Set(NoExit) then $
			exit, status=254
	endif else begin
		Pass = Pass + 1
	endelse
	
	if abs(max(beamd)) GT 1.0d-9 OR abs(min(beamd)) GT 1.0d-9 then begin
		print,'Error: log(Differnce between IDL and C methods) > -9'
		Fail = Fail + 1

		if NOT Keyword_Set(NoExit) then $
			exit, status=254
	endif else begin
		Pass = Pass + 1
	endelse
endfor

end
