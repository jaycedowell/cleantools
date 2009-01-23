pro cleantest

restore,'/home/jddowell/DataCapacitor/ALFALFA/1236+15/gridbf_1236+15a.sav'
data = double(reform(total(grid.d[*,*,48:92,12:56],2)/2.0))
help,data

build_beam_grid2, grid, 48, 12, 92, 56, beams, MapSize=51, FWHM=FWHM, PA=PA, /Single

; Part 1: IDL Only
start1 = systime(1)
alfa_clean7, data, beams, cleand_temp, MapRMS=3.2d0, NIter=50000L, Flux=3.2d0, $
	Sigma=1.0d0, Gain=1.0d-1, Range=[0L,1023L], AllowSum=0, Exit_Status=out2, $
	FWHM=FWHM, PA=PA, /Silent, Output=out3
end1 = systime(1)

wait,2

; Part 2: IDL+C
start2 = systime(1)
alfa_clean7, data, beams, cleand_temp, MapRMS=3.2d0, NIter=50000L, Flux=3.2d0, $
	Sigma=1.0d0, Gain=1.0d-1, Range=[0L,1023L], AllowSum=0, Exit_Status=out2, $
	FWHM=FWHM, PA=PA, /Silent, Output=out3, /CExts
end2 = systime(1)

print,'IDL:   '+strtrim(string(end1-start1),2)+" s"
print,'IDL+C: '+strtrim(string(end2-start2),2)+" s"

end

