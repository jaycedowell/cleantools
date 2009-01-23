pro beamtest

restore,'/home/jddowell/DataCapacitor/ALFALFA/1020+07/gridbf_1020+07a.sav'

; Part 1: IDL only
start1 = systime(1)
;build_beam_grid2, grid, 70, 70, 79, 79, beams, MapSize=51
end1 = systime(1)

DelVarX, beams

; Part 2: IDL+C @ 4 cores
start2 = systime(1)
build_beam_grid2, grid, 70, 70, 79, 79, beams, MapSize=51, /CExts
end2 = systime(1)

print,'IDL:   '+strtrim(string(end1-start1),2)+" s"
print,'IDL+C: '+strtrim(string(end2-start2),2)+" s"

end

