pro a1963_patch

common agcshare

inWidth = [21, 25, 31, 35, 41, 45, 51]
for w=0,(n_elements(inWidth)-1) do begin
	cWidth = inWidth[w]
	case cWidth of 
		21:	restore,agcdir+'../deconvolution/beams/a1963_21_180.sav'
		25:	restore,agcdir+'../deconvolution/beams/a1963_25_180.sav'
		31:	restore,agcdir+'../deconvolution/beams/a1963_31_180.sav'
		35:	restore,agcdir+'../deconvolution/beams/a1963_35_180.sav'
		41:	restore,agcdir+'../deconvolution/beams/a1963_41_180.sav'
		45:	restore,agcdir+'../deconvolution/beams/a1963_45_180.sav'
		51:	restore,agcdir+'../deconvolution/beams/a1963_51_180.sav'
		else:	restore,agcdir+'../deconvolution/beams/a1963_25_180.sav'
	endcase

	inBeams = [1, 4]
	for b=0,(n_elements(inBeams)-1) do begin
		if inBeams[b] EQ 0 then $
			continue
		new_beam = fixbeam(inBeams[b], Width=cWidth, InputBeams=ant_pat, /Verbose)
		ant_pat[inBeams[b],*,*] = new_beam
	endfor

	new_name = 'a1963_galfa_'+strtrim(string(cWidth),2)+'_180.sav'
	save,ant_pat,filename=new_name
endfor

end
