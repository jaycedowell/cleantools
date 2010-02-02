function fixbeam, b, Width=Width, Verbose=Verbose, InputBeams=InputBeams

common agcshare

if n_elements(Width) NE 1 then $
	Width = 25.0
Width2 = Width

if n_elements(InputBeams) EQ 0 then begin
	case Width of 
		21:	restore,agcdir+'../deconvolution/beams/a1963_21_180.sav'
		25:	restore,agcdir+'../deconvolution/beams/a1963_25_180.sav'
		31:	restore,agcdir+'../deconvolution/beams/a1963_31_180.sav'
		35:	restore,agcdir+'../deconvolution/beams/a1963_35_180.sav'
		41:	restore,agcdir+'../deconvolution/beams/a1963_41_180.sav'
		45:	restore,agcdir+'../deconvolution/beams/a1963_45_180.sav'
		51:	restore,agcdir+'../deconvolution/beams/a1963_51_180.sav'
		else:	restore,agcdir+'../deconvolution/beams/a1963_25_180.sav'
	endcase
endif else begin
	ant_pat = InputBeams
endelse

if Keyword_Set(Verbose) then $
	print,"% Loading data for central feed"
acen = reform( ant_pat[0,*,*] )
gcen = heiles(0, Width=Width2)

x_size = n_elements(acen[*,0])
y_size = n_elements(acen[0,*])
r = acen*0.0
t = acen*0.0
for i=0,(x_size-1) do begin
	for j=0,(y_size-1) do begin
		r[i,j] = sqrt( (i-x_size/2)^2.0 + (j-y_size/2)^2.0 ) * 0.05
		tmp = atan( float(j-x_size/2), float(i-y_size/2))
		if tmp LT 0 then $
			tmp += 2*!Pi
		t[i,j] = tmp*180/!PI
	endfor
endfor

toShift = fltarr(9)
for i=0,(n_elements(toShift)-1) do begin
	valid = where( abs(r - (5.4+0.1*i)) LE 0.01 )
	a_r = acen[valid]
	a_t = t[valid]
	g_r = gcen[valid]
	g_t = t[valid]
	
	a_order = sort(a_t)
	a_r = a_r[a_order] / max(a_r)
	a_t = a_t[a_order]
	g_order = sort(g_t)
	g_r = g_r[g_order] / max(g_r)
	g_t = g_t[g_order]
	

	lags = findgen(720)/2.0 - 180.0
	cc = c_correlate(a_r, g_r, lags)
	best = min(where( cc EQ max(CC)))
	toShift[i] = lags[best]
endfor
toShift = median(toShift)
if Keyword_Set(Verbose) then $
	print,'% Best Rotation: '+string(toShift, Format='(F5.2)')+' degrees'
gcen = rot(gcen, toShift)

toScale = fltarr(10)
for i=0,(n_elements(toScale)-1) do begin
	valid = where( abs(r - (5.0+0.2*i)) LE 0.01 )
	a_r = acen[valid]
	g_r = gcen[valid]
	
	toScale[i] = median(a_r / g_r)
endfor
if Keyword_Set(Verbose) then $
	print,'% Best Scaling: '+string(10.0*alog10(median(toScale)), Format='(F5.1)')+' dB'
toScale = median(toScale)
gcen *= toScale

if Keyword_Set(Verbose) then begin
	window,0,XSize=1024,YSize=512,Title='Central Feed Comparison: A1963 and GALFA 2004-01'
	tvimage, bytscl(acen, min=0, max=1e-4), Pos=[0,0,0.5,1]
	tvimage, bytscl(gcen, min=0, max=1e-4), Pos=[0.5,0,1,1]
endif

if Keyword_Set(Verbose) then $
	print,"% Loading data for feed #"+strtrim(string(b),2)
a1963 = reform( ant_pat[b,*,*] )
Width2 = Width
galfa = heiles(b, Width=Width2)

if Keyword_Set(Verbose) then $
	print,"% Rotating and scaling GALFA sidelobe"
galfa = rot(galfa, toShift)
galfa *= toScale

toReplace = where( galfa GT 1.1*a1963, numReplace )

final = a1963
mask = fix(final*0.0)+1
if toReplace[0] NE -1 then begin
	final[toReplace] = galfa[toReplace]
	mask[toReplace] = [0]
endif

if Keyword_Set(Verbose) then $
	print,"% "+strtrim(string(numReplace),2)+" pixels replaced ("+ $
		strtrim(string(1.0*numReplace/n_elements(final)*100.0),2)+"%)"

if Keyword_Set(Verbose) then begin
	window,1,XSize=1024,YSize=512,Title='Feed '+strtrim(string(b),2)+': Before and After'
	tvimage, bytscl(10.0*alog10(a1963/max(a1963)), min=-27.0, max=0.0), Pos=[0,0,0.5,1]
	tvimage, bytscl(10.0*alog10(final/max(final)), min=-27.0, max=0.0), Pos=[0.5,0,1,1]

	;window,2,XSize=512,YSize=512,Title='Feed '+strtrim(string(b),2)+': Replacement Mask'
	;tvimage, bytscl(mask), Pos=[0,0,1,1]
endif

return,final

end


