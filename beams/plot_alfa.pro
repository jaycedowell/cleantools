pro plot_alfa, file, PSFileBase=PSFileBase

restore,file

if n_elements(Scale) EQ 0 then Scale=1.0

beam = ant_pat
nS = n_elements(beam[0,*,0])
scale = nS / 25.0 / 60.0

nM = round(24.0 / 25.0 * nS)
map = fltarr(nM,nM)
map2 = fltarr(7, nM, nM)

azA = double(329.060)
zaA = double(384.005)
az = double([0.000, -164.530, -329.060, -164.530,  164.530, 329.060, 164.530])
za = double([0.000,  332.558,    0.000, -332.558, -332.558,   0.000, 332.558])

azR = az
zaR = za

ang = 19.0
for i=0,6 do begin
	ang2 = -ang*!DPi/180.0
	ang3 = ((i-5)*60.0 - ang)*!DPi/180.0

	;out = [[cos(ang2), sin(ang2)], [-sin(ang2), cos(ang2)]] # [az[i], za[i]]
	out = dblarr(2)
	if i NE 0 then begin
		out[0] = azA*cos(ang3)
		out[1] = zaA*sin(ang3)
	endif
	print,i,za[i]/60.0,(out[1]-za[0])/60.0/2.1
	azR[i] = out[0]
	zaR[i] = out[1]

	rBeam = reform(beam[i,*,*])
	;rBeam = ( rot(rBeam, -ang) )

	dX = round( scale*out[0] + nM/2 )
	dY = round( scale*out[1] + nM/2 )
	
	mlx = max( [0, (dX-nS/2)] )
	mux = min( [nM-1, (dX+nS/2)] )
	mly = max( [0, (dY-nS/2)] )
	muy = min( [nM-1, (dY+nS/2)] )

	blx = nS/2 - (dX - mlx)
	bux = nS/2 - (dX - mux)
	bly = nS/2 - (dY - mly)
	buy = nS/2 - (dY - muy)

	sBeam = rBeam[blx:bux, bly:buy]
	if i EQ 0 then $
		sBeam *= 1.27

	map[mlx:mux, mly:muy] = map[mlx:mux, mly:muy] + sBeam
	map2[i, mlx:mux, mly:muy] = sBeam
endfor

really_old_dev = !d.name
if n_elements(PSFileBase) NE 0 then begin
	file_name = PSFileBase+'-alfa.eps'

	set_plot,'ps'
	device, /Encapsulated, filename=file_name, /Color, /CMYK, $
		XSize=5, YSize=5, /Inches, Bits_Per_Pixel=8

	loadct, 39, /Silent
endif else begin
	set_plot,'x'
	window,2,XSize=512,YSize=512
endelse

map = map / max(map)
tvimage,255B-bytscl(map, min=0.1, max=1.0), $
	Pos=[0.13, 0.10, 0.98, 0.95]

plot, [0,0], /NoData, XRange=[-12,12], YRange=[-12,12], $
	XStyle=1, YStyle=1, $
	XTitle='Azimuth Offset [arc minutes]', $
	YTitle='Zenith Angle Offset [arc minutes]', $
	Pos=[0.13, 0.10, 0.98, 0.95], /NoErase, $
	XThick=3.0, YThick=3.0, Thick=3.0, CharThick=3.0, CharSize=1.01

for i=0,6 do begin
	;plots, azR[i]/60.0, zaR[i]/60.0, PSym=1
	xyouts, azR[i]/60.0+1.0, zaR[i]/60.0+0.5, string(i,Format='(I1)'), /Data, CharThick=3.0, CharSize=1.3
	if i EQ 0 then begin
		plots, [-12,12], zaR[i]/60.0*[1,1], LineStyle=3, Thick=3.0, Color=FSC_Color('Green')
	endif else begin
		plots, [-12,12], zaR[i]/60.0*[1,1], LineStyle=2, Thick=3.0
	endelse
endfor

plots, [-12,12],zaR[0]/60.0*[1,1]-(7+18/60.), LineStyle=3, Thick=3.0, Color=FSC_Color('Green')

colorbar, Range=[1.0, 0.1], Divisions=3, Minor=2, TickNames=string(reverse(findgen(4)*30+10), Format='(I3, "%")'), $
	Pos=[0.70, 0.86, 0.93, 0.89], Title='Sensitivity', CharThick=3.0, CharSize=1.01

if n_elements(PSFileBase) NE 0 then begin
	device, /close

	file_name = PSFileBase+'-response.eps'

	device, /Encapsulated, filename=file_name, /Color, /CMYK, $
		XSize=5, YSize=5, /Inches, Bits_Per_Pixel=8
endif else begin
	window,0,XSize=512,YSize=512
endelse

colors = ['Red', 'Forest Green', 'Blue', 'Saddle Brown', 'Magenta', 'Orange']
peak = 1.0
for i=0,6 do begin
	cut = fltarr(nM)
	for j=0,(nM-1) do begin
		cut[j] = total( reform(map2[i,*,j]), /Double )
	endfor
	x = 1/scale*(findgen(nM)-nM/2)/60.0
	if i EQ 0 then $
		peak = max(cut)
	cut = cut / peak	


	if i EQ 0 then begin
		plot, x, cut, XRange=[-12,12], YRange=[0,1.05], $
			XStyle=1, YStyle=1, $
			XTitle='Zenith Angle Offset [arc min]', $
			YTitle='Relative Total Feed Contribution', $
			Pos=[0.12, 0.10, 0.97, 0.95], /NoErase, $
			XThick=3.0, YThick=3.0, Thick=3.0, CharThick=3.0, CharSize=1.01
	endif else begin
		oplot,x,cut,Color=FSC_Color(colors[i-1]), Thick=3.0
	endelse
endfor

legend,string(findgen(7),Format='("Feed ", I1)'),LineStyle=replicate(0,7), Color=[0B, FSC_Color(colors)], $
	Thick=replicate(3.0,7), CharThick=3.0, CharSize=0.85

if n_elements(PSFileBase) NE 0 then begin
	device, /close
endif
set_plot,really_old_dev

end
