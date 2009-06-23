pro plot_beams, file, Scale=Scale, PSFileBase=PSFileBase

restore,file

if n_elements(Scale) EQ 0 then Scale=1.0

beam0 = reform( ant_pat[0,*,*] )
beam1 = reform( ant_pat[1,*,*] )
beam2 = reform( ant_pat[2,*,*] )
beam3 = reform( ant_pat[3,*,*] )
beam4 = reform( ant_pat[4,*,*] )
beam5 = reform( ant_pat[5,*,*] )
beam6 = reform( ant_pat[6,*,*] )

l0 = 10.0*alog10(beam0/max(beam0))
l1 = 10.0*alog10(beam1/max(beam1))
l2 = 10.0*alog10(beam2/max(beam2))
l3 = 10.0*alog10(beam3/max(beam3))
l4 = 10.0*alog10(beam4/max(beam4))
l5 = 10.0*alog10(beam5/max(beam5))
l6 = 10.0*alog10(beam6/max(beam6))

c_levels = [-21, -18, -15, -12, -9, -6, -3]
c_labels = [' -21 ', ' -18 ', ' -15 ', ' -12 ', ' -9 ', ' -6 ', ' -3 ']
xrange = [0, n_elements(l0[*,0])-1]
yrange = [0, n_elements(l0[0,*])-1]

really_old_dev = !d.name
if n_elements(PSFileBase) NE 0 then begin
	file_name = PSFileBase+'-beam0.eps'

	Scale = 1
	set_plot,'ps'
	device, /Encapsulated, filename=file_name, $
		XSize=5, YSize=5, /Inches, Bits_Per_Pixel=8

	loadct, 0, /Silent
	TVLCT, r, g, b, /Get
	TVLCT, Reverse(r), Reverse(g), Reverse(b)
endif else begin
	set_plot,'x'
	window,0,XSize=200*Scale,YSize=200*Scale
endelse
tvimage, 255B-bytscl( congrid(l0, 200*Scale,200*Scale), min=-24, max=0), Pos=[0.00,0.00,1.00,1.00]
contour, l0, levels=c_levels, c_labels=c_labels, xrange=xrange, yrange=yrange, /NoErase, Pos=[0.00,0.00,1.00,1.00], XStyle=5, YStyle=5, Thick=2.0, CharThick=2.0
; xyouts, 0.50, 0.95, 'Beam 0', /Norm, CharThick=2.0, CharSize=1.0, Align=0.5
plots, [50, 110], [50, 50], Thick=3.0
xyouts, 80, 60, '3''', /Data, CharThick=2.0, CharSize=1.0, Align=0.5


if n_elements(PSFileBase) NE 0 then begin
	device, /close

	file_name = PSFileBase+'-beam123456.eps'

	Scale = 1
	device, /Encapsulated, filename=file_name, $
		XSize=9, YSize=6, /Inches, Bits_Per_Pixel=8
endif else begin
	window,1,XSize=600*Scale,YSize=400*Scale
endelse
tvimage, 255B-bytscl( congrid(l2, 200*Scale,200*Scale), min=-24, max=0), Pos=[0.00,0.50,0.33,1.00]
tvimage, 255B-bytscl( congrid(l1, 200*Scale,200*Scale), min=-24, max=0), Pos=[0.33,0.50,0.67,1.00]
tvimage, 255B-bytscl( congrid(l6, 200*Scale,200*Scale), min=-24, max=0), Pos=[0.67,0.50,1.00,1.00]

tvimage, 255B-bytscl( congrid(l3, 200*Scale,200*Scale), min=-24, max=0), Pos=[0.00,0.00,0.33,0.50]
tvimage, 255B-bytscl( congrid(l4, 200*Scale,200*Scale), min=-24, max=0), Pos=[0.33,0.00,0.67,0.50]
tvimage, 255B-bytscl( congrid(l5, 200*Scale,200*Scale), min=-24, max=0), Pos=[0.67,0.00,1.00,0.50]

contour, l2, levels=c_levels, c_labels=c_labels, xrange=xrange, yrange=yrange, /NoErase, Pos=[0.00,0.50,0.33,1.00], XStyle=5, YStyle=5, Thick=2.0, CharThick=2.0
contour, l1, levels=c_levels, c_labels=c_labels, xrange=xrange, yrange=yrange, /NoErase, Pos=[0.33,0.50,0.67,1.00], XStyle=5, YStyle=5, Thick=2.0, CharThick=2.0

xs = xrange[1]/2
plots, [xs-40, xs+40], [-5, -5], Thick=3.0
xyouts, xs, 5, '4''', /Data, CharThick=2.0, CharSize=1.0, Align=0.5

contour, l6, levels=c_levels, c_labels=c_labels, xrange=xrange, yrange=yrange, /NoErase, Pos=[0.67,0.50,1.00,1.00], XStyle=5, YStyle=5, Thick=2.0, CharThick=2.0

contour, l3, levels=c_levels, c_labels=c_labels, xrange=xrange, yrange=yrange, /NoErase, Pos=[0.00,0.00,0.33,0.50], XStyle=5, YStyle=5, Thick=2.0, CharThick=2.0
contour, l4, levels=c_levels, c_labels=c_labels, xrange=xrange, yrange=yrange, /NoErase, Pos=[0.33,0.00,0.67,0.50], XStyle=5, YStyle=5, Thick=2.0, CharThick=2.0
contour, l5, levels=c_levels, c_labels=c_labels, xrange=xrange, yrange=yrange, /NoErase, Pos=[0.67,0.00,1.00,0.50], XStyle=5, YStyle=5, Thick=2.0, CharThick=2.0

xyouts, 0.17, 0.95, 'Beam 2', /Norm, CharThick=2.0, CharSize=1.0, Align=0.5
xyouts, 0.50, 0.95, 'Beam 1', /Norm, CharThick=2.0, CharSize=1.0, Align=0.5
xyouts, 0.83, 0.95, 'Beam 6', /Norm, CharThick=2.0, CharSize=1.0, Align=0.5

xyouts, 0.17, 0.05, 'Beam 3', /Norm, CharThick=2.0, CharSize=1.0, Align=0.5
xyouts, 0.50, 0.05, 'Beam 4', /Norm, CharThick=2.0, CharSize=1.0, Align=0.5
xyouts, 0.83, 0.05, 'Beam 5', /Norm, CharThick=2.0, CharSize=1.0, Align=0.5

if n_elements(PSFileBase) NE 0 then begin
	device, /close
endif
set_plot,really_old_dev

end