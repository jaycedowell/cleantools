function heiles, b, Width=Width

coeff = [0.0, 1.0, 2.0, 3.0, $
	-4.0, -3.0, -2.0, -1.0]
H = complexarr(7, 8)
C = complexarr(7, 8)
W = complexarr(7, 8)

; Feed 0
H[0,*] = [complex(0.0272, 0.0000), complex(-0.0036, 0.0024), complex(0.0006, 0.0006), complex(0.006, -0.0007), $
	complex(0.0023, 0.0000), complex(0.0006, 0.0007), complex(0.0006, -0.0006), complex(-0.0036, -0.0024)]
C[0,*] = [complex(5.7181, 0.0000), complex(-0.0085, 0.0383), complex(-0.1241, 0.0426), complex(0.0000, 0.0000), $
	complex(0.0000, 0.0000), complex(0.0000, 0.0000), complex(-0.1241, -0.0426), complex(-0.0085, -0.0383)]
W[0,*] = [complex(1.7305, 0.0000), complex(0.0072, 0.0355), complex(-0.0267, 0.0750), complex(0.0000, 0.0000), $
	complex(0.0000, 0.0000), complex(0.0000, 0.0000), complex(-0.0267, -0.0750), complex(0.0072, -0.0355)]

; Feed 1
H[1,*] = [complex(0.0487, 0.0000), complex(-0.0202, -0.0232), complex(0.0003, 0.0088), complex(-0.0016, -0.0041), $
	complex(-0.0016, 0.0000), complex(-0.0016, 0.0041), complex(0.0003, -0.0088), complex(-0.0202, 0.0232)]
C[1,*] = [complex(5.5117, 0.0000), complex(-0.0241, 0.2027), complex(0.0027, -0.0340), complex(0.0000, 0.0000), $
	complex(0.0000, 0.0000), complex(0.0000, 0.0000), complex(0.0027, 0.0340), complex(-0.0241, -0.2027)]
W[1,*] = [complex(2.0863, 0.0000), complex(-0.0442, -0.0899), complex(-0.0490, 0.1648), complex(0.0000, 0.0000), $
	complex(0.0000, 0.0000), complex(0.0000, 0.0000), complex(-0.0490, -0.1648), complex(-0.0442, 0.0899)]

; Feed 2
H[2,*] = [complex(0.0558, 0.0000), complex(-0.0374, 0.0060), complex(0.0141, -0.0059), complex(-0.0094, 0.0019), $
	complex(0.0096, 0.0000), complex(-0.0094, -0.0019), complex(0.0141, 0.0059), complex(-0.0374, -0.0060)]
C[2,*] = [complex(5.5310, 0.0000), complex(0.1441, 0.1763), complex(-0.1761, -0.0222), complex(0.0000, 0.0000), $
	complex(0.0000, 0.0000), complex(0.0000, 0.0000), complex(-0.1761, 0.0222), complex(0.1441, -0.1763)]
W[2,*] = [complex(2.0474, 0.0000), complex(-0.0999, 0.0215), complex(0.1120, -0.0039), complex(0.0000, 0.0000), $
	complex(0.0000, 0.0000), complex(0.0000, 0.0000), complex(0.1120, 0.0039), complex(-0.0999, -0.0215)]

; Feed 3
H[3,*] = [complex(0.0556, 0.0000), complex(-0.0184, 0.0329), complex(-0.0004, -0.0107), complex(-0.0019, 0.0111), $
	complex(-0.0086, 0.0000), complex(-0.0019, -0.0111), complex(-0.0004, 0.0107), complex(-0.0184, -0.0329)]
C[3,*] = [complex(5.4427, 0.0000), complex(0.1399, 0.0350), complex(-0.1273, 0.1049), complex(0.0000, 0.0000), $
	complex(0.0000, 0.0000), complex(0.0000, 0.0000), complex(-0.1273, -0.1049), complex(0.1399, -0.0350)]
W[3,*] = [complex(2.1206, 0.0000), complex(-0.0534, 0.0458), complex(-0.0686, -0.0366), complex(0.0000, 0.0000), $
	complex(0.0000, 0.0000), complex(0.0000, 0.0000), complex(-0.0686, 0.0366), complex(-0.0534, -0.0458)]

; Feed 4
H[4,*] = [complex(0.0487, 0.0000), complex(0.0202, 0.0232), complex(0.0003, 0.0088), complex(0.0016, 0.0041),$
	complex(-0.0016, 0.0000), complex(0.0016, -0.0041), complex(0.0003, -0.0088), complex(0.0202, -0.0232)]
C[4,*] = [complex(5.5117, 0.0000), complex(0.0241, -0.2027), complex(0.0027, -0.0340), complex(-0.0000, -0.0000), $
	complex(0.0000, 0.0000), complex(-0.0000, -0.0000), complex(0.0027, 0.0340), complex(0.0241, 0.2027)]
W[4,*] = [complex(2.0863, 0.0000), complex(0.0442, 0.0899), complex(-0.0490, 0.1648), complex(-0.0000, -0.0000), $
	complex(0.0000, 0.0000), complex(-0.0000, -0.0000), complex(-0.0490, -0.1648), complex(0.0442, -0.0899)]

; Feed 5
H[5,*] = [complex(0.0496, 0.0000), complex(0.0305, 0.0006), complex(0.0115, 0.0031), complex(0.0097, 0.0019), $
	complex(0.0094, 0.0000), complex(0.0097, -0.0019), complex(0.0115, -0.0031), complex(0.0305, -0.0006)]
C[5,*] = [complex(5.5929, 0.0000), complex(-0.1671, -0.0097), complex(-0.1557, 0.0741), complex(0.0000, 0.0000), $
	complex(0.0000, 0.0000), complex(0.0000, 0.0000), complex(-0.1557, -0.0741), complex(-0.1671, 0.0097)]
W[5,*] = [complex(2.0248, 0.0000), complex(0.0732, -0.0519), complex(0.1006, 0.0708), complex(0.0000, 0.0000), $
	complex(0.0000, 0.0000), complex(0.0000, 0.0000), complex(0.1006, -0.0708), complex(0.0732, 0.0519)]

; Feed 6
H[6,*] = [complex(0.0474, 0.0000), complex(0.0138, -0.0275), complex(-0.0017, -0.0093), complex(0.0005, -0.0021), $
	complex(0.0052, 0.0000), complex(0.0005, 0.0021), complex(-0.0017, 0.0093), complex(0.0138, 0.0275)]
C[6,*] = [complex(5.4416, 0.0000), complex(-0.0854, 0.0345), complex(-0.0291, 0.0699), complex(0.0000, 0.0000), $
	complex(0.0000, 0.0000), complex(0.0000, 0.0000), complex(-0.0291, -0.0699), complex(-0.0854, -0.0345)]
W[6,*] = [complex(2.0517, 0.0000), complex(0.0758, -0.1084), complex(-0.0982, -0.0506), complex(0.0000, 0.0000), $
	complex(0.0000, 0.0000), complex(0.0000, 0.0000), complex(-0.0982, 0.0506), complex(0.0758, 0.1084)]

if n_elements(Width) NE 1 then $
	Width = 25.0
Width = fix( Width / 0.05 )		; arc min -> pixels
if Width mod 2 EQ 0 then $
	Width++
beam = dblarr(Width,Width)		; new array (Width+.05)'x(Width+0.05)' @ 0.05'/px

useH = reform(H[b,*])
useC = reform(C[b,*])
useW = reform(W[b,*])

Sampling = 3.0
theta = findgen(360*Sampling)/Sampling*!Pi/180.0

HGT = complexarr(n_elements(theta))
CEN = complexarr(n_elements(theta))
WID = complexarr(n_elements(theta))
for f=0,(n_elements(coeff)-1) do begin
	HGT += useH[f]*exp(complex(0.0, 1.0)*coeff[f]*theta)
	CEN += useC[f]*exp(complex(0.0, 1.0)*coeff[f]*theta)
	WID += useW[f]*exp(complex(0.0, 1.0)*coeff[f]*theta)
endfor

x_size = n_elements(beam[*,0])
y_size = n_elements(beam[0,*])
for i=0,(x_size-1) do begin
	for j=0,(y_size-1) do begin
		r = sqrt( (i-x_size/2)^2.0 + (j-y_size/2)^2.0 ) * 0.05
		t = atan( float(j-x_size/2), float(i-y_size/2))
		if t LT 0 then $
			t += 2.0*!Pi
		useHGT = interpol(HGT, theta, t)
		useCEN = interpol(CEN, theta, t)
		useSWID = interpol(WID, theta, t) / 2.0 / sqrt(alog(2.0))

		beam[i,j] = useHGT*exp( -((r-useCEN)/useSWID)^2.0 )
	endfor
endfor

return, beam

end



