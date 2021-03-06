;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro figure2d
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pol_file='pol_products/pol_8_10_R_cstn_ESO.fits'
stn_file='pol_products/stn_8_10_R_cstn_ESO.fits'

pol=mrdfits(pol_file)
stn=mrdfits(stn_file)
intens=mrdfits('stacked/rebin/I_all_2dstack2_6514_2d_rebin.fits')
ss=size(pol)

intens1d=mrdfits('stacked/I_all_2dstack2_6514_1d.fits',0,h)
wlen=dindgen(ss[1])*sxpar(h,'CD1_1')+sxpar(h,'CRVAL1')
x=dindgen(ss[1])
xrange=where(wlen gt 4930. and wlen lt 5057.5)

bin=strsplit(pol_file,'_',/extract)
name=strsplit(pol_file,'/.',/extract)

stncut = 3.0


;stn=mrdfits('pola_stn_airm.fits')
;Pola=Pobs*sqrt(1. - (error/Pobs)^2.d)
nan=finite(pol,/nan)
pol[where(nan eq 1)]=0.
pol_sig=pol*0. -100
pol_sig[where(stn gt stncut and stn lt 5.)]=pol[where(stn gt stncut and stn lt 5.)]
;pol_sig[where(pol_sig gt 0.35)]=0.

pol_sig[0:9,*]=0.
pol_sig[15:24,*] = 0.
;stop
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

set_plot,'ps'
device,/encap,filename='fig2d_'+name[1]+'_'+string(stncut,format='(D3.1)')+'_.ps',/color

pos1=[0.15,0.4,0.9,0.65]
pos2=[0.15,0.65,0.9,0.9]

x=dindgen(ss[1])*0.636*fix(bin[2]) + 4930.-fix(bin[2]/2.) ;wlen(xrange[0]-fix(bin[2]/2));
y=dindgen(ss[1])/35.*0.25*13.*6

;sc=dindgen(10)/10.
sc=[0.0,.05,.1,.15,.2,.25,.3]
sci=dblarr(7,3)

sci[*,0]=sc
sci[*,1]=sc
sci[*,2]=sc

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

loadct,0
; this draws IMAGES ONLY.
cgloadct,8,clip=100,/reverse;, ncolors=100
tvimage,bytscl(pol_sig,min=0.,max=.5),POSITION=pos1
loadct,0
tvimage,bytscl(-intens,min=-0.004,max=0.0001),POSITION=pos2;,/noerase
loadct,0

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; this draws AXES ONLY.
plot,x,y,POSITION=pos1,/noerase,/nodata,thick=4,charthick=4,xthick=4,ythick=4,xtitle='Wavelength [Angstroms]',charsize=1.5,xstyle=1,yr=[0.001,19.99],ystyle=1;,xr=[4771.,4883.+5.*0.64]+ 156.360
;oplot,dindgen(100)*0. + 4972.09 ,dindgen(100),thick=5,color=1
oplot,dindgen(100)*0. + 4981.5,dindgen(100),thick=5,linestyle=2,color=1

xyouts,5010,5,'Polarization',charsize=1.2,charthick=4,/data
xyouts,5010,2,'S/N>3.0',charsize=1.2,charthick=4,/data

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

plot,x,y,POSITION=pos2,/noerase, /nodata,thick=4,charthick=4,xthick=4,ythick=4,charsize=1.5,xstyle=1,xtickformat='(a2)',yr=[0.001,19.99],ystyle=1, title='Bin size = '+string(bin[2])+' x '+string(bin[3])

xyouts,5010,2,'Intensity',charsize=1.2,charthick=4,/data
xyouts,4915,-10,'Arcseconds',charsize=1.5,charthick=4,/data,orientation=90

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

cgloadct,8,clip=100,/reverse
tvimage,bytscl(sci,min=0.,max=.5),POSITION=[0.15,0.12,0.9,0.2];,/noerase
loadct,0

plot,sc,sc/sc,POSITION=[0.15,0.12,0.9,0.2],/noerase,/nodata,thick=4,xtitle='P fraction',charsize=1.2,xstyle=1,ystyle=4,charthick=4,ytickformat='(a2)',xr=[0,.35]

device,/close
set_plot,'X'
;spawn,'ps2pdf'figure_airm',1
stop
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro figure1d
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; Prepare the total intensity line plot
tot1d=mrdfits('tot_i_calint_1d.fits',0,h)
ss=size(tot1d)
wlen=dindgen(ss[1])*sxpar(h,'CD1_1')+sxpar(h,'CRVAL1')
x=dindgen(ss[1])
xrange=where(wlen gt 4930. and wlen lt 5057.5)

; mask out the ly-a line 
mask = where((wlen lt 4960 and wlen gt 4700) or (wlen gt 5020 and wlen lt 5500))
; fit the background with a linear polynomial
coeffs = poly_fit(x[mask],tot1d[mask],1)
y = coeffs[0] + coeffs[1]*x
; subtract the linear fit from the spectrum
flat = tot1d - y
;plot,wlen,flat,xr=[4900,5050]

; prepare the polarization fraction points
pol_file1='pol_products/pol_8_R_cstn_Tin.fits'
stn_file1='pol_products/stn_8_R_cstn_Tin.fits'
pol_file2='pol_products/pol_8_R_cstn_ESO.fits'
stn_file2='pol_products/stn_8_R_cstn_ESO.fits'

bin=strsplit(pol_file1,'_',/extract)
name=strsplit(pol_file1,'/.',/extract)

pol_t=mrdfits(pol_file1)
nan=finite(pol_t,/nan)
pol_t[where(nan eq 1)]=0.

pol_e=mrdfits(pol_file2)
nan=finite(pol_e,/nan)
pol_e[where(nan eq 1)]=0.
		
stn_t=mrdfits(stn_file1)
stn_e=mrdfits(stn_file2)

stncut=3.0

;if bin[2] eq 't' then spec_type='sum' else spec_type='mean'
ss2=size(pol_t)
xp=dindgen(ss2[1])*0.636*fix(bin[2]) + wlen(xrange[0]-fix(bin[2]/2));
pol_t2=pol_t*0. - 100
pol_t2[where(stn_t gt stncut)]=pol_t[where(stn_t gt stncut)]

pol_e2=pol_e*0.-100
pol_e2[where(stn_e gt stncut)]=pol_e[where(stn_e gt stncut)]

;pol2[0:ss2[1]/3-1]=-100.
;pol2[ss2[1]*2/3+1:ss2[1]-1]=0.
;stop

; create plot of 1d total intensity overlaid with fractional polarization in bins
set_plot,'ps'
device,/encap,filename='fig1d_'+name[1]+'_'+string(stncut,format='(D3.1)')+'_comp.ps',/color
pos1=[0.15,0.2,0.85,0.85]
; plot the actual 1d total intensity spectrum
cgplot,wlen,flat,POSITION=pos1,/noerase,thick=4,charthick=4,xthick=4,ythick=4,xtitle='Wavelength [Angstroms]',charsize=1.5,xstyle=1,ystyle=4,xr=[4930,5057.5],yr=[-0.,0.28],title='Bin size = '+bin[2]+' pixels'

; overplot the bin lines so we can see what we've summed over
num_x = n_elements(xrange)/fix(bin[2])
for i =0,num_x-1 do begin
;	oplot,dindgen(100)*0. + wlen(xrange(i*bin[2])),dindgen(100),linestyle=5
;stop
endfor

; overplot the line center given z = 3.09
oplot,dindgen(200)*0. + 4981.5 ,dindgen(200),thick=5,linestyle=2;,color=cgcolor('blue')

; draw the left-hand y-axis -- intensity axis
axis,4930,yaxis=0,yrange=[0,1],ystyle=1,charthick=4,charsize=1.5,ythick=4,ytitle='Intensity [arbitrary units]',/save

; draw the right-hand y-axis -- polarization fraction axis
axis,5057.5,yaxis=1,yrange=[0,0.4],ystyle=1,charthick=4,charsize=1.5,ythick=4,ytitle='Polarization Fraction',color=cgcolor('dark green'),/save
;xyouts,5028,.35,'STN>'+string(stncut,format='(D3.1)'),charsize=1.2,charthick=4,/data,color=cgcolor('dark green')

; overplot the polarization fraction 
circsym

perr_t=pol_t2/stn_t
perr_e=pol_e2/stn_e

perr_t[where(perr_t lt 0.)]=0.
perr_t[where(perr_t gt 1.)]=0.

perr_e[where(perr_e lt 0.)]=0.
perr_e[where(perr_e gt 1.)]=0.

;oplot,xp,pol2,linestyle=6,psym=8,symsize=4,color=cgcolor('dark green')
oploterror,xp,pol_t2,perr_t,linestyle=6,psym=8,symsize=4,errthick=3,color=cgcolor('blue')
oploterror,xp,pol_e2,perr_e,linestyle=6,psym=8,symsize=4,errthick=3,color=cgcolor('dark green')
xyouts,5020,.35,'P calc method:',charsize=1.2,charthick=4,/data;,color=cgcolor('blue')
xyouts,5025,.32,'Tinbergen',charsize=1.2,charthick=4,/data,color=cgcolor('blue')
xyouts,5025,0.30,'ESO',charsize=1.2,charthick=4,/data,color=cgcolor('dark green')

device,/close
set_plot,'X'

stop
end

pro full_ccd

sci=mrdfits('N2_1/N2_1_00_bftcor_crc.fits',0,h)

set_plot,'ps'
device,/encap,filename='example_science.ps',/color

pos1=[0.03,0.1,0.78,0.9]
tvimage,bytscl(sci,min=0.,max=0.08),POSITION=pos1,/keep_aspect_ratio
xyouts,.80,.687,'Sky',charthick=3,charsize=1.2
xyouts,.80,.595,'Align Star 1',charthick=3,charsize=1.2
xyouts,.80,.50,'Sky',charthick=3,charsize=1.2
xyouts,.80,.41,'LAB1',charthick=3,charsize=1.2
xyouts,.80,.31,'Align Star 2',charthick=3,charsize=1.2
device,/close
set_plot,'X'

stop
end

pro total_int

tot=mrdfits('stacked/I_all_2dstack2_6514.fits')
sig=mrdfits('stacked/I_all_2dstack2_6514_sig.fits')


set_plot,'ps'
device,/encap,filename='fringy.ps',/color


pos1=[0.2,0.4,0.8,0.5]
pos2=[0.2,0.65,0.8,0.9]

tvimage,bytscl(-tot,min=-0.004,max=0.),POSITION=pos2,/keep_aspect_ratio
tvimage,bytscl(sig,min=0.,max=0.008),POSITION=pos1,/keep_aspect_ratio


device,/close
set_plot,'X'


end










