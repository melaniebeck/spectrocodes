pro figure

readcol,'data_lists/hmags_3d.list',hmag3d
readcol,'data_lists/hmags_20k.list',hmag20k
readcol,'data_lists/hmags_deep2.list',hmagdeep2

tdhst=histogram_ez_cla(hmag3d,binsize=0.5)
t0k=histogram_ez_cla(hmag20k,binsize=0.5)
tdeep2=histogram_ez_cla(hmagdeep2,binsize=0.5)

set_plot,'ps'
device,filename='hmags3.ps',/encap,/color

plot,tdhst[0,*],tdhst[1,*]/total(tdhst[1,*]),psym=10,yr=[0,0.4],thick=4,xthick=4,ythick=4,xtitle='H [AB magnitudes]',ytitle='Fraction',charsize=1.5,charthick=4,/nodata

loadct,13

oplot,t0k[0,*],t0k[1,*]/total(t0k[1,*]),psym=10,thick=4,color=cgcolor('red')
;oplot,tdeep2[0,*],tdeep2[1,*]/total(tdeep2[1,*]),psym=10,thick=4,color=cgcolor('red')
oplot,tdhst[0,*],tdhst[1,*]/total(tdhst[1,*]),psym=10,thick=4,color=cgcolor('blue')



; plot details for COS 20k data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
x=dindgen(2*n_elements(reform(t0k[0,*])))
y=dindgen(2*n_elements(reform(t0k[0,*])))

for j=0,n_elements(x)-1,2 do begin

   x[j]=t0k[0,j/2]-0.25
   x[j+1]=t0k[0,j/2]+0.25
   
   y[j]=t0k[1,j/2]/total(t0k[1,*])
   y[j+1]=t0k[1,j/2]/total(t0k[1,*])
endfor
polyfill,x,y,thick=3,color=cgcolor('red'),/line_fill,orientation=45

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; plot details for DEEP2 data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;x=dindgen(2*n_elements(reform(tdeep2[0,*])))
;y=dindgen(2*n_elements(reform(tdeep2[0,*])))

;for j=0,n_elements(x)-1,2 do begin
;   x[j]=tdeep2[0,j/2]-0.25
;   x[j+1]=tdeep2[0,j/2]+0.25
;   
;   y[j]=tdeep2[1,j/2]/total(tdeep2[1,*])
;   y[j+1]=tdeep2[1,j/2]/total(tdeep2[1,*])
;endfor
;polyfill,x,y,/line_fill,thick=3,orientation=-45,color=45

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; plot details for 3dhst data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
x=dindgen(2*n_elements(reform(tdhst[0,*])))
y=dindgen(2*n_elements(reform(tdhst[0,*])))

for j=0,n_elements(x)-1,2 do begin
   x[j]=tdhst[0,j/2]-0.25
   x[j+1]=tdhst[0,j/2]+0.25
   
   y[j]=tdhst[1,j/2]/total(tdhst[1,*])
   y[j+1]=tdhst[1,j/2]/total(tdhst[1,*])
endfor
polyfill,x,y,thick=3,color=cgcolor('blue'),/line_fill,orientation=-45

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; matching legend detail for 3DHST data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
xx=[24,25,25,24,24]-1
yy=[0.35,0.35,0.37,0.37,0.35]
oplot,xx,yy,thick=4,color=cgcolor('blue')
polyfill,xx,yy,thick=3,color=cgcolor('blue'),/line_fill,orientation=45

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; matching legend detail for COS 20k data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
xx=[24,25,25,24,24]-1
yy=[0.32,0.32,0.34,0.34,0.32]
oplot,xx,yy,thick=4,color=cgcolor('red')
polyfill,xx,yy,/line_fill,thick=3,orientation=45,color=cgcolor('red')

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; matching legend detail for DEEP2 data
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;xx=[24,25,25,24,24]-1
;yy=[0.29,0.29,0.31,0.31,0.29]
;oplot,xx,yy,thick=4,color=45
;polyfill,xx,yy,/line_fill,thick=3,orientation=45,color=45





xyouts,24.3,0.353,'3DHST members',charthick=4,charsize=1.2
xyouts,24.3,0.323,'20k members',charthick=4,charsize=1.2
;xyouts,24.3,0.293,'DEEP2 members',charthick=4,charsize=1.2

Device,/close
;makepdf,'mags',1
stop

end
