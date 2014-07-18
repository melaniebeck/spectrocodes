pro polcomp
;	spawn,'ls *_fitsub > fitsub.list'
;	spawn,'ls *_nofitsub > nofitsub.list'

	readcol,'fitsub.list',format='a', fitsubfiles
	readcol,'nofitsub.list',format='a', nofitsubfiles

	p_fs = fltarr(n_elements(fitsubfiles),25)
	p_nfs = fltarr(n_elements(fitsubfiles),25)
	pe_fs = fltarr(n_elements(fitsubfiles),25)
	pe_nfs = fltarr(n_elements(fitsubfiles),25)

	readcol,fitsubfiles[0], a0, b0
	p_fs[0,*]=a0
	pe_fs[0,*]=b0
	readcol,nofitsubfiles[0], a0, b0
	p_nfs[0,*] = a0
	pe_nfs[0,*] = b0
	for i=1,n_elements(fitsubfiles)-1 do begin
		readcol,fitsubfiles[i], a1, b1
		p_fs[i,*] = a1
		pe_fs[i,*] = b1
		readcol,nofitsubfiles[i], a1, b1 
		p_nfs[i,*] = a1
		pe_nfs[i,*] = b1
	endfor
	x=dindgen(25)
	set_plot,'PS'
	device,filename='polcomp_multidata.ps',/color,/encap
	plot,x, p_fs[0,*], psym=1,yr=[0,.2], xr=[5,20], ytit='Fractional Polarization', $
		 xtit='Bin Number (25 bins total)', $
	 	 tit='Pol comparison from 5 sets of reduced spectra', charsize=1
	oplot,x, p_fs[1,*], psym=2
	oplot,x, p_fs[2,*], psym=7
	oplot,x, p_fs[3,*], psym=4
	oplot,x, p_fs[4,*], psym=5
	oplot,x, p_fs[5,*], psym=6
	oplot,x, p_nfs[0,*], psym=1, color=cgcolor('red')
	oplot,x, p_nfs[1,*], psym=2, color=cgcolor('red')
	oplot,x, p_nfs[2,*], psym=7, color=cgcolor('red')
	oplot,x, p_nfs[3,*], psym=4, color=cgcolor('red')
	oplot,x, p_nfs[4,*], psym=5, color=cgcolor('red')
	oplot,x, p_nfs[5,*], psym=6, color=cgcolor('red')
	xyouts,10,.18,'Linear fit subtracted from 1D spectra',color=cgcolor('red')
	xyouts,10,.17,'No fit subtracted from 1D spectra'
	device,/close
	set_plot,'X'
stop
end
