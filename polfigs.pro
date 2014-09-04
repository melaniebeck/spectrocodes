;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro fig2d, MYSIG=mysig
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;pol_file='pol_products/pol_8_10_R_cstn_ESO.fits'
;stn_file='pol_products/stn_8_10_R_cstn_ESO.fits'

	;if KEYWORD_SET(MYSIG) then begin
	pol_file='pol_products/715/pol_715_mysig_8xV_bb8_resized.fits'
	stn_file='pol_products/715/MC_pstn95_715_mysig_8xV_bb8_resized.fits'
	err_file='pol_products/715/MC_sig95_715_mysig_8xV_bb8_resized.fits'
	;endif else begin
	;	pol_file='pol_products/'+date+'/pol_'+date+'_IRAF.fits'
	;	stn_file='pol_products/'+date+'/MC_pstn68_'+date+'_IRAF.fits'
	;endelse
	
	int_file='stacked/715/rebin/Iall_avg_2d_crop_rebin8xV_bb8_resized.fits'

	i00=mrdfits('stacked/715/rebin/I00_stack_715_2d_crop_rebin8xV_bb8_resized.fits')
	i00s=mrdfits('stacked/715/rebin/I00_stack_715_2d_MYsig_avg_rebin8xV_bb8_resized.fits')
	stn00 = i00/i00s

;	pol_file='pol_products/715/pol_vtbin_MYSIG_Iall_newavg_cut.fits'		
;	stn_file='pol_products/715/MC_pstn95_vtbin_MYSIG_Iall_newavg_cut.fits'

;	int_file='stacked/Iall_avg_crop.fits'
;	int_file='stacked/715/crop/Iall_stack_715_2d_sum2_crop.fits'

	int1d_file='stacked/715/Iall_stack_715_2d.fits'
	
	; read in appropriate files
	pol=mrdfits(pol_file)
	stn=mrdfits(stn_file)
	err=mrdfits(err_file)
	int=mrdfits(int_file)
	int1d=mrdfits(int1d_file,0,h)
	
	pstncut = 1.0
	stn00cut = 0.0
	ss=size(pol)
	ssi=size(int1d)

	; remove any nans
	nan=finite(pol,/nan)
	pol[where(nan eq 1)]=0.
	; ceate array where only pol above stn cut remains
	pol_sig=pol*0.-100
	err_sig = err*0.+1.	
	pol_sig[where(stn gt pstncut and stn00 gt stn00cut)]= $
		pol[where(stn gt pstncut and stn00 gt stn00cut)]
	pol_sig[where(pol_sig ge 0.5)]=0.

	; cut any bin in the error map that is greater than 0.4
	; P necessary to achieve S/N = 2  P=S/N*err=2*err
	err_sig[where(err lt 0.5)] = err[where(err lt 0.5)] 

	; page position of fractional polarization
	pos_p=[0.1,0.25,0.78,0.45]
	; position for the 2sigma error map on pol
	pos_e=[0.1,.45,0.78,0.65]
	; page position of total intensity image
	pos_i=[0.1,0.65,0.78,0.85]
	; position of the color bar for polarization
	pos_b=[0.8,0.25,0.85,0.85]

	; x axis for the the above plots
;	x = dindgen(ss[1])*sxpar(h,'CD1_1')*bin + 4930.+bin/2.	; for square binning
	wlen = dindgen(ssi[1])*sxpar(h,'CD1_1') + sxpar(h,'CRVAL1')- $
			sxpar(h,'CDPIX1')*sxpar(h,'CD1_1')				
	x = wlen[where(wlen gt 4930. and wlen lt 5032.)]			; for VT binning
	; y axis for the above plots 
	y = dindgen(ss[2])*0.25
	
	;sc=dindgen(10)/10.
	;sc=[0.,.05,.1,.15,.2,.25,.3]
	
	; crazy shit that claudia came up with to make the Frac Pol bar at the bottom
	;max_plevel=.12
	;sc=[0.,.02,.04,.06,.08,.1]
	max_plevel=.5
	sc=[0.,.05,.10,.15,.2,.25, .3,.35,.40,.45]
	sci=dblarr(n_elements(sc),3)
	sci[*,0]=sc
	sci[*,1]=sc
	sci[*,2]=sc
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; START THE 2D FIGURE
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	set_plot,'ps'
	device, filename='fig2d_bb8.ps',/color,/encap
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; DRAW 2D IMAGE SPECTRA
	loadct,0
	cgloadct,3,clip=100,/reverse;, ncolors=100
	pimg=bytscl(pol_sig,min=0.,max=max_plevel)
	tvimage,bytscl(pol_sig,min=0.,max=max_plevel),POSITION=pos_p;,/keep_aspect_ratio
	tvimage,bytscl(err_sig, min=0.,max=max_plevel), POSITION=pos_e
	loadct,0
	tvimage,bytscl(-int,min=-0.005,max=0.0009),POSITION=pos_i;,/noerase -- for 7814
;	totsm = gauss_smooth(int,.85,/edge_truncate)
;	tvimage,bytscl(-totsm,min=-0.0023,max=0.0001),POSITION=pos_i;,/keep_aspect_ratio
;	tvimage,bytscl(-int,min=-.4,max=-.01),POSITION=pos_i;,/keep_aspect_ratio;,/noerase  -- for 6514
;	tvimage,bytscl(-int,min=-0.003,max=0.0008),POSITION=pos_i;,/noerase -- for 714
	loadct,0
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; DRAW AXES FOR POLARIZATION PLOT
	plot, x, y, POSITION=pos_p, /noerase, $
		  /nodata, thick=4, charthick=4, xthick=4, ythick=4, ystyle=1, $
		  xtitle='Wavelength [Angstroms]', charsize=1.5, xstyle=1, yr=[0.001,17.99]
	oplot, dindgen(100)*0. + 4981.5,dindgen(100),thick=5,linestyle=2,color=1
	
	;xyouts,5005,2,'Polarization',charsize=1.2,charthick=4,/data
	;xyouts,5005,2,'S/N>2.0',charsize=1.2,charthick=4,/data 

	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; DRAW AXES FOR ERROR MAP PLOT
	plot, x, y, POSITION=pos_e, /noerase, /nodata, thick=4, charthick=4, $
			xthick=4, ythick=4, ystyle=1, charsize=1.5, xstyle=1, $
			yr=[0.001,17.99], xtickformat='(a2)'
	letter = "162B
	lettersigma = '!4'+string(letter) + '!X'

	;xyouts, 5005,2,'2'+lettersigma+' Error', charsize=1.3, charthick=4
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; DRAW AXES FOR INTENSITY PLOT
	plot, x, y, POSITION=pos_i, /noerase, title = title,  $
		  /nodata, thick=4, charthick=4, xthick=4, ythick=4, ystyle=1, $
		  charsize=1.5, xstyle=1, xtickformat='(a2)', yr=[0.001,17.99]
	
	;xyouts,5005,2,'Intensity',charsize=1.2,charthick=4,/data
	xyouts,4920,-18,'Arcseconds',charsize=1.5,charthick=4,/data,orientation=90
	
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; DRAW FRACTIONAL POLARIZATION BAR AT BOTTOM
	cgloadct,3,clip=100,/reverse 
	bar = replicate(1B, 10) # sc
	bar = bytscl(bar, min=0., max=max_plevel)
	tvimage, bar, position=pos_b
	loadct,0
	plots,[pos_b[0], pos_b[0], pos_b[2], pos_b[2], pos_b[0]], $
		  [pos_b[1], pos_b[3], pos_b[3], pos_b[1], pos_b[1]], thick=3, /NORMAL
	plot, [0,1], [0,1], /nodata, /noerase, position=pos_b, $ 
		  yr=[0.,sc[n_elements(sc)-1]], $
		  xstyle=4, ytickformat='(a2)', /data
	axis, yaxis=1, charthick=3, charsize=1.0, yr=[0.,max_plevel], $
		  ytickinterval=0.1,/save
	xyouts, 2.5, .35, 'P Fraction', charsize=1.5, charthick=3, $
		  orientation=270, /data
	;plot,transpose(sc),transpose(sc/sc),POSITION=pos_b,/noerase,/nodata,thick=4, $
	;	 xtitle='P fraction',charsize=1.2,xstyle=1,ystyle=4,charthick=4, $
	;	 ytickformat='(a2)',xr=[0,max_plevel]
	
	device,/close
	set_plot,'X'
	;spawn,'ps2pdf'figure_airm',1		[0.8,0.3,0.9,0.9]
stop
end


pro thetafig
	readcol,'spectropol.dat',p,perr,theta,terr

	ts=[7,13,19,26,34,44,54,61,69]
	;ts=[5,11,18,25,33,40,48,55,62,69]
	bs=[0,ts[0:n_elements(ts)-1]]

	x1=dindgen(71)*.25
	y2=dindgen(15)/100.

	j = ts[0]
	x = [0]
	; create an offsets array for plotting the pol in the middle of bins
	offsets = [0]
	; plot the red lines to delineate the boxes i summed over
	for i=0,n_elements(ts)-1 do begin
		x=[x,x1[j+1]]
		; if it's not the last red line, plot it
		offsets=[offsets,float(x[i+1]-x[i])/2.]
		if i le n_elements(ts)-2 then j+=bs[i+2]-bs[i+1] else j=ts[i]
	endfor
	; the first offset isn't really 0 so get rid of that
	offsets = offsets[1:n_elements(offsets)-1]

	pos_i=[0.2,0.51,0.9,0.9]
	pos_p=[0.2,0.1,0.9,0.5]

	set_plot,'PS'
	device,filename='pol_theta.ps',/color,/encap

	plot,x,y2,POSITION=pos_i,/noerase,/nodata,thick=4,xthick=4,ythick=4, $
		 charthick=4,charsize=1.25,xstyle=1,ystyle=1,yr=[-0.1,.4], $
		 xr=[0,17.5], yticklen=0.01, yticks=5, xticklen=0.025, $
		 ytit='Frac. P', xtickformat='(a2)' ;xtit='Arcseconds',
	oplot,x,dindgen(100)*0, linestyle=2,thick=3
	; plot the fractional polarization points using the x axis created during 
	; the plotting of the red lines
	;oplot,x+offsets,pol,psym=6,thick=3, symsize=1, color=cgcolor('blue')
	oploterror,x+offsets,p,perr,linestyle=6,psym=6,symsize=1,errthick=3, $
			   thick=3, color=cgcolor('blue')

	plot,x,y2,POSITION=pos_p,/noerase,/nodata,thick=4,xthick=4,ythick=4, $
		 charthick=4,charsize=1.25,xstyle=1,ystyle=1,yr=[-90.,270.], $
		 xr=[0,17.5], yticklen=0.01, yticks=5, xticklen=0.025, $
		 xtit='Arcseconds',ytit='Position Angle'
	oplot,x,dindgen(100)*0, linestyle=2,thick=3

	oploterror, x+offsets,theta*180./!pi,terr*180./!pi, $
			    linestyle=6,psym=6,symsize=1,errthick=3, $
			    thick=3, color=cgcolor('blue')
	
	device,/close
	set_plot,'X'

stop
	
end
	
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro fig1d, date, bin
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; 1D TOT INT top and bottom file names
	int_t_file = 'stacked/'+date+'/crop/Iall_stack_'+date+'_1d_top_crop_2.fits'
	int_b_file = 'stacked/'+date+'/crop/Iall_stack_'+date+'_1d_bot_crop_2.fits'
	; Full 1D TOT INT spectrum for wavelength purposes
	int1d_file = 'stacked/'+date+'/Iall_stack_'+date+'_1d_2.fits'
	; 1D POL top/bot file names
	pol_t_file = 'pol_products/'+date+'/pol_'+date+'_1d_top.fits'
	pol_b_file = 'pol_products/'+date+'/pol_'+date+'_1d_bot.fits'

	stn_t_file = 'pol_products/'+date+'/MC_polstn68_'+date+'_1d_top.fits'
	stn_b_file = 'pol_products/'+date+'/MC_polstn68_'+date+'_1d_top.fits'

	; READ IN DATA
	int_t = mrdfits(int_t_file)
	int_b = mrdfits(int_b_file)
	int1d = mrdfits(int1d_file,0,h)

;	pol_t = mrdfits(pol_t_file)
;	pol_b = mrdfits(pol_b_file)

;	stn_t = mrdfits(stn_t_file)
;	stn_b = mrdfits(stn_b_file)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;	SET PLOT PARAMETERS / NECESSITIES
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; CONSTRUCT WAVELENGTH RANGE FROM INT1D HEADER
	ss_int1d = size(int1d)
	wlen = dindgen(ss_int1d[1])*sxpar(h,'CD1_1')+sxpar(h,'CRVAL1')- $
		   sxpar(h,'CRPIX1')*sxpar(h,'CD1_1')
	xrange = where(wlen gt 4930. and wlen lt 5032.)
	; EXTRACT FILE NAME INFORMATION FOR LATER?	
;	name = strsplit(pol_t,'/.',/extract)

	; NOT SURE WHAT THIS IS YET...
;	ss_pol = size(pol_t)
;	xp = dindgen(ss_pol[1])*sxpar(h,'CD1_1')*bin + wlen(xrange[0]-bin/2)
	
	; SET POSITION ARGUMENTS
	pos_t = [0.15,0.55,0.85,0.9]
	pos_b = [0.15,0.2,0.85,0.545]	

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;	PREPARE DATA
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	; REMOVE POSSIBLE NANS FROM POLARIZATION DATA
;	nan = finite(pol_t,/nan)
;	pol_t[where(nan eq 1)] = 0.
;	nan = finite(pol_b,/nan)
;	pol_b[where(nan eq 1)] = 0.
		
	; SET SIGNAL-TO-NOISE CUT	
	pstncut = 3.0

	; APPLY S/N CUT TO POLARIZATION DATA
;	pol_t2 = pol_t*0. - 100
;	pol_t2[where(stn_t gt pstncut)] = pol_t[where(stn_t gt pstncut)]
	
;	pol_b2 = pol_b*0.-100
;	pol_b2[where(stn_b gt pstncut)] = pol_b[where(stn_b gt pstncut)]
	
	;pol2[0:ss2[1]/3-1]=-100.
	;pol2[ss2[1]*2/3+1:ss2[1]-1]=0.
	;stop

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;	BEGIN PLOTTING
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	set_plot,'PS'
	device,filename='fig1d_'+date+'_'+strtrim(bin,2)+'_t&b.ps',/color,/encap

	; PLOT THE TOP 1D TOT INT SPECTRA - SUPPRESS X-AXIS
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	plot, wlen[xrange], int_t, POSITION=pos_t, /noerase, thick=4, charthick=4, $
		    xthick=4, ythick=4, charsize=1.5, xstyle=1, ystyle=4, $
			xr=[4930,5032.], xtickformat='(a2)';, yr=[-0.,0.28]

	; draw the left-hand y-axis -- intensity axis
	axis, 4929.5, yaxis=0, yrange=[0,1], ystyle=1, charthick=4, charsize=1.5, $
		  ythick=4, /save;, yticks=4
	; overplot the line center given z = 3.09
	plots, [4981.5,4981.5], [0.,1.], thick=5, linestyle=2, /data
	
	; draw the right-hand y-axis -- polarization fraction axis
	axis, 5032.5,yaxis=1, yrange=[0,0.4], ystyle=1, charthick=4, charsize=1.5, $
		 ythick=4, color=cgcolor('dark green'), /save

	; CREATE LABEL FOR TOP SPECTRUM
	xyouts, 4935., 0.30, 'APERTURE A', charsize=1.2, charthick=4, /data

	; PLOT THE BOTTOM 1D TOT INT SPECTRA
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	plot, wlen[xrange], int_b, POSITION=pos_b, /noerase, thick=4, charthick=4, $
		    xthick=4, ythick=4, xtitle='Wavelength [Angstroms]', charsize=1.5, $
			xstyle=1, ystyle=4;, xr=[4930,5032.];, yr=[-0.,0.28]	
	
	; draw the left-hand y-axis -- intensity axis
	axis, 4930., yaxis=0, yrange=[0,1], ystyle=1, charthick=4, charsize=1.5, $
		  ythick=4, /save;, yticks=4
	; overplot the line center given z = 3.09
	plots, [4981.5,4981.5], [0.,1.], thick=5, linestyle=2, /data

	; draw the right-hand y-axis -- polarization fraction axis
	axis, 5032., yaxis=1, yrange=[0,0.4], ystyle=1, charthick=4, charsize=1.5, $
		 ythick=4, color=cgcolor('dark green'), /save

	; CREATE LABEL FOR BOTTOM SPECTRUM
	xyouts, 4935., 0.30, 'APERTURE B', charsize=1.2, charthick=4, /data

	; CREATE Y-AXIS LABELS
	xyouts, 4917., 0.09, 'Intensity [Arbitrary Units]', charsize=1.5, charthick=4, $
			orientation=90, /data
	xyouts, 5047., 0.12, 'Fractional Polarization', charsize=1.5, charthick=4, $
			orientation=90, /data, color=cgcolor('dark green')

	; overplot the bin lines so we can see what we've summed over
	num_x = n_elements(xrange)/bin
	for i =0,num_x-1 do begin
	;	oplot,dindgen(100)*0. + wlen(xrange(i*bin[2])),dindgen(100),linestyle=5
	;stop
	endfor
	
	; overplot the polarization fraction 
;	circsym
	
;	perr_t=pol_t2/stn_t
;	perr_e=pol_e2/stn_e
	
;	perr_t[where(perr_t lt 0.)]=0.
;	perr_t[where(perr_t gt 1.)]=0.
	
;	perr_e[where(perr_e lt 0.)]=0.
;	perr_e[where(perr_e gt 1.)]=0.
	
	;oplot,xp,pol2,linestyle=6,psym=8,symsize=4,color=cgcolor('dark green')
;	oploterror,xp,pol_t2,perr_t,linestyle=6,psym=8,symsize=4,errthick=3,color=cgcolor('blue')
;	oploterror,xp,pol_e2,perr_e,linestyle=6,psym=8,symsize=4,errthick=3,color=cgcolor('dark 	green')
;	xyouts,5020,.35,'P calc method:',charsize=1.2,charthick=4,/data;,color=cgcolor('blue')
;	xyouts,5025,.32,'Tinbergen',charsize=1.2,charthick=4,/data,color=cgcolor('blue')
;	xyouts,5025,0.30,'ESO',charsize=1.2,charthick=4,/data,color=cgcolor('dark green')
	
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

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro tot_int, date

	tot=mrdfits('stacked/Iall_avg_crop.fits')
	ss=size(tot)

	tot1d=mrdfits('stacked/715/Iall_stack_715_2d.fits',0,h)
	ss1=size(tot1d)
	wlen=dindgen(ss1[1])*sxpar(h,'CD1_1') + sxpar(h,'CRVAL1')- $
		 sxpar(h,'CRPIX1')*sxpar(h,'CD1_1')
	y=dindgen(ss[2])*0.25

	pos_i=[0.15,0.55,0.9,0.9]

	set_plot,'PS'
	device,filename='master_totint_8.1.14.ps',/color,/encap
	loadct,0
	;cgloadct,8,clip=100,/reverse;, ncolors=100
	; sigma = 0.85 corresponds to FWHM = 0.5"
	totsm = gauss_smooth(tot,.85,/edge_truncate)
	tvimage,bytscl(-totsm,min=-0.0066,max=0.00035),POSITION=pos_i,/keep_aspect_ratio


	plot,wlen[where(wlen gt 4930 and wlen lt 5032.)],y,POSITION=pos_i,/noerase, $
		 /nodata,thick=4,charthick=4,xthick=4,ythick=4, $
		 xtitle='Wavelength [Angstroms]',charsize=1.5,xstyle=1,yr=[0.001,17.5],ystyle=1
;	plots,[4930.5,4930.5],[0,20],/DATA,thick=5;,color=240
;	plots,[5031.5,5031.5],[0,20],/DATA,thick=5;,color=240

;	plots,[5034.,5034.],[2.75,7.],/DATA,thick=10;,color=240
;	xyouts,5036.,4.5,'B',/DATA,charthick=4,charsize=1.5
;	plots,[5034.,5034.],[8.75,13.5],/DATA,thick=10;,color=240
;	xyouts,5036.,10.5,'A',/DATA,charthick=4,charsize=1.5

	oplot,dindgen(100)*0. + 4981.5,dindgen(100),thick=5,linestyle=2,color=cgcolor('red')
	xyouts,4919,3,'Arcseconds',charsize=1.5,charthick=4,/data,orientation=90

;	xyouts,5010,2,'Intensity',charsize=1.2,charthick=4,/data


	device,/close
	set_plot,'X'
stop

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
