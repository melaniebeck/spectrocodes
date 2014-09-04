# this is 8 groups of 10 images each to create the ord and ext at 4 angles
!ls N*/2dspectra/N*slit1*00*wc_ssfit2_acorr3.fits > combine_temp.list
imcombine @combine_temp.list stacked/A00_stack_715_2d_sc.fits combine=average offsets=wcs reject=sigclip nlow=1 nhigh=1 sigmas=stacked/A00_stack_715_2d_sc_sig.fits

!ls N*/2dspectra/N*slit1*22*wc_ssfit2_acorr3.fits > combine_temp.list
imcombine @combine_temp.list stacked/A22_stack_715_2d_sc.fits combine=average offsets=wcs reject=sigclip nlow=1 nhigh=1 sigmas=stacked/A22_stack_715_2d_sc_sig.fits

!ls N*/2dspectra/N*slit1*45*wc_ssfit2_acorr3.fits > combine_temp.list
imcombine @combine_temp.list stacked/A45_stack_715_2d_sc.fits combine=average offsets=wcs reject=sigclip nlow=1 nhigh=1 sigmas=stacked/A45_stack_715_2d_sc_sig.fits

!ls N*/2dspectra/N*slit1*67*wc_ssfit2_acorr3.fits > combine_temp.list
imcombine @combine_temp.list stacked/A67_stack_715_2d_sc.fits combine=average offsets=wcs reject=sigclip nlow=1 nhigh=1 sigmas=stacked/A67_stack_715_2d_sc_sig.fits

!ls N*/2dspectra/N*slit2*00*wc_ssfit2_acorr3.fits > combine_temp.list
imcombine @combine_temp.list stacked/B00_stack_715_2d_sc.fits combine=average offsets=wcs reject=sigclip nlow=1 nhigh=1 sigmas=stacked/B00_stack_715_2d_sc_sig.fits

!ls N*/2dspectra/N*slit2*22*wc_ssfit2_acorr3.fits > combine_temp.list
imcombine @combine_temp.list stacked/B22_stack_715_2d_sc.fits combine=average offsets=wcs reject=sigclip nlow=1 nhigh=1 sigmas=stacked/B22_stack_715_2d_sc_sig.fits

!ls N*/2dspectra/N*slit2*45*wc_ssfit2_acorr3.fits > combine_temp.list
imcombine @combine_temp.list stacked/B45_stack_715_2d_sc.fits combine=average offsets=wcs reject=sigclip nlow=1 nhigh=1 sigmas=stacked/B45_stack_715_2d_sc_sig.fits

!ls N*/2dspectra/N*slit2*67*wc_ssfit2_acorr3.fits > combine_temp.list
imcombine @combine_temp.list stacked/B67_stack_715_2d_sc.fits combine=average offsets=wcs reject=sigclip nlow=1 nhigh=1 sigmas=stacked/B67_stack_715_2d_sc_sig.fits


# this is 4 groups of 20 images to create individual intensity images at each angle
!ls N*/2dspectra/N*slit[1,2]*00*wc_ssfit2_acorr3.fits > combine_temp.list
imcombine @combine_temp.list stacked/I00_stack_715_2d_sc.fits combine=average offsets=wcs reject=sigclip nlow=1 nhigh=1 sigmas=stacked/I00_stack_715_2d_sc_sig.fits

!ls N*/2dspectra/N*slit[1,2]*22*wc_ssfit2_acorr3.fits > combine_temp.list
imcombine @combine_temp.list stacked/I22_stack_715_2d_sc.fits combine=average offsets=wcs reject=sigclip nlow=1 nhigh=1 sigmas=stacked/I22_stack_715_2d_sc_sig.fits

!ls N*/2dspectra/N*slit[1,2]*45*wc_ssfit2_acorr3.fits > combine_temp.list
imcombine @combine_temp.list stacked/I45_stack_715_2d_sc.fits combine=average offsets=wcs reject=sigclip nlow=1 nhigh=1 sigmas=stacked/I45_stack_715_2d_sc_sig.fits

!ls N*/2dspectra/N*slit[1,2]*67*wc_ssfit2_acorr3.fits > combine_temp.list
imcombine @combine_temp.list stacked/I67_stack_715_2d_sc.fits combine=average offsets=wcs reject=sigclip nlow=1 nhigh=1 sigmas=stacked/I67_stack_715_2d_sc_sig.fits


# this is 1 group of 80 images 
!ls N*/2dspectra/N*slit[1,2]*wc_ssfit2_acorr3.fits > combine_temp.list
imcombine @combine_temp.list stacked/Iall_stack_715_2d_sc.fits combine=average offsets=wcs reject=sigclip nlow=1 nhigh=1 sigmas=stacked/Iall_stack_715_2d_sc_sig.fits


# this is 1 group of 4 Total Intensity images that were previously summed 
#!ls stacked/715/I[00,22,45,67]*_stack_715_2d_sc_sum2.fits > combine_temp.list
#imcombine @combine_temp.list stacked/Iall_stack_715_2d_sc_avgofsum.fits combine=average #offsets=wcs reject=sigclip #nlow=1 nhigh=1 sigmas=stacked/Iall_stack_715_2d_sc_sig_avgofsum.fits
