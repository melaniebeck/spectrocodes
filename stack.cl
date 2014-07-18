# this is 8 groups of 10 images each to create the ord and ext at 4 angles
!ls N*/1dspectra/N*slit1*00*wc_1d_ssfit2_acorr.fits > combine_temp.list
imcombine @combine_temp.list stacked/A00_stack_712.fits combine=average offsets=wcs reject=sigclip nlow=1 nhigh=1 sigmas=stacked/A00_stack_712_sig.fits

!ls N*/1dspectra/N*slit1*22*wc_1d_ssfit2_acorr.fits > combine_temp.list
imcombine @combine_temp.list stacked/A22_stack_712.fits combine=average offsets=wcs reject=sigclip nlow=1 nhigh=1 sigmas=stacked/A22_stack_712_sig.fits

!ls N*/1dspectra/N*slit1*45*wc_1d_ssfit2_acorr.fits > combine_temp.list
imcombine @combine_temp.list stacked/A45_stack_712.fits combine=average offsets=wcs reject=sigclip nlow=1 nhigh=1 sigmas=stacked/A45_stack_712_sig.fits

!ls N*/1dspectra/N*slit1*67*wc_1d_ssfit2_acorr.fits > combine_temp.list
imcombine @combine_temp.list stacked/A67_stack_712.fits combine=average offsets=wcs reject=sigclip nlow=1 nhigh=1 sigmas=stacked/A67_stack_712_sig.fits

!ls N*/1dspectra/N*slit2*00*wc_1d_ssfit2_acorr.fits > combine_temp.list
imcombine @combine_temp.list stacked/B00_stack_712.fits combine=average offsets=wcs reject=sigclip nlow=1 nhigh=1 sigmas=stacked/B00_stack_712_sig.fits

!ls N*/1dspectra/N*slit2*22*wc_1d_ssfit2_acorr.fits > combine_temp.list
imcombine @combine_temp.list stacked/B22_stack_712.fits combine=average offsets=wcs reject=sigclip nlow=1 nhigh=1 sigmas=stacked/B22_stack_712_sig.fits

!ls N*/1dspectra/N*slit2*45*wc_1d_ssfit2_acorr.fits > combine_temp.list
imcombine @combine_temp.list stacked/B45_stack_712.fits combine=average offsets=wcs reject=sigclip nlow=1 nhigh=1 sigmas=stacked/B45_stack_712_sig.fits

!ls N*/1dspectra/N*slit2*67*wc_1d_ssfit2_acorr.fits > combine_temp.list
imcombine @combine_temp.list stacked/B67_stack_712.fits combine=average offsets=wcs reject=sigclip nlow=1 nhigh=1 sigmas=stacked/B67_stack_712_sig.fits


# this is 4 groups of 20 images to create individual intensity images at each angle
!ls N*/1dspectra/N*slit[1,2]*00*wc_1d_ssfit2_acorr.fits > combine_temp.list
imcombine @combine_temp.list stacked/I00_stack_712.fits combine=average offsets=wcs reject=sigclip nlow=1 nhigh=1 sigmas=stacked/I00_stack_712_sig.fits

!ls N*/1dspectra/N*slit[1,2]*22*wc_1d_ssfit2_acorr.fits > combine_temp.list
imcombine @combine_temp.list stacked/I22_stack_712.fits combine=average offsets=wcs reject=sigclip nlow=1 nhigh=1 sigmas=stacked/I22_stack_712_sig.fits

!ls N*/1dspectra/N*slit[1,2]*45*wc_1d_ssfit2_acorr.fits > combine_temp.list
imcombine @combine_temp.list stacked/I45_stack_712.fits combine=average offsets=wcs reject=sigclip nlow=1 nhigh=1 sigmas=stacked/I45_stack_712_sig.fits

!ls N*/1dspectra/N*slit[1,2]*67*wc_1d_ssfit2_acorr.fits > combine_temp.list
imcombine @combine_temp.list stacked/I67_stack_712.fits combine=average offsets=wcs reject=sigclip nlow=1 nhigh=1 sigmas=stacked/I67_stack_712_sig.fits


# this is 1 group of 80 images 
!ls N*/1dspectra/N*slit[1,2]*wc_1d_ssfit2_acorr.fits > combine_temp.list
imcombine @combine_temp.list stacked/Iall_stack_712.fits combine=average offsets=wcs reject=sigclip nlow=1 nhigh=1 sigmas=stacked/Iall_stack_712_sig.fits
