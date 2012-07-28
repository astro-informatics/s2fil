#!/bin/sh

./bin/s2fil_filter_construct \
    -out data_filter/mf_Ttexture_Bcmb_Stmpl.fil \
    -nside 128 \
    -filename_dil data_in/dilation_unity.txt \
    -bkgnd_cmb data_in/wmap_lcdm_pl_model_wmap7baoh0_CAMB.dat \
    -noise_var 0.02e0 \
    -beam_fwhm 13.2e0 \
    -tmpl texture \
    -tmpl_param tmplparam01_texture.txt \
    -filter_heu .true. \
    -filter_type mf \
    -scale_type tmpl

#mean(Nobs)
#ans =
#     1.017326307754467e+04
#median(Nobs)
#ans =
#     8.885183105468750e+03
#std(Nobs)
#ans =
#     3.525765578333145e+03

# noise based on median
#n_var = (6.549/sqrt(median(Nobs)))^2
#n_var =
#   0.004827070021056
#n_std = (6.549/sqrt(median(Nobs)))
#n_std =
#   0.069477118687060

# noise based on mean
#n_var = (6.549/sqrt(mean(Nobs)))^2
#n_var =
#   0.004215894219296
#n_std = (6.549/sqrt(mean(Nobs)))
#n_std =
#   0.064929917752113