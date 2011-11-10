#!/bin/sh

bin/s2fil_filter_construct \
    -out data_filter/mf_Tbubble_Bcmb_Stmpl.fil \
    -nside 128 \
    -filename_dil data_in/dilation_unity.txt \
    -bkgnd_cmb data_in/wmap_lcdm_pl_model_wmap7baoh0_CAMB.dat \
    -noise_var 0.0048e0 \
    -beam_fwhm 13.2e0 \
    -tmpl bubble \
    -filter_heu .true. \
    -filter_type mf \
    -scale_type tmpl


# mean
#noise_var = (6.549/sqrt(mean(Nobs)))^2 = 0.004215894219296

# median
# (6.549/sqrt(median(Nobs)))^2 = 0.004827070021056