#!/bin/sh

bin/s2fil_filter_construct \
    -out data_filter/mf_Tbfly_Bcmb_Stmpl.fil \
    -nside 128 \
    -filename_dil data_in/dilation_radian.txt \
    -bkgnd_cmb data_in/wmap_lcdm_pl_model_yr1_v1.txt \
    -noise_var 0.05e0 \
    -beam_fwhm 13.0e0 \
    -tmpl butterfly \
    -filter_heu .true. \
    -filter_type mf \
    -scale_type tmpl


