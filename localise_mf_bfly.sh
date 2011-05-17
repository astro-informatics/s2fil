#!/bin/sh

# Localise compact sources.
bin/s2fil_localisation_thres \
    -field data_field/field_mf_Tbfly_Bcmb_Stmpl.fld \
    -filter data_filter/mf_Tbfly_Bcmb_Stmpl.fil \
    -nsigma 3.0e0 \
    -no_gamma_search .true. \
    -thres data_localise/thres_mf_bfly.cswt \
    -connected data_localise/conn_mf_bfly.cswt \
    -out data_localise/localisation_mf_bfly.txt 
#    > data_localise/loc_mf_bfly.txt

../comb-1.0/bin/comb_objgen \
    -inp data_localise/localisation_mf_bfly.txt \
    -out data_localise/bfly_obj_loc.fits \
    -tmpl butterfly -dil 0.1 -nside 128

map2gif -inp data_localise/bfly_obj_loc.fits \
    -out data_localise/bfly_obj_loc.gif -bar .true.

