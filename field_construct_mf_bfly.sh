#!/bin/sh

# Construct field.
bin/s2fil_field_construct \
    -filter data_filter/mf_Tbfly_Bcmb_Stmpl.fil \
    -data data_in/bfly_csky.fits \
    -filetype map \
    -field data_field/field_mf_Tbfly_Bcmb_Stmpl.fld \
    -ngamma 5 \
    -write_filter .false.

# Make sky map of wavelet coefficients.
../fastcswt/bin/cswt_tr2sky \
    -inp data_field/field_mf_Tbfly_Bcmb_Stmpl_tr01.cswt \
    -nside 128 \
    -interp .true. \
    -all .true.

# Make gif of wavelet coefficient sky map.
map2gif \
    -inp data_field/field_mf_Tbfly_Bcmb_Stmpl_tr01_sky_id01_ig01.fits \
    -out data_field/field_mf_Tbfly_Bcmb_Stmpl_tr01_sky_id01_ig01.gif  \
    -bar .true.

