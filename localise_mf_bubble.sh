#!/bin/sh

./bin/s2fil_axiloc \
    -inp data_in/bubble_csky.fits \
    -filter_data filter_data.txt \
    -nside 128 \
    -lmax 256 \
    -out data_axiloc/out_axiloc \
    -theta_filter_adj 10

../comb/bin/comb_objgen \
    -inp data_axiloc/out_axiloc_sources.txt \
    -out data_axiloc/out_axiloc_sources.fits \
    -tmpl bubble \
    -nside 128 -include_size .true.

map2gif_multi data_axiloc/*.fits
map2gif_sk_multi data_axiloc/out_axiloc_sources.fits
