#!/bin/sh

# Construct and plot Bubble optimal filter and associated data
#################################################################

# Construct filter.
./filter_construct_mf_bubble.sh

# Make fits map of sky.
s2_sky2map \
    -sky data_filter/mf_Tbubble_Bcmb_Stmpl_sky01.sky \
    -map data_filter/mf_Tbubble_Bcmb_Stmpl_sky01.fits \
    -beta_center .true.

# Map image of sky fits map.
map2gif -inp data_filter/mf_Tbubble_Bcmb_Stmpl_sky01.fits \
    -out data_filter/mf_Tbubble_Bcmb_Stmpl_sky01.gif -bar .true.

# Make fits map of tmpl.
s2_sky2map \
    -sky data_filter/mf_Tbubble_Bcmb_Stmpl_tmpl.sky \
    -map data_filter/mf_Tbubble_Bcmb_Stmpl_tmpl.fits \
    -dil1 1.0 \
    -beta_center .true.

# Map image of sky fits map.
map2gif -inp data_filter/mf_Tbubble_Bcmb_Stmpl_tmpl.fits \
    -out data_filter/mf_Tbubble_Bcmb_Stmpl_tmpl.gif -bar .true.




