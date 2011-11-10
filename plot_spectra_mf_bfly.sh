#!/bin/sh

# Plot Butterfly optimal filter and associated data spectra
#################################################################

# Produce cl files of filter sky and template (dilating template)

s2_sky2cl \
    -inp data_filter/mf_Tbfly_Bcmb_Stmpl_sky01.sky \
    -type sky \
    -out data_filter/mf_Tbfly_Bcmb_Stmpl_sky01.cl \
    -npres .false.

s2_sky2cl \
    -inp data_filter/mf_Tbfly_Bcmb_Stmpl_tmpl.sky \
    -type sky \
    -out data_filter/mf_Tbfly_Bcmb_Stmpl_tmpl_d01.cl \
    -beam data_filter/mf_Tbfly_Bcmb_Stmpl_beam.cl \
    -dil1 0.1 \
    -npres .false. 


# Produce ascii files of cl spectra for gnuplot input

s2_cl2ascii \
    -inp data_filter/mf_Tbfly_Bcmb_Stmpl_tmpl_d01.cl \
    -out data_filter/mf_Tbfly_Bcmb_Stmpl_tmpl_d01.cl.txt \
    -dl_status .false.

s2_cl2ascii \
    -inp data_filter/mf_Tbfly_Bcmb_Stmpl_sky01.cl \
    -out data_filter/mf_Tbfly_Bcmb_Stmpl_sky01.cl.txt \
    -dl_status .false.

s2_cl2ascii \
    -inp data_filter/mf_Tbfly_Bcmb_Stmpl_bkgnd.cl \
    -out data_filter/mf_Tbfly_Bcmb_Stmpl_bkgnd.cl.txt \
    -dl_status .false.


# Run gnuplot to produce plot

gnuplot gnuplot_mf_bfly.gp
evince mf_bfly_spectra.eps &
