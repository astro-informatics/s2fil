set terminal postscript color
set output 'mf_bubble_spectra.eps'

set xlabel 'Multipole'
show xlabel
set ylabel 'Cl'
show ylabel
set title 'BUBBLE'
show title

set logscale xy
set yrange[1e-14:1e4]

set output 'mf_bubble_spectra.eps'
plot [2:1000] 'data_filter/mf_Tbubble_Bcmb_Stmpl_bkgnd.cl.txt' \
	with linespoints pointtype 7 pointsize 0.5 lw 2 title "Background"

set output 'mf_bubble_spectra.eps'
replot 'data_filter/mf_Tbubble_Bcmb_Stmpl_sky01.cl.txt' \
	with linespoints pointtype 7 pointsize 0.5 lw 2 title "Filter"

set output 'mf_bubble_spectra.eps'
replot 'data_filter/mf_Tbubble_Bcmb_Stmpl_tmpl_d01.cl.txt' \
	with linespoints pointtype 7 pointsize 0.5 lw 2 title "Template"


