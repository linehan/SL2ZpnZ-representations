set datafile separator ","

#set terminal png size 1024, 768 


#plot 'g1_spectrum_disc.csv' using 4:2 pt 7 ps 0.5 lc rgb "red" #title "Principal series spectrum for the Cayley graph of SL2(F_p) with respect to generator set G1"

#set output "random_spectrum3.png"
#plot 'random_spectrum3.csv' using 4:2 pt 7 ps 0.5 lc rgb "red"

#set output "sle_princ.png"

#plot 'sle_growth_princ.csv' using 2:3 with linespoints lc rgb "red" #title "Principal series spectrum for the Cayley graph of SL2(F_p) with respect to generator set G2"

set terminal pdf size 11in,8.5in

set output "ME.pdf"
plot 'g1_spectrum_disc.csv' using 4:2 pt 7 ps 0.2 lc rgb "black" title "" #title "Principal series spectrum for the Cayley graph of SL2(F_p) with respect to generator set G1"

#set output "g2_evalues_113.pdf"
#plot 'g2_spectrum_113.csv' using 4:2 pt 7 ps 0.2 lc rgb "black" title "" #title "Principal series spectrum for the Cayley graph of SL2(F_p) with respect to generator set G1"


#set output "g2_spectrum_113.png"
#plot 'g2_spectrum_113.csv' using 4:2 pt 7 ps 0.5 lc rgb "red" title "Principal series spectrum for the Cayley graph of SL2(F_p) with respect to generator set G2"

#pause -1
