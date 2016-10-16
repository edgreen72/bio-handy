set terminal postscript eps color "Arial" 14;
#############################################################
#                CUSTOMIZE THESE OPTIONS                    #
#############################################################
# This template makes a map-damage like plot from the
# output file of pss-bam.pl run from merged reads, i.e.,
# using the -m option to pss-bam.pl
# Customize these 4 lines to control what the input file and
# output files are and how much damage was present in the
# reads. Then run:
# gnuplot THIS_FILE
# to generate the plot
#############################################################
set ytic 0.05;                        # tic marks on y-axis #
set output "JK469.pss.v0.6.eps";          # output filename #
datafile = "JK469.pss.dat";         # input pss-bam.pl file #
max_damage = 0.1           # high fraction of C->T in plot  #
#############################################################

set multiplot;
set size 0.5, 1;
set grid;
set style line 1 lt 1 lc 7 lw 0.5;
set style line 2 lt 1 lc rgb "#8b0000" lw 3;
set style line 3 lt 1 lc rgb "#228b22" lw 3;
set style line 4 lc rgb "green"; # A = green in Sanger
set style line 5 lc rgb "blue";  # C
set style line 6 lc rgb "black"; # G
set style line 7 lc rgb "red"; # T
set nokey;
unset xtic;
set tmargin 0.9;
set bmargin 0.5;
set style data histogram;
set style histogram rowstacked;
set style fill transparent solid 0.4 noborder;
set boxwidth 0.8;
set style rect fc lt -1 fs solid 0.35 noborder;
set obj rect from -1.0, graph 0 to 1.5, graph 1;
r = 14;
set origin 0, 0;

#set title "JK469 - Merged";

plot [-1.0:17.0][0:max_damage] datafile index 0 using ($2+$3+$4+$5) axes x1y2 ls 4,\
     '' index 0 using ($10+$11+$12+$13) axes x1y2 ls 6,\
     '' index 0 using ($6+$7+$8+$9) axes x1y2 ls 5,\
     '' index 0 using ($14+$15+$16+$17) axes x1y2 ls 7,\
     '' index 0 u ($1+2):($3==0?NaN:$3/($2+$3+$4+$5)) t "A->C" w l ls 1,\
     '' index 0 u ($1+2):($4==0?NaN:$4/($2+$3+$4+$5)) t "A->G" w l ls 1,\
     '' index 0 u ($1+2):($5==0?NaN:$5/($2+$3+$4+$5)) t "A->T" w l ls 1,\
     '' index 0 u ($1+2):($6==0?NaN:$6/($6+$7+$8+$9)) t "C->A" w l ls 1,\
     '' index 0 u ($1+2):($8==0?NaN:$8/($6+$7+$8+$9)) t "C->G" w l ls 1,\
     '' index 0 u ($1+2):($9==0?NaN:$9/($6+$7+$8+$9)) t "C->T" w l ls 2,\
     '' index 0 u ($1+2):($10==0?NaN:$10/($10+$11+$12+$13)) t "G->A" w l ls 3,\
     '' index 0 u ($1+2):($11==0?NaN:$11/($10+$11+$12+$13)) t "G->C" w l ls 1,\
     '' index 0 u ($1+2):($13==0?NaN:$13/($10+$11+$12+$13)) t "G->T" w l ls 1,\
     '' index 0 u ($1+2):($14==0?NaN:$14/($14+$15+$16+$17)) t "T->A" w l ls 1,\
     '' index 0 u ($1+2):($15==0?NaN:$15/($14+$15+$16+$17)) t "T->C" w l ls 1,\
     '' index 0 u ($1+2):($16==0?NaN:$16/($14+$15+$16+$17)) t "T->G" w l ls 1;

set origin 0.5, 0;
unset obj;
set obj rect from 14.5, graph 0 to 17, graph 1;
plot [-1.0:17.0][0:max_damage] datafile index 1 using ($2+$3+$4+$5) axes x1y2 ls 4,\
     '' index 1 using ($10+$11+$12+$13) axes x1y2 ls 6,\
     '' index 1 using ($6+$7+$8+$9) axes x1y2 ls 5,\
     '' index 1 using ($14+$15+$16+$17) axes x1y2 ls 7,\
     '' index 1 using ($3==0?NaN:$3/($2+$3+$4+$5)) t "A->C" w l ls 1,\
     '' index 1 using ($4==0?NaN:$4/($2+$3+$4+$5)) t "A->G" w l ls 1,\
     '' index 1 using ($5==0?NaN:$5/($2+$3+$4+$5)) t "A->T" w l ls 1,\
     '' index 1 using ($6==0?NaN:$6/($6+$7+$8+$9)) t "C->A" w l ls 1,\
     '' index 1 using ($8==0?NaN:$8/($6+$7+$8+$9)) t "C->G" w l ls 1,\
     '' index 1 using ($9==0?NaN:$9/($6+$7+$8+$9)) t "C->T" w l ls 2,\
     '' index 1 using ($10==0?NaN:$10/($10+$11+$12+$13)) t "G->A" w l ls 3,\
     '' index 1 using ($11==0?NaN:$11/($10+$11+$12+$13)) t "G->C" w l ls 1,\
     '' index 1 using ($13==0?NaN:$13/($10+$11+$12+$13)) t "G->T" w l ls 1,\
     '' index 1 using ($14==0?NaN:$14/($14+$15+$16+$17)) t "T->A" w l ls 1,\
     '' index 1 using ($15==0?NaN:$15/($14+$15+$16+$17)) t "T->C" w l ls 1,\
     '' index 1 using ($16==0?NaN:$16/($14+$15+$16+$17)) t "T->G" w l ls 1;
