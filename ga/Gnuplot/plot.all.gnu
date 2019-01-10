set term post eps
set xlabel "Generation"
set ylabel "Fitness"

set output "20runs.elite-on.eps"
set title "20 runs, Elitism off"
plot \
   "../Output/run.1/run.1.genstats" using 1:7 title "" w line,\
   "../Output/run.2/run.2.genstats" using 1:7 title "" w line,\
   "../Output/run.3/run.3.genstats" using 1:7 title "" w line,\
   "../Output/run.4/run.4.genstats" using 1:7 title "" w line,\
   "../Output/run.5/run.5.genstats" using 1:7 title "" w line,\
   "../Output/run.6/run.6.genstats" using 1:7 title "" w line,\
   "../Output/run.7/run.7.genstats" using 1:7 title "" w line,\
   "../Output/run.8/run.8.genstats" using 1:7 title "" w line,\
   "../Output/run.9/run.9.genstats" using 1:7 title "" w line,\
   "../Output/run.10/run.10.genstats" using 1:7 title "" w line,\
   "../Output/run.11/run.11.genstats" using 1:7 title "" w line,\
   "../Output/run.12/run.12.genstats" using 1:7 title "" w line,\
   "../Output/run.13/run.13.genstats" using 1:7 title "" w line,\
   "../Output/run.14/run.14.genstats" using 1:7 title "" w line,\
   "../Output/run.15/run.15.genstats" using 1:7 title "" w line,\
   "../Output/run.16/run.16.genstats" using 1:7 title "" w line,\
   "../Output/run.17/run.17.genstats" using 1:7 title "" w line,\
   "../Output/run.18/run.18.genstats" using 1:7 title "" w line,\
   "../Output/run.19/run.19.genstats" using 1:7 title "" w line,\
   "../Output/run.20/run.20.genstats" using 1:7 title "" w line

set nolabel
set output "20runs.elite-off.eps"
set title "20 runs, Elitism on"
plot \
   "../Output/run.21/run.21.genstats" using 1:7 title "" w line,\
   "../Output/run.22/run.22.genstats" using 1:7 title "" w line,\
   "../Output/run.23/run.23.genstats" using 1:7 title "" w line,\
   "../Output/run.24/run.24.genstats" using 1:7 title "" w line,\
   "../Output/run.25/run.25.genstats" using 1:7 title "" w line,\
   "../Output/run.26/run.26.genstats" using 1:7 title "" w line,\
   "../Output/run.27/run.27.genstats" using 1:7 title "" w line,\
   "../Output/run.28/run.28.genstats" using 1:7 title "" w line,\
   "../Output/run.29/run.29.genstats" using 1:7 title "" w line,\
   "../Output/run.30/run.30.genstats" using 1:7 title "" w line,\
   "../Output/run.31/run.31.genstats" using 1:7 title "" w line,\
   "../Output/run.32/run.32.genstats" using 1:7 title "" w line,\
   "../Output/run.33/run.33.genstats" using 1:7 title "" w line,\
   "../Output/run.34/run.34.genstats" using 1:7 title "" w line,\
   "../Output/run.35/run.35.genstats" using 1:7 title "" w line,\
   "../Output/run.36/run.36.genstats" using 1:7 title "" w line,\
   "../Output/run.37/run.37.genstats" using 1:7 title "" w line,\
   "../Output/run.38/run.38.genstats" using 1:7 title "" w line,\
   "../Output/run.39/run.39.genstats" using 1:7 title "" w line,\
   "../Output/run.40/run.40.genstats" using 1:7 title "" w line
