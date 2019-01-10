set term post eps
set output "run.1.genstats.eps"
set xlabel "Generation"
set ylabel "Fitness"
set title "Run 1"

plot \
   "../Output/run.1/run.1.genstats" using 1:7 title "Best fitness" w line,\
   "../Output/run.1/run.1.genstats" using 1:3 title "Average fitness" w line,\
   "../Output/run.1/run.1.genstats" using 1:5 title "Standard deviation" w line
