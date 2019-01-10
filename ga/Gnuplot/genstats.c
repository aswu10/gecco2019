#include <stdio.h>
#include <stdlib.h>

int main(int argc, char **argv)
   {
   int run_num;
   int display;
   char gnuplotfile[200];
   char datafile[500];
   char file_with_plot[200];
   char systemcmd[100];
   FILE *fp;

   if (argc < 3)
      {
      printf(" Usage: %s <run number> <display plot flag>\n", argv[0]);
      printf("    Creates an eps plot of best & average fitness");
      printf(" (and stdev) for specified run.\n");
      printf("    Assumes run data is in ../Output/run.<run_num>/ directory.\n");
      printf("    <display plot flag> 1 = open plot, 0 = don't open plot\n");
      return 0;
      }

   sscanf(argv[1], "%d", &run_num);
   printf(" Run number: %d\n", run_num);
   sscanf(argv[2], "%d", &display);
   printf(" Display plot flag: %d\n", display);

   // open file to print to
   sprintf(gnuplotfile, "genstats.gnu");
   fp = fopen(gnuplotfile, "w");

   // save path to input data file (file of data to plot)
   sprintf(datafile, "../Output/run.%d/run.%d.genstats",
           run_num, run_num);
   printf(" Datafile: %s\n", datafile);

   // name of output file with plot
   sprintf(file_with_plot, "run.%d.genstats.eps", run_num);

   // create gnuplot file
   fprintf(fp, "set term post eps\n");
   fprintf(fp, "set output \"%s\"\n", file_with_plot);
   fprintf(fp, "set xlabel \"Generation\"\n");
   fprintf(fp, "set ylabel \"Fitness\"\n");
   fprintf(fp, "set title \"Run %d\"\n\n", run_num);
   fprintf(fp, "plot \\\n");
   fprintf(fp, "   \"%s\" using 1:7 title \"Best fitness\" w line,\\\n",
           datafile);
   fprintf(fp, "   \"%s\" using 1:3 title \"Average fitness\" w line,\\\n",
           datafile);
   fprintf(fp, "   \"%s\" using 1:5 title \"Standard deviation\" w line\n",
           datafile);

   fclose(fp);

   system("gnuplot genstats.gnu");
   if (display)
      {
      sprintf(systemcmd, "open %s", file_with_plot);
      system(systemcmd);
      }

   return 0;
   }
