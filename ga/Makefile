CC = gcc
#CFLAGS = -ansi -fullwarn -O
LIBS = -lm -lc
OFILES = main.o ga.o params.o output.o pop.o init.o random.o stats.o reproduce.o util.o genops.o indv.o build.o parse.o futil.o
FXOFILES = fxtsp.o
HFILES = types.h global.h extern.h

ga:	$(OFILES) $(FXOFILES) $(HFILES)
	$(CC) $(CFLAGS) $(OFILES) $(FXOFILES) -o ga $(LIBS)

main.o:	main.c $(HFILES)
	$(CC) $(CFLAGS) -c main.c 

ga.o:	ga.c $(HFILES)
	$(CC) $(CFLAGS) -c ga.c 

params.o:	params.c params.h $(HFILES)
	$(CC) $(CFLAGS) -c params.c 

output.o:	output.c output.h $(HFILES)
	$(CC) $(CFLAGS) -c output.c 

pop.o:	pop.c pop.h $(HFILES)
	$(CC) $(CFLAGS) -c pop.c 

init.o:	init.c init.h $(HFILES)
	$(CC) $(CFLAGS) -c init.c 

random.o:	random.c random.h $(HFILES)
	$(CC) $(CFLAGS) -c random.c 

stats.o:	stats.c $(HFILES)
	$(CC) $(CFLAGS) -c stats.c 

reproduce.o:	reproduce.c $(HFILES)
	$(CC) $(CFLAGS) -c reproduce.c 

util.o:	util.c $(HFILES)
	$(CC) $(CFLAGS) -c util.c

genops.o:	genops.c $(HFILES)
	$(CC) $(CFLAGS) -c genops.c 

indv.o:	indv.c $(HFILES)
	$(CC) $(CFLAGS) -c indv.c 

build.o:	build.c $(HFILES)
	$(CC) $(CFLAGS) -c build.c

parse.o:	parse.c $(HFILES)
	$(CC) $(CFLAGS) -c parse.c

futil.o:	futil.c $(HFILES)
	$(CC) $(CFLAGS) -c futil.c

fxtsp.o:    fxtsp.c fxtsp.h $(HFILES)
	$(CC) $(CFLAGS) -c fxtsp.c 

fxtc.o:	fxtc.c fxtc.h $(HFILES)
	$(CC) $(CFLAGS) -c fxtc.c 
    