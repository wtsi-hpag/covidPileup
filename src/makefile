# Makefile for covidPileup 
CC= gcc
CFLAGS= -O4
LFLAGS= -lm -pthread -lz

SOURCES=covid_fastq.c covid_frequ.c covid_comms.c covidSNP.c covidPileup.c

OBJS = $(patsubst %.c,%.o,$(SOURCES)) fast.o
EXECS = $(patsubst %.c,%,$(SOURCES))
EXECS_BIN = $(patsubst %.c,covid-bin/%,$(SOURCES))
COMPILE = $(CC) $(CFLAGS) 


all:  cleanall iprint $(OBJS) executables clean oprint

executables:
	for exe in $(EXECS);  do $(COMPILE) -o $$exe $$exe.o fast.o $(LFLAGS); cp $$exe covid-bin/.; done

%.o: %.c fasta.h
	$(CC) $(CFLAGS)  -c $< -o $@

iprint:
	@echo '+++ Compiling All ... '

oprint:
	@echo 'All Done '


clean: 
	@echo '+++ Cleaning Up ... '
	@rm -f $(EXECS)
	@rm -f $(OBJS)
	@cp covid-bin/covidPileup .

cleanall: 
	@echo '+++ Cleaning All ... '
	@rm -f $(EXECS)
	@rm -f $(OBJS) fast.o
	@rm -f $(EXECS_BIN)
