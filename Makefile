#
# Makefile for REED-MPI
.SUFFIXES: .o .c
CXX = /usr/bin/mpicxx
FLAGS = -O3 -ansi
LIB = -lm 

OBJT = main.o domain.o spawning.o vverlet.o force.o util.o
TARG = REED_2007C_MPI

${TARG}: ${OBJT}
	${CXX} ${FLAGS} -o ${TARG} ${OBJT} ${LIB}

.cc.o:
	${CXX} ${FLAGS} -c $<
clean:
	rm *.o *.h~ *.cc~ ${TARG}
