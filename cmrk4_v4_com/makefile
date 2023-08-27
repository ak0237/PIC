COMPILER = gcc
FLAGS = -Wall -O3 -mtune=native

LIB = -lgsl -lgslcblas -lm

rps:
	@${COMPILER} ${FLAGS} src/rps.c src/op.c ${LIB} -o rps.out

dft:
	@${COMPILER} ${FLAGS} src/dft.c ${LIB} -o dft.out

avg-dft:
	@${COMPILER} ${FLAGS} src/avg-dft.c ${LIB} -o avg-dft.out

clean:
	@rm -f dat/*.dat
	@rm -f rps.out

clean-fig:
	@rm -f plt/*.png
	@rm -f plt/*.pdf
	@rm -f plt/*.mp4
