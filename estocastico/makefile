COMPILER = gcc
FLAGS = -Wall -O3 -mtune=native

LIB = -lgsl -lgslcblas -lm

ndm:
	@${COMPILER} ${FLAGS} src/ndm.c  ${LIB} -o ndm.out

clean:
	@rm -f rps.out
	@rm -f dat/*.dat

clean-fig:
	@rm -f plt/*.png
	@rm -f plt/*.pdf
	@rm -f plt/*.mp4
