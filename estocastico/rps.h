#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>

#define N  15000  /* number of individuals          */
#define L  1.0    /* grid size                      */
#define R  0.1    /* radius of interections - pre.  */
#define NG 500    /* number of generations          */
#define NF 10     /* number of output files         */
#define LS 0      /* 0-> nolog time || 1-> log time */

#define alpha -1.0
#define N_S 200
#define lambda_0 1.7
#define beta_0 -0.7
#define sigma 0.00447213595499957939
#define eta 1.0

struct DATA{
	int s;
	double  x;
	double  y;
};
