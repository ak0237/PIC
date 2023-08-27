#include <stdio.h>
#include "../rps.h"

void op(int l, double *p_1){
        int x, y;
        FILE *file;
        char name[100];

        sprintf(name, "dat/p_i-%d.dat", l);
        file= fopen(name, "w");
        for(x= 0; x< Nx; x++){
                for(y= 0; y< Ny; y++){
                        fprintf(file, "%.3e ", p_1[x*Ny+y]);
                }
                fprintf(file, "\n");
        }
        fclose(file);
}

