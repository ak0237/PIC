#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include "../rps.h"

#define r (tf/NF)

void op(int, double *);
double canetaazul(int, int, double *);
double canetaazul23(int, int, double *, double *);
double canetaazul4(int, int, double *, double *);

int main(int argc, char **argv){
	int i, j, l=0;
	double t= 0.0;
	double P_1 = 0.0;
	double a0= r;
	double *phi, *k_1, *k_2, *k_3, *k_4;
	phi = (double *) calloc(Nx*Ny, sizeof(double));
	k_1 = (double *) calloc(Nx*Ny, sizeof(double));
	k_2 = (double *) calloc(Nx*Ny, sizeof(double));
	k_3 = (double *) calloc(Nx*Ny, sizeof(double));
	k_4 = (double *) calloc(Nx*Ny, sizeof(double));

	// Definindo uma seed aleatoria	
	gsl_rng_default_seed= (argc == 2) ? atoi(argv[1]) : time(NULL) ;

	// Printando a seed
	printf("seed = %lu\n", gsl_rng_default_seed);

	gsl_rng *w= gsl_rng_alloc(gsl_rng_mt19937);
	
	// Populando a rede com valores aleatorios
	for(i= 0; i< Nx; i++){
		for(j= 0; j< Ny; j++){
			phi[i*Ny+j]= gsl_rng_uniform(w);
		}
	}
	gsl_rng_free(w);

	while(t < tf){
		
		// Loop k_1
		for(i= 0; i< Nx; i++){
			for(j= 0; j< Ny; j++){

				// Calculando integral de phi
				P_1 = canetaazul(i, j, phi);

				// Calculando k_1
                		k_1[(i)*Ny+j] = D * (phi[((i+1)%Nx)*Ny+j] + phi[((i-1+Nx)%Nx)*Ny+j] + phi[i*Ny+(j+1)%Ny] + phi[i*Ny+(j-1+Ny)%Ny] - 4.0*phi[i*Ny+j]) + mu*phi[i*Ny+j] - P_1*phi[i*Ny+j]/N_S;

			}
		}
	
		// Loop k_2
       		for(i= 0; i< Nx; i++){
			for(j= 0; j< Ny; j++){

				// Calculando integral de phi
				P_1 = canetaazul23(i, j, phi, k_1);

                		k_2[(i)*Ny+j] = D * (phi[((i+1)%Nx)*Ny+j] + 0.5*dt*k_1[((i+1)%Nx)*Ny+j]  + phi[((i-1+Nx)%Nx)*Ny+j] + 0.5*dt*k_1[((i-1+Nx)%Nx)*Ny+j]  + phi[i*Ny+(j+1)%Ny]  + 0.5*dt*k_1[i*Ny+(j+1)%Ny] + phi[i*Ny+(j-1+Ny)%Ny] + 0.5*dt*k_1[i*Ny+(j-1+Ny)%Ny] - 4.0*(phi[i*Ny+j] + 0.5*dt*k_1[i*Ny+j])) + mu*(phi[i*Ny+j] +0.5*dt*k_1[i*Ny+j]) - P_1*((phi[i*Ny+j] + 0.5*dt*k_1[i*Ny+j])/N_S);


			}
		} 

		// Calculando k_3
		for(i= 0; i< Nx; i++){
			for(j= 0; j< Ny; j++){
				P_1 = canetaazul23(i, j, phi, k_2);

                		k_3[(i)*Ny+j] = D * (phi[((i+1)%Nx)*Ny+j] + 0.5*dt*k_2[((i+1)%Nx)*Ny+j] + phi[((i-1+Nx)%Nx)*Ny+j] + 0.5*dt*k_2[((i-1+Nx)%Nx)*Ny+j] + phi[i*Ny+(j+1)%Ny] + 0.5*dt*k_2[i*Ny+(j+1)%Ny] + phi[i*Ny+(j-1+Ny)%Ny] + 0.5*dt*k_2[i*Ny+(j-1+Ny)%Ny] - (4.0*phi[i*Ny+j] + 0.5*dt*k_2[i*Ny+j])) + mu*(phi[i*Ny+j] + 0.5*dt*k_2[i*Ny+j]) - P_1*((phi[i*Ny+j] + 0.5*dt*k_2[i*Ny+j])/N_S);


			}
		}

		// Calculando k_4
        	for(i= 0; i< Nx; i++){
			for(j= 0; j< Ny; j++){
				P_1 = canetaazul4(i, j, phi, k_3);

                	k_4[(i)*Ny+j] = D * (phi[((i+1)%Nx)*Ny+j] + dt*k_3[((i+1)%Nx)*Ny+j]  + phi[((i-1+Nx)%Nx)*Ny+j] + dt*k_3[((i-1+Nx)%Nx)*Ny+j]  + phi[i*Ny+(j+1)%Ny]  + dt*k_3[i*Ny+(j+1)%Ny]  + phi[i*Ny+(j-1+Ny)%Ny] + dt*k_3[i*Ny+(j-1+Ny)%Ny] - 4.0*(phi[i*Ny+j] + dt*k_3[i*Ny+j])) + mu*(phi[i*Ny+j] + dt*k_3[i*Ny+j]) - P_1*((phi[i*Ny+j] + dt*k_3[i*Ny+j])/N_S);

                

			}
		}

		for(i=0; i<Nx; i++){
			for(j=0; j<Ny; j++){
				phi[(i)*Ny+j] = phi[i*Ny+j] + (dt/6.0) * (k_1[i*Ny+j] + 2.0*k_2[i*Ny+j] + 2.0*k_3[i*Ny+j] + k_4[i*Ny+j]);
			}
		}
		

		t+= dt;
		if(t >= a0){
			printf("%d %g\n", l, t);
			op(l++, phi);
			a0+= r;
		}
	}
	
	
	free(phi);
    	free(k_1);
    	free(k_2);
    	free(k_3);
    	free(k_4);
	return 0;
}

double canetaazul(int x, int y, double *p_1){
	int i, j;
	double P_1 = 0.0;
	for(i= 0; i< 20; i++){
		for(j= 0; j< 20; j++){
			if(sqrt((i-10)*(i-10)+(j-10)*(j-10)) < 10){
				P_1+= p_1[((x+i-10+Nx)%Nx)*Ny+((y+j-10+Ny)%Ny)];
				
			}
		}
	}
	return P_1;
}

double canetaazul23(int x, int y, double *p_1, double *k_n){
	int i, j;
	double P_1 = 0.0;
	for(i= 0; i< 20; i++){
		for(j= 0; j< 20; j++){
			if(sqrt((i-10)*(i-10)+(j-10)*(j-10)) < 10){
				P_1+= p_1[((x+i-10+Nx)%Nx)*Ny+((y+j-10+Ny)%Ny)] + 0.5*dt*k_n[((x+i-10+Nx)%Nx)*Ny+((y+j-10+Ny)%Ny)];  
			  //phi[((i+1)%Nx)*Ny+j] + 0.5*dt*k_1[((i+1)%Nx)*Ny+j]
				
			}
		}
	}
	return P_1;
}

double canetaazul4(int x, int y, double *p_1, double *k_n){
	int i, j;
	double P_1 = 0.0;
	for(i= 0; i< 20; i++){
		for(j= 0; j< 20; j++){
			if(sqrt((i-10)*(i-10)+(j-10)*(j-10)) < 10){
				P_1+= p_1[((x+i-10+Nx)%Nx)*Ny+((y+j-10+Ny)%Ny)] + dt*k_n[((x+i-10+Nx)%Nx)*Ny+((y+j-10+Ny)%Ny)];  
			  //phi[((i+1)%Nx)*Ny+j] + 0.5*dt*k_1[((i+1)%Nx)*Ny+j]
				
			}
		}
	}
	return P_1;
}
