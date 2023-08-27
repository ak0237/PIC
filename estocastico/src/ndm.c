#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>

#define NI 15000 // numero de individuos
#define TR 1.0 // tamanho da rede
#define R 0.1 // raio de interacao
#define NG 500 // numero de geracoes
#define NF 10 // numero de arquivos gerados

#define LS 0

#define alpha -1.0 // Se alpha/Ns for positivo, taxa de morte diminui. O inverso se for negativo
#define Ns 200 // Parametro de saturacao. Se positivo, penaliza reproducao em ambientes saturados
#define RpT 1.7 // Taxa de natalidade
#define MpT -0.7 // Taxa de mortalidade 
#define sigma 0.00447213595499957939
#define eta 1.0

struct INDV{
	int ezst; // 1 = individuo vivo | 0 = vazio ou individuo morto
	double x; // posicao x do individuo
	double y; // posicao y do individuo
};


#define r (double)NG/NF
#define q pow(NG, 1.0/(NF-1))

void ic(const gsl_rng *w, struct INDV *in){
	int c;
	
	// Criando 100 individuos em posicoes aleatorias
	for(c = 0; c < 100; c++){
		in[c].ezst = 1;
		in[c].x= gsl_rng_uniform(w)*TR;
		in[c].y= gsl_rng_uniform(w)*TR;
	}
}

void op(int t, struct INDV *in){
	int c;
	char name[100];
	FILE *file;
	
	sprintf(name, "dat/p-%d.dat", t);
	file = fopen(name, "w");
	for(c = 0; c < NI; c++){
		if(in[c].ezst != 0)
			fprintf(file, "%e %e %d\n", in[c].x, in[c].y, 1);
	}
	fclose(file);
}

double max(double x){
	if(x > 0){
		return x;
	}//else{
		return 0;
	//}
}

int main(int argc, char **argv){
	int i, j, n, t, N_R;
	int l= 0;
	double dx, dy, act, passo, angulo, p_rep;
	double c_rep, c_morte;
	struct INDV *in;
	
#if LS == 1
	double a0= 1.0;
#else
	double a0= r;
#endif
	
	in  = (struct INDV *) calloc(NI,    sizeof(struct INDV));
	
	// Definido a seed 
	gsl_rng_default_seed= (argc == 2) ? atoi(argv[1]) : time(NULL);
	
	// Criando um ponteiro para guardar um número aleatório
	gsl_rng *w= gsl_rng_alloc(gsl_rng_taus);
	
	// Definindo os parametros iniciais e exportando o primeiro arquivo
	ic(w, in);
	op(0, in);
	
	
	for(t = 1; t <= NG; t++){
		i = 0;
		while(i < NI){
			// selecionando individuo aleatoriamente
			do{
				n = gsl_rng_uniform(w)*NI;
			}while(in[n].ezst == 0);
			i++;
			
			
			
			// Verificando os individuos ao redor
			N_R= 0;
			for(j = 0; j < NI; j++){
				dx= fabs(in[j].x-in[n].x); // Distancia em x entre o individuo n e j
				if(dx > 0.5*TR){dx -= TR;} // Condicao periodica de contorno
				dy= fabs(in[j].y-in[n].y); // Distancia em y entre o individuo n e j
				if(dy > 0.5*TR){dy -= TR;} // Condicao periodica de contorno
				if(in[j].ezst == 1 && sqrt(dx*dx+dy*dy) < R){ // Verificando quantos individuos j estao pertos de n
					N_R++;
				}
			}
			
			// Definindo se n vai morrer, se reproduzir, ou nenhum 
			c_rep = max(RpT - N_R/Ns); // chance de reproducao
			c_morte  = max(MpT - alpha*N_R/Ns); // chance de morrer
			if(c_rep + c_morte > 0){
				p_rep= (c_rep)/(c_morte+c_rep); // Probabilidade de reproducao
				act= gsl_rng_uniform(w);
				if(act < p_rep){
					// Reproducao
					for(j = 0; j < NI; j++){
						if(in[j].ezst == 0 ){
							in[j].ezst = 1;
							in[j].x = in[n].x; 
							in[j].y = in[n].y; 
							j=NI+1;
						}
					}
				}else{
					// Morte
					in[n].ezst= 0;
				}
			}
			
			// Movimentando o individuo
			angulo = M_PI*(gsl_rng_uniform(w)-0.5)*eta; // Angulo
			passo = fabs(gsl_ran_gaussian(w, sigma)); // Tamanho do passo
			in[n].x += passo*cos(angulo); // Passo em x
			if(in[n].x > TR){in[n].x -= TR;} // Condicao periodica de contorno
			if(in[n].x < 0.0){in[n].x += TR;} // Condicao periodica de contorno
			in[n].y += passo*sin(angulo); // Passo em y
			if(in[n].y > TR){in[n].y -= TR;} // Condicao periodica de contorno
			if(in[n].y < 0.0){in[n].y += TR;} // Condicao periodica de contorno
		}
		// Contando quantos individuos existem por geracao
		int aa= 0;
		for(j= 0; j< NI; j++){
			if(in[j].ezst == 1){
				aa++;
			}
		}
		printf("%d %d\n", t, aa);
		
		// Registrando as mudancas de forma igualmente espaca no tempo
		if(t >= round(a0)){
#if LS == 1
			a0*= q;
#else 
			a0+= r;
#endif
			op(++l, in);

		}
	}
	
	gsl_rng_free(w);
	free(in);
	return 0;	
}









