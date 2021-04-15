/*
 Main program
 NTAILS->T
 n->N
 NPATHS->M
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <memory.h>
#include <math.h>
#include <float.h>
#include <ipps.h>
#include <iostream>
#include "engine.h"
#include "coder.h"
#include "bitio.h"
#include "hmm.h"
using namespace std;

//#include <intrin.h>  
//#include <Windows.h>  
#pragma comment(lib,"libeng.lib")
#pragma comment(lib,"ippi.lib")
/*
 * Declarations for local routines.
 */
void source(char type[]);
void compress(void);
void expand(void);
unsigned expandHMM(float p,float k);
void error_exit(char *message);

double p = 0.5, qs = 0.110028, qc = 0.0;
double ys[100001] = { 0.0 };
int cell = 100000;

//NTAILS->T
#define NTAILS	15U
unsigned src[BUF_SIZE], side[BUF_SIZE], error[BUF_SIZE];
//n->N
unsigned n = 1024, m = 0;


int	*O;		/* the observation sequence O[1..T]*/
int	*q; 	/* the state sequence q[1..T] */
HMM hmm; 	/* the HMM */
//NPATHS->M
#define NPATHS	2048U
struct DAC_node{
	unsigned low, high, code, ptr;
	double metric;
	double exmetric;
	double alpha[2];
	unsigned path[BUF_SIZE];
};

DAC_node node[NPATHS<<1]; 
unsigned nodeIdx[NPATHS<<1];

extern "C"{
	int *ivector(int nl, int nh);
	void free_ivector(int *v, int nl, int nh);
}

/*
 * This example program compresses an input string, sending
 * the output to a file.  It then expands the output file,
 * sending the decoded characters to the screen.
 */


void getDataCss(float k) {
	char nums[10] ;
	gcvt(k,2, nums);
	char path[20] = {};
	strcat(path, "data//");
	strcat(path, nums);
	strcat(path, ".txt\0");
	FILE * fp = fopen(path, "r");
	if (fp == NULL)
	{
		printf("file open error/n");
	}
	for (int i = 0; i <= 100000; i++)
		fscanf(fp, "%lf", &ys[i]);
	fclose(fp);
	printf("trying k=%lf , path=%s , data read ok\n",k,path);
}


int main(int argc, char *argv[]){

	srand(0);
	initLUT(p);

	// alloc HMM-related
	FILE *fp = fopen(argv[1], "rt");		
	if(fp==NULL){
		fprintf(stderr, "Error: File %s not found \n", argv[1]);
		exit (1);
	}
	ReadHMM(fp, &hmm);
	fclose(fp);
	O = ivector(1, n); /* alloc space for observation sequence O */
	q = ivector(1, n); /* alloc space for state sequence q */
	
	unsigned ntries = 100;
	float rate = 0.0f;
	unsigned nerr = 0;
	float r = 0.0f;
	// loop for trials
	for (unsigned i=0; i<ntries; i++){
		
		source("hmm");
		float k = atof(argv[2]);		
		while(1){// loop for rates
			getDataCss(k);
			initLUTdac(p, k);
			compress();
			nerr = expandHMM(p, k);
			if (nerr == 0) {
				//puts("Good!");				
				break;
			}
			k += 0.01000;
		}
		r += k;
		printf("CR %d, %f\n", m, (float)m/(float)n);
		rate += float(m)/ float(n);
		printf("\nk=%f , %d/%d decoding ok\n\n", k,i + 1,ntries);
		
	}
	printf("mean R : %lf\n", r/ntries );
	printf("(results of table II): %lf\n", rate/ntries);

	// free HMM-related 
	free_ivector(O, 1, n);
	free_ivector(q, 1, n);
	FreeHMM(&hmm);

	getchar();
	return 0;
}

void display(unsigned *s, unsigned n, char *name){
	puts(name);
	for(unsigned i=0; i<n; i++)		
		putchar(getbit(s, i)+0x30);
	putchar('\n');
	putchar('\n');
}
/*
 * This is the compress routine.  It shows the basic algorithm for
 * the compression programs used in this article.  First, an input
 * characters is loaded.  The modeling routines are called to
 * convert the character to a symbol, which has a high, low and
 * range.  Finally, the arithmetic coder module is called to
 * output the symbols to the bit stream.
 */
void source(char type[]){
	GenSequenceArray(&hmm, 0,n, O, q);//zy
	memset(src, 0, BUF_SIZE*sizeof(unsigned));
	memset(side, 0, BUF_SIZE*sizeof(unsigned));
	for(unsigned i=0; i<n; i++){
		byte c = (rand()<RAND_MAX*p);
		if(c)	
			putone(src, i);
		if(!memcmp(type, "bsc", 3))			
			c ^= (rand()<RAND_MAX*qs);
		else if(!memcmp(type, "hmm", 3))	
			c ^= (O[i+1]-1);
		if(c)	
			putone(side, i);
	}
	
}

void compress(void){
	//puts("Compressing...");
	m = 0;
	init_AE();
	for(unsigned i=0; i<n; i++){
		byte c = getbit(src, i);
		encode_symbol(m, c, i<(n-NTAILS));
    }
	flush_AE(m);
}

int compare(const void *a, const void *b){
	unsigned idxA = *((unsigned*)a), idxB = *((unsigned*)b);
	if(node[idxA].metric > node[idxB].metric)		return -1;
	else if(node[idxA].metric == node[idxB].metric)	return 0;  
	else							 				return 1;
}

/*
 * This expansion routine demonstrates the basic algorithm used for
 * decompression in this article.  It first goes to the modeling
 * module and gets the scale for the current context.  (Note that
 * the scale is fixed here, since this is not an adaptive model).
 * It then asks the arithmetic decoder to give a high and low
 * value for the current input number scaled to match the current
 * range.  Finally, it asks the modeling unit to convert the
 * high and low values to a symbol.
 */


 unsigned expandHMM(float p, float k){
	memset(node, 0, sizeof(node));   
	for(unsigned i=0; i<(NPATHS<<1); i++)	
		nodeIdx[i] = i;
	init_AD(node[0].low, node[0].high, node[0].code, node[0].ptr);
	unsigned nNode = 1, nSubNode;

	for(unsigned t=0; t<n; t++){
		bool dac = t<(n-NTAILS);
		byte s = getbit(side, t);
		nSubNode = nNode;
		for(unsigned j=0; j<nNode; j++){// for each active node
			unsigned jj = nodeIdx[j];   
			
			// forward algorithm
			double scale[2] = {0.0, 0.0}, alpha[2][2];
			if(!t){                                   
				for(int i=0; i<hmm.N; i++){
					alpha[0][i] = hmm.pi[i+1] * hmm.B[i+1][0+1];
					scale[0] += alpha[0][i];
					alpha[1][i] = hmm.pi[i+1] * hmm.B[i+1][1+1];
					scale[1] += alpha[1][i];
				}
			}else{
				for(int j=0; j<hmm.N; j++) {
					double sum = 0.0;
					for(int i=0; i<hmm.N; i++)	
						sum += node[jj].alpha[i] * hmm.A[i+1][j+1]; 
					alpha[0][j] = sum * hmm.B[j+1][0+1];
					scale[0] += alpha[0][j];             

					alpha[1][j] = sum * hmm.B[j+1][1+1];
					scale[1] += alpha[1][j];              
				}
			}
			for(int i=0; i<hmm.N; i++){ 
				alpha[0][i] /= scale[0]; 
				alpha[1][i] /= scale[1]; 
			}
			// end of forward algorithm

			// dac dec
			double weight0 = 0.0, weight1 = 0.0;
			byte c = get_symbol(node[jj].low, node[jj].high, node[jj].code, dac, &weight0, &weight1,k );
			if(c=='x'){				
				unsigned kk = nodeIdx[nSubNode];
				memcpy(node+kk, node+jj, sizeof(DAC_node));	// 
				if (s) {
					putone(node[jj].path, t);	// set the '1' branch to node[jj]	
				}								
				else {
					putone(node[kk].path, t);	// set the '1' branch to node[kk]	
				}			

				node[jj].metric += log(scale[0]);
				node[kk].metric += log(scale[1]);
	
				/*
				* If you want to use the forward algorithm and the initial spectrum as an overall measure
				* 
				if (t<n-32) {
					if (s) {
						node[jj].metric += log(weight1);
						node[kk].metric += log(weight0);
					}
					else {
						node[jj].metric += log(weight0);
						node[kk].metric += log(weight1);
					}
				}
				*/

				memcpy(node[jj].alpha, alpha[0], hmm.N*sizeof(double));
				memcpy(node[kk].alpha, alpha[1], hmm.N*sizeof(double));
				remove_symbol(node[jj].low, node[jj].high, node[jj].code, node[jj].ptr, s, dac);
				remove_symbol(node[kk].low, node[kk].high, node[kk].code, node[kk].ptr, !s, dac);
				nSubNode++;
			}else{
				if(c)	
					putone(node[jj].path, t);
				node[jj].metric += log(scale[c^s]);
				memcpy(node[jj].alpha, alpha[c^s], hmm.N*sizeof(double));
				remove_symbol(node[jj].low, node[jj].high, node[jj].code, node[jj].ptr, c, dac);
			}
		}

		// sort 
		if(nSubNode>1){	
			qsort(nodeIdx, nSubNode, sizeof(unsigned), compare);
			nNode = __min(NPATHS, nSubNode);		// select (<=NPATHS) paths 
		}
	}

	unsigned *dec = node[nodeIdx[0]].path;
	
	unsigned nerr = 0;
	for (unsigned i = 0; i < n; i++) {
		nerr += (getbit(dec, i) ^ getbit(src, i));
	}
	return nerr;

	
}


/*
 * A generic error routine.
 */
void error_exit( char *message )
{
    puts( message );
    exit( -1 );
}
