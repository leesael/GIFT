/*
* @file        Silenced-TF.cpp
* @author      Sejoon Oh (ohhenrie@snu.ac.kr), Seoul National University
* @author      Jungwoo Lee (muon9401@gmail.com), Seoul National University
* @author      Lee Sael (saellee@gmail.com), Seoul National University
* @version     1.0
* @date        2018-04-07
*
* GIFT: Guided and Interpretable Factorization for Tensors - Applications to Human Cancer Analytics
*
* This software is free of charge under research purposes.
* For commercial purposes, please contact the author.
*
* Description: the interpretable TF method which makes unintended entries as zeros according to mask matrices. 
* Usage:
*   - make stf
*   - ./stf input.txt mask.txt result/ 3 0.01 10 20 30
*/


/////    Header files     /////

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <armadillo>
#include <omp.h>
#define ARMA_USE_LAPACK
#define ARMA_USE_BLAS
using namespace std;
using namespace arma;
double lambda;
///////////////////////////////
char* InputPath;
char* ResultPath;
char* testPath;
char* MaskPath;
/////////////////////////////////                Variables                  /////////////////////////////////

//Related to input tensor
//Sparse tensor form : (i1, i2, ... , iN, value)
int order;										// Tensor order (e.g., 5)
int *dimensionality, max_dim;					// Tensor dimensionality (e.g., 100x100x100)
int Entries_N;									// Total number of observation entries
double *Entries;								// Containing all values of an input tensor
int *Index;										// Containing all indices of an input tensor (Vectorized form) 
int *WhereX, *CountX;					        // WhereX[n][I_n] contains all entries of a tensor X whose nth mode's index is I_n
double NormX;									// Norm of input tensor
int threadsN;									// Number of threads to be used
int samp;

int Core_N = 1;									// Total number of nonzeros in a core tensor
int *Core_size, Core_dim;						// Core tensor dimensionality
double *CoreTensor;								// Containing all values of a core tensor
int *CorePermu;									// Containing all indices of a core tensor (Vectorized form) 
double MaxCore;									// The maximum value of core tensor
int *WhereG, *CountG;							// WhereG[n][R_n] contains all entries of a core tensor G whose nth mode's index is R_n

double *FactorM; //Factor matrices in vectorized form


//Update-related variables
//Please see the paper and supplementary material for details
double Fit, pFit = -1;
double *Error_T;
double Error, RMSE, MAE, TEST_RMSE;
int *crows, rowcount;
double *tempCore;
int *Mul, *tempPermu;
//For-loop variables
int i, j, k, l;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


double frand(double x, double y) {//return the random value in (x,y) interval
	return ((y - x)*((double)rand() / RAND_MAX)) + x;
}

double abss(double x) {
	return x > 0 ? x : -x;
}

//[Input] Given tensor X
//[Output] Updated WhereX and CountX
//[Function] Assign all non-zeros to the corresponding rows, represented as WhereX and CountX
void assign_index() {
	int *tempX = (int *)malloc(sizeof(int)*max_dim*order);
	int pos = 0, i, j, k, l;
	for (i = 0; i < order; i++) {
		for (j = 0; j < dimensionality[i]; j++) {
			CountX[i*max_dim + j] = tempX[i*max_dim + j] = 0;
		}
	}
	for (i = 0; i < Entries_N; i++) {
		for (j = 0; j < order; j++) {
			k = Index[pos++];
			CountX[j*max_dim + k]++;
			tempX[j*max_dim + k]++;
		}
	}
	pos = 0;
	int now = 0;
	for (i = 0; i < order; i++) {
		pos = i*max_dim;
		for (j = 0; j < dimensionality[i]; j++) {
			k = CountX[pos];
			CountX[pos] = now;
			tempX[pos++] = now;
			now += k;
		}
		CountX[pos] = now;
		tempX[pos] = now;
	}
	pos = 0;
	for (i = 0; i < Entries_N; i++) {
		for (j = 0; j < order; j++) {
			k = Index[pos++];
			int now = tempX[j*max_dim + k];
			WhereX[now] = i;
			tempX[j*max_dim + k]++;
		}
	}
}

FILE *fmask;
int *Mask;
char tmp[1005];

//[Input] Metadata + input tensor as a sparse tensor format 
//[Output] Initialized core tensor G and factor matrices A^{(n)} (n=1...N)
//[Function] Getting all information about an input tensor X / Initialize all factor matrices and core tensor.
void Getting_Input() {
	FILE *fin = fopen(InputPath, "r");
	FILE *fin2 = fopen(InputPath, "r");
	
	double Timee = clock();

	printf("Reading Input...\n");

	dimensionality = (int *)malloc(sizeof(int)*order);

	for (i = 0; i < order; i++) {
		dimensionality[i] = 0;
		Core_N *= Core_size[i];
		if(Core_dim<Core_size[i]) Core_dim = Core_size[i];
	}

	int len = 0;
	while (fgets(tmp, 1005, fin)) {
		Entries_N++;
	}
	Index = (int *)malloc(sizeof(int)*Entries_N*order);
	Entries = (double *)malloc(sizeof(double)*Entries_N);
	int pos = 0;
	for (i = 0; i < Entries_N; i++) {
		fgets(tmp, 1005, fin2);
		len = strlen(tmp);
		int k = 0, idx = 0, flag = 0;
		double mul = 0.1, val = 0;
		for (j = 0; j < len; j++) {
			if (tmp[j] == ' ' || tmp[j] == '\t') {
				Index[pos++] = idx - 1;
				if (dimensionality[k] < idx) dimensionality[k] = idx;
				idx = 0;
				k++;
			}
			else if (tmp[j] >= '0' && tmp[j] <= '9') {
				if (flag == 1) {
					val += mul*(tmp[j] - '0');
					mul /= 10;
				}
				else idx = idx * 10 + tmp[j] - '0';
			}
			else if (tmp[j] == '.') {
				val += idx;
				flag = 1;
			}
		}
		if(flag==0) val = idx;
		Entries[i] = val;
		NormX += Entries[i] * Entries[i];
	}
	
	for (i = 0; i < order; i++) {
		if (max_dim < dimensionality[i]) max_dim = dimensionality[i];
	}
	max_dim++;
	WhereX = (int *)malloc(sizeof(int)*order*Entries_N);
	CountX = (int *)malloc(sizeof(int)*max_dim*order);

	NormX = sqrt(NormX);
	printf("Reading Done.\nTensor Order: %d\tSize: ", order);
	for (i = 0; i < order; i++) {
		if (i != order - 1) printf("%dx", dimensionality[i]);
		else printf("%d\t", dimensionality[i]);
	}
printf("]\tRank: [");
	for (i = 0; i < order; i++) {
		if (i != order - 1) printf("%dx", Core_size[i]);
		else printf("%d", Core_size[i]);
	}
	printf("]\tNNZ : %d\tThreads : %d\tNorm : %lf\tLambda : %lf\nInitialize...\n", Entries_N, omp_get_max_threads(), NormX, lambda);

	//INITIALIZE
	assign_index();
	FactorM = (double *)malloc(sizeof(double)*order*max_dim*Core_dim);
	Mask = (int *)malloc(sizeof(int)*order*max_dim*Core_dim);
	memset(Mask, 0, sizeof(int)*order*max_dim*Core_dim);
	fmask = fopen(MaskPath, "r");
	int xx, yy,zz;
	while(fscanf(fmask,"%d %d %d",&xx,&yy,&zz)!=EOF){
		Mask[(xx-1)*max_dim*Core_dim+(yy-1)*Core_dim+(zz-1)] = 1;
	}

	int mult = max_dim*Core_dim;
	for (i = 0; i < order; i++) {
		int row = dimensionality[i], col = Core_size[i];
		for (j = 0; j < row; j++) {
			for (k = 0; k < col; k++) {
				FactorM[i*mult + j*col + k] = frand(0, 1);
			}
		}
	}

	CoreTensor = (double *)malloc(sizeof(double)*Core_N);
	CorePermu = (int *)malloc(sizeof(int)*Core_N*order);
	pos = 0;
	for (i = 0; i < Core_N; i++) {
		CoreTensor[i] = frand(0, 1);
		if (i == 0) {
			for (j = 0; j < order; j++) CorePermu[j] = 0;
		}
		else {
			for (j = 0; j < order; j++) {
				CorePermu[i*order + j] = CorePermu[(i - 1)*order + j];
			}
			CorePermu[i*order + order - 1]++;  k = order - 1;
			while (CorePermu[i*order + k] >= Core_size[k]) {
				CorePermu[i*order + k] -= Core_size[k];
				CorePermu[i*order + k - 1]++; k--;
			}
		}
	}
	printf("Elapsed Time for Getting Input and INIT:\t%lf\n", (clock() - Timee) / CLOCKS_PER_SEC);
}


//[Input] Input tensor X, initialized core tensor G, and factor matrices A^{(n)} (n=1...N)  
//[Output] Updated factor matrices A^{(n)} (n=1...N)
//[Function] Update all factor matrices by a gradient-based ALS.
void Update_Factor_Matrices() {

	int mult = max_dim*Core_dim;
	for (i = 0; i < order; i++) { //Updating the ith Factor Matrix
		int row_size = dimensionality[i];
		int column_size = Core_size[i];
		int iter;
		rowcount = 0;
		crows = (int *)malloc(sizeof(int)*row_size);
		for (j = 0; j < row_size; j++) { //extracting rows who have non-zeros 
			int pos = i*max_dim + j, nnz = CountX[pos + 1] - CountX[pos];
			if (nnz != 0) {
				crows[rowcount++] = j;
			}
		}

#pragma omp parallel for schedule(dynamic) //in parallel
		for (iter = 0; iter < rowcount; iter++) {
			int j, k, l, ii, jj;
			j = crows[iter]; //Indicating that we are updating the jth row of ith factor matrix
			double *Delta = (double *)malloc(sizeof(double)*column_size);
			double *B = (double *)malloc(sizeof(double)*column_size*column_size);
			double *C = (double *)malloc(sizeof(double)*column_size);

			//Initialize B and C
			int pos = 0;
			for (k = 0; k < column_size; k++) {
				for (l = 0; l < column_size; l++) {
					B[pos] = 0;
					if (k == l){
						B[pos] = lambda;
					}
					pos++;
				}
				C[k] = 0;
			}

			pos = i*max_dim + j;
			int nnz = CountX[pos + 1] - CountX[pos];
			pos = CountX[pos];
			for (k = 0; k < nnz; k++) { //Updating Delta, B, and C
				int current_input_entry = WhereX[pos + k];
				int pre_val = current_input_entry*order;
				int *cach1 = (int *)malloc(sizeof(int)*order);
				for (l = 0; l < order; l++) cach1[l] = Index[pre_val++];
				for (l = 0; l < column_size; l++) Delta[l] = 0;
				for (l = 0; l < Core_N; l++) {
					int pre1 = l*order, pre2 = 0;
					int CorePos = CorePermu[pre1 + i];
					double res = CoreTensor[l];
					for (ii = 0; ii < order; ii++) {
						if (ii != i) {
							int mulrow = cach1[ii], mulcol = CorePermu[pre1];
							res *= FactorM[pre2 + mulrow*Core_size[ii] + mulcol];
						}
						pre1++;
						pre2 += mult;
					}
					Delta[CorePos] += res;
				}
				free(cach1);
				int now = 0;
				double Entry_val = Entries[current_input_entry];
				for (ii = 0; ii < column_size; ii++) {
					double cach = Delta[ii];
					for (jj = 0; jj < column_size; jj++) {
						B[now++] += cach * Delta[jj];
					}
					C[ii] += cach * Entry_val;
				}
			}
			free(Delta);
			//Getting the inverse matrix of [B+lambda*I]
			mat A = mat(column_size, column_size);
			pos = 0;
			for (k = 0; k < column_size; k++) {
				for (l = 0; l < column_size; l++) {
					A(k, l) = B[k*column_size + l];
				}
			}

			mat BB = inv_sympd(A);
			pos = 0;
			for (k = 0; k < column_size; k++) {
				for (l = 0; l < column_size; l++) {
					B[k*column_size + l] = BB(k, l);
				}
			}

			//Update the jth row of ith Factor Matrix 
			int cach = i*mult + j*column_size;
			for (k = 0; k < column_size; k++) {
				double res = 0;
				for (l = 0; l < column_size; l++) {
					res += C[l] * B[l*column_size + k];
				}
				if (Mask[i*mult+j*Core_dim+k] == 1) res=0;
				FactorM[cach + k] = res;
			}
			free(B);
			free(C);
		}
		free(crows);
	}
}


//[Input] Input tensor X, core tensor G, and factor matrices A^{(n)} (n=1...N)
//[Output] Fit = 1-||X-X'||/||X|| (Reconstruction error = ||X-X'||)
//[Function] Calculating fit and reconstruction error in a parallel way.
void Reconstruction() {
	RMSE = Error = 0;
	Error_T = (double *)malloc(sizeof(double)*Entries_N);
#pragma omp parallel for schedule(static)
	for (i = 0; i < Entries_N; i++) {
		Error_T[i] = Entries[i];
	}
	int mult = max_dim*Core_dim;
#pragma omp parallel for schedule(static)
	for (i = 0; i < Entries_N; i++) {
		int j, pre_val = i*order;
		double ans = 0;
		int *cach1 = (int *)malloc(sizeof(int)*order);
		for (j = 0; j < order; j++) cach1[j] = Index[pre_val++];
		for (j = 0; j < Core_N; j++) {
			double temp = CoreTensor[j];
			int k;
			int pos = j*order;
			int val = 0;
			for (k = 0; k < order; k++) {
				int mulrow = cach1[k], mulcol = CorePermu[pos++];
				temp *= FactorM[val + mulrow*Core_size[k] + mulcol];
				val += mult;
			}
			ans += temp;
		}
		free(cach1);
		Error_T[i] -= ans;
	}
#pragma omp parallel for schedule(static) reduction(+:Error)
	for (i = 0; i < Entries_N; i++) {
		Error += Error_T[i] * Error_T[i];
	}
	RMSE = sqrt(Error / Entries_N);
	if (NormX == 0) Fit = 1;
	else Fit = 1 - sqrt(Error) / NormX;
	free(Error_T);
}


//[Input] Input tensor X, initialized core tensor G, and initialized factor matrices A^{(n)} (n=1...N)
//[Output] Updated core tensor G and factor matrices A^{(n)} (n=1...N)
//[Function] Performing main algorithm which updates core tensor and factor matrices iteratively
void PTucker() {
	printf("PTucker START\n");

	double Stime = omp_get_wtime();
	int iter = 0;
	double avertime = 0;

	while (1) {

		double itertime = omp_get_wtime(), steptime;
		steptime = itertime;

		Update_Factor_Matrices();
		printf("Factor Time : %lf\n", omp_get_wtime() - steptime);
		steptime = omp_get_wtime();

		Reconstruction();
		printf("Recon Time : %lf\n", omp_get_wtime() - steptime);
		steptime = omp_get_wtime();

	
		avertime += omp_get_wtime() - itertime;
		printf("iter%d :      Fit : %lf\tTraining RMSE: %lf\tTest RMSE: %lf\tElapsed Time : %lf\n", ++iter, Fit, RMSE, TEST_RMSE, omp_get_wtime() - itertime);

		if (pFit!=-1 && abss(pFit-Fit)<=0.001) break;
		pFit = Fit;
	}

	avertime /= iter;

	printf("\nTotal update ended.\tFinal fit : %lf\tAverage iteration time : %lf\tTotal Elapsed time: %lf\n", Fit,avertime, omp_get_wtime() - Stime);
}


//[Input] Updated core tensor G and factor matrices A^{(n)} (n=1...N)
//[Output] core tensor G in sparse tensor format and factor matrices A^{(n)} (n=1...N) in full-dense matrix format
//[Function] Writing all factor matrices and core tensor in result path
void Print() {
	printf("\nWriting factor matrices and core tensor to file...\n");
	char temp[50];
	int pos = 0;
	int mult = max_dim*Core_dim;
	for (i = 0; i < order; i++) {
		sprintf(temp, "%s/FACTOR%d", ResultPath, i + 1);
		FILE *fin = fopen(temp, "w");
		for (j = 0; j < dimensionality[i]; j++) {
			for (k = 0; k < Core_size[i]; k++) {
				fprintf(fin, "%e\t", FactorM[i*mult + j*Core_size[i] + k]);
			}
			fprintf(fin, "\n");
		}
	}
	sprintf(temp, "%s/CORETENSOR", ResultPath);
	FILE *fcore = fopen(temp, "w");
	pos = 0;
	for (i = 0; i < Core_N; i++) {
		for (j = 0; j < order; j++) {
			fprintf(fcore, "%d\t", CorePermu[pos++] + 1);
		}
		fprintf(fcore, "%e\n", CoreTensor[i]);
	}
}


//[Input] Path of input tensor file, result directory, tensor order, tensor rank, and number of threads
//[Output] Core tensor G and factor matrices A^{(n)} (n=1,...,N)
//[Function] Performing P-Tucker for a given sparse tensor
int main(int argc, char* argv[]) {
	srand(time(NULL));
	if (argc !=0) {
		InputPath = argv[1];
		MaskPath = argv[2];
		ResultPath = argv[3];
		order = atoi(argv[4]);
		lambda = atof(argv[5]);
		Core_size = (int *)malloc(sizeof(int)*order);
		for (int i = 6; i <= 5 + order; i++) {
			Core_size[i - 6] = atoi(argv[i]);
		}

		Getting_Input();

		PTucker();

		Print();

	}

	else printf("Usage: ./stf [input_tensor_path] [mask_matrix_path] [result_directory_path] [tensor_order] [lambda] [tensor_rank]\n");
	return 0;

}
