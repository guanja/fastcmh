#include<stdlib.h>
#include<stdio.h>
#include "../../EasyGWAS/chi2.h"

// Number of tables
long long K;
// Number of observations per table
long long *Nt;
// Number of observations in positive class per table
long long *nt;
// Now some precomputed quantities to save time when computing the CMH test statistic
long long *x;
long long *Nt_nt; // Ni-ni for each of the K tables
long long *hypercorner_bnd; // max(ni,Ni-ni) for each of the K tables
double *gammat; // ni/Ni for each of the K tables
double *gammabint; // (ni/Ni)*(1-ni/Ni) for each of the K tables
double *f_vals, *g_vals, *betas;
double f_sum, g_sum, Tcmh_max_corner_l, Tcmh_max_corner_r, Tcmh_aux_corner;
long long *idx_betas_sorted;

unsigned long long *bitmasks;

extern double Chi2_sf(double,double); //Chi2 Survival function, from Dominik's EasyGWASCore

int qsort_cmp_betas ( const void *x, const void *y ){
	if ( betas[*((long long *)x)] < betas[*((long long *)y)] ) return (-1);
	else return 1;
}

double isprunable2(long long*x){
	// Variables for looping across 2^K cases
	unsigned long long idx_mask, max_mask;
	// Variables for computing minimum attainable p-value in each case
	long long k;
	double left_tail_num, right_tail_num, den;
	double aux1, aux2;

	// If for any of the K tables, its margin x is smaller than the maximum of n and N-n, then we cannot prune
	// the interval (we are not in the "top-right" hypercorner)
	for(k=0; k<K; k++) if(x[k] < hypercorner_bnd[k]) return 0;

	max_mask = 1<<K; Tcmh_max_corner_l = 0;
	for(idx_mask=1;idx_mask<max_mask;idx_mask++){
		left_tail_num = 0, right_tail_num = 0, den = 0;
		for(k=0; k<K; k++){
			// Skip if the bitmask contains a zero
			if(idx_mask&bitmasks[k]){
				aux1 = x[k]-Nt_nt[k]; aux2 = x[k]*gammat[k];
				left_tail_num += ((aux1 > 0) ? aux1 : 0) - aux2;
				right_tail_num += ((x[k] > nt[k]) ? nt[k] : x[k]) - aux2;
				den += x[k]*(1-((double)x[k])/Nt[k])*gammabint[k];
			}
		}
		left_tail_num *= left_tail_num; right_tail_num *= right_tail_num;
		if(den==0) Tcmh_aux_corner = 1;
		else Tcmh_aux_corner = ((left_tail_num > right_tail_num) ? left_tail_num : right_tail_num)/den;
		Tcmh_max_corner_l = (Tcmh_max_corner_l >= Tcmh_aux_corner) ? Tcmh_max_corner_l : Tcmh_aux_corner;
	}
	return Chi2_sf(Tcmh_max_corner_l,1);
}

double isprunable(long long *x){
	long long j,k;
	// If for any of the K tables, its margin x is smaller than the maximum of n and N-n, then we cannot prune
	// the interval (we are not in the "top-right" hypercorner)
	for(k=0; k<K; k++) if(x[k] < hypercorner_bnd[k]) return 0;

	//printf("a\n");
	// Compute the maximum value of left handside function
	for(j=0,k=0; k<K; k++){
		// Discard all dimensions for which x[k]==Nt[k], as they don't contribute to the function neither in the
		// numerator nor the denominator
		if(x[k] < Nt[k]){
			f_vals[j] = gammat[k]*(Nt[k]-x[k]);
			g_vals[j] = gammabint[k]*x[k]*(1-((double)x[k])/Nt[k]);
			betas[j] = g_vals[j]/f_vals[j];
			//printf("f_vals[%d]=%f\n",j,f_vals[j]);
			//printf("g_vals[%d]=%f\n",j,g_vals[j]);
			//printf("betas[%d]=%f\n",j,betas[j]);
			idx_betas_sorted[j] = j;
			j++;
		}
	}
	//printf("b,j=%d\n",j);
	qsort(idx_betas_sorted,j,sizeof(long long),qsort_cmp_betas); //Equivalent to argsort(betas[0:j])
	f_sum = 0; g_sum = 0; Tcmh_max_corner_l = 0;
	for(k=0; k<j; k++){
		//printf("%d,%d\n",k,idx_betas_sorted[k]);
		f_sum += f_vals[idx_betas_sorted[k]];
		g_sum += g_vals[idx_betas_sorted[k]];
		//printf("betas[idx_sorted[%d]]=%f\n",k,betas[idx_betas_sorted[k]]);
		Tcmh_aux_corner = (f_sum*f_sum)/g_sum;
		Tcmh_max_corner_l = (Tcmh_max_corner_l >= Tcmh_aux_corner) ? Tcmh_max_corner_l : Tcmh_aux_corner; //Tcmh_max_corner_l=max(Tcmh_max_corner_l,Tcmh_aux_corner)
	}
	//printf("c\n");
	// Compute the maximum value of right handside function
	for(j=0,k=0; k<K; k++){
		// Discard all dimensions for which x[k]==Nt[k], as they don't contribute to the function neither in the
		// numerator nor the denominator
		if(x[k] < Nt[k]){
			f_vals[j] = (1-gammat[k])*(Nt[k]-x[k]);
			//printf("f_vals[%d]=%f\n",j,f_vals[j]);
			//g_vals doesn't change, hence it does not need to be recomputed
			betas[j] = g_vals[j]/f_vals[j];
			//printf("betas[%d]=%f\n",j,betas[j]);
			idx_betas_sorted[j] = j;
			j++;
		}
	}
	//printf("d\n");
	qsort(idx_betas_sorted,j,sizeof(long long),qsort_cmp_betas); //Equivalent to argsort(betas[0:j])
	f_sum = 0; g_sum = 0; Tcmh_max_corner_r = 0;
	//printf("e,j=%d\n",j);
	for(k=0; k<j; k++){
		//printf("%d,%d\n",k,idx_betas_sorted[k]);
		f_sum += f_vals[idx_betas_sorted[k]];
		g_sum += g_vals[idx_betas_sorted[k]];
		//printf("betas[idx_sorted[%d]]=%f\n",k,betas[idx_betas_sorted[k]]);
		Tcmh_aux_corner = (f_sum*f_sum)/g_sum;
		Tcmh_max_corner_r = (Tcmh_max_corner_r >= Tcmh_aux_corner) ? Tcmh_max_corner_r : Tcmh_aux_corner; //Tcmh_max_corner_r=max(Tcmh_max_corner_r,Tcmh_aux_corner)
	}
	//printf("f\n");
	printf("l=%f\n",Tcmh_max_corner_l); printf("r=%f\n",Tcmh_max_corner_r);
	printf("max=%f\n",((Tcmh_max_corner_r >= Tcmh_max_corner_l) ? Tcmh_max_corner_r : Tcmh_max_corner_l));
	return Chi2_sf(((Tcmh_max_corner_r >= Tcmh_max_corner_l) ? Tcmh_max_corner_r : Tcmh_max_corner_l),1);
}

int main(int argc, char *argv[]){
	long long j,k;
	K = (argc-1)/3;
	//printf("K=%d\n",K);

	x = (long long *)calloc(K, sizeof(long long));
	if(!x){
		fprintf(stderr,"Error in function sis_init: couldn't allocate memory for array x\n");
		exit(1);
	}
	nt = (long long *)calloc(K, sizeof(long long));
	if(!nt){
		fprintf(stderr,"Error in function sis_init: couldn't allocate memory for array nt\n");
		exit(1);
	}
	Nt = (long long *)calloc(K, sizeof(long long));
	if(!Nt){
		fprintf(stderr,"Error in function sis_init: couldn't allocate memory for array Nt\n");
		exit(1);
	}
	for(k=0; k<K; k++){
		x[k] = atoll(argv[1+3*k]);
		nt[k] = atoll(argv[2+3*k]);
		Nt[k] = atoll(argv[3+3*k]);
	}


	Nt_nt = (long long *)calloc(K, sizeof(long long));
		if(!Nt_nt){
			fprintf(stderr,"Error in function sis_init: couldn't allocate memory for array Nt_nt\n");
			exit(1);
	}
	hypercorner_bnd = (long long *)calloc(K, sizeof(long long));
	if(!hypercorner_bnd){
		fprintf(stderr,"Error in function sis_init: couldn't allocate memory for array hypercorner_bnd\n");
		exit(1);
	}
	gammat = (double *)calloc(K, sizeof(double));
	if(!gammat){
		fprintf(stderr,"Error in function sis_init: couldn't allocate memory for array gammat\n");
		exit(1);
	}
	gammabint = (double *)calloc(K, sizeof(double));
	if(!gammabint){
		fprintf(stderr,"Error in function sis_init: couldn't allocate memory for array gammabint\n");
		exit(1);
	}
	// Fill them in
	for(j=0; j<K; j++){
		Nt_nt[j] = Nt[j]-nt[j];
		hypercorner_bnd[j] = (nt[j] > Nt_nt[j]) ? nt[j] : Nt_nt[j];
		gammat[j] = ((double)nt[j])/Nt[j];
		gammabint[j] = gammat[j]*(1-gammat[j]);
	}
	// Some other small vectors which will be used in the function isprunable to avoid recomputing things
	f_vals = (double *)calloc(K, sizeof(double));
	if(!f_vals){
		fprintf(stderr,"Error in function sis_init: couldn't allocate memory for array f_vals\n");
		exit(1);
	}
	g_vals = (double *)calloc(K, sizeof(double));
	if(!g_vals){
		fprintf(stderr,"Error in function sis_init: couldn't allocate memory for array g_vals\n");
		exit(1);
	}
	betas = (double *)calloc(K, sizeof(double));
	if(!betas){
		fprintf(stderr,"Error in function sis_init: couldn't allocate memory for array betas\n");
		exit(1);
	}
	idx_betas_sorted = (long long *)calloc(K, sizeof(long long));
	if(!idx_betas_sorted){
		fprintf(stderr,"Error in function sis_init: couldn't allocate memory for array idx_betas_sorted\n");
		exit(1);
	}

	bitmasks = (unsigned long long *)calloc(K, sizeof(unsigned long long));
	if(!bitmasks){
		fprintf(stderr,"Error in function sis_init: couldn't allocate memory for array bitmasks\n");
		exit(1);
	}
	bitmasks[0] = 1;
	for(j=1; j<K; j++){
		bitmasks[j] = 2*bitmasks[j-1];
	}

	for(j=0; j<K; j++){
		printf("%d,%d,%d\n",x[j],nt[j],Nt[j]);
	}

	printf("%e\n",isprunable(x));
	printf("%e\n",isprunable2(x));

	free(Nt); free(nt); free(Nt_nt); free(hypercorner_bnd); free(gammat); free(gammabint);
	free(f_vals); free(g_vals); free(betas); free(idx_betas_sorted); free(x); free(bitmasks);


}
