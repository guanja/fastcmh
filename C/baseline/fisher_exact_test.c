#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include<string.h>

long long n;
long long N;

// Array with all values of log(n!) in [0,N] pre-computed
double *loggamma;
// Logarithm of 1/binom(N,n). This terms appears for every evaluation of the hypergeometric
// PDF, so it makes sense to precompute it
double log_inv_binom_N_n;

/* Precompute values of log(x!) storing them in the array loggamma */
void loggamma_init(){
	long long x;
	// Allocate memory for log-gamma cache, raising error if it fails
	loggamma = (double *)malloc((N+1)*sizeof(double));
	if(!loggamma){
		fprintf(stderr,"Error in function loggamma_init: couldn't allocate memory for array loggamma\n");
		exit(1);
	}
	// Initialise cache with appropriate values
	for(x=0;x<=N;x++) loggamma[x] = lgamma(x+1);//Gamma(x) = (x-1)!
	// Initialise log_inv_binom_N_n
	log_inv_binom_N_n = loggamma[n] + loggamma[N-n] - loggamma[N];
}
#if 0
double hypergeom_pdf(long long a, long long x){
	return exp(log_inv_binom_N_n + loggamma[x] + loggamma[N-x] - (loggamma[a] + loggamma[n-a] + loggamma[x-a] + loggamma[(N-n)-(x-a)]));
}

double fisher_pval2(long long a, long long x){
	long long a_min, a_max, k;
	double p_left, p_right;
	a_min = ((n+x-N) > 0) ? (n+x-N) : 0;//max(0,n+x-N)
	a_max = (x < n) ? x : n;//min(x,n)
	for(k=a_min,p_left=0; k <= a; k++) p_left += hypergeom_pdf(k,x);
	for(k=a,p_right=0; k <= a_max; k++) p_right += hypergeom_pdf(k,x);
	return (p_left < p_right) ? p_left : p_right;//min(p_left,p_right)
}
#endif
/* Evaluate Fisher's exact test on a table with margins x, n and N and cell count a. Note that n and N are defined as global variables.
 * The p-value is defined as a two-tailed p-value which adds up the probabilities of all tables less or equally likely to occur than the
 * one we observed
 */
double fisher_pval(long long a, long long x){
	long long a_min, a_max, k;
	double p_left, p_right, pval;
	double pre_comp_xterms;

	// Compute the contribution of all terms depending on x but not on a
	pre_comp_xterms = loggamma[x] + loggamma[N-x];
	a_min = ((n+x-N) > 0) ? (n+x-N) : 0;//max(0,n+x-N)
	a_max = (x > n) ? n : x;//min(x,n)


	// The algorithm to compute the p-value proceeds as follows. We inspect simultaneously probability values on the left and right tails of the
	// hypergeometric distribution, "accepting" each time the one that is smaller. When that happens, we move the index in the appropriate direction,
	// that is, a_min++ if we "accept" a value on the left and a_max-- if we "accept" a value on the right. As soon as the "accepted" value is located
	// in index a, we know that we have already explored all values of the hypergeometric probability mass whose probabilities are smaller or equal
	// than the probability of a. The only tricky case occurs when the probability values on both sides happen to be identical. The way to handle
	// that case is by "accepting" both values simultaneously.
	pval = 0; //Accumulate probabilities in this variable
	while(a_min<a_max){
		p_left = exp(pre_comp_xterms + log_inv_binom_N_n - (loggamma[a_min] + loggamma[n-a_min] + loggamma[x-a_min] + loggamma[(N-n)-(x-a_min)]));
		p_right = exp(pre_comp_xterms + log_inv_binom_N_n - (loggamma[a_max] + loggamma[n-a_max] + loggamma[x-a_max] + loggamma[(N-n)-(x-a_max)]));
		if(p_left == p_right) {
			pval += (p_left+p_right);
			if((a==a_min) || (a==a_max)) return pval;
			a_min++; a_max--;
		}
		else if(p_left < p_right){
			pval += p_left;
			if(a==a_min) return pval;
			a_min++;
		}
		else{
			pval += p_right;
			if(a==a_max) return pval;
			a_max--;
		}
	}
	// If we get to this part of the code, it means is the mode of the distribution and, therefore, its associated p-value is 1
	return 1;
}


int main(int argc, char *argv[]){

	long long a,x;

	a = atoll(argv[1]); x = atoll(argv[2]); n = atoll(argv[3]); N = atoll(argv[4]);

	loggamma_init();

	//printf("%e\n",fisher_pval2(a,x));
	printf("%e\n",fisher_pval(a,x));

	free(loggamma);
}
