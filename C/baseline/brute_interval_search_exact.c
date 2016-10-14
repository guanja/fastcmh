#ifndef _brute_interval_search_exact_c_
#define _brute_interval_search_exact_c_

/* LIBRARY INCLUDES */
#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include<string.h>

/* CODE DEPENDENCIES */
#include"time_keeping.c"

/* MACROS */
#define ISPRUNABLE(X) ((X) > su2)
#define ISTESTABLE(X) (psi[(X)] <= delta)


/* CONSTANT DEFINES */
#define READ_BUF_SIZ 524288 //Size of buffer used to read files

/* -------------------------------------------- GLOBAL VARIABLES -----------------------------------------------------------*/

FILE *pval_out_file; // File to output the P-values of testable intervals
FILE *sig_int_file; // File to output the P-values of significant intervals
FILE *timing_file; // File to output information about the runtime of the algorithm
FILE *summary_file; // File to output varied information about the process (number of interval processed and so on)

// Number of observations, N; and midpoint of interval [0,N], floor(N/2)
long long N, N_over_2;
// Number of observations in positive class
long long n;
// Sequence length
long long L;
long long L_max;

// Current interval length
long long l;
// Number of testable intervals
long long m;
// Target FWER
double alpha;

// Region thresholds: Sigma_k = [sl1,sl2] U [su1,su2]
// We always have su1=N-sl2 and su2=N-sl1, but keeping each variable
// separate reduces the number of computations at the expense of a tiny
//vamount of memory
long long sl1, sl2, su1, su2;
// Current P-value threshold
double delta;
// Final corrected significance threshold
double delta_opt;
// Flag variable to keep track of the last change done to region Sigma_k
// If flag==1, the last change was sl1++ (shrink on extremes of the W)
// If flag==0, the last change was sl2-- (shrink on center of the W)
int flag;

// Array with all values of log(n!) in [0,N] pre-computed
double *loggamma;
// Logarithm of 1/binom(N,n). This terms appears for every evaluation of the hypergeometric
// PDF, so it makes sense to precompute it
double log_inv_binom_N_n;
// Array with all values of minimum attainable P-value in [0,N] pre-computed
double *psi;

// Vector of class labels
char *Y_tr;
// The original dataset matrix stored as a LxN matrix, where L is the sequence length
// and N the number of observed sequences. This exploits the row-major storage of C matrices.
char **X_tr;
// Another LxN matrix, storing the values of nodes in the layer above.
char **X_par;

// A L-dimensional vector storing the frequency of each interval as they are processed
long long *freq_par;
// A (N+1)-dimensional vector such that freq_cnt[j] = #intervals with x_{i}=j processed so far
long long *freq_cnt;

/* PROFILING VARIABLES */
long long n_intervals_processed;
long long n_pvalues_computed;
long long n_significant_intervals;

/* FUNCTION DECLARATIONS */
void loggamma_init();
void psi_init();
void get_N_n(char *);
void read_labels_file(char *);
void get_L(char *);
void read_dataset_file(char *, char *);

/* -------------------------------- INITIALISATION AND TERMINATION FUNCTIONS ----------------------------------------- */

/* Initialise the Westfall-Young permutation code
 * Input arguments are self-explanatory
 * */
void sis_init(char *X_filename, char *Y_filename, double target_fwer, long long l_max){
	int j; //Loop variable
	double tic,toc;//Internal ones, do not overwrite the ones in time_keeping.c

	// Compute total number of observations and number of observations in minority class
	tic = measureTime();
	get_N_n(Y_filename);
	get_L(X_filename);
	toc = measureTime();
	time_IO += toc-tic;

	// Store core constants
	N_over_2 = (N % 2) ? (N-1)/2 : N/2;//floor(N/2)
	alpha = target_fwer;
	L_max = l_max;
	// And initialise some others
	sl1 = 1; sl2 = N_over_2; su1 = N-sl2; su2 = N-sl1;
	flag = 1;
	delta = ((double) n)/N; //$\psi(1)=\frac{n}{N}$

	// Initialise cache for log(x!) and psi(x)
	loggamma_init(); psi_init();

	// Allocate space for class labels
	Y_tr = (char *)malloc(N*sizeof(char));
	if(!Y_tr){
		fprintf(stderr,"Error in function sis_init: couldn't allocate memory for array Y_tr\n");
		exit(1);
	}
	// And store them in memory from file
	tic = measureTime();
	read_labels_file(Y_filename);
	toc = measureTime();
	time_IO += toc-tic;

	// Ensure class 1 is the minority class
	if(n > (N/2)){
		for(j=0; j<N; j++) Y_tr[j] = !Y_tr[j];
		n = N-n;
	}

	// Initialise dataset matrix

	X_tr = (char **)malloc(L*sizeof(char *));
	if(!X_tr){
		fprintf(stderr,"Error in function sis_init: couldn't allocate memory for array X_tr\n");
		exit(1);
	}

	X_tr[0] = (char *)calloc(L*N, sizeof(char));
	if(!X_tr[0]){
		fprintf(stderr,"Error in function sis_init: couldn't allocate memory for array X_tr[0]\n");
		exit(1);
	}
	for(j=1; j<L; j++) X_tr[j] = X_tr[0] + j*N;
	// Same for parents
	X_par = (char **)malloc(L*sizeof(char *));
	if(!X_par){
		fprintf(stderr,"Error in function sis_init: couldn't allocate memory for array X_par\n");
		exit(1);
	}
	X_par[0] = (char *)calloc(L*N, sizeof(char));
	if(!X_par[0]){
		fprintf(stderr,"Error in function sis_init: couldn't allocate memory for array X_par[0]\n");
		exit(1);
	}
	for(j=1; j<L; j++) X_par[j] = X_par[0] + j*N;

	tic = measureTime();
	read_dataset_file(X_filename, X_tr[0]);
	toc = measureTime();
	time_IO += toc-tic;


	freq_par = (long long *)calloc(L, sizeof(long long));
	if(!freq_par){
		fprintf(stderr,"Error in function sis_init: couldn't allocate memory for array freq_par\n");
		exit(1);
	}
	freq_cnt = (long long *)calloc(N+1, sizeof(long long));
	if(!freq_cnt){
		fprintf(stderr,"Error in function sis_init: couldn't allocate memory for array freq_cnt\n");
		exit(1);
	}

}

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

/* Precompute minimum attainable P-values $\psi(x)$ for all x in [0,N] and store them in array psi */
void psi_init(){
	double xi1;
	long long x, x_init;
	// Allocate memory for psi, raising error if it fails
	psi = (double *)malloc((N+1)*sizeof(double));
	if(!psi){
		fprintf(stderr,"Error in function psi_and_xi1_init: couldn't allocate memory for array psi\n");
		exit(1);
	}

	/* Initialise caches with appropriate values */

	// First compute the left side of "the W", i.e. the range [0,n]
	psi[0] = 1;
	//In this range, the recursion $\psi(x)$=$\psi(x-1)$*[(n-(x-1))/(N-(x-1))] can be seen to be correct
	for(x=1; x<=n; x++) psi[x] = (((double)(n-(x-1)))/(N-(x-1)))*psi[x-1];

	// Now, start computing xi1 in the range [N-N_over_2,N] using another recursion, this time
	// starting in N
	// Note that we don't need to store all values, since this will be used only to initialise
	// psi[N_over_2]
	x_init = N-N_over_2;
	xi1 = 1;
	//In this range, the recursion $\xi_{1}(x)$=$\xi_{1}(x+1)$*[((x-1)-n)/(x+1)] can be seen to be correct
	for(x=(N-1); x>=x_init; x--) xi1 = (((double)((x+1)-n))/(x+1))*xi1;

	// Now, use the value of $\xi_{1}(N-N_over_2)$=xi1[0] to get $\psi(N_over_2)$=psi[N_over_2] using the
	// same recursion if N is odd, or simply copy the value of xi1[0] since $\xi_{1}(N-N_over_2)=\psi(N_over_2)$
	// if N is even
	if (N % 2) psi[N_over_2] = (((double)(x_init-n))/x_init)*xi1;
	else psi[N_over_2] = xi1;

	// Now, using psi[N_over_2] compute the right side of "the W", i.e. the range [n+1,N_over_2]
	// using the same recursion as for $\xi_{1}$
	for(x=(N_over_2-1); x > n; x--) psi[x] = (((double)((x+1)-n))/(x+1))*psi[x+1];

	// Finally, since $\psi(x)$ is symmetric around N_over_2, complete the right half by symmetry
	for(x=x_init; x<=N; x++) psi[x] = psi[N-x];

	// Correct minimum attainable P-value in some edge-cases
	if((N % 2)==0){
		if (n == (N/2)) for(x=1; x<N; x++) psi[x] *= 2;
		else psi[N/2] *= 2;
	}
}

/* Free all allocated memory and give some output for debugging purposes */
void sis_end(){
	// OUTPUT RESULTS
	// Execution time and peak memory consumption
	profileCode();

	// Free allocated memory
	free(loggamma); free(psi);
	free(Y_tr);
	free(X_tr[0]); free(X_par[0]);
	free(X_tr); free(X_par);
	free(freq_par); free(freq_cnt);

	// Close files
	if(pval_out_file) fclose(pval_out_file);
	fclose(sig_int_file); fclose(timing_file); fclose(summary_file);
}

/* ---------------------------------------FUNCTIONS TO FIND THE SIGNIFICANT INTERVALS-------------------------------- */

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

void process_first_layer_pvalues(){
	long long tau, j, a;
	char *X_tr_aux;
	double pval;
	// Process each length 1 interval
	for(tau=0; tau<L; tau++){
		// Compute number of 1s in the interval
		X_tr_aux = X_tr[tau];
		for(j=0; j<N; j++) freq_par[tau] += X_tr_aux[j];
		// Brute-force approach considers all intervals testable
		// Compute cell count
		#ifndef NO_SINGLE_FEATURES
		a = 0;
		for(j=0; j<N; j++) a += X_tr_aux[j]&&Y_tr[j];
		// Compute the p-value
		pval = fisher_pval(a,freq_par[tau]); n_pvalues_computed++;
		// Check if the P-value is significant
		if(pval_out_file) fprintf(pval_out_file,"%lld,%lld,%lld,%lld,%e\n",l+1,tau,a,freq_par[tau],pval);
		if(pval <= delta_opt) { fprintf(sig_int_file,"%lld,%lld,%lld,%lld,%e\n",l+1,tau,a,freq_par[tau],pval); n_significant_intervals++; }
		#endif
	}
}

void process_intervals_pvalues(){
	long long tau, j, a;
	char *X_tr_aux, *X_par_aux;
	double pval;
	for(l=1; l<L_max; l++){
		//printf("\tProcessing layer %lld...\n",l+1);
		for(tau=0; tau<(L-l); tau++){
			// Compute OR and frequency of the interval
			X_tr_aux = X_tr[tau+l]; X_par_aux = X_par[tau];
			for(j=0; j<N; j++) if((!X_par_aux[j]) && X_tr_aux[j]){ X_par_aux[j] = 1; freq_par[tau]++;}
			// Brute-force approach considers all intervals testable
			// Compute cell count
			a = 0;
			for(j=0; j<N; j++) a += X_par_aux[j]&&Y_tr[j];
			// Compute the P-value
			pval = fisher_pval(a,freq_par[tau]); n_pvalues_computed++;
			// Check if the P-value is significant
			if(pval_out_file) fprintf(pval_out_file,"%lld,%lld,%lld,%lld,%e\n",l+1,tau,a,freq_par[tau],pval);
			if(pval <= delta_opt) { fprintf(sig_int_file,"%lld,%lld,%lld,%lld,%e\n",l+1,tau,a,freq_par[tau],pval); n_significant_intervals++; }
		}
	}
	l = l - 1;
}

/* Wrapper function that encapsulates the functionality required to find significant intervals */
void find_significant_intervals(){
	// Give feedback to user
	#ifndef NO_VERBOSE
	printf("\n\nSCANNING DATASET FOR SIGNIFICANT INTERVALS...\n\n");
	#endif
	// Initialise current layer index and current number of computed p-values to 0
	l = 0; n_pvalues_computed = 0; n_significant_intervals = 0;
	// Clear the current layer frequency counters
	memset(freq_par,0,L*sizeof(long long));
	// Initialise the value of the OR vectors of current layer to original dataset
	memcpy(X_par[0],X_tr[0],L*N*sizeof(char));
	// Process the upper-most layer (i.e. the layer composed of length 1 intervals)
	process_first_layer_pvalues();
	// Process the rest of layers (i.e. intervals of length > 1) until the pruning naturally stops the execution
	process_intervals_pvalues();
	// Report number of significant intervals found
	fprintf(summary_file,"Number of significantly associated intervals found: %lld\n",n_significant_intervals);
}

/* -------------------FUNCTIONS TO FIND THE CORRECTED SIGNIFICANCE THRESHOLD -------------------------------------- */

void compute_significance_threshold(){
	if(L_max==0) L_max = L;
	for(l=0; l<L_max; l++) m += L-l;
	#ifdef NO_SINGLE_FEATURES
	m -= L; //If the intervals of length 1 are not to be tested, readjust the correction factor
	#endif
}

/* Function to give some feedback about the computation of the significance threshold */
void output_significance_threshold(){
	fprintf(summary_file,"DATASET CHARACTERISTICS:\n");
	fprintf(summary_file,"\tN = %lld, n = %lld, L = %lld\n",N,n,L);
	fprintf(summary_file,"RESULTS: \n");
	// Number of intervals processed, proportion of intervals pruned
	fprintf(summary_file,"Intervals processed: %lld (%f%% of total)\n",n_intervals_processed,((double)(200*n_intervals_processed))/(L*(L+1)));
	fprintf(summary_file,"Maximum testable interval length: %lld\n",l);
	if(L_max==0) fprintf(summary_file,"Maximum interval length to be processed: unlimited\n");
	else fprintf(summary_file,"Maximum interval length to be processed: %lld\n",L_max);
	fprintf(summary_file,"Number of testable intervals: %lld\n",m);
	fprintf(summary_file,"Corrected significance threshold at level %e: %e\n",alpha, delta_opt);

}

/* Wrapper function that encapsulates the functionality required to find the corrected significance threshold */
void compute_corrected_significance_threshold(){
	// Give feedback to user
	#ifndef NO_VERBOSE
	printf("COMPUTING CORRECTED SIGNIFICANCE THRESHOLD...\n");
	#endif
	// Initialise current layer index, current number of testable intervals
	l = 0; m = 0;
	// Compute number of intervals to be processed
	compute_significance_threshold();
	n_intervals_processed = m;
	// Set final corrected significance threshold
	delta_opt = alpha/m;
	// Print results to stdout
	output_significance_threshold();
}

/*--------------------------------------------------------- FILE I/O --------------------------------------------------------*/
/* Do a first scan of the file containing the class labels to compute the total number of observations, N,
 * and the total number of observations in the positive class, n
 * */
void get_N_n(char *labels_file){
	FILE *f_labels;//Stream with file containing class labels
	int n_read;//Number of chars read
	int i;// Iterator variable to be used in loops
	char char_to_int[256];//Array for converting chars to int fast
	char *read_buf, *read_buf_aux, *read_buf_end;//Buffer for reading from file and extra pointers for loops

	// Initialise both counters to 0 (the variables are defined as global variables in wy.c)
	N = 0; n = 0;

	//Try to open file, giving an error message if it fails
	if(!(f_labels = fopen(labels_file,"r"))){
		fprintf(stderr, "Error in function get_N_n when opening file %s\n",labels_file);
		exit(1);
	}

	//Try to allocate memory for the buffer, giving an error message if it fails
	read_buf = (char *)malloc(READ_BUF_SIZ*sizeof(char));
	if(!read_buf){
		fprintf(stderr,"Error in function get_N_n: couldn't allocate memory for array read_buf\n");
		exit(1);
	}

	//Initialize the char to int converter
	for(i=0;i<256;i++) char_to_int[i] = 127;
	// We only care about the chars '0' and '1'. Everything else is mapped into the same "bucket"
	char_to_int['0'] = 0; char_to_int['1'] = 1;

	// Read the entire file
	while(1){
		// Try to read READ_BUF_SIZ chars from the file containing the class labels
		n_read = fread(read_buf,sizeof(char),READ_BUF_SIZ,f_labels);
		// If the number of chars read, n_read_ is smaller than READ_BUF_SIZ, either the file ended
		// or there was an error. Check if it was the latter
		if((n_read < READ_BUF_SIZ) && !feof(f_labels)){
			fprintf(stderr,"Error in function get_N_n while reading the file %s\n",labels_file);
			exit(1);
		}
		// Process the n_read chars read from the file
		for(read_buf_aux=read_buf,read_buf_end=read_buf+n_read;read_buf_aux<read_buf_end;read_buf_aux++){
			//If the character is anything other than '0' or '1' go to process the next char
			if(char_to_int[*read_buf_aux] == 127) continue;
			N++;
			if(char_to_int[*read_buf_aux]) n++;
		}
		// Check if the file ended,. If yes, then exit the while loop
		if(feof(f_labels)) break;
	}

	//Close the file
	fclose(f_labels);

	//Free allocated memory
	free(read_buf);
}

void read_labels_file(char *labels_file){
	FILE *f_labels;//Stream with file containing class labels
	int n_read;//Number of chars read
	int i;// Iterator variable to be used in loops
	char char_to_int[256];//Array for converting chars to int fast
	char *read_buf, *read_buf_aux, *read_buf_end;//Buffer for reading from file and extra pointers for loops
	char *labels_aux = Y_tr;//Auxiliary pointer to array labels for increments

	//Try to open file, giving an error message if it fails
	if(!(f_labels = fopen(labels_file,"r"))){
		fprintf(stderr, "Error in function read_labels_file when opening file %s\n",labels_file);
		exit(1);
	}

	//Try to allocate memory for the buffer, giving an error message if it fails
	read_buf = (char *)malloc(READ_BUF_SIZ*sizeof(char));
	if(!read_buf){
		fprintf(stderr,"Error in function read_labels_file: couldn't allocate memory for array read_buf\n");
		exit(1);
	}

	//Initialize the char to int converter
	for(i=0;i<256;i++) char_to_int[i] = 127;
	// We only care about the chars '0' and '1'. Everything else is mapped into the same "bucket"
	char_to_int['0'] = 0; char_to_int['1'] = 1;

	// Read the entire file
	while(1){
		// Try to read READ_BUF_SIZ chars from the file containing the class labels
		n_read = fread(read_buf,sizeof(char),READ_BUF_SIZ,f_labels);
		// If the number of chars read, n_read_ is smaller than READ_BUF_SIZ, either the file ended
		// or there was an error. Check if it was the latter
		if((n_read < READ_BUF_SIZ) && !feof(f_labels)){
			fprintf(stderr,"Error in function read_labels_file while reading the file %s\n",labels_file);
			exit(1);
		}
		// Process the n_read chars read from the file
		for(read_buf_aux=read_buf,read_buf_end=read_buf+n_read;read_buf_aux<read_buf_end;read_buf_aux++){
			//If the character is anything other than '0' or '1' go to process the next char
			if(char_to_int[*read_buf_aux] == 127) continue;
			*labels_aux++ = char_to_int[*read_buf_aux];
		}
		// Check if the file ended,. If yes, then exit the while loop
		if(feof(f_labels)) break;
	}

	// Sanity check to see if we successfully read the correct number of labels
	i = labels_aux-Y_tr;
	if(i != N){
		fprintf(stderr,"Error in function read_labels_file: incorrect number of labels read. Read %d, correct number %lld\n",i,N);
		exit(1);
	}

	//Close the file
	fclose(f_labels);

	//Free allocated memory
	free(read_buf);
}

void get_L(char *filename){
	FILE *f_dat = ((FILE*)0);
	int i, j, n_read;
	char char_to_int[256];//Array for converting chars to int fast
	char *read_buf, *read_buf_aux, *read_buf_end;//Buffer for reading from file and extra pointers for loops

	//Try to open file, giving an error message if it fails
	if(!(f_dat = fopen(filename,"r"))){
		fprintf(stderr, "Error in function get_L when opening file %s\n",filename);
		exit(1);
	}
	//Try to allocate memory for the buffer, giving an error message if it fails
	read_buf = (char *)malloc(READ_BUF_SIZ*sizeof(char));
	if(!read_buf){
		fprintf(stderr,"Error in function get_L: couldn't allocate memory for array read_buf\n");
		exit(1);
	}

	//Initialize the char to int converter
	for(i=0;i<256;i++) char_to_int[i] = 0;
	// We only care about the chars '0' and '1'. Everything else is mapped into the same "bucket"
	char_to_int['\n'] = 1;

	// Read the entire file, counting the number of lines
	L = 0;
	while(1){
		// Try to read READ_BUF_SIZ chars from the file containing the class labels
		n_read = fread(read_buf,sizeof(char),READ_BUF_SIZ,f_dat);
		// If the number of chars read, n_read_ is smaller than READ_BUF_SIZ, either the file ended
		// or there was an error. Check if it was the latter
		if((n_read < READ_BUF_SIZ) && !feof(f_dat)){
			fprintf(stderr,"Error in function get_L while reading the file %s\n",filename);
			exit(1);
		}
		// Process the n_read chars read from the file
		for(read_buf_aux=read_buf,read_buf_end=read_buf+n_read;read_buf_aux<read_buf_end;read_buf_aux++) if(char_to_int[*read_buf_aux]) L++;
		// Check if the file ended,. If yes, then exit the while loop
		if(feof(f_dat)) break;
	}

	// Close file
	fclose(f_dat);
	// Free allocated memory
	free(read_buf);
}


void read_dataset_file(char *filename, char *ptr){
	FILE *f_dat = ((FILE*)0);
	int i, j, n_read;
	char char_to_int[256];//Array for converting chars to int fast
	char *read_buf, *read_buf_aux, *read_buf_end;//Buffer for reading from file and extra pointers for loops

	//Try to open file, giving an error message if it fails
	if(!(f_dat = fopen(filename,"r"))){
		fprintf(stderr, "Error in function read_dataset_file when opening file %s\n",filename);
		exit(1);
	}
	//Try to allocate memory for the buffer, giving an error message if it fails
	read_buf = (char *)malloc(READ_BUF_SIZ*sizeof(char));
	if(!read_buf){
		fprintf(stderr,"Error in function read_dataset_file: couldn't allocate memory for array read_buf\n");
		exit(1);
	}

	//Initialize the char to int converter
	for(i=0;i<256;i++) char_to_int[i] = 127;
	// We only care about the chars '0' and '1'. Everything else is mapped into the same "bucket"
	char_to_int['0'] = 0; char_to_int['1'] = 1;

	// Read the entire file
	while(1){
		// Try to read READ_BUF_SIZ chars from the file containing the class labels
		n_read = fread(read_buf,sizeof(char),READ_BUF_SIZ,f_dat);
		// If the number of chars read, n_read_ is smaller than READ_BUF_SIZ, either the file ended
		// or there was an error. Check if it was the latter
		if((n_read < READ_BUF_SIZ) && !feof(f_dat)){
			fprintf(stderr,"Error in function read_dataset_file while reading the file %s\n",filename);
			exit(1);
		}
		// Process the n_read chars read from the file
		for(read_buf_aux=read_buf,read_buf_end=read_buf+n_read;read_buf_aux<read_buf_end;read_buf_aux++){
			//If the character is anything other than '0' or '1' go to process the next char
			if(char_to_int[*read_buf_aux] == 127) continue;
			*ptr++ = char_to_int[*read_buf_aux];
		}
		// Check if the file ended,. If yes, then exit the while loop
		if(feof(f_dat)) break;
	}

	// Close file
	fclose(f_dat);
	// Free allocated memory
	free(read_buf);
}

/* ----------------------------------------------------ENTRY POINT---------------------------------------------------------- */

int main(int argc, char *argv[]){
	char *tmp_filename;
	char *R_command;
    char i;
    char n_fixed_args = 0;
    char idx_fixed_args[5] = {-1,-1,-1,-1,-1};
    char idx_postprocessing_folder = -1, idx_out_pvals_file = -1;

    // Process input
    for(i=1; i<argc; i++){
        if((strcmp(argv[i],"-postprocessing_folder")==0) || (strcmp(argv[i],"-pp")==0)) idx_postprocessing_folder = ++i;
        else if((strcmp(argv[i],"-pval_file")==0)) idx_out_pvals_file = ++i;
        else idx_fixed_args[n_fixed_args++] = i;
    }

	// Check input
	if(n_fixed_args != 5){
		printf("ERROR: INCORRECT SYNTAX!\n");
		printf("\tUSAGE: ./program_name X_file Y_file alpha L_max base_filename [-postprocessing_folder path_to_pyfiltering.py] [-pval_file all_pvals_file]\n");
    }

	// Get time when program started to run
	t_init = measureTime();

	// INITIALISATION
	tic = measureTime();

	// First allocate memory for filename holder
	tmp_filename = (char *)malloc((strlen(argv[idx_fixed_args[4]])+512)*sizeof(char));

	// Create a file to report significant intervals
	strcpy(tmp_filename,argv[idx_fixed_args[4]]); strcat(tmp_filename,"_sigints.csv");
	if(!(sig_int_file = fopen(tmp_filename,"w"))){
		fprintf(stderr, "Error in function main when opening file %s\n",tmp_filename);
		exit(1);
	}
	// If the file was successfully create, write the file header
	fprintf(sig_int_file,"l,tau,a,x,P-value\n");

	// Create a file to report runtime information
	strcpy(tmp_filename,argv[idx_fixed_args[4]]); strcat(tmp_filename,"_timing.txt");
	if(!(timing_file = fopen(tmp_filename,"w"))){
		fprintf(stderr, "Error in function main when opening file %s\n",tmp_filename);
		exit(1);
	}

	// Create a file to report varied information such as number of intervals processed, etc.
	strcpy(tmp_filename,argv[idx_fixed_args[4]]); strcat(tmp_filename,"_summary.txt");
	if(!(summary_file = fopen(tmp_filename,"w"))){
		fprintf(stderr, "Error in function main when opening file %s\n",tmp_filename);
		exit(1);
	}

	// Free filename holder
	free(tmp_filename);

	// If the optional argument was used, create a file to output all testable P-values
	pval_out_file = ((FILE *)0);
	if(idx_out_pvals_file != -1){
		// Actually try to create file, if that fails raise an error
		if(!(pval_out_file = fopen(argv[idx_out_pvals_file],"w"))){
			fprintf(stderr, "Error in function main when opening file %s\n",argv[idx_out_pvals_file]);
			exit(1);
		}
		// If the file was successfully create, write the file header
		fprintf(pval_out_file,"l,tau,a,x,P-value\n");
	}

	// Initialise main variables
	sis_init(argv[idx_fixed_args[0]], argv[idx_fixed_args[1]], atof(argv[idx_fixed_args[2]]), atoll(argv[idx_fixed_args[3]]));

	toc = measureTime();
	time_initialisation = toc-tic;

	// MAIN FUNCTIONALITY
	tic = measureTime();
	compute_corrected_significance_threshold();
	toc = measureTime();
	time_comp_threshold = toc-tic;
	tic = measureTime();
	find_significant_intervals();
	toc = measureTime();
	time_comp_significant_intervals = toc-tic;

	// Get time when program finished
	t_end = measureTime();

	// Produce output and free memory
	sis_end();

	// Call R postprocessing script
	R_command = (char *)malloc((strlen(argv[idx_fixed_args[4]])+512)*sizeof(char));
	strcpy(R_command,"python "); 
    if(idx_postprocessing_folder != -1) strcat(R_command,argv[idx_postprocessing_folder]);
    strcat(R_command,"pyfiltering.py "); strcat(R_command,argv[idx_fixed_args[4]]);
	strcat(R_command,"_sigints.csv "); strcat(R_command,argv[idx_fixed_args[4]]); strcat(R_command,"_sigints_filtered.corrected.csv");
	system(R_command);
	free(R_command);

	// Return
	exit(0);
}

#endif
