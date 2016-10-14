#ifndef _significant_interval_search_wy_chi_c_
#define _significant_interval_search_wy_chi_c_

/* LIBRARY INCLUDES */
#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include<string.h>

/* CODE DEPENDENCIES */
#include"time_keeping.c"
#include "../../EasyGWAS/chi2.h"

/* MACROS */
#define ISPRUNABLE(X) ((X) > su2)
#define ISTESTABLE(X) (psi[(X)] <= delta)


/* CONSTANT DEFINES */
#define READ_BUF_SIZ 524288 //Size of buffer used to read files
#define STORE_PERMUTATIONS_FILE 0 //Set the flag to 1 if the generated permutions should be written to a file for reuse/debugging

/* -------------------------------------------- GLOBAL VARIABLES -----------------------------------------------------------*/

FILE *pval_out_file; // File to output the P-values of testable intervals
FILE *sig_int_file; // File to output the P-values of significant intervals
FILE *timing_file; // File to output information about the runtime of the algorithm
FILE *summary_file; // File to output varied information about the process (number of interval processed and so on)
FILE *minpvals_file; // File to output the minimum p-values per permutation

// Number of observations, N; and midpoint of interval [0,N], floor(N/2)
long long N, N_over_2;
// Number of observations in positive class
long long n;
// Sequence length
long long L;
long long L_max;

// Current interval length
long long l;
// Current FWER
double FWER;
// FWER at corrected significance threshold
double FWER_opt;
// Target FWER
double alpha;

// Region thresholds: Sigma_k = [sl1,sl2] U [su1,su2]
// We always have su1=N-sl2 and su2=N-sl1, but keeping each variable
// separate reduces the number of computations at the expense of a tiny
// amount of memory
long long sl1, sl2, su1, su2;
// Current P-value threshold
double delta;
// Final corrected significance threshold
double delta_opt;
// Flag variable to keep track of the last change done to region Sigma_k
// If flag==1, the last change was sl1++ (shrink on extremes of the W)
// If flag==0, the last change was sl2-- (shrink on center of the W)
int flag;

// Array with all values of minimum attainable P-value in [0,N] pre-computed
double *psi;
// And some constants which are useful to precompute
double class_ratio, class_ratio_bin;

// Vector of class labels
char *Y_tr;
// The original dataset matrix stored as a LxN matrix, where L is the sequence length
// and N the number of observed sequences. This exploits the row-major storage of C matrices.
char **X_tr;
// Another LxN matrix, storing the values of nodes in the layer above.
char **X_par;

// Variables for Westfall-Young permutation
// Number of permutations
long long J;
//Random generator seed for reproducibility
time_t seed;
// Matrix of size J x (# non-empty transactions)
char **Y_tr_perm;

// A L-dimensional vector storing the frequency of each interval as they are processed
long long *freq_par;

// Queue of testable intervals in the layer below
long long *testable_queue;
long long testable_queue_front;
long long testable_queue_length;

// Auxiliary variable to keep track of current layer
long long last_tau;

// Minimum P-value for each permutation
double *min_pval;

// Cell-count counter (table entry (X=1,Y=1))
// While finding significance threshold J-dimensional vector (one cell count per permutation)
long long *a_cnt;

/* PROFILING VARIABLES */
long long n_intervals_processed;
long long n_pvalues_computed;
long long n_significant_intervals;

/* FUNCTION DECLARATIONS */
void psi_init();
void get_N_n(char *);
void read_labels_file(char *);
void get_L(char *);
void read_dataset_file(char *, char *);
void randperm(char *, char *, int, FILE *);
int doublecomp(const void*,const void*);
extern double Chi2_sf(double,double); //Chi2 Survival function, from Dominik's EasyGWASCore

/* -------------------------------- INITIALISATION AND TERMINATION FUNCTIONS ----------------------------------------- */

/* Initialise the Westfall-Young permutation code
 * Input arguments are self-explanatory
 * */
void siswy_init(char *X_filename, char *Y_filename, double target_fwer, long long n_perm, long long l_max){
	long long i,j;// Iterator variables to be used in loops
	FILE *f_perm = ((FILE*)0);
	char *perm_buffer;
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
	delta = psi[1];

	// Allocate space for class labels
	Y_tr = (char *)malloc(N*sizeof(char));
	if(!Y_tr){
		fprintf(stderr,"Error in function siswy_init: couldn't allocate memory for array Y_tr\n");
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

	// Initialise cache for psi(x)
	psi_init();
	// Initialise corrected significance threshold
	delta = psi[sl1];

	// Initialise dataset matrix
	X_tr = (char **)malloc(L*sizeof(char *));
	if(!X_tr){
		fprintf(stderr,"Error in function siswy_init: couldn't allocate memory for array X_tr\n");
		exit(1);
	}

	X_tr[0] = (char *)calloc(L*N, sizeof(char));
	if(!X_tr[0]){
		fprintf(stderr,"Error in function siswy_init: couldn't allocate memory for array X_tr[0]\n");
		exit(1);
	}
	for(j=1; j<L; j++) X_tr[j] = X_tr[0] + j*N;
	// Same for parents
	X_par = (char **)malloc(L*sizeof(char *));
	if(!X_par){
		fprintf(stderr,"Error in function siswy_init: couldn't allocate memory for array X_par\n");
		exit(1);
	}
	X_par[0] = (char *)calloc(L*N, sizeof(char));
	if(!X_par[0]){
		fprintf(stderr,"Error in function siswy_init: couldn't allocate memory for array X_par[0]\n");
		exit(1);
	}
	for(j=1; j<L; j++) X_par[j] = X_par[0] + j*N;

	tic = measureTime();
	read_dataset_file(X_filename, X_tr[0]);
	toc = measureTime();
	time_IO += toc-tic;

	// Allocate memory for several vectors
	testable_queue = (long long *)calloc(L, sizeof(long long));
	if(!testable_queue){
		fprintf(stderr,"Error in function siswy_init: couldn't allocate memory for array testable_queue\n");
		exit(1);
	}
	freq_par = (long long *)calloc(L, sizeof(long long));
	if(!freq_par){
		fprintf(stderr,"Error in function sis_init: couldn't allocate memory for array freq_par\n");
		exit(1);
	}

	// Initialise Westfall-Young permutation-related variables
	// Save core constants as global variables
	J = n_perm;

	/* Allocate memory for the matrix of permuted labels */
	// First allocate memory for the J row pointers
	Y_tr_perm = (char **)malloc(N*sizeof(char *));
	if(!Y_tr_perm){
		fprintf(stderr,"Error in function siswy_init: couldn't allocate memory for array Y_tr_perm\n");
		exit(1);
	}
	// Now allocate memory for a contiguous block of J*(# non-empty transactions) chars
	Y_tr_perm[0] = (char *)malloc(J*N*sizeof(char));
	if(!Y_tr_perm[0]){
		fprintf(stderr,"Error in function siswy_init: couldn't allocate memory for array Y_tr_perm[0]\n");
		exit(1);
	}
	// And make each row pointer point to the appropriate position (labels_perm[0] is
	// already correctly set)
	for(i=1;i<N;i++) Y_tr_perm[i] = Y_tr_perm[0] + i*J;

	// Initialize random number generator
	seed = time(NULL);
	srand(seed);

	/* Compute the random permutations of the class labels */

	// Create file for storing the permutations (only if STORE_PERMUTATIONS_FILE==1)
	if(STORE_PERMUTATIONS_FILE){
		f_perm = fopen("permutations.csv","w");
		// Check if file was correctly created, otherwise raise an error
		if(!f_perm){
			fprintf(stderr,"Error in function siswy_init: file %s couldn't be created\n","permutations.csv");
			exit(1);
		}
	}

	// Do the J permutations themselves, storing them in labels_perm and
	// saving them to the file f_perm if necessary
	perm_buffer = (char *)malloc(N*sizeof(char));
	for(j=0;j<J;j++) {
		randperm(perm_buffer, Y_tr, N, f_perm);
		// Dump contents of buffer into destination, skipping values corresponding to empty observations/transactions
		for(i=0;i<N;i++) Y_tr_perm[i][j] = perm_buffer[i];
	}
	free(perm_buffer);

	// Close the file f_perm (if necessary)
	if(f_perm) fclose(f_perm);

	// Allocate memory for minimum p-values, raising error if it fails
	min_pval = (double *)malloc(J*sizeof(double));
	if(!min_pval){
		fprintf(stderr,"Error in function siswy_init: couldn't allocate memory for array min_pval\n");
		exit(1);
	}
	// Initialise all p-values to 1
	for(j=0; j<J; j++) min_pval[j] = 1;

	// Allocate memory for cell counts, raising an error if it fails
	a_cnt = (long long *)calloc(J, sizeof(long long));
	if(!a_cnt){
		fprintf(stderr,"Error in function siswy_init: couldn't allocate memory for array a_cnt\n");
		exit(1);
	}

}

/* Precompute minimum attainable P-values $\psi(x)$ for all x in [0,N] and store them in array psi */
void psi_init(){
	double num, den;
	long long x;

	// Allocate memory for psi, raising error if it fails
	psi = (double *)malloc((N+1)*sizeof(double));
	if(!psi){
		fprintf(stderr,"Error in function psi_init: couldn't allocate memory for array psi\n");
		exit(1);
	}

	/* Initialise caches with appropriate values */

	// Precompute some useful constants
	class_ratio = ((double)n)/N; class_ratio_bin = class_ratio*(1-class_ratio);

	// First compute the left side of "the W", i.e. the range [0,n]
	psi[0] = 1;
	for(x=1; x<=n; x++) {
		num = x*(1-class_ratio); num *= num;
		den = x*(1-((double)x)/N)*class_ratio_bin;
		psi[x] = Chi2_sf(num/den,1);
	}

	// Now, compute the minimum attainable p-values in the range [N-N_over_2,N]
	for(x=(n+1); x<=N_over_2; x++) {
		num = n*(1-((double)x)/N); num *= num;
		den = x*(1-((double)x)/N)*class_ratio_bin;
		psi[x] = Chi2_sf(num/den,1);
	}

	// Finally, since $\psi(x)$ is symmetric around N_over_2, complete the right half by symmetry
	for(x=(N_over_2+1); x<=N; x++) psi[x] = psi[N-x];
}

/* Free all allocated memory and give some output for debugging purposes */
void siswy_end(){
	//long long j;

	// Output execution time and peak memory consumption
	profileCode();

	// Free allocated memory
	free(psi);
	free(Y_tr);
	free(X_tr[0]); free(X_par[0]);
	free(X_tr); free(X_par);
	free(freq_par); free(testable_queue);
	free(Y_tr_perm[0]); free(Y_tr_perm);
	free(min_pval); free(a_cnt);

	// Close files
	if(pval_out_file) fclose(pval_out_file);
	fclose(sig_int_file); fclose(timing_file); fclose(summary_file); fclose(minpvals_file);
}

/* ---------------------------------------FUNCTIONS TO FIND THE SIGNIFICANT INTERVALS-------------------------------- */

double compute_pval(long long a, long long x){
	double aux, num, den;
	aux = ((double)x)/N;
	num = a-n*aux; num = pow(num,2);
	den = x*(1-aux)*class_ratio_bin;
	return Chi2_sf(num/den,1);
}

void process_first_layer_pvalues(){
	long long tau, j, queue_idx, a;
	char *X_tr_aux;
	double pval;
	// Process each length 1 interval
	for(tau=0; tau<L; tau++){
		// Compute number of 1s in the interval
		X_tr_aux = X_tr[tau];
		for(j=0; j<N; j++) freq_par[tau] += X_tr_aux[j];
		// If the interval is testable...
		// Update frequency-buckets and number of testable intervals
		#ifndef NO_SINGLE_FEATURES
		if(ISTESTABLE(freq_par[tau])){
			// Compute cell count
			a = 0;
			for(j=0; j<N; j++) if(X_tr_aux[j]) a += Y_tr[j];
			// Compute the p-value
			pval = compute_pval(a,freq_par[tau]); n_pvalues_computed++;
			// Check if the P-value is significant
			if(pval_out_file) fprintf(pval_out_file,"%lld,%lld,%lld,%lld,%e\n",l+1,tau,a,freq_par[tau],pval);
			if(pval <= delta_opt) { fprintf(sig_int_file,"%lld,%lld,%lld,%lld,%e\n",l+1,tau,a,freq_par[tau],pval); n_significant_intervals++; }
		}
		#endif
		// If either the current interval or the previous one are prunable (i.e. have more than su2 ones)
		// then do NOT append the left-child to the testable queue (i.e. prune it)
		if((tau==0) || ISPRUNABLE(freq_par[tau]) || ISPRUNABLE(freq_par[tau-1])) continue;
		// Compute index of array position to be used, wrapping around if necessary
		queue_idx = testable_queue_front + testable_queue_length;
		queue_idx = (queue_idx < L) ? queue_idx : queue_idx - L;
		// Actually append children to testable queue
		testable_queue[queue_idx] = tau-1;
		// Update queue length
		testable_queue_length++;
	}
}

void process_intervals_pvalues(){
	long long tau, j, queue_idx, a;
	char *X_tr_aux, *X_par_aux;
	double pval;
	// While testable-interval queue is not empty, continue to process intervals
	while(testable_queue_length){
		// Pop a testable interval from the queue
		tau = testable_queue[testable_queue_front];
		testable_queue_front = (testable_queue_front<(L-1)) ? testable_queue_front + 1: 0;
		testable_queue_length--;
		// Check if we have started processing a new layer by detecting non-monotonicity in tau
		if(tau < last_tau) {
			l++;
			#ifndef NO_VERBOSE
			printf("\tProcessing layer %lld...\n",l+1);
			#endif
		}
		if((L_max>0) && ((l+1) > L_max)) {
			#ifndef NO_VERBOSE
			printf("\tMaximum interval length achieved at l=%lld. Stopping enumeration...\n",l+1);
			#endif
			break;
		}
		last_tau = tau;
		// In this case, the testable region does not change, so we don't need to check if the interval
		// has to be pruned now. If it was appended to the queue, then it has to be processed for sure
		// Compute OR and frequency of the interval
		X_tr_aux = X_tr[tau+l]; X_par_aux = X_par[tau];
		for(j=0; j<N; j++){
			if((!X_par_aux[j]) && X_tr_aux[j]){
				X_par_aux[j] = 1; freq_par[tau]++;
			}
		}
		// If the interval is testable, increase counter of testable items and frequency-buckets and
		// check if the corrected significance threshold must be reduced
		if(ISTESTABLE(freq_par[tau])){
			// Compute cell count
			a = 0;
			for(j=0; j<N; j++) if(X_par_aux[j]) a += Y_tr[j];
			// Compute the P-value
			pval = compute_pval(a,freq_par[tau]); n_pvalues_computed++;
			// Check if the P-value is significant
			if(pval_out_file) fprintf(pval_out_file,"%lld,%lld,%lld,%lld,%e\n",l+1,tau,a,freq_par[tau],pval);
			if(pval <= delta_opt) { fprintf(sig_int_file,"%lld,%lld,%lld,%lld,%e\n",l+1,tau,a,freq_par[tau],pval); n_significant_intervals++; }
		}
		// If either the current interval or the previous one are prunable (i.e. have more than su2 ones)
		// then do NOT append the left-child to the testable queue (i.e. prune it)
		if((tau==0) || ISPRUNABLE(freq_par[tau]) || ISPRUNABLE(freq_par[tau-1])) continue;
		// Compute index of array position to be used, wrapping around if necessary
		queue_idx = testable_queue_front + testable_queue_length;
		queue_idx = (queue_idx < L) ? queue_idx : queue_idx - L;
		// Actually append children to testable queue
		testable_queue[queue_idx] = tau-1;
		// Update queue length
		testable_queue_length++;
	}
}

/* Wrapper function that encapsulates the functionality required to find significant intervals */
void find_significant_intervals(){
	long long j;
	// Give feedback to user
	#ifndef NO_VERBOSE
	printf("\n\nSCANNING DATASET FOR SIGNIFICANT INTERVALS...\n\n");
	#endif
	// Initialise the queue as empty
	testable_queue_front = 0; testable_queue_length = 0;
	// Initialise current layer index and current number of computed p-values to 0
	l = 0; n_pvalues_computed = 0; n_significant_intervals = 0;
	// Clear the current layer frequency counters
	memset(freq_par,0,L*sizeof(long long));
	// Initialise the value of the OR vectors of current layer to original dataset
	memcpy(X_par[0],X_tr[0],L*N*sizeof(char));
	// Process the upper-most layer (i.e. the layer composed of length 1 intervals)
	#ifndef NO_VERBOSE
	printf("\tProcessing layer %lld...\n",l+1);
	#endif
	process_first_layer_pvalues();
	// Artificially initialise last_tau to L-1 to ensure the first iteration of process_intervals()
	// increases the number of layers processed if the testable queue is non-empty
	last_tau = L-1;
	// Process the rest of layers (i.e. intervals of length > 1) until the pruning naturally stops the execution
	process_intervals_pvalues();
	// Report number of significant intervals found
	fprintf(summary_file,"Number of significantly associated intervals found: %lld\n",n_significant_intervals);
}

/* -------------------FUNCTIONS TO FIND THE CORRECTED SIGNIFICANCE THRESHOLD -------------------------------------- */

/* Decrease the minimum p-value threshold one level
 * Main operations that need to be performed are:
 * 1) Figure out whether we have to shrink "the W" on the left side or the right side, that is, if the current region
 *    is Sigma_{k} = [sl1,sl2] U [N-sl2,N-sl1], we need to figure out if Sigma_{k+1} is of the form
 *    Sigma_{k+1} = [sl1+1,sl2] U [N-sl2,N-sl1-1] (shrink left side) or Sigma_{k+1} = [sl1,sl2-1] U [N-sl2+1,N-sl1-1]
 *    (shrink right side). This is done with help of a binary flag that remembers which of the two types of region
 *    change happened the last time the threshold was decreased.
 * 2) Update variables sl1, sl2 and delta accordingly
 * 4) Recompute the number of testable items by removing those being excluded from the testable region due to the
 *    threshold change
 * */
void decrease_threshold(){
	int j; //Loop iterator
	int false_positives; //Number of false positives (a false positive occurs if min_pval[j] <= delta)
	if(flag){ // Flag==1 means the last call to decrease_threshold() shrunk "the W" on the left side
		sl1++; su2--; // Shrink Sigma_k on extremes of the W
		// Check what the new case will be
		if (psi[sl1] >= psi[sl2]) delta = psi[sl1];
		else{ delta = psi[sl2]; flag = 0; }
	}else{ // Flag==0 means the last call to decrease_threshold() shrunk "the W" on the right side
		sl2--; su1++; // Shrink Sigma_k on center of the W
		// Check what the new case will be
		if (psi[sl1] >= psi[sl2]){ delta = psi[sl1]; flag = 1; }
		else delta = psi[sl2];
		//No need to update LCM minimum support in this case, since sl1 remains the same
	}
	// Recompute FWER from scratch
	false_positives = 0;
	for(j=0; j<J; j++) false_positives += (min_pval[j]<=delta) ? 1 : 0;
	FWER = ((double)false_positives)/J;
}

void process_first_layer_threshold(){
	long long tau, queue_idx, i, j;
	char *X_tr_aux;
	char *Y_tr_perm_aux;
	double Tval, pval, aux, num_precomp, den_precomp;
	// Process each length 1 interval
	for(tau=0; tau<L; tau++){
		n_intervals_processed++;
		// Compute number of 1s in the interval
		X_tr_aux = X_tr[tau];
		for(j=0; j<N; j++) freq_par[tau] += X_tr_aux[j];
		// If the interval is testable...
		// Update frequency-buckets and number of testable intervals
		#ifndef NO_SINGLE_FEATURES
		if(ISTESTABLE(freq_par[tau])){
			// Precompute common parts of Chi2 statistic
			aux = ((double)freq_par[tau])/N; num_precomp = -n*aux; den_precomp = freq_par[tau]*(1-aux)*class_ratio_bin;
			if(den_precomp>0){//If den_precomp==0, then pval=1 and the min_pvals will never change
				// Compute cell-counts for all permutations
				for(i=0; i<N; i++){
					if(X_tr_aux[i]) for(j=0, Y_tr_perm_aux = Y_tr_perm[i]; j<J; j++) a_cnt[j] += Y_tr_perm_aux[j];
				}
				// Check if we have a new minimum P-value for some of the permutations
				for(j=0; j<J; j++){
					Tval = a_cnt[j]+num_precomp; Tval *= Tval; Tval /= den_precomp;
					pval = Chi2_sf(Tval,1); a_cnt[j] = 0;
					if(pval < min_pval[j]){
						// Increase FWER only if the previous minimum p-value is above current threshold
						if((pval <= delta) && (min_pval[j] > delta)) FWER += ((double)1)/J;
						// Update minimum p-value
						min_pval[j] = pval;
					}
				}
				// Update threshold until FWER condition is satisfied again
				while(FWER > alpha) decrease_threshold();
			}
		}
		#endif
		// If either the current interval or the previous one are prunable (i.e. have more than su2 ones)
		// then do NOT append the left-child to the testable queue (i.e. prune it)
		if((tau==0) || ISPRUNABLE(freq_par[tau]) || ISPRUNABLE(freq_par[tau-1])) continue;
		// Compute index of array position to be used, wrapping around if necessary
		queue_idx = testable_queue_front + testable_queue_length;
		queue_idx = (queue_idx < L) ? queue_idx : queue_idx - L;
		// Actually append children to testable queue
		testable_queue[queue_idx] = tau-1;
		// Update queue length
		testable_queue_length++;
	}
}

void process_intervals_threshold(){
	long long tau, queue_idx, i, j;
	char *X_tr_aux, *X_par_aux;
	char *Y_tr_perm_aux;
	double Tval, pval, aux, num_precomp, den_precomp;
	// While testable-interval queue is not empty, continue to process intervals
	while(testable_queue_length){
		// Pop a testable interval from the queue
		tau = testable_queue[testable_queue_front];
		testable_queue_front = (testable_queue_front<(L-1)) ? testable_queue_front + 1: 0;
		testable_queue_length--;
		// Check if we have started processing a new layer by detecting non-monotonicity in tau
		if(tau < last_tau) {
			l++;
			#ifndef NO_VERBOSE
			printf("\tProcessing layer %lld...\n",l+1);
			#endif
		}
		if((L_max>0) && ((l+1) > L_max)) {
			#ifndef NO_VERBOSE
			printf("\tMaximum interval length achieved at l=%lld. Stopping enumeration...\n",l+1);
			#endif
			break;
		}
		last_tau = tau;
		// Check any of the two parents is prunable, stop processing. Notice that this check is necessary
		// even if the current interval was appended to the testable queue, because the threshold and
		// testability regions might have been modified between the time in which the current interval
		// was appended to the queue and the time in which it is being processed
		if(ISPRUNABLE(freq_par[tau]) || ISPRUNABLE(freq_par[tau+1])) continue;
		n_intervals_processed++;
		// Compute OR and frequency of the interval
		X_tr_aux = X_tr[tau+l]; X_par_aux = X_par[tau];
		for(j=0; j<N; j++) if((!X_par_aux[j]) && X_tr_aux[j]){ X_par_aux[j] = 1; freq_par[tau]++;}
		// If the interval is testable, increase counter of testable items and frequency-buckets and
		// check if the corrected significance threshold must be reduced
		if(ISTESTABLE(freq_par[tau])){
			// Precompute common parts of Chi2 statistic
			aux = ((double)freq_par[tau])/N; num_precomp = -n*aux; den_precomp = freq_par[tau]*(1-aux)*class_ratio_bin;
			if(den_precomp>0){//If den_precomp==0, then pval=1 and the min_pvals will never change
				// Compute cell-counts for all permutations
				for(i=0; i<N; i++){
					if(X_par_aux[i]) for(j=0, Y_tr_perm_aux = Y_tr_perm[i]; j<J; j++) a_cnt[j] += Y_tr_perm_aux[j];
				}
				// Check if we have a new minimum P-value for some of the permutations
				for(j=0; j<J; j++){
					Tval = a_cnt[j]+num_precomp; Tval *= Tval; Tval /= den_precomp;
					pval = Chi2_sf(Tval,1); a_cnt[j] = 0;
					if(pval < min_pval[j]){
						// Increase FWER only if the previous minimum p-value is above current threshold
						if((pval <= delta) && (min_pval[j] > delta)) FWER += ((double)1)/J;
						// Update minimum p-value
						min_pval[j] = pval;
					}
				}
				// Update threshold until FWER condition is satisfied again
				while(FWER > alpha) decrease_threshold();
			}
		}
		// If either the current interval or the previous one are prunable (i.e. have more than su2 ones)
		// then do NOT append the left-child to the testable queue (i.e. prune it)
		if((tau==0) || ISPRUNABLE(freq_par[tau]) || ISPRUNABLE(freq_par[tau-1])) continue;
		// Compute index of array position to be used, wrapping around if necessary
		queue_idx = testable_queue_front + testable_queue_length;
		queue_idx = (queue_idx < L) ? queue_idx : queue_idx - L;
		// Actually append children to testable queue
		testable_queue[queue_idx] = tau-1;
		// Update queue length
		testable_queue_length++;
	}
}

/* Function to give some feedback about the computation of the significance threshold */
void output_significance_threshold(){
	long long j;
	fprintf(summary_file,"DATASET CHARACTERISTICS:\n");
	fprintf(summary_file,"\tN = %lld, n = %lld, L = %lld\n",N,n,L);
	fprintf(summary_file,"RESULTS: \n");
	fprintf(summary_file,"Intervals processed: %lld (%f%% of total)\n",n_intervals_processed,((double)(200*n_intervals_processed))/(L*(L+1)));
	fprintf(summary_file,"Maximum testable interval length: %lld\n",l+1);
	if(L_max==0) fprintf(summary_file,"Maximum interval length to be processed: unlimited\n");
	else fprintf(summary_file,"Maximum interval length to be processed: %lld\n",L_max);
	fprintf(summary_file,"Resultant testability region: [%lld,%lld] U [%lld,%lld]\n",sl1,sl2,su1,su2);
	fprintf(summary_file,"Resultant testability threshold: %e\n",delta);
	fprintf(summary_file,"\t FWER at testability threshold: %e\n",FWER);
	fprintf(summary_file,"Corrected significance threshold at level %e: %e\n",alpha, delta_opt);
	fprintf(summary_file,"\t FWER at corrected significance threshold: %e\n",FWER_opt);
	fprintf(minpvals_file,"MINIMUM P-VALS (%lld PERMUTATIONS)\n",J);
	for(j=0;j<(J-1);j++) fprintf(minpvals_file,"%e,",min_pval[j]);
	fprintf(minpvals_file,"%e\n",min_pval[J-1]);
}

/* Wrapper function that encapsulates the functionality required to find the corrected significance threshold */
void compute_corrected_significance_threshold(){
	long long j, idx_max;
	// Give feedback to user
	#ifndef NO_VERBOSE
	printf("COMPUTING CORRECTED SIGNIFICANCE THRESHOLD...\n");
	#endif
	// Initialise the queue as empty
	testable_queue_front = 0; testable_queue_length = 0;
	// Initialise current layer index, current number of testable intervals and current number of intervals processed to 0
	l = 0; FWER = 0; n_intervals_processed = 0;
	// Initialise the value of the OR vectors of current layer to original dataset
	memcpy(X_par[0],X_tr[0],L*N*sizeof(char));
	// Process the upper-most layer (i.e. the layer composed of length 1 intervals)
	#ifndef NO_VERBOSE
	printf("\tProcessing layer %lld...\n",l+1);
	#endif
	process_first_layer_threshold();
	// Artificially initialise last_tau to L-1 to ensure the first iteration of process_intervals()
	// increases the number of layers processed if the testable queue is non-empty
	last_tau = L-1;
	// Process the rest of layers (i.e. intervals of length > 1) until the pruning naturally stops the execution
	process_intervals_threshold();
	// Sort p-values
	qsort(min_pval,J,sizeof(double),doublecomp);
	// Tentative index to corrected significance threshold
	idx_max = (int)floor(alpha*J)-1; delta_opt = min_pval[idx_max];
	// Check and correct (if necessary) boundary cases
	if(delta_opt==min_pval[idx_max+1]){
		while(min_pval[--idx_max]==delta_opt);
		delta_opt = min_pval[idx_max];
	}
	FWER_opt = floor(idx_max+1)/J;
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

/* -----------------FUNCTIONS TO SAMPLE RANDOM INTEGERS AND GENERATE RANDOM PERMUTATIONS--------------------------- */

/* Sample a random integer uniformly distributed the range [0,x)
 * Default random number generator samples integers uniformly in the range [0,RAND_MAX)
 * To sample in the range [0,x) with x < RAND_MAX, one can think of sampling an integer
 * rnd from [0,RAND_MAX) and then returning rnd % x. However, this may generate a non-uniform
 * distribution unless RAND_MAX is a multiple of x. To fix that, we combine that basic idea
 * with rejection sampling, rejecting any value rnd greater or equal than RAND_MAX - RAND_MAX *x.
 * This is the same as ensuring that the "effective" RAND_MAX is an exact multiple of x, so that
 * returning rnd % x leads to a uniform sampling scheme in the range [0,x)
 * The seed of the random number generator should have been initialised externally
 * */
static int rand_int(int x){
	int rnd;
	int limit = RAND_MAX - RAND_MAX % x;

	do{
		rnd = rand();
	}while(rnd >= limit);
	return rnd % x;
}

void randperm(char *buffer, char *src, int l, FILE *f_perm){
	int i,j; // Variables for looping and swapping
	char tmp; // Temp int for swapping
	//Array to store the permutation and temp int for swapping (only needed if f_perm is not NULL)
	int *perm_array, tmp_int;

	// First of all, copy the original array in the buffer
	for(i=0;i<l;i++) buffer[i] = src[i];

	// If the generated permutation is to be stored in a file, initialise memory for perm_array and int_to_str
	if(f_perm){
		perm_array = (int *)malloc(l*sizeof(int));
		if(!perm_array){
			fprintf(stderr,"Error in function randperm: couldn't allocate memory for array perm_array\n");
			exit(1);
		}
		// Initialise indices of the permutation to the identity permutation
		for(i=0;i<l;i++) perm_array[i] = i;
	}

	// Fisher-Yates algorithm
	for(i = l-1; i > 0; i--){
		// Sample a random integer in [0,i]
		j = rand_int(i + 1);
		// Swap dest[j] and dest[i]
		tmp = buffer[j];
		buffer[j] = buffer[i];
		buffer[i] = tmp;
		// If the generated permutation has to be written to the stream f_perm, keep
		// track of the indices as well
		if(f_perm){
			tmp_int = perm_array[j];
			perm_array[j] = perm_array[i];
			perm_array[i] = tmp_int;
		}
	}
	// If the generated permutation is to be stored in a file, write perm_array to the stream
	if(f_perm){
		// Write the first l-1 indices with comma as a delimiter
		for(i=0;i<(l-1);i++) fprintf(f_perm,"%d,",perm_array[i]);
		// For the last index, change the comma by a newline char
		fprintf(f_perm,"%d\n",perm_array[i]);
		// Free allocated memory
		free(perm_array);
	}
}

// Comparison function used by quicksort implementation in C (increasing sorting order)
int doublecomp(const void* elem1, const void* elem2){
    if(*(const double*)elem1 < *(const double*)elem2)
        return -1;
    return *(const double*)elem1 > *(const double*)elem2;
}

/* ----------------------------------------------------ENTRY POINT---------------------------------------------------------- */

int main(int argc, char *argv[]){
	char *tmp_filename;
	char *R_command;
    char i;
    char n_fixed_args = 0;
    char idx_fixed_args[6] = {-1,-1,-1,-1,-1,-1};
    char idx_postprocessing_folder = -1, idx_out_pvals_file = -1;


    // Process input
    for(i=1; i<argc; i++){
        if((strcmp(argv[i],"-postprocessing_folder")==0) || (strcmp(argv[i],"-pp")==0)) idx_postprocessing_folder = ++i;
        else if((strcmp(argv[i],"-pval_file")==0)) idx_out_pvals_file = ++i;
        else idx_fixed_args[n_fixed_args++] = i;
    }

	// Check input
	if(n_fixed_args != 6){
		printf("ERROR: INCORRECT SYNTAX!\n");
		printf("\tUSAGE: ./program_name X_file Y_file alpha n_perm L_max base_filename [-postprocessing_folder path_to_pyfiltering.py] [-pval_file all_pvals_file]\n");
		exit(1);
	}

	// Get time when program started to run
	t_init = measureTime();

	// INITIALISATION
	tic = measureTime();

	// First allocate memory for filename holder
	tmp_filename = (char *)malloc((strlen(argv[idx_fixed_args[5]])+512)*sizeof(char));

	// Create a file to report significant intervals
	strcpy(tmp_filename,argv[idx_fixed_args[5]]); strcat(tmp_filename,"_sigints.csv");
	if(!(sig_int_file = fopen(tmp_filename,"w"))){
		fprintf(stderr, "Error in function main when opening file %s\n",tmp_filename);
		exit(1);
	}
	// If the file was successfully create, write the file header
	fprintf(sig_int_file,"l,tau,a,x,P-value\n");

	// Create a file to report runtime information
	strcpy(tmp_filename,argv[idx_fixed_args[5]]); strcat(tmp_filename,"_timing.txt");
	if(!(timing_file = fopen(tmp_filename,"w"))){
		fprintf(stderr, "Error in function main when opening file %s\n",tmp_filename);
		exit(1);
	}

	// Create a file to report varied information such as number of intervals processed, etc.
	strcpy(tmp_filename,argv[idx_fixed_args[5]]); strcat(tmp_filename,"_summary.txt");
	if(!(summary_file = fopen(tmp_filename,"w"))){
		fprintf(stderr, "Error in function main when opening file %s\n",tmp_filename);
		exit(1);
	}

	// Create a file to store the minimum p-value for each permutation
	strcpy(tmp_filename,argv[idx_fixed_args[5]]); strcat(tmp_filename,"_minpvals.txt");
	if(!(minpvals_file = fopen(tmp_filename,"w"))){
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
	siswy_init(argv[idx_fixed_args[0]], argv[idx_fixed_args[1]], atof(argv[idx_fixed_args[2]]), atoll(argv[idx_fixed_args[3]]), atoll(argv[idx_fixed_args[4]]));

	toc = measureTime();
	time_initialisation = toc-tic;

	// Main functionality
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
	siswy_end();

	// Call R postprocessing script
	R_command = (char *)malloc((strlen(argv[idx_fixed_args[5]])+512)*sizeof(char));
	strcpy(R_command,"python "); 
    if(idx_postprocessing_folder != -1) strcat(R_command,argv[idx_postprocessing_folder]);
    strcat(R_command,"pyfiltering.py "); strcat(R_command,argv[idx_fixed_args[5]]);
	strcat(R_command,"_sigints.csv "); strcat(R_command,argv[idx_fixed_args[5]]); strcat(R_command,"_sigints_filtered.corrected.csv");
	system(R_command);
	free(R_command);

	// Return
	exit(0);
}



#endif
