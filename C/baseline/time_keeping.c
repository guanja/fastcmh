#ifndef _time_keeping_c_
#define _time_keeping_c_

/* LIBRARY INCLUDES FOR MEASURING EXECUTION TIME */
#include <time.h> //Already included in original LCM source code above
#include <sys/time.h>
#include <sys/resource.h>
/* END OF LIBRARY INCLUDES FOR MEASURING EXECUTION TIME */

/* GLOBAL VARIABLES (TIME SPENT) */
double t_init=0, t_end=0;
double time_initialisation=0;
double time_IO=0;
double time_comp_threshold=0;
double time_comp_significant_intervals=0;
double tic,toc;
extern FILE *timing_file; // File to output information about the runtime of the algorithm, defined in the significant_interval_search_* C files

// Measure running time
double measureTime(){
  struct rusage t;
  struct timeval tv;
  getrusage(RUSAGE_SELF, &t);
  tv = t.ru_utime;
  return tv.tv_sec + (double)tv.tv_usec * 1e-6;
}

// Measure peak memory usage
size_t measurePeakMemory(){
  struct rusage t;
  getrusage(RUSAGE_SELF, &t);
  return (size_t)t.ru_maxrss;
}

// Display execution time and memory consumption
void profileCode(){
	size_t peak_memory;
	fprintf(timing_file,"CODE PROFILING\n");
	fprintf(timing_file,"Total execution time: %f (s).\n",t_end-t_init);
	fprintf(timing_file,"\tInitialisation time: %f (s).\n",time_initialisation);
	fprintf(timing_file,"\t\tFile I/O time: %f (s).\n",time_IO);
	fprintf(timing_file,"\tTime to compute corrected significance threshold: %f (s).\n",time_comp_threshold);
	fprintf(timing_file,"\tTime to find significant intervals: %f (s).\n",time_comp_significant_intervals);
	peak_memory = measurePeakMemory();
	fprintf(timing_file,"Peak memory usage: %zu (bytes)\n",peak_memory);
}

#endif
