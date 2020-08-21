all: single meta baseline
all_nv: single_nv meta_nv baseline_nv
single: non_wy wy
single_nv: non_wy_nv wy_nv
meta: meta_non_wy
meta_nv: meta_non_wy_nv
baseline: baseline_p
baseline_nv: baseline_p_nv
non_wy:
	gcc -O3 -o compiled/significant_interval_search_exact single_analysis/non_wy/significant_interval_search_exact.c -lm
	gcc -O3 -o compiled/significant_interval_search_chi single_analysis/non_wy/significant_interval_search_chi.c -lm
non_wy_nv:
	gcc -DNO_VERBOSE -O3 -o compiled/significant_interval_search_exact single_analysis/non_wy/significant_interval_search_exact.c -lm
	gcc -DNO_VERBOSE -O3 -o compiled/significant_interval_search_chi single_analysis/non_wy/significant_interval_search_chi.c -lm
wy:
	gcc -O3 -o compiled/significant_interval_search_wy_exact single_analysis/wy/significant_interval_search_wy_exact.c -lm
	gcc -O3 -o compiled/significant_interval_search_wy_chi single_analysis/wy/significant_interval_search_wy_chi.c -lm
wy_nv:
	gcc -DNO_VERBOSE -O3 -o compiled/significant_interval_search_wy_exact single_analysis/wy/significant_interval_search_wy_exact.c -lm
	gcc -DNO_VERBOSE -O3 -o compiled/significant_interval_search_wy_chi single_analysis/wy/significant_interval_search_wy_chi.c -lm
meta_non_wy:
	gcc -g -o compiled/significant_interval_search_meta_cmh meta_analysis/non_wy/significant_interval_search_meta_cmh.c -lm
	gcc -O3 -o compiled/significant_interval_search_meta_fisher meta_analysis/non_wy/significant_interval_search_meta_fisher.c -lm
meta_non_wy_nv:
	gcc -DNO_VERBOSE -O3 -o compiled/significant_interval_search_meta_cmh meta_analysis/non_wy/significant_interval_search_meta_cmh.c -lm
	gcc -DNO_VERBOSE -O3 -o compiled/significant_interval_search_meta_fisher meta_analysis/non_wy/significant_interval_search_meta_fisher.c -lm
baseline_p:
	gcc -O3 -o compiled/brute_interval_search_exact baseline/brute_interval_search_exact.c -lm
	gcc -O3 -o compiled/brute_interval_search_meta_cmh baseline/brute_interval_search_meta_cmh.c -lm
	gcc -O3 -o compiled/significant_interval_search_meta_cmh_no_heuristic baseline/significant_interval_search_meta_cmh_no_heuristic.c -lm
baseline_p_nv:
	gcc -DNO_VERBOSE -O3 -o compiled/brute_interval_search_exact baseline/brute_interval_search_exact.c -lm
	gcc -DNO_VERBOSE -O3 -o compiled/brute_interval_search_meta_cmh baseline/brute_interval_search_meta_cmh.c -lm
	gcc -DNO_VERBOSE -O3 -o compiled/significant_interval_search_meta_cmh_no_heuristic baseline/significant_interval_search_meta_cmh_no_heuristic.c
