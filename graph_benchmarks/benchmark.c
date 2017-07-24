#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "graph.h"
#include  "degree_sort.h"
#include "cluster_coeff.h"
#include "null_sort.h"
#include "jaccard.h"
#include "omp.h"
#define NELEMS(_array) (sizeof (_array)/sizeof *(_array))

typedef struct state {
	double *(*benchfunc)(struct Graph *);
	void (*sortfunc)(struct Graph *);

	char *graphfilename;
	struct Graph *graph;
} state_t;

static const struct {
	const char *name;
	double *(*const func)(struct Graph *);
} benchmarks[] = {
	{"cluster_coeff", cluster_coeff},
	{"jaccard", jaccard_all_pairs}
};

static const struct {
	const char *name;
	void (*const func)(struct Graph *);
} sorts[] = {
	{"none", null_sort},
	{"degree", degree_sort}
};

static void process_args(state_t *state, int argc, char **argv) {
	bool needhelp;

	needhelp = false;

	for(int argi = 1; argi < argc; argi++) {
		if(!state->benchfunc) {
			for(size_t benchi = 0; benchi < NELEMS(benchmarks); benchi++)
				if(strcmp(argv[argi], benchmarks[benchi].name) == 0)
					state->benchfunc = benchmarks[benchi].func;
		} else if(!state->sortfunc) {
			for(size_t sorti = 0; sorti < NELEMS(sorts); sorti++)
				if(strcmp(argv[argi], sorts[sorti].name) == 0)
					state->sortfunc = sorts[sorti].func;
		} else if(!state->graphfilename) {
			state->graphfilename = argv[argi];
		} else needhelp = true;
	}

	// Sanity checks
	needhelp |= !state->benchfunc;
	needhelp |= !state->sortfunc;
	needhelp |= !state->graphfilename;

	if(needhelp) {
		fprintf(stderr, "usage: %s benchmark-type sort-type graph-file\n"
			"\tbenchmark-type: cluster_coeff\n"
			"\t     sort-type: degree | none\n",
			argv[0]);
		exit(EXIT_FAILURE);
	}
}

static void load_graph(state_t *state) {
	FILE *f;

	// Open the graph file
	if(f = fopen(state->graphfilename, "rb"), !f) {
		fprintf(stderr, "error: cannot open '%s' for reading", state->graphfilename);
		exit(EXIT_FAILURE);
	}

	// Parse the graph from the file
	state->graph = readGraph(f);

	// Clean up
	fclose(f);
}

static void run_benchmark(state_t *state) {
	double *outputs;
	struct timespec t0, t1;
	
	//clock_gettime(CLOCK_MONOTONIC, &t0);
	state->sortfunc(state->graph); // don't factor sort into exec time
	for(int x = 16; x <= 257; x+=16)
	{
		omp_set_num_threads(x);
		clock_gettime(CLOCK_MONOTONIC, &t0);
	
		outputs = state->benchfunc(state->graph);

		clock_gettime(CLOCK_MONOTONIC, &t1);

//	for(int i = 0; i < state->graph->V; i++)
//		printf("%i: %.6f\n", i, outputs[i]);
		if(outputs)
		{
			free(outputs);
		}
	

		printf("OpenMP %d: %.6f\n", x, (t1.tv_sec + t1.tv_nsec/1e9) - (t0.tv_sec + t0.tv_nsec/1e9));
	}
}

int main(int argc, char **argv) {
	state_t *state = &(state_t) {
		.benchfunc = NULL,
		.sortfunc = NULL,

		.graphfilename = NULL,
		.graph = NULL
	};

	process_args(state, argc, argv);

	load_graph(state);

	run_benchmark(state);

	return 0;
}

