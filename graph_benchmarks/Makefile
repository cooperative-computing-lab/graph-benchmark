ALL_CFLAGS := -std=c99 -D_POSIX_C_SOURCE=200809L -fopenmp -g -lmemkind $(CFLAGS)

BENCHMARK_CSRC := benchmark.c cluster_coeff.c graph.c null_sort.c jaccard.c degree_sort.c
BENCHMARK_OBJS := $(addprefix obj/, $(BENCHMARK_CSRC:.c=.o))

.PHONY: clean
.PRECIOUS: %/ obj/%.o

all: bin/benchmark

bin/benchmark: $(BENCHMARK_OBJS) | bin/
	$(CC) $(ALL_CFLAGS) -o $@ $^

%/:
	mkdir -p $@

obj/%.o: %.c | obj/
	$(CC) $(ALL_CFLAGS) -c -o $@ $<

clean:
	rm -rf bin obj ./results/jac* ./results/clust*

