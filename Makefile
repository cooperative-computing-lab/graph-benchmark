STINGER = ./stinger
BENCH = ./graph_benchmarks/bin

all: $(STINGER) $(BENCH)

./stinger:
	git clone git@github.com:stingergraph/stinger.git
	cp ./graph_benchmarks/stinger_addons/adamic_adar_test.* ./stinger/src/tests/adamic_adar_test
	cd ./stinger && mkdir ./build && cd ./build && ccmake ..
	cd ./stinger/build && make

./graph_benchmarks/bin:
	cd ./graph_benchmarks && make

build: $(STINGER) $(BENCH)

stinger: $(STINGER)

bench: $(STINGER) $(BENCH)
	export KMP_AFFINITY=verbose,granularity=core,compact
	export OMP_NUM_THREADS=256
	cd ./graph_benchmarks && ./run.sh

adar: $(STINGER) $(BENCH)
	export KMP_AFFINITY=verbose,granularity=core,compact
	export OMP_NUM_THREADS=256
	./stinger/build/bin/stinger_adamic_adar_test ./graph_benchmarks/soc-Slashdot0811.txt

clean:
	rm -rf stinger
	cd ./graph_benchmarks && make clean

.PHONY: clean
