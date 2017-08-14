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

bench:
	cd ./graph_benchmarks && ./run.sh

clean:
	rm -rf stinger
	cd ./graph_benchmarks && make clean

.PHONY: clean
