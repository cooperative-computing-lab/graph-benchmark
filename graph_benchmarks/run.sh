#!/bin/bash

./bin/benchmark jaccard none soc-Slashdot0811.txt > results/jacSlashdotNone

./bin/benchmark jaccard degree soc-Slashdot0811.txt > results/jacSlashdotDegree

./bin/benchmark cluster_coeff none soc-Slashdot0811.txt > results/clustSlashdotNone

./bin/benchmark cluster_coeff degree soc-Slashdot0811.txt > results/clustSlashdotDegree

