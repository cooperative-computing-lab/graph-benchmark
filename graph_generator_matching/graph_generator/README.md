Scale-free graph generator using similar techniques to R-MAT Graph500 generator. The main differences from the Graph500 generator are the addition of noise to the process for smoother distributions of edges and extending the idea to generate bipartite graphs.

The following command will generate a bipartite graph:

`python generate.py <scale_l> <scale_r> <edge_factor> -o <outfile> --cover --weight --visual`

The size of the graph is determined by the first three mandatory arguments:
* scale_l : Scale of the left set of vertices, N, where N = 2<sup>scale_l</sup>
* scale_r : Scale of the right set of vertices, P, where P = 2<sup>scale_r</sup>
* edge_factor : Determines number of edges, M, where M = edge_factor * N

The remaining optional arguments determine the following:
* -o, -\-outfile : Name of file to write graph/s to. Will be preceded by 'adj\_' or 'edge\_' depending on the type of graph to be written.
* -c, -\-cover : The sets of vertices will contain all elements from 1 to N, without missing elements.
* -w, -\-weight : The graphs will have weighted edges.
* -v, -\-visual : Pyplot will be utilized to visualize the degree distribution of the generated graph.
* -r, -\-random : Whether or not to generate random probabilities within R-MAT.

**Notes about -\-cover argument**: If the -\-cover argument is not used, the size of the generated graph may fluctuate greatly around your original input parameters due to the stochasticity and duplication of edges. Using it will generate a graph closer to your input parameters, however, it may cause more possible matches than are desirable.

**Adjusting the R-MAT probabilities**: Adjustments to A, B, C, D, and a random noise variable within the graph generator greatly affect the distribution of edges within the graph. The default values are:
* A = 0.57
* B = 0.19
* C = 0.19
* D = 0.05
* noise = 0.05

Increasing the noise leads to smoother distributions. A good way to find which values fit your needs is by utilizing the -\-visual argument to see the curve.