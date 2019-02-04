import graph
import time
import sys
import argparse

def main( args ):

	# Setup arguments necessary for graph
	scale_l = args.scale_l
	scale_r = args.scale_r
	edge_factor = args.edge_factor
	weighted = args.weighted
	covered = args.covered
	visual = args.visual
	rand_probs = args.rand_probs
	if args.output is not None:
		adj_outfile = 'adj_' + args.output
		edge_outfile = 'edge_' + args.output
	else:
		adj_outfile = 'adj_%dx%dx%d.csv' % (scale_l, scale_r, edge_factor)
		edge_outfile = 'edge_%dx%dx%d.csv' % (scale_l, scale_r, edge_factor)

	print('\nCovered :', covered)
	print('Random :', rand_probs)
	print('Weighted :', weighted, '\n')


	#############################################################################
	############################ Create and Generate ############################
	#############################################################################

	# Create graph object with input parameters
	bipart = graph.Graph(scale_l, scale_r, edge_factor, weighted=weighted)

	# Generate bipartite graph with option of coverting all vertices or not
	bipart.generate_bipartite(covered=covered, rand_probs=rand_probs)

	# Get the stats from the graph that was generated
	bipart.get_stats()

	# Write the graph to a file in the form of adjacency list or edge list
	bipart.write_adj_list(adj_outfile)
	#bipart.write_edge_list(edge_outfile)

	# Produce visuals of the distribution of the generated graph
	if args.visual is True:
		bipart.plot_distribution()
		bipart.plot_histogram()

	#############################################################################
	#############################################################################
	#############################################################################

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description = 'This script generates \
					scale-free graphs of user specified scale. The graphs can be \
					represented as adjacency lists or edge lists. These graphs may then \
					be written to files of your specification.')
	parser.add_argument('scale_l', help='Scale of the left set of vertices, N,\
											where N = 2^scale_l.', type=int)
	parser.add_argument('scale_r', help='Scale of the right set of vertices, P,\
											where P = 2^scale_r.', type=int)
	parser.add_argument('edge_factor', help='Determines number of edges, M, where\
											 M = edge_factor * N.', type=int)
	parser.add_argument('-o', '--output', action='store',
										 dest='output', type=str, metavar='',
										 help='Name of file to write graph to. Will be appended to\
										 edge_ or adj_ depending on type of graph.')
	parser.add_argument('-c', '--cover', action='store_true',
											dest='covered', help='Whether or not the vertex sets\
											should cover all elements from [-N, -1] and [1, P]. NOTE:\
											Generates graphs that tend to be consistent with input\
											parameters, but forces a higher percentage of matches\
											than without.')
	parser.add_argument('-w', '--weight', action='store_true', dest='weighted',
											help='Whether or not the graph to be generated should \
											have weighted edges.')
	parser.add_argument('-v', '--visual', action='store_true', dest='visual',
											help='Uses pyplot to visualize the degree distibution of\
											the generated graph.')
	parser.add_argument('-r', '--random', action='store_true', dest='rand_probs',
											help='Randomly determines a, b, c, d probabilities during\
											generation.')

	args = parser.parse_args()

	main( args )
