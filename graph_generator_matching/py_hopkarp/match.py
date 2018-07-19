import hop_karp
import argparse

def main( args ):
	G = hop_karp.Hopcroft_Karp()
	G.read_graph( args.filename )
	num_matches = G.bipartite_match()
	print 'Total number of matches:', num_matches
	#G.print_matches()


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'This script performs bipartite\
																	 matching on a graph read from an input file.\
																	 The input file should be csv containing an\
																	 adjacency list.')
	parser.add_argument('filename', help='Name of file containing adjacency list.')
	args = parser.parse_args()

	main( args )
