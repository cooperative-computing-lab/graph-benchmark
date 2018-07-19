import numpy as np
import random
import matplotlib.pyplot as plt
from collections import defaultdict
import csv


class Graph:

	adj_list = None
	edge_list = None

	# Graph constructor:
	# Creates Graph object with specified vertex scales and edge factor.
	# scale_l : number of vertices in left, N, where N = 2^scale_l
	# scale_r : number of vertices in right, P, where P = 2^scale_r
	# edge_factor : number of edges, M, where M = edge_factor * N
	# weighted (Optional) : Indicates whether or not the graph is weighted
	def __init__(self, scale_l, scale_r, edge_factor, weighted):
		self.scale_l = scale_l
		self.scale_r = scale_r
		self.edge_factor = edge_factor
		self.weighted = weighted

		
	# Generates a graph with number of vertices as 2^scale, respectively. The
	# total number of edges will be approx. edge_factor * 2^scale. 
	# Returns: edge list as numpy array of dimensions 3 x edgeCount
	def generate_graph( self, rand_probs ):
		self.bipartite = False
		N = 2 ** self.scale_l
		M = self.edge_factor * N
		noise = 0.08

		#A, B, C, and D determine the shape of the degree distribution.
		#noise determines the amount of smoothing to be applied to distribution.
		if rand_probs:
			A, B, C, D = rand_probs()
		else:
			A = 0.57; B = 0.19; C = 0.19; D = 0.05

		i = np.zeros(M)
		j = np.zeros(M)

		for c in range(self.scale_l):
			ab = A + B
			c_norm = C/(1 - (A + B))
			a_norm = A/(A + B)

			#Apply noise at each level to remove oscillations
			u = random.uniform( -noise, noise )
			A = A - ((2 * u * A)/(A + D))
			B = B + u
			C = C + u
			D = D - ((2 * u * D)/(A + D))

			i_bit = np.random.rand(1, M) > ab
			j_bit = np.greater( np.random.rand(1, M), 
												((c_norm * i_bit) + (a_norm * ~i_bit)) )
			i = i + 2**c * i_bit
			j = j + 2**c * j_bit

		if self.weighted:
			edge_list = np.ones( (3, M) )
			edge_list[2, :] = np.random.rand(1, M)
		else:
			edge_list = np.ones( (2, M) )

		edge_list[0, :] = i
		edge_list[1, :] = j

		p = np.random.permutation(N)
		edge_list[0:2, :] = p[edge_list[0:2, :].astype(int)]
		p = np.random.permutation(M)
		edge_list = edge_list[:, p]

		edge_list = remove_self(edge_list)
		edge_list = remove_duplicates(edge_list)

		self.edge_list = edge_list


	# Generates a bipartite graph with number of vertices in left and right set
	# determined by parameters scale_l and scale_r. The total number of edges
	# will be approx. edge_factor * 2^scale_l. The vertices in the left set are
	# negative to easily differentiate the two sets.
	#
	# covered: Indicates whether or not vertices from 1..(2^scale) should all
	#					 exist in the graph. Ex: [1, 3, 2, 4] vs. [1, 3, 4]
	# Returns: edge list as numpy array of dimensions 3 x edgeCount
	def generate_bipartite( self, covered, rand_probs ):
		self.bipartite = True
		N = 2 ** self.scale_l
		P = 2 ** self.scale_r
		M = ((self.edge_factor-2) * N) / 2
		noise = 0.05

		#These determine the shape of the degree distribution for left set.
		#noise determines the amount of smoothing to be applied to distribution
		#For a somewhat realistic baseline use [0.57, 0.19, 0.19, 0.05]
		if rand_probs:
			A_I, B_I, C_I, D_I = get_probs()
			A_O, B_O, C_O, D_O = get_probs()
		else:
			A_I = 0.57; B_I = 0.19; C_I = 0.19; D_I = 0.05
			A_O = 0.57; B_O = 0.19; C_O = 0.19; D_O = 0.05

		i = np.zeros(M)
		j = np.zeros(M)

		# Generate vertices present in left set
		for c in range(self.scale_l):
			ab_I = A_I + B_I
			c_norm_I = C_I/(1 - (A_I + B_I))
			a_norm_I = A_I/(A_I + B_I)

			#Apply noise at each level to remove oscillations
			u = random.uniform( -noise, noise )
			A_I = A_I - ((2 * u * A_I)/(A_I + D_I))
			B_I = B_I + u
			C_I = C_I + u
			D_I = D_I - ((2 * u * D_I)/(A_I + D_I))

			i_bit = np.random.rand(1, M) > ab_I
			j_bit = np.greater( np.random.rand(1, M),
												((c_norm_I * i_bit) + (a_norm_I * ~i_bit)) )
			i = i + 2**(c) * i_bit
			j = j + 2**(c) * j_bit

		ij = np.append(i, j, axis=1)

		a = np.zeros(M)
		b = np.zeros(M)

		# Generate vertices present in right set
		for d in range(self.scale_r):
			ab_O = A_O + B_O
			c_norm_O = C_O/(1 - (A_O + B_O))
			a_norm_O = A_O/(A_O + B_O)

			#Apply noise at each level to remove oscillations
			u = random.uniform( -noise, noise )
			A_O = A_O - ((2 * u * A_O)/(A_O + D_O))
			B_O = B_O + u
			C_O = C_O + u
			D_O = D_O - ((2 * u * D_O)/(A_O + D_O))

			a_bit = np.random.rand(1, M) > ab_O
			b_bit = np.greater( np.random.rand(1, M), 
												((c_norm_O * a_bit) + (a_norm_O * ~a_bit)) )
			a = a + 2**(d) * a_bit
			b = b + 2**(d) * b_bit

		ab = np.append(a, b, axis=1)

		if self.weighted:
			edge_list = np.ones( (3, 2*M) )
			edge_list[2, :] = np.random.rand(1, 2*M)
		else:
			edge_list = np.ones( (2, 2*M) )

		edge_list[0, :] = ij
		edge_list[1, :] = ab
		
		# Make sure to cover all elements from 1 to 2^scale for both sets
		if covered:
			edge_list = cover_sets( edge_list, N, P, self.weighted )

		# Permute the vertices and edges
		perm_N = np.random.permutation(N)
		edge_list[0, :] = -1 - perm_N[ edge_list[0, :].astype(int) ]
		perm_P = np.random.permutation(P)
		edge_list[1, :] = 1 + perm_P[ edge_list[1, :].astype(int) ]
		perm = np.random.permutation(edge_list[0].size)
		edge_list = edge_list[:, perm].astype(int)

		self.edge_list = edge_list



	# Writes adjacency list to a file in the form a csv:
	# N,P,M
	# -1,24,23
	# -32,6
	# -31,19,4,31
	# 
	# adj_list : Adjacency list to be written
	# left : number of vertices in the left set
	# right : number of vertices in the right set
	# edges : number of edges in the graph
	# filename : name of the file to write
	def write_adj_list( self, filename ):
		if self.adj_list is None:
			self.adj_list = to_adj_list(self.edge_list)

		with open(filename, 'w') as outfile:
			stats = [self.left, self.right, self.edges]
			writer = csv.writer(outfile, lineterminator='\n')
			writer.writerow(stats)
			for i in self.adj_list:
				writer.writerow([i] + self.adj_list[i])


	# Writes edge list to a file in the form of csv:
	# N,P,M
	# -16,11
	# -4,8
	# -25,15
	#
	# edge_list : Edge list to be written
	# left : number of vertices in the left set
	# right : number of vertices in the right set
	# edges : number of edges in the graph
	# filename : name of the file to write
	def write_edge_list( self, filename ):
		if self.edge_list is None:
			exit('Graph has not been generated!\n')

		with open(filename, 'w') as outfile:
			stats = [self.left, self.right, self.edges]
			writer = csv.writer(outfile, lineterminator='\n')
			writer.writerow(stats)
			for i in range(self.edges):
				#outfile.write(str(int(self.edge_list[0][i])) + ' ' 
				#							+ str(int(self.edge_list[1][i])))
				#outfile.write(' \n')
				writer.writerow(self.edge_list[:,i])


	# Prints the general statistics of the generated graph such as counts of
	# edges, vertices, and degree
	def get_stats(self):
		if self.edge_list is None:
			exit('Graph has not been generated!\n')

		count_O, count_I = degree( self.edge_list )
		self.edges = self.edge_list[0].size

		if self.bipartite == True:
			self.left = count_O.size
			self.right = count_I.size
			self.vertex_count = self.left + self.right
			self.ave_degree = 2*self.edges / float(self.vertex_count)
			self.left_degree = self.edges / float(self.left)
			self.right_degree = self.edges / float(self.right)
			print 'Average degree of entire graph:', self.ave_degree
			print 'Average degree of left:', self.left_degree
			print 'Average degree of right:', self.right_degree, '\n'
			print 'Total number of vertices:', self.vertex_count
			print 'Number of vertices in left:', self.left
			print 'Number of vertices in right:', self.right
			print 'Number of edges:', self.edges, '\n'
		else:
			self.vertex_count = count_O.size
			self.ave_degree = 2*self.edges / float(count_O.size + count_I.size)
			self.in_degree = self.edges / float(count_I.size)
			self.out_degree = self.edges / float(count_O.size)
			print 'Average degree of entire graph', self.ave_degree
			print 'Average outdegree:', self.out_degree
			print 'Average indegree:', self.in_degree
			print 'Total number of vertices:', self.vertex_count
			print 'Number of edges:', self.edges, '\n'


	# Plots a histogram of degree distribution
	def plot_histogram( self ):
		count_O, count_I = degree( self.edge_list )
		if self.bipartite == True:
			log_max_l = np.log10(np.max(count_O))
			logbins_l = np.logspace(0, log_max_l, num=500)
			log_max_r = np.log10(np.max(count_I))
			logbins_r = np.logspace(0, log_max_r, num=500)
			fig, (deg_l, deg_r) = plt.subplots(1, 2, figsize=(10, 4.5))

			deg_l.hist( count_O, bins=logbins_l, log=True, align='left' )
			deg_l.set_title('Left Set')
			deg_l.set_xlabel('Degree')
			deg_l.set_ylabel('Frequency')
			deg_l.set_xscale('log')
			deg_l.set_yscale('log')

			deg_r.hist( count_I, bins=logbins_r, log=True, align='left' )
			deg_r.set_title('Right Set')
			deg_r.set_xlabel('Degree')
			deg_r.set_ylabel('Frequency')
			deg_r.set_xscale('log')
			deg_r.set_yscale('log')

			fig.suptitle('Degree Distributions of Bipartite Graph')

		else:
			log_max = np.log10(np.max(count_O))
			logbins = np.logspace(0, log_max, num=500)

			plt.hist( count_O, bins=logbins, log=True, align='left' )
			plt.title('Out Degree Distribution')
			plt.xlabel('Degree')
			plt.ylabel('Frequency')
			plt.xscale('log')
			plt.yscale('log')

		plt.show()


	# Plots the degree distribution in a scatter plot
	def plot_distribution( self ):
		if self.bipartite == True:
			count_l, count_r = degree( self.edge_list )
			unique_l, frequency_l = np.unique(count_l, return_counts=True)
			unique_r, frequency_r = np.unique(count_r, return_counts=True)
			fig, (deg_l, deg_r) = plt.subplots(1, 2, figsize=(10, 4.5))

			deg_l.scatter( unique_l, frequency_l, s=0.6 )
			deg_l.set_title( 'Left Set' )
			deg_l.set_xlabel('Degree')
			deg_l.set_ylabel('Frequency')
			deg_l.set_xscale('log')
			deg_l.set_yscale('log')

			deg_r.scatter( unique_r, frequency_r, s=0.6 )
			deg_r.set_title( 'Right Set' )
			deg_r.set_xlabel('Degree')
			deg_r.set_ylabel('Frequency')
			deg_r.set_xscale('log')
			deg_r.set_yscale('log')

			fig.suptitle('Degree Distributions of Bipartite Graph')

		else:
			count_O, count_I = degree(self.edge_list)
			unique_x, frequency_x = np.unique(count_O, return_counts=True)
			fig, ax = plt.subplots(figsize=(5.35, 4.13))
			color = np.asarray([[0, 0.2, 1]])
			lab_size = 15
			tick_w = 0.45
			min_tick = 7
			maj_tick = 15

			ax.scatter( unique_x, frequency_x, s=14, marker='.', c=color)
			ax.set_xlabel('Degree (d)', fontsize=lab_size)
			ax.set_ylabel('Frequency', fontsize=lab_size)
			ax.set_xscale('log')
			ax.set_yscale('log')
			ax.tick_params(axis='y', which='major', direction='in', labelsize=lab_size,
										 length=maj_tick, width=tick_w)
			ax.tick_params(axis='x', which='major', direction='in', labelsize=lab_size,
										 length=maj_tick, width=tick_w)
			ax.tick_params(axis='y', which='minor', direction='in', labelsize=lab_size,
										 length=min_tick, width=tick_w )
			ax.tick_params(axis='x', which='minor', direction='in', labelsize=lab_size,
										 length=min_tick, width=tick_w)
			ax.xaxis.set_ticks_position('both')
			ax.yaxis.set_ticks_position('both')
			ax.set_ylim(0, 10**5)
			ax.set_xlim(0, 10**4)
			#fig.suptitle('Degree Distribution of Graph')
			fig.savefig('../Visuals/temp.png', dpi=fig.dpi)

		plt.tight_layout()
		plt.show()


############################# HELPER FUNCTIONS #############################

# All vertices in the set from 1..(2^scale_x) should have at minimum one edge
# in the edgelist. This covers both the left and right sets without modifying
# the degree distribution much.
def cover_sets( edge_list, N, P, weighted ):
	left_cover = np.random.permutation(N).reshape((1, N))
	left_cover = np.append(left_cover, np.random.randint(1, P, (1, N)), axis=0)
	
	right_cover = np.random.randint(0, N, (1, P))
	all_right = np.random.permutation(P).reshape((1,P))
	right_cover = np.append(right_cover, all_right, axis=0)
	
	if weighted:
		left_cover = np.append(left_cover, np.random.rand(1, N), axis=0)
		right_cover = np.append(right_cover, np.random.rand(1, P), axis=0)

	edge_list = np.append(edge_list, left_cover, axis=1)
	edge_list = np.append(edge_list, right_cover, axis=1)
	return edge_list


def remove_duplicates( edge_list ):
	edge_list = np.transpose(edge_list[0:2, :])
	edge_list = np.ascontiguousarray(edge_list)
	unduped = np.unique(edge_list.view([('', edge_list.dtype)]
																		 * edge_list.shape[1]))
	unduped = unduped.view(edge_list.dtype). \
						reshape((unduped.shape[0], edge_list.shape[1]))
	edge_list = np.transpose(unduped)
	return edge_list


# Removes self edges in an edge list
def remove_self( edge_list ):
	edge_list = edge_list[:, edge_list[0,:] != edge_list[1,:]]
	return edge_list


# Returns an edgelist that has been made undirected
def make_undirected( edge_list ):
	flipped_el = np.asarray([edge_list[1,:], edge_list[0,:], edge_list[2,:]])
	edge_list = np.append(edge_list, flipped_el, axis=1)
	return edge_list


# Returns adjacency representation of edgelist
# NOTE: By using sets, duplicate edges are eliminated from the graph
# edge_list : edge list to be converted to adjacency list
def to_adj_list( edge_list ):
	adj_list = defaultdict(set)
	for i in edge_list.T:
		adj_list[int(i[0])].add(int(i[1]))

	for i in adj_list:
		adj_list[i] = list(adj_list[i])

	return adj_list


# Returns the total number of edges in an adjacency list
def count_edges( adj_list ):
	count = sum(len(v) for v in adj_list.itervalues())
	return count


# Returns the outdegree of the left set and indegree of the right set.
def degree( edge_list ):
	unique_i, counts_i = np.unique(edge_list[0,:], return_counts=True)
	unique_j, counts_j = np.unique(edge_list[1,:], return_counts=True)
	return counts_i, counts_j

# Returns random probabilities for the R-MAT algorithm.
def get_probs():
	a = random.uniform(0.45, 0.6)
	b = c = random.uniform(0.1, 0.2)
	d = 1 - (a + b + c)
	return a, b, c, d