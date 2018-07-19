from collections import defaultdict, deque
import sys

INT_MAX = sys.maxint

class Hopcroft_Karp :

	def read_graph( self, filename ) :
		adj_list = defaultdict(list)

		with open(filename, 'r') as infile :
			stats = infile.readline()
			self.left = int(stats.split(',')[0])
			self.right = int(stats.split(',')[1])
			self.edges = int(stats.split(',')[2])

			for line in infile :
				line = line.rstrip('\n')
				left_vert = int(line.split(',')[0])
				right_verts = map(int, line.split(',')[1:])
				adj_list[left_vert] = right_verts

		self.adj_list = adj_list


	# Find maximum cardinality matching of a bipartite graph (U,V,E).
	# The output is a triple (M,A,B) where M is a
	# dictionary mapping members of V to their matches in U, A is the part
	# of the maximum independent set in U, and B is the part of the MIS in V.
	# The same object may occur in both U and V, and is treated as two
	# distinct vertices if this happens.
	# graph : dictionary mapping members of U to a list of their neighbors in V
	def bipartite_match( self ):
		# pairU[u] stores pair of u in matching where u
		# is a vertex on left side of Bipartite Graph.
		# If u doesn't have any pair, then pairU[u] is 0
		self.pairU = [0] * (self.left + 1)

		# pairV[v] stores pair of v in matching. If v
		# doesn't have any pair, then pairU[v] is NIL
		self.pairV = [0] * (self.right + 1)

		# dist[u] stores distance of left side vertices
		# dist[u] is one more than dist[u'] if u is next
		# to u' in augmenting path
		self.dist = [0] * (self.left + 1)

		result = 0

		while self.bfs() :
			for i in range(1, self.left + 1) :
				if self.pairU[i] == 0 and self.dfs(i) :
					result+=1

		return result


	def bfs( self ) :
		Q = deque()

		for i in range(1, self.left + 1) :
			if self.pairU[i] == 0 :
				self.dist[i] = 0
				Q.append(i)
			else :
				 self.dist[i] = INT_MAX

		self.dist[0] = INT_MAX
	 
		while Q :
			u = Q.popleft()

			if self.dist[u] < self.dist[0] :
				tmp = self.adj_list[-u]

				for i in range(0, len(tmp)) :
					v = tmp[i]
				
					if self.dist[self.pairV[v]] == INT_MAX :
						self.dist[self.pairV[v]] = self.dist[u] + 1
						Q.append(self.pairV[v])


		return self.dist[0] != INT_MAX


	def dfs( self, u ) :
		if( u != 0 ) :
			tmp = self.adj_list[-u]

			for i in range(0, len(tmp)) :
				v = tmp[i]
				if( self.dist[self.pairV[v]] == self.dist[u] + 1 ) :
					if( self.dfs( self.pairV[v] ) == 1 ) :
						self.pairV[v] = u
						self.pairU[u] = v
						return True

			self.dist[u] = INT_MAX
			return False
		return True


	# Utility function to print the matchings found from the matching algorithm.
	# matching : dictionary with matching[v] containing vertex matched to v
	def print_matching(self):
		for key in self.adj_list :
			print "%d\t:" % key,
			print "\t%d" % self.pairU[-key]
		print "\n"