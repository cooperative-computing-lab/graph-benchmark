#ifndef GRAPH_H
#define GRAPH_H

// A structure to represent an adjacency list node
struct AdjListNode
{
      int dest;
      struct AdjListNode* next;
};

	  // A structure to represent an adjacency list
struct AdjList
{
      struct AdjListNode *head;  // pointer to head node of list
};

	       
//	A structure to represent a graph. A graph is an array of adjacency lists.
//      Size of array will be V (number of vertices in graph)
struct Graph
{
      int V;
      struct AdjList* array;
};

struct Graph *createGraph(int V);
struct Graph *createGraph_hbw(int V);
struct Graph *readGraph(FILE *f);
void addEdge(struct Graph *graph, int src, int dst);
void printGraph(struct Graph *graph);

#endif

