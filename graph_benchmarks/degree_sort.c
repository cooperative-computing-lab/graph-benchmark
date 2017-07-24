/*
 Luigi Grazioso
 April 17th, 2017
 */

#include <stdio.h>
#include <stdlib.h>
#include "graph.h"
#include "degree_sort.h"
#include <hbwmalloc.h>
void degree_sort(struct Graph* graph){
    int v;
    int *bins;
    int *degrees;
    int threshold;
    double percentage = .75;			    //*****
    int quantity = graph->V*percentage;		    //*****
    printf("quantity: %d \n", quantity);
    struct Graph* sorted_graph;
    int *hbwNodes;

    bins = malloc(graph->V*sizeof *bins); // Max possible degree is |N| - 1
    degrees = malloc(graph->V*sizeof *degrees);
    hbwNodes = malloc(graph->V*sizeof *hbwNodes);

    for(int x = 0; x < graph->V; x++)
    {
        bins[x] = 0;
	hbwNodes[x]=0;
    }

    //Calculates the degree of each vertex and stores it in an array
    for (v = 0; v < graph->V; ++v){
	degrees[v] = 0;
	struct AdjListNode* pCrawl = graph->array[v].head;
	while (pCrawl) {
	    degrees[v]++;
	    pCrawl = pCrawl->next;
	}
    }
	printf("calculated degrees:\n", quantity);

    // Perform a counting sort to find the threshold for selecting the top quantity nodes
    for(int nodei = 0; nodei < graph->V; nodei++)
        bins[degrees[nodei]]++;

    for(int nodei = 1; nodei < graph->V; nodei++)
        bins[nodei] += bins[nodei - 1];

    threshold = graph->V; // You must have this many nodes to be in the HBW memory

    for(int nodei = 0; nodei < graph->V; nodei++) {
        if(graph->V - bins[nodei] < quantity) {
            threshold = nodei;
            break;
        }
    }

    printf("threshold: %i\n", threshold);

    for(int nodei = 0; nodei < graph->V; nodei++) {
        hbwNodes[nodei] = degrees[nodei] >= threshold;
//        printf("node %i with degree %i does %sgo to HBW\n",
//		nodei, degrees[nodei], hbwNodes[nodei] ? "" : "NOT ");
    }

	printf("sorted: %p \n", graph);
     graph = copytohbw(graph, hbwNodes, quantity);
      printf("cpy to hbw %p \n", graph);
    // free(graph);
    //*graph = *sorted_graph;

    free(bins);
    free(degrees);
    free(hbwNodes);
}

struct Graph* copytohbw(struct Graph* graph, int hbwNodes[], int quantity){
    int v;
    struct Graph* sorted_graph = createGraph_hbw(graph->V);

    for (v = 0; v < graph->V; ++v){
	struct AdjListNode* pCrawl = graph->array[v].head;
	//printf("\n Adjacency list of vertex %d\n head ", graph->array[v].Node);
	while (pCrawl){
	    if (checkLoc(hbwNodes, v, quantity) == 1){
		addEdgehbw(sorted_graph, v, pCrawl->dest);
	    } else {
		addEdgelbw(sorted_graph, v, pCrawl->dest);
	    }
	    pCrawl = pCrawl->next;
	}
	//printf("\n");
    }

    //printGraph(sorted_graph);

    return sorted_graph;
}

int checkLoc(int* hbwNodes, int node, int quantity){
    int i;
	return (hbwNodes[node] == 1);
}


// A utility function to create a new adjacency list node in hbwm
struct AdjListNode* newAdjListNodehbw(int dest){
    struct AdjListNode* newNode =
    (struct AdjListNode*) hbw_malloc(sizeof(struct AdjListNode));
    newNode->dest = dest;
    newNode->next = NULL;
    return newNode;
}

// A utility function to create a new adjacency list node in lbwm
struct AdjListNode* newAdjListNodelbw(int dest){
    struct AdjListNode* newNode =
    (struct AdjListNode*) malloc(sizeof(struct AdjListNode));
    newNode->dest = dest;
    newNode->next = NULL;
    return newNode;
}

// Adds an edge to an undirected graph hbwm
void addEdgehbw(struct Graph* graph, int src, int dest){

    // Add an edge from src to dest.  A new node is added to the adjacency
    // list of src.  The node is added at the begining
    struct AdjListNode* newNode = newAdjListNodehbw(dest);
    newNode->next = graph->array[src].head;
    graph->array[src].head = newNode;
}

// Adds an edge to an undirected graph lbwm
void addEdgelbw(struct Graph* graph, int src, int dest){

    // Add an edge from src to dest.  A new node is added to the adjacency
    // list of src.  The node is added at the begining
    struct AdjListNode* newNode = newAdjListNodelbw(dest);
    newNode->next = graph->array[src].head;
    graph->array[src].head = newNode;
}


// Driver program to test above functions
// 
/*
int main()
{
int V = 5;
struct Graph* graph = createGraph(V);
struct Graph* sorted_graph;
addEdge(graph, 0, 1);
addEdge(graph, 0, 4);
addEdge(graph, 1, 2);
addEdge(graph, 1, 3);
addEdge(graph, 1, 4);
addEdge(graph, 2, 3);
addEdge(graph, 3, 4);
// print the adjacency list representation of the above graph
printGraph(graph);
printf("SORTING\n");
sorted_graph = degree_sort(graph);
printGraph(sorted_graph);

return 0;
}
*/
