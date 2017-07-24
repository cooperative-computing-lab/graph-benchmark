#ifndef DEG_SORT_H
#define DEG_SORT_H



void degree_sort(struct Graph* graph);

struct Graph* copytohbw(struct Graph* graph, int hbwNodes[], int quantity);

int checkLoc(int hbwNodes[], int node, int quantity);

void addEdgehbw(struct Graph* graph, int src, int dest);

void addEdgelbw(struct Graph* graph, int src, int dest);

struct AdjListNode* newAdjListNodehbw(int dest);

struct AdjListNode* newAdjListNodelbw(int dest);

#endif
