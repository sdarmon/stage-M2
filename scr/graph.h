#ifndef DEF_GRAPH
#define DEF_GRAPH

using namespace std;
#include <string>
#include <vector>



// ===========================================================================
//                                  Structures
// ===========================================================================

/* 
 * Noeud d'un graphe 
 * Il est défini par une valeur val et par un label (string).
 * On rajoute aussi un argument weight pour lui donner un poids (variable)
 */
class Node {
public:
    int val, weight;
    string label;
    Node(int v, int w, string l);
};


/* 
 * Noeud dans la liste d'adjacence d'un graphe. 
 * Il est défini par une valeur val et par un label (ici deux caractères correspondant au type de l'arête).
 * On rajoute aussi un argument weight pour lui donner un poids (variable)
 */
class Neighbor {
public:
    int val, weight;
    char label[2];
    Neighbor(int v, int w, char* l);
};

/* 
 * Liste des noeuds adajacents à un certain noeud donné. 
 * Il est défini seulement par un vecteur de Neighbor.
 * Je pense qu'ici il serait préférable de ne pas utiliser une classe...
 */
class LstNode {
public:
    vector<Neighbor> adjVec;
    LstNode();
    LstNode(vector<Neighbor>& A);
};


/* 
 * Arête d'un graphe. 
 * Il est défini par une valeur de départ (start), d'arrivé (end) et un label (ici deux caractères correspondant au type de l'arête).
 * On rajoute aussi un argument weight pour lui donner un poids (variable)
 */
class Edge {
public:
    int start, end, weight;
    char label[2];
    Edge(int s, int e, int w, char* l);
};



/* 
 * Arête d'un graphe. 
 * Il est défini par une valeur de départ (start), d'arrivé (end) et un label (ici deux caractères correspondant au type de l'arête).
 * On rajoute aussi un argument weight pour lui donner un poids (variable)
 */
class Graph
{
    public:

    //---------------
    // Variables
    //---------------

    int N;  // number of nodes in the graph
    int M; //number of edges in the graph
    vector<LstNode> Adj;                //adjacency list as a vector of adjNodes
    vector<Node> Vertices;         //the vertices list as a vector of embNodes

    //---------------
    // Constructor
    //---------------

    Graph(vector<Node>& vertices, vector<Edge>& edges);

    //---------------
    // Methodes
    //---------------

    void add(Edge &e);
    void add(int start,int end, int weight, char* labelEdge);
    void add(Node &v);
    void add(int val,int weight, string labelEdge);

    Neighbor* link(int start, int end);

    vector<Neighbor>* Neighbors(int n);
    void weighing();
    int BFSCount(vector<int> &rayons, int acc,vector<Neighbor*> &aVoir,vector<int> &vu);
    void weighingANode(int source, int rayon);
    void weighingAllNodes(int rayon);
};


//Other functions
int count_nb_lines( FILE* file );
void read_node_file( FILE* node_file, vector<Node>& seqs);
void read_edge_file( FILE *edge_file, vector<Edge>& edges );

void printGraphVertices(Graph& G);
void printVertices(vector<Node>& V);
void printGraphEdges(Graph& G);
void printEdges(vector<Edge>& E);

#endif //def_graph