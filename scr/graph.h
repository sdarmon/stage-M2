#ifndef DEF_GRAPH
#define DEF_GRAPH

using namespace std;
#include <string>
#include <vector>


class Node {
public:
    int val, weight;
    string label;
    Node(int v, int w, string l);
};

class Neighbor {
public:
    int val, weight;
    char label[2];
    Neighbor(int v, int w, char* l);
};


class LstNode {
public:
    vector<Neighbor> adjVec;
    LstNode();
    LstNode(vector<Neighbor>& A);
};



class Edge {
public:
    int start, end, weight;
    char label[2];
    Edge(int s, int e, int w, char* l);
};




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