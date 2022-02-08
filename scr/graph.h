#ifndef DEF_GRAPH
#define DEF_GRAPH

using namespace std;
#include <string>
#include <vector>


class embNode {
public:
    int val, weight;
    char* label;
    embNode(int v, int w, char* l);
};



class adjNode {
public:
    vector<embNode> adjVec;
    adjNode();
    adjNode(vector<embNode>& A);
};



class edge {
public:
    int start, end, weight;
    char label[2];
    edge(int s, int e, int w, char* l);
};




class Graph
{
    public:

    int N;  // number of nodes in the graph
    int M; //number of edges in the graph
    vector<adjNode> Adj;                //adjacency list as a vector of adjNodes
    vector<embNode> Vertices;         //the vertices list as a vector of embNodes
    void add(edge e);
    Graph(vector<embNode>& vertices, vector<edge>& edges);
    private:

    void insEdgeAdjList(int start,int end, int weight, char* labelEdge);
};

#endif //def_graph