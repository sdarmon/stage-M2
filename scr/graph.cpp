// ===========================================================================
//                                  Constructors
// ===========================================================================
/* All the data is coded with int data type. For graph with more than 2147483648
* (> 2x10^9) nodes or edges, you will get trouble!
* Instead of int, use unsigned int (for the same weight, extend this upper bound
* to 4x10^9) or long long int (for twice the weight, extention to 2^63 ~ 10^19).
*
*/

// ===========================================================================
//                               Include Libraries
// ===========================================================================
#include <list>
#include <vector>
#include <iostream>
#include <string>
#include <stdio.h>
#include <cstring>

// ===========================================================================
//                             Include Project Files
// ===========================================================================
#include "graph.h"


// ===========================================================================
//                             Declare Used Namespaces
// ===========================================================================
using namespace std;
#define MAX 1024

// ===========================================================================
//                                  Structures
// ===========================================================================

embNode::embNode(int v, int w, char* l){
        val = v;
        weight = w;
        label = l;
    }

// stores adjacency list items
adjNode::adjNode() {
        adjVec.clear();
    }
adjNode::adjNode(vector<embNode>& A) {
        adjVec = A;
    }


// structure to store edges
edge::edge(int s, int e, int w, char* l) {
        start = s;
        end = e;
        weight = w;
        strcpy(label,l);
    }

// ===========================================================================
//                                  Class
// ===========================================================================



void Graph::insEdgeAdjList(int start,int end, int weight, char* labelEdge)   {
        embNode node(end,weight,labelEdge);
       (Adj[start].adjVec).push_back(node);
    }


Graph::Graph(vector<embNode>& vertices, vector<edge>& edges)  {

        N = vertices.size();
        M = edges.size();
        Vertices = vertices;
        Adj.clear();

        // initialize the adjNode for all vertices
        for (int i = 0; i < vertices.size(); i++){
            adjNode adj;
            Adj.push_back(adj);
            (Adj.back().adjVec).clear();
        }

        // construct directed graph by adding edges to it
        for (unsigned i = 0; i < edges.size(); i++)  {
            int start = edges[i].start;
            int end = edges[i].end;
            int weight = edges[i].weight;
            char* label = edges[i].label;

            insEdgeAdjList(start, end, weight, label);
             }
    }


//================================================================
//                  Public methodes
//================================================================

void Graph::add( edge e) {
    embNode node(e.end,e.weight,e.label);
   (Adj[e.start].adjVec).push_back(node);
   M++;
    }



//================================================================
//                  Reading functions
//================================================================

int count_nb_lines( FILE* file )
{
  int number_of_lines = 0;
  char ch;

  
  while (true) {
        ch=fgetc(file);
        if (ch == (int)'\n') {
            number_of_lines++;
        } 
        if (ch == EOF){
            break;
        }
    }
  
  // Set the cursor back to the begining of the file.
  rewind(file);
  // Don't care if the last line has a '\n' or not. We over-estimate it.
  return number_of_lines; 
}

void read_node_file( FILE* node_file, vector<embNode>& seqs)
{
    seqs.clear();
    char* buffer = new char[100 * MAX];
    char* seq;

    seqs.reserve(count_nb_lines(node_file));
    int compt = 0;
    while ( fgets(buffer, 100 * MAX, node_file) != NULL )
    {
        char* p;

        if (strlen(buffer) == 100 * MAX)
        {  
          p = strtok(buffer, "\t\n");
          fprintf(stdout, "ERROR: node %s with sequence larger than %d!", p, 100 * MAX);
          exit(0);
        }
          
        // Node label
        p = strtok( buffer, "\t\n" );
              
        // Node seq
        p = strtok( NULL, "\t\n"  );
        seq = new char[strlen(p) + 1];
        strcpy( seq, p );
        embNode node(compt,0,seq);
        seqs.push_back(node);
        compt++;
    }

    delete [] buffer;
}

void read_edge_file( FILE *edge_file, vector<edge>& edges )
{
  char* buffer = new char[100 * MAX];
  char* u = new char[MAX];
  char* v = new char[MAX];

  edges.clear();
  edges.reserve(count_nb_lines(edge_file));
  while ( fgets(buffer, 100 * MAX, edge_file) != NULL )
  {
    char* p;

    // outgoing node
    p = strtok( buffer, "\t\n" );
    strcpy( u, p );

    // incoming node
    p = strtok( NULL, "\t\n" );
    strcpy( v, p );

    // edge label
    p = strtok( NULL, "\t\n" );
  
    edge e(atoi(u),atoi(v),0,p);
    edges.push_back(e);

  }
  
  
  delete [] buffer;
  delete [] u;
  delete [] v;
  
}




//================================================================
//                  Test part
//================================================================



// print all adjacent vertices of given vertex
void printGraphVertices(Graph& G)
{
    vector<embNode> V = G.Vertices;
    for (vector<embNode>::iterator it = V.begin(); it != V.end(); ++it) {
        cout << it->val << "\t" << (string)it->label << "\n";
    }
    cout << endl;
}

void printVertices(vector<embNode>& V)
{
    for (vector<embNode>::iterator it = V.begin(); it != V.end(); ++it) {
        cout << it->val << "\t" << (string)it->label << "\n";
    }
    cout << endl;
}

void printGraphEdges(Graph& G)
{
    vector<adjNode> Adj = G.Adj;
    int N = G.N;
    for (int i = 0; i<N ; i++ ) {
        for (vector<embNode>::iterator it = Adj[i].adjVec.begin(); it != Adj[i].adjVec.end(); ++it) {
            cout << i << " to " << it->val << " of type " << (string)it->label << " and of weight " << it->weight << endl;
        }
    }
    cout << endl;
}

void printEdges(vector<edge>& E)
{
    int M = E.size();
    edge* e;
    for (int i = 0; i<M ; i++ ) {
        *e = E[i];
        cout << e->start << " to " << e->end << " of type " << (string)e->label << " and of weight " << e->weight << endl;
    }
    cout << endl;
}

// graph implementation
int main()
{
    vector<edge> E;
    vector<embNode> V;
    FILE * edges;
    FILE * nodes;

    edges = fopen("../data/graph.edges","r");
    nodes= fopen("../data/graph.nodes","r");

    read_edge_file(edges,E);
    read_node_file(nodes,V);

    printVertices(V);
    printEdges(E);

    Graph G(V,E);

    printGraphVertices(G);
    printGraphEdges(G);

    cout << "Il y a " << G.M << " arrêtes" << endl;

    char l[2] = {'F','R'};
    edge e(0,1,2,l);
    G.add(e);
    cout << "Il y a " << G.M << " arrêtes" << endl;


    printGraphEdges(G);
    return 0;
}