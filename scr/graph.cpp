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
#include <fstream>
#include <string>
#include <stdio.h>
#include <cstring>
#include <algorithm>

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

Node::Node(int v, int w, string l){
        val = v;
        weight = w;
        label=(string)l;
    }

Neighbor::Neighbor(int v, int w, char* l){
        val = v;
        weight = w;
        strcpy(label,l);
    }

// stores adjacency list items
LstNode::LstNode() {
        adjVec.clear();
    }
LstNode::LstNode(vector<Neighbor>& A) {
        adjVec = A;
    }


// structure to store edges
Edge::Edge(int s, int e, int w, char* l) {
        start = s;
        end = e;
        weight = w;
        strcpy(label,l);
    }

// ===========================================================================
//                                  Class
// ===========================================================================



Graph::Graph(vector<Node>& vertices, vector<Edge>& edges)  {

        N = vertices.size();
        M = 0;
        Vertices = vertices;
        Adj.clear();

        // initialize the LstNode for all vertices
        for (int i = 0; i < vertices.size(); i++){
            LstNode adj;
            Adj.push_back(adj);
            (Adj.back().adjVec).clear();
        }

        // construct directed graph by adding edges to it
        for (unsigned i = 0; i < edges.size(); i++)  {
            add(edges[i].start, edges[i].end, edges[i].weight, edges[i].label);
             }
    }


//================================================================
//                  Public methodes  
//================================================================

void Graph::add( Edge &e) {
    Neighbor node(e.end,e.weight,e.label);
   (Adj[e.start].adjVec).push_back(node);
   M++;
    }

void Graph::add(int start,int end, int weight, char* labelEdge)   {
        Neighbor node(end,weight,labelEdge);
       (Adj[start].adjVec).push_back(node);
       M++;
    }

void Graph::add(Node &v) {
   Vertices.push_back(v);
       LstNode A;
       Adj.push_back(A);
   N++;
    }

void Graph::add(int val,int weight, string label)   {
        Node node(val,weight,label);
       Vertices.push_back(node);
       LstNode A;
       Adj.push_back(A);
       N++;
    }

Neighbor* Graph::link(int start, int end) {
    if (start >= N) {
        return NULL;
    }
    for (int index = 0; index < Adj[start].adjVec.size(); index++){
        if (Adj[start].adjVec[index].val == end) {
            return &Adj[start].adjVec[index];
        }
    }
    return NULL;
}

vector<Neighbor>* Graph::Neighbors(int n){
    return &(Adj[n].adjVec);
}

void Graph::weighing(){
    for (int index = 0; index < N; index++){
        for (vector<Neighbor>::iterator it = Neighbors(index)->begin(); it != Neighbors(index)->end(); ++it){
            it->weight = Vertices[it->val].label.size();
        }
    }
}




// /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ 
//   Petite approximation ici, je suppose que dès lors que l'on 
//   voit un sommet par le sens foward, on ne le recroisera pas
//   par le sens reverse! Cela semble peu probable d'arriver et
//   l'impact semble peu important (car BFS) mais on sous-estime
//   la vraie valeur.
// /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ 


int Graph::BFSCount(vector<int> &rayons, int acc,vector<Neighbor*> &aVoir,vector<int> &vu){
    if (aVoir.size() == 0){
        return acc;
    }
    int rayon = rayons.front();
    rayons.erase(rayons.begin());
    Neighbor* node = aVoir.front();
    aVoir.erase(aVoir.begin());
    if (find(vu.begin(),vu.end(),node->val) != vu.end()){
        //Cas où le sommet a été vu par le BFS
        return BFSCount(rayons,acc,aVoir,vu);
    }
    vu.push_back(node->val);
    if (node->weight - 40 <= rayon){
        //Cas où le sommet est bien dans la boule de rayon
        for (vector<Neighbor>::iterator it = Neighbors(node->val)->begin(); it != Neighbors(node->val)->end(); ++it){
            //On boucle sur ses voisins
            if (it->label[0] == node->label[1] && find(vu.begin(),vu.end(),it->val) == vu.end()){ 
                //Cas où l'arrêt est bien valide et sommet non vu avant, ce voisin est rajouté dans la file des visites
                aVoir.push_back(&(*it));
                rayons.push_back(rayon-node->weight+40);
            }
        }
        return BFSCount(rayons,acc+1,aVoir,vu); //On traite les cas suivants, en prennant en compte le sommet
    }

    return BFSCount(rayons,acc,aVoir,vu); //On continue sans prendre en compte le sommet.
}

void Graph::weighingANode(int source, int rayon) {
    vector<Neighbor*> aVoir;
    vector<int> vu;
    vector<int> rayons;
    rayons.clear();
    aVoir.clear();
    vu.clear();
    vu.push_back(source);
    for (vector<Neighbor>::iterator it = Neighbors(source)->begin(); it != Neighbors(source)->end(); ++it){
        aVoir.push_back(&(*it));
        rayons.push_back(rayon);
    }
    Vertices[source].weight = BFSCount(rayons,1,aVoir,vu);
}

void Graph::weighingAllNodes(int rayon) {
    for (vector<Node>::iterator it = Vertices.begin(); it != Vertices.end(); ++it){
        weighingANode(it->val, rayon);
    }
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

void read_node_file( FILE* node_file, vector<Node>& seqs)
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
        Node node(compt,0,seq);
        seqs.push_back(node);
        compt++;
    }

    delete [] buffer;
}

void read_edge_file( FILE *edge_file, vector<Edge>& edges )
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

    // Edge label
    p = strtok( NULL, "\t\n" );
  
    Edge e(atoi(u),atoi(v),0,p);
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
    vector<Node> V = G.Vertices;
    for (vector<Node>::iterator it = V.begin(); it != V.end(); ++it) {
        cout << it->val << "\t" << (string)it->label << "\t" << it->weight << "\n";
    }
    cout << endl;
}

void printGraphVertices(Graph& G,ofstream& output)
{
    vector<Node> V = G.Vertices;
    for (vector<Node>::iterator it = V.begin(); it != V.end(); ++it) {
        output << it->val << "\t" <<(string)it->label << "\t" << it->weight << "\n";
    }
}


void printVertices(vector<Node>& V)
{
    for (vector<Node>::iterator it = V.begin(); it != V.end(); ++it) {
        cout << it->val << "\t" << (string)it->label << "\t" << it->weight << "\n";
    }
    cout << endl;
}

void printGraphEdges(Graph& G)
{
    vector<LstNode> Adj = G.Adj;
    int N = G.N;
    for (int i = 0; i<N ; i++ ) {
        for (vector<Neighbor>::iterator it = Adj[i].adjVec.begin(); it != Adj[i].adjVec.end(); ++it) {
            cout << i << " to " << it->val << " of type " << (string)it->label << " and of weight " << it->weight << endl;
        }
    }
    cout << endl;
}

void printEdges(vector<Edge>& E)
{
    int M = E.size();
    Edge* e;
    for (int i = 0; i<M ; i++ ) {
        *e = E[i];
        cout << e->start << " to " << e->end << " of type " << (string)e->label << " and of weight " << e->weight << endl;
    }
    cout << endl;
}




// graph implementation
int main(int argc, char** argv)
{
    if (argc!=3 and argc!=5){
        cout << "Expected use of this program: \n\n\t" <<argv[0] << " file.nodes file.edges -o output.txt\n" << endl;
        return 0;
    }

    vector<Edge> E;
    vector<Node> V;

    char* nodesPath = argv[1];
    char* edgesPath = argv[2];


    FILE * edges;
    FILE * nodes;

    edges = fopen(edgesPath,"r");
    nodes= fopen(nodesPath,"r");

    read_edge_file(edges,E);
    read_node_file(nodes,V);

    Graph G(V,E);

    G.weighing();
    G.weighingAllNodes(200);

    if(argc == 5){
        ofstream output;
        output.open(argv[4]);
        printGraphVertices(G,output);
        output.close();
    } else{
        printGraphVertices(G);
    }
    return 0;

}

    // printVertices(V);
    // printEdges(E);

    // printGraphVertices(G);
    // printGraphEdges(G);

    // cout << "Il y a " << G.M << " arrêtes" << endl;

    // char l[2] = {'F','R'};
    // Edge e(0,1,2,l);
    // G.add(e);
    // cout << "Il y a " << G.M << " arrêtes" << endl;


    // printGraphEdges(G);


    // Neighbor* p;
    // p =G.link(2,6);
    // if (p != NULL){
    //     cout << "Edge found! weight: " << p->weight << "  type: " << (string)p->label << endl;
    // } else {
    //     cout << "Edge not found..." << endl;
    // }

    // p = G.link(20,3);
    // if (p != NULL){
    //     cout << "Edge found! weight: " << p->weight << "  type: " << (string)p->label << endl;
    // } else {
    //     cout << "Edge not found..." << endl;
    // }

    // p = G.link(1,6);
    // if (p != NULL){
    //     cout << "Edge found! weight: " << p->weight << "  type: " << (string)p->label << endl;
    // } else {
    //     cout << "Edge not found..." << endl;
    // }



//We do not allow to remove anything for the moment, maybe it will be mandatory for some future implementations....

// bool Graph::remove(int valNode){
//     for (int index = 0, index < N, index++){
//         if (Vertices[index].val == valNode){
//             pop(index);
//             return 1;
//         }
//         return 0;
//     }  
// }

// bool Graph::remove(Node n){
//     for (int index = 0, index < N, index++){
//         if (Vertices[index] == n){
//             pop(index);
//             return 1;
//         }
//         return 0;
//     }  
// }

// bool Graph::remove(int startEdge, int endEdge){
//     for (int index = 0, index < Adj[startEdge].size(), index++){
//         if (Vertices[index] == n){
//             pop(index);
//             return 1;
//         }
//         return 0;
//     }  
// }

// bool Graph::remove(Edge e){
//     int start, end;
//     for (int index = 0, index < N, index++){
//         if (Vertices[index] == n){
//             pop(index);
//             return 1;
//         }
//         return 0;
//     }  
// }

// Node Graph::pop(int index){
//     int e;
//     Node n(Vertices[index]);
//     Vertices.erase(Vertices.begin()+index);
//     e = Adj[index].size();
//     Adj.erase(Adj.begin()+index);
//     N--;
//     M= M-e;
//     return n;
// }