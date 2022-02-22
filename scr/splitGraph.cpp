/*
 * Ce programme permet de couper un graphe en deux par rapport
 * à une lecture dans le sens foward et dans le sens reverse
 */

#include <list>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <cstring>
#include <algorithm>
#include "graph.h"
#define MAX 1024

//Lit un fichier contenant les sommets du graphe et les ajoute au vecteur seqs (réalisé par Pierre Peterlongo et Vincent Lacroix)
void read_node_file_4args( FILE* node_file, vector<Node>& seqs)
{
    seqs.clear();
    char* buffer = new char[100 * MAX];
    char* seq;
    char* u = new char[MAX];
    char* v = new char[MAX];

    seqs.reserve(count_nb_lines(node_file));

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
        strcpy( u, p );

        // Node seq
        p = strtok( NULL, "\t\n"  );
        seq = new char[strlen(p) + 1];
        strcpy( seq, p );

        // Node weight
        p = strtok( NULL, "\t\n" );
        strcpy( v, p );

        // Node weight
        p = strtok( NULL, "\t\n" );

        if (atoi(v)){
            Node node(atoi(u),atoi(v),seq);
            seqs.push_back(node);
        }
    }

    delete [] buffer;
    delete [] u;
    delete [] v;
}



// graph implementation
int main(int argc, char** argv)
{
    if (argc!=5){
        cout << "Expected use of this program: \n\n\t" <<argv[0] << " originalFile.nodes file.node file.edges outputPrefixe\n" << endl;
        return 0;
    }

    vector<Edge> E;
    vector<Edge> foward;
    vector<Edge> reverse;
    vector<Node> V;
    vector<Node> Vvu;

    char* nodesPath = argv[1];
    char* nodesVuPath = argv[2];
    char* edgesPath = argv[3];

    FILE * edges;
    FILE * nodes;
    FILE * nodesVu;

    edges = fopen(edgesPath,"r");
    nodes= fopen(nodesPath,"r");
    nodesVu= fopen(nodesVuPath,"r");

    read_edge_file(edges,E);
    read_node_file(nodes,V);
    read_node_file_4args(nodesVu,Vvu);

    Graph G(V,E);


    foward.clear();
    reverse.clear();

    if (Vvu.empty()){
        cout << "Pas de références!" << endl;
        return 0;
    }

    Node node = Vvu.front();
    cout << "Le noeud de départ est le " << node.val << " de séquence " << node.label << endl;
    vector<Neighbor*> aVoir;
    vector<int> vu;
    aVoir.clear();
    vu.clear();
    vu.push_back(node.val);
    for (vector<Neighbor>::iterator it = G.Neighbors(node.val)->begin(); it != G.Neighbors(node.val)->end(); ++it){
        if (it->label[0] == 'F'){
            aVoir.push_back(&(*it));
            Edge edge(node.val,it->val,0,it->label);
            foward.push_back(edge);
        }
    }
    G.BFS(1,foward,aVoir,vu);

    aVoir.clear();
    vu.clear();
    vu.push_back(node.val);
    for (vector<Neighbor>::iterator it = G.Neighbors(node.val)->begin(); it != G.Neighbors(node.val)->end(); ++it){
        if (it->label[0] == 'R'){
            aVoir.push_back(&(*it));
            char r[2];
            comp(it->label,r);
            Edge edge(it->val,node.val,0,r);
            reverse.push_back(edge);
        }
    }
    G.BFS(0,reverse,aVoir,vu);




    ofstream edges_foward;
    ofstream edges_reverse;
    edges_foward.open((string)argv[4]+(string)"_foward.edges");
    printEdges(foward,edges_foward);
    edges_foward.close();
    edges_reverse.open((string)argv[4]+(string)"_reverse.edges");
    printEdges(reverse,edges_reverse);
    edges_reverse.close();
    return 0;

}