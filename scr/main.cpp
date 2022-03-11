/*
 * Ce programme permet de loader un graphe de De Bruijn
 * et de calculer les poids des arêtes et des sommets à
 * un rayon donné en entrée.
 */

#include <list>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <algorithm>
#include "graph.h"


// graph implementation
int main(int argc, char** argv)
{
    if (argc!=4 and argc!=6){
        cout << "Expected use of this program: \n\n\t" <<argv[0] << " file.nodes file.edges radius -o output.txt\n" << endl;
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
    G.weighingAllNodes(atoi(argv[3]));

    if(argc == 6){
        ofstream output;
        output.open(argv[5]);
        printGraphVertices(G,output);
        output.close();
    } else{
        printGraphVertices(G);
    }
    return 0;

}