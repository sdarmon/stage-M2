/*
 * Ce programme permet de faire une agglomeration des composantes
 * des graphes.
 */

#include <list>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <set>
#include <iterator>
#include "graph.h"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wstring-compare"
#define MAX 1024


//Fonction permettant d'enregistrer une composante du graphe `G` dans un fichier
void save_comp(Graph &G, vector<int> &compo, string outputPrefix, int rang){
    ofstream output;
    output.open(outputPrefix+"/processing/comp"+to_string(rang)+".txt");
    for (vector<int>::iterator it = compo.begin(); it != compo.end(); it++){
        output << *it << "\t" << G.Vertices[*it].label << "\t" << G.Vertices[*it].weight  << "\n";
    }
    return;
}

int main(int argc, char** argv) {
    int c,dis;
    if (argc == 4) {
        c = 0;
        dis = 1;
    } else if (argc == 6 and argv[4] == "-d"){
        int dis = atoi(argv[5]);
        c= 0;
    }else if (argc == 6 and argv[4] == "-c"){
        int c = atoi(argv[5]);
        dis = 1;
    }else if (argc == 8 and argv[4] == "-c"){
        int c = atoi(argv[5]);
        dis = atoi(argv[7]);
    }else if (argc == 8 and argv[4] == "-d"){
        int c = atoi(argv[7]);
        dis = atoi(argv[5]);
    } else {
        cout << "Expected use of this program: \n\n\t" << argv[0]
             << " file.nodes file.edges node_id [-c value] [-d dis] \n"
             << "The option c allows to give a threshold to the weight of visited vertices (default 0) and the option d gives the "
                "depth of neighbor gathered (default 1)" << endl;
        return 0;
    }

    //On charge le graphe
    vector <Edge> E;
    vector <Node> V;

    char *nodesPath = argv[1];
    char *edgesPath = argv[2];
    int node_id = argv[3];
    ifstream edges(edgesPath, std::ios::binary);
    ifstream nodes(nodesPath, std::ios::binary);

    read_edge_file(edges, E);
    read_node_file_weighted(nodes, V);

    Graph G(V, E);

    cout << "Graphe chargé et construit" << endl;

}