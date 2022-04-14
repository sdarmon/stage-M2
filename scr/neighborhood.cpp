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
    int threshold,dis;
    if (argc == 4) {
        threshold = 0;
        dis = 1;
    } else if (argc == 6 and argv[4] == "-d"){
        int dis = atoi(argv[5]);
        threshold= 0;
    }else if (argc == 6 and argv[4] == "-c"){
        int threshold = atoi(argv[5]);
        dis = 1;
    }else if (argc == 8 and argv[4] == "-c"){
        int threshold = atoi(argv[5]);
        dis = atoi(argv[7]);
    }else if (argc == 8 and argv[4] == "-d"){
        int threshold = atoi(argv[7]);
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
    vector<Neighbor*> aVoir;
    aVoir.clear();
    vector<int> pronf;
    pronf.clear();
    Neighbor* node ;
    int p ;
    vector<int> vu(G.N,0);
    vu[node_id]=-1; //On le marque comme sommet d'entrée

    //On ajoute les voisins à visiter
    for (vector<Neighbor>::iterator it = G.Neighbors(node_id)->begin(); it != G.Neighbors(node_id)->end(); ++it) {
        if (G.Vertices[it->val].weight >= threshold) {
            aVoir.push_back(&(*it));
            pronf.push_back(1);
        }
    }

    //On fait un BFS
    while (aVoir.size() != 0){ //Cas de terminaison, on a terminé le BFS
        node = aVoir.front();
        aVoir.erase(aVoir.begin());
        p = pronf.front();
        pronf.erase(pronf.begin());
        if (vu[node->val]){ //Cas où le sommet a été vu par le BFS
            continue;
        }
        vu[node->val]=p;

        for (vector<Neighbor>::iterator it = G.Neighbors(node->val)->begin(); it != G.Neighbors(node->val)->end(); ++it){
            //On boucle sur ses voisins
            if (Vertices[it->val].weight >= threshold and vu[it->val]==0){
                //Cas où l'arrêt est bien valide et sommet non vu avant, ce voisin est rajouté dans la file des visites
                aVoir.push_back(&(*it));
                pronf.push_back(p+1);
            }
        }
    }

    for (int i = 0; i < G.N ; i++){
        if (vu[i]){
            cout << i << "\t" << G.Vertices[i].label << "\t" << ((vu[i] < 0) ? 0 : vu[i]) << "\t" << G.Vertices[i].weight << endl;
        }
    }

}