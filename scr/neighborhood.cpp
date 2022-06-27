/*
 *
 */

#include <list>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <set>
#include "graph.h"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wstring-compare"
#define MAX 1024


//Fonction permettant d'enregistrer une composante du graphe `G` dans un fichier
void save_comp(Graph &G, vector<int> &vu, string outputPrefix, int rang){
    ofstream output;
    output.open(outputPrefix+".nodes");
    for (int i = 0; i < G.N ; i++){
        if (vu[i]){
            output << i << "\t" << G.Vertices[i].label << "\t" << ((vu[i] < 0) ? 0 : vu[i]) << "\t" << G.Vertices[i].weight << endl;
        }
    }
    return;
}

int main(int argc, char** argv) {
    int threshold,dis;
    string prefixOutput;
    if (argc == 4) {
        prefixOutput = "./";
        threshold = 0;
        dis = 1;
    }else if (argc == 6 and (string)argv[4] == "-o"){
        prefixOutput = argv[5];
        threshold = 0;
        dis = 1;
    }else if (argc == 6 and (string)argv[4] == "-c"){
        prefixOutput = "./";
        threshold = atoi(argv[5]);
        dis = 1;
    }else if (argc == 6 and (string)argv[4] == "-d"){
        prefixOutput = "./";
        threshold = 0;
        dis = atoi(argv[5]);
    }else if (argc == 8 and (string)argv[4] == "-o" and (string)argv[6] == "-d"){
        prefixOutput = argv[5];
        dis = atoi(argv[7]);
        threshold = 0;
    }else if (argc == 8 and (string)argv[4] == "-c" and (string)argv[6] == "-d"){
        threshold = atoi(argv[5]);
        dis = atoi(argv[7]);
        prefixOutput = "./";
    }else if (argc == 10 and (string)argv[4] == "-o" and (string)argv[6] == "-c" and (string)argv[8] == "-d"){
        prefixOutput = argv[5];
        threshold = atoi(argv[7]);
        dis = atoi(argv[9]);
    } else {
        cout << "Expected use of this program: \n\n\t" << argv[0]
             << " file.nodes file.edges node_id/input.txt [-o prefixOutput] [-c value] [-d dis] \n"
             << "The option c allows to give a threshold to the weight of visited vertices (default 0) and the option d gives the "
                "depth of neighbor gathered (default 1)" << endl;
        cout << argc << " arguments have been used " ;
        for(int i =0; i<argc ; i++){
            cout << "arg" << i << ": " << argv[i] << "\t" ;
        }
        cout << endl;
        return 0;
    }

    //On charge le graphe
    vector <Edge> E;
    vector <Node> V;

    char *nodesPath = argv[1];
    char *edgesPath = argv[2];
    vector<int> nodes_id;
    nodes_id.clear();
    string arg = argv[3];
    if (arg.find('.') != string::npos){
        char *idPath =argv[3];
        ifstream ids(idPath, std::ios::binary);
        string line;
        while (getline(ids, line)){
            nodes_id.push_back(stoi(line));
        }
        ids.close();
    } else {
        nodes_id.push_back(atoi(argv[3]));
    }
    ifstream edges(edgesPath, std::ios::binary);
    ifstream nodes(nodesPath, std::ios::binary);

    read_edge_file(edges, E);
    read_node_file_weighted(nodes, V);

    Graph G(V, E);
    vector<Neighbor*> aVoir;
    vector<Edge> E2;
    vector<int> pronf;
    Neighbor* node ;
    int node_id;
    int p ;
    vector<int> vu;
    int taille = nodes_id.size();
    for(int rang = 0; rang <taille; rang++) {

        aVoir.clear();
        E2.clear();
        pronf.clear();
        for (int i = 0; i<G.N; i++){
            vu.push_back(0);
        }
        node_id = nodes_id.front();
        nodes_id.erase(nodes_id.begin());
        vu[node_id] = -1; //On le marque comme sommet d'entrée

        //On ajoute les voisins à visiter
        for (vector<Neighbor>::iterator it = G.Neighbors(node_id)->begin(); it != G.Neighbors(node_id)->end(); ++it) {
            if (G.Vertices[it->val].weight >= threshold) {
                aVoir.push_back(&(*it));
                pronf.push_back(1);
                E2.push_back(Edge(node_id,it->val,0,it->label));
            }
        }

        //On fait un BFS
        while (aVoir.size() != 0) { //Cas de terminaison, on a terminé le BFS
            node = aVoir.front();
            aVoir.erase(aVoir.begin());
            p = pronf.front();
            pronf.erase(pronf.begin());
            if (vu[node->val] != 0) { //Cas où le sommet a été vu par le BFS
                continue;
            }
            vu[node->val] = p;

            for (vector<Neighbor>::iterator it = G.Neighbors(node->val)->begin();
                 it != G.Neighbors(node->val)->end(); ++it) {
                //On boucle sur ses voisins
                if (G.Vertices[it->val].weight >= threshold and vu[it->val] == 0 and p <= dis) {
                    //Cas où l'arrêt est bien valide et sommet non vu avant, ce voisin est rajouté dans la file des visites
                    aVoir.push_back(&(*it));
                    pronf.push_back(p + 1);
                    E2.push_back(Edge(node->val,it->val,0,it->label));
                }
            }
        }
        save_comp(G, vu, prefixOutput, rang);
        ofstream edges;
        edges.open(prefixOutput+".edges");
        printEdges(E2,edges);
    }
}