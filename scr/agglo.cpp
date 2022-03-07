/*
 * Ce programme permet de faire une aglomération des composantes
 * des graphes.
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

// Permet de déterminer s'il existe un chemin de la composante i à la composante j, et s'il existe, renvoie la distance.
int chemin(int i, int j, Graph &G, vector<vector<int>> &components){
    int step = 1;
    vector<int> comp1;
    vector<int> comp2;
    vector<Neighbor*> voisin1;
    voisin1.clear();
    comp1.clear();
    comp2.clear();

    for(int k= 0; k< components[i].size(); k++){
        comp1.push_back(components[i][k]);
        for (vector<Neighbor>::iterator it = G.Neighbors(components[i][k])->begin(); it != G.Neighbors(components[i][k])->end(); ++it){
            if (find(comp1.begin(),comp1.end(), it->val) != comp1.end()){ 
                continue;
            }
            voisin1.push_back(&(*it));
        }
    }
    for(int k= 0; k< components[j].size(); k++){
        if (find(comp1.begin(),comp1.end(),components[j][k]) != comp1.end()){ 
            return 0;
        }
        comp2.push_back(components[j][k]);
    }

    int modif = 1;
    int size = voisin1.size();
    Neighbor* node;
    while (modif ) { //Tant qu'un sommet a été rajouté à la couche précédente, on regarde tous les sommets de cette couche là dans le BFS.
        modif = 0;
        size = voisin1.size();
        for (int k = 0; k<size; ++k){ //On boucle sur les sommets de la couche du BFS uniquement
            node = voisin1.front();
            voisin1.erase(voisin1.begin());
            if (find(comp1.begin(),comp1.end(), node->val) != comp1.end()){ // Sommet déjà présent
                continue;
            }
            if (find(comp2.begin(),comp2.end(),node->val) != comp2.end()){ // Sommet également présent dans la composante d'arrivée
                return step;
            }
            //Sinon on rajoute ce sommet comme vu et on ajoute ses voisins à traiter dans le BFS pour la prochaine couche
            comp1.push_back(node->val);
            modif = 1;
            for (vector<Neighbor>::iterator it = G.Neighbors(node->val)->begin(); it != G.Neighbors(node->val)->end(); ++it){ 
                if (it->label[0] != node->label[1]){ 
                    continue;
                }
                voisin1.push_back(&(*it));
            }
        }
        step ++; //On passe à la couche suivante
    }
    //Aucun nouveau sommet n'a été rajouté à traiter... Il n'y a pas de chemin de la composante i à j.
    return 0;
}

int indexMax(Graph &G, vector<int> &vu_total){
    int index;
    int val = -1;
    for (int i = 0; i<G.N; i++){
        if (vu_total[i] == 0 && G.Vertices[i].weight > val){
            index = i;
            val = G.Vertices[i].weight;
        }
    }
    return index;
}
int main(int argc, char** argv){
    if (argc!=6){
        cout << "Expected use of this program: \n\n\t" <<argv[0] << " file.node file.edges -c value outputPrefixe\n" << endl;
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
    read_node_file_weighted(nodes,V);

    Graph G(V,E);

    cout << "Graphe chargé et construit" << endl;

    int index;
    vector<int> vu_total;
    for (int i=0; i<G.N ; i++){
        vu_total.push_back(0);
    }

    vector<vector<int>> components;
    components.clear();

    int threshold;
    threshold =atoi(argv[4]);
    int m = 0;
    int val;
    vector<Neighbor*> aVoir;

    index = indexMax(G,vu_total);
    //cout << index << " de poids " << G.Vertices[index].weight << endl;
    vu_total[index]= 1;
    while (G.Vertices[index].weight >= threshold and m < 30) //On cherche les composantes
    {
        vector<int> compo;
        compo.clear();
        aVoir.clear();
        aVoir.push_back(val);
        m++;
        G.BFS_func(threshold ,aVoir,compo);

        components.push_back(compo);
        for (int i = 1; i< compo.size(); i++){
            vu_total[compo[i]] = 1;
        }
        index = indexMax(G,vu_total);
        vu_total[index]= 1;
        //cout << index << " de poids " << G.Vertices[index].weight << endl;

    }

    cout << "Fin de la recherche de composantes." << endl;

    //Maintenant on construit nouveau graphe contracté
    vector<Edge> E2;
    vector<Node> V2;
    E2.clear();
    V2.clear();
    int c;
    for (int i = 0; i < components.size(); i++){
        V2.push_back(Node(i,components[i].size(),""));
        for (int j = 0; j< components.size(); j++){
            if (i == j){
                continue;
            }
            c = chemin(i,j,G,components);
            if (c>=0) {                
                char arete[2] = {'N','N'};
                E2.push_back(Edge(i,j,c,arete));
            }
        }
    }

    Graph H(V2,E2);

    cout << "Graphe aggloméré construit" << endl;

    ofstream outputNodes;
    outputNodes.open((string) argv[5]+".nodes");
    printGraphVertices(H,outputNodes);
    outputNodes.close();

    ofstream outputEdges;
    outputEdges.open( (string) argv[5]+".edges");
    printEdges(E2,outputEdges);
    outputEdges.close();
}