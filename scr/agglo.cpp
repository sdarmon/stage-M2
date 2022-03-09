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
void chemin(int i, vector<int> &endings, Graph &G, vector<vector<int>> &components, int profondeurMax){
    int step = 1;
    vector<int> comp1;
    vector<Neighbor*> voisin1;
    voisin1.clear();
    comp1.clear();
    vector<vector<int>> vu;
    vu.clear();

    for (int j = 0; j<components.size(); j++){
        vector<int> vu_row(G.N,0);
        vu.push_back(vu_row);
    }

    for(int k= 0; k< components[i].size(); k++){
        vu[i][components[i][k]] = 1;
        comp1.push_back(components[i][k]);
        //cout << components[i][k] << endl;
        for (vector<Neighbor>::iterator it = G.Neighbors(components[i][k])->begin(); it != G.Neighbors(components[i][k])->end(); ++it){
            voisin1.push_back(&(*it));
        }
    }
    //On regarde si la j-eme composante s'intersecte avec la comp1
    for (int j = 0; j<components.size(); j++){
        for(int k= 0; k< components[j].size(); k++){
            vu[j][components[j][k]]=1;
            if (vu[i][components[j][k]]){ 
                endings[j] = 0;
                continue;
            }
        }
    }
    cout << "\tData initialisée" << endl;
    int modif = 1;
    int trouve = 0;
    int size;
    Neighbor* node;

    //cout << "Initialisation chemin ok" << endl;
    while (modif && step < profondeurMax) { //Tant qu'un sommet a été rajouté à la couche précédente, on regarde tous les sommets de cette couche là dans le BFS.
        modif = 0;
        trouve = 0;
        size = voisin1.size();
        for (int k = 0; k<size; ++k){ //On boucle sur les sommets de la couche du BFS uniquement

            node = voisin1.front();
            voisin1.erase(voisin1.begin());
            if (vu[i][node->val]){ // Sommet déjà présent
                continue;
            }
            // Sommet également présent dans la composante d'arrivée
            for (int j = 0; j<components.size(); j++){
                if (vu[j][node->val]){ 
                    if (endings[j] == -1){ //On garde toujours le plus court chemin
                        endings[j] = step;
                    }
                    vu[i][node->val] = 1;
                    trouve = 1;
                    continue;
                }
            }
            if (not trouve) {
                //Sinon on rajoute ce sommet comme vu et on ajoute ses voisins à traiter dans le BFS pour la prochaine couche
                comp1.push_back(node->val);
                vu[i][node->val] = 1;
                modif = 1;
                for (vector<Neighbor>::iterator it = G.Neighbors(node->val)->begin(); it != G.Neighbors(node->val)->end(); ++it){ 
                    if (it->label[0] != node->label[1]){ 
                        continue;
                    }
                    voisin1.push_back(&(*it));
                }
            }
        }
        step ++; //On passe à la couche suivante
    }
    //Aucun nouveau sommet n'a été rajouté à traiter... Il n'y a plus de chemins à visiter...
    return;
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

void initVec(vector<int> &vec, int n){
    vec.clear();
    for(int i = 0; i<n; i++){
        vec.push_back(-1);
    }
    return;
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
    int sizeMax;
    vector<int> vu(G.N,0);
    sizeMax = 2000000;
    vector<Neighbor*> aVoir;

    index = indexMax(G,vu_total);
    val = G.Vertices[index].weight;
    //cout << index << " de poids " << G.Vertices[index].weight << endl;
    vu_total[index]= 1;

    cout << "Début de la recherche des composantes" << endl;
    while (val >= threshold and m < 50) //On cherche les composantes
    {
        for (int i = 0; i<G.N; i++){
            vu[i] = 0;
        }
        vector<int> compo;
        compo.clear();
        aVoir.clear();
        vu[index] = 1;
        for (vector<Neighbor>::iterator it = G.Neighbors(index)->begin(); it != G.Neighbors(index)->end(); ++it){
            aVoir.push_back(&(*it));
        }
        m++;
        G.BFS_func(threshold,sizeMax, aVoir,vu); 


        for (int i = 0; i<G.N; i++){
            if (vu[i]){
                compo.push_back(i);
            }
        }
        components.push_back(compo);
        for (int i = 1; i< compo.size(); i++){
            vu_total[compo[i]] = 1;
        }

        cout << "Composante trouvée de départ " << index << " et de poids " << G.Vertices[index].weight << " et de taille " << compo.size() << endl;
        index = indexMax(G,vu_total);
        val = G.Vertices[index].weight;
        vu_total[index]= 1;

    }

    cout << "Fin de la recherche de composantes." << endl;

    //Maintenant on construit nouveau graphe contracté
    vector<Edge> E2;
    vector<Node> V2;
    vector<int> endings;
    E2.clear();
    V2.clear();
    int c;
    for (int i = 0; i < components.size(); i++){
        cout << "Calcul des chemins partant de " << i << endl;
        initVec(endings, components.size());
        chemin(i,endings,G,components,5);

        V2.push_back(Node(i,components[i].size(),""));

        for (int j = 0; j< components.size(); j++){
            if (endings[j] >= 0){             
                char arete[2] = {'N','N'};
                E2.push_back(Edge(i,j,endings[j],arete));
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