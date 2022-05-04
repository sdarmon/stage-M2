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
#include <queue>
#include <map>
#include "graph.h"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wstring-compare"
#define MAX 1024


//Fonction permettant de récupérer l'index du sommet de `G` ayant le poids le plus élevé et qui n'a pas été déjà vu dans
//`vu_total` (i.e. `vu_total[index] == 0`)
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

//Initialisation d'un vecteur avec la valeur -1. Il me semble que vector<int> vec(-1,n) doit fonctionner, à corriger!!
void initVec(vector<int> &vec, int n){
    vec.clear();
    for(int i = 0; i<n; i++){
        vec.push_back(-1);
    }
    return;
}

//Fonction permettant d'enregistrer une composante du graphe `G` dans un fichier
void save_comp(Graph &G, vector<int> &compo, string outputPrefix, int rang){
    ofstream output;
    output.open(outputPrefix+"/processing/comp"+to_string(rang)+".txt");
    for (vector<int>::iterator it = compo.begin(); it != compo.end(); it++){
        output << *it << "\t" << G.Vertices[*it].label << "\t" << G.Vertices[*it].weight  << "\n";
    }
    return;
}

//Fonction indiquant si le setA est un sous-ensemble du setB. Le setA et le setB doivent être triés.
int subset(vector<int> &setA, vector<int> &setB){
    if (setA.size() > setB.size()){
        return 0;
    }
    int posi = 0;
    int posj = 0;
    while(posi< setA.size() and posj < setB.size()) {
        if (setA[posi] == setB[posj]){
            posi++;
            posj++;
        } else if (setA[posi] < setB[posj]){ //Cas où le i-ème élément n'a pas été trouvé dans le setB
            return 0;
        }
        else {
            posj++;
        }
    }
    return posi == setA.size();
}


int main(int argc, char** argv) {
    if (argc != 8 and argc != 6) {
        cout << "Expected use of this program: \n\n\t" << argv[0]
             << " file.nodes file.edges -c value [-k kmer] outputPrefix \n" << endl;
        return 0;
    }

    //On charge le graphe
    vector <Edge> E;
    vector <Node> V;

    char *nodesPath = argv[1];
    char *edgesPath = argv[2];
    string outputPrefix;

    ifstream edges(edgesPath, std::ios::binary);
    ifstream nodes(nodesPath, std::ios::binary);

    read_edge_file(edges, E);
    read_node_file_weighted(nodes, V);

    Graph G(V, E);
    if(argc >= 7 and argv[5][1]=='k' ){
        G.kmer = stoi(argv[6]);
        outputPrefix = argv[7];
    } else {
        outputPrefix = argv[5];
    }

    cout << "Graphe chargé et construit" << endl;

    int index;
    //Construction du vector `vu_total`
    vector<int> vu_total(G.N, 0);

    vector <vector<int>> components;
    components.clear();

    //Ici le critère choisi pour définir une composante est juste d'être le sous-graphe connexe maximum du graphe G
    //auquel on a enlevé tous les sommets de poids inférieur au threshold `threshold`.
    int threshold;
    threshold = atoi(argv[4]);
    int m = 0;
    int compt;
    vector<int> vu(G.N, 0);
    queue < Neighbor * > aVoir;
    vector<int> indexTrie;
    indexTrie.clear();
    for (int sommet = 0; sommet < G.N ; sommet++){
        if (G.Vertices[sommet].weight >= threshold){
            indexTrie.push_back(sommet);
        }
    }
    cout << "Début du tri des " << indexTrie.size() << " sommets en fonction de leur poids" << endl;
    sort(indexTrie.begin(),indexTrie.end(),G);
    index = indexTrie.back();
    indexTrie.pop_back();
    vu_total[index] = 1;

    cout << "Début de la recherche des composantes" << endl;
    while (!indexTrie.empty()) //On cherche des composantes tant qu'il existe encore un sommet vérifiant
        //le critère et dans la limite des 100 composantes.
        //Condition m < 100 enlevée
    {
        //On réalise alors un BFS depuis le sommet `index`
        for (int i = 0; i < G.N; i++) {
            vu[i] = 0;
        }
        vector<int> compo;
        compo.clear();
        //aVoir étant une queue de BFS; est toujours censé être vide ici. Pour la sanité des algos qui suivent on la
        //vide au cas où...
        while (!aVoir.empty()){
            aVoir.pop();
        }
        vu[index] = 1;
        for (vector<Neighbor>::iterator it = G.Neighbors(index)->begin(); it != G.Neighbors(index)->end(); ++it) {
            if (G.Vertices[it->val].weight >= threshold) {
                aVoir.push(&(*it));
            }
        }
        G.BFS_func(threshold, aVoir, vu);

        //Puis on ajoute les sommets trouvés à la composante
        compt=0;
        for (int i = 0; i < G.N; i++) {
            if (vu[i]) {
                compt++;
                compo.push_back(i);
                vu_total[i] = 1;
            }
        }
        if (compt > 1){
            save_comp(G, compo, outputPrefix, m); //On enregistre la composante sur l'ordinateur
            components.push_back(compo); //Et on ajoute la composante au vecteur de composantes
            cout << "Composante trouvée de départ " << index << " et de poids " << G.Vertices[index].weight
                 << " et de taille " << compo.size() << endl;
        }

        //Finalement, on recommence la boucle while
        m++;
        while (!indexTrie.empty() and vu_total[indexTrie.back()]){
            indexTrie.pop_back();
        }
        if (!indexTrie.empty()){
            index = indexTrie.back();
            indexTrie.pop_back();
            vu_total[index] = 1;
        }

    }

    cout << "Fin de la recherche de composantes." << endl;
}