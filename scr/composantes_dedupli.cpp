/*
 * Ce programme permet de faire une agglomeration des composantes
 * des graphes.
 */

#include <list>
#include <vector>
#include <iostream>
#include <fstream>
//#include <string>
#include <set>
#include <map>
#include <queue>
#include "graph.h"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wstring-compare"
#define MAX 1024

//Initialisation d'un vecteur avec la valeur -1. Il me semble que vector<int> vec(-1,n) doit fonctionner, à corriger!!
void initVec(vector<int> &vec, int n){
    vec.clear();
    for(int i = 0; i<n; i++){
        vec.push_back(-1);
    }
    return;
}

bool foundedEdge(int i ,int j,vector <Edge> &E2,int limite){
    //Ici il me semble que l'on peut faire une recherche dichotomique mais je n'en suis pas sûr à 100%
    for (int k = 0; k< limite; k++){
        if (E2[k].start == i and E2[k].end == j) {
            return true;
        }
    }
    return false;
}

int reverse(int n, int N){
    if (n>= N/2){
        return n-N/2;
    } else {
        return n+N/2;
    }
}

struct indexDic {
    int key;
    int weight;
    friend bool operator< (indexDic const& lhs, indexDic const& rhs) {
        return (lhs.weight < rhs.weight);
    }
};

int main(int argc, char** argv) {
    if (argc != 8) {
        cout << "Expected use of this program: \n\n\t" << argv[0]
             << " file.nodes file.edges -k kmer compoPrefix nbComp outputPrefix \n" << endl;
        return 0;
    }

    //On charge le graphe
    vector <Edge> E;
    vector <Node> V;

    string nodesPath = argv[1];
    string edgesPath = argv[2];
    string outputPrefix;
    string compoPrefix;
    int nbComp;
    ifstream edges(edgesPath, std::ios::binary);
    ifstream nodes(nodesPath, std::ios::binary);

    read_edge_file(edges, E);
    read_node_file_weighted(nodes, V);

    Graph G(V, E);
    G.kmer = stoi(argv[4]);
    compoPrefix = argv[5];
    nbComp = stoi(argv[6]);
    outputPrefix = argv[7];


    cout << "Graphe chargé et construit" << endl;

    cout << "Début construction graphe aggloméré" << endl;
    vector <Edge> E3;
    vector <Node> V3;
    E3.clear();
    V3.clear();
    vector<int> seen(G.N, 0);
    vector<int> correspondingVertex(G.N, 0);
    //seen est le tableau de correspondance entre les anciens sommets et les nouveaux.


    //Après, on fait un BFS à partir de chaque sommet dans la composante afin d'obtenir les voisins du périmètre.
    set<int> vu2;
    set<int> sonSet;
    typedef map<pair<int,int>,string> dicChem;
    dicChem areteFF;
    int limiteAretes;
    vector<int> sons;
    vector<int> depthSons;
    queue<int> depth;
    vector<string> labelSons;
    queue<string> labels;
    vector<Neighbor*> aretes;
    queue < Neighbor * > aVoir;
    int index = 0;
    int compteurDeBoucle = 0;
    char aretFF[3] = {'F', 'F'};
    // char aretRR[3] = {'R', 'R'};
    vector <int> comp;
    string line;
    string file_name;
    int compt;
    //On boucle sur les composantes
    cout << "Début des calculs des arêtes des comp " << nbComp<<endl;
    for (int component = 0; component<nbComp; component++){
        comp.clear();
        file_name = compoPrefix+to_string(component)+".txt";
        ifstream file(file_name, std::ios::binary);
        while (getline(file, line)){
            compt=stoi(line.substr(0,line.find('\t')));
            comp.push_back(compt);
            seen[compt] = component + 1;
        }
        areteFF.clear();
        sonSet.clear();
        //On boucle sur les sommets de la composante
        cout << "Pré-calcul pour la composante " << compteurDeBoucle << " de taille " << comp.size() << endl;
        compteurDeBoucle++;
        for (vector<int>::iterator it = comp.begin(); it != comp.end(); ++it) {
            if (G.Vertices[(*it)].label.size() < G.kmer) {
                continue;
            }
            //Cas sommet en périphérie
            sons.clear();
            aretes.clear();
            depthSons.clear();
            labelSons.clear();
            while (!aVoir.empty()) { // N'est jamais censé arriver !
                aVoir.pop();
            }
            while (!depth.empty()) { // N'est jamais censé arriver !
                depth.pop();
            }
            while (!labels.empty()) { // N'est jamais censé arriver !
                labels.pop();
            }
            vu2.clear();
            vu2.insert(reverse(*it, G.N)); //On voit bien le sommet duquel on part

            //On fait un BFS à partir de chaque sommet afin de savoir quels sommets du périmètre sont
            //atteignables à partir de chaque sommet de la composante
            //Cas en partant du forward
            for (vector<Neighbor>::iterator voisin = G.Neighbors(reverse(*it, G.N))->begin();
                 voisin != G.Neighbors(reverse(*it, G.N))->end(); ++voisin) {
                aVoir.push(&(*voisin));
                depth.push(1);
                string aux = G.Vertices[reverse(*it, G.N)].label;
                labels.push(aux);
            }
            G.BFS_comp(seen, vu2, aVoir, depth, labels, sons, depthSons, aretes, labelSons);
            for (int i = 0; i < sons.size(); i++) {
                sonSet.insert(reverse(sons[i], G.N)); //Cette opération permet de récupérer les sommets complémentaires
            }
        }
        for (set<int>::iterator it=sonSet.begin(); it != sonSet.end(); ++it){
            sons.clear();
            aretes.clear();
            depthSons.clear();
            labelSons.clear();
            while (!aVoir.empty()) { // N'est jamais censé arriver !
                aVoir.pop();
            }
            while (!depth.empty()) { // N'est jamais censé arriver !
                depth.pop();
            }
            while (!labels.empty()) { // N'est jamais censé arriver !
                labels.pop();
            }
            vu2.clear();
            vu2.insert((*it)); //On voit bien le sommet duquel on part
            for (vector<Neighbor>::iterator voisin = G.Neighbors((*it))->begin();
                 voisin != G.Neighbors((*it))->end(); ++voisin) {
                aVoir.push(&(*voisin));
                depth.push(1);
                labels.push("");
            }
            G.BFS_comp_GraphDedupli(seen, vu2, aVoir, depth, labels, sons, depthSons, aretes, labelSons);
            for (int i = 0; i < sons.size(); i++){
                areteFF[make_pair((*it), sons[i])] = labelSons[i];
            }
        } //On termine de traiter tous les sommets de la composante
        cout << "BFS sur tous les sommets terminés" << endl;


        //On peut donc passer la construction du graphe. Commençons par les sommets en péri.
        for (set<int>::iterator setIt = sonSet.begin(); setIt != sonSet.end(); ++setIt) {
            V3.push_back(Node(index, G.Vertices[(*setIt)].weight, G.Vertices[(*setIt)].label));
            correspondingVertex[(*setIt)] = index;
            index++;
            seen[(*setIt)] = -compteurDeBoucle;
        }

        //Puis les arêtes au sein de la composante :
        for (dicChem::iterator itDic = areteFF.begin(); itDic != areteFF.end(); ++itDic) {
            if (itDic->second.size()==0){ //Ce cas là arrive lorsque les sommets i et j sont adjacents
                continue;
            }
            V3.push_back(Node(index, 0, itDic->second));
            E3.push_back(Edge(correspondingVertex[itDic->first.first], index, 0, aretFF));
            E3.push_back(Edge(index, correspondingVertex[itDic->first.second], 0, aretFF));
            index++;
        }
    } //On vient de terminer cette composante !

    //Maintenant que les composantes ont bien été ajouté, on s'occupe des sommets restants
    for (int i = 0; i < G.N; i++) {
        if (seen[i] == 0 and G.Vertices[i].label.size()>=G.kmer) {
            V3.push_back(Node(index, G.Vertices[i].weight, G.Vertices[i].label));
            correspondingVertex[i] = index;
            index++;
        } else if (seen[i] < 0) {
            continue; //Cas déjà traité précédemment
        } else {
            correspondingVertex[i] = -1; //Cas sommet de comp ou de label trop court
        }
    }
    cout << "Sommets restants ajoutés" << endl;
    limiteAretes = E3.size();
    //Maintenant on s'occupe des arêtes
    for (int i = 0; i<G.N ; i++) {
        if (seen[i] == 0 and correspondingVertex[i]>=0){//Cas sommet non vu et valide
            for(vector<Neighbor>::iterator node = G.Neighbors(i)->begin(); node != G.Neighbors(i)->end(); ++node){
                if (seen[node->val] == 0 and correspondingVertex[node->val]>=0){ //Cas voisin valide
                    E3.push_back(Edge(correspondingVertex[i],correspondingVertex[node->val],0,node->label));
                } else if (seen[node->val] < 0 and correspondingVertex[node->val]>=0
                and not foundedEdge(seen[i],seen[node->val],E3,limiteAretes)){ //Cas voisin valide
                        E3.push_back(Edge(correspondingVertex[i],correspondingVertex[node->val],0,node->label));
                    }
            }
        } else if (seen[i] < 0 ){ //Cas sommet de comp en péri
            for (vector<Neighbor>::iterator node = G.Neighbors(i)->begin(); node != G.Neighbors(i)->end(); ++node) {
                if (seen[node->val] == 0 and correspondingVertex[node->val]>=0) { //Voisin hors comp et valide
                    E3.push_back(Edge(correspondingVertex[i], correspondingVertex[node->val], 0, node->label));
                } else if (seen[node->val] < 0 and correspondingVertex[node->val]>=0 and seen[i] != seen[node->val] and
                    not foundedEdge(seen[i],seen[node->val],E3,limiteAretes)) { //Voisin hors comp et valide
                    E3.push_back(Edge(correspondingVertex[i], correspondingVertex[node->val], 0, node->label));
                }
            }
        }
    }
    cout << "Arêtes restantes ajoutées" << endl;

    //Enfin, on enregistre le graphe créé.
    ofstream outputNodes2;
    outputNodes2.open(outputPrefix + "/clean.nodes");
    printVerticesBcalm(V3, outputNodes2);
    outputNodes2.close();

    ofstream outputEdges2;
    outputEdges2.open(outputPrefix + "/clean.edges");
    printEdgesBcalm(E3, outputEdges2);
    outputEdges2.close();

    //On récupère aussi l'abondance qui est nécessaire pour kissplice
    ifstream ab(edgesPath.substr(0,edgesPath.size()-12)+".abundance", std::ios::binary);

    vector <double> A;
    read_abundance_file(ab,A);
    vector<double> A3(V3.size(),0.0);
    for (int i = 0; i<A.size();i++){
        if(correspondingVertex[i] >= 0){
            A3[correspondingVertex[i]]+= A[i];
        }
    }

    ofstream outputAb;
    outputAb.open(outputPrefix + "/clean.abundance");
    printAbundance(A3,outputAb);
    cout << "Graphe enregistré" << endl;

}