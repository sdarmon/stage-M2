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
#include <queue>
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

struct indexDic {
    int key;
    int weight;
    friend bool operator< (indexDic const& lhs, indexDic const& rhs) {
        return (lhs.weight < rhs.weight);
    }
};

int main(int argc, char** argv) {
    if (argc != 8 and argc != 6) {
        cout << "Expected use of this program: \n\n\t" << argv[0]
             << " file.nodes file.edges -c value [-k kmer] outputPrefix \n" << endl;
        return 0;
    }

    //On charge le graphe
    vector <Edge> E;
    vector <Node> V;

    string nodesPath = argv[1];
    string edgesPath = argv[2];
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
    set<int> setVu;
    queue < Neighbor * > aVoir;

    cout << "Début du tas min" << endl;
    priority_queue<indexDic> indexTrie;
    for (int sommet = 0; sommet < G.N ; sommet++){
        if (G.Vertices[sommet].weight >= max(threshold,10)){
            indexTrie.push(indexDic{sommet,G.Vertices[sommet].weight});
        }
    } //Attention ici il y a une opti possible : on est en N log N mais c'est clairement possible de faire en N
    index = indexTrie.top().key;
    indexTrie.pop();
    vu_total[index] = 1;

    cout << "Début de la recherche des composantes" << endl;
    while (!indexTrie.empty()) //On cherche des composantes tant qu'il existe encore un sommet vérifiant
        //le critère et dans la limite des 100 composantes.
        //Condition m < 100 enlevée
    {
        //On réalise alors un BFS depuis le sommet `index`
        for (set<int>::iterator it = setVu.begin(); it != setVu.end() ; it++) {
            vu[(*it)] = 0;
        }
        setVu.clear();
        vector<int> compo;
        compo.clear();
        //aVoir étant une queue de BFS; est toujours censé être vide ici. Pour la sanité des algos qui suivent on la
        //vide au cas où...
        while (!aVoir.empty()){
            aVoir.pop();
        }
        vu[index] = 1;
        setVu.insert(index);
        for (vector<Neighbor>::iterator it = G.Neighbors(index)->begin(); it != G.Neighbors(index)->end(); ++it) {
            if (G.Vertices[it->val].weight >= threshold) {
                aVoir.push(&(*it));
            }
        }
        G.BFS_func(threshold, aVoir, vu, setVu);

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
            m++;
            components.push_back(compo); //Et on ajoute la composante au vecteur de composantes
            cout << "Composante trouvée de départ " << index << " et de poids " << G.Vertices[index].weight
                 << " et de taille " << compo.size() << endl;
        }

        //Finalement, on recommence la boucle while
        while (!indexTrie.empty() and vu_total[indexTrie.top().key]){
            indexTrie.pop();
        }
        if (!indexTrie.empty()){
            index = indexTrie.top().key;
            indexTrie.pop();
            vu_total[index] = 1;
        } else {
            break;
        }

    }

    cout << "Fin de la recherche de composantes." << endl;


    cout << "Début construction graphe aggloméré" << endl;
    vector <Edge> E3;
    vector <Node> V3;
    E3.clear();
    V3.clear();
    vector<int> seen(G.N, 0);
    vector<int> correspondingVertex(G.N, 0);
    compt = 0;
    //seen est le tableau de correspondance entre les anciens sommets et les nouveaux.
    //Pour les sommets des composantes on les flag avec des indices positifs
    for (vector < vector < int >> ::iterator comp = components.begin(); comp != components.end();
    ++comp){
        compt++;
        for (vector<int>::iterator it = comp->begin(); it != comp->end(); ++it) {
            seen[(*it)] = compt;
        }
    }

    //Après, on fait un BFS à partir de chaque sommet dans la composante afin d'obtenir les voisins du périmètre.
    set<int> vu2;
    vector<vector<int>> neighborsPeri;
    vector<vector<Neighbor*>> neighborsInF;
    vector<vector<Neighbor*>> neighborsInR;
    vector<vector<Neighbor*>> neighborsAretes;
    vector<int> limiteAretes;
    vector<int> indexation;
    vector<vector<int>> setVoisinInF;
    vector<vector<int>> setVoisinInR;
    index = 0;
    int pos;
    int modif;
    int compteurDeBoucle = 0;
    //On boucle sur les composantes
    for (vector < vector < int >> ::iterator comp = components.begin(); comp != components.end();
    ++comp){
        setVoisinInF.clear();
        setVoisinInR.clear();
        neighborsPeri.clear(); //Vecteur qui contiendra l'ensemble des voisins de chaque sommet en périmètre.
        neighborsAretes.clear();
        neighborsInF.clear();
        neighborsInR.clear();
        indexation.clear();
        limiteAretes.clear();
        //On boucle sur les sommets de la composante
        cout << "Pré-calcul pour la composante "<<compteurDeBoucle << endl;
        compteurDeBoucle++;
        for (vector<int>::iterator it = comp->begin(); it != comp->end(); ++it) {
            if (seen[(*it)] == 0) { //Cas sommet en périphérie
                vector<int> sons;
                sons.clear();
                vector<Neighbor*> aretes;
                aretes.clear();
                vector<int> inF;
                vector<int> inR;
                inF.clear();
                inR.clear();
                vector<Neighbor*> nodeInF;
                vector<Neighbor*> nodeInR;
                nodeInF.clear();
                nodeInR.clear();
                while (!aVoir.empty()){
                    aVoir.pop();
                }
                vu2.clear();
                vu2.insert((*it)); //On voit bien le sommet duquel on part

                //On fait un BFS à partir de chaque sommet afin de savoir quels sommets du périmètre sont
                //atteignables à partir de chaque sommet du périmètre
                //Cas en partant du foward
                for (vector<Neighbor>::iterator voisin = G.Neighbors((*it))->begin();
                     voisin != G.Neighbors((*it))->end(); ++voisin) {
                    if (voisin->label[0] == 'F') {
                        aVoir.push(&(*voisin));
                    } else {
                        pos= distance(inF.begin(), upper_bound(inF.begin(),inF.end(),voisin->val));
                        inF.insert(inF.begin()+pos,voisin->val);
                        nodeInF.insert(nodeInF.begin()+ pos,&(*voisin));
                    }
                }
                G.BFS_comp(seen, vu2, aVoir, sons, aretes);

                limiteAretes.push_back(aretes.size()); //On fait ça pour garder en mémoire le fait que ce ne sont
                //pas les mêmes arêtes

                //Cas en partant du reverse
                for (vector<Neighbor>::iterator voisin = G.Neighbors((*it))->begin();
                     voisin != G.Neighbors((*it))->end(); ++voisin) {
                    if (voisin->label[0] == 'R') {
                        aVoir.push(&(*voisin));
                    } else {
                        pos= distance(inR.begin(), upper_bound(inR.begin(),inR.end(),voisin->val));
                        inR.insert(inR.begin()+pos,voisin->val);
                        nodeInR.insert(nodeInR.begin()+ pos,&(*voisin));
                    }
                }
                G.BFS_comp(seen, vu2, aVoir, sons, aretes);

                neighborsPeri.push_back(sons);
                neighborsAretes.push_back(aretes);
                neighborsInF.push_back(nodeInF);
                neighborsInR.push_back(nodeInR);
                setVoisinInF.push_back(inF);
                setVoisinInR.push_back(inR);
                indexation.push_back((*it)); //On garde en mémoire la correspondance entre l'index et la position
                //du sommet
            }
        } //On termine de traiter tous les sommets de la composante
        cout << "BFS sur tous les sommets terminés" << endl;


        //On peut donc passer la construction du graphe. Commençons par les sommets.
        for (int i = 0; i < indexation.size(); i++) {
            V3.push_back(Node(index,G.Vertices[indexation[i]].weight,G.Vertices[indexation[i]].label));
            correspondingVertex[indexation[i]] = index;
            index++;
            seen[indexation[i]] = -1;
        }

        //Puis les arêtes au sein de la composante :
        for (int i = 0; i < indexation.size(); i++) {
            for (int j = 0; j < neighborsPeri[i].size(); j++) {
                if (seen[neighborsPeri[i][j]] < 0) {
                    if (j<limiteAretes[i]){
                        char aret[3] = {'F', neighborsAretes[i][j]->label[1]};
                        E3.push_back(Edge(correspondingVertex[indexation[i]],
                                          correspondingVertex[neighborsPeri[i][j]], 0, aret
                        ));
                    } else {
                        char aret[3] = {'R', neighborsAretes[i][j]->label[1]};
                        E3.push_back(Edge(correspondingVertex[indexation[i]],
                                          correspondingVertex[neighborsPeri[i][j]], 0, aret
                        ));
                    }
                }
            }
        }

    } //On vient de terminer cette composante !

    //Maintenant que les composantes ont bien été ajouté, on s'occupe des sommets restants
    for (int i = 0; i < G.N; i++) {
        if (seen[i] == 0 and G.Vertices[i].label.size()>G.kmer-1) {
            V3.push_back(Node(index, G.Vertices[i].weight, G.Vertices[i].label));
            correspondingVertex[i] = index;
            index++;
        } else if (seen[i] < 0) {
            continue; //Cas déjà traité précédemment
        } else {
            correspondingVertex[i] = -1; //Cas sommet fusionné ou de label trop court
        }
    }
    cout << "Sommets restants ajoutés" << endl;

    //Maintenant on s'occupe des arêtes
    for (int i = 0; i<G.N ; i++) {
        if (seen[i] == 0 and correspondingVertex[i]>=0){//Cas sommet hors comp et valide
            for(vector<Neighbor>::iterator node = G.Neighbors(i)->begin(); node != G.Neighbors(i)->end(); ++node){
                if (correspondingVertex[node->val]>=0){ //Cas voisin valide
                    E3.push_back(Edge(correspondingVertex[i],correspondingVertex[node->val],0,node->label));
                }
            }
        } else if (seen[i] < 0 ){
            for (vector<Neighbor>::iterator node = G.Neighbors(i)->begin(); node != G.Neighbors(i)->end(); ++node) {
                if (seen[node->val] == 0 and correspondingVertex[node->val]>=0) { //Voisin hors comp et valide
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