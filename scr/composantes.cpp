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
#include <map>
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
    if (argc != 9 and argc != 7) {
        cout << "Expected use of this program: \n\n\t" << argv[0]
             << " file.nodes file.edges -c value [-k kmer] outputPrefix [-with/without]\n" << endl;
        return 0;
    }

    //On charge le graphe
    vector <Edge> E;
    vector <Node> V;

    string nodesPath = argv[1];
    string edgesPath = argv[2];
    string outputPrefix;
    int with;
    ifstream edges(edgesPath, std::ios::binary);
    ifstream nodes(nodesPath, std::ios::binary);

    read_edge_file(edges, E);
    read_node_file_weighted(nodes, V);

    Graph G(V, E);
    if(argc >= 8 and argv[5][1]=='k' ){
        G.kmer = stoi(argv[6]);
        outputPrefix = argv[7];
        if (((string)argv[8]).size() == 5){
            with = 1;
        } else {
            with = 0;
        }
    } else {
        outputPrefix = argv[5];
        if (((string)argv[6]).size() == 5){
            with = 1;
        } else {
            with = 0;
        }
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
        if (G.Vertices[sommet].weight >= threshold){
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
    set<int> sonSet;
    typedef map<pair<int,int>,pair<int,int>> dic;
    typedef map<pair<int,int>,string> dicChem;
    dic BulF_FF; //la clé est le couple (i,j) trouvé et la valeur associé est le couple (s,taille)
    dic BulF_FR;
    dic BulF_RF;
    dic BulF_RR;
    dic BulR_FF; //la clé est le couple (i,j) trouvé et la valeur associé est le couple (s,taille)
    dic BulR_FR;
    dic BulR_RF;
    dic BulR_RR;
    dicChem areteFF;
    dicChem areteFR;
    dicChem areteRF;
    dicChem areteRR;
    int limiteAretes;
    vector<int> sons;
    vector<int> depthSons;
    queue<int> depth;
    vector<string> labelSons;
    queue<string> labels;
    vector<Neighbor*> aretes;
    index = 0;
    int pos;
    int modif;
    int sizeLabel;
    int I;
    int J;
    int compteurDeBoucle = 0;
    char aretFF[3] = {'F', 'F'};
    char aretFR[3] = {'F', 'R'};
    char aretRF[3] = {'R', 'F'};
    char aretRR[3] = {'R', 'R'};
    //On boucle sur les composantes
    for (vector < vector < int >> ::iterator comp = components.begin(); comp != components.end();
    ++comp){
        BulF_FF.clear();
        BulF_RF.clear();
        BulF_FR.clear();
        BulF_RR.clear();
        BulR_FF.clear();
        BulR_RF.clear();
        BulR_FR.clear();
        BulR_RR.clear();
        areteFF.clear();
        areteFR.clear();
        areteRR.clear();
        areteRF.clear();
        sonSet.clear();
        //On boucle sur les sommets de la composante
        if (comp->size() > 5000) {
        cout << "Composante " << compteurDeBoucle << " passée car trop grosses (de taille " << comp->size() << ")" << endl; }
        compteurDeBoucle++;
        for (vector<int>::iterator it = comp->begin(); it != comp->end(); ++it){
            seen[(*it)] = 0;
        }
        else {
            cout << "Pré-calcul pour la composante " << compteurDeBoucle << " de taille " << comp->size() << endl;
            compteurDeBoucle++;
            for (vector<int>::iterator it = comp->begin(); it != comp->end(); ++it) {
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
                vu2.insert((*it)); //On voit bien le sommet duquel on part

                //On fait un BFS à partir de chaque sommet afin de savoir quels sommets du périmètre sont
                //atteignables à partir de chaque sommet de la composante
                //Cas en partant du forward
                for (vector<Neighbor>::iterator voisin = G.Neighbors((*it))->begin();
                     voisin != G.Neighbors((*it))->end(); ++voisin) {
                    if (voisin->label[0] == 'F') {
                        aVoir.push(&(*voisin));
                        depth.push(1);
                        string aux = G.Vertices[(*it)].label;
                        labels.push(aux);
                    }
                }
                G.BFS_comp(seen, vu2, aVoir, depth, labels, sons, depthSons, aretes, labelSons);
                for (int i = 0; i < sons.size(); i++) {
                    sonSet.insert(sons[i]);
                }
                if (!with) {
                    for (int i = 0; i < sons.size(); i++) {
                        for (int j = i + 1; j < sons.size(); j++) {
                            if (sons[i] < sons[j]) {
                                I = i;
                                J = j;
                            } else {
                                I = j;
                                J = i;
                            }
                            if (aretes[I]->label[1] == 'F' and aretes[J]->label[1] == 'F') {
                                if (BulF_FF.find(make_pair(sons[I], sons[J])) == BulF_FF.end() or
                                    BulF_FF[make_pair(sons[I], sons[J])].second > depthSons[I] + depthSons[J]) {
                                    BulF_FF[make_pair(sons[I], sons[J])] = make_pair((*it),
                                                                                     depthSons[I] + depthSons[J]);
                                }
                            } else if (aretes[I]->label[1] == 'F' and aretes[J]->label[1] == 'R') {
                                if (BulF_FR.find(make_pair(sons[I], sons[J])) == BulF_FR.end() or
                                    BulF_FR[make_pair(sons[I], sons[J])].second > depthSons[I] + depthSons[J]) {
                                    BulF_FR[make_pair(sons[I], sons[J])] = make_pair((*it),
                                                                                     depthSons[I] + depthSons[J]);
                                }

                            } else if (aretes[I]->label[1] == 'R' and aretes[J]->label[1] == 'F') {
                                if (BulF_RF.find(make_pair(sons[I], sons[J])) == BulF_RF.end() or
                                    BulF_RF[make_pair(sons[I], sons[J])].second > depthSons[I] + depthSons[J]) {
                                    BulF_RF[make_pair(sons[I], sons[J])] = make_pair((*it),
                                                                                     depthSons[I] + depthSons[J]);
                                }

                            } else {
                                if (BulF_RR.find(make_pair(sons[I], sons[J])) == BulF_RR.end() or
                                    BulF_RR[make_pair(sons[I], sons[J])].second > depthSons[I] + depthSons[J]) {
                                    BulF_RR[make_pair(sons[I], sons[J])] = make_pair((*it),
                                                                                     depthSons[I] + depthSons[J]);
                                }
                            }
                        }
                    }
                }

                limiteAretes = aretes.size(); //On fait ça pour garder en mémoire le fait que ce ne sont
                //pas les mêmes arêtes

                //Cas en partant du reverse
                for (vector<Neighbor>::iterator voisin = G.Neighbors((*it))->begin();
                     voisin != G.Neighbors((*it))->end(); ++voisin) {
                    if (voisin->label[0] == 'R') {
                        aVoir.push(&(*voisin));
                        depth.push(1);
                        labels.push("");
                    }
                }
                G.BFS_comp(seen, vu2, aVoir, depth, labels, sons, depthSons, aretes, labelSons);

                for (int i = limiteAretes; i < sons.size(); i++) {
                    sonSet.insert(sons[i]);
                }
                if (!with) {
                    for (int i = limiteAretes; i < sons.size(); i++) {
                        for (int j = i + 1; j < sons.size(); j++) {
                            if (sons[i] < sons[j]) {
                                I = i;
                                J = j;
                            } else {
                                I = j;
                                J = i;
                            }
                            if (aretes[I]->label[1] == 'F' and aretes[J]->label[1] == 'F') {
                                if (BulR_FF.find(make_pair(sons[I], sons[J])) == BulR_FF.end() or
                                    BulR_FF[make_pair(sons[I], sons[J])].second > depthSons[I] + depthSons[J]) {
                                    BulR_FF[make_pair(sons[I], sons[J])] = make_pair((*it),
                                                                                     depthSons[I] + depthSons[J]);
                                }
                            } else if (aretes[I]->label[1] == 'F' and aretes[J]->label[1] == 'R') {
                                if (BulR_FR.find(make_pair(sons[I], sons[J])) == BulR_FR.end() or
                                    BulR_FR[make_pair(sons[I], sons[J])].second > depthSons[I] + depthSons[J]) {
                                    BulR_FR[make_pair(sons[I], sons[J])] = make_pair((*it),
                                                                                     depthSons[I] + depthSons[J]);
                                }

                            } else if (aretes[I]->label[1] == 'R' and aretes[J]->label[1] == 'F') {
                                if (BulR_RF.find(make_pair(sons[I], sons[J])) == BulR_RF.end() or
                                    BulR_RF[make_pair(sons[I], sons[J])].second > depthSons[I] + depthSons[J]) {
                                    BulR_RF[make_pair(sons[I], sons[J])] = make_pair((*it),
                                                                                     depthSons[I] + depthSons[J]);
                                }

                            } else {
                                if (BulR_RR.find(make_pair(sons[I], sons[J])) == BulR_RR.end() or
                                    BulR_RR[make_pair(sons[I], sons[J])].second > depthSons[I] + depthSons[J]) {
                                    BulR_RR[make_pair(sons[I], sons[J])] = make_pair((*it),
                                                                                     depthSons[I] + depthSons[J]);
                                }
                            }
                        }
                    }
                }
                if (with) {
                    for (int i = 0; i < limiteAretes; i++) {
                        for (int j = limiteAretes; j < sons.size(); j++) {
                            string aux = reverse_complement(labelSons[j]) + labelSons[i];
                            if (aretes[i]->label[1] == 'F' and
                                aretes[j]->label[1] == 'F') { //attention, ici la seconde lettre doit être comprise
                                //comme le complément ! (Voir cahier, ce n'est pas forcément évident)
                                if (areteRF.find(make_pair(sons[j], sons[i])) == areteRF.end() or
                                    areteRF[make_pair(sons[j], sons[i])].size() > aux.size()) {
                                    areteRF[make_pair(sons[j], sons[i])] = aux;
                                }
                            } else if (aretes[i]->label[1] == 'F' and aretes[j]->label[1] == 'R') {
                                if (areteFF.find(make_pair(sons[j], sons[i])) == areteFF.end() or
                                    areteFF[make_pair(sons[j], sons[i])].size() > aux.size()) {
                                    areteFF[make_pair(sons[j], sons[i])] = aux;
                                }
                            } else if (aretes[i]->label[1] == 'R' and aretes[j]->label[1] == 'R') {
                                if (areteFR.find(make_pair(sons[j], sons[i])) == areteFR.end() or
                                    areteFR[make_pair(sons[j], sons[i])].size() > aux.size()) {
                                    areteFR[make_pair(sons[j], sons[i])] = aux;
                                }
                            } else {
                                if (areteRR.find(make_pair(sons[j], sons[i])) == areteRR.end() or
                                    areteRR[make_pair(sons[j], sons[i])].size() > aux.size()) {
                                    areteRR[make_pair(sons[j], sons[i])] = aux;
                                }
                            }
                        }
                    }
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

            //On continue par les sommets à dédoubler et les arêtes associées :


            if (!with) {
                for (dic::iterator itDic = BulF_FF.begin(); itDic != BulF_FF.end(); ++itDic) {
                    V3.push_back(Node(index, 0, G.Vertices[itDic->second.first].label));
                    E3.push_back(Edge(index, correspondingVertex[itDic->first.first], 0, aretFF));
                    E3.push_back(Edge(index, correspondingVertex[itDic->first.second], 0, aretFF));
                    E3.push_back(Edge(correspondingVertex[itDic->first.first], index, 0, aretRR));
                    E3.push_back(Edge(correspondingVertex[itDic->first.second], index, 0, aretRR));
                    index++;
                }
                for (dic::iterator itDic = BulF_FR.begin(); itDic != BulF_FR.end(); ++itDic) {
                    V3.push_back(Node(index, 0, G.Vertices[itDic->second.first].label));
                    E3.push_back(Edge(index, correspondingVertex[itDic->first.first], 0, aretFF));
                    E3.push_back(Edge(index, correspondingVertex[itDic->first.second], 0, aretFR));
                    E3.push_back(Edge(correspondingVertex[itDic->first.first], index, 0, aretRR));
                    E3.push_back(Edge(correspondingVertex[itDic->first.second], index, 0, aretFR));
                    index++;
                }
                for (dic::iterator itDic = BulF_RF.begin(); itDic != BulF_RF.end(); ++itDic) {
                    V3.push_back(Node(index, 0, G.Vertices[itDic->second.first].label));
                    E3.push_back(Edge(index, correspondingVertex[itDic->first.first], 0, aretFR));
                    E3.push_back(Edge(index, correspondingVertex[itDic->first.second], 0, aretFF));
                    E3.push_back(Edge(correspondingVertex[itDic->first.first], index, 0, aretFR));
                    E3.push_back(Edge(correspondingVertex[itDic->first.second], index, 0, aretRR));
                    index++;
                }
                for (dic::iterator itDic = BulF_RR.begin(); itDic != BulF_RR.end(); ++itDic) {
                    V3.push_back(Node(index, 0, G.Vertices[itDic->second.first].label));
                    E3.push_back(Edge(index, correspondingVertex[itDic->first.first], 0, aretFR));
                    E3.push_back(Edge(index, correspondingVertex[itDic->first.second], 0, aretFR));
                    E3.push_back(Edge(correspondingVertex[itDic->first.first], index, 0, aretFR));
                    E3.push_back(Edge(correspondingVertex[itDic->first.second], index, 0, aretFR));
                    index++;
                }
                for (dic::iterator itDic = BulR_FF.begin(); itDic != BulR_FF.end(); ++itDic) {
                    V3.push_back(Node(index, 0, G.Vertices[itDic->second.first].label));
                    E3.push_back(Edge(index, correspondingVertex[itDic->first.first], 0, aretRF));
                    E3.push_back(Edge(index, correspondingVertex[itDic->first.second], 0, aretRF));
                    E3.push_back(Edge(correspondingVertex[itDic->first.first], index, 0, aretRF));
                    E3.push_back(Edge(correspondingVertex[itDic->first.second], index, 0, aretRF));
                    index++;
                }
                for (dic::iterator itDic = BulR_FR.begin(); itDic != BulR_FR.end(); ++itDic) {
                    V3.push_back(Node(index, 0, G.Vertices[itDic->second.first].label));
                    E3.push_back(Edge(index, correspondingVertex[itDic->first.first], 0, aretRF));
                    E3.push_back(Edge(index, correspondingVertex[itDic->first.second], 0, aretRR));
                    E3.push_back(Edge(correspondingVertex[itDic->first.first], index, 0, aretRF));
                    E3.push_back(Edge(correspondingVertex[itDic->first.second], index, 0, aretFF));
                    index++;
                }
                for (dic::iterator itDic = BulR_RF.begin(); itDic != BulR_RF.end(); ++itDic) {
                    V3.push_back(Node(index, 0, G.Vertices[itDic->second.first].label));
                    E3.push_back(Edge(index, correspondingVertex[itDic->first.first], 0, aretRR));
                    E3.push_back(Edge(index, correspondingVertex[itDic->first.second], 0, aretRF));
                    E3.push_back(Edge(correspondingVertex[itDic->first.first], index, 0, aretFF));
                    E3.push_back(Edge(correspondingVertex[itDic->first.second], index, 0, aretRF));
                    index++;
                }
                for (dic::iterator itDic = BulR_RR.begin(); itDic != BulR_RR.end(); ++itDic) {
                    V3.push_back(Node(index, 0, G.Vertices[itDic->second.first].label));
                    E3.push_back(Edge(index, correspondingVertex[itDic->first.first], 0, aretRR));
                    E3.push_back(Edge(index, correspondingVertex[itDic->first.second], 0, aretRR));
                    E3.push_back(Edge(correspondingVertex[itDic->first.first], index, 0, aretFF));
                    E3.push_back(Edge(correspondingVertex[itDic->first.second], index, 0, aretFF));
                    index++;
                }
            }

            //Puis les arêtes au sein de la composante :
            if (with) {
                for (dicChem::iterator itDic = areteFF.begin(); itDic != areteFF.end(); ++itDic) {
                    V3.push_back(Node(index, 0, itDic->second));
                    E3.push_back(Edge(correspondingVertex[itDic->first.first], index, 0, aretFF));
                    E3.push_back(Edge(index, correspondingVertex[itDic->first.second], 0, aretFF));
                    E3.push_back(Edge(correspondingVertex[itDic->first.second], index, 0, aretRR));
                    E3.push_back(Edge(index, correspondingVertex[itDic->first.first], 0, aretRR));
                    index++;
                }
                for (dicChem::iterator itDic = areteFR.begin(); itDic != areteFR.end(); ++itDic) {
                    V3.push_back(Node(index, 0, itDic->second));
                    E3.push_back(Edge(correspondingVertex[itDic->first.first], index, 0, aretFF));
                    E3.push_back(Edge(index, correspondingVertex[itDic->first.second], 0, aretFR));
                    E3.push_back(Edge(correspondingVertex[itDic->first.second], index, 0, aretFR));
                    E3.push_back(Edge(index, correspondingVertex[itDic->first.first], 0, aretRR));
                    index++;
                }
                for (dicChem::iterator itDic = areteRF.begin(); itDic != areteRF.end(); ++itDic) {
                    V3.push_back(Node(index, 0, itDic->second));
                    E3.push_back(Edge(correspondingVertex[itDic->first.first], index, 0, aretRF));
                    E3.push_back(Edge(index, correspondingVertex[itDic->first.second], 0, aretFF));
                    E3.push_back(Edge(correspondingVertex[itDic->first.second], index, 0, aretRR));
                    E3.push_back(Edge(index, correspondingVertex[itDic->first.first], 0, aretRF));
                    index++;
                }
                for (dicChem::iterator itDic = areteRR.begin(); itDic != areteRR.end(); ++itDic) {
                    V3.push_back(Node(index, 0, itDic->second));
                    E3.push_back(Edge(correspondingVertex[itDic->first.first], index, 0, aretRF));
                    E3.push_back(Edge(index, correspondingVertex[itDic->first.second], 0, aretFR));
                    E3.push_back(Edge(correspondingVertex[itDic->first.second], index, 0, aretFR));
                    E3.push_back(Edge(index, correspondingVertex[itDic->first.first], 0, aretRF));
                    index++;
                }
            }
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

    //Maintenant on s'occupe des arêtes
    for (int i = 0; i<G.N ; i++) {
        if (seen[i] == 0 and correspondingVertex[i]>=0){//Cas sommet non vu et valide
            for(vector<Neighbor>::iterator node = G.Neighbors(i)->begin(); node != G.Neighbors(i)->end(); ++node){
                if (correspondingVertex[node->val]>=0){ //Cas voisin valide
                    E3.push_back(Edge(correspondingVertex[i],correspondingVertex[node->val],0,node->label));
                }
            }
        } else if (seen[i] < 0 ){ //Cas sommet de comp en péri
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