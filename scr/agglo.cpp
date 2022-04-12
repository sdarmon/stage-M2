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
#include "graph.h"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wstring-compare"
#define MAX 1024

// Permet de déterminer s'il existe un chemin de la composante i à la composante j, et s'il existe, renvoie la distance.

void chemin_local(int i, vector<int> &endings, Graph &G, vector<vector<int>> &components, int profondeurMax){
    vector<int> comp1;
    vector<Neighbor*> voisin1;
    vector<int> pronf;
    pronf.clear();
    voisin1.clear();
    comp1.clear();

    //Étape d'initialisation:
    for(int k= 0; k< components[i].size(); k++){
        comp1.push_back(components[i][k]); //On récupère les sommets de la composante que l'on souhaite étudier

        //On récupère les voisins de ces sommets
        for (vector<Neighbor>::iterator it = G.Neighbors(components[i][k])->begin(); it != G.Neighbors(components[i][k])->end(); ++it){
            if (find(comp1.begin(),comp1.end(), it->val) == comp1.end()){ 
                voisin1.push_back(&(*it));
                pronf.push_back(G.Vertices[it->val].label.size()-40);
            }
        }
    }
    //On regarde si la j-eme composante s'intersecte avec la comp1 (possible en fonction du critère de la composante
    //choisi
    for (int j = 0; j<components.size(); j++){
        for(int k= 0; k< components[j].size(); k++){
            if (find(comp1.begin(),comp1.end(),components[j][k]) != comp1.end()){ 
                endings[j] = 0;
                continue;
            }
        }
    }

    int step = 1; //On effectue le BFS par couches successives : on traite à chaque boucle les sommets à distance `step`
    //de la composante
    int modif = 1; //Booléen indiquant s'il y a eu une modification à la `step`-ème étape du BFS
    int size;
    Neighbor* node;
    int pro; //Variable indiquant à quelle profondeur (en termes de nucléotides) de la composante

    while (modif && step < profondeurMax) { //Tant qu'un sommet a été rajouté à la couche précédente, on regarde tous les sommets de cette couche-là dans le BFS.
        modif = 0;
        size = voisin1.size();
        for (int k = 0; k<size; ++k){ //On boucle sur les sommets de la couche du BFS uniquement
            node = voisin1.front();
            pro = pronf.front();
            voisin1.erase(voisin1.begin());
            pronf.erase(pronf.begin());
            if (find(comp1.begin(),comp1.end(), node->val) != comp1.end()){ // Sommet déjà présent
                continue;
            }
            // Sommet également présent dans la composante d'arrivée
            for (int j = 0; j<components.size(); j++){
                if (find(components[j].begin(),components[j].end(),node->val) != components[j].end()){ 
                    if (endings[j] == -1){ //On garde toujours le plus court chemin
                        endings[j] = pro;
                    } else {
                        endings[j] = min(endings[j],pro);
                    }
                    continue;
                }
            }
            //Sinon on rajoute ce sommet comme vu et on ajoute ses voisins à traiter dans le BFS pour la prochaine couche
            comp1.push_back(node->val);
            modif = 1;
            for (vector<Neighbor>::iterator it = G.Neighbors(node->val)->begin(); it != G.Neighbors(node->val)->end(); ++it){
                if (it->label[0] != node->label[1]){
                    continue;
                }
                voisin1.push_back(&(*it));
                pronf.push_back(pro+G.Vertices[it->val].label.size()-40);
            }
        }
        step ++; //On passe à la couche suivante
    }
    //Aucun nouveau sommet n'a été rajouté à traiter... Il n'y a plus de chemins à visiter...
    return;
}

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

int main(int argc, char** argv) {
    if (argc != 8 and argc != 9) {
        cout << "Expected use of this program: \n\n\t" << argv[0]
             << " file.nodes file.edges -c value -d dis outputPrefix [-clean]\n" << endl;
        return 0;
    }

    //On charge le graphe
    vector <Edge> E;
    vector <Node> V;

    char *nodesPath = argv[1];
    char *edgesPath = argv[2];
    int dis = atoi(argv[6]);

    ifstream edges(edgesPath, std::ios::binary);
    ifstream nodes(nodesPath, std::ios::binary);

    read_edge_file(edges, E);
    read_node_file_weighted(nodes, V);

    Graph G(V, E);

    cout << "Graphe chargé et construit" << endl;

    int index;
    //Construction du vector `vu_total`
    vector<int> vu_total(G.N, 0);

    vector <vector<int>> components;
    components.clear();

    //Ici le critère choisi pour définir une composante est juste d'être le sous-graphe connexe maximum du graphe G
    //auquel on a enlevé tous les sommets de poids inférieur au seuil `threshold`.
    int threshold;
    threshold = atoi(argv[4]);
    int m = 0;
    int weight;
    vector<int> vu(G.N, 0);
    vector < Neighbor * > aVoir;

    index = indexMax(G, vu_total);
    weight = G.Vertices[index].weight;
    vu_total[index] = 1;

    cout << "Début de la recherche des composantes" << endl;
    while (weight >= threshold and m < 100) //On cherche des composantes tant qu'il existe encore un sommet vérifiant
        //le critère et dans la limite des 100 composantes.
    {
        //On réalise alors un BFS depuis le sommet/graine `index`
        for (int i = 0; i < G.N; i++) {
            vu[i] = 0;
        }
        vector<int> compo;
        compo.clear();
        aVoir.clear();
        vu[index] = 1;
        for (vector<Neighbor>::iterator it = G.Neighbors(index)->begin(); it != G.Neighbors(index)->end(); ++it) {
            if (G.Vertices[it->val].weight >= threshold) {
                aVoir.push_back(&(*it));
            }
        }
        G.BFS_func(threshold, aVoir, vu);

        //Puis on ajoute les sommets trouvés à la composante
        for (int i = 0; i < G.N; i++) {
            if (vu[i]) {
                compo.push_back(i);
                vu_total[i] = 1;
            }
        }
        save_comp(G, compo, argv[7], m); //On enregistre la composante sur l'ordinateur
        components.push_back(compo); //Et on ajoute la composante au vecteur de composantes
        cout << "Composante trouvée de départ " << index << " et de poids " << G.Vertices[index].weight
             << " et de taille " << compo.size() << endl;

        //Finalement, on recommence la boucle while
        m++;
        index = indexMax(G, vu_total);
        weight = G.Vertices[index].weight;
        vu_total[index] = 1;

    }

    cout << "Fin de la recherche de composantes." << endl;
    if (argc < 9 ) {
    //Maintenant on construit nouveau graphe contracté
    vector <Edge> E2;
    vector <Node> V2;
    vector<int> endings;
    E2.clear();
    V2.clear();

    for (int i = 0; i < components.size(); i++) {//On boucle sur les nouveaux sommets i.e. les composantes
        cout << "Calcul des chemins partant de " << i << endl;
        initVec(endings, components.size());
        chemin_local(i, endings, G, components, dis); //On calcule s'il existe des chemins entre les composantes

        V2.push_back(Node(i, components[i].size(), ""));//On rajoute le sommet au graphe

        for (int j = 0; j < components.size(); j++) {//Puis on rajoute les arêtes au graphe
            if (endings[j] >= 0) {
                char arete[3] = {'N', 'N'};
                E2.push_back(Edge(i, j, endings[j], arete));
            }
        }
    }

    Graph H(V2, E2);

    cout << "Graphe aggloméré construit" << endl;

    //Enfin, on enregistre le graphe créé.
    ofstream outputNodes;
    outputNodes.open(string(argv[7]) + "/agglo.nodes");
    printGraphVertices(H, outputNodes);
    outputNodes.close();

    ofstream outputEdges;
    outputEdges.open(string(argv[7]) + "/agglo.edges");
    printEdges(E2, outputEdges);
    outputEdges.close();
}
    if (argc == 9) {
        cout << "Début construction graphe aggloméré" << endl;
        vector <Edge> E3;
        vector <Node> V3;
        E3.clear();
        V3.clear();
        vector<int> seen(G.N, 0);
        vector<int> correspondingVertex(G.N, 0);
        int compt = 0;
        //seen est le tableau de correspondance entre les anciens sommets et les nouveaux.
        //Pour les sommets des composantes on les flag avec des indices positifs
        for (vector < vector < int >> ::iterator comp = components.begin(); comp != components.end();
        ++comp){
            compt++;
            for (vector<int>::iterator it = comp->begin(); it != comp->end(); ++it) {
                seen[(*it)] = compt;
            }
        }

        //Ensuite, on regarde qui est sur le périmètre, on leur donne un indice négatif.
        int out;
        for (vector < vector < int >> ::iterator comp = components.begin(); comp != components.end();
        ++comp){
            for (vector<int>::iterator it = comp->begin(); it != comp->end(); ++it) {
                out = 0;
                for(vector<Neighbor>::iterator it2 = G.Adj[(*it)].adjVec.begin(); it2 != G.Adj[(*it)].adjVec.end(); ++it2){
                    if (seen[it2->val] == 0){
                        out = 1;
                        break;
                    }
                }
                if (out){
                    seen[(*it)] = min(-seen[(*it)],seen[(*it)]); //Pour être sûr de prendre une valeur négative;
                    // normalement le min n'est pas censé servir car composantes distantes d'au moins 2 deux à deux
                }
            }
        }
        cout << "Sommets en périmètre trouvés" << endl;
        //Après, on fait un BFS à partir de chaque sommet dans la composante afin d'obtenir les voisins du périmètre.
        set<int> vu;
        vector<vector<int>> neighborsPeri;
        vector<vector<Neighbor*>> neighborsAretes;
        vector<int> limiteAretes;
        vector<int> indexation;
        vector<int> fusion;
        set<int> inter;
        vector<set<int>> setVoisin;
        index = 0;
        int modif;
        int compteurDeBoucle = 0;
        //On boucle sur les composantes
        for (vector < vector < int >> ::iterator comp = components.begin(); comp != components.end();
        ++comp){
            neighborsPeri.clear(); //Vecteur qui contiendra l'ensemble des voisins de chaque sommet en périmètre.
            neighborsAretes.clear();
            indexation.clear();
            limiteAretes.clear();
            //On boucles sur les sommets de la composante
            cout << "Pré-calcul pour la composante "<<compteurDeBoucle << endl;
            compteurDeBoucle++;
            for (vector<int>::iterator it = comp->begin(); it != comp->end(); ++it) {
                if (seen[(*it)] < 0) { //Cas sommet en périphérie
                    vector<int> sons;
                    sons.clear();
                    vector<Neighbor*> aretes;
                    aretes.clear();
                    aVoir.clear();
                    vu.clear();
                    vu.insert((*it)); //On voit bien le sommet duquel on part

                    //On fait un BFS à partir de chaque sommet afin de savoir quels sommets du périmètre sont
                    //atteignables à partir de chaque sommet du périmètre
                    //Cas en partant du foward
                    for (vector<Neighbor>::iterator voisin = G.Neighbors((*it))->begin();
                         voisin != G.Neighbors((*it))->end(); ++voisin) {
                        if (voisin->label[0] == 'F') {
                            aVoir.push_back(&(*voisin));
                        }
                    }
                    G.BFS_comp(seen, vu, aVoir, sons, aretes);

                    limiteAretes.push_back(aretes.size()); //On fait ça pour garder en mémoire le fait que ce ne sont
                    //pas les mêmes arêtes

                    //Cas en partant du reverse
                    for (vector<Neighbor>::iterator voisin = G.Neighbors((*it))->begin();
                         voisin != G.Neighbors((*it))->end(); ++voisin) {
                        if (voisin->label[0] == 'R') {
                            aVoir.push_back(&(*voisin));
                        }
                    }
                    G.BFS_comp(seen, vu, aVoir, sons, aretes);

                    neighborsPeri.push_back(sons);
                    neighborsAretes.push_back(aretes);
                    indexation.push_back((*it)); //On garde en mémoire la correspondance entre l'index et la position
                    //du sommet
                }
            } //On termine de traiter tous les sommets de la composante
            cout << "BFS sur tous les sommets terminés" << endl;

            //On trie par cardinal décroissant le vecteur neighborsPeri, un tri bulle suffit car algo suivant en n²
            //aussi
            for (int i = indexation.size() - 1; i > 0; i--) {
                modif = 0;
                for (int j = 0; j < i; j++) {
                    if (neighborsPeri[j].size() < neighborsPeri[j+1].size()) {
                        swap(neighborsPeri[j + 1], neighborsPeri[j]);
                        swap(neighborsAretes[j + 1], neighborsAretes[j]);
                        swap(indexation[j + 1], indexation[j]);
                        modif = 1;
                    }
                }
                if (!modif) {
                    break;
                }
            }
            cout << "Ensembles des " << indexation.size() << " voisins triés" << endl;
            //Ici, on ne veut garder que la plus grande antichaine de sommets (en terme d'inclusion des voisinages)
            //Et alors, on veut fusionner les sommets comparables.
            fusion.clear(); //Ce vecteur contiendra -1 si le sommet ne doit pas fusionner dans un autre sommet, et
            //i s'il correspond à un sous-ensemble du i-ème sommet

            //On récupère les voisins hors composante
            setVoisin.clear();
            for (int i = 0; i < indexation.size(); i++) {
                set<int> setA;
                setA.clear();
                for (vector<int>::iterator el =neighborsPeri[i].begin(); el!= neighborsPeri[i].end(); ++el ){
                    if (seen[(*el)] == 0){
                        setA.insert((*el));
                    }
                }
                setVoisin.push_back(setA);
            }

            //Ainsi on peut parcourir les sommets pour former l'antichaine
            for (int i = 0; i < indexation.size(); i++) {
                fusion.push_back(-1);
                inter.clear();
                //Maintenant, on regarde si un précédent cas correspond à un sur-ensemble de notre i-ème cas
                for (int j = 0; j < i; j++) {
                    if (fusion[j] < 0 and setVoisin[i].size() <=
                                                  setVoisin[j].size()) { //Cas où le sommet n'est pas déjà marqué comme à fusionner avec un autre sommet
                        //Pour savoir si l'un est inclus dans l'autre, on fait l'intersection des deux puis on
                        //vérifie si le cardinal de l'intersection correspond à celui de l'un des deux sommets
                        set_intersection(setVoisin[i].begin(), setVoisin[i].end(),
                                         setVoisin[j].begin(), setVoisin[j].end(), inserter(inter, inter.begin()));
                        if (inter.size() == setVoisin[i].size()) { //Cas i inclus dans j; on rappelle que le
                            //vecteur est trié par cardinal décroissant.
                            fusion[i] = j;
                            seen[indexation[i]] = max(-seen[indexation[i]],seen[indexation[i]]); //On retire le sommet comme étant en péri
                            break;
                        }
                    }
                }
            }
            //Ici on est sûr que fusion indique bien les sommets étant dans sous-cas des autres.
            cout << "Antichaine maximal de l'ensemble de voisins faite" << endl;
            //On peut donc passer la construction du graphe. Commençons par les sommets.
            for (int i = 0; i < indexation.size(); i++) {
                if (fusion[i] < 0) {
                    V3.push_back(Node(index,G.Vertices[indexation[i]].weight,G.Vertices[indexation[i]].label));
                    correspondingVertex[indexation[i]] = index;
                    index++;
                } else {
                    correspondingVertex[indexation[i]] = correspondingVertex[indexation[fusion[i]]];
                }
            }

            //Puis les arêtes au sein de la composante :
            for (int i = 0; i < indexation.size(); i++) {
                if (seen[indexation[i]] < 0) {
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
            }
        } //On vient de terminer cette composante !

        //Maintenant que les composantes ont bien été ajouté, on s'occupe des sommets restants
        for (int i = 0; i < G.N; i++) {
            if (seen[i] == 0) {
                V3.push_back(Node(index, G.Vertices[i].weight, G.Vertices[i].label));
                correspondingVertex[i] = index;
                index++;
            }
        }
        cout << "Sommets restants ajoutés" << endl;

        //Maintenant on s'occupe des arêtes
        for (int i = 0; i<G.N ; i++) {
            if (seen[i] == 0){
                for(vector<Neighbor>::iterator node = G.Neighbors(i)->begin(); node != G.Neighbors(i)->end(); ++node){
                    E3.push_back(Edge(correspondingVertex[i],correspondingVertex[node->val],0,node->label));
                }
            } else {
                for (vector<Neighbor>::iterator node = G.Neighbors(i)->begin(); node != G.Neighbors(i)->end(); ++node) {
                    if (seen[node->val] < 0) { //Ici on regarde que les négatifs car tous les positifs ont été fusionné
                        //avec des sur-ensembles de leurs voisins, donc normalement c'est bon, on n'a pas d'arêtes en double!
                        E3.push_back(Edge(correspondingVertex[i], correspondingVertex[node->val], 0, node->label));
                    }
                }
            }
        }
        cout << "Arêtes restantes ajoutées" << endl;

        //Enfin, on enregistre le graphe créé.
        ofstream outputNodes2;
        outputNodes2.open(string(argv[7]) + "/clean.nodes");
        printVerticesBcalm(V3, outputNodes2);
        outputNodes2.close();

        ofstream outputEdges2;
        outputEdges2.open(string(argv[7]) + "/clean.edges");
        printEdgesBcalm(E3, outputEdges2);
        outputEdges2.close();
        cout << "Graphe enregistré" << endl;
    }
}