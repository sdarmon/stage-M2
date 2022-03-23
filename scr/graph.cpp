// ===========================================================================
// ===========================================================================
//                               Classe de graphes
// ===========================================================================
// ===========================================================================


/* 
 * Ce fichier permet de manipuler les graphes pondérés avec des classes adaptées 
 * Pour voir les definitions des classes, allez dans le fichier graph.h.
 * Ce fichier contient seulement les fonctions et méthodes liées aux classes
 * définies dans le fichier graph.h.
 */

/* 
 * WARNING
 * All the data is coded with int data type. For graph with more than 2147483648
 * (> 2x10^9) nodes or edges, you will get trouble!
 * Instead of int, use unsigned int (for the same weight, extend this upper bound
 * to 4x10^9) or long long int (for twice the weight, extension to 2^63 ~ 10^19).
 *
 */

// ===========================================================================
//                               Include Libraries
// ===========================================================================
#include <list>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <cstring>
#include <algorithm>
#include <sstream>

// ===========================================================================
//                             Include Project Files
// ===========================================================================
#include "graph.h"


// ===========================================================================
//                             Declare Used Namespaces
// ===========================================================================
using namespace std;
#define MAX 1024


// ===========================================================================
//                                  Constructors
// ===========================================================================


Node::Node(int v, int w, string l){
        val = v;
        weight = w;
        label=(string)l;
    }

Neighbor::Neighbor(int v, int w, char* l){
        val = v;
        weight = w;
        strcpy(label,l);
    }

LstNode::LstNode() { // Je crois que ce constructeur ne fonctionne pas, à éviter
        adjVec.clear();
    }

LstNode::LstNode(vector<Neighbor>& A) {
        adjVec = A;
    }


Edge::Edge(int s, int e, int w, char* l) {
        start = s;
        end = e;
        weight = w;
        strcpy(label,l);
    }

/* Un graphe est défini par ses sommets (vecteur de Nodes) et
 * ses arêtes (des listes d'adjacence i.e. un vecteur de LstNode).
 * On tient également en mémoire le nombre de sommets (N) et le 
 * nombre d'arête M.
 */
Graph::Graph(vector<Node>& vertices, vector<Edge>& edges)  {
        N = vertices.size();
        M = 0;
        Vertices = vertices;
        Adj.clear();

        // initialize the LstNode for all vertices
        for (int i = 0; i < vertices.size(); i++){
            LstNode adj;
            Adj.push_back(adj);
            (Adj.back().adjVec).clear();
        }

        // construct directed graph by adding edges to it
        for (unsigned i = 0; i < edges.size(); i++)  {
            add(edges[i].start, edges[i].end, edges[i].weight, edges[i].label);
             }
    }

 bool Graph::operator() (int i,int j) { return (i<j);}
//================================================================
//                  Public methods
//================================================================

//Ajoute une arête (de type Edge) au graphe.
void Graph::add( Edge &e) {
    Neighbor node(e.end,e.weight,e.label);
   (Adj[e.start].adjVec).push_back(node);
   M++;
    }

//Ajoute une arête (définie par ses quatre paramètres) au graphe.
void Graph::add(int start,int end, int weight, char* labelEdge)   {
        Neighbor node(end,weight,labelEdge);
       (Adj[start].adjVec).push_back(node);
       M++;
    }

//Ajoute un sommet (de type Node) au graphe.
void Graph::add(Node &v) {
   Vertices.push_back(v);
       LstNode A;
       Adj.push_back(A);
   N++;
    }

//Ajoute un sommet (défini par ses trois paramètres) au graphe.
void Graph::add(int val,int weight, string label)   {
        Node node(val,weight,label);
       Vertices.push_back(node);
       LstNode A;
       Adj.push_back(A);
       N++;
    }

/* Permet de savoir s'il existe une arête du sommet start au sommet end.
 * Si c'est le cas, revoit un pointeur vers l'arête dans la liste d'adjacence.
 * Sinon, revoit le pointeur NULL.
 */
Neighbor* Graph::link(int start, int end) {
    if (start >= N) {
        return NULL;
    }
    for (int index = 0; index < Adj[start].adjVec.size(); index++){
        if (Adj[start].adjVec[index].val == end) {
            return &Adj[start].adjVec[index];
        }
    }
    return NULL;
}

/* WARNING
 * Dans ce model de graphe, on suppose que l'on ne peut pas retirer des
 * sommets et donc la valeur d'un sommet est tout simplement son indice
 * dans le vecteur des sommets. La fonction suivante sert juste à
 * renvoyer la liste d'adjacence du n-ème sommet du graphe.
 * En cas d'ajout de la possibilité de retirer des sommets, cette
 * fonction sera bien à préciser car ne renverra pas forcément les
 * voisins d'un sommet de valeur donnée!!!!
 */ 
vector<Neighbor>* Graph::Neighbors(int n){
    return &(Adj[n].adjVec);
}


//Donne comme poids à toutes les arêtes la taille du label du sommet d'arrivée
void Graph::weighing(){
    for (int index = 0; index < N; index++){
        for (vector<Neighbor>::iterator it = Neighbors(index)->begin(); it != Neighbors(index)->end(); ++it){
            it->weight = Vertices[it->val].label.size();
        }
    }
}

//Calcule l'arête complémentaire de s.
void comp(char* s, char* r){
    if (s[0] == 'F')
    {
        r[1]='R';
    } else {
        r[1] = 'F';
    }
    if (s[1] == 'F')
    {
        r[0]='R';
    } else {
        r[0] = 'F';
    }
}

// BFS qui stocke les arêtes vues dans e. Le premier argument permet de dire dans quel sens on va
void Graph::BFS(int r, vector<Edge>& e ,vector<Neighbor*> &aVoir,vector<int> &vu){
    if (aVoir.size() == 0){ //Cas de terminaison, on a terminé le BFS
        return;
    }
    Neighbor* node = aVoir.front();
    aVoir.erase(aVoir.begin());
    if (find(vu.begin(),vu.end(),node->val) != vu.end()){ //Cas où le sommet a été vu par le BFS
        return BFS(r,e,aVoir,vu);
    }
    vu.push_back(node->val);
    for (vector<Neighbor>::iterator it = Neighbors(node->val)->begin(); it != Neighbors(node->val)->end(); ++it){
        //On boucle sur ses voisins
        if (it->label[0] == node->label[1] && find(vu.begin(),vu.end(),it->val) == vu.end()){ 
            //Cas où l'arrêt est bien valide et sommet non vu avant, ce voisin est rajouté dans la file des visites
            aVoir.push_back(&(*it));
            if (r){
                Edge edge(node->val,it->val,0,it->label);
                e.push_back(edge);
            } else {
                char R[3];
                comp(it->label,R);
                Edge edge(it->val,node->val,0,R);
                e.push_back(edge);
            }
        }
    }
    return BFS(r,e,aVoir,vu); //Sinon, on continue sans prendre en compte le sommet.
}



// BFS qui teste si les sommets vérifient bien une condition (par une fonction)
void Graph::BFS_func(int threshold ,vector<Neighbor*> &aVoir,vector<int> &vu){
    Neighbor* node ;
    while (aVoir.size() != 0){ //Cas de terminaison, on a terminé le BFS
        node = aVoir.front();
        aVoir.erase(aVoir.begin());
        if (vu[node->val]){ //Cas où le sommet a été vu par le BFS
            continue;
        }
        vu[node->val]=1;

        for (vector<Neighbor>::iterator it = Neighbors(node->val)->begin(); it != Neighbors(node->val)->end(); ++it){
            //On boucle sur ses voisins
            if (Vertices[it->val].weight >= threshold and vu[it->val]==0){
                //Cas où l'arrêt est bien valide et sommet non vu avant, ce voisin est rajouté dans la file des visites
            }
        }
    }
    return;
}


/* /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ 
 *   Petite approximation ici, je suppose que dès lors que l'on 
 *   voit un sommet par le sens forward, on ne le recroisera pas
 *   par le sens reverse ! Cela semble peu probable d'arriver et
 *   l'impact semble peu important (car BFS) mais on sous-estime
 *   la vraie valeur.
 *  /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ /!\ 
 *
 * Fonction inductive permettant d'effectuer un BFS, en mettant à jour
 * les données pour les appels inductifs suivants. De manière générale,
 * cette fonction sert à compter le nombre de sommets dans une boule de
 * rayon fixé, autour d'un sommet donné. 
 */
int Graph::BFSCount(vector<int> &rayons, int acc,vector<Neighbor*> &aVoir,vector<int> &vu){
    if (aVoir.size() == 0){ //Cas de terminaison, on a terminé le BFS
        return acc;
    }
    int rayon = rayons.front();
    rayons.erase(rayons.begin());
    Neighbor* node = aVoir.front();
    aVoir.erase(aVoir.begin());
    if (find(vu.begin(),vu.end(),node->val) != vu.end()){ //Cas où le sommet a été vu par le BFS
        return BFSCount(rayons,acc,aVoir,vu);
    }
    vu.push_back(node->val);
    if (node->weight <= rayon+40){ //Cas où le sommet est bien dans la boule
        for (vector<Neighbor>::iterator it = Neighbors(node->val)->begin(); it != Neighbors(node->val)->end(); ++it){
            //On boucle sur ses voisins
            if (it->label[0] == node->label[1] && find(vu.begin(),vu.end(),it->val) == vu.end()){ 
                //Cas où l'arrêt est bien valide et sommet non vu avant, ce voisin est rajouté dans la file des visites
                aVoir.push_back(&(*it));
                rayons.push_back(rayon+40-node->weight);
            }
        }
        return BFSCount(rayons,acc+1,aVoir,vu); //On traite les cas suivants, en prenant en compte le sommet
    }

    return BFSCount(rayons,acc,aVoir,vu); //Sinon, on continue sans prendre en compte le sommet.
}

//Permet de donner un poids à un sommet, correspondant aux nombres de sommets présents dans un rayon donné.
void Graph::weighingANode(int source, int rayon) {
    vector<Neighbor*> aVoir;
    vector<int> vu;
    vector<int> rayons;
    rayons.clear();
    aVoir.clear();
    vu.clear();
    vu.push_back(source);
    for (vector<Neighbor>::iterator it = Neighbors(source)->begin(); it != Neighbors(source)->end(); ++it){
        aVoir.push_back(&(*it));
        rayons.push_back(rayon);
    }
    Vertices[source].weight = BFSCount(rayons,1,aVoir,vu);
}

//Permet de donner un poids à tous les sommets du graphe.
void Graph::weighingAllNodes(int rayon) {
    for (vector<Node>::iterator it = Vertices.begin(); it != Vertices.end(); ++it){
        weighingANode(it->val, rayon);
    }
}

//================================================================
//                  Reading functions
//================================================================


//Lit un fichier contenant les sommets du graphe et les ajoute au vecteur seqs (réalisé par Pierre Peterlongo et Vincent Lacroix)
void read_node_file( ifstream &node_file, vector<Node>& seqs)
{
    seqs.clear();
    string line;
    int u,v;
    string p;
    int compt;
    while (getline(node_file, line)) {
        istringstream ss(line);
        string substr;
        compt=0;
        while (getline(ss, substr, '\t')) {
            if (compt == 0) {
                u = stoi(substr);
            }
            else if (compt == 1){
                p=substr;
            }
            compt++;
        }
        Node node(u,0,p);
        seqs.push_back(node);
    }
}

//Lit un fichier contenant les sommets du graphe et les ajoute au vecteur seqs
void read_node_file_weighted( ifstream &node_file,vector<Node>& seqs)
{
    seqs.clear();
    string line;
    int u,v;
    string p;
    int compt;
    while (getline(node_file, line)) {
        istringstream ss(line);
        string substr;
        compt=0;
        while (getline(ss, substr, '\t')) {
            if (compt == 0) {
                u = stoi(substr);
            }
            else if (compt == 1){
                p=substr;
            } else if (compt == 2){
                v = stoi(substr);
            }
            compt++;
        }
        Node node(u,v,p);
        seqs.push_back(node);
    }
}


//Lit un fichier d'arêtes et les ajoute au vecteur edges
void read_edge_file( ifstream &edge_file, vector<Edge>& edges ) {
    edges.clear();
    string line;
    int u,v;
    char p[3];
    int compt;
    while (getline(edge_file, line)) {
        istringstream ss(line);
        string substr;
        compt=0;
        while (getline(ss, substr, '\t')) {
            if (compt == 0) {
                u = stoi(substr);
            }
            else if (compt == 1){
                v = stoi(substr);
            } else if (compt == 2){
                strcpy( p, substr.c_str() );
            }
            compt++;
        }
        Edge e(u,v,0,p);
        edges.push_back(e);
    }
}


//================================================================
//                  Fonctions d'affichage
//================================================================



//Affiche les sommets d'un graphe
void printGraphVertices(Graph& G)
{
    vector<Node> V = G.Vertices;
    for (vector<Node>::iterator it = V.begin(); it != V.end(); ++it) {
        cout << it->val << "\t" << (string)it->label << "\t" << it->weight << "\n";
    }
    cout << endl;
}

//Affiche les sommets d'un graphe dans une autre sortie (dans un fichier par exemple)
void printGraphVertices(Graph& G,ofstream& output)
{
    vector<Node> V = G.Vertices;
    for (vector<Node>::iterator it = V.begin(); it != V.end(); ++it) {
        output << it->val << "\t" <<(string)it->label << "\t" << it->weight << "\n";
    }
}

//Affiche les sommets d'un graphe à partir d'un vecteur de Nodes
void printVertices(vector<Node>& V)
{
    for (vector<Node>::iterator it = V.begin(); it != V.end(); ++it) {
        cout << it->val << "\t" << (string)it->label << "\t" << it->weight << "\n";
    }
    cout << endl;
}

//Affiche les arêtes d'un graphe
void printGraphEdges(Graph& G)
{
    vector<LstNode> Adj = G.Adj;
    int N = G.N;
    for (int i = 0; i<N ; i++ ) {
        for (vector<Neighbor>::iterator it = Adj[i].adjVec.begin(); it != Adj[i].adjVec.end(); ++it) {
            cout << i << " to " << it->val << " of type " << (string)it->label << " and of weight " << it->weight << endl;
        }
    }
    cout << endl;
}

//Affiche les arêtes d'un graphe à partir d'un vecteur d'Edges
void printEdges(vector<Edge>& E)
{
    int M = E.size();
    Edge* e;
    for (int i = 0; i<M ; i++ ) {
        *e = E[i];
        cout << e->start << " to " << e->end << " of type " << (string)e->label << " and of weight " << e->weight << endl;
    }
    cout << endl;
}

//Affiche les arêtes d'un graphe à partir d'un vecteur d'Edges
void printEdges(vector<Edge>& E,ofstream& output)
{
    for (vector<Edge>::iterator it = E.begin(); it != E.end(); ++it) {
        output << it->start << "\t" << it->end << "\t" << (string)it->label << "\t" << it->weight << "\n";
    }
}
