/*
 * Ce programme permet de loader un graphe de De Bruijn
 * et de calculer les poids des arêtes et des sommets à
 * un rayon donné en entrée.
 */

#include <list>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "graph.h"


// graph implementation
int main(int argc, char** argv)
{
    if (argc!=4 and argc!=6 and argc!=8){
        cout << "Expected use of this program: \n\n\t" <<argv[0] << " file.nodes file.edges radius -k kmer -o output.txt\n" << endl;
        return 0;
    }

    vector<Edge> E;
    vector<Node> V;

    char* nodesPath = argv[1];
    char* edgesPath = argv[2];

    ifstream edges(edgesPath, std::ios::binary);
    ifstream nodes(nodesPath, std::ios::binary);

    read_edge_file(edges,E);
    read_node_file(nodes,V);


    vector<Edge> E_duppli;
    int n = V.size();
    for (int i = 0; i<n ; i++){
        V.push_back(Node(n+i,V[i].weight, reverse_complement(V[i].label)));
    }
    char aretFF[3] = {'F', 'F'};
    char aretFR[3] = {'F', 'R'};
    char aretRF[3] = {'R', 'F'};
    char aretRR[3] = {'R', 'R'};
    for (int i = 0; i<n ; i++){
        if (E[i].label[0] == 'F' and E[i].label[1] == 'F'){
            E_duppli.push_back(Edge(E[i].start,E[i].end,E[i].weight,aretFF));
        } else if (E[i].label[0] == 'R' and E[i].label[1] == 'R'){
            E_duppli.push_back(Edge(n+E[i].start,n+E[i].end,E[i].weight,aretRR));
        } else if (E[i].label[0] == 'F' and E[i].label[1] == 'R'){
            E_duppli.push_back(Edge(E[i].start,n+E[i].end,E[i].weight,aretFR));
        } else if (E[i].label[0] == 'R' and E[i].label[1] == 'F'){
            E_duppli.push_back(Edge(n+E[i].start,E[i].end,E[i].weight,aretRF));
        }
    }


    Graph G(V,E_duppli);
    if(argc >= 6 and argv[4][1]=='k' ){
        G.kmer = stoi(argv[5]);
    }

    G.weighing();
    G.weighingAllNodesGraphDuppli(atoi(argv[3]));

    if(argc == 6 and argv[4][1]=='o'){
        ofstream output;
        output.open(argv[5]);
        printGraphVertices(G,output);
        output.close();
    } else if(argc == 8 and argv[6][1]=='o'){
        ofstream output;
        output.open(argv[7]);
        printGraphVertices(G,output);
        output.close();
    } else {
        printGraphVertices(G);
    }
    return 0;

}