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
    for (int i = 0; i<V.size() ; i++){
        if (E[i].label[0] == 'F' and E[i].label[1] == 'F'){
            E_duppli.push_back(Edge(E[i].start,E[i].end,E[i].weight,aretFF));
        } else if (E[i].label[0] == 'R' and E[i].label[1] == 'R'){
            E_duppli.push_back(Edge(n+E[i].start,n+E[i].end,E[i].weight,aretFF));
        } else if (E[i].label[0] == 'F' and E[i].label[1] == 'R'){
            E_duppli.push_back(Edge(E[i].start,n+E[i].end,E[i].weight,aretFF));
        } else if (E[i].label[0] == 'R' and E[i].label[1] == 'F'){
            E_duppli.push_back(Edge(n+E[i].start,E[i].end,E[i].weight,aretFF));
        }
    }


    Graph G(V,E_duppli);
    if(argc >= 6 and argv[4][1]=='k' ){
        G.kmer = stoi(argv[5]);
    }

    G.weighing();
    G.weighingAllNodesGraphDuppli(atoi(argv[3]));

    if(argc == 6 and argv[4][1]=='o'){
        string prefix = argv[5];
        ofstream output;
        output.open(prefix);
        ofstream outputE;
        outputE.open(prefix.substr(0,prefix.size()-5)+"edges");
        printGraphVertices(G,output);
        printEdgesBcalm(E_duppli,outputE);
        outputE.close();
        output.close();
    } else if(argc == 8 and argv[6][1]=='o'){
        string prefix = argv[7];
        ofstream output;
        output.open(prefix);
        ofstream outputE;
        outputE.open(prefix.substr(0,prefix.size()-5)+"edges");
        printGraphVertices(G,output);
        printEdgesBcalm(E_duppli,outputE);
        outputE.close();
        output.close();
    } else {
        printGraphVertices(G);
    }
    return 0;

}