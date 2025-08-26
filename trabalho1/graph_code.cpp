#include <bits/stdc++.h>
#include <chrono>

using namespace std;

struct graph {
    /*Edges pairs*/
    vector <vector <int>> graph_edges;

    /*graph's info*/
    int n = 0; //number of vectors
    int m = graph_edges.size(); //number of edges
    int G_min, G_max, G_med, Med_g; //maximum, minimum, medium and mediana of the degrees
    int dt; //execution time to create the structure

    /*Setting type of structure*/
    int ListOrMat = 0;

    /*Creating the basics structures*/
    vector <vector <int>> matrix; // matrix
    vector<int> linklist; //linked list

    void start() {
        /*matrix estructure*/ 
        if (ListOrMat){
            /*creating marix nxn with 0's*/

            //OBS VER COM PROFESSOR SE É MELHOR COMEÇAR NO 0 OU NO 1
            for (int i=0; i<=n; i++){
                vector <int> support;
                for (int j=0; j<=n; j++){support.push_back(0);}
                matrix.push_back(support);
            }

            /*Placing edges*/
            for (auto item : graph_edges){
                matrix[item[0]][item[1]] = 1;
                matrix[item[1]][item[0]] = 1;
            }
        }

        /*Linked list estructure*/
        else {
            int nada;
        }
        
    }

    /*Creating output grafics*/
    void print(){
        // if is a matrix
        if (ListOrMat){
            for (auto line : matrix){
                cout << "|  ";
                for (auto item : line) {
                    cout << item << "   ";
                }
                cout << "|\n";
            }
        }
    }
};


int main() {

    /*getting the number of lines*/
    int nlines; cin >> nlines;

    /*creating the a vector of vectors to keep all edges information*/
    vector <vector <int>> edges;

    /*Looking for the biggest number to see the number of vertex*/
    int biggest = 0;

    /*getting all the edges of the graph*/
    for (int i=0; i<nlines; i++){
        int a, b; cin >> a >> b;
        vector <int> line = {a, b};
        edges.push_back(line);
        if (a > biggest) biggest = a;
        if (b > biggest) biggest = b;
    }

    /*
    TO CREATE A GRAPH USING A MATRIX DATA STRUCTURE, YOU SHOULD CALL:
        graph.ListOrMat = 1;
    */


    graph test;
    test.ListOrMat = 1;
    test.graph_edges = edges;
    test.n = biggest;

    test.start();
    test.print();


/*
    //printing the lines
    cout << '\n';
    for (int i=0; i<nlines; i++){
        cout << edges[i][0] << "~>" << edges[i][1] << '\n';
    }

    return 0;
*/

}
