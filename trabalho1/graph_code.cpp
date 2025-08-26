#include <bits/stdc++.h>
#include <chrono>

using namespace std;

struct graph {
    /*Edges pairs*/
    vector <vector <int>> graph_edges;

    /*graph's info*/
    int n = 0; //number of vectors
    int m = graph_edges.size(); //number of edges
    int G_min = 0, G_max = 0, G_med = 0, Med_g = 0; //maximum, minimum, medium and median of the degrees
    double dt = 0; //execution time to create the structure. Only not 0 when start() executed

    /*Setting type of structure*/
    int ListOrMat = 0;

    /*Creating the basics structures*/
    vector <vector <int>> matrix; // matrix
    vector<int> linklist; //linked list

    //Support structures and variables
    vector <int> G_list;

    void start() {
        /*matrix estructure*/ 
        if (ListOrMat){
            auto start_time = chrono::high_resolution_clock::now(); //getting initial time
            /*creating marix nxn with 0's*/
            for (int i=0; i<=n; i++){
                vector <int> support;
                for (int j=0; j<=n; j++){support.push_back(0);}
                matrix.push_back(support);
                G_list.push_back(0);
            }

            /*Placing edges*/
            for (auto item : graph_edges){
                matrix[item[0]][item[1]] = 1;
                matrix[item[1]][item[0]] = 1;
                G_list[item[0]] += 1;
                G_list[item[1]] += 1;
            }
            auto end_time = chrono::high_resolution_clock::now(); //getting ending time
            chrono::duration<double,std::milli> duration = end_time - start_time;
            dt = duration.count(); //em ms
        }

        /*Linked list estructure*/
        else {
            auto start_time = chrono::high_resolution_clock::now(); //getting initial time
            int nada;
            auto end_time = chrono::high_resolution_clock::now(); //getting ending time
            chrono::duration<double,std::milli> duration = end_time - start_time;
            dt = duration.count(); //em ms
        }
    }

    void getinfo() {
        G_min = n; //seting G_min for the max value (the biggest degree a vertex can have is n-1, that's why I settle it n)
        for (int i=1; i<=n; i++){
            int value = G_list[i];
            if (value < G_min) G_min = value; //getting lowest degree
            if (value > G_max) G_max = value; //getting highest degree
            G_med += value;
        }
        G_med = G_med / n; //getting medium degree
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

    //calling comands
};

//-------------------------------------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------------------------------------
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
    test.getinfo();
    test.print();
    cout << "\n" << test.dt << " ms\n";
    cout << test.G_max << ' ' << test.G_min << ' ' << test.G_med << "\n";

/*
    //printing the lines
    cout << '\n';
    for (int i=0; i<nlines; i++){
        cout << edges[i][0] << "~>" << edges[i][1] << '\n';
    }

    return 0;
*/

}
