#include <bits/stdc++.h>
#include <fstream>
#include <chrono>

using namespace std;

int main(){

    ifstream infile("data.txt"); //opening the data file
    
    int m; //number of edges
    infile >> m;

    vector<vector<int>> edges; //vector to get the edges

    //getting the edges
    for (int i=0; i<m; i++){ 
        int u, v;
        infile >> u >> v;
        edges.push_back({u, v});
    }

    infile.close(); //closing the data file

    //printing the vector edges
    for (auto &edge : edges) {
            cout << "[" << edge[0] << "," << edge[1] << "] ";
        }
        cout << endl;

    return 0;
}