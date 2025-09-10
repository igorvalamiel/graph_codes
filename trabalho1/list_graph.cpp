#include <bits/stdc++.h>
#include <fstream>
#include <string>
#include <chrono>
#include <windows.h> 
#include <psapi.h> //to get memory info

using namespace std;

string printMemoryUsage() {
    PROCESS_MEMORY_COUNTERS_EX pmc;
    GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));
    SIZE_T virtualMemUsedByMe = pmc.PrivateUsage;
    return "Memória utilizada pelo processo: " + to_string(virtualMemUsedByMe / 1024) + " KB\n";
}

/* Creating node structure */
struct node {
    int vertex;
    node* next;
};

/* Creating linked-list structure*/
struct graph {
    /*Edges pairs*/
    vector <vector <int>> graph_edges;

    /*graph's info*/
    int n = 0; //number of vectors
    int m = 0; //number of edges
    int G_min = 0, G_max = 0, Medi_g = 0; //maximum, minimum, medium and median of the degrees
    double G_med = 0;
    double dt = 0; //execution time to create the structure. Only not 0 when start() executed
    int diam = -1;
    string mem_graph;

    //-----------------------------------------------------------------------------------------------------------------------
    /* Starting the list */
    void start() {
        vector <node*> Linklist(m+1); //creating the linked-list

        //placing all the edges
        for (auto item : graph_edges){
            
            int a = item[0], b = item[1];

            // creating edge a -> b
            node* auxA = new node; auxA->vertex = b;
            auxA->next = Linklist[a];
            Linklist[a] = auxA;

            // creating edge b -> a
            node* auxB = new node; auxB->vertex = a;
            auxB->next = Linklist[b];
            Linklist[b] = auxB;
        }
        mem_graph = printMemoryUsage();
    }

    //-----------------------------------------------------------------------------------------------------------------------
    /*Executing the other functions to work properly*/
    graph(const vector<vector<int>>& edges, int num_vertex, int num_edges){
        graph_edges = edges;
        n = num_vertex;
        m = (int)graph_edges.size();

        //As soon as the structure graph is called, all these functions are also called
        start();
        cout << "Start ok\n";
        /*getinfo();
        cout << "getinfo ok\n";
        ConctComp();
        cout << "CC ok\n";
        diameter();
        cout << "diameter ok\n";*/
    }

};


int main(){
    //opening the data file
    ifstream infile("data.txt");

    //getting the number of lines
    int nlines; infile >> nlines;

    //creating the a vector of vectors to keep all edges information
    vector <vector <int>> edges;

    //getting n and m
    int n = nlines;
    int m = 0;

    //stopping point
    int last1, last2;

    //getting all the edges of the graph
    while (true){
        int a, b; infile >> a >> b;
        if (a == last1 && b == last2){break;}
        else {
            vector <int> line = {a, b};
            edges.push_back(line);
            last1 = a;
            last2 = b;
            m++;
        }
    }

    //closing the data file
    infile.close();


    graph test(edges, n, m);

    cout << test.n << '\n';
    cout << test.m << '\n';

    return 0;
}