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
    return "Memoria utilizada pelo processo: " + to_string(virtualMemUsedByMe / 1024) + " KB\n";
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

    vector <int> G_list; //getting the degrees of each vertex

    //-----------------------------------------------------------------------------------------------------------------------
    /* Starting the list */
    void start() {
        vector <node*> Linklist(m+1); //creating the linked-list
        for (int i=0; i<=m; i++) {G_list.push_back(0);} //initiating the G_list
        
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

            // adding a degree to a and b
            G_list[a]++;
            G_list[b]++;
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
        getinfo();
        cout << "getinfo ok\n";
        /*ConctComp();
        cout << "CC ok\n";
        diameter();
        cout << "diameter ok\n";*/
    }

    //-----------------------------------------------------------------------------------------------------------------------
    /*Getting all the information needed*/
    void getinfo() {
        G_min = n; //seting G_min for the max value (the biggest degree a vertex can have is n-1, that's why I settle it n)
        for (int i=1; i<=n; i++){
            double value = G_list[i];
            if (value < G_min) G_min = value; //getting lowest degree
            if (value > G_max) G_max = value; //getting highest degree
            G_med += value;
        }

        G_med = G_med / (double) n; //getting the medium degree

        //creating a copy of G_list to find the median
        vector <int> Copy_G_list; 
        for (int i=1; i<=n; i++) {Copy_G_list.push_back(G_list[i]);}

        sort(Copy_G_list.begin(), Copy_G_list.end()); //sorting the copy list
        
        //getting the median
        if (n % 2 == 1) {Medi_g = Copy_G_list[(n/2)+1];} //if the number of vertexes are even
        else {Medi_g = (Copy_G_list[(n/2)-1] + Copy_G_list[n/2]) / 2;} //if the number of vertexes are odd
    }

    //-----------------------------------------------------------------------------------------------------------------------
    //Implementing BFS
    vector <vector <int>> BFS(int s, bool diam_detect = false){

        auto start_time = chrono::high_resolution_clock::now(); //getting initial time

        vector <bool> visit_stats(n+1, 0); //creating a vector to mark if the vertex was already visited
        queue <int> Q; //creating the queue for getting the next item to be visited

        vector <int> parent(n+1, 0); //vector to register the parent of each vertex
        vector <int> level(n+1, 0); //vector to register the level of each vertex

        visit_stats[s] = 1; //marking s as visited
        Q.push(s); //placing s in the queue

        while (Q.size() > 0){ //While there is any item on the queue
            int v = Q.front(); //getting the head
            Q.pop(); //deleting the head
        }

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

    //Output model (it should appear in another file just like that)
    cout << "\nNumero de vertices: " << test.n << '\n';
    cout << "Numero de arestas: " << test.m << '\n';
    cout << "Grau minimo: " << test.G_min << '\n';
    cout << "Grau maximo: " << test.G_max << '\n';
    cout << "Grau medio: " << test.G_med << '\n';
    cout << "Mediana de grau: " << test.Medi_g << '\n';
    cout << test.mem_graph << '\n';
    //if (test.diam < 0){cout << "Diametro do Grafo: infinito\n";}
    //else {cout << "Diametro do Grafo: " << test.diam << "\n";}
    //cout << "Componentes Conexas (" << test.quantCC << " CC's)\n";

    return 0;
}