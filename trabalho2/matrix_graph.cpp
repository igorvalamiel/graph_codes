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

struct graph {
    /*Edges pairs*/
    vector <vector <float>> graph_edges;

    /*graph's info*/
    int n = 0; //number of vectors
    int m = 0; //number of edges
    int G_min = 0, G_max = 0, Medi_g = 0; //maximum, minimum, medium and median of the degrees
    double G_med = 0;
    double dt = 0; //execution time to create the structure. Only not 0 when start() executed
    int diam = -1;
    bool w = false; //verifying if the graph is weightened
    string mem_graph;

    /*Creating the basics structures*/
    vector <vector <bool>> matrix; // matrix
    vector <vector <int>> CC; // conected components
    vector <vector <int>> sizesCC; //sizes of each CC
    int quantCC = 0; // quantity of CC

    //Support structures and variables
    vector <int> G_list;

    //-----------------------------------------------------------------------------------------------------------------------
    void start() {
        /*matrix estructure*/ 
        /*creating marix nxn with 0's*/
        for (int i=0; i<=n; i++){
            vector <bool> support;
            for (int j=0; j<=n; j++){support.push_back('.');}
            matrix.push_back(support);
            G_list.push_back(0);
        }

        /*Placing edges*/
        for (auto item : graph_edges){
            matrix[item[0]][item[1]] = item[2];
            matrix[item[1]][item[0]] = item[2];
            G_list[item[0]] += 1;
            G_list[item[1]] += 1;
        }
        mem_graph = printMemoryUsage();
    }

    //-----------------------------------------------------------------------------------------------------------------------
    void getinfo() {
        G_min = n; //seting G_min for the max value (the biggest degree a vertex can have is n-1, that's why I settle it n)
        for (int i=1; i<=n; i++){
            double value = G_list[i];
            if (value < G_min) G_min = value; //getting lowest degree
            if (value > G_max) G_max = value; //getting highest degree
            G_med += value;
        }
        
        G_med = G_med / (double) n; //getting medium degree
        
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

        while (Q.size() > 0){ //while there is any item on the queue
            int v = Q.front(); //getting the head
            Q.pop(); //deleting the head

            for (int i=1; i<=n; i++){ //the matrix representation uses matrix[v][i] to say if i is a neighbor of v
                if (matrix[v][i] != 0){ //if they are neighbors
                    if (visit_stats[i] == 0){ //if not visited yet
                        visit_stats[i] = 1; //mark as visited
                        parent[i] = v; //getting parent
                        level[i] = level[v] + 1; //setting level
                        Q.push(i); //placing the neighbor in the queue
                    }
                }
            }
        }

        vector <vector <int>> ret;
        for (int i=0; i<=n; i++){
            vector <int> aux = {parent[i], level[i]};
            ret.push_back(aux);
        }

        auto end_time = chrono::high_resolution_clock::now(); //getting ending time
        chrono::duration<double,std::milli> duration = end_time - start_time;
        dt = duration.count(); //em ms

        createFile("BFS", ret, diam_detect, dt);

        return ret;
    }

    //-----------------------------------------------------------------------------------------------------------------------
    //implemneting DFS
    vector <vector <int>> DFS(int s, bool diam_detect = false){

        auto start_time = chrono::high_resolution_clock::now(); //getting initial time

        vector <bool> visit_stats(n+1, 0); //creating a vector to mark if the vertex was already visited
        stack <int> P; //creating the stack for getting the next item to be visited

        vector <int> parent(n+1, 0); //vector to register the parent of each vertex
        vector <int> level(n+1, 0); //vector to register the level of each vertex

        P.push(s); //adding s to the stack

        while (!P.empty()){
            int u = P.top(); //getting the highest element
            P.pop(); //removing the highest element
            if (visit_stats[u] == 0){ //verifying if u was already visited
                visit_stats[u] = 1; //marking u as visited
                for (int j=n; j>=1; j--){ //looking for the next neighbor
                    if (matrix[u][j] != 0){
                        if (visit_stats[j] == 0){ //if the neighbor wasn't visited
                            parent[j] = u; //setting parent
                            level[j] = level[u] + 1; //setting level
                            P.push(j); //putting neighbor in the stack
                        }
                    }
                }
            }
        }

        vector <vector <int>> ret;
        for (int i=0; i<=n; i++){
            vector <int> aux = {parent[i], level[i]};
            ret.push_back(aux);
        }

        auto end_time = chrono::high_resolution_clock::now(); //getting ending time
        chrono::duration<double,std::milli> duration = end_time - start_time;
        dt = duration.count(); //em ms

        createFile("DFS", ret, diam_detect, dt);

        return ret;
    }

    //-----------------------------------------------------------------------------------------------------------------------
    /*Getting the distance between the vertex a & b (obs: the distance between two vertex )*/
    int dist(int a, int b){
        vector <vector <int>> bfs_res = BFS(a, true); //creating a vector to receive the BFS values
        return bfs_res[b][1];
    }

    //-----------------------------------------------------------------------------------------------------------------------
    /*Getting the diameter of the graph*/
    void diameter(bool complete = 0){
        if (complete){
            int big = 0; //setting the counter
            for (int i = 1; i<=n; i++){
                vector <vector <int>> l = BFS(1, true); //doing the BFS
                for (auto j : l) {
                    if (j[1] > big) {big = j[1];} //finding the biggest distance
                }
            }
        } else if (quantCC == 1) {
            int big = 0; //setting the counter
            vector <vector <int>> l = BFS(1, true); //doing the BFS
            for (auto i : l) {
                if (i[1] > big) {big = i[1];} //finding the biggest distance
            }
            int double_diam = 0; //setting the counter
            vector <vector <int>> l2 = BFS(big, true); //doing the BFS
            for (auto i : l) {
                if (i[1] > double_diam) {double_diam = i[1];} //finding the biggest distance
            }
            diam = double_diam;
        }
    }

    //-----------------------------------------------------------------------------------------------------------------------
    /*Getting all conected components*/
    void ConctComp() {
        //Making atributtes empty
        CC.clear();
        sizesCC.clear();
        quantCC = 0;

        //Placing the first item because the vertex 0 doesn't exist
        CC.push_back({});       
        sizesCC.push_back({0, 0}); //{size, idCC}

        //Marking first vertex as visited
        vector<int> visited(n+1, 0);

        for (int start = 1; start <= n; start++) {
            if (!visited[start]) {
                //Checking if found a new component
                quantCC++;
                CC.push_back({});
                sizesCC.push_back({0, quantCC});

                //creating the stack to get the Conected Components
                stack<int> P;
                P.push(start);
                visited[start] = 1;

                while (!P.empty()) { //While there is no more vertex in the CC
                    int u = P.top();
                    P.pop();

                    //Add u to current component
                    CC.back().push_back(u);
                    sizesCC.back()[0]++;

                    //Explore neighbors
                    for (int v = 1; v <= n; v++) {
                        if (matrix[u][v] != 0 && !visited[v]) {
                            visited[v] = 1;
                            P.push(v);
                        }
                    }
                }
            }
        }

        //sorting the vector to print in decreasing order
        sort(sizesCC.begin(), sizesCC.end());
    }

    //-----------------------------------------------------------------------------------------------------------------------
    /*Creating output grafics*/
    void print(){
        for (int i=0; i<n; i++){
            cout << "|  ";
            for (int j=0; j<n; j++) {
                cout << matrix[i][j] << "   ";
            }
            cout << "|\n";
        }
    }

    //-----------------------------------------------------------------------------------------------------------------------
    /*Executing the other functions to work properly*/
    graph(const vector<vector<float>>& edges, int num_vertex, int num_edges, bool weightened){
        graph_edges = edges;
        n = num_vertex;
        m = (int)graph_edges.size();
        w = weightened;

        //As soon as the structure graph is called, all these functions are also called
        start();
        cout << "Start ok\n";
        getinfo();
        cout << "getinfo ok\n";
        ConctComp();
        cout << "CC ok\n";
        diameter();
        cout << "diameter ok\n";
    }

    //-------------------------------------------------------------------------------------------------------------------------
    /*Creating a function to create and/or modify a file*/
    void createFile(string name, vector <vector <int>> s, bool get_diam, int t){

        if (!get_diam) {
            if (name == "BFS"){
                ofstream testFile("bfs_output.txt", std::ios::app);
                testFile << "BFS ~   ";

                testFile << "Levels: [ ";
                for (auto par : s){
                    testFile << par[1] << ' ';
                } testFile << "]    ";
                testFile << "|   Parents: [ ";
                for (auto par : s){
                    testFile << par[0] << ' ';
                } testFile << "]";
                testFile << "   |   Runtime: " << t << "ms\n";

                testFile.close();
            } else {
                ofstream testFile("dfs_output.txt", std::ios::app);
                testFile << "DFS ~   ";

                testFile << "Levels: [ ";
                for (auto par : s){
                    testFile << par[1] << ' ';
                } testFile << "]    ";
                testFile << "|   Parents: [ ";
                for (auto par : s){
                    testFile << par[0] << ' ';
                } testFile << "]";
                testFile << "   |   Runtime: " << t << "ms\n";

                testFile.close();
            }
        }
    }
};

//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------
int main() {

    //opening the data file
    ifstream infile("data.txt");

    //getting the number of lines
    int nlines; infile >> nlines;
    cout << nlines << '\n';

    //creating the a vector of vectors to keep all edges information
    vector <vector <float>> edges;
    string line;
    bool weightened = false;

    //getting n and m
    int n = nlines;
    int m = 0;

    //stopping point
    int last1, last2, lastw;

    //getting all the edges of the graph
    while (true){
        getline(infile, line);
        stringstream ss(line);

        float a, b, w; ss >> a >> b;

        if (ss >> w) weightened = true;

        if (weightened) {
            if (a == last1 && b == last2 && w == lastw){break;}
            else {
                vector <float> line = {a, b, w};
                edges.push_back(line);
                last1 = a;
                last2 = b;
                lastw = w;
                m++;
            }
        } else {
            if (a == last1 && b == last2){break;}
            else {
                vector <float> line = {a, b};
                edges.push_back(line);
                last1 = a;
                last2 = b;
                m++;
            }
        }
        
    }

    //closing the data file
    infile.close();

    //opening the output_data file
    ofstream outD("out_data.txt", std::ios::app);

    graph test(edges, n, m, weightened);

    test.print();
    cout << endl;

    //Output model (it should appear in another file just like that)
    cout << "\nNumero de vertices: " << test.n << '\n';
    cout << "Numero de arestas: " << test.m << '\n';
    cout << "Grau minimo: " << test.G_min << '\n';
    cout << "Grau maximo: " << test.G_max << '\n';
    cout << "Grau medio: " << test.G_med << '\n';
    cout << "Mediana de grau: " << test.Medi_g << '\n';
    cout << test.mem_graph << '\n';
    if (test.diam < 0){cout << "Diametro do Grafo: infinito\n";}
    else {cout << "Diametro do Grafo: " << test.diam << "\n";}
    cout << "Componentes Conexas (" << test.quantCC << " CC's)\n";

    for (int i=test.quantCC; i>0; i--) {
        vector <int> vecCC = test.sizesCC[i];
        cout << "CC " << vecCC[1] << ": (" << vecCC[0] << " vertices) ~ [ ";
        for (auto item : test.CC[vecCC[1]]){cout << item << " ";}
        cout << "]\n";
    }

    /*
    vector <vector <int>> a = test.BFS(1);
    cout << "Pai[10] = " << a[10][0] << '\n';
    cout << "Pai[20] = " << a[20][0] << '\n';
    cout << "Pai[30] = " << a[30][0] << '\n';
    */

    cout << "=================================================\n";
    outD.close();

    return 0;
}