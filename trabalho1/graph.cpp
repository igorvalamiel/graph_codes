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
    node* back;
    node* next;
};

/* Creating graph structure */
struct graph {
    /*Getting graph type [ matrix | list ]*/
    // if graph_type == 1 ~> List
    // if graph_type == 0 ~> Matrix
    bool graph_type = 1;

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

    /*Creating the basics structures*/
    vector <vector <bool>> matrix; // matrix
    vector <vector <int>> CC; // conected components
    vector <vector <int>> sizesCC; //sizes of each CC
    int quantCC = 0; // quantity of CC
    vector <int> G_list; //getting the degrees of each vertex
    vector <node*> Linklist; //creating the linked-list
    vector <node*> TailLL; //creating a vector to keep all the last itens of the linked list

    //-----------------------------------------------------------------------------------------------------------------------
    /*Executing the other functions to work properly*/
    graph(const vector<vector<int>>& edges, int num_vertex, int num_edges, bool gt = 1){
        graph_edges = edges;
        n = num_vertex;
        m = (int)graph_edges.size();
        graph_type = gt;

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

    //-----------------------------------------------------------------------------------------------------------------------
    /* Starting the graph */
    void start() {
        if (graph_type) {
            //initiating the G_list
            for (int i=0; i<n; i++) {
                G_list.push_back(0);
                node* aux = new node;
                aux->vertex = i;
                aux->next = nullptr;
                aux->back = nullptr;
                Linklist.push_back(aux);
                TailLL.push_back(aux);
            }
            
            //placing all the edges
            for (auto item : graph_edges){
                
                int a = item[0], b = item[1];

                // creating edge a -> b
                node* auxA = new node; auxA->vertex = b;
                auxA->next = Linklist[a];
                if (Linklist[a] != nullptr) Linklist[a]->back = auxA;
                Linklist[a] = auxA;

                // creating edge b -> a
                node* auxB = new node; auxB->vertex = a;
                auxB->next = Linklist[b];
                if (Linklist[b] != nullptr) Linklist[b]->back = auxB;
                Linklist[b] = auxB;

                // adding a degree to a and b
                G_list[a]++;
                G_list[b]++;
            }
            mem_graph = printMemoryUsage();
        } else {
            /*matrix estructure*/ 
            /*creating marix nxn with 0's*/
            for (int i=0; i<=n; i++){
                vector <bool> support;
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
            mem_graph = printMemoryUsage();
        }
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
        for (int i=0; i<n; i++) {Copy_G_list.push_back(G_list[i]);}

        sort(Copy_G_list.begin(), Copy_G_list.end()); //sorting the copy list
        
        //getting the median
        if (n % 2 == 1) {Medi_g = Copy_G_list[(n/2)+1];} //if the number of vertexes are even
        else {Medi_g = (Copy_G_list[(n/2)-1] + Copy_G_list[n/2]) / 2;} //if the number of vertexes are odd
    }

    //-----------------------------------------------------------------------------------------------------------------------
    //Implementing BFS
    vector <vector <int>> BFS(int s, bool diam_detect = false){

        if (graph_type){
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
                node* aux = new node; aux = Linklist[v]; //creating a auxiliar node
                for (int i=1; i<=G_list[v]; i++){ //for each node neighbor
                    int v_aux = aux->vertex; //getting vertex number
                    if (!visit_stats[v_aux]) { //if not visited
                        visit_stats[v_aux] = 1; //mark as visited
                        parent[v_aux] = v; //getting parent
                        level[v_aux] = level[v] + 1; //setting level
                        Q.push(v_aux); //placing in the queue
                    }
                    aux = aux->next; //getting next neighbor
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
        } else {
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
    }

    //-----------------------------------------------------------------------------------------------------------------------
    //Implementing DFS
    vector <vector <int>> DFS(int s, bool diam_detect = false){

        if (graph_type) {
            auto start_time = chrono::high_resolution_clock::now(); //getting initial time

            vector <bool> visit_stats(n+1, 0); //creating a vector to mark if the vertex was already visited
            stack <int> P; //creating the queue for getting the next item to be visited

            vector <int> parent(n+1, 0); //vector to register the parent of each vertex
            vector <int> level(n+1, 0); //vector to register the level of each vertex

            P.push(s); //placing s in the queue

            while (!P.empty()){
                int v = P.top();
                P.pop();
                if (!visit_stats[v]){
                    visit_stats[v] = 1;

                    // Coleta todos os vizinhos de v
                    vector<int> neighbors;
                    node* aux = Linklist[v];
                    while (aux != nullptr) {
                        neighbors.push_back(aux->vertex);
                        aux = aux->next;
                    }

                    // Ordena do MAIOR para o MENOR para garantir ordem crescente de exploração
                    sort(neighbors.rbegin(), neighbors.rend());

                    for (int w : neighbors) {
                        P.push(w);
                        if (!visit_stats[w]) {
                            parent[w] = v;
                            level[w] = level[v] + 1;
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
        } else {
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
    }

    //-----------------------------------------------------------------------------------------------------------------------
    /*Getting the distance between the vertex a & b (obs: the distance between two vertex )*/
    int dist(int a, int b){
        vector <vector <int>> bfs_res = BFS(a, true); //creating a vector to receive the BFS values
        return bfs_res[b][1];   
    }

    //-----------------------------------------------------------------------------------------------------------------------
    /*Getting the diameter of the graph*/
    void diameter(){
        if (quantCC == 1) {
            int highest_level = 0; //setting the counter
            vector <vector <int>> l = BFS(1, true); //doing the BFS
            for (auto i : l) {
                if (i[1] > highest_level) {highest_level = i[1];} //finding the biggest distance
            }
            diam = highest_level;
        }
    }

    //-----------------------------------------------------------------------------------------------------------------------
    /*Getting all connected components*/
    void ConctComp() {
        //Making atributtes empty
        CC.clear();
        sizesCC.clear();
        quantCC = 0;

        //Placing the first item because the vertex 0 doesn't exist
        CC.push_back({});
        sizesCC.push_back({0, 0});

        //Marking first vertex as visited
        vector<bool> visited(n+1, 0);

        if (graph_type) {
            for (int start = 1; start <= n; start++) {
                if (!visited[start]) {
                    quantCC++;
                    vector<int> CC_itens;
                    int ctng_CC = 0;

                    queue<int> Q;
                    Q.push(start);
                    visited[start] = true;

                    while (!Q.empty()) {
                        int v = Q.front(); Q.pop();
                        CC_itens.push_back(v);
                        ctng_CC++;
                        node* aux = Linklist[v];
                        while (aux != nullptr) {
                            int w = aux->vertex;
                            if (!visited[w]) {
                                visited[w] = true;
                                Q.push(w);
                            }
                            aux = aux->next;
                        }
                    }

                    CC.push_back(CC_itens);
                    sizesCC.push_back({ctng_CC, quantCC});
                }
            }
        } else {
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
        }

        sort(sizesCC.begin(), sizesCC.end());
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

    //-------------------------------------------------------------------------------------------------------------------------
    /*Creating output grafics*/
    void print(){
        if (graph_type) {
            node* aux;
            int cnt = 0;
            for (auto line : Linklist){
                aux = line;
                cout << cnt << " => ";
                while (aux != nullptr){
                    cout << aux->vertex << " ~> ";
                aux = aux->next;
                }
                cnt++;
                cout << '\n';
            }
        } else {
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

//-------------------------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------

int main() {

    //opening the data file
    ifstream infile("data.txt");

    //getting the number of lines
    int nlines; infile >> nlines;
    cout << nlines << '\n';

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

    //opening the output_data file
    ofstream outD("out_data.txt", std::ios::app);

    graph testL(edges, n, m);
    graph testM(edges, n, m, 0);

    //1st Question
    outD << "Questão 1\n";
    outD << "Lista 1: " << testL.mem_graph;
    outD << "Matriz 1: " << testM.mem_graph << '\n';

    outD << "=================================================\n";
    outD.close();

    return 0;
}
