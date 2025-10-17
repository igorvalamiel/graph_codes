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
    float weight;
    node* back;
    node* next;
};

/* Creating graph structure */
struct graph {
    bool graph_type = 1; // 1 = lista, 0 = matriz

    vector<vector<float>> graph_edges;
    int n = 0; //number of vertices
    int m = 0; //number of edges
    int G_min = 0, G_max = 0, Medi_g = 0;
    double G_med = 0;
    double dt = 0;
    int diam = -1;
    string mem_graph;

    vector<vector<bool>> matrix;
    vector<vector<int>> CC;
    vector<vector<int>> sizesCC;
    int quantCC = 0;
    vector<int> G_list;
    vector<node*> Linklist;
    vector<node*> TailLL;
    vector<vector<float>> weight_matrix; // matriz de pesos
    float inf = numeric_limits<float>::infinity();

//-----------------------------------------------------------------------------------------------------------------------

    graph(const vector<vector<float>>& edges, int num_vertex, int num_edges, bool gt = 1) {
        graph_edges = edges;
        n = num_vertex;
        m = (int)graph_edges.size();
        graph_type = gt;

        start();
        cout << "Start ok\n";
        getinfo();
        cout << "getinfo ok\n";
        
        // print grafo no final
        printGraph();
    }
//-----------------------------------------------------------------------------------------------------------------------
/* Starting the graph */
    void start() {
        if (graph_type) {
            //initiating the G_list
            for (int i=0; i<=n; i++) {
                G_list.push_back(0); //adding a vertex in the degree's list
                node* aux = new node;
                aux->vertex = i;
                aux->next = nullptr;
                aux->back = nullptr;
                Linklist.push_back(nullptr);
                TailLL.push_back(nullptr);
            }      
            //placing all the edges
            for (auto item : graph_edges){   
                int a = (int)item[0], b = (int)item[1]; 
                float w = item[2]; //including weight
               
                // creating edge a -> b
                node* auxA = new node;
                auxA->vertex = b;
                auxA->next = Linklist[a];
                auxA->weight = w;
                if (Linklist[a] != nullptr) Linklist[a]->back = auxA;
                Linklist[a] = auxA;
                
                // creating edge b -> a
                node* auxB = new node;
                auxB->vertex = a;
                auxB->next = Linklist[b];
                auxB->weight = w;
                if (Linklist[b] != nullptr) Linklist[b]->back = auxB;
                Linklist[b] = auxB;

                // adding a degree to a and b
                G_list[a]++;
                G_list[b]++;
            }
            mem_graph = printMemoryUsage();
        } else {
            /*matrix estructure*/ 
            /*creating matrix nxn with 0's*/
            for (int i=0; i<=n; i++){
                vector <bool> support;
                vector <float> support_weight;
                for (int j=0; j<=n; j++){support.push_back(0);support_weight.push_back(0.0);}
                matrix.push_back(support);
                weight_matrix.push_back(support_weight);
                G_list.push_back(0);
            }
            /*Placing edges*/
            for (auto item : graph_edges){
                matrix[item[0]][item[1]] = 1;
                matrix[item[1]][item[0]] = 1;
                weight_matrix[item[0]][item[1]] = item[2];
                weight_matrix[item[1]][item[0]] = item[2];
                G_list[item[0]] += 1;
                G_list[item[1]] += 1;
            }
            mem_graph = printMemoryUsage();
        }
    }
//-----------------------------------------------------------------------------------------------------------------------
    void getinfo() {
        G_min = INT_MAX;
        G_max = INT_MIN;
        G_med = 0;

        for (int i = 1; i <= n; i++) {
            int value = G_list[i];
            G_min = min(G_min, value);
            G_max = max(G_max, value);
            G_med += value;
        }
        G_med /= (double)n;

        vector<int> Copy_G_list(G_list.begin() + 1, G_list.end());
        sort(Copy_G_list.begin(), Copy_G_list.end());

        if (n % 2 == 1)
            Medi_g = Copy_G_list[n / 2];
        else
            Medi_g = (Copy_G_list[n / 2 - 1] + Copy_G_list[n / 2]) / 2;
    }

//-----------------------------------------------------------------------------------------------------------------------
    void printGraph() {
    
    if (graph_type) {
        for (int i = 0; i <= n; i++) {
            cout << i << " -> bao ";
            node* current = Linklist[i];
            while (current != nullptr) {
                cout << current->vertex << "(" << current->weight << ") ";
                current = current->next;
            }
            cout << "\n";
        }
    } else {
        for (int i = 1; i <= n; i++){
                cout << "|  ";
                for (int j = 1; j <= n; j++) { 
                    if (matrix[i][j]) { 
                        cout << weight_matrix[i][j] << "  ";
                    }
                    else{
                        cout << "." << "   ";
                    }
                }
                cout << "|\n";
            }
        
        }
         
    }


//-----------------------------------------------------------------------------------------------------------------------

    /*Implementing Dijkstra with vector (it returns the spanning tree)*/
    vector<vector<float>> Dijkstra_vector(int s){
        float inf = numeric_limits<float>::infinity();
        vector<float> dist(n + 1, inf); // distancia entre s e cada vértice i
        vector<int> parent(n + 1, -1);
        vector<int> discovered;
        vector<int> unexplored;
        vector<int> level(n+1, -1);
        vector<int> explored(n+1,-1);
        
        for (int i = 0; i <= n; i++){
            unexplored.push_back(i);
        }

        dist[s] = 0;
        level[s]=0;
        discovered.push_back(s);
        unexplored[s] = -1;

        if (graph_type){

            // seleciona o vértice no conjunto de descobertos
            while(!discovered.empty()){
                int u = discovered[0];
                
                node* current = Linklist[u];
                while(current != nullptr){
                    int v = current->vertex;
                    float w = current->weight;

// mais facil ir por indice
                    if (unexplored[v] == v){ // verificando se tá nos inesplorados
                        discovered.push_back(v); // adicionando aos descobertos
                        unexplored[v] = -1; //marcando v como "não inesplorado"
                    }
                    if (explored[v] == -1) { //aqui eu preciso verificar se v tá nos explorados
                        if (dist[v] > dist[u] + w){
                            dist[v] = dist[u] + w;
                            parent[v] = u;
                            level[v] = level[u]+ 1;
                        }
                    }

                    current = current->next;
                }
            
                discovered.erase(discovered.begin());
                explored[u] = u;
            }
        } else {
            // seleciona o vértice no conjunto de descobertos
            while(!discovered.empty()){
                int u = discovered[0];

                for (int v = 0; v <= n; v++){
                if (matrix[u][v]){
                    float w = weight_matrix[u][v];

                    // só adiciona v em discovered se estiver em unexplored
                    if (unexplored[v] == v){
                        discovered.push_back(v);
                        unexplored[v] = -1;
                    }

                    // só atualiza a distancia se ainda não estiver explorado pra n dar xabu
                    if (explored[v] == -1){
                        if (dist[v] > dist[u] + w){
                            dist[v] = dist[u] + w;
                            parent[v] = u;
                            level[v] = level[u]+ 1;

                        }
                    }
                }
            }

            // aqui o u já foi explorado
            discovered.erase(discovered.begin());
            explored[u] = u;
            }

            }
        
    vector<vector<float>> tree(n + 1, vector<float>(3));
    for (int i = 0; i <= n; i++){
        tree[i][0] = parent[i];
        tree[i][1] = dist[i];
        tree[i][2] = level[i];
    }
            
    return tree;

    }
    void printTree(const vector<vector<float>>& tree) {
    cout << "Vértice\tPai\tDistância\tNível\n";
    for (int i = 0; i < tree.size(); i++) {
        cout << i << "\t" 
             << tree[i][0] << "\t" 
             << tree[i][1] << "\t\t" 
             << tree[i][2] << "\n";
    }
}


};
    


//-------------------------------------------------------------------------------------------------------------------------
int main() {
    // Abrir arquivo
    ifstream infile("C:/Users/Julia/Desktop/teste_grafo.txt");

    // Ler número de vértices
    int n;
    infile >> n;

    // Criar vetor de arestas
    vector<vector<float>> edges; 
    float a, b, w;

    // Ler todas as arestas do arquivo
    while (infile >> a >> b >> w) {
        edges.push_back({(float)a, (float)b, w});
    }
    infile.close();

    // Número de arestas
    int m = (int)edges.size();

    graph test(edges, n, m,0);

    vector<vector<float>> tree = test.Dijkstra_vector(1);
    test.printTree(tree);




    // Salvar memória usada
    ofstream outD("out_data.txt", ios::app);
    outD << "Memória grafo: " << test.mem_graph;
    outD.close();

    return 0;
}


