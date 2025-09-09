#include <bits/stdc++.h>
#include <fstream>
#include <string>
#include <chrono>
#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#endif

using namespace std;

string printMemoryUsage() {
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS_EX pmc;
    GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc));
    SIZE_T virtualMemUsedByMe = pmc.PrivateUsage;
    return "Memória utilizada pelo processo: " + to_string(virtualMemUsedByMe / 1024) + " KB\n";
#else
    return "Memória: unavailable\n";
#endif
}

struct graph {
    vector<vector<int>> graph_edges; 
    int n = 0; // número de vértices
    int m = 0; // número de arestas
    int G_min = 0, G_max = 0, G_med = 0, Medi_g = 0;
    double dt = 0;
    int diam = 0;
    string mem_graph;

    vector<vector<int>> matrix;   // matriz de adjacência
    vector<int> G_list;           // graus

    // componentes conexas
    vector<vector<int>> CC;
    vector<pair<int,int>> sizesCC;
    int quantCC = 0;

    // construtor
    graph(const vector<vector<int>>& edges, int num_vertex, bool diam_exact = false) {
        graph_edges = edges;
        n = num_vertex;
        m = (int)graph_edges.size();
        start();
        getinfo();
        ConctComp();
        diameter(diam_exact);
    }

    void start() {
        // cria a matriz nxn já zerada de uma vez
        matrix.assign(n+1, vector<int>(n+1, 0));
        G_list.assign(n+1, 0);

        // insere arestas
        for (auto &e : graph_edges) {
            int a = e[0], b = e[1];
            matrix[a][b] = 1;
            matrix[b][a] = 1;
            G_list[a]++; G_list[b]++;
        }
        mem_graph = printMemoryUsage();
    }

    void getinfo() {
        if (n <= 0) return;
        G_min = INT_MAX; G_max = 0;
        long long soma = 0;
        vector<int> graus; graus.reserve(n);
        for (int i=1; i<=n; i++) {
            int g = G_list[i];
            graus.push_back(g);
            soma += g;
            G_min = min(G_min, g);
            G_max = max(G_max, g);
        }
        G_med = (int)(soma / n);
        sort(graus.begin(), graus.end());
        if (n % 2 == 1) Medi_g = graus[n/2];
        else Medi_g = (graus[n/2 - 1] + graus[n/2]) / 2;
    }

    // BFS na matriz
    vector<pair<int,int>> BFS(int s, bool diam_detec = false) {
        auto t0 = chrono::high_resolution_clock::now();
        vector<int> vis(n+1,0), parent(n+1,0), level(n+1,-1);
        queue<int> Q;
        vis[s]=1; level[s]=0; Q.push(s);
        while(!Q.empty()) {
            int v=Q.front(); Q.pop();
            for (int i=1;i<=n;i++) {
                if (matrix[v][i] && !vis[i]) {
                    vis[i]=1;
                    parent[i]=v;
                    level[i]=level[v]+1;
                    Q.push(i);
                }
            }
        }
        vector<pair<int,int>> ret(n+1);
        for (int i=0;i<=n;i++) ret[i]={parent[i], level[i]};
        auto t1 = chrono::high_resolution_clock::now();
        dt = chrono::duration<double,milli>(t1-t0).count();
        return ret;
    }

    // DFS na matriz
    vector<pair<int,int>> DFS(int s) {
        auto t0 = chrono::high_resolution_clock::now();
        vector<int> vis(n+1,0), parent(n+1,0), level(n+1,-1);
        stack<int> st; st.push(s); level[s]=0;
        while(!st.empty()) {
            int u=st.top(); st.pop();
            if(!vis[u]) {
                vis[u]=1;
                for(int j=n;j>=1;j--) {
                    if(matrix[u][j] && !vis[j]) {
                        parent[j]=u;
                        level[j]=level[u]+1;
                        st.push(j);
                    }
                }
            }
        }
        vector<pair<int,int>> ret(n+1);
        for (int i=0;i<=n;i++) ret[i]={parent[i], level[i]};
        auto t1 = chrono::high_resolution_clock::now();
        dt = chrono::duration<double,milli>(t1-t0).count();
        return ret;
    }

    // distância entre a e b
    int dist(int a,int b) {
        auto bfs=BFS(a,true);
        return bfs[b].second;
    }

    // diâmetro: exato (muito lento) ou approx (2 BFS)
    void diameter(bool exact=false) {
        if(!exact) {
            int start=1;
            auto bfs1=BFS(start,true);
            int u=max_element(bfs1.begin(), bfs1.end(),
                              [](auto&a,auto&b){return a.second<b.second;})->first;
            auto bfs2=BFS(u,true);
            diam=0;
            for(auto &p:bfs2) diam=max(diam,p.second);
        } else {
            diam=0;
            for(int i=1;i<=n;i++){
                auto bfs=BFS(i,true);
                for(auto &p:bfs) diam=max(diam,p.second);
            }
        }
    }

    void ConctComp() {
        CC.clear(); sizesCC.clear(); quantCC=0;
        vector<int> vis(n+1,0);
        for(int s=1;s<=n;s++){
            if(!vis[s]) {
                quantCC++;
                CC.push_back({});
                sizesCC.push_back({0,quantCC});
                stack<int> st; st.push(s); vis[s]=1;
                while(!st.empty()){
                    int u=st.top(); st.pop();
                    CC.back().push_back(u);
                    sizesCC.back().first++;
                    for(int v=1;v<=n;v++){
                        if(matrix[u][v] && !vis[v]){
                            vis[v]=1; st.push(v);
                        }
                    }
                }
            }
        }
        sort(sizesCC.begin(),sizesCC.end(),
             [](auto&a,auto&b){return a.first>b.first;});
    }
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    ifstream infile("data.txt");
    int nlines; infile>>nlines;
    vector<vector<int>> edges; edges.reserve(nlines);
    int biggest=0;
    for(int i=0;i<nlines;i++){
        int a,b; infile>>a>>b;
        edges.push_back({a,b});
        biggest=max({biggest,a,b});
    }
    infile.close();

    ofstream outD("out_data.txt",ios::app);
    graph G(edges,biggest,false); // false = diâmetro approx
    
    //Output model (it should appear in another file just like that)
    outD<<"\nNumero de vertices: "<<G.n<<"\n";
    outD<<"Numero de arestas: "<<G.m<<"\n";
    outD<<"Grau minimo: "<<G.G_min<<"\n";
    outD<<"Grau maximo: "<<G.G_max<<"\n";
    outD<<"Grau medio: "<<G.G_med<<"\n";
    outD<<"Mediana de grau: "<<G.Medi_g<<"\n";
    outD<<G.mem_graph<<"\n";
    outD<<"Diametro do Grafo: "<<G.diam<<"\n";
    outD<<"Componentes Conexas ("<<G.quantCC<<" CC's)\n";

    for (auto [tam, idx] : G.sizesCC) {
        outD << "CC " << idx+1 << ": (" << tam << " vertices) ~ [ ";
        for (auto v : G.CC[idx]) outD << v << " ";
        outD << "]\n";
    }


    outD<<"=================================================\n";
    outD.close();
}
