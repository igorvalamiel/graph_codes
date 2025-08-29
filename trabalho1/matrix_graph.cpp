#include <bits/stdc++.h>
#include <chrono>

using namespace std;

struct graph {
    /*Edges pairs*/
    vector <vector <int>> graph_edges;

    /*graph's info*/
    int n = 0; //number of vectors
    int m = graph_edges.size(); //number of edges
    int G_min = 0, G_max = 0, G_med = 0, Medi_g = 0; //maximum, minimum, medium and median of the degrees
    double dt = 0; //execution time to create the structure. Only not 0 when start() executed

    /*Creating the basics structures*/
    vector <vector <int>> matrix; // matrix
    vector <vector <int>> CC; // conected components
    vector <int> sizesCC; //sizes of each CC
    int quantCC = CC.size(); // quantity of CC

    //Support structures and variables
    vector <int> G_list;

    //-----------------------------------------------------------------------------------------------------------------------
    void start() {
        /*matrix estructure*/ 
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

    //-----------------------------------------------------------------------------------------------------------------------
    void getinfo() {
        G_min = n; //seting G_min for the max value (the biggest degree a vertex can have is n-1, that's why I settle it n)
        for (int i=1; i<=n; i++){
            int value = G_list[i];
            if (value < G_min) G_min = value; //getting lowest degree
            if (value > G_max) G_max = value; //getting highest degree
            G_med += value;
        }
        
        G_med = G_med / n; //getting medium degree
        
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
    vector <vector <int>> BFS(int s){
        vector <int> visit_stats(n+1, 0); //creating a vector to mark if the vertex was already visited
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
        return ret;
    }

    //-----------------------------------------------------------------------------------------------------------------------
    //implemneting DFS
    vector <vector <int>> DFS(int s){
        vector <int> visit_stats(n+1, 0); //creating a vector to mark if the vertex was already visited
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
        return ret;
    }

    //-----------------------------------------------------------------------------------------------------------------------
    /*Getting the distance between the vertex a & b (obs: the distance between two vertex )*/
    int dist(int a, int b){
        vector <vector <int>> bfs_res = BFS(a); //creating a vector to receive the BFS values
        return bfs_res[b][1];
    }

    //-----------------------------------------------------------------------------------------------------------------------
    /*Getting the diameter of the graph*/
    int diameter(){
        int max = 0; //setting auxiliar variable to find the maximum (minimun) distances
        for (int i=1; i<=n; i++){ //running a BFS for each vertex
            vector <vector <int>> bfs_res = BFS(i);
            for (vector <int> v : bfs_res){ //finding the longest path in the BFS
                if (v[1] > max) max = v[1];
            }
        }
        return max;
    }

    //-----------------------------------------------------------------------------------------------------------------------
    /*Getting all conected components*/
    void ConctComp(){
        vector <int> visit_stats(n+1, 0); //creating a vector to mark if the vertex was already visited
        vector <int> cc_list(n+1, 0); //creating a vector
        int idCC = 0;

        for (int i=0; i<=n; i++){
            if (visit_stats[i] == 0){
                idCC++;
                stack <int> P;
                visit_stats[i] = 1;
                cc_list[i] = idCC;

                while (!P.empty()){
                    int t = P.top();
                    P.pop();
                    for (int j=n; j>=1; j--){
                        if (matrix[t][j] != 0){
                            if (visit_stats[j] == 0){
                                P.push(j);
                                visit_stats[j] = 1;
                                cc_list[j] = idCC;
                            }
                        }
                    }
                }
            }
        }

        quantCC = idCC;
        for (int k=0; k<=idCC; k++) {
            sizesCC.push_back(0);
            vector <int> aux;
            CC.push_back(aux);
        }
        for (int l=0; l<=idCC; l++) {
            sizesCC[cc_list[l]]++;
            CC[cc_list[l]].push_back(l);
        }
    }

    //-----------------------------------------------------------------------------------------------------------------------
    /*Creating output grafics*/
    void print(){
        for (auto line : matrix){
            cout << "|  ";
            for (auto item : line) {
                cout << item << "   ";
            }
            cout << "|\n";
        }
    }
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


    graph test;
    test.graph_edges = edges;
    test.n = biggest;

    
    test.start();
    test.getinfo();
    test.print();
    test.ConctComp();
    cout << "\nExecution Time: " << test.dt << " ms\n";
    cout << "Max G: " << test.G_max << "\nMin G: " << test.G_min << "\n";
    cout << "Medium G: " << test.G_med << "\nMedian G: " << test.Medi_g << "\n"; 
    /*
    for (int i=1; i<=5; i++){
        vector <vector <int>> dfs = test.BFS(i);
        cout << '\n';
        for (int i=1; i<=5; i++){
            cout << "vertex:" << i << "(level " << dfs[i][1] << ") ~> " << dfs[i][0] << '\n';
        }
        cout << "-----------------------------\n";
    }
    */

    for (int i=1; i<=5; i++){
        cout << test.dist(1,i) << ' ';
    }
    cout << '\n';

    cout << test.diameter() << "\n\n";

    cout << test.quantCC << '\n';
    for (int i=1; i<=test.quantCC; i++){
        cout << "CC " << i << " ~> Size: " << test.sizesCC[i] << '\n';
        for (auto item : test.CC[i]){
            cout << item << " ";
        } cout << '\n';
    }

    return 0;
}
