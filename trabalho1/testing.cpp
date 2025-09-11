#include <bits/stdc++.h>
#include <fstream>
#include <string>
#include <chrono>

using namespace std;

vector<int> DFS(vector<vector<bool>> graph, int vertex, float lenght){

    vector<int> marks(lenght, -1);
    stack<int> vertexstack;
    marks[vertex - 1] = 0;
    vertexstack.push(vertex - 1);
    while (vertexstack.empty() == false){
        int actual = vertexstack.top();
        vertexstack.pop();
        for (int i = 0; i < lenght; i++){
            if (graph[actual][i] != 0){
                if (marks[i] == -1){
                    marks[i] = marks[actual] + 1;
                    vertexstack.push(i);
                }
            }
        }
    }
    return marks;

}

int main(){
    
    return 0;
}
