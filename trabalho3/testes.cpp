#include <bits/stdc++.h>

using namespace std;

void printV(vector <vector <int>> v) {
    for (auto i : v) {
        cout << "{ ";
        for (auto j : i) {
            cout << j << ' ';
        } cout << "}, ";
    }
    cout << '\n';
}

int main() {
    vector <vector <int>> edges;

    edges.push_back({2,1});
    edges.push_back({3,2});
    edges.push_back({0,-3});
    edges.push_back({4,0});

    printV(edges);

    sort(edges.begin(), edges.end());

    printV(edges);

    return 0;
}