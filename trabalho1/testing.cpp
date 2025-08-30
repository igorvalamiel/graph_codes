#include <bits/stdc++.h>
#include <fstream>
#include <chrono>

using namespace std;

int main(){

    ifstream infile("data.txt");
    string line;

    while (getline(infile, line)){
        cout << line << endl;
    }

    return 0;
}