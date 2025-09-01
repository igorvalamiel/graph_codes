#include <bits/stdc++.h>
#include <fstream>
#include <string>
#include <chrono>

using namespace std;

void createFile(int n){

    string fileName = "testFile" + to_string(n) + ".txt";

    ofstream testFile(fileName);

    if (testFile.is_open()){
        testFile << "Testando essa merda aqui.\n";
        testFile << "I'm working on this code for " << n << " days";
    }
}

int main(){
    
    vector<int> l = {1,2,3,4};

    for (auto i : l){
        createFile(i);
    }

    return 0;
}
