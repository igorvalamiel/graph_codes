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

int main(){
    //opening the data file
    ifstream infile("data.txt");
}