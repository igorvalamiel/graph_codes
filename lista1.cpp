#include <bits/stdc++.h>
#include <chrono>
using namespace std;

//funcao criada para printar o vetor
void print_vec(vector <int> item) {
    cout << '[';
    for (int i : item) {
        cout << i << ", ";
    }
    cout << "]\n";
}

//funcao recursiva criada para realizar as permutacoes
void permute(vector <int> v, int ini, int end) {
    //verificando se o inicio eh igual o final
    if (ini == end) {
        print_vec(v);
        return;
    }

    //trocando as posicoes e fazendo a recursao
    for (int i=ini; i<=end; i++){
        swap(v[ini], v[i]);
        permute(v, ini+1, end);
        swap(v[ini], v[i]);
    }

}

//codigo principal
int main() {
    int a;
    cin >> a;
    vector <int> V = {};

    //criando o vetor com todos os valores de 1 a 'a'
    for (int i=1; i<=a; i++) {
        V.push_back(i);
    }
    
    //começando a contagem do tempo
    auto start = chrono::high_resolution_clock::now();
    //chamando a funcao de permutacao
    permute(V, 0, a-1);
    //terminando a contagem do tempo
    auto end = chrono::high_resolution_clock::now();

    //calculando o tempo total de execucao
    chrono::duration<double,std::milli> duration = end - start;
    cout << "\nO tempo total de execucao foi: " << duration.count() << "ms" << endl;

    return 0;
}
