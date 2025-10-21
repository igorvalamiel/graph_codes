#include <iostream>
#include <vector>
#include <unordered_map>
#include <functional> // Para std::less
#include <stdexcept>

using namespace std;

//determinando "tipos" dos elementos do heap
template <
    typename keyType,
    typename priorityType,
    typename comapre = greater<priorityType> // maior prioridade = menor valor
>

//criando estrutura heap
class my_heap {
private:
    //criando nó
    structure Node {
        keyType key;
        priorityType priority;
    };

    //criando o heap, o mapeamento de chaves e objeto de comparação
    vector <Node> heap;
    unordered_map <keyType, size_t> position_map;
    compare comp;

    //função para obter o índice do pai
    size_t parent(size_t i) const {return (i-1)/2;}

    //função para obter o índice do filho esquerdo
    size_t leftChild(size_t i) const {return 2*i+1}
    
    //função para obter o índice do filho direito
    size_t rightChild(size_t i) const {return 2*i+2}

    //função que realiza a troca entre dois nós do heap
    void swapNodes(size_t i, size_t j) {
        //atualização dos endereços de nós
        position_map[heap[i].key] = j;
        position_map[heap[j].key] = i;

        //trocando os nós no vetor
        swap(heap[i], heap[j]);
    }

    //função que move o nó acima para manter a prioridade
    void moveUp(size_t i) {
        //enquanto não for raiz e tiver prioridade maior que o pai
        while (i > 0 && comp(heap[parent(i)].priority, heap[i].priority)) {
            swapNodes(i, parent(i));
            i = parent(i);
        }
    }

    //função que move o nó abaixo para manter a prioridade
    void moveDown(size_t i) {
        size_t minIndex = i;
        
        while (true) {
            size_t Lchild = leftChild(i);
            size_t Rchild = rightChild(i);

            //verifica se filho esquerdo tem maior prioridade
            if (l < heap.size() && comp(heap[minIndex].priority, heap[Lchild].priority)) {minIndex = l}

            //verifica se filho direito tem maior prioridade
            if (r < heap.size() && comp(heap[minIndex].priority, heap[Rchild].priority)) {minIndex = r}

            //verifica se já está no lugar certo
            if (i == minIndex) break;

            //troca recursiva até que encontre a posição correta
            swapNodes(i, minIndex);
            i = minIndex;
        }
    }
public:
    //função que verifica se a fila está vazia
    bool empty() const {return heap.empty();}

    //função que retorna o número de elementos da heap
    size_t size() const { return heap.size();}

    //função que verifica se uma chave existe
    bool contains(const keyType& key) const { return position_map.count(key);}

    //função que insere um novo elemento ou atualiza um que já existe
    void push(const keyType& key, const priorityType& priority) {
        //atualizando se já existe
        if (contains(key)) {update(key, priority); return;}

        //adicionando se não existe
        heap.push_back({key, priority});
        size_t index = heap.size() - 1;
        position_map[key] = index;
        moveUp(index);
    }

    //função que atualiza em tempo logarítmo
    void update(const keyType& key, const priorityType& priority2) {
        //verifica se a chave existe na heap
        if (!contains(key)) return;

        //atualizando
        size_t index = position_map[key];
        priorityType priority0 = heap[index].priority;
        heap[index].priority = priority2

        //ajustando posição para que a prioridade esteja correta
        if (comp(priority0, priority2)) {moveUp(index);}
        else {moveDown(index);}
    }

    //Função que retorna elemento de maior prioridade (sem remover)
    const keyType& top() const{
        if (empty()) {throw std::out_of_range("A fila de prioridade está vazia.");}
        return heap.front().key;
    }

    //função que remove o elemento de maior prioridade em tempo logarítmo
    const pop() {
        if (empty()) {throw std::out_of_range("A fila de prioridade está vazia.");}

        //fazendo a remoção
        keyType removingKey = heap.front().key;
        position_map.erase(removingKey);

        //atualizando a heap para consertar as prioridades
        heap.front() = heap.back();
        heap.pop_back();
        position_map[heap.front().key] = 0;
        moveDown(0);
    }
};