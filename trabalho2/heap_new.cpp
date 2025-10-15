#include <iostream>
#include <vector>
#include <unordered_map>
#include <functional> // Para std::less
#include <stdexcept>

// TKey: Tipo da chave/identificador do nó (ex: int para um vértice, string para um nome)
// TPriority: Tipo da prioridade (ex: int, double)
// Compare: Functor de comparação para definir a ordem (padrão para min-heap)
template<
    typename TKey, 
    typename TPriority, 
    typename Compare = std::greater<TPriority> // Padrão: maior prioridade = menor valor
>
class priory_queue_update {
private:
    // Estrutura interna para armazenar os nós no heap
    struct Node {
        TKey key;
        TPriority priority;
    };

    std::vector<Node> heap;                     // O heap em si, armazenado em um vetor
    std::unordered_map<TKey, size_t> position_map; // Mapeia a chave para seu índice no vetor 'heap'
    Compare comp;                               // Objeto de comparação

    // --- Funções Auxiliares do Heap ---

    // Obtém o índice do pai
    size_t parent(size_t i) const {
        return (i - 1) / 2;
    }

    // Obtém o índice do filho esquerdo
    size_t leftChild(size_t i) const {
        return 2 * i + 1;
    }

    // Obtém o índice do filho direito
    size_t rightChild(size_t i) const {
        return 2 * i + 2;
    }

    // Troca dois nós no heap e atualiza seus mapeamentos de posição
    void swapNodes(size_t i, size_t j) {
        // Atualiza o mapa de posições
        position_map[heap[i].key] = j;
        position_map[heap[j].key] = i;

        // Troca os nós no vetor
        std::swap(heap[i], heap[j]);
    }

    // Move um nó para cima para manter a propriedade do heap (usado no push/update)
    void siftUp(size_t i) {
        // Enquanto o nó não for a raiz e tiver maior prioridade que seu pai
        while (i > 0 && comp(heap[parent(i)].priority, heap[i].priority)) {
            swapNodes(i, parent(i));
            i = parent(i);
        }
    }

    // Move um nó para baixo para manter a propriedade do heap (usado no pop/update)
    void siftDown(size_t i) {
        size_t minIndex = i;

        while (true) {
            size_t l = leftChild(i);
            size_t r = rightChild(i);
            
            // Verifica se o filho esquerdo existe e tem maior prioridade
            if (l < heap.size() && comp(heap[minIndex].priority, heap[l].priority)) {
                minIndex = l;
            }
            
            // Verifica se o filho direito existe e tem maior prioridade
            if (r < heap.size() && comp(heap[minIndex].priority, heap[r].priority)) {
                minIndex = r;
            }
            
            // Se o nó atual já tem a maior prioridade, está no lugar certo
            if (i == minIndex) {
                break;
            }

            // Caso contrário, troca com o filho de maior prioridade e continua
            swapNodes(i, minIndex);
            i = minIndex;
        }
    }

public:
    // --- Funções Públicas ---

    // Verifica se a fila está vazia
    bool empty() const {
        return heap.empty();
    }

    // Retorna o número de elementos na fila
    size_t size() const {
        return heap.size();
    }
    
    // Verifica se uma chave existe na fila
    bool contains(const TKey& key) const {
        return position_map.count(key);
    }

    // Insere um novo elemento ou atualiza um existente
    void push(const TKey& key, const TPriority& priority) {
        if (contains(key)) {
            update(key, priority);
            return;
        }

        // Adiciona o novo nó no final do heap
        heap.push_back({key, priority});
        size_t index = heap.size() - 1;
        position_map[key] = index;

        // Move o nó para cima para sua posição correta
        siftUp(index);
    }

    // Atualiza a prioridade de uma chave existente
    // A chave DEVE existir na fila
    void update(const TKey& key, const TPriority& new_priority) {
        if (!contains(key)) {
            // Ou poderia lançar uma exceção:
            // throw std::out_of_range("A chave não existe na fila de prioridade.");
            return;
        }

        size_t index = position_map[key];
        TPriority old_priority = heap[index].priority;
        heap[index].priority = new_priority;

        // Se a nova prioridade é "melhor" (menor para um min-heap), tenta subir o nó.
        // Se for "pior", tenta descer o nó.
        if (comp(old_priority, new_priority)) {
            siftUp(index);
        } else {
            siftDown(index);
        }
    }

    // Retorna o elemento de maior prioridade sem removê-lo
    const TKey& top() const {
        if (empty()) {
            throw std::out_of_range("A fila de prioridade está vazia.");
        }
        return heap.front().key;
    }

    // Remove e retorna o elemento de maior prioridade
    void pop() {
        if (empty()) {
            throw std::out_of_range("A fila de prioridade está vazia.");
        }
        
        // Remove a chave do elemento raiz do mapa
        TKey key_to_remove = heap.front().key;
        position_map.erase(key_to_remove);

        // Se for o último elemento, apenas o remova
        if (heap.size() == 1) {
            heap.pop_back();
            return;
        }

        // Move o último elemento para a raiz
        heap.front() = heap.back();
        heap.pop_back();
        
        // Atualiza a posição do nó que foi movido para a raiz
        position_map[heap.front().key] = 0;

        // Restaura a propriedade do heap movendo a nova raiz para baixo
        siftDown(0);
    }
};

