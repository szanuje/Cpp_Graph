#pragma once

#include <cstdint>
#include <vector>
#include <optional>
#include <limits>
#include <stack>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <map>
#include <iomanip>
#include <cmath>

using namespace std;

template <typename V, typename E>
class Graph
{
public:

    // DFS IT

    class IteratorDFS
    {
    private:
        IteratorDFS(Graph<V, E> *graph, std::size_t curr_vector_idx);
        IteratorDFS();

    public:
        IteratorDFS(const IteratorDFS&) = default;
        IteratorDFS& operator=(const IteratorDFS&) = default;
        IteratorDFS(IteratorDFS&&) = default;
        IteratorDFS& operator=(IteratorDFS&&) = default;

        friend class Graph<V, E>;

    private:
        Graph<V,E> *graph_ptr_;
        std::size_t idx_;
        vector<std::size_t> v;
        vector<bool> visited;
        std::size_t curr;
        stack<std::size_t> s;

    public:
        bool operator==(const IteratorDFS &vi) const;
        bool operator!=(const IteratorDFS &vi) const;
        IteratorDFS& operator++();
        IteratorDFS operator++(int);
        V& operator*() const;
        V* operator->() const;
    };

    // BFS IT

    class IteratorBFS
    {
    private:
        IteratorBFS(Graph<V, E> *graph, std::size_t curr_vector_idx);
        IteratorBFS();

    public:
        IteratorBFS(const IteratorBFS&) = default;
        IteratorBFS& operator=(const IteratorBFS&) = default;
        IteratorBFS(IteratorBFS&&) = default;
        IteratorBFS& operator=(IteratorBFS&&) = default;

        friend class Graph<V, E>;

    private:
        Graph<V,E> *graph_ptr_;
        std::size_t idx_;
        vector<std::size_t> v;
        vector<bool> visited;
        std::size_t curr;

    public:
        bool operator==(const IteratorBFS &vi) const;
        bool operator!=(const IteratorBFS &vi) const;
        IteratorBFS& operator++();
        IteratorBFS operator++(int);
        V& operator*() const;
        V* operator->() const;
    };

    // VERTICES IT

    class VerticesIterator
    {
    private:
        VerticesIterator(Graph<V, E> *graph, std::size_t current_vertex_id);

    public:
        VerticesIterator(const VerticesIterator&) = default;
        VerticesIterator& operator=(const VerticesIterator&) = default;
        VerticesIterator(VerticesIterator&&) = default;
        VerticesIterator& operator=(VerticesIterator&&) = default;

        friend class Graph<V, E>;

    private:
        Graph<V,E> *graph_ptr_;
        std::size_t idx_;

    public:
        bool operator==(const VerticesIterator &vi) const;
        bool operator!=(const VerticesIterator &vi) const;
        VerticesIterator& operator++();
        VerticesIterator operator++(int);
        V& operator*() const;
        V* operator->() const;
    };

    // EDGES IT

    class EdgesIterator
    {
    private:
        EdgesIterator(Graph<V, E> *graph, std::size_t nm_row, std::size_t nm_col);

    public:
        EdgesIterator(const EdgesIterator&) = default;
        EdgesIterator& operator=(const EdgesIterator&) = default;
        EdgesIterator(EdgesIterator&&) = default;
        EdgesIterator& operator=(EdgesIterator&&) = default;

        friend class Graph<V, E>;

    private:
        Graph<V,E> *graph_ptr_;
        std::size_t Vidx_;
        std::size_t Hidx_;

    public:
        bool operator==(const EdgesIterator &ei) const;
        bool operator!=(const EdgesIterator &ei) const;
        EdgesIterator& operator++();
        EdgesIterator operator++(int);
        E& operator*() const;
        E* operator->() const;
    };

public:
    Graph() = default;

private:
    vector<vector<optional<E>>> edges_;
    vector<V> vertices_;
    int num_of_edges_ = 0;

public:

    VerticesIterator insertVertex(const V &vertex_data);
    std::pair<EdgesIterator, bool> insertEdge(std::size_t vertex1_id, std::size_t vertex2_id, E edge);
    bool removeVertex(std::size_t vertex_id);
    bool removeEdge(std::size_t vertex1_id, std::size_t vertex2_id);
    bool edgeExist(std::size_t vertex1_id, std::size_t vertex2_id) const;
    std::size_t nrOfVertices() const;
    std::size_t nrOfEdges() const;
    void printNeighborhoodMatrix() const;
    VerticesIterator vertex(std::size_t vertex_id);
    EdgesIterator edge(std::size_t vertex1_id, std::size_t vertex2_id);
    VerticesIterator begin();
    VerticesIterator end();
    VerticesIterator beginVertices();
    VerticesIterator endVertices();
    EdgesIterator beginEdges();
    EdgesIterator endEdges();

    vector<std::size_t> BFS(std::size_t index) const;
    vector<std::size_t> DFS(std::size_t index) const;

    IteratorBFS beginBFS(std::size_t index);
    IteratorBFS endBFS();

    IteratorDFS beginDFS(std::size_t index);
    IteratorDFS endDFS();

    std::size_t findIndexOfVertice(const V &vertice);
    V getVertice(std::size_t idx) const { return vertices_[idx]; }
    double getDistance(std::size_t first, std::size_t second);

    std::pair<double, vector<std::size_t>> dijkstra(std::size_t begin, std::size_t end);
};

template <typename V, typename E>
typename Graph<V,E>::VerticesIterator Graph<V,E>::insertVertex(const V &vertex_data) {
    vertices_.push_back(vertex_data);
    vector<optional<E>> vec(vertices_.size());
    fill(vec.begin(), vec.end(), optional<E>());
    edges_.push_back(vec);
    if(edges_.size() > 1) {
        int k = edges_.size() - 2;
        while(k >= 0) {
            edges_[k].emplace_back(optional<E>());
            k--;
        }
    }
    return VerticesIterator(this, vertices_.size() - 1);
}

template <typename V, typename E>
std::pair<typename Graph<V,E>::EdgesIterator, bool> Graph<V,E>::insertEdge(std::size_t vertex1_id, std::size_t vertex2_id, E edge) {
    if(edges_[vertex1_id][vertex2_id]) return std::pair(EdgesIterator(this, vertex1_id, vertex2_id), false);
    edges_[vertex1_id][vertex2_id].emplace(edge);
    num_of_edges_++;
    return std::pair(EdgesIterator(this, vertex1_id, vertex2_id), true);
}

template <typename V, typename E>
bool Graph<V,E>::removeVertex(std::size_t vertex_id) {
    if(vertex_id > vertices_.size()) return false;
    else {
        vertices_.erase(vertices_.begin() + vertex_id);
        for(auto &i : edges_) {
            if(i[vertex_id]) num_of_edges_--;
            i.erase(i.begin() + vertex_id);
        }
        for(const auto &i : edges_[vertex_id]) {
            if(i) num_of_edges_--;
        }
        edges_.erase(edges_.begin() + vertex_id);
        return true;
    }
}

template <typename V, typename E>
bool Graph<V,E>::removeEdge(std::size_t vertex1_id, std::size_t vertex2_id) {
    if(vertex1_id > edges_.size() || vertex2_id > edges_.size() || !edges_[vertex1_id][vertex2_id]) return false;
    else {
        edges_[vertex1_id][vertex2_id].reset();
        num_of_edges_--;
        return true;
    }
}

template <typename V, typename E>
bool Graph<V,E>::edgeExist(std::size_t vertex1_id, std::size_t vertex2_id) const {
    return edges_[vertex1_id][vertex2_id];
}

template <typename V, typename E>
std::size_t Graph<V,E>::nrOfVertices() const {
    return vertices_.size();
}

template <typename V, typename E>
std::size_t Graph<V,E>::nrOfEdges() const {
    return num_of_edges_;
}

template <typename V, typename E>
void Graph<V,E>::printNeighborhoodMatrix() const {
    for(const auto &e : edges_) {
        for(const auto &i : e) {
            if(!i) cout << left << setw(6) << "0 ";
            else cout << left << setw(6) << i.value();
        }
        cout << endl;
    }
}

template <typename V, typename E>
typename Graph<V,E>::VerticesIterator Graph<V,E>::vertex(std::size_t vertex_id) {
    return VerticesIterator(*this, vertex_id);
}

template <typename V, typename E>
typename Graph<V,E>::EdgesIterator Graph<V,E>::edge(std::size_t vertex1_id, std::size_t vertex2_id) {
    return EdgesIterator(*this, vertex1_id, vertex2_id);
}

template <typename V, typename E>
typename Graph<V,E>::VerticesIterator Graph<V,E>::begin() { return beginVertices(); }

template <typename V, typename E>
typename Graph<V,E>::VerticesIterator Graph<V,E>::end() { return endVertices(); }

template <typename V, typename E>
typename Graph<V,E>::VerticesIterator Graph<V,E>::beginVertices() {
    return VerticesIterator(this, 0);
}

template <typename V, typename E>
typename Graph<V,E>::VerticesIterator Graph<V,E>::endVertices() {
    return VerticesIterator(this, vertices_.size());
}

template <typename V, typename E>
typename Graph<V,E>::EdgesIterator Graph<V,E>::beginEdges() {
    for(std::size_t i = 0; i < edges_.size(); i++) {
        for(std::size_t j = 0; j < edges_.size(); j++) {
            if(edges_[i][j]) {
                return EdgesIterator(this, i, j);
            }
        }
    }
    return EdgesIterator(this, edges_.size(), 0);
}

template <typename V, typename E>
typename Graph<V,E>::EdgesIterator Graph<V,E>::endEdges() {
    return EdgesIterator(this, edges_.size(), 0);
}

// VERTICES IT

template <typename V, typename E>
Graph<V,E>::VerticesIterator::VerticesIterator(Graph<V, E> *graph, std::size_t current_vertex_id) : graph_ptr_(graph), idx_(current_vertex_id) {}

template <typename V, typename E>
bool Graph<V,E>::VerticesIterator::operator==(const VerticesIterator &vi) const { return this->idx_ == vi.idx_; }

template <typename V, typename E>
bool Graph<V,E>::VerticesIterator::operator!=(const VerticesIterator &vi) const { return !(this->idx_ == vi.idx_); }

template <typename V, typename E>
typename Graph<V,E>::VerticesIterator& Graph<V,E>::VerticesIterator::operator++() {
    if(idx_ != graph_ptr_->nrOfVertices()) idx_++;
    return *this;
}

template <typename V, typename E>
typename Graph<V,E>::VerticesIterator Graph<V,E>::VerticesIterator::operator++(int) {
    VerticesIterator temp(*this);
    ++(*this);
    return temp;
}

template <typename V, typename E>
V& Graph<V,E>::VerticesIterator::operator*() const {
    return graph_ptr_->vertices_[idx_];
}

template <typename V, typename E>
V* Graph<V,E>::VerticesIterator::operator->() const {
    return &graph_ptr_->vertices_[idx_];
}

// EDGES IT

template <typename V, typename E>
Graph<V,E>::EdgesIterator::EdgesIterator(Graph<V, E> *graph, std::size_t nm_row, std::size_t nm_col) : graph_ptr_(graph), Vidx_(nm_row), Hidx_(nm_col) {}

template <typename V, typename E>
bool Graph<V,E>::EdgesIterator::operator==(const EdgesIterator &ei) const { return this->Hidx_ == ei.Hidx_ && this->Vidx_ == ei.Vidx_; }

template <typename V, typename E>
bool Graph<V,E>::EdgesIterator::operator!=(const EdgesIterator &ei) const { return this->Hidx_ != ei.Hidx_ || this->Vidx_ != ei.Vidx_; }

template <typename V, typename E>
typename Graph<V,E>::EdgesIterator& Graph<V,E>::EdgesIterator::operator++() {
    while(Vidx_ != graph_ptr_->edges_.size()) {
        //        cout << "[" << Vidx_  << "]" << "[" << Hidx_ << "]";
        if(Hidx_ == graph_ptr_->edges_.size() - 1) {
            Vidx_++;
            Hidx_ = 0;
        } else Hidx_++;
        if(Vidx_ == graph_ptr_->edges_.size()) break;
        if(graph_ptr_->edges_[Vidx_][Hidx_]) break;
    }
    return *this;
}

template <typename V, typename E>
typename Graph<V,E>::EdgesIterator Graph<V,E>::EdgesIterator::operator++(int) {
    EdgesIterator temp(*this);
    ++(*this);
    return temp;
}

template <typename V, typename E>
E& Graph<V,E>::EdgesIterator::operator*() const {
    return graph_ptr_->edges_[Vidx_][Hidx_].value();
}

template <typename V, typename E>
E* Graph<V,E>::EdgesIterator::operator->() const {
    return &graph_ptr_->edges_[Vidx_][Hidx_];
}


template <typename V, typename E>
std::size_t Graph<V,E>::findIndexOfVertice(const V &vertice) {
    for(std::size_t i = 0; i < vertices_.size(); i++) {
        if(vertices_[i] == vertice) return i;
    }
    return 0;
}

// TRAVERSALS

template <typename V, typename E>
vector<std::size_t> Graph<V,E>::BFS(std::size_t index) const {
    bool visited[vertices_.size()];
    memset(visited, false, sizeof(visited));
    vector<std::size_t> v;
    v.push_back(index);
    visited[index] = true;

    for(std::size_t i = 0; i < v.size(); i++){

        //cout << vertices_[ v[i] ] << ", ";

        for(std::size_t j = 0; j < edges_.size(); j++) {

            if(edges_[v[i]][j] && !visited[j]) {
                v.push_back(j);
                visited[j] = true;

            }
        }
    }

    return v;
}

template <typename V, typename E>
vector<std::size_t> Graph<V,E>::DFS(std::size_t index) const {
    bool visited[vertices_.size()];
    memset(visited, false, sizeof(visited));
    vector<std::size_t> v;
    stack<std::size_t> s;
    v.push_back(index);
    visited[index] = true;

    while(true) {

        while(true) {

            std::size_t size_before = v.size();
            std::size_t size_after = v.size();
            //cout << vertices_[ v[i] ] << ", ";
            bool go = false;

            for(std::size_t j = 0; j < edges_.size(); j++) {

                if(edges_[v[v.size() - 1]][j] && !visited[j]) {

                    if(!go) {
                        v.push_back(j);
                        visited[j] = true;
                        go = true;
                        size_after++;
                    }

                    if(go) {
                        s.push(j);
                    }
                }
            }

            if(size_before == size_after) break;
        }

        while(true) {
            if(s.empty()) return v;
            if(visited[s.top()]) s.pop();
            if(!s.empty()) {
                if(!visited[s.top()]) {
                    v.push_back(s.top());
                    visited[s.top()] = true;
                    s.pop();
                    break;
                }
            }
        }
    }
}

// BFS IT

template <typename V, typename E>
Graph<V,E>::IteratorBFS::IteratorBFS(Graph<V, E> *graph, std::size_t curr_vector_idx) :
    graph_ptr_(graph), curr(curr_vector_idx) {}

template <typename V, typename E>
bool Graph<V,E>::IteratorBFS::operator==(const IteratorBFS &vi) const { return this->curr == vi.curr; }

template <typename V, typename E>
bool Graph<V,E>::IteratorBFS::operator!=(const IteratorBFS &vi) const { return this->curr != vi.curr; }

template <typename V, typename E>
typename Graph<V,E>::IteratorBFS& Graph<V,E>::IteratorBFS::operator++() {

    //cout << this->graph_ptr_->vertices_[ v[curr] ] << "x, ";

    for(std::size_t j = 0; j < graph_ptr_->edges_.size(); j++) {

        if(graph_ptr_->edges_[v[curr]][j] && !visited[j]) {
            v.push_back(j);
            visited[j] = true;

        }
    }

    curr++;
    if(curr < v.size()) {
        idx_ = v[curr];
    } else curr = this->graph_ptr_->vertices_.size();

    return *this;
}

template <typename V, typename E>
typename Graph<V,E>::IteratorBFS Graph<V,E>::IteratorBFS::operator++(int) {
    IteratorBFS temp(*this);
    ++(*this);
    return temp;
}

template <typename V, typename E>
V& Graph<V,E>::IteratorBFS::operator*() const {
    return graph_ptr_->vertices_[idx_];
}

template <typename V, typename E>
V* Graph<V,E>::IteratorBFS::operator->() const {
    return &graph_ptr_->vertices_[idx_];
}

template <typename V, typename E>
typename Graph<V,E>::IteratorBFS Graph<V,E>::beginBFS(std::size_t index) {
    IteratorBFS b(this, 0);
    b.visited.resize(vertices_.size(), false);
    b.v.push_back(index);
    b.idx_ = index;
    b.visited[index] = true;
    return b;
}

template <typename V, typename E>
typename Graph<V,E>::IteratorBFS Graph<V,E>::endBFS() {
    IteratorBFS b(this, vertices_.size());
    return b;
}

// DFS IT

template <typename V, typename E>
Graph<V,E>::IteratorDFS::IteratorDFS(Graph<V, E> *graph, std::size_t curr_vector_idx) :
    graph_ptr_(graph), curr(curr_vector_idx) {}

template <typename V, typename E>
bool Graph<V,E>::IteratorDFS::operator==(const IteratorDFS &vi) const { return this->curr == vi.curr; }

template <typename V, typename E>
bool Graph<V,E>::IteratorDFS::operator!=(const IteratorDFS &vi) const { return this->curr != vi.curr; }

template <typename V, typename E>
typename Graph<V,E>::IteratorDFS& Graph<V,E>::IteratorDFS::operator++() {

    bool go = false;

    bool repeat = true;
    while(repeat) {
        std::size_t status_before = v.size();
        std::size_t status_after = v.size();

        while(true) {
            for(std::size_t j = 0; j < graph_ptr_->edges_.size(); j++) {

                if(graph_ptr_->edges_[v[v.size() - 1]][j] && !visited[j]) {

                    if(!go) {
                        v.push_back(j);
                        status_after++;
                        visited[j] = true;
                        go = true;
                    }

                    if(go) {
                        s.push(j);
                    }
                }
            }

            if(status_before == status_after) {
                break;
            }
            else {
                curr++;
                idx_ = v[curr];
                return *this;
            }
        }

        while(true) {
            if(s.empty()) {
                repeat = false;
                curr = graph_ptr_->vertices_.size();
                break;
            }
            if(visited[s.top()]) s.pop();

            if(!s.empty()) {
                if(!visited[s.top()]) {
                    v.push_back(s.top());
                    visited[s.top()] = true;
                    s.pop();
                    curr++;
                    idx_ = v[curr];
                    return *this;
                }
            }
        }
    }

    curr++;
    if(s.empty()) curr = graph_ptr_->vertices_.size();

    return *this;
}

template <typename V, typename E>
typename Graph<V,E>::IteratorDFS Graph<V,E>::IteratorDFS::operator++(int) {
    IteratorDFS temp(*this);
    ++(*this);
    return temp;
}

template <typename V, typename E>
V& Graph<V,E>::IteratorDFS::operator*() const {
    return graph_ptr_->vertices_[idx_];
}

template <typename V, typename E>
V* Graph<V,E>::IteratorDFS::operator->() const {
    return &graph_ptr_->vertices_[idx_];
}

template <typename V, typename E>
typename Graph<V,E>::IteratorDFS Graph<V,E>::beginDFS(std::size_t index) {
    IteratorDFS b(this, 0);
    b.visited.resize(vertices_.size(), false);
    b.v.push_back(index);
    b.idx_ = index;
    b.visited[index] = true;
    return b;
}

template <typename V, typename E>
typename Graph<V,E>::IteratorDFS Graph<V,E>::endDFS() {
    IteratorDFS b(this, vertices_.size());
    return b;
}

template <typename V, typename E>
double Graph<V,E>::getDistance(std::size_t first, std::size_t second){
    return edges_[first][second];
}

// DIJKSTRA

template <typename V, typename E>
std::pair<double, vector<std::size_t>> Graph<V,E>::dijkstra(std::size_t begin, std::size_t end) {

    vector<bool> visited;
    visited.assign(vertices_.size(), false);

    vector<double> sums;
    for(std::size_t i = 0; i < vertices_.size(); i++) {
        sums.push_back(std::numeric_limits<double>::infinity());
    }

    vector<std::size_t> result;
    double sum = 0;

    stack<std::size_t> remain;

    remain.push(begin);
    sums[begin] = 0;

    while(!remain.empty()) {

        visited[remain.top()] = true;
        std::size_t current = remain.top();
        remain.pop();

        for(std::size_t i = 0; i < edges_.size(); i++) {
            if(edges_[current][i] && !visited[i]) {
                if(sums[current] + edges_[current][i].value() < sums[i]) {
                    sums[i] = sums[current] + edges_[current][i].value();
                }
                remain.push(i);
            }
        }
    }

    if(isinf(sums[end])) {
        return std::make_pair(0, result);
    }

    std::size_t lowest = begin;
    visited.assign(visited.size(), false);

    result.push_back(lowest);
    visited[lowest] = true;

    while(lowest != end) {

        std::size_t current_lowest;
        double current = std::numeric_limits<double>::infinity();

        for(std::size_t i = 0; i < edges_.size(); i++) {
            if(edges_[lowest][i] && !visited[i]) {
                if(i == end) {
                    current_lowest = i;
                    break;
                }
                if(sums[i] < current) {
                    current = sums[i];
                    current_lowest = i;
                }
            }
        }
        lowest = current_lowest;
        result.push_back(lowest);
        visited[lowest] = true;

    }

    sum = sums[result[result.size() - 1]];

    return std::make_pair(sum, result);
}
