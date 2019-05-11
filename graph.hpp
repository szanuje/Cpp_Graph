#pragma once

#include <cstdint>
#include <vector>
#include <optional>
#include <limits>
#include <stack>
#include <iostream>
#include <algorithm>
#include <cstring>

using namespace std;

template <typename V, typename E>
class Graph
{
public:
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
        unsigned int idx_;

    public:
        bool operator==(const VerticesIterator &vi) const;
        bool operator!=(const VerticesIterator &vi) const;
        VerticesIterator& operator++();
        VerticesIterator operator++(int);
        V& operator*() const;
        V* operator->() const;
    };

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
        unsigned int Vidx_;
        unsigned int Hidx_;

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

    vector<unsigned int> BFS(std::size_t index) const;
    vector<unsigned int> DFS(std::size_t index) const;
    unsigned findIndexOfVertice(const V &vertice);
    V getVertice(unsigned int idx) const { return vertices_[idx]; }
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
            if(!i) cout << "0 ";
            else cout << i.value() << " ";
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
    for(unsigned int i = 0; i < edges_.size(); i++) {
        for(unsigned int j = 0; j < edges_.size(); j++) {
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
unsigned int Graph<V,E>::findIndexOfVertice(const V &vertice) {
    for(unsigned int i = 0; i < vertices_.size(); i++) {
        if(vertices_[i] == vertice) return i;
    }
    return 0;
}

// TRAVERSALS

template <typename V, typename E>
vector<unsigned int> Graph<V,E>::BFS(std::size_t index) const {
    bool visited[vertices_.size()];
    memset(visited, false, sizeof(visited));
    vector<unsigned int> v;
    v.push_back(index);
    visited[index] = true;

    for(unsigned int i = 0; i < v.size(); i++){

        //cout << vertices_[ v[i] ] << ", ";

        for(unsigned int j = 0; j < edges_.size(); j++) {

            if(edges_[v[i]][j] && !visited[j]) {
                v.push_back(j);
                visited[j] = true;

            }
        }
    }

    return v;
}

template <typename V, typename E>
vector<unsigned int> Graph<V,E>::DFS(std::size_t index) const {
    bool visited[vertices_.size()];
    memset(visited, false, sizeof(visited));
    vector<unsigned int> v;
    stack<unsigned int> s;
    v.push_back(index);
    visited[index] = true;

    while(v.size() != vertices_.size()) {

        for(unsigned int i = v.size() - 1; i < v.size(); i++) {

            //cout << vertices_[ v[i] ] << ", ";
            bool go = false;

            for(unsigned int j = 0; j < edges_.size(); j++) {

                if(edges_[v[i]][j] && !visited[j]) {

                    if(!go) {
                        v.push_back(j);
                        visited[j] = true;
                        go = true;
                    }

                    if(go) {
                        s.push(j);
                    }
                }
            }
        }

        while(true) {
            if(s.empty()) break;
            if(visited[s.top()]) s.pop();
            if(!s.empty()) {
                if(!visited[s.top()]) break;
            }
        }

        if(!s.empty()) {
            v.push_back(s.top());
            visited[s.top()] = true;
            s.pop();
        }
    }

    return v;
}
