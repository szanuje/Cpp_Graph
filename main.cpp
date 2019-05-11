#include <iostream>
#include <functional>
#include "graph.hpp"

using namespace std;

template <typename V, typename E>
void BFS(const Graph<V,E> &graph, std::size_t start_idx, std::function<void(const V&)> func) {

    cout << "BFS: ";
    vector<unsigned int> v = graph.BFS(start_idx);
    for(unsigned int i = 0; i < v.size(); i++) {
        func(graph.getVertice(v[i]));
    }
}

template <typename V, typename E>
void DFS(const Graph<V,E> &graph, std::size_t start_idx, std::function<void(const V&)> func) {

    cout << "DFS: ";
    vector<unsigned int> v = graph.DFS(start_idx);
    for(unsigned int i = 0; i < v.size(); i++) {
        func(graph.getVertice(v[i]));
    }
}

int main()
{
    Graph<std::string, double> g;

    for(std::size_t i = 0u; i < 6u; ++i)
    {
        g.insertVertex("data " + std::to_string(i));
    }

    for(std::size_t i = 0u; i < g.nrOfVertices(); ++i)
    {
        for(std::size_t j = 0u; j < g.nrOfVertices(); ++j)
        {
            if((i + j) & 1u || i & 1u)
            {
                g.insertEdge(i, j, ((i != j) ? (i+j)/2. : 1.));
            }
        }
    }
    g.insertEdge(0, 2, 4.);

    std::cout << (g.removeVertex(1) ? "Udalo sie" : "Nie udalo sie") << std::endl;
    std::cout << (g.removeEdge(2, 2) ? "Udalo sie" : "Nie udalo sie") << std::endl;
    std::cout << (g.removeEdge(2, 3) ? "Udalo sie" : "Nie udalo sie") << std::endl;
    std::cout << (g.removeEdge(4, 3) ? "Udalo sie" : "Nie udalo sie") << std::endl;
    std::cout << "Nr of vertices: " << g.nrOfVertices() << std::endl;
    std::cout << "Nr of edges: " << g.nrOfEdges() << std::endl;
    std::cout << std::endl;

    g.printNeighborhoodMatrix();
    std::cout << std::endl;
    std::cout << "Vertices data:" << std::endl;
    for(auto v_it = g.beginVertices(); v_it != g.endVertices(); ++v_it)
    {
        std::cout << *v_it << ", ";
    }
    std::cout << std::endl << std::endl;
    std::cout << "Edges data:" << std::endl;
    for(auto e_it = g.beginEdges(); e_it != g.endEdges(); ++e_it)
    {
        std::cout << *e_it << ", ";
    }

    cout << endl << endl;

    BFS<std::string, double>(g, 3u, [](const std::string &v) -> decltype(auto) {
        std::cout << v << ", ";
    });

    cout << endl;

    DFS<std::string, double>(g, 3u, [](const std::string &v) -> decltype(auto) {
        std::cout << v << ", ";
    });

    cout << endl << endl;

    return 0;
}
