// dydaktykafais@outlook.com

#include <iostream>
#include <functional>
#include "graph.hpp"

using namespace std;

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

    DFS<std::string, double>(g, 3u, [](const std::string &v) -> decltype(auto) {
        std::cout << v << ", ";
    });


    cout << endl;

    //    Graph<std::string, int> g;
    //    g.insertVertex("V1");
    //    g.insertVertex("V2");
    //    g.insertVertex("V3");
    //    g.insertVertex("V4");

    //    g.insertEdge(0, 0, 1);
    //    g.insertEdge(0, 1, 2);
    //    g.insertEdge(1, 2, 3);
    //    g.insertEdge(2, 2, 4);
    //    g.insertEdge(3, 2, 5);
    //    g.insertEdge(3, 0, 6);
    //    g.insertEdge(0, 3, 7);
    //    g.insertEdge(1, 3, 8);

    //    std::cout << (g.removeEdge(0, 1) ? "Udalo sie" : "Nie udalo sie") << std::endl;
    //    std::cout << (g.removeEdge(1, 0) ? "Udalo sie" : "Nie udalo sie") << std::endl;
    //    std::cout << (g.removeVertex(2) ? "Udalo sie" : "Nie udalo sie") << std::endl;
    //    std::cout << (g.removeVertex(5) ? "Udalo sie" : "Nie udalo sie") << std::endl;
    //    std::cout << "Nr of vertices: " << g.nrOfVertices() << std::endl;
    //    std::cout << "Nr of edges: " << g.nrOfEdges() << std::endl;

    //    g.printNeighborhoodMatrix();

    //    cout << "Vertices: ";
    //    for(const auto &i : g) cout << i << ", ";
    //    cout << endl;

    //    cout << "Edges: ";
    //    for(auto i = g.beginEdges(); i != g.endEdges(); ++i) {
    //        cout << *i << ", ";
    //    }
    //    cout << endl;

    //    g.insertVertex("a1");
    //    g.insertVertex("a2");
    //    g.insertVertex("a3");
    //    g.insertVertex("a4");
    //    g.insertEdge(4, 0, 10);
    //    g.insertEdge(4, 1, 20);
    //    g.insertEdge(1, 5, 30);
    //    g.insertEdge(2, 7, 40);
    //    g.insertEdge(3, 7, 50);
    //    g.insertEdge(3, 4, 60);
    //    g.insertEdge(6, 3, 70);
    //    g.insertEdge(1, 6, 80);

    //    g.removeEdge(4,0);

    //    cout << endl;
    //    std::cout << "Nr of vertices: " << g.nrOfVertices() << std::endl;
    //    std::cout << "Nr of edges: " << g.nrOfEdges() << std::endl;
    //    g.printNeighborhoodMatrix();

    //    cout << "Vertices: ";
    //    for(const auto &i : g) cout << i << ", ";
    //    cout << endl;

    //    cout << "Edges: ";
    //    for(auto i = g.beginEdges(); i != g.endEdges(); ++i) {
    //        cout << *i << ", ";
    //    }
    //    cout << endl;

    return 0;
}
