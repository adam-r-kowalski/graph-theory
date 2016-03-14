/* file: main.cpp
 * description: entry point of the program
 */

#include "header/cli.h"

void sandbox();

int main(int argc, char *argv[]) {
  cli::parse(argc, argv);

  sandbox();

  return 0;
}

void sandbox() {
  // This function shows examples of using the api

  // make a matrix
  matrix<double> adj =
    {{0, 0, 0, 3, 2, 0, 0}
    ,{0, 0, 3, 0, 0, 0, 0}
    ,{0, 3, 0, 6, 0, 0, 3}
    ,{3, 0, 6, 0, 5, 4, 0}
    ,{2, 0, 0, 5, 0, 2, 0}
    ,{0, 0, 0, 4, 2, 0, 2}
    ,{0, 0, 3, 0, 0, 2, 0}};
  cout << "adjacency matrix" << endl;
  println(adj);
  newline();

  // get the degree matrix
  matrix<double> deg = degree(adj);
  cout << "degree matrix" << endl;
  println(deg);
  newline();

  // get the incidence matrix
  matrix<double> inc = incidence(adj);
  cout << "incidence matrix" << endl;
  println(inc);
  newline();

  // get the laplacian matrix
  matrix<double> lap = laplacian(adj);
  cout << "laplacian matrix" << endl;
  println(lap);
  newline();

  // get the algebraic connectivity
  double ac = algebraicConnectivity(adj);
  cout << "algebraic connectivity = " << ac;
  newline(2);

  // get the eigenvalues
  vector<gsl_complex> evs = eigenvalues(lap);
  println(evs);
  newline();

  // get the euler circuit
  vector<int> ec = eulerCircuit(adj);
  cout << "euler circuit = ";
  println(ec);
  newline();

  // get the minimal spanning tree
  vector<pair<int,int>> mst = minimalSpanningTree(adj);
  double cost = costOfRoute(adj, mst);
  matrix<double> mstAdjacency = adjacency(adj, mst);

  cout << "minimal spanning tree adjacency matrix" << endl;
  println(mstAdjacency);
  newline();

  cout << "total cost of minimal spanning tree = " << cost;
  newline(2);

  // get the hamiltonian circuit
  hamiltonianCircuit(adj);
}

