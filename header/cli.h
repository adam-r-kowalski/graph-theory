#pragma once

#include "linear.h"
#include <sstream>
#include "print.h"

namespace cli {
  const string argumentError = "error: wrong number of arguments supplied";
  const string missingFunctionError = "error: no function by that name";
  const string argumentTypeError = "error: wrong type of argument supplied";

  void parse(int argc, char *argv[]);
  double getScalar(const string &input);
  vector<double> getVector(const string &input);
  matrix<double> getMatrix(const string &input);
  vector<string> splitLines(const string &input);

  void parse(int argc, char *argv[]) {
    vector<string> args;
    for (int i = 1; i < argc; i++)
      args.push_back(argv[i]);

    if (args.size() == 0)
      return;

    try {
      string func = args[0];
      if (func == "magnitude") {
        if (args.size() != 2)
          cerr << argumentError << endl;
        println(magnitude(getVector(args[1])));
      }
      else if (func == "dot") {
        if (args.size() != 3)
          cerr << argumentError << endl;
        println(dot(getVector(args[1]), getVector(args[2])));
      }
      else if (func == "scaleVector") {
        if (args.size() != 3)
          cerr << argumentError << endl;
        println(scale(getScalar(args[1]), getVector(args[2])));
      }
      else if (func == "lower") {
        if (args.size() != 2)
          cerr << argumentError << endl;

        println(lower(getMatrix(args[1])));
      }
      else if (func == "upper") {
        if (args.size() != 2)
          cerr << argumentError << endl;
        println(upper(getMatrix(args[1])));
      }
      else if (func == "column") {
        if (args.size() != 3)
          cerr << argumentError << endl;
        println(column(getScalar(args[1]), getMatrix(args[2])));
      }
      else if (func == "vectorMatrixInnerProduct") {
        if (args.size() != 3)
          cerr << argumentError << endl;
        println(innerProduct<double>(getVector(args[1]), getMatrix(args[2])));
      }
      else if (func == "matrixVectorInnerProduct") {
        if (args.size() != 3)
          cerr << argumentError << endl;
        println(innerProduct<double>(getMatrix(args[1]), getVector(args[2])));
      }
      else if (func == "innerProduct") {
        if (args.size() != 3)
          cerr << argumentError << endl;
        println(innerProduct<double>(getMatrix(args[1]), getMatrix(args[2])));
      }
      else if (func == "outerProduct") {
        if (args.size() != 3)
          cerr << argumentError << endl;
        println(outerProduct<double>(getMatrix(args[1]), getMatrix(args[2])));
      }
      else if (func == "diagonal") {
        if (args.size() != 2)
          cerr << argumentError << endl;
        println(diagonal(getMatrix(args[1])));
      }
      else if (func == "minorDiagonal") {
        if (args.size() != 2)
          cerr << argumentError << endl;
        println(minorDiagonal(getMatrix(args[1])));
      }
      else if (func == "trace") {
        if (args.size() != 2)
          cerr << argumentError << endl;
        println(trace(getMatrix(args[1])));
      }
      else if (func == "scaleMatrix") {
        if (args.size() != 3)
          cerr << argumentError << endl;
        println(scale(getScalar(args[1]), getMatrix(args[2])));
      }
      else if (func == "subtract") {
        if (args.size() != 3)
          cerr << argumentError << endl;
        println(subtract(getMatrix(args[1]), getMatrix(args[2])));
      }
      else if (func == "transpose") {
        if (args.size() != 2)
          cerr << argumentError << endl;
        println(transpose<double>(getMatrix(args[1])));
      }
      else if (func == "adjacencyFromLaplacian") {
        if (args.size() != 2)
          cerr << argumentError << endl;
        println(adjacency(getMatrix(args[1])));
      }
      else if (func == "complement") {
        if (args.size() != 2)
          cerr << argumentError << endl;
        println(complement(getMatrix(args[1])));
      }
      else if (func == "degree") {
        if (args.size() != 2)
          cerr << argumentError << endl;
        println(degree(getMatrix(args[1])));
      }
      else if (func == "laplacian") {
        if (args.size() != 2)
          cerr << argumentError << endl;
        println(laplacian(getMatrix(args[1])));
      }
      else if (func == "incidence") {
        if (args.size() != 2)
          cerr << argumentError << endl;
        println(incidence(getMatrix(args[1])));
      }
      else if (func == "algebraicConnectivity") {
        if (args.size() != 2)
          cerr << argumentError << endl;
        println(algebraicConnectivity(getMatrix(args[1])));
      }
      else if (func == "eigenvalues") {
        if (args.size() != 2)
          cerr << argumentError << endl;
        println(eigenvalues(getMatrix(args[1])));
      }
      else if (func == "determinant") {
        if (args.size() != 2)
          cerr << argumentError << endl;
        println(determinant(getMatrix(args[1])));
      }
      else if (func == "connected") {
        if (args.size() != 2)
          cerr << argumentError << endl;
        println(connected(getMatrix(args[1])));
      }
      else if (func == "hasEulerCircuit") {
        if (args.size() != 2)
          cerr << argumentError << endl;
        println(hasEulerCircuit(getMatrix(args[1])));
      }
      else if (func == "eulerCircuit") {
        if (args.size() != 2)
          cerr << argumentError << endl;
        println(eulerCircuit(getMatrix(args[1])));
      }
      else if (func == "minimalSpanningTree") {
        if (args.size() != 2)
          cerr << argumentError << endl;
        auto mat = getMatrix(args[1]);
        println(adjacency(mat, minimalSpanningTree(mat)));
      }
      else if (func == "spectralGap") {
        if (args.size() != 2)
          cerr << argumentError << endl;
        println(spectralGap(getMatrix(args[1])));
      }
      else if (func == "spectralRadius") {
        if (args.size() != 2)
          cerr << argumentError << endl;
        println(spectralRadius(getMatrix(args[1])));
      }
      else if (func == "connectedComponents") {
        if (args.size() != 2)
          cerr << argumentError << endl;
        println(connectedComponents(getMatrix(args[1])));
      }
      else if (func == "numSpanningTrees") {
        if (args.size() != 2)
          cerr << argumentError << endl;
        println(numSpanningTrees(getMatrix(args[1])));
      }
      else if (func == "isSymmetric") {
        if (args.size() != 2)
          cerr << argumentError << endl;
        println(isSymmetric(getMatrix(args[1])));
      }
      else if (func == "isSquare") {
        if (args.size() != 2)
          cerr << argumentError << endl;
        println(isSquare(getMatrix(args[1])));
      }
      else if (func == "isJagged") {
        if (args.size() != 2)
          cerr << argumentError << endl;
        println(isJagged(getMatrix(args[1])));
      }
      else if (func == "hamiltonianCircuit") {
        if (args.size() != 2)
          cerr << argumentError << endl;
        println(hamiltonianCircuit(getMatrix(args[1])));
      }
      else {
        cerr << missingFunctionError << endl;
      }
    } catch (const exception &e) {
       cerr << "error: " << e.what() << endl;
    }

    return;
  }

  vector<double> getVector(const string &input) {
    vector<double> elements;
    string word;

    istringstream iss(input);
    while (iss >> word) {
      double number;
      stringstream ss(word);
      ss >> number;
      if (ss.fail()) {
        cerr << argumentTypeError << endl;
        return {};
      }

      elements.push_back(number);
    }

    if (elements.size() == 0) {
      cout << endl << argumentError << endl;
      return {};
    }

    return elements;
  }

  double getScalar(const string &input) {
    double number;
    string word;
    istringstream iss(input);
    while (iss >> word) {
      stringstream ss(word);
      ss >> number;
      if (ss.fail()) {
        cerr << argumentTypeError << endl;
        return 0;
      }
      return number;
    }
    return 0;
  }

  matrix<double> getMatrix(const string &input) {
    auto lines = splitLines(input);
    matrix<double> ret;
    vector<double> elements;
    string word;
    double number;

    for (auto &l : lines) {
      elements = {};

      istringstream iss(l);
      while (iss >> word) {
        stringstream ss(word);
        ss >> number;
        if (ss.fail()) {
          cerr << argumentTypeError << endl;
          return {};
        }

        elements.push_back(number);
      }

      if (elements.size() == 0) {
        cout << endl << argumentError << endl;
        return {};
      }

      ret.push_back(elements);
    }

    return ret;
  }

  vector<string> splitLines(const string &input) {
    string temp = input;
    vector<string> ret;

    while (true) {
      size_t pos = temp.find('\\');
      if (pos == string::npos) {
        ret.push_back(temp);
        return ret;
      } else {
        string first = temp.substr(0, pos);
        string second = temp.substr(pos + 2, temp.size() - 1);
        ret.push_back(first);
        temp = second;
      }
    }

    return {};
  }
}
