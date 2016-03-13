#pragma once

#include "linear.h"
#include <sstream>
#include "print.h"

namespace cli {
  const string argumentError = "error: wrong number of arguments supplied";
  const string missingFunctionError = "error: no function by that name";
  const string argumentTypeError = "error: wrong type of argument supplied";

  void parse(const vector<string> &args);
  double getScalar(const string &input);
  vector<double> getVector(const string &input);
  matrix<double> getMatrix(const string &input);

  void parse(const vector<string> &args) {
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
    else if (func == "scale") {
      if (args.size() != 3)
        cerr << argumentError << endl;
      println(scale(getScalar(args[1]), getVector(args[2])));
    }
    else if (func == "lower") {
      if (args.size() != 2)
        cerr << argumentError << endl;
      println(lower(getMatrix(args[1])));
    }
    else {
      cerr << missingFunctionError << endl;
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
      return getVector(input);
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
    matrix<double> ret;
    vector<double> elements;
    string line;
    string word;
    unsigned long rowLength = 0;

    stringstream linestream(input);
    do {
      elements = {};

      getline(linestream, line);
      istringstream iss(line);
      while (iss >> word) {
        double number;
        stringstream ss(word);
        ss >> number;
        if (ss.fail()) {
          cout << argumentError << endl;
          return {};
        }

        elements.push_back(number);
      }

      if (elements.size() != 0) {
        if (rowLength == 0)
          rowLength = elements.size();
        else if (rowLength != elements.size()) {
          cout << argumentError << endl;
          return {};
        }
        ret.push_back(elements);
      }
    } while (elements.size() != 0);

    return ret;
  }
}

