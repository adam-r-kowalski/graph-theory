/* file: print.h
 *  description: provides utitlies for viewing a variety of data structures
 */

#pragma once

#include <iostream>
#include <utility>
#include <iomanip>

#include "linear.h"

#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */

template <class V>
void print(const V &value);

template <class V>
void print(const vector<V> &vec);

template <class A, class B>
void print(const pair<A,B> &p);

template <class M>
void print(const matrix<M> &mat);

void print(const vector<gsl_complex> &vec);

void print(const gsl_complex &complex);

void print(const map<string, vector<double>> &vecs);

void print(const map<string, matrix<double>> &mats);

template <class V>
void println(const V &value);

void newline(unsigned long times = 1);

template <class V>
void print(const V &value) {
  cout << value;
}

template <class V>
void print(const vector<V> &vec) {
  unsigned long size = vec.size();
  print("|");
  for (unsigned long i = 0; i < size; i++) {
    print(vec[i]);
    if (i != size - 1) {
      print(", ");
    }
  }
  print("|");
}

template <class A, class B>
void print(const pair<A,B> &p) {
  print("(");
  print(p.first);
  print(", ");
  print(p.second);
  print(")");
};

template <class M>
void print(const matrix<M> &mat) {
  unsigned long rows = mat.size();
  unsigned long columns = mat[0].size();

  for (unsigned long r = 0; r < rows; r++) {
    if (r != 0) { newline(); }
    print("|");
    for (unsigned long c = 0; c < columns; c++) {
      print(mat[r][c]);
      if (c != columns - 1) {
        print(", ");
      }
    }
    print("|");
  }
}

void print(const vector<gsl_complex> &vec) {
  cout << "EigenValues" << endl;
  for (unsigned long i = 0; i < vec.size(); i++) {
    cout << "\u03BB" << i + 1
    << " = " << GSL_REAL(vec[i])
    << " + " << GSL_IMAG(vec[i]) << "i" << endl;
  }
}

void print(const gsl_complex &complex) {
  cout << GSL_REAL(complex) << " + " << GSL_IMAG(complex) << "i" << endl;
}

void print(const map<string, vector<double>> &vecs) {
  for (auto &v : vecs) {
    cout << v.first << " = ";
    println(v.second);
  }
}

void print(const map<string, matrix<double>> &mats) {
  for (auto &m : mats) {
    cout << m.first << endl;
    println(m.second);
    newline();
  }
}

template <class V>
void println(const V &value) {
  print(value);
  newline();
}

void newline(unsigned long times) {
  for (unsigned long i = 0; i < times; i++) {
    cout << endl;
  }
}

