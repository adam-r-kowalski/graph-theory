/* file: main.cpp
 * description: entry point of the program
 */

#include "header/cli.h"

int main(int argc, char *argv[]) {
  /*
  vector<string> args;

  for (int i = 1; i < argc; i++)
    args.push_back(argv[i]);

  if (args.size() == 0)
    return 0;

  cli::parse(args);
  */

  matrix<double> a =
    {{1, 2, 3}
    ,{3, 2, 1}
    ,{1, 2, 3}};

  HamiltonianCircuit(a);

  return 0;
}
