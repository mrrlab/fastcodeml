#include <iostream>
#include "Newick.h"

using namespace std;

int main(int argc, char **argv) {
  if (argc < 2) {
    cerr << "Please specify file name\n";
    return 1;
  }
  Newick tree(0);
  tree.readFile(argv[1]);
  tree.printTreeAnnotated(cout);
}
