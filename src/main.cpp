#include <iostream>
#include "Interface.hpp"

int main(int argc, char **argv) {
   try {
       Interface I(argc, argv);
   } catch(ProjectError &e) {
        std::cerr << e.what() << std::endl;
        return e.value();
    }
   return 0;
}
