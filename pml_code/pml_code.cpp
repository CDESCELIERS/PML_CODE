//
//  pml_code.cpp
//

#include "FEM_problem.hpp"

using namespace FEM;

/* -------------------------
 
 MAIN FUNCTION
 
 --------------------------- */

int main(int argc, char** argv)
{

    FEM_problem MY_PROBLEM;
    
    string CONFIG_FILE(argv[1]);
    
    if (MY_PROBLEM.read_cfg(CONFIG_FILE) < 0) return -1;
    
    MY_PROBLEM.solve();
    
    return 0;
}



