//
//  FEM_problem.hpp
//


#ifndef FEM_problem_hpp
#define FEM_problem_hpp

#include "PML_solver.hpp"

/* -------------------------
 
 CLASS FEM_problem
 
 --------------------------- */

namespace FEM
{
    
    class FEM_problem {
        
    private:
        int problem_type ;
        long Nsamples;
        
        FEM::PML_solver problem;
        
        int  solve_determinist();
        int  solve_random();  // TO DO
        
    public:
        int  solve();
        int  read_cfg(std::string CONFIG_FILE);
        
    };
}

#endif
