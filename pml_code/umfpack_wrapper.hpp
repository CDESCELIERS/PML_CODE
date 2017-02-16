//
//  umfpack_wrapper.hpp
//  pml_code
//
//  Created by Christophe Desceliers on 13/02/2017.
//  Copyright © 2017 Christophe Desceliers. All rights reserved.
//

#ifndef umfpack_wrapper_hpp
#define umfpack_wrapper_hpp

#include <armadillo>


namespace FEM
{
    class umfpack_wrapper {
        
    private:
        double * Values_x  ;
        double * Values_y  ;
        long * OuterStart  ;
        long * InnerIndices;
        void *Numeric;
        long n_cols;
        long n_rows;
        long n_nz ;
        
    public:
        ~umfpack_wrapper();
        
        umfpack_wrapper(arma::sp_cx_mat &MAT);
        
        umfpack_wrapper(arma::sp_mat &MAT);
        
        arma::cx_vec solve(arma::cx_vec &F);
        
        arma::vec solve(arma::vec &F);
        
    };

}


#endif
