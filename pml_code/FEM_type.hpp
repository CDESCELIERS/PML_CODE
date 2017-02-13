//
//  FEM_type.hpp
//  pml_code
//
//  Created by Christophe Desceliers on 13/02/2017.
//  Copyright © 2017 Christophe Desceliers. All rights reserved.
//

#include <eigen3/Eigen/Sparse>

#ifndef FEM_type_hpp
#define FM_type_hpp


namespace FEM
{
    
    typedef Eigen::Triplet<double> T;
    
    typedef std::complex<double> dcomplex;
    
    typedef std::map<unsigned, Eigen::Vector4d> Materials;
    
    typedef Eigen::SparseMatrix<double,Eigen::ColMajor> SP_mat;
    
    typedef Eigen::SparseMatrix<dcomplex,Eigen::ColMajor> SP_mat_complex;
}


#endif
