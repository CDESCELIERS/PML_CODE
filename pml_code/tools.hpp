//
//  tools.hpp
//

#ifndef tools_hpp
#define tools_hpp

#include <string>
#include "FEM_type.hpp"


/* -------------------------
 
 TOOLS
 
 --------------------------- */

void start_chrono();

void chrono_message(std::string message);

void printProgress (double percentage);

void save_sp2ascii(FEM::SP_mat_complex &MAT, std::string filename);

void save_sp2ascii(FEM::SP_mat &MAT, std::string filename);

void save_vec2ascii(Eigen::VectorXcd &VEC, std::string filename);

void save_vec2ascii(Eigen::VectorXd  &VEC, std::string filename);

void save_mat2ascii(Eigen::Matrix<unsigned,Eigen::Dynamic, 7> &MAT, std::string filename);

void save_mat2ascii(Eigen::Matrix<double,2, Eigen::Dynamic> &MAT, std::string filename);

void save_int2ascii(int value, std::string filename);

#endif
