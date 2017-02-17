//
//  tools.hpp
//

#ifndef tools_hpp
#define tools_hpp

#include <string>

/* -------------------------
 
 TOOLS
 
 --------------------------- */

void start_chrono();

void chrono_message(std::string message);

void printProgress (double percentage);

void save_sp2ascii(SP_mat_complex &MAT, string filename);

void save_sp2ascii(SP_mat &MAT, string filename);

void save_vec2ascii(VectorXcd &VEC, string filename);

void save_vec2ascii(VectorXd  &VEC, string filename);

void save_mat2ascii(Matrix<unsigned,Dynamic, 7> &MAT, string filename);

void save_mat2ascii(Matrix<double,2, Dynamic> &MAT, string filename);

void save_int2ascii(int value, string filename);

#endif
