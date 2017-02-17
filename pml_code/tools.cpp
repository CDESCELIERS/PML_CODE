//
//  tools.cpp
//

#include "tools.hpp"
#include <chrono>
#include <iostream>

using namespace std;

/* -------------------------
 
 TOOLS
 
 --------------------------- */

chrono::high_resolution_clock::time_point t_start;

void start_chrono()
{
    t_start = chrono::high_resolution_clock::now();
    
}

void chrono_message(string message)
{
    chrono::high_resolution_clock::time_point t_now = chrono::high_resolution_clock::now();
    
    cout << endl <<  "(" << std::chrono::duration<double, std::milli>(t_now-t_start).count()/1000 << " secs) "<< message;
}

void printProgress (double percentage)
{
    string PBSTR("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||");
    int PBWIDTH= 60;
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR.c_str(), rpad, "");
    fflush (stdout);
}




void FEM::save_sp2ascii(SP_mat &MAT, string filename)
{
    ofstream fid;
    fid.open (filename,ios::out);
    cout << filename << " " << MAT.nonZeros() <<  " non zeros" << endl;
    for (int k=0; k<MAT.outerSize(); ++k)
        for (SP_mat::InnerIterator it(MAT,k); it; ++it)
            fid<< it.value() << " " <<  it.row()+1 << " " << it.col()+1<< endl;
    fid.close();
}

void FEM::save_sp2ascii(SP_mat_complex &MAT, string filename)
{
    ofstream fid;
    fid.open (filename,ios::out);
    cout << filename << " " << MAT.nonZeros() <<  " non zeros" << endl;
    for (int k=0; k<MAT.outerSize(); ++k)
        for (SP_mat_complex::InnerIterator it(MAT,k); it; ++it)
            fid<< it.value().real() << " " <<  it.value().imag() << " " << it.row()+1 << " " << it.col()+1<< endl;
    fid.close();
}

void FEM::save_vec2ascii(VectorXcd &VEC, string filename)
{
    ofstream fid;
    fid.open (filename,ios::out);
    cout << filename << " " << VEC.nonZeros() <<  " non zeros" << endl;
    for (int k=0; k<VEC.size(); ++k)
        fid<< VEC(k).real() << " " <<  VEC(k).imag() << endl ;
    fid.close();
}

void FEM::save_vec2ascii(VectorXd &VEC, string filename)
{
    ofstream fid;
    fid.open (filename,ios::out);
    cout << filename << " " << VEC.nonZeros() <<  " non zeros" << endl;
    for (int k=0; k<VEC.size(); ++k)
        fid<< VEC(k) << endl ;
    fid.close();
}

void FEM::save_mat2ascii(Matrix<unsigned,Dynamic,7> &MAT, string filename)
{
    ofstream fid;
    fid.open (filename,ios::out);
    cout << filename << " " << MAT.nonZeros() <<  " non zeros" << endl;
    for (int k=0; k<MAT.rows(); ++k)
    {
        for(int j=0; j<MAT.cols();++j)
            fid<< MAT(k,j) << " " ;
        fid << endl;
    }
    fid.close();
}
void FEM::save_mat2ascii(Matrix<double,2, Dynamic> &MAT, string filename)
{
    ofstream fid;
    fid.open (filename,ios::out);
    cout << filename << " " << MAT.nonZeros() <<  " non zeros" << endl;
    for (int k=0; k<MAT.rows(); ++k)
    {
        for(int j=0; j<MAT.cols();++j)
            fid<< MAT(k,j) << " " ;
        fid << endl;
    }
    fid.close();
}

void FEM::save_int2ascii(int value, string filename)
{
    ofstream fid;
    fid.open (filename,ios::out);
    cout << filename << " 1 non zero" << endl;
    fid<< value <<  endl;    fid.close();
}


