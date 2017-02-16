//
//  FEM.cpp
//

#include "PML_solver.hpp"
#include "tools.hpp"
#include "Tokenizer.hpp"
#include <iostream>
#include <fstream>
#include <algorithm>
#include "umfpack_wrapper.hpp"

#include <armadillo>

using namespace std;
using namespace arma;

using namespace FEM;


/* -------------------------
 
 CLASS PML_solver
 
 --------------------------- */

void FEM::save_sp2ascii(sp_mat &MAT, string filename)
{
    ofstream fid;
    fid.open (filename,ios::out);
    cout << filename << " " << MAT.n_nonzero << " entries" << endl;
    for (sp_mat::const_iterator it = MAT.begin(); it != MAT.end() ;it++)
        fid<< *it << " " <<  it.row()+1 << " " << it.col()+1<< endl;
    fid.close();
}

void FEM::save_sp2ascii(sp_cx_mat &MAT, string filename)
{
    ofstream fid;
    fid.open (filename,ios::out);
    cout << filename << " " << MAT.n_nonzero << " entries" << endl;
    for (sp_cx_mat::const_iterator it = MAT.begin(); it != MAT.end() ;it++)
        fid<< real((*it)) << " " << imag((*it))<< " " <<  it.row()+1 << " " << it.col()+1<< endl;
    fid.close();
}

void FEM::save_vec2ascii(cx_vec &VEC, string filename)
{
    ofstream fid;
    fid.open (filename,ios::out);
    cout << filename << " " << VEC.n_elem <<  " elements" << endl;
    for (int k=0; k<VEC.n_elem; ++k)
        fid<< real(VEC(k)) << " " <<  imag(VEC(k)) << endl ;
    fid.close();
}

void FEM::save_vec2ascii(vec &VEC, string filename)
{
    ofstream fid;
    fid.open (filename,ios::out);
    cout << filename << " " << VEC.n_elem <<  " elements" << endl;
    for (int k=0; k<VEC.n_elem; ++k)
        fid<< VEC(k) << endl ;
    fid.close();
}

void FEM::save_mat2ascii(umat &MAT, string filename)
{
    ofstream fid;
    fid.open (filename,ios::out);
    cout << filename << " " << MAT.n_elem <<  " entries" << endl;
    for (int k=0; k<MAT.n_rows; ++k)
    {
        for(int j=0; j<MAT.n_cols;++j)
            fid<< MAT(k,j) << " " ;
        fid << endl;
    }
    fid.close();
}
void FEM::save_mat2ascii(mat &MAT, string filename)
{
    ofstream fid;
    fid.open (filename,ios::out);
    cout << filename << " " << MAT.n_elem <<  " entries" << endl;
    for (int k=0; k<MAT.n_rows; ++k)
    {
        for(int j=0; j<MAT.n_cols;++j)
            fid<< MAT(k,j) << " " ;
        fid << endl;
    }
    fid.close();
}

void FEM::save_int2ascii(sword value, string filename)
{
    ofstream fid;
    fid.open (filename,ios::out);
    cout << filename << " 1 non zero" << endl;
    fid<< value <<  endl;    fid.close();
}


string PML_solver::get_mshfile()
{
    return settings.msh_file;
}

void PML_solver::set_settings(config_settings new_settings)
{
    settings = new_settings;
}

int PML_solver::read_mesh()
{
    const int element_type = 9;
    const int BC_type = 8;
    
    std::string buffer;
    ifstream input;
    Tokenizer Tokens ;
    int NumberNode=0;
    int NumberElement=0;
    uvec Check_ref;
    
    input.open ( settings.msh_file.c_str ( ) );
    while (!input.eof())
    {
        getline(input, buffer);
        if (Tokenizer(buffer).next() == "$Nodes")
        {
            getline(input, buffer);
            
            NumberNode = stoi(Tokenizer(buffer).next() );
            Table_Noeuds =  mat(2,NumberNode,fill::zeros);
            while (!input.eof())
            {
                getline(input, buffer);
                if (Tokenizer(buffer).next() == "$EndNodes")
                    break;
                else
                {
                    Tokens.set(buffer);
                    int num = stoi(Tokens.next());
                    for (int k=0; k<2; k++)
                        Table_Noeuds(k,num-1) = stod(Tokens.next());
                }
            }
            cout <<  endl << "Number of Nodes = " << NumberNode << endl ;
            Table_Noeuds.resize (2, NumberNode );
            Table_BC = uvec(NumberNode,fill::zeros);
            Check_ref = uvec(NumberNode,fill::zeros);
        }
        else if (Tokenizer(buffer).next() == "$Elements")
        {
            getline(input, buffer);
            NumberElement = stoi(Tokenizer(buffer).next() );
            Table_Elements =  umat(NumberElement,7,fill::zeros);
            NumberElement = 0;
            while (!input.eof())
            {
                getline(input, buffer);
                if (Tokenizer(buffer).next() == "$EndElements")
                    break;
                else
                {
                    Tokens.set(buffer);
                    string token = Tokens.next();
                    
                    int type = stoi(Tokens.next());
                    
                    int nb_tag = stoi(Tokens.next()) ;
                    
                    stoi(Tokens.next());
                    uword tag = (uword) stoi(Tokens.next());
                    
                    for (int k=2; k<nb_tag; k++)
                        stoi(Tokens.next());
                    
                    if (type == element_type)
                    {
                        Table_Elements(NumberElement,0) = tag;
                        for (int k=1; k<7; k++)
                        {
                            int num =  stoi(Tokens.next());
                            Table_Elements(NumberElement,k) = num;
                            Check_ref(num-1) = 1 ;
                        }
                        
                        NumberElement++;
                    }
                    
                    if ( type == BC_type && settings.BC_list.count(tag) != 0 )
                    {
                        for (int k=0; k<3; k++)
                        {
                            int num = stoi(Tokens.next());
                            Table_BC(num-1) = 3;
                        }
                    }
                    
                }
            }
            cout << "Number of Elements = " << NumberElement << endl ;
            Table_Elements.resize ( NumberElement,  7 );
        }
        
    }
    
    for (int k=0;k<NumberNode;k++)
    {
        if (Check_ref(k) != 1 )
        {
            cout << "UNREFERENCED NODE " << k+1 << endl;
            return -1;
        }
    }
    
    int nb_bc = 0;
    for (int k=0;k<NumberNode;k++)
    {
        if (Table_BC(k) != 0 ) nb_bc++;
    }
    cout << "Number of fixed nodes = " << nb_bc << endl <<endl ;
    
    input.close();
    
    set_xmid_ytop();
    sub_element_orientation_check();
    
    return 0;
};

void PML_solver::sub_element_orientation_check()
{
    uword Ne = Table_Elements.n_rows;
    
    for(int e=0; e < Ne; e++)
    {
        vec X1, X2, X3;
        uword num1, num2, num3;
        num1 = Table_Elements(e,1);
        num2 = Table_Elements(e,2);
        num3 = Table_Elements(e,3);
        
        X1 = Table_Noeuds.col(num1-1);
        X2 = Table_Noeuds.col(num2-1);
        X3 = Table_Noeuds.col(num3-1);
        
        double test =(X2(0)-X1(0))*(X3(1)-X1(1)) - (X2(1)-X1(1))*(X3(0)-X1(0)) ;
        if (test < 0)
        {
            uword N;
            N = Table_Elements(e,3);
            Table_Elements(e,3) = Table_Elements(e,2);
            Table_Elements(e,2) = N;
            N = Table_Elements(e,6);
            Table_Elements(e,6) = Table_Elements(e,4);
            Table_Elements(e,4) = N;
        }
    }
}


void PML_solver::set_xmid_ytop()
{
    long Nn = Table_Noeuds.n_cols;
    
    double x_left  = 1000;
    double x_right = -1000;
    ytop = -1000;
    for (uword k=0; k<Nn; k++)
    {
        double y = Table_Noeuds(1,k) ;
        double x = Table_Noeuds(0,k) ;
        if ( y > ytop ) ytop = y ;
        if (x_left  > x) x_left  = x;
        if (x_right < x) x_right = x;
    }
    
    x_mid = (x_left + x_right)/2.f;
}

int PML_solver::sub_element_T6()
{
    uword Ne = Table_Elements.n_rows;
    uword Nn = Table_Noeuds.n_cols;
    
    uword Np = 6;
    
    double a = 0.445948490915965;
    double b = 0.091576213509771;
    
    mat Liste_Yp =   { { b, 1-2*b, b    , a    , a, 1-2*a },
                       { b, b    , 1-2*b, 1-2*a, a, a     } } ;
    
    double P1 = 0.11169079483905;
    double P2 = 0.0549758718227661;
    
    rowvec Liste_gammap = { P2, P2, P2, P1, P1, P1 } ;
    
    mat M_I2 = eye<mat>(2,2);
    mat M_I3 = eye<mat>(3,3);
    M_I3(2,2) = sqrt(2);
    
    double s2 = sqrt(2)/2.f;
    double s3 = 0.5;
    double s4 = 0.25;
    
    mat M_P =
    {   { 1,  0,  0, 0 } ,
        { 0,  0,  0, 1 } ,
        { 0, s2, s2, 0 }
    } ;
    
    mat C_coeff =
    {
        { 1,  0,  0,  0,  0,  0 } ,
        { 1,  1,  0,  1,  0,  0 } ,
        { 1,  0,  1,  0,  1,  0 } ,
        { 1, s3,  0, s4,  0,  0 } ,
        { 1, s3, s3, s4, s4, s4 } ,
        { 1,  0, s3,  0, s4,  0 }
    };
    
    
    vec V(6), V_coeff ;
    
    V.zeros();
    V(0) = 1;
    V_coeff = arma::solve(C_coeff, V);
    
    double a1 = V_coeff(0);
    double b1 = V_coeff(1);
    double c1 = V_coeff(2);
    double d1 = V_coeff(3);
    double e1 = V_coeff(4);
    double f1 = V_coeff(5);
    
    V.zeros();
    V(1) = 1;
    V_coeff = arma::solve(C_coeff, V);
    
    double a2 = V_coeff(0);
    double b2 = V_coeff(1);
    double c2 = V_coeff(2);
    double d2 = V_coeff(3);
    double e2 = V_coeff(4);
    double f2 = V_coeff(5);
    
    V.zeros();
    V(2) = 1;
    V_coeff = arma::solve(C_coeff, V);
    
    double a3 = V_coeff(0);
    double b3 = V_coeff(1);
    double c3 = V_coeff(2);
    double d3 = V_coeff(3);
    double e3 = V_coeff(4);
    double f3 = V_coeff(5);
    
    V.zeros();
    V(3) = 1;
    V_coeff = arma::solve(C_coeff, V);
    
    double a4 = V_coeff(0);
    double b4 = V_coeff(1);
    double c4 = V_coeff(2);
    double d4 = V_coeff(3);
    double e4 = V_coeff(4);
    double f4 = V_coeff(5);
    
    V.zeros();
    V(4) = 1;
    V_coeff = arma::solve(C_coeff, V);
    
    double a5 = V_coeff(0);
    double b5 = V_coeff(1);
    double c5 = V_coeff(2);
    double d5 = V_coeff(3);
    double e5 = V_coeff(4);
    double f5 = V_coeff(5);
    
    V.zeros();
    V(5) = 1;
    V_coeff = arma::solve(C_coeff, V);
    
    double a6 = V_coeff(0);
    double b6 = V_coeff(1);
    double c6 = V_coeff(2);
    double d6 = V_coeff(3);
    double e6 = V_coeff(4);
    double f6 = V_coeff(5);
    
    cube T_K2, T_K, T_M, T_D, T_D0, T_H, T_H1_tilde, T_K1_tilde, T_K_tilde, T_M_tilde, T_D_tilde, T_K1, T_D1;
    ucube T_I_11, T_J_11, T_I_22, T_J_22, T_I_12, T_J_12, T_I_21, T_J_21;
    
    T_I_11.zeros(12,12,Ne);
    T_J_11.zeros(12,12,Ne);
    T_I_22.zeros(18,18,Ne);
    T_J_22.zeros(18,18,Ne);
    
    T_I_12.zeros(12,18,Ne);
    T_J_12.zeros(12,18,Ne);
    T_I_21.zeros(18,12,Ne);
    T_J_21.zeros(18,12,Ne);
    
    T_K2.zeros(12,12,Ne);
    T_K.zeros(12,12,Ne);
    T_M.zeros(12,12,Ne);
    T_D.zeros(12,12,Ne);
    T_D0.zeros(12,12,Ne);
    T_H.zeros(12,12,Ne);
    T_H1_tilde.zeros(12,18,Ne);
    T_K1_tilde.zeros(12,18,Ne);
    T_K_tilde.zeros(18,18,Ne);
    T_M_tilde.zeros(18,18,Ne);
    T_D_tilde.zeros(18,18,Ne);
    T_K1.zeros(18,12,Ne);
    T_D1.zeros(18,12,Ne);
    
    for(int e = 0; e<Ne; e++)
    {
        uword n1 = Table_Elements(e,1);
        uword n2 = Table_Elements(e,2);
        uword n3 = Table_Elements(e,3);
        uword n4 = Table_Elements(e,4);
        uword n5 = Table_Elements(e,5);
        uword n6 = Table_Elements(e,6);
        
        vec Xe1 = Table_Noeuds.col(n1-1);
        vec Xe2 = Table_Noeuds.col(n2-1);
        vec Xe3 = Table_Noeuds.col(n3-1);
        vec Xe4 = Table_Noeuds.col(n4-1);
        vec Xe5 = Table_Noeuds.col(n5-1);
        vec Xe6 = Table_Noeuds.col(n6-1);
        
        vec Xe =
        {
            Xe1(0) ,
            Xe1(1) ,
            Xe2(0) ,
            Xe2(1) ,
            Xe3(0) ,
            Xe3(1) ,
            Xe4(0) ,
            Xe4(1) ,
            Xe5(0) ,
            Xe5(1) ,
            Xe6(0) ,
            Xe6(1)
        } ;
        
        mat C_coeff2 =
        {
            { 1, Xe1(0), Xe1(1), Xe1(0)*Xe1(0), Xe1(1)*Xe1(1), Xe1(0)*Xe1(1) } ,
            { 1, Xe2(0), Xe2(1), Xe2(0)*Xe2(0), Xe2(1)*Xe2(1), Xe2(0)*Xe2(1) } ,
            { 1, Xe3(0), Xe3(1), Xe3(0)*Xe3(0), Xe3(1)*Xe3(1), Xe3(0)*Xe3(1) } ,
            { 1, Xe4(0), Xe4(1), Xe4(0)*Xe4(0), Xe4(1)*Xe4(1), Xe4(0)*Xe4(1) } ,
            { 1, Xe5(0), Xe5(1), Xe5(0)*Xe5(0), Xe5(1)*Xe5(1), Xe5(0)*Xe5(1) } ,
            { 1, Xe6(0), Xe6(1), Xe6(0)*Xe6(0), Xe6(1)*Xe6(1), Xe6(0)*Xe6(1) }
        };
        
        vec V_coeff_tile;
        
        V.zeros(6);
        V(0) = 1;
        V_coeff_tile = arma::solve(C_coeff2, V);
        
        double atilde_1 = V_coeff_tile(0);
        double btilde_1 = V_coeff_tile(1);
        double ctilde_1 = V_coeff_tile(2);
        double dtilde_1 = V_coeff_tile(3);
        double etilde_1 = V_coeff_tile(4);
        double ftilde_1 = V_coeff_tile(5);
        
        
        V.zeros(6);
        V(1) = 1;
        V_coeff_tile = arma::solve(C_coeff2, V);
        
        double atilde_2 = V_coeff_tile(0);
        double btilde_2 = V_coeff_tile(1);
        double ctilde_2 = V_coeff_tile(2);
        double dtilde_2 = V_coeff_tile(3);
        double etilde_2 = V_coeff_tile(4);
        double ftilde_2 = V_coeff_tile(5);
        
        V.zeros(6);
        V(2) = 1;
        V_coeff_tile = arma::solve(C_coeff2, V);
        
        double atilde_3 = V_coeff_tile(0);
        double btilde_3 = V_coeff_tile(1);
        double ctilde_3 = V_coeff_tile(2);
        double dtilde_3 = V_coeff_tile(3);
        double etilde_3 = V_coeff_tile(4);
        double ftilde_3 = V_coeff_tile(5);
        
        V.zeros(6);
        V(3) = 1;
        V_coeff_tile = arma::solve(C_coeff2, V);
        
        double atilde_4 = V_coeff_tile(0);
        double btilde_4 = V_coeff_tile(1);
        double ctilde_4 = V_coeff_tile(2);
        double dtilde_4 = V_coeff_tile(3);
        double etilde_4 = V_coeff_tile(4);
        double ftilde_4 = V_coeff_tile(5);
        
        V.zeros(6);
        V(4) = 1;
        V_coeff_tile = arma::solve(C_coeff2, V);
        
        double atilde_5 = V_coeff_tile(0);
        double btilde_5 = V_coeff_tile(1);
        double ctilde_5 = V_coeff_tile(2);
        double dtilde_5 = V_coeff_tile(3);
        double etilde_5 = V_coeff_tile(4);
        double ftilde_5 = V_coeff_tile(5);
        
        V.zeros(6);
        V(5) = 1;
        V_coeff_tile = arma::solve(C_coeff2, V);
        
        double atilde_6 = V_coeff_tile(0);
        double btilde_6 = V_coeff_tile(1);
        double ctilde_6 = V_coeff_tile(2);
        double dtilde_6 = V_coeff_tile(3);
        double etilde_6 = V_coeff_tile(4);
        double ftilde_6 = V_coeff_tile(5);
        
        mat M_K2e = zeros<mat>(12,12);
        mat M_Ke  = zeros<mat>(12,12);
        mat M_Me  = zeros<mat>(12,12);
        
        mat M_De  = zeros<mat>(12,12);
        mat M_He  = zeros<mat>(12,12);
        mat M_D0e = zeros<mat>(12,12);
        
        mat M_Me_tilde = zeros<mat>(18,18);
        mat M_De_tilde = zeros<mat>(18,18);
        mat M_Ke_tilde = zeros<mat>(18,18);
        
        mat M_K1e_tilde = zeros<mat>(12,18);
        mat M_H1e_tilde = zeros<mat>(12,18);
        
        mat M_K1e = zeros<mat>(18,12);
        mat M_D1e = zeros<mat>(18,12);
        
        for(int p = 0; p<Np; p++)
        {
            
            vec Yp = Liste_Yp.col(p);
            double gammap = Liste_gammap(p);
            
            double E1Yp = a1 + b1*Yp(0) + c1*Yp(1) + d1*Yp(0)*Yp(0) + e1*Yp(1)*Yp(1) + f1*Yp(0)*Yp(1);
            double E2Yp = a2 + b2*Yp(0) + c2*Yp(1) + d2*Yp(0)*Yp(0) + e2*Yp(1)*Yp(1) + f2*Yp(0)*Yp(1);
            double E3Yp = a3 + b3*Yp(0) + c3*Yp(1) + d3*Yp(0)*Yp(0) + e3*Yp(1)*Yp(1) + f3*Yp(0)*Yp(1);
            double E4Yp = a4 + b4*Yp(0) + c4*Yp(1) + d4*Yp(0)*Yp(0) + e4*Yp(1)*Yp(1) + f4*Yp(0)*Yp(1);
            double E5Yp = a5 + b5*Yp(0) + c5*Yp(1) + d5*Yp(0)*Yp(0) + e5*Yp(1)*Yp(1) + f5*Yp(0)*Yp(1);
            double E6Yp = a6 + b6*Yp(0) + c6*Yp(1) + d6*Yp(0)*Yp(0) + e6*Yp(1)*Yp(1) + f6*Yp(0)*Yp(1);
            
            mat M_EYp(2,12);
            
            M_EYp.cols(0 , 1) = E1Yp*M_I2;
            M_EYp.cols(2 , 3) = E2Yp*M_I2;
            M_EYp.cols(4 , 5) = E3Yp*M_I2;
            M_EYp.cols(6 , 7) = E4Yp*M_I2;
            M_EYp.cols(8 , 9) = E5Yp*M_I2;
            M_EYp.cols(10,11) = E6Yp*M_I2 ;
            
            vec Xp = M_EYp * Xe ;
            
            mat M_E_PML_Xp(3,18);
            mat M_EXp(2,12);
            
            {
                double E1Xp = atilde_1 + btilde_1*Xp(0) + ctilde_1*Xp(1)
                + dtilde_1*Xp(0)*Xp(0) + etilde_1*Xp(1)*Xp(1) + ftilde_1*Xp(0)*Xp(1);
                
                double E2Xp = atilde_2 + btilde_2*Xp(0) + ctilde_2*Xp(1)
                + dtilde_2*Xp(0)*Xp(0) + etilde_2*Xp(1)*Xp(1) + ftilde_2*Xp(0)*Xp(1);
                
                double E3Xp = atilde_3 + btilde_3*Xp(0) + ctilde_3*Xp(1)
                + dtilde_3*Xp(0)*Xp(0) + etilde_3*Xp(1)*Xp(1) + ftilde_3*Xp(0)*Xp(1);
                
                double E4Xp = atilde_4 + btilde_4*Xp(0) + ctilde_4*Xp(1)
                + dtilde_4*Xp(0)*Xp(0) + etilde_4*Xp(1)*Xp(1) + ftilde_4*Xp(0)*Xp(1);
                
                double E5Xp = atilde_5 + btilde_5*Xp(0) + ctilde_5*Xp(1)
                + dtilde_5*Xp(0)*Xp(0) + etilde_5*Xp(1)*Xp(1) + ftilde_5*Xp(0)*Xp(1);
                
                double E6Xp = atilde_6 + btilde_6*Xp(0) + ctilde_6*Xp(1)
                + dtilde_6*Xp(0)*Xp(0) + etilde_6*Xp(1)*Xp(1) + ftilde_6*Xp(0)*Xp(1);
                
                M_E_PML_Xp.cols(0,  2) = E1Xp*M_I3;
                M_E_PML_Xp.cols(3 ,  5) = E2Xp*M_I3;
                M_E_PML_Xp.cols(6 ,  8) = E3Xp*M_I3;
                M_E_PML_Xp.cols(9 , 11) = E4Xp*M_I3;
                M_E_PML_Xp.cols(12 , 14) = E5Xp*M_I3;
                M_E_PML_Xp.cols(15 , 17) = E6Xp*M_I3 ;
                
                M_EXp.cols(0 , 1) = E1Xp*M_I2;
                M_EXp.cols(2 , 3) = E2Xp*M_I2;
                M_EXp.cols(4 , 5) = E3Xp*M_I2;
                M_EXp.cols(6 , 7) = E4Xp*M_I2;
                M_EXp.cols(8 , 9) = E5Xp*M_I2;
                M_EXp.cols(10,11) = E6Xp*M_I2 ;
            }
            
            mat M_Ce(4,12);
            {
                double A1 = btilde_1 + 2.f*dtilde_1*Xp(0) + ftilde_1*Xp(1);
                double A2 = btilde_2 + 2.f*dtilde_2*Xp(0) + ftilde_2*Xp(1);
                double A3 = btilde_3 + 2.f*dtilde_3*Xp(0) + ftilde_3*Xp(1);
                double A4 = btilde_4 + 2.f*dtilde_4*Xp(0) + ftilde_4*Xp(1);
                double A5 = btilde_5 + 2.f*dtilde_5*Xp(0) + ftilde_5*Xp(1);
                double A6 = btilde_6 + 2.f*dtilde_6*Xp(0) + ftilde_6*Xp(1);
                
                double B1 = ctilde_1 + 2.f*etilde_1*Xp(1) + ftilde_1*Xp(0);
                double B2 = ctilde_2 + 2.f*etilde_2*Xp(1) + ftilde_2*Xp(0);
                double B3 = ctilde_3 + 2.f*etilde_3*Xp(1) + ftilde_3*Xp(0);
                double B4 = ctilde_4 + 2.f*etilde_4*Xp(1) + ftilde_4*Xp(0);
                double B5 = ctilde_5 + 2.f*etilde_5*Xp(1) + ftilde_5*Xp(0);
                double B6 = ctilde_6 + 2.f*etilde_6*Xp(1) + ftilde_6*Xp(0);
                
                M_Ce.submat(span(0, 1), span(0 , 1)) = A1*M_I2;
                M_Ce.submat(span(0, 1), span(2 , 3)) = A2*M_I2;
                M_Ce.submat(span(0, 1), span(4 , 5)) = A3*M_I2;
                M_Ce.submat(span(0, 1), span(6 , 7)) = A4*M_I2;
                M_Ce.submat(span(0, 1), span(8 , 9)) = A5*M_I2;
                M_Ce.submat(span(0, 1), span(10,11)) = A6*M_I2 ;
                
                M_Ce.submat(span(2, 3), span(0 , 1)) = B1*M_I2;
                M_Ce.submat(span(2, 3), span(2 , 3)) = B2*M_I2;
                M_Ce.submat(span(2, 3), span(4 , 5)) = B3*M_I2;
                M_Ce.submat(span(2, 3), span(6 , 7)) = B4*M_I2;
                M_Ce.submat(span(2, 3), span(8 , 9)) = B5*M_I2;
                M_Ce.submat(span(2, 3), span(10,11)) = B6*M_I2 ;
            }
        
            mat M_Btilde  = M_P * M_Ce ;
            
            double JeYp;
            {
                double A1 = b1 + 2.f*d1*Yp(0) + f1*Yp(1);
                double A2 = b2 + 2.f*d2*Yp(0) + f2*Yp(1);
                double A3 = b3 + 2.f*d3*Yp(0) + f3*Yp(1);
                double A4 = b4 + 2.f*d4*Yp(0) + f4*Yp(1);
                double A5 = b5 + 2.f*d5*Yp(0) + f5*Yp(1);
                double A6 = b6 + 2.f*d6*Yp(0) + f6*Yp(1);
                
                double B1 = c1 + 2.f*e1*Yp(1) + f1*Yp(0);
                double B2 = c2 + 2.f*e2*Yp(1) + f2*Yp(0);
                double B3 = c3 + 2.f*e3*Yp(1) + f3*Yp(0);
                double B4 = c4 + 2.f*e4*Yp(1) + f4*Yp(0);
                double B5 = c5 + 2.f*e5*Yp(1) + f5*Yp(0);
                double B6 = c6 + 2.f*e6*Yp(1) + f6*Yp(0);
                
                vec Phi1 = A1*Xe1 + A2*Xe2 +  A3*Xe3  + A4*Xe4 + A5*Xe5 + A6*Xe6;
                vec Phi2 = B1*Xe1 + B2*Xe2 +  B3*Xe3  + B4*Xe4 + B5*Xe5 + B6*Xe6;
                
                mat Mat_J;
                Mat_J = join_horiz(Phi1, Phi2) ;
                
                JeYp = det(Mat_J);
                if (JeYp < 0)
                {
                    cout << endl << "Jacobian is negative" << endl ;
                    return -1 ;
                }
            }
            
            uword num_mat = Table_Elements(e,0);
            
            double E,v,rho,Xi;
            
            vec parameters = settings.Mat[num_mat];
            E   = parameters(0);
            v   = parameters(1);
            rho = parameters(2);
            Xi  = parameters(3);
            
            double mu1_PML = 0;
            double mu2_PML = 0;
            
            if ( Xp(1) <= settings.ybot_pml)
            {
                mu2_PML = settings.a_PML * pow(settings.ybot_pml - Xp(1),3);
                mu2_PML = min( mu2_PML , settings.max_mu_PML);
            }
            
            if ( abs(Xp(0) - x_mid) >= settings.x_pml)
            {
                mu1_PML = settings.a_PML * pow(abs(Xp(0)-x_mid) - settings.x_pml ,3);
                mu1_PML = min( mu1_PML , settings.max_mu_PML);
            }
            
            mat M_Phat =
            {
                { mu2_PML, 0         , 0         , 0       },
                { 0      , 0         , 0         , mu1_PML },
                { 0      , s2*mu2_PML, s2*mu1_PML, 0       }
            };
            
            mat M_PDEpsilon =
            {
                { mu1_PML, 0          , 0         , 0       },
                { 0      , 0          , 0         , mu2_PML },
                { 0      , s2*mu1_PML , s2*mu2_PML, 0       }
            };
            
            mat  M_Btilde_DEpsilon = M_PDEpsilon * M_Ce ;
            mat  M_Btildehat = M_Phat * M_Ce ;
            
            double Lambda = (E*v)/(1.f+v)/(1.f-2.f*v);
            double Mu = 0.5*E/(1.f+v);
            
            mat M_Ae =
            {
                { Lambda+2*Mu, Lambda      , 0      },
                { Lambda     , Lambda+2*Mu , 0      },
                { 0          , 0           , 4.f*Mu }
            };
            
            M_Ke  = M_Ke + gammap *  M_Btilde.t() * M_Ae * M_Btilde * JeYp ;
            
            M_K2e  = M_K2e + mu1_PML * mu2_PML * gammap * rho * M_EXp.t() * M_EXp * JeYp;
            
            M_Me  = M_Me + gammap * rho * M_EXp.t() * M_EXp * JeYp;
            
            M_De  = M_De +  (mu1_PML + mu2_PML) * gammap * rho * M_EXp.t() * M_EXp * JeYp;
            
            M_He =  M_He +  gammap * M_Btildehat.t() * M_Ae * M_Btilde * JeYp;
            
            M_H1e_tilde = M_H1e_tilde + gammap * M_Btildehat.t() * M_Ae * M_E_PML_Xp * JeYp ;
            
            M_D0e  = M_D0e + 2.f * Xi * sqrt(E * rho) * gammap  * M_EXp.t() * M_EXp * JeYp;
            
            M_Ke_tilde = M_Ke_tilde + gammap * mu1_PML * mu2_PML * M_E_PML_Xp.t() * M_E_PML_Xp * JeYp;
            
            M_De_tilde = M_De_tilde + gammap * (mu1_PML + mu2_PML) * M_E_PML_Xp.t() * M_E_PML_Xp * JeYp;
            
            M_Me_tilde = M_Me_tilde + gammap * M_E_PML_Xp.t() * M_E_PML_Xp * JeYp;
            
            M_K1e =  M_K1e + gammap *  mu1_PML * mu2_PML * M_E_PML_Xp.t() * M_Btilde * JeYp ;
            
            M_K1e_tilde =  M_K1e_tilde +  gammap * M_Btilde.t() * M_Ae * M_E_PML_Xp * JeYp ;
            
            M_D1e = M_D1e + gammap * M_E_PML_Xp.t() * M_Btilde_DEpsilon * JeYp;
            
        }
        
        urowvec T_num_noeuds = Table_Elements.row(e).subvec(1,6);
        
        uvec T_Nek =
        {
            (T_num_noeuds(0) - 1)*2, 1+(T_num_noeuds(0) - 1)*2,
            (T_num_noeuds(1) - 1)*2, 1+(T_num_noeuds(1) - 1)*2,
            (T_num_noeuds(2) - 1)*2, 1+(T_num_noeuds(2) - 1)*2,
            (T_num_noeuds(3) - 1)*2, 1+(T_num_noeuds(3) - 1)*2,
            (T_num_noeuds(4) - 1)*2, 1+(T_num_noeuds(4) - 1)*2,
            (T_num_noeuds(5) - 1)*2, 1+(T_num_noeuds(5) - 1)*2
        };
        
        uvec T_Nek_PML =
        {
            (T_num_noeuds(0) - 1)*3, 1+(T_num_noeuds(0) - 1)*3, 2+(T_num_noeuds(0) - 1)*3,
            (T_num_noeuds(1) - 1)*3, 1+(T_num_noeuds(1) - 1)*3, 2+(T_num_noeuds(1) - 1)*3,
            (T_num_noeuds(2) - 1)*3, 1+(T_num_noeuds(2) - 1)*3, 2+(T_num_noeuds(2) - 1)*3,
            (T_num_noeuds(3) - 1)*3, 1+(T_num_noeuds(3) - 1)*3, 2+(T_num_noeuds(3) - 1)*3,
            (T_num_noeuds(4) - 1)*3, 1+(T_num_noeuds(4) - 1)*3, 2+(T_num_noeuds(4) - 1)*3,
            (T_num_noeuds(5) - 1)*3, 1+(T_num_noeuds(5) - 1)*3, 2+(T_num_noeuds(5) - 1)*3
        };
        
        T_I_11.slice(e) = repmat( T_Nek     ,  1, 12);
        T_J_11.slice(e) = repmat( T_Nek.t() , 12,  1);
        
        T_I_22.slice(e) = repmat( T_Nek_PML     ,  1, 18);
        T_J_22.slice(e) = repmat( T_Nek_PML.t() , 18,  1);
        
        T_I_12.slice(e) = repmat( T_Nek         ,  1, 18);
        T_J_12.slice(e) = repmat( T_Nek_PML.t() , 12,  1);
        
        T_I_21.slice(e) = repmat( T_Nek_PML ,  1, 12);
        T_J_21.slice(e) = repmat( T_Nek.t() , 18,  1);
        
        T_K2.slice(e) = M_K2e;
        T_M .slice(e) = M_Me ;
        T_D .slice(e) = M_De ;
        T_K .slice(e) = M_Ke ;
        T_H .slice(e) = M_He ;
        T_D0.slice(e) = M_D0e;
        
        T_H1_tilde.slice(e) = M_H1e_tilde;
        T_K1_tilde.slice(e) = M_K1e_tilde;
        
        T_M_tilde.slice(e) = M_Me_tilde;
        T_K_tilde.slice(e) = M_Ke_tilde;
        T_D_tilde.slice(e) = M_De_tilde;
        
        T_K1.slice(e) = M_K1e;
        T_D1.slice(e) = M_D1e;
        
    }
    urowvec vec_T_I_11(T_I_11.memptr(),Ne*12*12,false);
    urowvec vec_T_J_11(T_J_11.memptr(),Ne*12*12,false);
    
    urowvec vec_T_I_22(T_I_22.memptr(),Ne*18*18,false);
    urowvec vec_T_J_22(T_J_22.memptr(),Ne*18*18,false);
    
    urowvec vec_T_I_12(T_I_12.memptr(),Ne*12*18,false);
    urowvec vec_T_J_12(T_J_12.memptr(),Ne*12*18,false);
    
    urowvec vec_T_I_21(T_I_21.memptr(),Ne*18*12,false);
    urowvec vec_T_J_21(T_J_21.memptr(),Ne*18*12,false);
    
    umat Location_11 = join_vert(vec_T_I_11, vec_T_J_11) ;
    umat Location_22 = join_vert(vec_T_I_22, vec_T_J_22) ;
    umat Location_12 = join_vert(vec_T_I_12, vec_T_J_12) ;
    umat Location_21 = join_vert(vec_T_I_21, vec_T_J_21) ;
    
    vec vec_K2(T_K2.memptr(),Ne*12*12,false);
    vec vec_M (T_M .memptr(),Ne*12*12,false);
    vec vec_D( T_D .memptr(),Ne*12*12,false);
    vec vec_D0(T_D0.memptr(),Ne*12*12,false);
    vec vec_K (T_K .memptr(),Ne*12*12,false);
    vec vec_H (T_K2.memptr(),Ne*12*12,false);
    
    vec vec_H1_tilde(T_H1_tilde.memptr(),Ne*12*18,false);
    vec vec_K1_tilde(T_K1_tilde.memptr(),Ne*12*18,false);
    
    vec vec_K_tilde(T_K_tilde.memptr(),Ne*18*18,false);
    vec vec_M_tilde(T_M_tilde.memptr(),Ne*18*18,false);
    vec vec_D_tilde(T_D_tilde.memptr(),Ne*18*18,false);
    
    vec vec_K1(T_K1.memptr(),Ne*18*12,false);
    vec vec_D1(T_D1.memptr(),Ne*18*12,false);
    
    M_K2 = sp_mat( true, Location_11, vec_K2, 2*Nn, 2*Nn);
    M_M  = sp_mat( true, Location_11, vec_M , 2*Nn, 2*Nn);
    M_D  = sp_mat( true, Location_11, vec_D , 2*Nn, 2*Nn);
    M_D0 = sp_mat( true, Location_11, vec_D0, 2*Nn, 2*Nn);
    M_K  = sp_mat( true, Location_11, vec_K , 2*Nn, 2*Nn);
    M_H  = sp_mat( true, Location_11, vec_H , 2*Nn, 2*Nn);
    
    M_H1_tilde  = sp_mat( true, Location_12, vec_H1_tilde , 2*Nn, 3*Nn);
    M_K1_tilde  = sp_mat( true, Location_12, vec_K1_tilde , 2*Nn, 3*Nn);
    
    M_K_tilde  = sp_mat( true, Location_22, vec_K_tilde , 3*Nn, 3*Nn);
    M_M_tilde  = sp_mat( true, Location_22, vec_M_tilde , 3*Nn, 3*Nn);
    M_D_tilde  = sp_mat( true, Location_22, vec_D_tilde , 3*Nn, 3*Nn);
    
    M_K1  = sp_mat( true, Location_21, vec_K1 , 3*Nn, 2*Nn);
    M_D1  = sp_mat( true, Location_21, vec_D1 , 3*Nn, 2*Nn);
    
    return 0;
}

void PML_solver::sub_fixed_dof()
{
    uword Nn = Table_Noeuds.n_cols;
    
    num_ddl_u = 0;
    num_ddl_E = 0;
    
    vec T_IDDL_U(2*Nn,fill::zeros);
    vec T_IDDL_E(3*Nn,fill::zeros);
    urowvec T_I_U(2*Nn,fill::zeros);
    urowvec T_J_U(2*Nn,fill::zeros);
    urowvec T_I_E(3*Nn,fill::zeros);
    urowvec T_J_E(3*Nn,fill::zeros);
    
    
    for(int k =0; k< Nn; k++)
    {
        if(Table_BC(k) == 0 )
        {
            T_IDDL_U(num_ddl_u) = 1.0;
            T_I_U(num_ddl_u) = 2*k ;
            T_J_U(num_ddl_u) = num_ddl_u;
            num_ddl_u++ ;
            
            T_IDDL_U(num_ddl_u) = 1.0;
            T_I_U(num_ddl_u) = 2*k+1 ;
            T_J_U(num_ddl_u) = num_ddl_u;
            num_ddl_u++ ;

        }
        
        if(Table_Noeuds(1,k) <= settings.ybot_pml|| abs(Table_Noeuds(0,k)-x_mid) >= settings.x_pml )
        {
            if (abs(Table_Noeuds(0,k)-x_mid) >= settings.x_pml )
            {
                T_IDDL_E(num_ddl_E) = 1.0;
                T_I_E(num_ddl_E) = 3*k ;
                T_J_E(num_ddl_E) = num_ddl_E;
                num_ddl_E++ ;
            }
            
            if (Table_Noeuds(1,k) <= settings.ybot_pml )
            {
                T_IDDL_E(num_ddl_E) = 1.0;
                T_I_E(num_ddl_E) = 3*k+1 ;
                T_J_E(num_ddl_E) = num_ddl_E;
                num_ddl_E++ ;
            }
            
            T_IDDL_E(num_ddl_E) = 1.0;
            T_I_E(num_ddl_E) = 3*k+2 ;
            T_J_E(num_ddl_E) = num_ddl_E;
            num_ddl_E++ ;
        }
        
    }
    
    T_IDDL_U.resize(num_ddl_u);
    T_I_U.resize(num_ddl_u);
    T_J_U.resize(num_ddl_u);
    
    
    T_IDDL_E.resize(num_ddl_E);
    T_I_E.resize(num_ddl_E);
    T_J_E.resize(num_ddl_E);

    umat Location_U = join_vert(T_I_U, T_J_U) ;
    M_DDL_U = sp_mat( Location_U, T_IDDL_U, 2*Nn, num_ddl_u);
    
    umat Location_E = join_vert(T_I_E, T_J_E) ;
    M_DDL_E = sp_mat( Location_E, T_IDDL_E, 3*Nn, num_ddl_E);
    
    M_K2 = M_DDL_U.t() * M_K2 * M_DDL_U;
    M_M  = M_DDL_U.t() * M_M  * M_DDL_U;
    M_D  = M_DDL_U.t() * M_D  * M_DDL_U;
    M_K  = M_DDL_U.t() * M_K  * M_DDL_U;
    M_D0 = M_DDL_U.t() * M_D0 * M_DDL_U;
    M_H  = M_DDL_U.t() * M_H  * M_DDL_U;
    
    M_H1_tilde = M_DDL_U.t() * M_H1_tilde * M_DDL_E;
    M_K1_tilde = M_DDL_U.t() * M_K1_tilde * M_DDL_E;
    M_K_tilde = M_DDL_E.t() * M_K_tilde * M_DDL_E;
    M_M_tilde = M_DDL_E.t() * M_M_tilde * M_DDL_E;
    M_D_tilde = M_DDL_E.t() * M_D_tilde * M_DDL_E;
    
    M_K1 = M_DDL_E.t() * M_K1 * M_DDL_U ;
    M_D1 = M_DDL_E.t() * M_D1 * M_DDL_U ;
    
    V_F = M_DDL_U.t() * V_F;
}

void PML_solver::set_external_force_parameters()
{
    double dist_essieux = 3.f;
    double Vtrain = settings.speed/3.6f;
    double ff = Vtrain/dist_essieux;
    f_1 = 10;
    omega = 2.f * M_PI * ff;
    double m = 4.f ;
    dt = 0.5f * M_PI / omega / m;
}

void PML_solver::set_frequency(double w)
{
    omega = w;
    double m = 4.f ;
    dt = 0.5f * M_PI / omega / m;
}

void PML_solver::set_Newmark()
{
    double ff = omega / M_PI / 2.f ;
    Force_time.zeros(settings.Niter);
    Force_time(0) = 0.f ;
    for (int k=1; k<settings.Niter; k++)
    {
        double t_k = dt * k ;
        Force_time(k) = f_1 * sin(omega * t_k ) * exp( - settings.a_force * pow( t_k * ff - 1.f , 2) ) ;
    }
    gamma = 0.5;
    beta = 0.25;
    Alpha_1 = 1.f/beta/(dt*dt);
    Alpha_2 = -1.f/beta/dt;
    Alpha_3 = -0.5f/beta;
    Alpha_4 = gamma/beta/dt;
    Alpha_5 = -gamma/beta;
    Alpha_6 = dt*(1.f-0.5f * gamma/beta);
}

void PML_solver::set_unit_force()
{
    uword Nn = Table_Noeuds.n_cols;
    ddl_force=0;
    {
        double dist = 1000 ;
        
        for (uword k=0; k<Nn; k++)
        {
            double x = Table_Noeuds(0,k) - x_mid ;
            double y = Table_Noeuds(1,k) - ytop;
            double bid =  x*x + y*y;
            if ( bid < dist)
            {
                dist = bid;
                ddl_force = 2*k+1;
            }
        }
    }

    V_F.zeros(2*Nn) ;
    V_F(ddl_force) = 1.0;
}

void PML_solver::sub_construction_global_Indice()
{
    vec T_Ind1 = ones<vec>(num_ddl_u) ;
    vec T_Ind2 = ones<vec>(num_ddl_E) ;
    urowvec T_I_1 = regspace<urowvec>(0,num_ddl_u-1) ;
    urowvec T_J_1 = regspace<urowvec>(0,num_ddl_u-1)  ;
    urowvec T_I_2 = regspace<urowvec>(0,num_ddl_E-1)  ;
    urowvec T_J_2 = regspace<urowvec>(num_ddl_u,num_ddl_u+num_ddl_E-1)  ;
    
    umat Location_1 = join_vert(T_I_1, T_J_1) ;
    umat Location_2 = join_vert(T_I_2, T_J_2) ;
    
    M_Ind1 = sp_mat( Location_1, T_Ind1, num_ddl_u, num_ddl_u+num_ddl_E);
    M_Ind2 = sp_mat( Location_2, T_Ind2, num_ddl_E, num_ddl_u+num_ddl_E);
}

void PML_solver::set_matrix_system()
{
    uword Icase = settings.Icase;
    
    if (Icase == 4 || Icase == 5 || Icase == 6)
        set_Newmark();
    
    if (Icase == 1 || Icase == 3)
    {
        M_K_dyn = sp_cx_mat( M_K -(omega * omega) * M_M,  omega * M_D0 ) ;
        
        if (Icase == 3)
            M_K_dyn += omega *  sp_cx_mat( sp_mat(num_ddl_u,num_ddl_u), M_D) ;
        
        V_F_pml_complex = cx_vec( V_F , zeros<vec>(num_ddl_u));
    }
    
    if (Icase == 4 || Icase == 6)
    {
        M_K_eff = Alpha_1 * M_M + M_K + Alpha_4 * M_D0 ;
        
        if (Icase == 6)
            M_K_eff += Alpha_4 * M_D ;
        
        V_F_pml_real = V_F ;
    }
    
    if (Icase == 2 || Icase == 5)
    {
        sub_construction_global_Indice();
        
        M_Mglob = M_Ind1.t() * M_M * M_Ind1
                                            + M_Ind2.t()  * M_M_tilde * M_Ind2 ;
        
        M_Dglob = M_Ind1.t()  * (M_D + M_D0) * M_Ind1
                + M_Ind2.t()  * M_D1         * M_Ind1 + M_Ind2.t() * M_D_tilde * M_Ind2 ;
        
        M_Kglob = M_Ind1.t()  * (M_K+M_K2) * M_Ind1 + M_Ind1.t() * M_K1_tilde * M_Ind2
                + M_Ind2.t()  * M_K1       * M_Ind1 + M_Ind2.t() * M_K_tilde  * M_Ind2 ;
        
        M_Hglob = M_Ind1.t()  * M_H * M_Ind1 + M_Ind1.t()  * M_H1_tilde * M_Ind2;
        
        if (Icase == 2)
        {
            M_K_dyn =   sp_cx_mat( M_Kglob - (omega * omega) * M_Mglob ,  omega * M_Dglob  -1.f/omega * M_Hglob );
            V_F_pml_complex =  cx_vec ( M_Ind1.t() * V_F, vec(num_ddl_u + num_ddl_E, fill::zeros) ) ;
            
            if (settings.debug & DEBUG_PRINT_GLOBAL_MATRICES)
            {
                save_sp2ascii(M_K_dyn, "matlab/DATA_K_dyn.ascii");
                save_sp2ascii(M_DDL_U, "matlab/DATA_M_DDL_U.ascii");
                save_sp2ascii(M_DDL_E, "matlab/DATA_M_DDL_E.ascii");
                
                save_vec2ascii(V_F_pml_complex, "matlab/DATA_V_F_pml_complex.ascii");
                
                save_mat2ascii(  Table_Noeuds, "matlab/Table_Noeuds.ascii");
                save_mat2ascii(Table_Elements, "matlab/Table_Elements.ascii");
                
                save_int2ascii(num_ddl_u, "matlab/num_ddl_u.ascii");
                save_int2ascii(num_ddl_E, "matlab/num_ddl_E.ascii");
            }

        }
        else
        {
            M_K_eff = Alpha_1 * M_Mglob + Alpha_4 * M_Dglob  + M_Kglob + 0.5 * dt *  M_Hglob;
            V_F_pml_real = M_Ind1.t() * V_F;
        }
    }
}

uword PML_solver::get_Icase()
{
    return settings.Icase;
}

double PML_solver::get_final_time()
{
    return settings.Niter * dt ;
}

int PML_solver::solve()
{
    uword Icase = settings.Icase ;
    
    if (Icase == 1 || Icase == 2 || Icase == 3)
    {
        cx_vec  U_pml;
        umfpack_wrapper solver(M_K_dyn);
        U_pml = solver.solve(V_F_pml_complex);
        
        if (Icase == 2)
            Ustock_pml_complex = cx_vec(M_DDL_U * M_Ind1 * real(U_pml), M_DDL_U * M_Ind1 * imag(U_pml)) ;
        
        else if (Icase == 1|| Icase == 3)
            Ustock_pml_complex = cx_vec(M_DDL_U * real(U_pml), M_DDL_U * imag(U_pml));
    }
    
    if (Icase == 4 || Icase == 5 || Icase == 6)
    {
        
        long Nn = Table_Noeuds.n_cols;
        Ustock_pml_real = zeros<vec>( 2*Nn , settings.Niter);
        
        long NDDL = M_K_eff.n_rows;
        vec Feff, U_pml = zeros<vec>(NDDL), V_pml = zeros<vec>(NDDL), A_pml = zeros<vec>(NDDL) ;
        
        double fnew = 0.f, fold = 0.f ;
        
        umfpack_wrapper solver(M_K_eff);
        
        for (int k=1; k<settings.Niter; k++)
        {
            printProgress( (double)k/(double)settings.Niter);
            
            fold = fnew ;
            fnew = Force_time(k);
            if (Icase == 4)
                Feff    = (fnew-fold) * V_F_pml_real - M_M * ( Alpha_2 * V_pml + Alpha_3 * A_pml) ;
            
            if (Icase == 5)
                Feff    = (fnew-fold) * V_F_pml_real - M_Mglob * ( Alpha_2 * V_pml + Alpha_3 * A_pml)
                - M_Dglob * ( Alpha_5 * V_pml + Alpha_6 * A_pml)
                - dt * M_Hglob * U_pml ;
            
            if (Icase == 6)
                Feff    = (fnew-fold) * V_F_pml_real - M_M * ( Alpha_2 * V_pml + Alpha_3 * A_pml)
                - (M_D  + M_D0) * ( Alpha_5 * V_pml + Alpha_6 * A_pml)  ;
            
            vec DU = solver.solve(Feff);
            
            vec DV = Alpha_4*DU + Alpha_5*V_pml + Alpha_6*A_pml;
            vec DA = Alpha_1*DU + Alpha_2*V_pml + Alpha_3*A_pml;
            
            U_pml = U_pml + DU;
            V_pml = V_pml + DV;
            A_pml = A_pml + DA;
            
            
            if (Icase == 5)
                Ustock_pml_real.col(k) = M_DDL_U * M_Ind1 * U_pml;
            
            else if (Icase == 4|| Icase == 6)
                Ustock_pml_real.col(k) = M_DDL_U * U_pml;
        }
    }
    
    cout << endl ;
    return 0;
}

void PML_solver::save_vtk()
{
    if (settings.problem_type == DETERMINIST)
        save_vtk_determinist();
    
    if (settings.problem_type == RANDOM)
        save_vtk_random();
}

void PML_solver::save_vtk_determinist()
{
    uword Ne = Table_Elements.n_rows;
    uword Nn = Table_Noeuds.n_cols;
    
    uword  Icase = settings.Icase ;
    
    if (Icase == 4)
    {
        ofstream output;
        string filename("RESULTS_ref_transient.bin");
        
        chrono_message("Ecriture du fichier " +filename +"\n");
        
        output.open (filename,ios::out | ios::binary| ios::trunc);
        output.write( (char*)Ustock_pml_real.memptr(),sizeof(double) * Ustock_pml_real.size() );
        output.close();
    }
    
    mat Ustock_ref_real;
    
    if (Icase == 5 || Icase == 6)
    {
        long Nb = Ustock_pml_real.size() ;
        Ustock_ref_real.zeros(Ustock_pml_real.n_rows,Ustock_pml_real.n_cols);
        string filename("RESULTS_ref_transient.bin");
        std::ifstream input( filename, std::ios::binary  );
        input.read((char*)Ustock_ref_real.memptr(),Nb*sizeof(double));
        input.close();
    }
    //----------------- Sauvegarde des fichier VTK ----------------
    
    if (Icase == 1 || Icase == 2 || Icase == 3)
    {
        chrono_message("Ecriture du fichier " + settings.FILEOUT +".vtk\n") ;
        
        ofstream fid;
        fid.open (settings.FILEOUT+".vtk",ios::out);
        
        fid << "# vtk DataFile Version 2.0" << endl ;
        
        if (Icase == 3)
            fid << "Resultats etude CAL en frequence " << endl ;
        
        if (Icase == 2)
            fid << "Resultats etude PML en frequence " << endl ;
        
        if (Icase == 1)
            fid << "Resultats etude Référence en frequence " << endl ;
        
        fid << "ASCII" << endl << endl ;
        fid << "DATASET UNSTRUCTURED_GRID" << endl ;
        fid << "POINTS " << Nn << " double" << endl ;
        for (int j=0; j< Nn; j++)
            fid << Table_Noeuds(0,j) << " " << Table_Noeuds(1,j) << " 0" << endl ;
        
        fid << endl ;
        fid << "CELLS " << Ne << " " << 7*Ne << endl ;
        for(int e=0; e<Ne; e++)
            fid << "6 " << Table_Elements(e,1)-1 << " " << Table_Elements(e,2)-1 << " "
            << Table_Elements(e,3)-1 << " " << Table_Elements(e,4)-1 << " "
            << Table_Elements(e,5)-1 << " " << Table_Elements(e,6)-1 << endl ;
        fid << endl ;
        fid << "CELL_TYPES " << Ne << endl ;
        for(int e=0; e<Ne; e++)
            fid << "22" << endl ;
        
        fid << endl ;
        
        fid << "POINT_DATA " << Nn << endl ;
        fid << "VECTORS Disp double" << endl ;
        for (int k=0; k< Nn; k++)
            fid << "    " << real(Ustock_pml_complex(2*k,0)) << " " << real(Ustock_pml_complex(2*k+1,0)) << " 0 " << endl ;
        
        fid << endl ;
        
        fid.close();
    }
    
    if (Icase == 4 || Icase == 5 || Icase == 6)
    {
        chrono_message("Ecriture des fichiers " + settings.FILEOUT  + "XYZ.vtk\n");
        
        for (int k=0; k<settings.Niter; k++)
        {
            printProgress( (double)k/(double)settings.Niter);
            
            char filename[100];
            
            sprintf (filename, "%s%d.vtk", settings.FILEOUT.c_str(), k );
            
            if (k < 1000)
                sprintf (filename, "%s0%d.vtk",settings.FILEOUT.c_str(),  k);
            
            if (k < 100)
                sprintf (filename, "%s00%d.vtk",settings.FILEOUT.c_str(),  k );
            
            if (k < 10 )
                sprintf (filename, "%s000%d.vtk", settings.FILEOUT.c_str(), k );
            
            ofstream fid;
            fid.open (filename,ios::out);
            
            fid << "# vtk DataFile Version 2.0" << endl ;
            
            if (Icase == 6)
                fid << "Resultats etude CAL en transitoire " << endl ;
            
            if (Icase == 5)
                fid << "Resultats etude PML en transitoire " << endl ;
            
            if (Icase == 4)
                fid << "Resultats etude Référence en transitoire " << endl ;
            
            fid << "ASCII" << endl << endl ;
            fid << "DATASET UNSTRUCTURED_GRID" << endl ;
            fid << "POINTS " << Nn << " double" << endl ;
            for (int j=0; j< Nn; j++)
                fid << Table_Noeuds(0,j) << " " << Table_Noeuds(1,j) << " 0" << endl ;
            
            fid << endl ;
            fid << "CELLS " << Ne << " " << 7*Ne << endl ;
            for(int e=0; e<Ne; e++)
                fid << "6 " << Table_Elements(e,1)-1 << " " << Table_Elements(e,2)-1 << " "
                << Table_Elements(e,3)-1 << " " << Table_Elements(e,4)-1 << " "
                << Table_Elements(e,5)-1 << " " << Table_Elements(e,6)-1 << endl ;
            fid << endl ;
            fid << "CELL_TYPES " << Ne << endl ;
            for(int e=0; e<Ne; e++)
                fid << "22" << endl ;
            
            fid << endl ;
            
            fid << "POINT_DATA " << Nn << endl ;
            fid << "VECTORS Disp double" << endl ;
            for (int j=0; j< Nn; j++)
                fid << "    " << Ustock_pml_real(2*j,k) << " " << Ustock_pml_real(2*j+1,k) << " 0 " << endl ;
            
            fid << endl ;
            
            if (Icase == 5 || Icase == 6) {
                Ustock_ref_real = abs(Ustock_pml_real - Ustock_ref_real);
                double normalization = Ustock_ref_real.max();
                Ustock_ref_real = Ustock_ref_real/normalization;
                fid << "VECTORS Difference double" << endl;
                for (int j=0; j< Nn; j++)
                    fid << "    " << Ustock_ref_real(2*j,k) << " " << Ustock_ref_real(2*j+1,k) << " 0 " << endl ;
                fid << endl ;
            }
            
            fid.close();
            
        }
        cout << endl;
    }
}

void PML_solver::sample_random_settings()
{
    // Generate new samples of "settings"
    //
    // TO DO
}

void PML_solver::save_vtk_random()
{
    long Ne = Table_Elements.n_rows;
    long Nn = Table_Noeuds.n_cols;
    
    uword  Icase = settings.Icase ;
    
    // TO DO
}






