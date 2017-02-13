//
//  FEM.cpp
//

#include "PML_solver.hpp"
#include "tools.hpp"
#include "Tokenizer.hpp"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include "umfpack_wrapper.hpp"

using namespace std;
using namespace Eigen;

using namespace FEM;


/* -------------------------
 
 CLASS PML_solver
 
 --------------------------- */



void FEM::save_sp2matlab(SP_mat &MAT, string filename)
{
    ofstream fid;
    fid.open (filename,ios::out);
    cout << filename << " " << MAT.nonZeros() <<  " non zeros" << endl;
    for (int k=0; k<MAT.outerSize(); ++k)
        for (SP_mat::InnerIterator it(MAT,k); it; ++it)
            fid<< it.value() << " " <<  it.row()+1 << " " << it.col()+1<< endl;
    fid.close();
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
    Matrix<unsigned,Dynamic, 1> Check_ref;
    
    input.open ( settings.msh_file.c_str ( ) );
    while (!input.eof())
    {
        getline(input, buffer);
        if (Tokenizer(buffer).next() == "$Nodes")
        {
            getline(input, buffer);
            
            NumberNode = stoi(Tokenizer(buffer).next() );
            Table_Noeuds =  Matrix<double,2,Dynamic>::Zero(2,NumberNode);
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
            Table_Noeuds.conservativeResize (2, NumberNode );
            Table_BC = Matrix<unsigned,Dynamic, 1>::Zero(NumberNode,1);
            Check_ref = Matrix<unsigned,Dynamic, 1>::Zero(NumberNode,1);
        }
        else if (Tokenizer(buffer).next() == "$Elements")
        {
            getline(input, buffer);
            NumberElement = stoi(Tokenizer(buffer).next() );
            Table_Elements =  Matrix<unsigned,Dynamic,7>::Zero(NumberElement,7);
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
                    unsigned tag = (unsigned) stoi(Tokens.next());
                    
                    for (int k=2; k<nb_tag; k++)
                        stoi(Tokens.next());
                    
                    if (type == element_type)
                    {
                        Table_Elements(NumberElement,0) = tag;
                        //cout << tag << endl;
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
            Table_Elements.conservativeResize ( NumberElement,  7 );
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
    long Ne = Table_Elements.rows();
    
    for(int e=0; e < Ne; e++)
    {
        Vector2d X1, X2, X3;
        unsigned num1, num2, num3;
        num1 = Table_Elements(e,1);
        num2 = Table_Elements(e,2);
        num3 = Table_Elements(e,3);
        
        X1 = Table_Noeuds.col(num1-1);
        X2 = Table_Noeuds.col(num2-1);
        X3 = Table_Noeuds.col(num3-1);
        
        double test =(X2(0)-X1(0))*(X3(1)-X1(1)) - (X2(1)-X1(1))*(X3(0)-X1(0)) ;
        if (test < 0)
        {
            unsigned N;
            N = Table_Elements(e,3);
            Table_Elements(e,3) = Table_Elements(e,2);
            Table_Elements(e,2) = N;
            N = Table_Elements(e,6);
            Table_Elements(e,6) = Table_Elements(e,4);
            Table_Elements(e,4) = N;
            //cout << "Reorientation of element " << e << endl;
        }
    }
}


void PML_solver::set_xmid_ytop()
{
    long Nn = Table_Noeuds.cols();
    
    double x_left  = 1000;
    double x_right = -1000;
    ytop = -1000;
    for (unsigned k=0; k<Nn; k++)
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
    long Ne = Table_Elements.rows();
    long Nn = Table_Noeuds.cols();
    
    unsigned Np = 6;
    
    double a = 0.445948490915965;
    double b = 0.091576213509771;
    
    Matrix<double,2,6> Liste_Yp;
    Liste_Yp <<  b, 1-2*b, b    , a    , a, 1-2*a ,
    b, b    , 1-2*b, 1-2*a, a, a     ;
    
    double P1 = 0.11169079483905;
    double P2 = 0.0549758718227661;
    
    Matrix<double,1,6> Liste_gammap;
    Liste_gammap <<  P2, P2, P2, P1, P1, P1 ;
    
    Matrix2d M_I2 = Matrix2d::Identity();
    
    Matrix3d M_I3= Matrix3d::Identity();
    M_I3(2,2) = sqrt(2);
    
    double s2 = sqrt(2)/2.f;
    double s3 = 0.5;
    double s4 = 0.25;
    
    Matrix<double,3,4>  M_P;
    M_P <<
    1,  0,  0, 0 ,
    0,  0,  0, 1 ,
    0, s2, s2, 0 ;
    
    Matrix<double,6,6>   C_coeff;
    C_coeff <<
    1,  0,  0,  0,  0,  0 ,
    1,  1,  0,  1,  0,  0 ,
    1,  0,  1,  0,  1,  0 ,
    1, s3,  0, s4,  0,  0 ,
    1, s3, s3, s4, s4, s4 ,
    1,  0, s3,  0, s4,  0 ;
    
    ColPivHouseholderQR<Matrix<double,6,6>> dec_C_coeff(C_coeff);
    
    VectorXd V(6);
    
    V.setZero();
    V(0) = 1;
    VectorXd V_coeff = dec_C_coeff.solve(V);
    
    double a1 = V_coeff(0);
    double b1 = V_coeff(1);
    double c1 = V_coeff(2);
    double d1 = V_coeff(3);
    double e1 = V_coeff(4);
    double f1 = V_coeff(5);
    
    V.setZero();
    V(1) = 1;
    V_coeff = dec_C_coeff.solve(V);
    
    double a2 = V_coeff(0);
    double b2 = V_coeff(1);
    double c2 = V_coeff(2);
    double d2 = V_coeff(3);
    double e2 = V_coeff(4);
    double f2 = V_coeff(5);
    
    V.setZero();
    V(2) = 1;
    V_coeff = dec_C_coeff.solve(V);
    
    double a3 = V_coeff(0);
    double b3 = V_coeff(1);
    double c3 = V_coeff(2);
    double d3 = V_coeff(3);
    double e3 = V_coeff(4);
    double f3 = V_coeff(5);
    
    V.setZero();
    V(3) = 1;
    V_coeff = dec_C_coeff.solve(V);
    
    double a4 = V_coeff(0);
    double b4 = V_coeff(1);
    double c4 = V_coeff(2);
    double d4 = V_coeff(3);
    double e4 = V_coeff(4);
    double f4 = V_coeff(5);
    
    V.setZero();
    V(4) = 1;
    V_coeff = dec_C_coeff.solve(V);
    
    double a5 = V_coeff(0);
    double b5 = V_coeff(1);
    double c5 = V_coeff(2);
    double d5 = V_coeff(3);
    double e5 = V_coeff(4);
    double f5 = V_coeff(5);
    
    V.setZero();
    V(5) = 1;
    
    V_coeff = dec_C_coeff.solve(V);
    
    double a6 = V_coeff(0);
    double b6 = V_coeff(1);
    double c6 = V_coeff(2);
    double d6 = V_coeff(3);
    double e6 = V_coeff(4);
    double f6 = V_coeff(5);
    
    vector<T> T_K2, T_K, T_M, T_D, T_D0, T_H, T_H1_tilde, T_K1_tilde, T_K_tilde, T_M_tilde, T_D_tilde, T_K1, T_D1;
    
    T_K2.reserve(12*12*Ne);
    T_K.reserve(12*12*Ne);
    T_M.reserve(12*12*Ne);
    T_D.reserve(12*12*Ne);
    T_D0.reserve(12*12*Ne);
    T_H.reserve(12*12*Ne);
    T_H1_tilde.reserve(12*18*Ne);
    T_K1_tilde.reserve(12*18*Ne);
    T_K_tilde.reserve(18*18*Ne);
    T_M_tilde.reserve(18*18*Ne);
    T_D_tilde.reserve(18*18*Ne);
    T_K1.reserve(12*18*Ne);
    T_D1.reserve(12*18*Ne);
    
    for(int e = 0; e<Ne; e++)
    {
        unsigned n1 = Table_Elements(e,1);
        unsigned n2 = Table_Elements(e,2);
        unsigned n3 = Table_Elements(e,3);
        unsigned n4 = Table_Elements(e,4);
        unsigned n5 = Table_Elements(e,5);
        unsigned n6 = Table_Elements(e,6);
        
        Vector2d Xe1 = Table_Noeuds.col(n1-1);
        Vector2d Xe2 = Table_Noeuds.col(n2-1);
        Vector2d Xe3 = Table_Noeuds.col(n3-1);
        Vector2d Xe4 = Table_Noeuds.col(n4-1);
        Vector2d Xe5 = Table_Noeuds.col(n5-1);
        Vector2d Xe6 = Table_Noeuds.col(n6-1);
        
        VectorXd Xe(12);
        Xe  << Xe1 , Xe2 , Xe3 , Xe4 , Xe5 , Xe6 ;
        
        Matrix<double,6,6> C_coeff2;
        C_coeff2 <<
        1, Xe1(0), Xe1(1), Xe1(0)*Xe1(0), Xe1(1)*Xe1(1), Xe1(0)*Xe1(1) ,
        1, Xe2(0), Xe2(1), Xe2(0)*Xe2(0), Xe2(1)*Xe2(1), Xe2(0)*Xe2(1) ,
        1, Xe3(0), Xe3(1), Xe3(0)*Xe3(0), Xe3(1)*Xe3(1), Xe3(0)*Xe3(1) ,
        1, Xe4(0), Xe4(1), Xe4(0)*Xe4(0), Xe4(1)*Xe4(1), Xe4(0)*Xe4(1) ,
        1, Xe5(0), Xe5(1), Xe5(0)*Xe5(0), Xe5(1)*Xe5(1), Xe5(0)*Xe5(1) ,
        1, Xe6(0), Xe6(1), Xe6(0)*Xe6(0), Xe6(1)*Xe6(1), Xe6(0)*Xe6(1) ;
        
        
        ColPivHouseholderQR<Matrix<double,6,6>> dec_C_coeff2(C_coeff2);
        
        VectorXd V_coeff_tile(6);
        
        V.setZero();
        V(0) = 1;
        V_coeff_tile = dec_C_coeff2.solve(V);
        
        double atilde_1 = V_coeff_tile(0);
        double btilde_1 = V_coeff_tile(1);
        double ctilde_1 = V_coeff_tile(2);
        double dtilde_1 = V_coeff_tile(3);
        double etilde_1 = V_coeff_tile(4);
        double ftilde_1 = V_coeff_tile(5);
        
        
        V.setZero();
        V(1) = 1;
        V_coeff_tile = dec_C_coeff2.solve(V);
        
        double atilde_2 = V_coeff_tile(0);
        double btilde_2 = V_coeff_tile(1);
        double ctilde_2 = V_coeff_tile(2);
        double dtilde_2 = V_coeff_tile(3);
        double etilde_2 = V_coeff_tile(4);
        double ftilde_2 = V_coeff_tile(5);
        
        V.setZero();
        V(2) = 1;
        V_coeff_tile = dec_C_coeff2.solve(V);
        
        double atilde_3 = V_coeff_tile(0);
        double btilde_3 = V_coeff_tile(1);
        double ctilde_3 = V_coeff_tile(2);
        double dtilde_3 = V_coeff_tile(3);
        double etilde_3 = V_coeff_tile(4);
        double ftilde_3 = V_coeff_tile(5);
        
        V.setZero();
        V(3) = 1;
        V_coeff_tile = dec_C_coeff2.solve(V);
        
        double atilde_4 = V_coeff_tile(0);
        double btilde_4 = V_coeff_tile(1);
        double ctilde_4 = V_coeff_tile(2);
        double dtilde_4 = V_coeff_tile(3);
        double etilde_4 = V_coeff_tile(4);
        double ftilde_4 = V_coeff_tile(5);
        
        V.setZero();
        V(4) = 1;
        V_coeff_tile = dec_C_coeff2.solve(V);
        
        double atilde_5 = V_coeff_tile(0);
        double btilde_5 = V_coeff_tile(1);
        double ctilde_5 = V_coeff_tile(2);
        double dtilde_5 = V_coeff_tile(3);
        double etilde_5 = V_coeff_tile(4);
        double ftilde_5 = V_coeff_tile(5);
        
        V.setZero();
        V(5) = 1;
        V_coeff_tile = dec_C_coeff2.solve(V);
        
        double atilde_6 = V_coeff_tile(0);
        double btilde_6 = V_coeff_tile(1);
        double ctilde_6 = V_coeff_tile(2);
        double dtilde_6 = V_coeff_tile(3);
        double etilde_6 = V_coeff_tile(4);
        double ftilde_6 = V_coeff_tile(5);
        
        MatrixXd M_K2e=MatrixXd::Zero(12,12);
        MatrixXd M_Ke=MatrixXd::Zero(12,12);
        MatrixXd M_Me=MatrixXd::Zero(12,12);
        
        MatrixXd M_De=MatrixXd::Zero(12,12);
        MatrixXd M_He=MatrixXd::Zero(12,12);
        MatrixXd M_D0e=MatrixXd::Zero(12,12);
        
        MatrixXd M_Me_tilde=MatrixXd::Zero(18,18);
        MatrixXd M_De_tilde=MatrixXd::Zero(18,18);
        MatrixXd M_Ke_tilde=MatrixXd::Zero(18,18);
        
        MatrixXd M_K1e_tilde=MatrixXd::Zero(12,18);
        MatrixXd M_H1e_tilde=MatrixXd::Zero(12,18);
        
        MatrixXd M_K1e=MatrixXd::Zero(18,12);
        MatrixXd M_D1e=MatrixXd::Zero(18,12);
        
        for(int p = 0; p<Np; p++)
        {
            
            Vector2d Yp = Liste_Yp.col(p);
            double gammap = Liste_gammap(p);
            
            Matrix<double,2,12> M_EYp;
            
            {
                double E1Yp = a1 + b1*Yp(0) + c1*Yp(1) + d1*Yp(0)*Yp(0) + e1*Yp(1)*Yp(1) + f1*Yp(0)*Yp(1);
                double E2Yp = a2 + b2*Yp(0) + c2*Yp(1) + d2*Yp(0)*Yp(0) + e2*Yp(1)*Yp(1) + f2*Yp(0)*Yp(1);
                double E3Yp = a3 + b3*Yp(0) + c3*Yp(1) + d3*Yp(0)*Yp(0) + e3*Yp(1)*Yp(1) + f3*Yp(0)*Yp(1);
                double E4Yp = a4 + b4*Yp(0) + c4*Yp(1) + d4*Yp(0)*Yp(0) + e4*Yp(1)*Yp(1) + f4*Yp(0)*Yp(1);
                double E5Yp = a5 + b5*Yp(0) + c5*Yp(1) + d5*Yp(0)*Yp(0) + e5*Yp(1)*Yp(1) + f5*Yp(0)*Yp(1);
                double E6Yp = a6 + b6*Yp(0) + c6*Yp(1) + d6*Yp(0)*Yp(0) + e6*Yp(1)*Yp(1) + f6*Yp(0)*Yp(1);
                
                M_EYp <<  M_I2*E1Yp , M_I2*E2Yp, M_I2*E3Yp, M_I2*E4Yp, M_I2*E5Yp, M_I2*E6Yp ;
                
            }
            
            Vector2d Xp = M_EYp * Xe ;
            
            Matrix<double,3,18> M_E_PML_Xp;
            Matrix<double,2,12> M_EXp;
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
                
                M_EXp << E1Xp*M_I2, E2Xp*M_I2, E3Xp*M_I2, E4Xp*M_I2, E5Xp*M_I2, E6Xp*M_I2 ;
                
                M_E_PML_Xp << E1Xp*M_I3, E2Xp*M_I3, E3Xp*M_I3, E4Xp*M_I3, E5Xp*M_I3, E6Xp*M_I3 ;
                
            }
            
            
            Matrix<double,4,12> M_Ce;
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
                
                Matrix<double,2,12> M_Ce1;
                M_Ce1 << A1*M_I2, A2*M_I2, A3*M_I2, A4*M_I2, A5*M_I2, A6*M_I2;
                
                Matrix<double,2,12> M_Ce2;
                M_Ce2 << B1*M_I2, B2*M_I2, B3*M_I2, B4*M_I2, B5*M_I2, B6*M_I2;
                
                M_Ce << M_Ce1 , M_Ce2 ;
                
            }
            
            
            Matrix<double,3,12> M_Btilde  = M_P * M_Ce ;
            
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
                
                Vector2d Phi1 = A1*Xe1 + A2*Xe2 +  A3*Xe3  + A4*Xe4 + A5*Xe5 + A6*Xe6;
                Vector2d Phi2 = B1*Xe1 + B2*Xe2 +  B3*Xe3  + B4*Xe4 + B5*Xe5 + B6*Xe6;
                JeYp = (MatrixXd(2,2) << Phi1, Phi2).finished().determinant() ;
                
                if (JeYp < 0)
                {
                    cout << endl << "Jacobian is negative" << endl ;
                    return -1 ;
                }
            }
            
            unsigned num_mat = Table_Elements(e,0);
            
            double E,v,rho,Xi;
            
            Vector4d parameters = settings.Mat[num_mat];
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
            
            Matrix<double,3,4> M_Phat;
            M_Phat <<
            mu2_PML, 0         , 0         , 0       ,
            0      , 0         , 0         , mu1_PML ,
            0      , s2*mu2_PML, s2*mu1_PML, 0       ;
            
            Matrix<double,3,4> M_PDEpsilon ;
            M_PDEpsilon <<
            mu1_PML, 0          , 0         , 0       ,
            0      , 0          , 0         , mu2_PML ,
            0      , s2*mu1_PML , s2*mu2_PML, 0       ;
            
            Matrix<double,3,12>  M_Btilde_DEpsilon = M_PDEpsilon * M_Ce ;
            Matrix<double,3,12>  M_Btildehat = M_Phat * M_Ce ;
            
            double Lambda = (E*v)/(1.f+v)/(1.f-2.f*v);
            double Mu = 0.5*E/(1.f+v);
            
            Matrix3d M_Ae;
            M_Ae << Lambda+2*Mu, Lambda      ,0       ,
            Lambda     , Lambda+2*Mu , 0      ,
            0          , 0           , 4.f*Mu ;
            
            M_Ke  = M_Ke + gammap *  M_Btilde.transpose() * M_Ae * M_Btilde * JeYp ;
            
            M_K2e  = M_K2e + mu1_PML * mu2_PML * gammap * rho * M_EXp.transpose() * M_EXp * JeYp;
            
            M_Me  = M_Me + gammap * rho * M_EXp.transpose() * M_EXp * JeYp;
            
            M_De  = M_De +  (mu1_PML + mu2_PML) * gammap * rho * M_EXp.transpose() * M_EXp * JeYp;
            
            M_He =  M_He +  gammap * M_Btildehat.transpose() * M_Ae * M_Btilde * JeYp;
            
            M_H1e_tilde = M_H1e_tilde + gammap * M_Btildehat.transpose() * M_Ae * M_E_PML_Xp * JeYp ;
            
            M_D0e  = M_D0e + 2.f * Xi * sqrt(E * rho) * gammap  * M_EXp.transpose() * M_EXp * JeYp;
            
            M_Ke_tilde = M_Ke_tilde + gammap * mu1_PML * mu2_PML * M_E_PML_Xp.transpose() * M_E_PML_Xp * JeYp;
            
            M_De_tilde = M_De_tilde + gammap * (mu1_PML + mu2_PML) * M_E_PML_Xp.transpose() * M_E_PML_Xp * JeYp;
            
            M_Me_tilde = M_Me_tilde + gammap * M_E_PML_Xp.transpose() * M_E_PML_Xp * JeYp;
            
            M_K1e =  M_K1e + gammap *  mu1_PML * mu2_PML * M_E_PML_Xp.transpose() * M_Btilde * JeYp ;
            
            M_K1e_tilde =  M_K1e_tilde +  gammap * M_Btilde.transpose() * M_Ae * M_E_PML_Xp * JeYp ;
            
            M_D1e = M_D1e + gammap * M_E_PML_Xp.transpose() * M_Btilde_DEpsilon * JeYp;
            
        }
        
        Matrix<unsigned,6,1> T_num_noeuds(Table_Elements.block<1,6>(e,1));
        
        Matrix<int,12,1> T_Nek;
        T_Nek <<
        (T_num_noeuds(0) - 1)*2, 1+(T_num_noeuds(0) - 1)*2,
        (T_num_noeuds(1) - 1)*2, 1+(T_num_noeuds(1) - 1)*2,
        (T_num_noeuds(2) - 1)*2, 1+(T_num_noeuds(2) - 1)*2,
        (T_num_noeuds(3) - 1)*2, 1+(T_num_noeuds(3) - 1)*2,
        (T_num_noeuds(4) - 1)*2, 1+(T_num_noeuds(4) - 1)*2,
        (T_num_noeuds(5) - 1)*2, 1+(T_num_noeuds(5) - 1)*2;
        
        Matrix<int,18,1> T_Nek_PML;
        T_Nek_PML <<
        (T_num_noeuds(0) - 1)*3, 1+(T_num_noeuds(0) - 1)*3, 2+(T_num_noeuds(0) - 1)*3,
        (T_num_noeuds(1) - 1)*3, 1+(T_num_noeuds(1) - 1)*3, 2+(T_num_noeuds(1) - 1)*3,
        (T_num_noeuds(2) - 1)*3, 1+(T_num_noeuds(2) - 1)*3, 2+(T_num_noeuds(2) - 1)*3,
        (T_num_noeuds(3) - 1)*3, 1+(T_num_noeuds(3) - 1)*3, 2+(T_num_noeuds(3) - 1)*3,
        (T_num_noeuds(4) - 1)*3, 1+(T_num_noeuds(4) - 1)*3, 2+(T_num_noeuds(4) - 1)*3,
        (T_num_noeuds(5) - 1)*3, 1+(T_num_noeuds(5) - 1)*3, 2+(T_num_noeuds(5) - 1)*3;
        
        for (int i=0; i<12; i++)
        {
            int I = T_Nek(i);
            
            for (int j=0; j<12; j++)
            {
                int J = T_Nek(j);
                
                if (M_K2e (i,j) != 0)
                    T_K2.push_back(T(I,J,M_K2e(i,j)));
                
                if (M_Me (i,j) != 0)
                    T_M.push_back (T(I,J,M_Me (i,j)));
                
                if (M_De (i,j) != 0)
                    T_D.push_back (T(I,J,M_De (i,j)));
                
                if (M_Ke (i,j) != 0)
                    T_K.push_back (T(I,J,M_Ke (i,j)));
                
                if (M_He (i,j) != 0)
                    T_H.push_back (T(I,J,M_He (i,j)));
                
                if (M_D0e (i,j) != 0)
                    T_D0.push_back(T(I,J,M_D0e(i,j)));
                
            }
            
            for (int j=0; j<18; j++)
            {
                
                int J = T_Nek_PML(j);
                
                if (M_H1e_tilde (i,j) != 0)
                    T_H1_tilde.push_back (T(I,J,M_H1e_tilde (i,j)));
                
                if (M_K1e_tilde (i,j) != 0)
                    T_K1_tilde.push_back (T(I,J,M_K1e_tilde (i,j)));
            }
        }
        
        for (int i=0; i<18; i++)
        {
            int I = T_Nek_PML(i);
            
            for (int j=0; j<18; j++)
            {
                int J = T_Nek_PML(j);
                
                if (M_Ke_tilde (i,j) != 0)
                    T_K_tilde.push_back (T(I,J,M_Ke_tilde (i,j)));
                
                if (M_Me_tilde (i,j) != 0)
                    T_M_tilde.push_back (T(I,J,M_Me_tilde (i,j)));
                
                if (M_De_tilde (i,j) != 0)
                    T_D_tilde.push_back (T(I,J,M_De_tilde (i,j)));
            }
            
            for (int j=0; j<12; j++)
            {
                int J = T_Nek(j);
                
                if (M_K1e (i,j) != 0)
                    T_K1.push_back (T(I,J,M_K1e (i,j)));
                
                if (M_D1e (i,j) != 0)
                    T_D1.push_back (T(I,J,M_D1e (i,j)));
            }
        }
    }
    
    M_K2.resize(Nn*2,Nn*2);
    M_M.resize(Nn*2,Nn*2);
    M_D.resize(Nn*2,Nn*2);
    M_D0.resize(Nn*2,Nn*2);
    M_K.resize(Nn*2,Nn*2);
    M_H.resize(Nn*2,Nn*2);
    
    M_H1_tilde.resize(Nn*2,Nn*3);
    M_K1_tilde.resize(Nn*2,Nn*3);
    
    M_K_tilde.resize(Nn*3,Nn*3);
    M_M_tilde.resize(Nn*3,Nn*3);
    M_D_tilde.resize(Nn*3,Nn*3);
    
    M_K1.resize(Nn*3,Nn*2);
    M_D1.resize(Nn*3,Nn*2);
    
    M_K2.setFromTriplets(T_K2.begin(), T_K2.end());
    M_M .setFromTriplets(T_M .begin(), T_M .end());
    M_D .setFromTriplets(T_D .begin(), T_D .end());
    M_D0.setFromTriplets(T_D0.begin(), T_D0.end());
    M_K .setFromTriplets(T_K .begin(), T_K .end());
    M_H .setFromTriplets(T_H .begin(), T_H .end());
    
    M_H1_tilde.setFromTriplets(T_H1_tilde.begin(), T_H1_tilde.end());
    M_K1_tilde.setFromTriplets(T_K1_tilde.begin(), T_K1_tilde.end());
    
    M_K_tilde.setFromTriplets(T_K_tilde .begin(), T_K_tilde.end());
    M_M_tilde.setFromTriplets(T_M_tilde .begin(), T_M_tilde.end());
    M_D_tilde.setFromTriplets(T_D_tilde.begin(), T_D_tilde.end());
    
    M_K1.setFromTriplets(T_K1.begin(), T_K1.end());
    M_D1.setFromTriplets(T_D1.begin(), T_D1.end());
    
    return 0;
}

void PML_solver::sub_fixed_dof()
{
    long Nn = Table_Noeuds.cols();
    
    num_ddl_u = 0;
    num_ddl_E = 0;
    
    vector<T> T_IDDL_U;
    T_IDDL_U.reserve(2*Nn);
    
    vector<T> T_IDDL_E;
    T_IDDL_E.reserve(3*Nn);
    
    for(int k =0; k< Nn; k++)
    {
        if(Table_BC(k) == 0 )
        {
            T_IDDL_U.push_back(T(2*k  ,num_ddl_u++  ,1.f));
            T_IDDL_U.push_back(T(2*k+1,num_ddl_u++,1.f));
        }
        
        if(Table_Noeuds(1,k) <= settings.ybot_pml|| abs(Table_Noeuds(0,k)-x_mid) >= settings.x_pml )
        {
            if (abs(Table_Noeuds(0,k)-x_mid) >= settings.x_pml )
                T_IDDL_E.push_back(T(3*k,num_ddl_E++,1.f));
            
            if (Table_Noeuds(1,k) <= settings.ybot_pml )
                T_IDDL_E.push_back(T(3*k+1,num_ddl_E++,1.f));
            
            T_IDDL_E.push_back(T(3*k+2,num_ddl_E++,1.f));
        }
        
    }
    
    M_DDL_U.resize(2*Nn ,num_ddl_u);
    M_DDL_U.setFromTriplets(T_IDDL_U.begin(), T_IDDL_U.end());
    
    SP_mat M_DDL_E(3*Nn ,num_ddl_E);
    M_DDL_E.setFromTriplets(T_IDDL_E.begin(), T_IDDL_E.end());
    
    M_K2 = (M_DDL_U.transpose() * M_K2 * M_DDL_U).eval();
    M_M  = (M_DDL_U.transpose() * M_M  * M_DDL_U).eval();
    M_D  = (M_DDL_U.transpose() * M_D  * M_DDL_U).eval();
    M_K  = (M_DDL_U.transpose() * M_K  * M_DDL_U).eval();
    M_D0 = (M_DDL_U.transpose() * M_D0 * M_DDL_U).eval();
    M_H  = (M_DDL_U.transpose() * M_H  * M_DDL_U).eval();
    
    M_H1_tilde = (M_DDL_U.transpose() * M_H1_tilde * M_DDL_E ).eval();
    M_K1_tilde = (M_DDL_U.transpose() * M_K1_tilde * M_DDL_E ).eval();
    M_K_tilde = (M_DDL_E.transpose() * M_K_tilde * M_DDL_E ).eval();
    M_M_tilde = (M_DDL_E.transpose() * M_M_tilde * M_DDL_E ).eval();
    M_D_tilde = (M_DDL_E.transpose() * M_D_tilde * M_DDL_E ).eval();
    
    
    M_K1 = (M_DDL_E.transpose() * M_K1 * M_DDL_U ).eval();
    M_D1 = (M_DDL_E.transpose() * M_D1 * M_DDL_U ).eval();
    
    
    V_F = (M_DDL_U.transpose() * V_F).eval();

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
    Force_time = VectorXd::Zero(settings.Niter);
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
    long Nn = Table_Noeuds.cols();
    ddl_force=0;
    {
        double dist = 1000 ;
        
        for (unsigned k=0; k<Nn; k++)
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

    V_F = VectorXd::Zero(2*Nn) ;
    V_F(ddl_force) = 1.F;
}

void PML_solver::sub_construction_global_Indice()
{
    vector<T> T_Ind1, T_Ind2;
    
    T_Ind1.reserve(num_ddl_u);
    for (int j=0; j< num_ddl_u; j++)
        T_Ind1.push_back(T(j,j,1.f));
    
    T_Ind2.reserve(num_ddl_E);
    for (int j=0; j< num_ddl_E; j++)
        T_Ind2.push_back(T(j,num_ddl_u+j,1.f));
    
    
    M_Ind1.resize(num_ddl_u,num_ddl_u+num_ddl_E);
    M_Ind1.setFromTriplets(T_Ind1.begin(), T_Ind1.end());
    
    
    M_Ind2.resize(num_ddl_E,num_ddl_u+num_ddl_E);
    M_Ind2.setFromTriplets(T_Ind2.begin(), T_Ind2.end());
}

void PML_solver::set_matrix_system()
{
    int Icase = settings.Icase;
    
    if (Icase == 4 || Icase == 5 || Icase == 6)
        set_Newmark();
    
    if (Icase == 1 || Icase == 3)
    {
        M_K_dyn = dcomplex(1,0) * ( M_K -(omega * omega) * M_M ) + dcomplex(0,1) * omega * M_D0 ;
        
        if (Icase == 3)
            M_K_dyn += dcomplex(0,1) * omega * M_D ;
        
        V_F_pml_complex = dcomplex(1,0) * V_F;
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
        
        M_Mglob = M_Ind1.transpose() * M_M * M_Ind1
                                                  + M_Ind2.transpose()  * M_M_tilde * M_Ind2 ;
        
        M_Dglob = M_Ind1.transpose()  * (M_D + M_D0) * M_Ind1
                + M_Ind2.transpose()  * M_D1         * M_Ind1 + M_Ind2.transpose() * M_D_tilde * M_Ind2 ;
        
        M_Kglob = M_Ind1.transpose()  * (M_K+M_K2) * M_Ind1 + M_Ind1.transpose() * M_K1_tilde * M_Ind2
                + M_Ind2.transpose()  * M_K1       * M_Ind1 + M_Ind2.transpose() * M_K_tilde  * M_Ind2 ;
        
        M_Hglob = M_Ind1.transpose()  * M_H * M_Ind1 + M_Ind1.transpose()  * M_H1_tilde * M_Ind2;
        
        if (Icase == 2)
        {
            M_K_dyn =   dcomplex(1,0) * ( M_Kglob - (omega * omega) * M_Mglob )
            + dcomplex(0,1) * ( omega * M_Dglob  -1.f/omega * M_Hglob );
            
            V_F_pml_complex =  M_Ind1.transpose() * V_F ;
        }
        else
        {
            M_K_eff = Alpha_1 * M_Mglob + Alpha_4 * M_Dglob  + M_Kglob + 0.5 * dt *  M_Hglob;
            V_F_pml_real = M_Ind1.transpose() * V_F;
        }
        
        
    }
}

unsigned PML_solver::get_Icase()
{
    return settings.Icase;
}

double PML_solver::get_final_time()
{
    return settings.Niter * dt ;
}

int PML_solver::solve()
{
    int Icase = settings.Icase ;
    
    if (Icase == 1 || Icase == 2 || Icase == 3)
    {
        M_K_dyn.makeCompressed();
        VectorXcd  U_pml;
        umfpack_wrapper solver(M_K_dyn);
        U_pml = solver.solve(V_F_pml_complex);
        
        if (Icase == 2)
            Ustock_pml_complex = M_DDL_U * M_Ind1 * U_pml;
        
        else if (Icase == 1|| Icase == 3)
            Ustock_pml_complex = M_DDL_U * U_pml;
        
    }
    
    if (Icase == 4 || Icase == 5 || Icase == 6)
    {
        
        long Nn = Table_Noeuds.cols();
        Ustock_pml_real = MatrixXd::Zero( 2*Nn , settings.Niter);
        
        long NDDL = M_K_eff.rows();
        VectorXd Feff, U_pml = VectorXd::Zero(NDDL), V_pml = VectorXd::Zero(NDDL), A_pml = VectorXd::Zero(NDDL) ;
        
        double fnew = 0.f, fold = 0.f ;
        
        M_K_eff.makeCompressed();
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
            
            VectorXd DU = solver.solve(Feff);
            
            VectorXd DV = Alpha_4*DU + Alpha_5*V_pml + Alpha_6*A_pml;
            VectorXd DA = Alpha_1*DU + Alpha_2*V_pml + Alpha_3*A_pml;
            
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
    long Ne = Table_Elements.rows();
    long Nn = Table_Noeuds.cols();
    
    int  Icase = settings.Icase ;
    
    if (Icase == 4)
    {
        ofstream output;
        string filename("RESULTS_ref_transient.bin");
        
        chrono_message("Ecriture du fichier " +filename +"\n");
        
        output.open (filename,ios::out | ios::binary| ios::trunc);
        output.write( (char*)Ustock_pml_real.data(),sizeof(double) * Ustock_pml_real.size() );
        output.close();
    }
    
    MatrixXd Ustock_ref_real;
    
    if (Icase == 5 || Icase == 6)
    {
        long Nb = Ustock_pml_real.size() ;
        Ustock_ref_real = MatrixXd::Zero(Ustock_pml_real.rows(),Ustock_pml_real.cols());
        string filename("RESULTS_ref_transient.bin");
        std::ifstream input( filename, std::ios::binary  );
        input.read((char*)Ustock_ref_real.data(),Nb*sizeof(double));
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
                Ustock_ref_real = (Ustock_pml_real - Ustock_ref_real).array().abs();
                double normalization = Ustock_ref_real.maxCoeff();
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
    long Ne = Table_Elements.rows();
    long Nn = Table_Noeuds.cols();
    
    int  Icase = settings.Icase ;
    
    // TO DO
}






