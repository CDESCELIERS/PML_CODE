//
//  FEM_solver.hpp
//

#ifndef PML_SOLVER_hpp
#define PML_SOLVER_hpp

#include "FEM_type.hpp"
#include <set>

using namespace std;
using namespace Eigen;

/* -------------------------
 
 CLASS PML_solver
 
 --------------------------- */

namespace FEM
{
    static const int DEBUG_PRINT_GLOBAL_MATRICES = 1;
    static const int RANDOM = 1;
    static const int DETERMINIST = 0;
    
    
    void save_sp2ascii(SP_mat_complex &MAT, string filename);
    void save_sp2ascii(SP_mat &MAT, string filename);
    void save_vec2ascii(VectorXcd &VEC, string filename);
    void save_vec2ascii(VectorXd  &VEC, string filename);
    void save_mat2ascii(Matrix<unsigned,Dynamic, 7> &MAT, string filename);
    void save_mat2ascii(Matrix<double,2, Dynamic> &MAT, string filename);
    void save_int2ascii(int value, string filename);
    
    
    struct config_settings {
        unsigned Icase ;
        string FILEOUT ;
        double a_force ;
        unsigned Niter ;
        double speed  ;
        Materials Mat;
        double ybot_pml;
        double a_PML;
        double max_mu_PML;
        string msh_file;
        set<unsigned> BC_list;
        double x_pml;
        int problem_type;
        long Nsamples;
        int debug;
    } ;
    
    class PML_solver {
        
    private:
        
        config_settings settings;
        
        int debug;
        
        double x_mid;
        double ytop;
        
        SP_mat M_K;
        SP_mat M_M;
        SP_mat M_D;
        SP_mat M_D0;
        SP_mat M_H;
        SP_mat M_H1_tilde;
        SP_mat M_K1_tilde;
        SP_mat M_K_tilde;
        SP_mat M_M_tilde;
        SP_mat M_D_tilde;
        SP_mat M_K2;
        SP_mat M_D1;
        SP_mat M_K1;
        
        VectorXd V_F;
        int      num_ddl_u;
        int      num_ddl_E;
        SP_mat   M_DDL_U;
        SP_mat   M_DDL_E;
        SP_mat   M_Ind1;
        SP_mat   M_Ind2;
        
        double f_1;
        double omega;
        double dt;
        double gamma;
        double beta;
        double Alpha_1;
        double Alpha_2;
        double Alpha_3;
        double Alpha_4;
        double Alpha_5;
        double Alpha_6;
        VectorXd Force_time ;
        
        unsigned ddl_force;
        
        SP_mat M_Mglob;
        SP_mat M_Dglob;
        SP_mat M_Kglob;
        SP_mat M_Hglob ;
        
        SP_mat_complex M_K_dyn ;
        SP_mat M_K_eff ;
        
        Matrix<dcomplex, Dynamic, 1> V_F_pml_complex;
        VectorXd V_F_pml_real;
          
        Matrix<dcomplex, Dynamic, 1> Ustock_pml_complex;
        MatrixXd Ustock_pml_real;
        
        Matrix<unsigned,Dynamic,7>   Table_Elements;
        Matrix<double,2,Dynamic>     Table_Noeuds;
        Matrix<unsigned,Dynamic, 1>  Table_BC;
        
        void set_xmid_ytop();
        void sub_construction_global_Indice() ;
        void sub_element_orientation_check();
        void set_Newmark();
        
        void save_vtk_determinist();
        void save_vtk_random();   // TO DO
        
    public:
        
        PML_solver();
        
        void set_settings(config_settings settings);
        int read_mesh();
        int  sub_element_T6();
        void sub_fixed_dof() ;
        void set_external_force_parameters();
        void set_frequency(double w);
        void set_unit_force();
        void set_matrix_system();
        unsigned get_Icase();
        double get_final_time();
        int solve();
        void save_vtk();
        string get_mshfile();
        
        void sample_random_settings(); // TO DO
    };
    
    
    
}
#endif 
