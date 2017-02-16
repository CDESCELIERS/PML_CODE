//
//  FEM_solver.hpp
//

#ifndef PML_SOLVER_hpp
#define PML_SOLVER_hpp

#include <armadillo>
#include <map>
#include <set>

using namespace std;
using namespace arma;

/* -------------------------
 
 CLASS PML_solver
 
 --------------------------- */

namespace FEM
{
    static const int RANDOM = 1;
    static const int DETERMINIST = 0;
    static const int DEBUG_PRINT_GLOBAL_MATRICES = 1;
    
    typedef std::map<uword, arma::vec> Materials;
    
    void save_sp2ascii(sp_mat &MAT, string filename);
    void save_sp2ascii(sp_cx_mat &MAT, string filename);
    void save_vec2ascii(cx_vec &VEC, string filename);
    void save_vec2ascii(vec &VEC, string filename);
    void save_mat2ascii(umat &MAT, string filename);
    void save_mat2ascii(mat &MAT, string filename);
    void save_int2ascii(sword value, string filename);

    
    struct config_settings {
        uword Icase ;
        string FILEOUT ;
        double a_force ;
        uword Niter ;
        double speed  ;
        Materials Mat;
        double ybot_pml;
        double a_PML;
        double max_mu_PML;
        string msh_file;
        set<uword> BC_list;
        double x_pml;
        int problem_type;
        long Nsamples;
        uword debug;
    } ;
    
    class PML_solver {
        
    private:
        
        config_settings settings;
        
        double x_mid;
        double ytop;
        
        sp_mat M_K;
        sp_mat M_M;
        sp_mat M_D;
        sp_mat M_D0;
        sp_mat M_H;
        sp_mat M_H1_tilde;
        sp_mat M_K1_tilde;
        sp_mat M_K_tilde;
        sp_mat M_M_tilde;
        sp_mat M_D_tilde;
        sp_mat M_K2;
        sp_mat M_D1;
        sp_mat M_K1;
        
        vec V_F;
        int      num_ddl_u;
        int      num_ddl_E;
        sp_mat   M_DDL_U;
        sp_mat   M_DDL_E;
        sp_mat   M_Ind1;
        sp_mat   M_Ind2;
        
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
        vec Force_time ;
        
        uword ddl_force;
        
        sp_mat M_Mglob;
        sp_mat M_Dglob;
        sp_mat M_Kglob;
        sp_mat M_Hglob ;
        
        sp_cx_mat M_K_dyn ;
        sp_mat M_K_eff ;
        
        cx_vec V_F_pml_complex;
        vec V_F_pml_real;
          
        cx_vec Ustock_pml_complex;
        mat Ustock_pml_real;
        
        umat  Table_Elements;
        mat   Table_Noeuds;
        umat  Table_BC;
        
        void set_xmid_ytop();
        void sub_construction_global_Indice() ;
        void sub_element_orientation_check();
        void set_Newmark();
        
        void save_vtk_determinist();
        void save_vtk_random();   // TO DO
        
    public:
        
        void set_settings(config_settings settings);
        int read_mesh();
        int  sub_element_T6();
        void sub_fixed_dof() ;
        void set_external_force_parameters();
        void set_frequency(double w);
        void set_unit_force();
        void set_matrix_system();
        uword get_Icase();
        double get_final_time();
        int solve();
        void save_vtk();
        string get_mshfile();
        
        void sample_random_settings(); // TO DO
    };
    
    
    
}
#endif 
