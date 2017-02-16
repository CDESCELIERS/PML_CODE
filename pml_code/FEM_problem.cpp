//
//  FEM_problem.cpp
//

#include "FEM_problem.hpp"
#include "tools.hpp"
#include <libconfig.h++>
#include <iostream>

using namespace libconfig;
using namespace FEM;

/* -------------------------
 
 CLASS FEM_problem
 
 --------------------------- */

int FEM_problem::read_cfg(string CONFIG_FILE)
{
    Config cfg;
    struct config_settings settings;
    
    try
    {
        cfg.readFile(CONFIG_FILE.c_str());
    }
    catch(const FileIOException &fioex)
    {
        cout << endl << "I/O error while reading configuration file" << endl << endl;
        return -1 ;
    }
    catch(const ParseException &fioex)
    {
        cout << endl << "Parsing error while reading configuration file" << endl << endl;
        return -1 ;
    }
    
    cout << endl ;
    try {
        unsigned val = (unsigned) cfg.lookup("Icase");
        settings.Icase = val ;
    }
    catch(const SettingNotFoundException &nfex) {
        cout << endl <<"Value of 'Icase' is missing" << endl << endl;
        return -1;
    }
    catch(const SettingTypeException &nfex) {
        cout << endl <<"Wrong type for 'Icase'" << endl << endl;
        return -1;
    }
    
    try {
        string val  = cfg.lookup("vtk");
        settings.FILEOUT = val ;
    }
    catch(const SettingNotFoundException &nfex) {
        cout << endl <<"Value of 'vtk' is missing" << endl << endl;
        return -1;
    }
    catch(const SettingTypeException &nfex) {
        cout << endl <<"Wrong type for 'vtk'" << endl<< endl;
        return -1;
    }
    
    try {
        double val  = (double)cfg.lookup("max_mu");
        settings.max_mu_PML = val ;
    }
    catch(const SettingNotFoundException &nfex) {
        cout << endl <<"Value of 'max_mu' is missing" << endl << endl;
        return -1;
    }
    catch(const SettingTypeException &nfex) {
        cout << endl <<"Wrong type for 'max_mu'" << endl << endl;
        return -1;
    }
    
    try {
        double val  = (double)cfg.lookup("a_PML");
        settings.a_PML = val ;
    }
    catch(const SettingNotFoundException &nfex) {
        cout << endl <<"Value of 'a_PML' is missing" << endl << endl;
        return -1;
    }
    catch(const SettingTypeException &nfex) {
        cout << endl <<"Wrong type for 'a_PML'" << endl << endl;
        return -1;
    }
    
    try {
        double val  = (double)cfg.lookup("a_force");
        settings.a_force = val ;
    }
    catch(const SettingNotFoundException &nfex) {
        cout << endl <<"Value of 'a_force' is missing" << endl << endl;
        return -1;
    }
    catch(const SettingTypeException &nfex) {
        cout << endl <<"Wrong type for 'a_force'" << endl << endl;
        return -1;
    }
    
    try {
        double val  = (double)cfg.lookup("x_pml");
        settings.x_pml = val ;
    }
    catch(const SettingNotFoundException &nfex) {
        cout << endl <<"Value of 'x_pml' is missing" << endl << endl;
        return -1;
    }
    catch(const SettingTypeException &nfex) {
        cout << endl <<"Wrong type for 'x_pml' is" << endl << endl;
        return -1;
    }
    
    try {
        double val  = (double)cfg.lookup("ybot_pml");
        settings.ybot_pml = val ;
    }
    catch(const SettingNotFoundException &nfex) {
        cout << endl <<"Value of 'ybot_pml' is missing" << endl << endl;
        return -1;
    }
    catch(const SettingTypeException &nfex) {
        cout << endl <<"Wrong type for 'ybot_pml'" << endl << endl;
        return -1;
    }
    
    try {
        unsigned val  = (unsigned)cfg.lookup("Niter");
        settings.Niter = val ;
    }
    catch(const SettingNotFoundException &nfex) {
        cout << endl <<"Value of 'Niter' is missing" << endl << endl;
        return -1;
    }
    catch(const SettingTypeException &nfex) {
        cout << endl <<"Wrong type for 'Niter'" << endl << endl;
        return -1;
    }
    
    try {
        double val  = (double)cfg.lookup("speed");
        settings.speed = val ;
    }
    catch(const SettingNotFoundException &nfex) {
        cout << endl <<"Value of 'speed' is missing" << endl << endl;
        return -1;
    }
    catch(const SettingTypeException &nfex) {
        cout << endl << "Wrong type for 'speed'" << endl << endl;
        return -1;
    }
    
    try {
        string val = cfg.lookup("msh");
        settings.msh_file = val ;
    }
    catch(const SettingNotFoundException &nfex) {
        cout << endl <<"Value of 'mshfile' is missing" << endl << endl;
        return -1;
    }
    catch(const SettingTypeException &nfex) {
        cout << endl <<"Wrong type for 'mshfile'" << endl << endl;
        return -1;
    }
    
    try {
        int val = cfg.lookup("problem_type");
        problem_type = val ;
        settings.problem_type = val;
    }
    catch(const SettingNotFoundException &nfex) {
        cout << endl <<"Value of 'mshfile' is problem_type" << endl << endl;
        return -1;
    }
    catch(const SettingTypeException &nfex) {
        cout << endl <<"Wrong type for 'problem_type'" << endl << endl;
        return -1;
    }
    
    try {
        int val = cfg.lookup("Nsamples");
        Nsamples = (long)val ;
        settings.Nsamples = (long)val;
    }
    catch(const SettingNotFoundException &nfex) {
        cout << endl <<"Value of 'mshfile' is Nsamples" << endl << endl;
        return -1;
    }
    catch(const SettingTypeException &nfex) {
        cout << endl <<"Wrong type for 'Nsamples'" << endl << endl;
        return -1;
    }
    
    try {
        int val = cfg.lookup("DEBUG");
        settings.debug = val;
    }
    catch(const SettingNotFoundException &nfex) {
        cout << endl <<"Value of 'mshfile' is Nsamples" << endl << endl;
        return -1;
    }
    catch(const SettingTypeException &nfex) {
        cout << endl <<"Wrong type for 'Nsamples'" << endl << endl;
        return -1;
    }

    
    const Setting& root = cfg.getRoot();
    
    try
    {
        const Setting &materials = root["materials"];
        int count = materials.getLength();
        for(int i = 0; i < count; ++i)
        {
            const Setting &material = materials[i];
            double E=0.0, v=0.0, rho=0.0, Xi=0.0;
            int case_mat=0;
            if(!(material.lookupValue("E", E)
                 && material.lookupValue("v", v)
                 && material.lookupValue("case_mat", case_mat)
                 && material.lookupValue("rho", rho)
                 && material.lookupValue("Xi", Xi))) {
                cout << endl << "Wrong settings for material " << i+1 << endl << endl;
                return -1;
            }
            
            vec parameters =
            {
                E, v, rho, Xi
            };
            
            settings.Mat[case_mat] = parameters ;
            
        }
    }
    catch(const SettingNotFoundException &nfex)
    {
        cout << endl << "No materials seetings" << endl << endl;
        return -1;
    }
    
    try
    {
        const Setting &bcs = root["BC"];
        int count = bcs.getLength();
        for(int i = 0; i < count; ++i)
        {
            const Setting &bc = bcs[i];
            int type=0, num=0;
            if(!(bc.lookupValue("type", type)
                 && bc.lookupValue("num", num))) {
                cout << endl << "Wrong settings for boundary conditions " << i+1 << endl << endl;
                return -1;
            }
            
            if (type == 1) {
                // TO DO IF WE EVER NEED TO BLOCK ONLY X
            }
            else if (type == 2){
                // TO DO IF WE EVER NEED TO BLOCK ONLY Y
            }
            else if (type == 3)
                settings.BC_list.insert(num);
            
        }
    }
    catch(const SettingNotFoundException &nfex)
    {
        cout << endl << "No materials seetings" << endl << endl;
        return -1;
    }
    
    problem.set_settings(settings);
    
    return 0;
}

int FEM_problem::solve()
{
    if (problem_type == DETERMINIST)
        return solve_determinist();
    
    if (problem_type == RANDOM)
        return solve_random();
    
    return 0;
}

int FEM_problem::solve_determinist()
{
    start_chrono();
    cout <<  "(0.0 secs) Lecture du maillage " << problem.get_mshfile() << endl;
    
    if (problem.read_mesh() < 0) return -1;
    
    chrono_message("Construction de la force\n");
    
    problem.set_external_force_parameters();
    
    problem.set_unit_force();
    
    chrono_message("Construction des matrices élémentaires et assemblage\n");
    
    if (problem.sub_element_T6() < 0 ) return -1;
    
    chrono_message("Elimination des ddl nuls\n");
    
    problem.sub_fixed_dof();
    
    chrono_message("Construction du système matriciel\n");
    
    problem.set_matrix_system();
    
    chrono_message("Resolution");
    
    unsigned Icase = problem.get_Icase();
    
    if (Icase == 4 || Icase == 5 || Icase == 6)
        cout << " (T_f = "<< problem.get_final_time() << " secs)" << endl << endl;
    
    if (problem.solve() < 0 ) return -1;
    
    problem.save_vtk();
    
    chrono_message("Fin de la simulation \n\n");

    return 0;
}

int FEM_problem::solve_random()
{
    start_chrono();
    cout << endl << endl <<  "(0.0 secs) Lecture du maillage " << problem.get_mshfile() << endl;
    
    if (problem.read_mesh() < 0) return -1;
    
    chrono_message("Construction de la force\n");
    
    problem.set_external_force_parameters();
    
    problem.set_unit_force();
    
    chrono_message("Début de la simulation de MonteCarlo\n");
    
    for (int k=0; k<Nsamples; k++) 
    {
        printProgress( (double)k/(double)Nsamples);

        problem.sample_random_settings();
        
        if (problem.sub_element_T6() < 0 ) return -1;
        
        problem.sub_fixed_dof();
    
        problem.set_matrix_system();
    
        problem.solve() ;
    }
    
    problem.save_vtk();
    
    chrono_message("Fin de la simulation \n\n");

    return 0;
}
