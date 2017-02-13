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


