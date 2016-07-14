/****************************************************
**  dislocMMC
**
**  A Metropolis Monte Carlo algorithm to model the
**  evolution of structure in defective sp2 bonded
**  carbon systems
**
**  The code calls the LAMMPS atomistic simulation
**  program for structural optimisation with a
**  reactive force-field
**
**  T.Trevethan 2016
*****************************************************/

#include "setup.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>

Setup::Setup(std::string& setupfile)
{
    //set default values
    m_cutoff = 1.80;
    m_pstep = 1;
    m_maxs = 9999;
    m_type = 1;
    m_temp = 300.0;
    m_etol = 0.001;

    //open data input file for inital structure
    std::ifstream insetup(setupfile.c_str());
    if (!insetup) {
        std::cerr << "Error: set-up file could not be opened" << std::endl;
        exit(1);
    }

    //read from file (enhanced xyz file format)
    std::string data;
    std::stringstream ss; //use stringstream

    std::getline(insetup,data);
    m_lmpcommand = data;
    std::getline(insetup,data);
    ss >> m_type >> m_maxs >> m_pstep >> m_cutoff >> m_temp >> m_etol;

}

Setup::~Setup()
{
    //dtor
}
