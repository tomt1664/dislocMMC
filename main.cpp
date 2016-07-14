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

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include "bond.h"
#include "config.h"
#include "atom.h"
#include "polygon.h"
#include "setup.h"

int main()
{
    std::cout << "dislocMMC v.0.60" << std::endl;

    //initial configuration
    std::string inconfig = "data.in";

    //options input file
    std::string inoptions = "setup.in";

    //name output files
    std::string polyfile = "poly.xyz"; // 5-fold and 7-fold polygon coordinates: xyz movie
    std::string outfile = "traj.xyz"; // full system coordinates: xyz movie
    std::string bondfile = "bonds.xyz"; // active bond site coordinates: xyz movie
    std::string enfile = "energy.out"; // optimised system energy per step

    // polygon and active bond storage vectors
    std::vector<Polygon> polys;
    std::vector<Bond> bonds;

    // get the simulation set-up and options from input file
    Setup setup(inoptions);

    // construct initial configuration from input file
    Config config(inconfig);

    // do structure analysis
    std::cout << "Initial structure analysis ..." << std::endl;
    long numbnd = config.analyse(setup.coff(),bonds);
    std::cout << "Found " << numbnd << " active bonds" << std::endl;


    //write initial configuration and features to output
    config.write(outfile);
    config.writePoly(polyfile);
    config.writeBonds(bondfile);

    //perform initial relaxation
    std::cout << "First optimisation ..." << std::endl;
    double en = config.relax(setup.lmp());
    std::cout << "Complete. Final energy: " << en << " eV" << std::endl;

    //rotate bond (2)
    std::cout << "Rotate bond" << std::endl;
    config.rotate(1);
    std::cout << "Optimise ..." << std::endl;
    en = config.relax(setup.lmp());
    std::cout << "Complete. Final energy: " << en << " eV" << std::endl;

    std::cout << "Update structure analysis ..." << std::endl;
    numbnd = config.update(setup.coff(),1);
    std::cout << "Found " << numbnd << " active bonds" << std::endl;

    //write new configuration and features to output
    config.write(outfile);
    config.writePoly(polyfile);
    config.writeBonds(bondfile);

    return 0;
}
