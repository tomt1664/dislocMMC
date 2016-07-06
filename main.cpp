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
*****************************************************/

/*

Objects: configurations, atoms, rings
set-ups, lammps parameters

*/

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include "bond.h"
#include "config.h"
#include "atom.h"
#include "polygon.h"

int main()
{
    std::cout << "dislocMMC v.0.48" << std::endl;

    //initial configuration
    std::string inconfig = "data.in";

    //output files
    std::string polyfile = "poly.xyz";
    std::string outfile = "traj.xyz";
    std::string bondfile = "bonds.xyz";

    std::vector<Polygon> polys;
    std::vector<Bond> bonds;

    Config config(inconfig);

//    long natoms = config.nat();

    double cutof = 1.8;

 //   long nbnd = config.analyse(cutof,bonds);


    int berror = config.getbondlist(cutof);
    long nbnd = config.getn2();
    std::cout << "numbnd " << nbnd << std::endl;
    nbnd = config.getn3();
    std::cout << "numn3 " << nbnd << std::endl;
    nbnd = config.getn4();
    std::cout << "numn4 " << nbnd << std::endl;
    nbnd = config.getn5();
    std::cout << "numn5 " << nbnd << std::endl;
    nbnd = config.getPents(polys);
    std::cout << "numn pents " << nbnd << " " << polys.size() << std::endl;
    polys.clear();
    nbnd = config.getHepts(polys);
    std::cout << "numn hepts " << nbnd << " " << polys.size() << std::endl;
    nbnd = config.getBonds(bonds);

    std::cout << "numn bonds " << nbnd << std::endl;


    config.write(outfile);
    config.writePoly(polyfile);
    config.writeBonds(bondfile);

    double en = config.relax();
    std::cout << " en " << en << std::endl;

    nbnd = config.update(cutof,1);

    std::cout << "nbnd update " << nbnd << std::endl;

    return 0;
}
