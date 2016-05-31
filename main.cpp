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
#include "config.h"
#include "atom.h"
#include "optim.h"

using namespace std;

int main()
{
    cout << "dislocMMC v.0.22" << endl;

    Config config("data.in");


    return 0;
}
