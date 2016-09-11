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
#include <random>
#include "bond.h"
#include "config.h"
#include "atom.h"
#include "polygon.h"
#include "setup.h"

int main()
{
    std::cout << "dislocMMC v.0.60" << std::endl;

    //initial configuration file
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
    std::vector<Bond> bondc;

    // get the simulation set-up and options and LAMMPS command from the input file
    Setup setup(inoptions);

    //initialise Mersenne Twister with seed
    std::mt19937 generator (setup.seed());

    // construct initial configuration from input file
    Config config(inconfig);

    // perform initial structure analysis
    std::cout << "Initial structure analysis ..." << std::endl;
    long numbnd = config.analyse(setup.coff(),bonds,bondc);
    std::cout << "Found " << numbnd << " active bonds" << std::endl;

    //write initial configuration and features (rotatable bonds and polygons) to output
    config.write(outfile);
    config.writePoly(polyfile);
    config.writeBonds(bondfile);

    //perform initial relaxation
    std::cout << "First optimisation ..." << std::endl;
    double en = config.relax(setup.lmp());
    std::cout << "Complete. Final energy: " << en << " eV" << std::endl;

    //do downhill energy minimisation
    int ibnd = -1;
    if(setup.type() == 1) {
        for(int istep = 0; istep < setup.maxs(); istep++) {
            //record each bond rotated configuration
            std::vector<Config> configR;
            configR.clear();
            //loop over all rotatable atom pairs - optimise and find lowest energy
            double minen = 0.0;
            ibnd = -1;
            for(int nb = 0; nb < numbnd; nb++) {
                //copy configuration
                Config tconfig = config;
                //rotate bond
                tconfig.rotate(nb);
                std::cout << "Optimise configuration " << nb << " ..." << std::endl;
                en = tconfig.relax(setup.lmp());
                std::cout << "Complete. Final energy: " << en << " eV" << std::endl;
                configR.push_back(tconfig);
                //get the lowest energy structure
                if(en < minen) {
                    minen = en;
                    ibnd = nb;
                }
            }
            std::cout << "Selected rotation: " << ibnd << std::endl;
            //copy lowest energy rotation into config
            config = configR[ibnd];

            numbnd = config.update(setup.coff(),ibnd);
            std::cout << "Update. Bonds: " << numbnd << std::endl;

            //write new configuration and features to output
            config.write(outfile);
//            config.writePoly(polyfile);
            config.writeBonds(bondfile);
            config.writeEn(enfile);

            //update bond table
        }
    } else if(setup.type() == 2)
    //do Metropolis-Hastings Monte Carlo
    {
        for(int istep = 0; istep < setup.maxs(); istep++) {
            //pick a rotatable bond at random
            double rnum = generator()*1.0/(generator.max()*1.0);
            int nb = 0;
            for(int ib = 0; ib < numbnd; ib++) {
                double rscale = rnum*numbnd;
                if(rscale > ib*1.0 && rscale <= (ib+1)*1.0) {
                    std::cout << "r " << rnum << " rs " << rscale << " ib " << ib << std::endl;
                    nb = ib;
                    break;
                }
            }
            //copy configuration
            Config tconfig = config;
            //rotate bond
            tconfig.rotate(nb);
            std::cout << "Optimise configuration " << nb << " ..." << std::endl;
            double en2 = tconfig.relax(setup.lmp());
            std::cout << "Complete. Final energy: " << en << " eV" << std::endl;

            //if energy is lower, accept move
            if(en2 <= en) {
                std::cout << "Energy change: " << (en2 - en) << " eV. Move accepted." << std::endl;
                config = tconfig;
                numbnd = config.update(setup.coff(),nb);
                std::cout << "Update. Bonds: " << numbnd << std::endl;
                //write new configuration and features to output
                config.write(outfile);
//                config.writePoly(polyfile);
                config.writeBonds(bondfile);
                config.writeEn(enfile);
            } else {
            //if energy is higher, calculate the probability from the Boltzman distribution
                double sienergy = (en2 - en)*1.602e-19;
                double boltz = exp(-sienergy/(1.38e-23*setup.tmp()));
                std::cout << "Energy change: " << (en2 - en) << " eV. Boltz Prob: " << boltz << std::endl;
                //get new random  number and determine move
                double rnum = generator()*1.0/(generator.max()*1.0);
                if(boltz > rnum) {
                    std::cout << "Move Accepted." << std::endl;
                    config = tconfig;
                    numbnd = config.update(setup.coff(),nb);
                    std::cout << "Update. Bonds: " << numbnd << std::endl;
                    config.write(outfile);
//                    config.writePoly(polyfile);
                    config.writeBonds(bondfile);
                    config.writeEn(enfile);
                } else {
                    std::cout << "Move Rejected." << std::endl;
                }
            }
        }

    } else if(setup.type() == 3)
    //do double climb
    {
        config.dclimb(0);
        config.relax(setup.lmp());
        config.write(outfile);
    } else if(setup.type() == 4)
    //do single climb
    {
        config.sclimb(2);
        config.relax(setup.lmp());
        config.write(outfile);
        config.analyse(setup.coff(),bonds,bondc);
        std::vector<Polygon> octs;
        long n8 = config.getOcts(octs);
        std::cout << " 8s: " << n8 << std::endl;

        config.shuffle(0,-1);
        config.write(outfile);
        config.relax(setup.lmp());
        config.write(outfile);

        config.shclimb(0);
        config.relax(setup.lmp());
        config.write(outfile);

    } else if(setup.type() == 5)
    //do single climb on shuffle
    {
        config.shclimb(0);
        config.relax(setup.lmp());
        config.write(outfile);
        config.analyse(setup.coff(),bonds,bondc);
        std::vector<Polygon> octs;
        long n8 = config.getOcts(octs);
        std::cout << " 8s: " << n8 << std::endl;
    } else {
        std::cout << "Invalid type option" << std::endl;
        exit(1);
    }

    std::cout << "Complete. Exiting ..." << std::endl;

    return 0;
}
