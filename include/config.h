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

//class for the atomic configuration of the system and boundary conditions

#ifndef CONFIG_H
#define CONFIG_H

#include <vector>
#include <string>

#include "atom.h"
#include "bond.h"
#include "polygon.h"

class Config
{
    public:
        Config(const std::string& infile); //construct config object with input file
        Config(std::vector<Atom>& inconf, double nbounds[5]); //constrct config object with atom vector
        virtual ~Config();
        int getbondlist(double cutoff); // get the bond list
        long getn2(); //find all bonds
        long getn3(); //find all 3 sequentially bonded atoms
        long getn4(); // find all 4 sequentially bonded atoms
        long getn5(); // find all 5 sequentially bonded atoms

        long getPents(std::vector<Polygon>& pents); //get a list of all the pentagons in the system (return number)
        long getHepts(std::vector<Polygon>& hepts); //get a list of all the heptagons in the system (return number)
        long getBonds(std::vector<Bond>& bonds); //get a list of all the bond rotation pairs (return number)

        long analyse(double cutoff, std::vector<Bond>& bondl);
        long update(double cutoff, long ibnd);

        int rotate(long ibnd); //rotate the bond pair of index ibnd (return success)
        long getCoords(std::vector<Atom>& coords); //get the atomic coordinates
        long nat() { return m_nat; } //get the number of atoms in the configuration
        int ipm(int n); //modular function for pentagon
        int ihm(int n); //modular function for heptagon

        double relax(const std::string& lmp_command);  //optimise the configuration by calling LAMMPS and return the final energy

        void write(const std::string& configfile); //write the system coordintes to configfile in XYZ format
        void writePoly(const std::string& polyfile); //write pentagon and hepatgon coordinates to file
        void writeBonds(const std::string& bondfile); //write the active bond coordinates to file
        void writeEn(const std::string& enfile); //write the configuration system energy to file

    private:
        long m_nat; //number of atoms
        std::vector<Atom> m_coords; //atomic coordinates

        //working arrays for structure scanning
        std::vector<long> m_b1;  //atomic bond table
        std::vector<long> m_b2;
        std::vector<long> m_b3;
        std::vector<long> m_b4;
        std::vector<int> m_nb; //number of bonds (per atom)

        std::vector<long> m_n21;  // n2 vectors with bond list
        std::vector<long> m_n22;
        std::vector<long> m_n31;  // n3 list with all 3 sequentially bonded atoms
        std::vector<long> m_n32;
        std::vector<long> m_n33;
        std::vector<long> m_n41;  // n4 list with all 4 sequentially bonded atoms
        std::vector<long> m_n42;
        std::vector<long> m_n43;
        std::vector<long> m_n44;
        std::vector<long> m_n51;  //  n5 list with all 5 sequentially bonded atom
        std::vector<long> m_n52;
        std::vector<long> m_n53;
        std::vector<long> m_n54;
        std::vector<long> m_n55;

        std::vector<Polygon> m_pent; //5 fold and 7 fold rings
        std::vector<Polygon> m_hept;
        std::vector<Bond> m_bonds; //list of all the potential bond rotation pairs (5-7 intersects)
        double m_en; //configuration energy
        double bounds [6]; //periodic boundary conditions
};

#endif // CONFIG_H
