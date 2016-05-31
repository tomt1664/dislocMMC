//class for atomic configuration of the system

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
        Config(std::string infile); //construct object with input file
        virtual ~Config();
        int getbondlist(); // get the bond list
        long getn2(double cutoff); //find all bonds
        long getn3(); //find all 3 sequentially bonded atoms
        long getn4(); // find all 4 sequentially bonded atoms
        long getn5(); // find all 5 sequentially bonded atoms

        long getPents(std::vector<Polygon>& pents); //get a list of all the pentagons in the system (return number)
        long getHepts(std::vector<Polygon>& hepts); //get a list of all the heptagons in the system (return number)
        long getBonds(std::vector<Bond>& bonds); //get a list of all the bond rotation pairs (return number)
        long rotate(long ibnd); //rotate the bond pair of index ibnd (return success)
        long getCoords(std::vector<Atom>& coords); //get the atomic coordinates
        long nat() { return m_nat; } //get the number of atoms in the configuration
    private:
        long m_nat; //number of atoms
        std::vector<Atom> m_coords; //atomic coordinates
        double bounds [5]; //periodic boundary conditions

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

};

#endif // CONFIG_H
