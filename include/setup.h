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

#ifndef SETUP_H
#define SETUP_H

#include <string>

class Setup //class to read in and store the simultion set-up options
{
    public:
        Setup(std::string& setupfile); // read-file constructor
        virtual ~Setup();

        double tmp() { return m_temp; }
        double coff() { return m_cutoff; }
        int type() { return m_type; }
        int pstep() { return m_pstep; }
        long maxs() { return m_maxs; }
        double etol() { return m_etol; }
        std::string lmp() { return m_lmpcommand; }

    private:
        double m_temp; // effective temperature for MMC algorithm
        double m_cutoff; // nearest neighbour distance cutoff for bond table formation
        int m_type; // simulation type: 1 = minimisation, 2 = MMC
        int m_pstep; // printing interval for output
        long m_maxs; // maximum number of simulation steps
        double m_etol; // tolerance for exiting energy minimisation
        std::string m_lmpcommand; //command for LAMMPS call
};

#endif // SETUP_H
