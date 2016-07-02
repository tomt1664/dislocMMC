#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <vector>
#include "config.h"
#include "polygon.h"

// input constructor
Config::Config(std::string& infile)
{
    //open data input file for inital structure
    std::ifstream initstruct(infile.c_str());
    if (!initstruct)
    {
        std::cerr << "Error: input structure could not be opened" << std::endl;
        exit(1);
    }

    //read from file (enhanced xyz file format)
    std::string data;
    std::stringstream ss; //use stringstream
    std::string attype;
    double x,y,z; //atom coordinates

    std::getline(initstruct,data);
    ss << data;
    ss >> m_nat;
    std::cout << "Input structure: " << m_nat << " atoms" << std::endl;
    ss.clear();
    ss.str(""); //clear stringsteam
    std::getline(initstruct,data);  //get PBC
    ss << data;
    ss >> bounds[0] >> bounds[1] >> bounds[2] >> bounds[3] >> bounds[4] >> bounds[5];
    ss.clear();
    ss.str(""); //clear stringsteam


    for(long iat = 0; iat < m_nat; iat++)
    {
        if(initstruct.eof())
        {
            std::cerr << "Error: incomplete input configuration" << std::endl;
            exit(1);
        }
        std::getline(initstruct,data);
        ss << data;
        ss >> attype >> x >> y >> z;
        ss.clear();
        ss.str(""); //clear stringsteam
        Atom atom(x,y,z);
        m_coords.push_back(atom);
    }
    std::cout << "Read in " << m_coords.size() << " atoms" << std::endl;
}

Config::~Config() //default destructor
{
    //dtor
}

//determine the bond list for the configuration: return 1 if any overcoordination
int Config::getbondlist(double cutoff) // get the bond list
{
    int overc = 0;
    //double loop over all atoms in the configuration
    for(int i=0; i < m_nat; i++)
    {
        long bnd[4] = {0,0,0,0};
        int ib = 0;
        for(int ii=0; ii < m_nat; ii++)
        {
            if(i != ii)
            {
                double bd = m_coords[i].dist(m_coords[ii]);
                if(bd < cutoff) //if separation within cutoff then create bond
                {
                    if(ib > 3)
                    {
                        overc = 1;
                    } else
                    {
                        bnd[ib] = ii+1;
                        ib++;
                    }
                }
            }
        }
        //add bond list to vector
        m_b1.push_back(bnd[0]);
        m_b2.push_back(bnd[1]);
        m_b3.push_back(bnd[2]);
        m_b4.push_back(bnd[3]);
        m_nb.push_back(ib);

    }
    return overc;
}

//get the list of all bonds and return the number of bonds
long Config::getn2()
{
    long numbnd = 0;
    for(int i=0; i < m_nat; i++)
    {
        if(m_b1[i] > i+1)
        {
            numbnd++;
            m_n21.push_back(i);
            m_n22.push_back(m_b1[i]-1);
        }
        if(m_b2[i] > i+1)
        {
            numbnd++;
            m_n21.push_back(i);
            m_n22.push_back(m_b2[i]-1);
        }
        if(m_b3[i] > i+1)
        {
            numbnd++;
            m_n21.push_back(i);
            m_n22.push_back(m_b3[i]-1);
        }
        if(m_b4[i] > i+1)
        {
            numbnd++;
            m_n21.push_back(i);
            m_n22.push_back(m_b4[i]-1);
        }
    }
    return numbnd;
}

//get the list of all angles (n3) and return the number
long Config::getn3()
{
    long numn3 = 0;
    long numn2 = m_n21.size();
    for(int i=0; i < numn2; i++)
    {
        for(int j = i+1; j < numn2; j++)
        {
            if(m_n21[i] == m_n21[j])
            {
                m_n31.push_back(m_n22[i]);
                m_n32.push_back(m_n21[i]);
                m_n33.push_back(m_n22[j]);
                numn3++;
            }
            else if(m_n22[i] == m_n21[j])
            {
                m_n31.push_back(m_n21[i]);
                m_n32.push_back(m_n22[i]);
                m_n33.push_back(m_n22[j]);
                numn3++;
            }
            else if(m_n22[i] == m_n22[j])
            {
                m_n31.push_back(m_n21[i]);
                m_n32.push_back(m_n22[i]);
                m_n33.push_back(m_n21[j]);
                numn3++;
            }
            else if(m_n21[i] == m_n22[j])
            {
                m_n31.push_back(m_n22[i]);
                m_n32.push_back(m_n21[i]);
                m_n33.push_back(m_n21[j]);
                numn3++;
            }

        }
    }
    return numn3;
}

//get the list of all torsions (n4) and return the number
long Config::getn4()
{
    long numn4 = 0;
    long numn3 = m_n31.size();
    for(int i=0; i < numn3; i++)
    {
        for(int j = i+1; j < numn3; j++)
        {
            if(m_n32[i] == m_n31[j] && m_n33[i] == m_n32[j])
            {
                m_n41.push_back(m_n31[i]);
                m_n42.push_back(m_n32[i]);
                m_n43.push_back(m_n33[i]);
                m_n44.push_back(m_n33[j]);
                numn4++;
            }
            else if(m_n31[i] == m_n32[j] && m_n32[i] == m_n33[j])
            {
                m_n41.push_back(m_n31[j]);
                m_n42.push_back(m_n31[i]);
                m_n43.push_back(m_n32[i]);
                m_n44.push_back(m_n33[i]);
                numn4++;
            }
            else if(m_n32[i] == m_n33[j] && m_n33[i] == m_n32[j])
            {
                m_n41.push_back(m_n31[i]);
                m_n42.push_back(m_n32[i]);
                m_n43.push_back(m_n33[i]);
                m_n44.push_back(m_n31[j]);
                numn4++;
            }
            else if(m_n31[i] == m_n32[j] && m_n32[i] == m_n31[j])
            {
                m_n41.push_back(m_n33[j]);
                m_n42.push_back(m_n32[j]);
                m_n43.push_back(m_n31[j]);
                m_n44.push_back(m_n33[i]);
                numn4++;
            }
            else if(m_n32[i] == m_n31[j] && m_n31[i] == m_n32[j])
            {
                m_n41.push_back(m_n33[j]);
                m_n42.push_back(m_n31[i]);
                m_n43.push_back(m_n32[i]);
                m_n44.push_back(m_n33[i]);
                numn4++;
            }
            else if(m_n33[i] == m_n32[j] && m_n32[i] == m_n33[j])
            {
                m_n41.push_back(m_n31[i]);
                m_n42.push_back(m_n32[i]);
                m_n43.push_back(m_n33[i]);
                m_n44.push_back(m_n31[j]);
                numn4++;
            }
            else if(m_n32[i] == m_n33[j] && m_n31[i] == m_n32[j])
            {
                m_n41.push_back(m_n31[j]);
                m_n42.push_back(m_n32[j]);
                m_n43.push_back(m_n33[j]);
                m_n44.push_back(m_n33[i]);
                numn4++;
            }
            else if(m_n33[i] == m_n32[j] && m_n32[i] == m_n31[j])
            {
                m_n41.push_back(m_n31[i]);
                m_n42.push_back(m_n32[i]);
                m_n43.push_back(m_n33[i]);
                m_n44.push_back(m_n33[j]);
                numn4++;
            }
        }
    }
    return numn4;
}

//get the list of all n5 and return the number
long Config::getn5()
{
    long numn5 = 0;
    long numn4 = m_n41.size();
    for(int i=0; i < numn4; i++)
    {
        for(int j = i+1; j < numn4; j++)
        {
            if(m_n42[i] == m_n41[j] && m_n43[i] == m_n42[j] && m_n44[i] == m_n43[j])
            {
                m_n51.push_back(m_n41[i]);
                m_n52.push_back(m_n42[i]);
                m_n53.push_back(m_n43[i]);
                m_n54.push_back(m_n44[i]);
                m_n55.push_back(m_n44[j]);
                numn5++;
            }
            else if(m_n41[i] == m_n42[j] && m_n42[i] == m_n43[j] && m_n43[i] == m_n44[j])
            {
                m_n51.push_back(m_n41[j]);
                m_n52.push_back(m_n41[i]);
                m_n53.push_back(m_n42[i]);
                m_n54.push_back(m_n43[i]);
                m_n55.push_back(m_n44[i]);
                numn5++;
            }
            else if(m_n42[i] == m_n44[j] && m_n43[i] == m_n43[j] && m_n44[i] == m_n42[j])
            {
                m_n51.push_back(m_n41[i]);
                m_n52.push_back(m_n42[i]);
                m_n53.push_back(m_n43[i]);
                m_n54.push_back(m_n44[i]);
                m_n55.push_back(m_n41[j]);
                numn5++;
            }
            else if(m_n41[i] == m_n43[j] && m_n42[i] == m_n42[j] && m_n43[i] == m_n41[j])
            {
                m_n51.push_back(m_n44[j]);
                m_n52.push_back(m_n41[i]);
                m_n53.push_back(m_n42[i]);
                m_n54.push_back(m_n43[i]);
                m_n55.push_back(m_n44[i]);
                numn5++;
            }
        }
    }
    return numn5;
}

//get a list of all the pentagons in the system
long Config::getPents(std::vector<Polygon>& pents)
{
    long np = 0;
    long numn4 = m_n41.size();
    long numn3 = m_n31.size();
    for(int i=0; i < numn4; i++)
    {
        for(int j = 0; j < numn3; j++)
        {
            if(m_n41[i] == m_n31[j] && m_n44[i] == m_n33[j])
            {
                std::vector<long> indx;
                std::vector<Atom> atoms;
                indx.push_back(m_n41[i]);
                indx.push_back(m_n42[i]);
                indx.push_back(m_n43[i]);
                indx.push_back(m_n44[i]);
                indx.push_back(m_n32[j]);

                for(int n=0; n < 5; n++)
                {
                    atoms.push_back(m_coords[indx[n]]);
                }

                Polygon tpent(atoms,indx);

                long psum = tpent.sum();
                long psize = m_pent.size();
                int padd = 1;
                for(long ps=0; ps < psize; ps++)
                {
                    long psum2 = m_pent[ps].sum();
                    if(psum == psum2)
                    {
                        padd = 0;
                    }

                }
                if(padd)
                {
                    m_pent.push_back(tpent);
                    pents.push_back(tpent);
                    np++;
                }
            }
            else if(m_n44[i] == m_n31[j] && m_n41[i] == m_n33[j])
            {
                std::cout << "l2" << std::endl;


                std::vector<long> indx;
                std::vector<Atom> atoms;
                indx.push_back(m_n41[i]);
                indx.push_back(m_n42[i]);
                indx.push_back(m_n43[i]);
                indx.push_back(m_n44[i]);
                indx.push_back(m_n32[j]);

                for(int n=0; n < 5; n++)
                {
                    atoms.push_back(m_coords[indx[n]]);
                }

                Polygon tpent(atoms,indx);
                long psum = tpent.sum();
                long psize = m_pent.size();
                int padd = 1;
                for(long ps=0; ps < psize; ps++)
                {
                    long psum2 = m_pent[ps].sum();
                    if(psum == psum2)
                    {
                        padd = 0;
                    }
                }
                if(padd)
                {
                    m_pent.push_back(tpent);
                    pents.push_back(tpent);
                    np++;
                }
            }
        }
    }
    return np;
}

//get a list of all the heptagons in the system
long Config::getHepts(std::vector<Polygon>& hepts)
{
    long nh = 0;
    long numn5 = m_n51.size();
    long numn4 = m_n41.size();
    for(int i=0; i < numn5; i++)
    {
        for(int j = 0; j < numn4; j++)
        {
            if(m_n51[i] == m_n41[j] && m_n55[i] == m_n44[j] && m_n52[i] != m_n42[j] && m_n54[i] != m_n43[j])
            {
                std::vector<long> indx;
                std::vector<Atom> atoms;
                atoms.clear();
                indx.clear();
                indx.push_back(m_n51[i]);
                indx.push_back(m_n52[i]);
                indx.push_back(m_n53[i]);
                indx.push_back(m_n54[i]);
                indx.push_back(m_n55[i]);
                indx.push_back(m_n43[j]);
                indx.push_back(m_n42[j]);

                for(int n=0; n < 7; n++)
                {
                    atoms.push_back(m_coords[indx[n]]);
                }

                Polygon thept(atoms,indx);

                long hsum = thept.sum();
                long hsize = m_hept.size();
                int hadd = 1;
                for(long hs=0; hs < hsize; hs++)
                {
                    long hsum2 = m_hept[hs].sum();
                    if(hsum == hsum2)
                    {
                        hadd = 0;
                    }

                }
                if(hadd)
                {
                    m_hept.push_back(thept);
                    hepts.push_back(thept);
                    nh++;
                }
            }

            else if(m_n55[i] == m_n41[j] && m_n51[i] == m_n44[j] && m_n52[i] != m_n43[j] && m_n54[i] != m_n42[j])
            {
                std::vector<long> indx;
                std::vector<Atom> atoms;
                atoms.clear();
                indx.clear();
                indx.push_back(m_n51[i]);
                indx.push_back(m_n52[i]);
                indx.push_back(m_n53[i]);
                indx.push_back(m_n54[i]);
                indx.push_back(m_n55[i]);
                indx.push_back(m_n42[j]);
                indx.push_back(m_n43[j]);

                for(int n=0; n < 7; n++)
                {
                    atoms.push_back(m_coords[indx[n]]);
                }

                Polygon thept(atoms,indx);

                long hsum = thept.sum();
                long hsize = m_hept.size();
                int hadd = 1;
                for(long hs=0; hs < hsize; hs++)
                {
                    long hsum2 = m_hept[hs].sum();
                    if(hsum == hsum2)
                    {
                        hadd = 0;
                    }

                }
                if(hadd)
                {
                    m_hept.push_back(thept);
                    hepts.push_back(thept);
                    nh++;
                }
            }
        }
    }
    return nh;
}

//get a list of all bonds in a heptagon where it intercepts with a pentagon
long Config::getBonds(std::vector<Bond>& bonds)
{
    long nb = 0;
    long nump = m_pent.size();
    long numh = m_hept.size();
    for(int i=0; i < nump; i++)
    {
        std::vector<long> pindx;
        m_pent[i].index(pindx);
        for(int j = 0; j < numh; j++)
        {
            std::vector<long> hindx;
            m_hept[j].index(hindx);
                for(int n=0; n < 5; n++)
                {
                    for(int m=0; m < 7; m++)
                    {
                        if(pindx[n] == hindx[m])
                        {
                            long indx1 = pindx[n];
                            long indx2;
//                            std::cout << indx1 << std::endl;
                            if(hindx[ihm(m+1)] == pindx[ipm(n+1)])
                            {
                                indx2 = hindx[ihm(m-1)];
                            } else if(hindx[ihm(m-1)] == pindx[ipm(n-1)])
                            {
                                indx2 = hindx[ihm(m+1)];
                            } else if(hindx[ihm(m+1)] == pindx[ipm(n-1)])
                            {
                                indx2 = hindx[ihm(m-1)];
                            } else if(hindx[ihm(m-1)] == pindx[ipm(n+1)])
                            {
                                indx2 = hindx[ihm(m+1)];
                            } else
                            {
                                std::cout << "Error in getBonds" << std::endl;
                                exit(1);
                            }
                            Bond bond(m_coords[indx1],m_coords[indx2],indx1,indx2);
                            m_bonds.push_back(bond);
                            bonds.push_back(bond);
                            nb++;
                        }
                    }
                }
        }
    }
    return nb;
}

//write the configuration to file (in enhanced XYZ format)
void Config::write(std::string& configfile)
{
    //open data input file for inital structure
    std::ofstream writeconfig(configfile.c_str(),std::ios::app);

    writeconfig << m_nat << std::endl;
    for(int i=0; i < 6; i++)
    {
        writeconfig << bounds[i] << " ";
    }
    writeconfig << std::endl;

    //write coordinates
    for(long i=0; i < m_nat; i++)
    {
        writeconfig << "C " << m_coords[i].getx() << " " << m_coords[i].gety() << " " << m_coords[i].getz() << std::endl;
    }
    writeconfig.close();
}

//write the polygons to file (in XYZ format)
void Config::writePoly(std::string& polyfile)
{
    //open data input file for inital structure
    std::ofstream writepoly(polyfile.c_str(),std::ios::app);

    int npoly = m_pent.size()*5 + m_hept.size()*7;
    writepoly << npoly << std::endl;
    writepoly << " " << std::endl;

    //write coordinates
    for(unsigned i=0; i < m_pent.size(); i++)
    {
        std::vector<long> pindx;
        pindx.clear();
        m_pent[i].index(pindx);
        for(int j=0; j < 5; j++)
        {
            writepoly << "P " << m_coords[pindx[j]].getx() << " " << m_coords[pindx[j]].gety() << " " << m_coords[pindx[j]].getz() << std::endl;
        }
    }

    //write coordinates
    for(unsigned i=0; i < m_hept.size(); i++)
    {
        std::vector<long> hindx;
        hindx.clear();
        m_hept[i].index(hindx);
        for(int j=0; j < 7; j++)
        {
            writepoly << "P " << m_coords[hindx[j]].getx() << " " << m_coords[hindx[j]].gety() << " " << m_coords[hindx[j]].getz() << std::endl;
        }
    }
    writepoly.close();
}

//write the polygons to file (in XYZ format)
void Config::writeBonds(std::string& bondfile)
{
    //open data input file for inital structure
    std::ofstream writebond(bondfile.c_str(),std::ios::app);

    int nbnd = m_bonds.size()*2;
    writebond << nbnd << std::endl;
    writebond << " " << std::endl;

    //write coordinates
    for(unsigned i=0; i < m_bonds.size(); i++)
    {
        Atom at1;
        Atom at2;
        m_bonds[i].coords(at1,at2);

        writebond << "B " << at1.getx() << " " << at1.gety() << " " << at1.getz() << std::endl;
        writebond << "B " << at2.getx() << " " << at2.gety() << " " << at2.getz() << std::endl;
    }
    writebond.close();
}

//modular function for pentagons
int Config::ipm(int n)
{
    if(n > 4)
    {
        n -= 5;
    } else if(n < 0)
    {
        n += 5;
    }
    return n;
}

//modular function for hexagons
int Config::ihm(int n)
{
    if(n > 6)
    {
        n -= 7;
    } else if(n < 0)
    {
        n += 7;
    }
    return n;
}

//optimise the configuration by calling LAMMPS and return the final energy
double Config::relax()
{
    double energy = 0.0;
    //write the LAMMPS data file
    std::ofstream lammps("data.rot",std::ios::trunc);
    lammps << "LAMMPS dislocMMC optimisation" << std::endl;
    lammps << " " << std::endl;
    lammps << "      " << m_nat << " atoms" << std::endl;
    lammps << "          0 bonds" << std::endl;
    lammps << "          0 angles" << std::endl;
    lammps << "          0 dihedrals" << std::endl;
    lammps << "          0 impropers" << std::endl;
    lammps << " " << std::endl;
    lammps << "          1 atom types" << std::endl;
    lammps << "          0 bond types" << std::endl;
    lammps << "          0 angle types" << std::endl;
    lammps << "          0 dihedral types" << std::endl;
    lammps << "          0 improper types" << std::endl;
    lammps << " " << std::endl;
    lammps << bounds[0] << " " << bounds[1] << "   xlo xhi" << std::endl;
    lammps << bounds[2] << " " << bounds[3] << "   ylo yhi" << std::endl;
    lammps << bounds[4] << " " << bounds[5] << "   zlo zhi" << std::endl;
    lammps << " " << std::endl;
    lammps << " Masses" << std::endl;
    lammps << " " << std::endl;
    lammps << "          1   12.0000" << std::endl;
    lammps << " " << std::endl;
    lammps << " Atoms" << std::endl;
    lammps << " " << std::endl;
    for(long i=0; i < m_nat; i++)
    {
        lammps << i+1 << "    1 " << m_coords[i].getx() << " " << m_coords[i].gety() << " " << m_coords[i].getz() << std::endl;
    }
    lammps.close();

    //call and run LAMMPS
    system("lmp_serial < in.rot > lmp.out");

    //open output file and retrieve energy value
    std::ifstream logfile("log.lammps");
    std::string line;
    std::stringstream ss; //use stringstream
    int count = 0;
    while(std::getline(logfile, line))
    {
        count++;
        size_t found = line.find("next-to-last");
        if(found != std::string::npos)
        {
            std::getline(logfile, line);

            std::string attype;
            double en_init,en_pen,en_final; //potential energies during the optimisation

            ss << line;
            ss >> en_init >> en_pen >> en_final;
            ss.clear();
            ss.str(""); //clear stringsteam

            std::cout << "line "<< en_final << std::endl;
        }
    }

    m_en = en_final;
    return en_final;
}
