#ifndef _libscf_src_molecule_h_
#define _libscf_src_molecule_h_

#include <string>
#include <vector>
#include <libint2.hpp>


namespace scf
{   
    using libint2::Shell;

    struct coordinate {
        double x;
        double y;
        double z;
    };

    class Molecule
    {
    public:
        std::vector<int> atom_numbers;
        std::vector<coordinate> atom_position;
        std::vector<Shell> shells;
        int nelectron;
        
        Molecule(std::vector<int> _atom_numbers, std::vector<coordinate> _atom_position , std::vector<Shell> _shells);
        ~Molecule() { 
            std::cout << "Molecule was destructed"  << std::endl; 
        }
    };

    std::shared_ptr<Molecule> get_molecule_from_xyz(std::string filename, std::string basisset_name);
    
}

#endif
