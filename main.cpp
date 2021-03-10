#include "src/test.hpp"
#include "src/molecule.hpp"
#include "src/iwl.hpp"
#include "src/basisset.hpp"
#include "src/hf.hpp"
#include <libint2.hpp>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>


int main() {
    hello_world_func();
    std::shared_ptr<scf::Molecule> mol_ptr = scf::get_molecule_from_xyz("c6h12", "sto-3g");
    
    scf::HF hf(mol_ptr);
    
    hf.init();
    hf.compute();
}





