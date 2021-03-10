#ifndef _libscf_src_hf_h_
#define _libscf_src_hf_h_

#include <libint2.hpp>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "src/molecule.hpp"
#include "src/iwl.hpp"

namespace scf {
using real_t = libint2::scalar_type;
typedef Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
using libint2::Shell;
using libint2::Operator;
using libint2::Engine;

int compute_nao(const std::vector<Shell>& shells);
int max_l(const std::vector<libint2::Shell>& shells);
std::vector<size_t> map_shell_to_basis_function(const std::vector<libint2::Shell>& shells);
std::vector<size_t> map_basis_to_shell_function(const std::vector<libint2::Shell>& shells);
size_t max_nprim(const std::vector<libint2::Shell>& shells);

class HF {
public:
    std::shared_ptr<Molecule> molecule;
    int nao;
    Matrix S;
    Matrix T;
    Matrix V;
    Matrix H;
    Matrix D;
    Matrix G;
    IWL iwl_buf;
    bool computed;
    std::vector<size_t> shell2bf;
    std::vector<size_t> bf2shell;
    double nuc_rep_energy;

    HF(std::shared_ptr<Molecule> molecule_);
    ~HF();
    void init();
    void compute(const int maxiter = 100, const double conv = 1e-12);
    void compute_1body_ints(Matrix& M, const std::vector<Shell>& shells, Operator optype, 
        const std::vector<int>& atom_numbers = std::vector<int>(), const std::vector<coordinate>& atom_positions = std::vector<coordinate>());
    void compute_2body_fock(Matrix& M, const std::vector<Shell>& shells, const Matrix& D);
};
}

#endif