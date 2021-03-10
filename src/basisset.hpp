#ifndef _libscf_src_basisset_h_
#define _libscf_src_basisset_h_

#include <libint2.hpp>
#include <vector>

namespace scf {
using libint2::Shell;

int compute_nao(const std::vector<Shell>& shells);
int max_l(const std::vector<libint2::Shell>& shells);
size_t max_nprim(const std::vector<libint2::Shell>& shells);
std::vector<size_t> map_shell_to_basis_function(const std::vector<libint2::Shell>& shells);
std::vector<size_t> map_basis_to_shell_function(const std::vector<libint2::Shell>& shells);
void normalize(libint2::Shell& shell);

class IntegralIterator {
public:
    struct Integral {
        int i;
        int j;
        int k;
        int l;
    };

    Integral current;
    Shell* shells;
    Shell ishell;
    Shell jshell;
    Shell kshell;
    Shell lshell;
    int shell_size;
    int isize, jsize, ksize, lsize;
    int imax, jmax, kmax, lmax;
    bool done;

    IntegralIterator(std::vector<Shell>& shells_);

    void first();
    void next();
    bool is_done() { return done; }

    int i() const { return current.i; }
    int j() const { return current.j; }
    int k() const { return current.k; }
    int l() const { return current.l; }
};
}

#endif