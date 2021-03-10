#include <libint2.hpp>
#include <vector>
#include "src/basisset.hpp"


namespace scf {
using libint2::Shell;

int compute_nao(const std::vector<Shell>& shells) {
    int count = 0;
    for(Shell shell: shells) {
        for(auto c: shell.contr) {
            count += c.cartesian_size();
        }
    }
    return count;
}

int max_l(const std::vector<libint2::Shell>& shells) {
    int l = 0;
    for (auto shell: shells) 
        for (auto c: shell.contr) 
            l = std::max(c.l, l);
    return l;
}

std::vector<size_t> map_shell_to_basis_function(const std::vector<libint2::Shell>& shells) {
    std::vector<size_t> result;
    result.reserve(shells.size());

    size_t n = 0;
    for (auto shell: shells) {
        result.push_back(n);
        n += shell.size();
    }

    return result;
}

std::vector<size_t> map_basis_to_shell_function(const std::vector<libint2::Shell>& shells) {
    std::vector<size_t> result;
    result.reserve(compute_nao(shells));
    size_t n = 0;
    for (auto shell: shells) {
        for(int i = 0; i < shell.size(); i++) {
            result.push_back(n);
        }
        n++;
    }
    return result;
}

size_t max_nprim(const std::vector<libint2::Shell>& shells) {
    size_t n = 0;
    for (auto shell: shells) {
        n = std::max(shell.nprim(), n);
    }
    return n;
}

void normalize(libint2::Shell& shell) {
    using libint2::math::df_Kminus1;
    using std::pow;
    const auto sqrt_Pi_cubed = double{5.56832799683170784528481798212};;
    const auto np = shell.nprim();
    for (auto& c : shell.contr) {
        assert(c.l <= 15); // due to df_Kminus1[] a 64-bit integer type; kinda ridiculous restriction anyway
        for (auto p=0ul; p!=np; ++p) {
            assert(shell.alpha[p] >=0);
            if (shell.alpha[p] != 0) {
                const auto two_alpha = 2 * shell.alpha[p];
                const auto two_alpha_to_am32 = pow(two_alpha, c.l+1) * sqrt(two_alpha);
                auto normalization_factor = sqrt(pow(2, c.l) * two_alpha_to_am32/(sqrt_Pi_cubed * df_Kminus1[2*c.l]));
                c.coeff[p] *= normalization_factor;
            }
        }
        // compute the self-overlap of the , scale coefficients by its inverse square root
        double norm{0};
        for (auto p=0ul; p!=np; ++p) {
            for (decltype(p) q=0ul; q<=p; ++q) {
                auto gamma = shell.alpha[p] + shell.alpha[q];
                norm += (p==q ? 1 : 2) * df_Kminus1[2*c.l] * sqrt_Pi_cubed * c.coeff[p] * c.coeff[q] /
                        (pow(2,c.l) * pow(gamma,c.l+1) * sqrt(gamma));
            }
        }
        auto normalization_factor = 1 / sqrt(norm);
        for (auto p = 0ul; p!=np; ++p) {
            c.coeff[p] *= normalization_factor;
        }
    }

    shell.max_ln_coeff.resize(np);
    for(auto p=0ul; p!=np; ++p) {
        double max_ln_c = -std::numeric_limits<double>::max();
        for (auto& c: shell.contr) {
            max_ln_c = std::max(max_ln_c, std::log(std::abs(c.coeff[p])));
        }
        shell.max_ln_coeff[p] = max_ln_c;
    }
}

IntegralIterator::IntegralIterator(std::vector<Shell>& shells_){
    current = Integral{0,0,0,0};
    shells = shells_.data();
    shell_size = shells_.size();
    done = false;
}

void IntegralIterator::first() {
    current.i = 0;
    current.j = 0;
    current.k = 0;
    current.l = 0;

    isize = shells[current.i].size();
    jsize = shells[current.j].size();
    ksize = shells[current.k].size();
    lsize = shells[current.l].size();

    imax = shell_size-1; jmax = 0; kmax = 0; lmax = 0;
}

void IntegralIterator::next() {
    current.l++;
    if (current.l > lmax){
        current.l = 0;
        current.k++;
        if (current.k > kmax) {
            current.k = 0;
            current.j++;
            if (current.j > jmax) {
                current.j = 0;
                current.i++;
                if (current.i > imax) {
                    done = true;
                    return;
                }
                jmax = current.i;
                ishell = shells[current.i];
                isize = ishell.size();
            }
            kmax = current.i;
            jshell = shells[current.j];
            jsize = jshell.size();
        }
        lmax = (current.i == current.k) ? current.j : current.k;
        kshell = shells[current.k];
        ksize = kshell.size();
    }
    lshell = shells[current.l];
    lsize = lshell.size();
}

}