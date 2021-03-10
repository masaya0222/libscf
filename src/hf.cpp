#include <libint2.hpp>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <chrono>

#include "src/molecule.hpp"
#include "src/hf.hpp"
#include "src/basisset.hpp"

namespace scf {

HF::HF(std::shared_ptr<Molecule> moleucle_) {
    molecule = moleucle_;
    nao = compute_nao(molecule->shells);
    shell2bf = map_shell_to_basis_function(molecule->shells);
    bf2shell = map_basis_to_shell_function(molecule->shells);
    computed = false;
    S = Matrix::Zero(nao, nao);
    T = Matrix::Zero(nao, nao);
    V = Matrix::Zero(nao, nao);
    H = Matrix::Zero(nao, nao);
    D = Matrix::Zero(nao, nao);
    G = Matrix::Zero(nao, nao);
    
    iwl_buf.init_write();
}

HF::~HF() { 
    libint2::finalize(); 
    std::cout << "hf was destructed" << std::endl;
}

void HF::init() {
    // compute the nuclear repulsion energy 
    double enuc = 0.0;
    for (int i = 0; i < molecule->atom_numbers.size(); i++) {
        for (int j = i + 1; j < molecule->atom_numbers.size(); j++) {
            double xij = molecule->atom_position[i].x - molecule->atom_position[j].x;
            double yij = molecule->atom_position[i].y - molecule->atom_position[j].y;
            double zij = molecule->atom_position[i].z - molecule->atom_position[j].z;
            double r2 = xij*xij + yij*yij + zij*zij;
            double r = sqrt(r2);
            enuc += molecule->atom_numbers[i] * molecule->atom_numbers[j] / r;
        }
    }
    nuc_rep_energy = enuc;
    std::cout << "\tNuclear repulsion energy = " << enuc << std::endl;

    // compute 1-e integrals                                                                                             
    // initializes the Libint integrals library ... now ready to compute
    libint2::initialize();

    // compute overlap integrals
    compute_1body_ints(S, molecule->shells, Operator::overlap);
    std::cout << "\n\tOverlap Integrals:" << std::endl;
    std::cout << S << std::endl;

    // compute kinetic-energy integrals
    compute_1body_ints(T, molecule->shells, Operator::kinetic);
    std::cout << "\n\tKinetic-Energy Integrals:" << std::endl;
    std::cout << T << std::endl;

    // compute nuclear-attraction integrals
    compute_1body_ints(V, molecule->shells, Operator::nuclear, molecule->atom_numbers, molecule->atom_position);
    std::cout << "\n\tNuclear Attraction Integrals:" << std::endl;
    std::cout << V << std::endl;

    // Core Hamiltonian H = T + V
    H = T + V;
    std::cout << "\n\tCore Hamiltonian:" << std::endl;
    std::cout << H << std::endl;

    // build initial-guess density
    // solve H C = e S C
    Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(H, S);
    auto eps = gen_eig_solver.eigenvalues();
    auto C = gen_eig_solver.eigenvectors();
    std::cout << "\n\tInitial C Matrix:" << std::endl;
    std::cout << C << std::endl;

    // compute density, D = C(occ) . C(occ)T
    auto C_occ = C.leftCols(molecule->nelectron/2);
    D = C_occ * C_occ.transpose();

    std::cout << "\n\tInitial Density Matrix:" << std::endl;
    std::cout << D << std::endl;

    Engine engine(Operator::coulomb, max_nprim(molecule->shells), max_l(molecule->shells), 0);
    const auto& buf = engine.results();

    IntegralIterator iter(molecule->shells);
    const auto tstart = std::chrono::high_resolution_clock::now();
    for (iter.first(); !iter.done; iter.next()) {
        engine.compute(iter.shells[iter.current.i], iter.shells[iter.current.j], iter.shells[iter.current.k], iter.shells[iter.current.l]);
        const auto* buf_1234 = buf[0];
        if (buf_1234 == nullptr)
            continue;
        iwl_buf.write(shell2bf[iter.current.i], shell2bf[iter.current.j], shell2bf[iter.current.k], shell2bf[iter.current.l], iter.isize, iter.jsize, iter.ksize, iter.lsize, (double*)buf_1234);
    }
    
    iwl_buf.flush(1);
    const auto tstop = std::chrono::high_resolution_clock::now();
    const std::chrono::duration<double> time_elapsed = tstop - tstart;
    std::cout << "computed 2 electron integral" << std::endl;
    printf("takes %10.5lf (s)\n",time_elapsed.count());
    iwl_buf.close_buf();
    iwl_buf.init_read();
}

void HF::compute(const int maxiter, const double conv) {
    int iteration = 0;
    double rmsd = 0.0;
    double ediff = 0.0;
    double ehf = 0.0;
    do {
        const auto tstart = std::chrono::high_resolution_clock::now();
        ++iteration;

        // Save a copy of the energy and the density
        double ehf_last = ehf;
        Matrix D_last = D;

        // build a new Fock matrix
        Matrix F = H;
        compute_2body_fock(G, molecule->shells, D);
        
        F += G;

        if (iteration == 1) {
            std::cout << "\n\tFock Matrix:" << std::endl;
            std::cout << F << std::endl;
        }
        
        // solve F C = e S C
        Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(F, S);

        auto eps = gen_eig_solver.eigenvalues();
        auto C = gen_eig_solver.eigenvectors();

        // compute density, D = C(occ) . C(occ)T
        auto C_occ = C.leftCols(molecule->nelectron/2);
        D = C_occ * C_occ.transpose();

        // compute HF energy
        ehf = 0.0;
        for (int i = 0; i < nao; i++) 
            for (int j = 0; j < nao; j++) 
                ehf += D(i, j) * (H(i, j) + F(i, j));
        // compute difference with last iteration
        ediff = ehf - ehf_last;
        rmsd = (D - D_last).norm();

        const auto tstop = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> time_elapsed = tstop - tstart;

        if (iteration == 1)
            std::cout << "\n\n Iter        E(elec)              E(tot)               Delta(E)             RMS(D)         Time(s)\n";
        printf(" %02d %20.12f %20.12f %20.12f %20.12f %10.5lf\n", iteration, ehf, ehf + nuc_rep_energy, ediff, rmsd, time_elapsed.count());
    } while(((fabs(ediff) > conv) || (fabs(rmsd) > conv)) && (iteration < maxiter));

    printf("** Hartree-Fock energy = %20.12f\n", ehf + nuc_rep_energy);

}

void HF::compute_1body_ints(Matrix& M, const std::vector<Shell>& shells, Operator optype,
        const std::vector<int>& atom_numbers, const std::vector<coordinate>& atom_positions){
    // construct the overlap integrals engine
    Engine engine(optype, max_nprim(shells), max_l(shells), 0);
    // nuclear attraction ints engine needs to know where the charges sit ...
    // the nuclei are charges in this case; in QM/MM there will also classical charges
    if (optype == Operator::nuclear) {
        std::vector<std::pair<real_t, std::array<real_t, 3>>> q;
        for (int i = 0; i < atom_numbers.size(); i++) {
            q.push_back({ (real_t)atom_numbers[i], {{atom_positions[i].x, atom_positions[i].y, atom_positions[i].z}} });
        }
        engine.set_params(q);
    }

    // buf[0] points to the target shell set after every call to engine.compute()                                                                                            
    const auto& buf = engine.results();
    // loop over unique shell pairs, {s1,s2} such that s1 >= s2                                                               
    // this is due to the permutational symmetry of the real integrals over Hermitian operators: (1|2) = (2|1)
    for (int s1 = 0; s1!=shells.size(); ++s1) {
        size_t bf1 = shell2bf[s1];
        size_t n1 = shells[s1].size();
        for (int s2 = 0; s2 <= s1; ++s2) {
            size_t bf2 = shell2bf[s2];
            size_t n2 = shells[s2].size();

            // compute shell pair
            engine.compute(shells[s1], shells[s2]);

            // "map" buffer to a const Eigen Matrix, and copy it to the corresponding blocks of the result
            Eigen::Map<const Matrix> buf_mat(buf[0], n1, n2);
            M.block(bf1, bf2, n1, n2) = buf_mat;
            if (s1 != s2) {
                // if s1 > s2, copy {s1, s2} to the corresponding {s2,s1} block, note the transpose!
                M.block(bf2, bf1, n2, n1) = buf_mat.transpose();
            }
        }
    }
}

void HF::compute_2body_fock(Matrix& M, const std::vector<Shell>& shells, const Matrix& D) {
    int idx, idx_;
    int s1, s2, s3, s4;
    short int p, q, r, s;
    double value, value_scale_by_deg;
    short int* lblptr;
    double* valptr;
    
    double s12_deg, s34_deg, s12_34_deg, s1234_deg;
    G.setZero();
    do
    {
        iwl_buf.fetch();
        idx_ = iwl_buf.index();
        lblptr = iwl_buf.labels();
        valptr = iwl_buf.values();
        
        while (idx_ < iwl_buf.buffer_count())
        {
            idx = 4 * idx_;
            p = lblptr[idx++];
            q = lblptr[idx++];
            r = lblptr[idx++];
            s = lblptr[idx++];

            value = valptr[idx_];
            
            s1 = bf2shell[p]; s2 = bf2shell[q];
            s3 = bf2shell[r]; s4 = bf2shell[s];
            s12_deg = (s1 == s2) ? 1.0 : 2.0;
            s34_deg = (s3 == s4) ? 1.0 : 2.0;
            s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1.0 : 2.0) : 2.0;
            s1234_deg = s12_deg * s34_deg * s12_34_deg;
            
            value_scale_by_deg = value * s1234_deg;
            
            G(p, q) += D(r, s) * value_scale_by_deg;
            G(r, s) += D(p, q) * value_scale_by_deg;
            G(p, r) -= 0.25 * D(q, s) * value_scale_by_deg;
            G(q, s) -= 0.25 * D(p, r) * value_scale_by_deg;
            G(p, s) -= 0.25 * D(q, r) * value_scale_by_deg;
            G(q, r) -= 0.25 * D(p, s) * value_scale_by_deg;
            
            idx_++;
        }
    } while (iwl_buf.last_buffer() == 0);

    iwl_buf.close_buf();
    iwl_buf.init_read();

    Matrix Gt = G.transpose();
    G = 0.5 * (G + Gt);
}

}