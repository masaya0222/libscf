#include "src/molecule.hpp"
#include "src/basisset.hpp"
#include "src/to_atomic_num.hpp"

#include <libint2.hpp>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>

namespace scf {
Molecule::Molecule(std::vector<int> _atom_numbers, std::vector<coordinate> _atom_position, std::vector<libint2::Shell> _shells)
    : atom_numbers(_atom_numbers), atom_position(_atom_position), shells(_shells)
{
    nelectron = std::accumulate(atom_numbers.begin(), atom_numbers.end(), 0);
    std::cout << "Molecule constructor" << std::endl;
}

std::vector<std::string> split(std::string str, std::string del) {
    int first = str.find_first_not_of(del);
    int last = str.find_first_of(del, first+1);

    std::vector<std::string> result;

    while(first < str.size()) {
        std::string subStr(str, first, last - first);
        result.push_back(subStr);

        first = str.find_first_not_of(del, last+1);
        last = str.find_first_of(del, first+1);

        if (last == std::string::npos) {
        last = str.size();
        }
    }
    return result;
}

    
std::shared_ptr<Molecule> get_molecule_from_xyz(std::string filename, std::string basisset_name) {
    filename += ".xyz"; // change filename to xyz filename
    std::ifstream inputFile;
    std::string input_line_buffer;
    
    std::string read_directory = "../";
    inputFile.open(read_directory + filename);

    if(inputFile.fail()) {
        std::cout << "could not access input file" << std::endl;
        exit(1);
    }

    std::vector<int> atoms;
    std::vector<coordinate> atom_pos;
    std::vector<libint2::Shell> shells;
    int count_line = 0;
    while (!inputFile.eof())
    {
        std::getline(inputFile, input_line_buffer);
        if (count_line == 0) {
            atoms.reserve(std::stoi(input_line_buffer));
        } else if (count_line == 1 ){
            ; // inputfile のdescriotionが入る
        } else {
            if (input_line_buffer.empty()) {
                continue;
            }
            std::vector<std::string> v = split(input_line_buffer, " /t");
            
            if (v.size() != 4) {
                std::cout << "input file does not match format " << std::endl;
                exit(1);
            }
            atoms.push_back(label2atomicnumber[v[0]]);
            atom_pos.push_back(coordinate{ stod(v[1]), stod(v[2]) ,stod(v[3]) });
        }
        count_line++;
    }

    std::transform(basisset_name.begin(), basisset_name.end(), basisset_name.begin(), ::tolower);
    std::ifstream basisfile;
    std::string basisset_dir = "../basis/";
    if (basisset_name == "sto-3g" || basisset_name == "sto3g") {
        basisfile.open(basisset_dir + "STO-3G");
    } else {
        std::cout << "can't use '" << basisset_name << "' basis set in this program" << std::endl;
        exit(1);
    }
    if (basisfile.fail()) {
        std::cout << "can't open '" << basisset_name << "' basis set in this program" << std::endl;
    }
    std::vector<std::string> basisset_data;
    int max_atom_num = *std::max_element(atoms.begin(), atoms.end());
    std::string buffer;
    int count = 0;
    while (!basisfile.eof()) {
        std::getline(basisfile, input_line_buffer);
        if (input_line_buffer == "****") {
            basisset_data.push_back(buffer);
            if (count++ == max_atom_num)
                break;
            buffer.clear();
            continue;
        }
        buffer += " " + input_line_buffer;
    }

    buffer.clear();
    
    for(int i = 0; i < atoms.size(); i++) {
        buffer = basisset_data[atoms[i]];
        std::vector<std::string> v = split(buffer, " \t\n");
        std::string angular_str;
        int contr_num ;
        for(int j = 0; j < v.size(); ) {
            if (j == 0 || j == 1) {
                j++;
                continue;
            }
            
            angular_str = v[j++];
            contr_num = std::stoi(v[j++]);
            j++; // S 3 1.00 の1.00部分 なんのため？
            
            libint2::svector<double> alpha(contr_num,0);
            libint2::svector<double> coeff(contr_num,0);
            libint2::svector<double> max_ln_coeff(contr_num,0);
            int c = 0;
            for(auto k : angular_str ){
                int l = (k == 'S') ? 0 : (k == 'P') ? 1 : (k == 'D') ? 2 : (k == 'F') ? 3 : -1;
                for(int l = 0; l < contr_num; l++) {
                    
                    v[j][v[j].size() -4] = 'e'; // Dをeの表示に変える
                    alpha[l] = std::stod(v[j++]);
                    j += c;
                    v[j][v[j].size() -4] = 'e'; // Dをeの表示に変える
                    coeff[l] = std::stod(v[j]);
                    max_ln_coeff[l] = std::log(coeff[l]);
                    j += (angular_str.size()-c) ; // 次のalphaの最初にしておく

                }
                libint2::Shell shell;
                libint2::Shell::Contraction contr;
                contr.l = l;
                contr.pure = false; // 必ずpure関数でない
                contr.coeff = coeff;
                shell.alpha = alpha;
                shell.max_ln_coeff = max_ln_coeff;
                shell.contr = {contr};
                shell.O = { atom_pos[i].x, atom_pos[i].y, atom_pos[i].z };
                // shell.renorm();
                normalize(shell);
                shells.push_back(shell);
                j -= contr_num*(angular_str.size()+1);
                c++;
            }
            j += contr_num*(angular_str.size()+1);
        }
    }

    std::shared_ptr<Molecule> mol(new Molecule( atoms, atom_pos, shells ));
    return mol;
}
}