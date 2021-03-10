#ifndef _to_atomic_num_h
#define _to_atomic_num_h

#include <map>
#include <string>
#include <utility>

std::map<std::string, int> label2atomicnumber = 
{
    std::pair<std::string, int>{ "H", 1},
    std::pair<std::string, int>{"He", 2},
    std::pair<std::string, int>{"Li", 3},
    std::pair<std::string, int>{"Be", 4},
    std::pair<std::string, int>{ "B", 5},
    std::pair<std::string, int>{ "C", 6},
    std::pair<std::string, int>{ "N", 7},
    std::pair<std::string, int>{ "O", 8},
    std::pair<std::string, int>{ "F", 9},
    std::pair<std::string, int>{"Ne", 10},
};

#endif