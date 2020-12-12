#ifndef SETTINGS_HPP
#define SETTINGS_HPP

#include <string>
#include <fstream>
#include <vector>

#include "particle.hpp"

struct PtStat;

struct Settings{
    double D; 	//slit height m
    double L; 	//slit length m
    double m; 	//mass kg
    double Q;	//charge C
    double phi_val; 	//Scale for potential curve
    double R; 	//reflection coefficient
    size_t pt_num;
    std::string out_file;
    std::string energy_fname;
    std::string phi_fname;

    void SaveSettings(FILE* savefile) const ;
};

std::ifstream load_helper(const std::string& filename);
void SaveResult(const std::string& filename, const std::vector<PtStat>& stat,
                const Settings& opt);


#endif //SETTINGS_HPP
