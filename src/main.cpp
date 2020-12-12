﻿#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <optional>

#include "spline.hpp"
#include "settings.hpp"
#include "particle.hpp"


int main(){

    Settings opt;
    opt.D = 10e-6;
    opt.L = 0.04e-2;
    opt.m = 9.1e-31;
    opt.Q = -1.602e-19;
    opt.phi_val = 2.0;
    opt.pt_num = 1000000;
    opt.R = 0.5;
    opt.out_file = "statistics.txt";
    opt.energy_fname = "Maxwel_T_10.00eV.txt";
    opt.phi_fname = "Phi_shape_L_0.04cm.txt";

    //download energy distribution and make generator
    std::ifstream energy_dist_file = load_helper(opt.energy_fname);
    Spline energy(energy_dist_file);
    Spline E_gen = energy.GetGenerator();
    //download mfp dependence on energy
    //Spline mfp_dep(load_helper(""));
    //download potetnial curve
    std::ifstream phi_shape_file = load_helper(opt.phi_fname);
    Spline potential(phi_shape_file);
    potential.Scale(opt.phi_val);

    std::mt19937 rnd_gen(42);
    std::vector<PtStat> particle_statistics;
    for(size_t i=0; i<opt.pt_num; i++){
        auto res = trace_particle(rnd_gen, opt, E_gen, potential);
        if(res){
            particle_statistics.push_back(res.value());
        }
    }

    SaveResult(opt.out_file, particle_statistics, opt);


    return 0;
}
