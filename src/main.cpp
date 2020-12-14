#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <optional>
#include <omp.h>

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
    opt.R = 0.7;
    opt.dt = 5e-14;
    opt.num_threads = 6;
    opt.out_file = "statistics.txt";
    opt.energy_fname = "Maxwel_T_10.00eV.txt";
    opt.phi_fname = "Phi_shape_L_0.04cm.txt";
    opt.mfp_fname = "mfp_P_5.00_Pa.txt";

    //download energy distribution and make generator
    std::ifstream energy_dist_file = load_helper(opt.energy_fname);
    Spline energy(energy_dist_file);
    Spline E_gen = energy.GetGenerator();
    //download mfp dependence on energy
    std::ifstream mfp_file = load_helper(opt.mfp_fname);
    Spline mfp(mfp_file);
    //download potetnial curve
    std::ifstream phi_shape_file = load_helper(opt.phi_fname);
    Spline potential(phi_shape_file);
    potential.Scale(opt.phi_val);

    std::vector<PtStat> particle_statistics;
    omp_set_num_threads(static_cast<int>(opt.num_threads));
    std::vector<size_t> thread_load(opt.num_threads, opt.pt_num/opt.num_threads);
    thread_load.back() += opt.pt_num%opt.num_threads;
    #pragma omp parallel
    {
        size_t tid = static_cast<uint>(omp_get_thread_num());
        std::mt19937 rnd_gen(static_cast<uint>(time(0))+tid);
        size_t counter = 0;
        while(counter<thread_load[tid]){
            auto res = trace_particle(rnd_gen, opt, E_gen, potential, mfp);
            #pragma omp critical
            if(res) particle_statistics.push_back(res.value());
            #pragma omp master
            {
                if(counter%(thread_load[tid]/10)==0){
                    std::cout << (100*counter)/thread_load[tid] << "%" << std::endl;
                }
            }
            counter++;
        }
    }
    SaveResult(opt.out_file, particle_statistics, opt);

    return 0;
}
