#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <optional>
#include <omp.h>
#include <lyra/lyra.hpp>

#include "spline.hpp"
#include "settings.hpp"
#include "particle.hpp"


int main(int argc, const char** argv){

    std::string out_file = "stat_at_slit_end.txt";
    std::string energy_file = "Maxwel_T_10.00eV.txt";
    std::string phi_file = "Phi_shape_L_0.04cm.txt";
    std::string free_path_file = "mfp_P_5.00_Pa.txt";
    size_t num_threads = 6u;
    size_t particle_num = 1000000u;
    double slit_width = 10e-6;
    double slit_length = 4e-4;
    double phi_multiplier = 1;
    double refl_coeff = 0.7;
    double time_step = 5e-14;
    double pt_charge = -1.602e-19;
    double pt_mass = 9.1e-31;
    bool show_help = false;

    auto cli = lyra::cli()
            | lyra::help(show_help)
            | lyra::opt(out_file, "output file")["-o"]["--output"]
            ("Output file name [default = stat_at_slit_end.txt]")
            | lyra::opt(energy_file, "energy file")["-e"]["--energy"]
            ("File with energy distribution [def = Maxwel_T_10.txt]")
            | lyra::opt(phi_file, "phi file")["--phi"]
            ("File with potential shape (will be multiplied!) [def = Phi_shape_L_0.txt]")
            | lyra::opt(free_path_file, "mean free path file")["--mfp"]
            ("File with mfp dependence on the energy [def = mfp_P_5.00_Pa.txt]")
            | lyra::opt(num_threads, "num_threads")["-t"]["--threads"]
            ("Number of threads will be used [def = 6]")
            | lyra::opt(particle_num, "particle number")["--pt_num"]
            ("Number of particles will be ejected [def = 1e6]")
            | lyra::opt(slit_width, "slit width")["-w"]["--width"]
            ("Simulated slit width in m [def = 10e-6]")
            | lyra::opt(slit_length, "slit length")["-l"]["--length"]
            ("Simulated slit length in m [def = 4e-4")
            | lyra::opt(phi_multiplier, "phi_multiplier")["--phi_m"]
            ("Positive valu with which potential curve will be scaled [def = 1]")
            | lyra::opt(refl_coeff, "reflection coefficient")["-r"]["--reflection"]
            ("Particle reflection from walls coefficient [def = 0.7]")
            | lyra::opt(time_step, "time step")["--dt"]
            ("Time step of the particle mover in s [def = 5e-14]")
            | lyra::opt(pt_charge, "particle charge")["-q"]["--charge"]
            ("Simulated particle charge in C [def = 1.602e-19]")
            | lyra::opt(pt_mass, "particle mass")["-m"]["--mass"]
            ("Simulated particle mass in kg [def = 9.1e-31]");
    auto cmd_parse = cli.parse({argc, argv});
    if(show_help){
        std::cout << cli;
        return 0;
    }
    if(!cmd_parse){
        std::cerr << cmd_parse.errorMessage() << '\n';
    }

    Settings opt;
    opt.D = slit_width;
    opt.L = slit_length;
    opt.m = pt_mass;
    opt.Q = pt_charge;
    opt.phi_val = phi_multiplier;
    opt.pt_num = particle_num;
    opt.R = refl_coeff;
    opt.dt = time_step;
    opt.num_threads = num_threads;
    opt.out_file = out_file;
    opt.energy_fname = energy_file;
    opt.phi_fname = phi_file;
    opt.mfp_fname = free_path_file;

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
