#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <optional>
#include <random>
#include <vector>
#include <array>

#include "settings.hpp"
#include "spline.hpp"

struct Settings;

typedef std::array<double, 3> TVec;

struct PtStat{
    TVec V;
    double E;
    size_t surf_col_count;
    size_t vol_col_count;
};

std::optional<PtStat> trace_particle(std::mt19937& rnd_gen, const Settings& opt,
   const Spline& E_gen, const Spline& potential, const Spline& mfp);

bool push_particle(TVec& pos, TVec& V, const TVec& step, std::mt19937& rnd_gen,
                   const double mean_free_path);

bool reflect_from_wall(TVec& pos, const TVec& step, const double wall_y,
                       TVec& V, std::mt19937& rnd_gen, const double mean_free_path,
                       size_t& surf_col_count);

TVec get_random_velocity(std::mt19937& rnd_gen, std::pair<double, double> cos_span,
                         const double Vscale);
double get_energy_eV(const TVec& V, const Settings& opt);
#endif //PARTICLE_HPP
