#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include <optional>
#include <random>
#include <vector>

#include "settings.hpp"
#include "spline.hpp"

struct Settings;

struct PtStat{
    double Vx;
    double Vy;
    double Vz;
    double E;
    size_t surf_col_count;
};

std::optional<PtStat> trace_particle(std::mt19937& rnd_gen, const Settings& opt,
   const Spline& E_gen, const Spline& potential, const Spline& mfp);



#endif //PARTICLE_HPP
