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
   const Spline& E_gen, const Spline& potential);

void push_particle(double& x, double& y, const double x_step, const double y_step);
void reflect_from_wall(double& x, double& y, const double dy, const double wall_y,
                       const double Vx, double& Vy, size_t& surf_col_count);


#endif //PARTICLE_HPP
