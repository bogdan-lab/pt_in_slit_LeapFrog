#include "particle.hpp"


std::optional<PtStat> trace_particle(std::mt19937& rnd_gen, const Settings& opt,
          const Spline& E_gen, const Spline& potential, const Spline& mfp){
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    double energy = E_gen(dist(rnd_gen));
    if(energy<fabs(potential(0) - potential(opt.L))) return std::nullopt; //no chance particle will reach!
    double Vscale = sqrt(2*energy*1.602e-19/opt.m);
    TVec V = get_random_velocity(rnd_gen, std::make_pair(0.0, 1.0), Vscale);
    //Vz is not traced since in this direction slit is infinite
    TVec pos = {0.0, opt.D/2, 0.0};

    auto F = [&](double cur_x){
        double field = -(potential(cur_x+1e-6)-potential(cur_x))/1e-6;
        return opt.Q*field;
    };

    //LeapFrog
    double Vx_prev = V[0] - 0.5*opt.dt*F(pos[0])/opt.m;
    size_t surf_col_count = 0;
    size_t vol_col_count = 0;
    while(pos[0]>=0){
        V[0] = Vx_prev + opt.dt*F(pos[0])/opt.m;
        energy = get_energy_eV(V, opt);
        double mean_free_path = mfp(energy);
        TVec step = {V[0]*opt.dt, V[1]*opt.dt, V[2]*opt.dt};
        double y_next = pos[1]+step[1];
        if(y_next>opt.D){
            if(dist(rnd_gen)>opt.R) break; 	//died on wall
            if(reflect_from_wall(pos, step, opt.D, V, rnd_gen, mean_free_path, surf_col_count)){
                if(pos[0]>=0) Vx_prev= V[0] - 0.5*opt.dt*F(pos[0])/opt.m;
                vol_col_count++;
                continue;
            }
        }
        else if(y_next<0){
            if(dist(rnd_gen)>opt.R) break; 	//died on wall
            if(reflect_from_wall(pos, step, 0, V, rnd_gen, mean_free_path, surf_col_count)){
                if(pos[0]>=0) Vx_prev = V[0] - 0.5*opt.dt*F(pos[0])/opt.m;
                vol_col_count++;
                continue;
            }
        }
        else{
            if(push_particle(pos, V, step, rnd_gen, mean_free_path)){
                if(pos[0]>=0) Vx_prev = V[0] - 0.5*opt.dt*F(pos[0])/opt.m;
                vol_col_count++;
                continue;
            }
        }
        Vx_prev = V[0];
        if(pos[0]>=opt.L){
            //successfull particle
            PtStat info;
            info.V = V;
            info.surf_col_count = surf_col_count;
            info.vol_col_count = vol_col_count;
            info.E = get_energy_eV(V, opt);
            return info;
        }
    }
    return std::nullopt;
}


bool push_particle(TVec& pos, TVec& V, const TVec& step, std::mt19937& rnd_gen,
                   const double mean_free_path){
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    double dl = sqrt(step[0]*step[0] + step[1]*step[1] + step[2]*step[2]);
    double l_free = mean_free_path*log(1.0/(1.0-dist(rnd_gen)));
    if(l_free < dl){
        //collide in gas
        double Vscale = sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
        double t_before = l_free/Vscale;
        pos[0] += V[0]*t_before;
        pos[1] += V[1]*t_before;
        pos[2] += V[2]*t_before;
        V = get_random_velocity(rnd_gen, std::make_pair(-1.0, 1.0), Vscale);
        return true;
    }
    pos[0] += step[0];
    pos[1] += step[1];
    return false;
}


bool reflect_from_wall(TVec& pos, const TVec& step, const double wall_y,
                       TVec& V, std::mt19937& rnd_gen, const double mean_free_path,
                       size_t& surf_col_count){
    double dy_before = wall_y - pos[1];
    double dt_before = fabs(dy_before/V[1]);
    if(push_particle(pos, V, {V[0]*dt_before, dy_before, 0.0}, rnd_gen, mean_free_path)) return true;
    V[1] *=-1;
    surf_col_count++;
    double dy_after = dy_before - step[1];
    double dt_after = fabs(dy_after/V[1]);
    return push_particle(pos, V, {V[0]*dt_after, dy_after, 0.0}, rnd_gen, mean_free_path);
}


TVec get_random_velocity(std::mt19937& rnd_gen, std::pair<double, double> cos_span,
                                const double Vscale){
    std::uniform_real_distribution<double> dist(cos_span.first, cos_span.second);
    double cos_alpha = dist(rnd_gen);
    double sin_alpha = sqrt(1-cos_alpha*cos_alpha);
    double phi = dist(rnd_gen)*2*M_PI;
    TVec V ={cos_alpha, sin(phi)*sin_alpha, cos(phi)*cos_alpha};
    std::for_each(V.begin(), V.end(), [&](double& d){d*=Vscale;});
    return V;
}

double get_energy_eV(const TVec& V, const Settings& opt){
    return 0.5*opt.m*(V[0]*V[0]+V[1]*V[1]+V[2]*V[2])/1.602e-19;
}
