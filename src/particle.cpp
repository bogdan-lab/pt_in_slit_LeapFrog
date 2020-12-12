#include "particle.hpp"


std::optional<PtStat> trace_particle(std::mt19937& rnd_gen, const Settings& opt,
          const Spline& E_gen, const Spline& potential, const Spline& mfp){
    std::uniform_real_distribution<double> dist(0, 1.0);
    double cos_alpha = dist(rnd_gen);
    double phi = dist(rnd_gen)*2*M_PI;
    double V = sqrt(2*E_gen(dist(rnd_gen))*1.602e-19/opt.m);
    double Vx = V*cos_alpha;
    double sin_alpha = sqrt(1-cos_alpha*cos_alpha);
    double Vy = V*sin(phi)*sin_alpha;
    double Vz = V*cos(phi)*sin_alpha;
    //Vz is not traced since in this direction slit is infinite
    double dt = opt.D/(Vy*10);
    double x = 0;
    double y = opt.D/2;
    std::vector<double> x_trace = {x};
    std::vector<double> y_trace = {y};

    auto F = [&](double cur_x){
        double field = -(potential(cur_x+1e-6)-potential(cur_x))/1e-6;
        return opt.Q*field;
    };

    //LeapFrog
    double Vx_prev = Vx - 0.5*dt*F(x)/opt.m;
    size_t surf_col_count = 0;
    while(x>=0){
        Vx = Vx_prev + dt*F(x)/opt.m;
        double dx = Vx*dt;
        double dy = Vy*dt;
        double dz = Vz*dt;
        double dl = sqrt(dx*dx + dy*dy + dz*dz);
        x += dx;
        x_trace.push_back(x);
        y+= dy;
        //WITH REFLECTION I NEED TO TAKE INTO ACCOUNT THAT PARTICLE MOVEMENT
        //AWAY FROM THE WALL NOW I SIMPLY THROW IT AWAY!
        if(y>opt.D){
            y = opt.D;
            if(dist(rnd_gen)>opt.R) return std::nullopt; 	//did not reflect
            Vy *= -1;
            surf_col_count++;
        }
        if(y<0){
            y = 0.0;
            if(dist(rnd_gen)>opt.R) return std::nullopt; 	//did not reflect
            Vy *=-1;
            surf_col_count++;
        }
        y_trace.push_back(y);
        Vx_prev = Vx;
        if(x>=opt.L){
            //successfull particle
            PtStat info;
            info.Vx = Vx;
            info.Vy = Vy;
            info.Vz = Vz;
            info.surf_col_count = surf_col_count;
            info.E = 0.5*opt.m*(Vx*Vx+Vy*Vy+Vz*Vz)/1.602e-19; 	//eV
            return info;
        }
    }
    return std::nullopt;
}

