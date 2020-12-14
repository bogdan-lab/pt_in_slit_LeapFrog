#include "particle.hpp"


std::optional<PtStat> trace_particle(std::mt19937& rnd_gen, const Settings& opt,
          const Spline& E_gen, const Spline& potential){
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
        //double dz = Vz*dt;
        //double dl = sqrt(dx*dx + dy*dy + dz*dz);
        if(y+dy>opt.D){
            if(dist(rnd_gen)>opt.R) break;
            reflect_from_wall(x, y, dy, opt.D, Vx, Vy, surf_col_count);
        }
        else if(y+dy<0){
            if(dist(rnd_gen)>opt.R) break;
            reflect_from_wall(x, y, dy, 0, Vx, Vy, surf_col_count);
        }
        else{
            push_particle(x, y, dx, dy);
        }
        Vx_prev = Vx;
        x_trace.push_back(x);
        y_trace.push_back(y);
        if(x>=opt.L){
            //successfull particle
            PtStat info;
            info.Vx = Vx;
            info.Vy = Vy;
            info.Vz = Vz;
            info.surf_col_count = surf_col_count;
            info.E = 0.5*opt.m*(Vx*Vx+Vy*Vy+Vz*Vz)/1.602e-19; 	//eV
            break;
            //return info;
        }
    }

    FILE* out_f = fopen("test_trace.txt", "w");
    if(!out_f){
        fprintf(stderr, "CANNOT OPEN TEST TRACE FILE\n");
        exit(EXIT_FAILURE);
    }
    fprintf(out_f, "#x, cm\ty, cm\n");
    for(size_t i=0; i<x_trace.size(); i++){
        fprintf(out_f, "%.10e\t%.10e\n", x_trace[i], y_trace[i]);
    }
    fclose(out_f);

    return std::nullopt;
}


void push_particle(double& x, double& y, const double x_step, const double y_step){
    x += x_step;
    y += y_step;
}


void reflect_from_wall(double& x, double& y, const double dy, const double wall_y,
                       const double Vx,double& Vy, size_t& surf_col_count){
    double dy_before = wall_y - y;
    double dt_before = fabs(dy_before/Vy);
    push_particle(x, y, Vx*dt_before, dy_before);
    Vy *=-1;
    double dy_after = dy_before - dy;
    double dt_after = fabs(dy_after/Vy);
    push_particle(x, y, Vx*dt_after, dy_after);
    surf_col_count++;
}
