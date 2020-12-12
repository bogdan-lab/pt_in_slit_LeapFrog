#include <iostream>
#include <random>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

#include "spline.hpp"

struct Settings{
    double D; 	//slit height m
    double L; 	//slit length m
    double m; 	//mass kg
    double Q;	//charge C
    double E; 	//energy eV
    double phi_val;
    std::string out_file;
};

std::ifstream load_helper(const std::string& filename){
    std::ifstream file(filename);
    if(!file.is_open()){
        fprintf(stderr, "Cannot open file %s \n", filename.c_str());
        exit(EXIT_FAILURE);
    }
    return file;
}

int main(){

    Settings opt;
    opt.D = 10e-6;
    opt.L = 0.04e-2;
    opt.m = 9.1e-31;
    opt.Q = -1.602e-19;
    opt.E = 10.0;
    opt.phi_val = 10.0;
    opt.out_file = "test_trace.txt";

    //download energy distribution and make generator
    //download mfp dependence on energy
    //Spline mfp_dep(load_helper(""));
    //download potetnial curve
    std::ifstream phi_shape_file = load_helper("Phi_shape_L_0.04cm.txt");
    Spline potential(phi_shape_file);
    potential.Scale(opt.phi_val);

    std::mt19937 rnd_gen(42);
    std::uniform_real_distribution<double> dist(0, 1.0);

    double cos_alpha =1/sqrt(2); //dist(rnd_gen);
    double phi = 0.5*M_PI; //dist(rnd_gen)*2*M_PI;
    double V = sqrt(2*opt.E*1.602e-19/opt.m);
    double Vx = V*cos_alpha;
    double sin_alpha = sqrt(1-cos_alpha*cos_alpha);
    double Vy = V*sin(phi)*sin_alpha;
    double Vz = V*cos(phi)*sin_alpha;
    //Vz is not traced since in this direction slit is infinite
    double dt = opt.D/(Vy*100);
    double x = 0;
    double y = opt.D/2;
    std::vector<double> x_trace = {x};
    std::vector<double> y_trace = {y};

    auto F = [&](double cur_x){
        double field = -(potential(cur_x+1e-6)-potential(cur_x))/1e-6;
        return opt.Q*field;
    };

    std::vector<double> force;
    std::vector<double> energy = {0.5*opt.m*V*V};

    //LeapFrog
    double Vx_prev = Vx - 0.5*dt*F(x)/opt.m;
    while(x<opt.L && x>=0){
        force.push_back(F(x));
        Vx = Vx_prev + dt*F(x)/opt.m;
        x += Vx*dt;
        x_trace.push_back(x);
        y+= Vy*dt;
        if(y>opt.D){
            y = opt.D;
            Vy *= -1;
        }
        if(y<0){
            y = 0.0;
            Vy *=-1;
        }
        y_trace.push_back(y);
        Vx_prev = Vx;
        energy.push_back(0.5*opt.m*(Vx*Vx+Vy*Vy+Vz*Vz));
    }

    FILE* f_out = fopen(opt.out_file.c_str(), "w");
    if(!f_out){
        fprintf(stderr, "Cannot open file %s\n", opt.out_file.c_str());
        exit(EXIT_FAILURE);
    }
    fprintf(f_out, "#x, m\ty, m\tF, N\tE, J\n");
    for(size_t i=0; i<x_trace.size(); i++){
        fprintf(f_out, "%.10e\t%.10e\t%.10e\t%.10e\n", x_trace[i], y_trace[i], force[i], energy[i]);
    }
    fclose(f_out);

    return 0;
}
