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
    std::string out_file;
};

int main(){

    Settings opt;
    opt.D = 10e-6;
    opt.L = 0.04e-2;
    opt.m = 9.1e-31;
    opt.Q = 1.602e-19;
    opt.E = 10.0;
    opt.out_file = "test_trace.txt";

    std::mt19937 rnd_gen(42);
    std::uniform_real_distribution<double> dist(0, 0.5*M_PI);

    double alpha = dist(rnd_gen);
    double V = sqrt(2*opt.E*opt.Q/opt.m);
    double Vx = V*cos(alpha);
    double Vy = V*sin(alpha);
    double x = 0;
    double y = opt.D/2;
    std::vector<double> x_trace = {x};
    std::vector<double> y_trace = {y};
    while(x<opt.L && x>=0){
        double t_wall_coll = [&](){
            if(Vy>0) return (opt.D-y)/Vy;
            return -y/Vy;
        }();
        x += Vx*t_wall_coll;
        x_trace.push_back(x);
        y+= Vy*t_wall_coll;
        y_trace.push_back(y);
        Vy *= -1;
    }

    FILE* f_out = fopen(opt.out_file.c_str(), "w");
    if(!f_out){
        fprintf(stderr, "Cannot open file %s\n", opt.out_file.c_str());
        exit(EXIT_FAILURE);
    }
    fprintf(f_out, "#x, m\ty, m\n");
    for(size_t i=0; i<x_trace.size(); i++){
        fprintf(f_out, "%.10e\t%.10e\n", x_trace[i], y_trace[i]);
    }
    fclose(f_out);

    return 0;
}
