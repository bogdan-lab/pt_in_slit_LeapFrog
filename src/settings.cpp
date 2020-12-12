#include "settings.hpp"
#include "particle.hpp"

void Settings::SaveSettings(FILE* savefile) const {
    fprintf(savefile, "#Slit D[m]             = %.6e\n", D);
    fprintf(savefile, "#Slit L[m]             = %.6e\n", L);
    fprintf(savefile, "#Particle m[kg]        = %.6e\n", m);
    fprintf(savefile, "#Particle Q[C]         = %.6e\n", Q);
    fprintf(savefile, "#Slit abs(phi_val[V])  = %.6e\n", phi_val);
    fprintf(savefile, "#Slit ReflectionCoeff  = %.6e\n", R);
    fprintf(savefile, "#Particle Number       = %zu\n", pt_num);
    fprintf(savefile, "#Particle energy dist  = %s\n", energy_fname.c_str());
    fprintf(savefile, "#Slit potential shape  = %s\n", phi_fname.c_str());
}

std::ifstream load_helper(const std::string& filename){
    std::ifstream file(filename);
    if(!file.is_open()){
        fprintf(stderr, "Cannot open file %s \n", filename.c_str());
        exit(EXIT_FAILURE);
    }
    return file;
}


void SaveResult(const std::string& filename, const std::vector<PtStat>& stat,
                const Settings& opt){
    FILE* f_out = fopen(filename.c_str(), "w");
    if(!f_out){
        fprintf(stderr, "Cannot open file %s\n", filename.c_str());
        exit(EXIT_FAILURE);
    }
    opt.SaveSettings(f_out);
    fprintf(f_out, "#Vx, m/s\tVy, m/s\tVz, m/s\tE, eV\tsurface collisions\n");
    for(size_t i=0; i<stat.size(); i++){
        fprintf(f_out, "%.10e\t%.10e\t%.10e\t%.10e\t%zu\n",
         stat[i].Vx, stat[i].Vy, stat[i].Vz, stat[i].E, stat[i].surf_col_count);
    }
    fclose(f_out);
}

