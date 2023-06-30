#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <complex>
#include "distances.hpp"
#include "neighbours.hpp"

void print_positions(const std::vector<std::vector<double>>& positions){
    unsigned long size = positions.size();
    unsigned long length = positions[0].size();
    for(int i = 0; i < size; ++i){
        std::cout << i << ") ";
        for(int j = 0; j < length; ++j){
            std::cout << positions[i][j] << " ";
        }std::cout << std::endl;
    }
}

std::string PHASE = "fluid";
int num = 1;
int step = 50;

int main(int argc, const char * argv[]) {
    
    std::string phase_dir {PHASE + "s/" + PHASE + std::to_string(num) + "/"};
    std::string phase_step_dir {phase_dir + "step" + std::to_string(step) + "/"};
    
    std::string data_path {"/Users/georgesmyridis/Desktop/Trading/ModSim/StructDet/data/"};
    std::string positions_path {data_path + "positions/" + phase_dir + "pos" + std::to_string(step) + ".dat"};
    std::string dist_mat_path {data_path + "distances/" + phase_step_dir};
    std::string sanns_pos_path {data_path + "sanns_pos/" + phase_step_dir};
    std::string sanns_angles_path {data_path + "sanns_data/"};

    //LOAD POSITIONS AND BOX PARAMETERS
    std::vector<double> box = {};
    std::vector<std::vector<double>> positions = {};
    std::cout << positions_path << std::endl;
    read_positions(positions, box, positions_path, false);
    unsigned long num_particles = positions.size();

    //CALCULATE DISTANCE MATRIX
    std::vector<std::vector<double>> dist_mat = distance_matrix(positions, box);

    //GET NEAREST NEIGHBOURS
    for(int idx = 0; idx < positions.size(); ++idx){
        std::vector<NbData> sanns = getSANNs(idx, dist_mat);
        calculate_theta_phi(sanns, positions, box);
        save_sanss_positions(sanns_pos_path, idx, sanns, positions, box);
        save_sanns_angles(sanns_angles_path, phase_step_dir, false, idx, sanns);
    }

    for(int i = 0; i < num_particles; ++i){
        symmetrize_neighbours(sanns_angles_path, phase_step_dir, i);
    }

    
    

    
    return 0;
}


