#include "distances.hpp"

void read_positions(std::vector<std::vector<double>>& positions, std::vector<double>& box, std::string filepath, bool print){
    //IMPORT POSITIONS AND BOX SPECS
    std::ifstream infile(filepath);
    if(!infile.is_open()){
        std::cerr << "Error: Unable to open file." << std::endl;
        return;
    }
    double N{0}, Lx{0}, Ly{0}, Lz{0};
    double x{0}, y{0}, z{0}, d{0};
    infile >> N >> Lx >> Ly >> Lz;
    box = {Lx, Ly, Lz, N};
    for(int i = 0; i < N; ++i){
        infile >> x >> y >> z >> d;
        positions.push_back({x,y,z,d});
    }
    
    if(print){
        //PRINT BOX
        for(int i = 0; i < 4; ++i){
            std::cout << box[i] << " ";
        }std::cout << std::endl;
        //PRINT POSITIONS
        for(int i = 0; i < box[3]; ++i){
            std::cout << i << ") ";
            for(int j = 0; j < 4; ++j){
                std::cout << positions[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }
    
}

double distance(std::vector<double>& vector1, std::vector<double>& vector2, std::vector<double>& box){
    /*
     Calculates distance between two particles taking into account the periodic boundary conditions.
     */
    double minD {0}, dist2 {0}, factor {0};
    unsigned long ndim {vector1.size() - 1};
    for (auto i {0}; i < ndim; ++i){
        minD = vector1[i] - vector2[i];
        factor = (int)(2.0 * minD /  box[i]);
        minD -= factor * box[i];
        dist2 += minD * minD;
    }
    return std::sqrt(dist2);
}

void print_distance_matrix(const std::vector<std::vector<double>>& dist_mat){
    /*
     Prints the distance matrix.
     */
    unsigned long size = dist_mat.size();
    for(int i = 0; i < size; ++i){
        for(int j = 0; j < size; ++j){
            std::cout << dist_mat[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

std::vector<std::vector<double>> distance_matrix(std::vector<std::vector<double>>& positions, std::vector<double>& box){
    /*
     /Calculates the distance matrix for the particles. i.e. it calculates all the distances among the particles.
     */
    
    std::vector<std::vector<double>> dist_mat = {};
    unsigned long num_particles {positions.size()};
    
    for(int particle1 = 0; particle1 < num_particles; ++particle1){
        std::vector<double> dist_from_particle1 = {};
        for(int particle2 = 0; particle2 < num_particles; ++particle2){
            double dist = distance(positions[particle1], positions[particle2], box);
            dist_from_particle1.push_back(dist);
        }
        dist_mat.push_back(dist_from_particle1);
    }
    return dist_mat;
}

void write_distance_matrix(std::vector<std::vector<double>>& dist_mat, std::string filepath){
    /*
     Saves the distance matrix in a dat file with filepath pathname.
     */
    
    std::ofstream outfile(filepath);
    if (!outfile.is_open()) {
        std::cerr << "Error: Unable to open output file" << std::endl;
        return;
    }
    
    unsigned long num_particles = dist_mat.size();
    for(int i = 0; i < num_particles; ++i){
        for(int j = 0; j < num_particles; ++j){
            outfile << dist_mat[i][j] << " ";
        }
        outfile << std::endl;
    }
}

std::vector<std::vector<double>> read_distance_matrix(std::string filepath, unsigned long size){
    /*
     Reads distance matrix from dat file with pathname = filepath.
     */
    
    std::ifstream infile(filepath);
    if(!infile.is_open()){
        std::cerr << "Error: could not open file." << std::endl;
        return {{0}};
    }
    
    std::vector<std::vector<double>> dist_mat = {};
    double dist {0};
    for(int i = 0; i < size; ++i){
        std::vector<double> row = {};
        for(int j = 0; j < size; ++j){
            infile >> dist;
            row.push_back(dist);
        }
        dist_mat.push_back(row);
    }
    return dist_mat;
}


std::vector<std::vector<double>> distance_matrix_from_positions(std::string positions_path, std::string dist_mat_path){
    
    //Load positions and box parameters.
    std::vector<std::vector<double>> positions {};
    std::vector<double> box = {};
    read_positions(positions, box, positions_path, false);
    unsigned long num_particles = positions.size();
    
    //Calculate & Save the distance matrix.
    std::vector<std::vector<double>> dist_mat = distance_matrix(positions, box);
    write_distance_matrix(dist_mat, dist_mat_path);
    
    //Read the distance matrix.
    std::vector<std::vector<double>> loaded_dist_mat = read_distance_matrix(dist_mat_path, num_particles);
    
    //Check if the loaded is the same as the calculated.
    for(int i = 0; i < num_particles; ++i){
        for(int j = 0; j < num_particles; ++j){
            assert(dist_mat[i][j] == loaded_dist_mat[i][j]);
        }
    }
    
    return dist_mat;
}


