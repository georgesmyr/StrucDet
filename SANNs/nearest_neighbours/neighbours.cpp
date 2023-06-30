#include "neighbours.hpp"


NbData::NbData(int identity, int from_identity, double dist)
: id{identity}, distance{dist}, from_id{from_identity}, theta{0}, phi{0}{
}

NbData::NbData(int identity, int from_identity, double dist, double theta_angle, double phi_angle)
: id{identity}, distance{dist}, from_id{from_identity}, theta{theta_angle}, phi{phi_angle}{
}

void NbData::theta_phi(const std::vector<std::vector<double>> &positions, const std::vector<double>& box){
    /*
     DESCRIPTION: Calculates the spherical coordinates in respect to the particle in reference.
     */
    std::vector<double> dist_vec = {0, 0, 0};
    for(int i = 0; i < 3; ++i){
        dist_vec[i] = positions[id][i] - positions[from_id][i];
        if(dist_vec[i] > box[i] / 2){
            dist_vec[i] -= box[i];
        }else if(dist_vec[i] < - box[i] / 2){
            dist_vec[i] += box[i];
        }
    }
    phi = std::atan2(dist_vec[1], dist_vec[0]);
    theta = std::acos(dist_vec[2] / distance);
}

void calculate_theta_phi(std::vector<NbData>& sanns,const std::vector<std::vector<double>>& positions, const std::vector<double>& box){
    /*
     DESCRIPTION: Calculates theta and phi angles for a vector of SANNs.
     */
    unsigned long size = sanns.size();
    for(int i = 0; i < size; ++i){
        sanns[i].theta_phi(positions, box);
    }
}

int NbLess(const NbData& nb1, const NbData& nb2){
    /*
     Order of neighbours.
     */
    return nb1.distance < nb2.distance;
}

void print_neighbours(const std::vector<NbData>& vec){
    /*
     Prints out the neighbours in the list.
     */
    unsigned long N = vec.size();
    for(auto i = 0; i < N; ++i){
        printf("(%d, %d, %f, %f, %f)\n", vec[i].from_id, vec[i].id, vec[i].distance, vec[i].theta, vec[i].phi);
    }
    printf("-----------------\n");
}


std::vector<NbData> getSANNs(int idx, std::vector<std::vector<double>>& dist_mat){
    /*
     DESCRIPTION: Finds the SANNs of particle idx.
     */
    
    //Create vector with neighbours
    std::vector<NbData> neighbours = {};
    unsigned long num_particles = dist_mat.size();
    unsigned long count = num_particles - 1;
    for(int i = 0; i < num_particles; ++i){
        if(i == idx){continue;}
        neighbours.push_back(NbData(i, idx, dist_mat[idx][i]));
    }
    
    //Sort neighbours based on distance
    sort(neighbours.begin(), neighbours.end(), NbLess);
    
    double distanceSum {0};
    // We define it here so that after loop is over, j = 3.
    int j = 0;
    for(j = 0; j < 3; ++j){
        distanceSum += neighbours[j].distance;
    }
    
    // Set SANN radius to distanceSum / (j - 2)
    double radius = distanceSum;
    
    while((j < count) and (radius > neighbours[j].distance)){
        
        distanceSum += neighbours[j].distance;
        radius = distanceSum / (j - 2);
        
        ++j;
    }
    if(j == count){
        return {{}};
    }
    std::vector<NbData> sanns = {};
    for(int k = 0; k < j; ++k){
        sanns.push_back(neighbours[k]);
    }
    return sanns;
}

void save_sanns_angles(std::string path1, std::string path2 , bool csv_ext, int idx, std::vector<NbData>& sanns){
    /*
     DESCRIPTION: Saves SANNs, i.e. the id of the neighbours, the distance,
     and the phi and theta angles of the distance vector.
     */
    std::string folder {"dats/"};
    std::string extension {".dat"};
    std::string delimeter {" "};
    if(csv_ext){
        folder = "csvs/";
        extension = ".csv";
        delimeter = ",";
    }
    std::string path = path1 + folder + path2;
    
    std::ofstream outfile(path + "sanns_" + std::to_string(idx) + extension);
    if(!outfile.is_open()){
        std::cerr << "Error: File not found when saving sanns angles." << std::endl;
        return;
    }
    if(csv_ext){
        outfile << "id,distance,theta,phi" << std::endl;
    }
    for(int i = 0; i < sanns.size(); ++i){
        outfile << sanns[i].id << delimeter << sanns[i].distance << delimeter << sanns[i].theta << delimeter << sanns[i].phi << std::endl;
    }
    outfile.close();
    
}

std::vector<NbData> read_sanns(std::string path1, std::string path2, bool csv_ext, int idx){
    /*
     DESCRIPTION: Reads the SANNs from a file, i.e. id, distance, theta and phi angles.
     */
    
    std::string folder {"dats/"};
    std::string extension {".dat"};
    std::string delimeter {" "};
    if(csv_ext){
        folder = "csvs/";
        extension = ".csv";
        delimeter = ",";
    }
    std::string path = path1 + folder + path2;

    std::ifstream infile(path + "sanns_" + std::to_string(idx) + extension);
    if(!infile.is_open()){
        std::cerr << "Error: File not found when reading sanns." << std::endl;
    }
    
    std::vector<NbData> sanns = {};
    int id = 0;
    double distance{0}, theta{0}, phi{0};
    
    std::string line {""};
    while(std::getline(infile, line)){
        std::istringstream iss(line);
        iss >> id;
        iss >> distance >> theta >> phi;
//        printf("(%d, %f, %f, %f)\n", id, distance, theta, phi);
        NbData nb = NbData(id, idx, distance, theta, phi);
        sanns.push_back(nb);
    }
    infile.close();
    return sanns;
}


void symmetrize_neighbours(std::string path1, std::string path2, int idx){
    /*
     DESCRIPTION: Symmetrises neighbours, i.e. if j is in the neighbours of i,
     and i is not in the neighbours of j, we remove j from the neighbours of i.
     */
    
    //LOAD SANNS FOR PARTICLE idx
    std::vector<NbData> sanns = read_sanns(path1, path2, false, idx);
    std::vector<NbData> sanns_new = {};
    
    //EXTRACT ids OF NEIGHBOURS OF PARTICLE idx
    std::vector<int> sanns_ids = {};
    for(const NbData& nb : sanns){
        sanns_ids.push_back(nb.id);
    }

    //FOR EACH NEIGHBOUR i ITS NEIGHBOURS. IF idx IS NOT A NEIGHBOUR OF i,
    //DELETE i from THE NEIGHBOURS OF idx.
    for(int i : sanns_ids){
        
        //Load the neighbours of particle i
        std::vector<NbData> sanns_i = read_sanns(path1, path2, false, i);
        
        //Look for idx in the neighbours of i
        bool found = false;
        for(int j = 0; j < sanns_i.size(); ++j){
            //If found, stop the loop and turn the flag found to true.
            if(sanns_i[j].id == idx){
                found = true;
                break;
            }
        }

        //If idx is not found in the neighbours of i, remove i from the neighbours of idx.
        //i.e. save all the neighbours of idx except i.
        if(!found){
            sanns_new = {};
            for(const NbData& nb : sanns){
                if(nb.id != i){
                    sanns_new.push_back(nb);
                }
            }
            
            sanns = {};
            for(NbData& nb : sanns_new){
                sanns.push_back(nb);
            }
        }
    }
    
    save_sanns_angles(path1, path2, true, idx, sanns);
}

void save_sanss_positions(std::string path, int idx, std::vector<NbData>& sanns,
                          std::vector<std::vector<double>>& positions, std::vector<double>& box){
    /*
     DESCRIPTION: saves the positions of all particles, and colours the SANNs of the particle idx
     */
    std::vector<int> sanns_ids = {};
    for(int i = 0; i < sanns.size(); ++i){
        sanns_ids.push_back(sanns[i].id);
    }
    
    std::vector<std::vector<double>> colored_pos = positions;
    for(int i = 0; i < colored_pos.size(); ++i){
        auto it = std::find(sanns_ids.begin(), sanns_ids.end(), i);
        if (it == sanns_ids.end()){
            colored_pos[i].push_back(0.11);
            colored_pos[i].push_back(0.249);
            colored_pos[i].push_back(0.249);
        }
        else{
            colored_pos[i].push_back(1);
            colored_pos[i].push_back(0);
            colored_pos[i].push_back(0);
        }
                     
    }
    colored_pos[idx][4] = 0;
    colored_pos[idx][5] = 1;
    colored_pos[idx][6] = 0;
    
    std::ofstream outfile(path + "sanns_pos" + std::to_string(idx) + ".dat");
    if(!outfile.is_open()){
        std::cerr << "Error: file not found." << std:: endl;
        return;
    }
    outfile << box[3] << std::endl;
    for(int d = 0; d < box.size() - 1; ++d){
        outfile << box[d] << std::endl;
    }
    for(int i = 0; i < colored_pos.size(); ++i){
        for(int j = 0; j < colored_pos[0].size(); ++j){
            outfile << colored_pos[i][j] << " ";
        }
        outfile << std::endl;
    }

    outfile.close();
}

