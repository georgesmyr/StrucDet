//
//  neighbours.hpp
//  nearest_neighbours
//
//  Created by Γιώργος Σμυρίδης on 10/6/23.
//

#ifndef neighbours_hpp
#define neighbours_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>

class NbData{
public:
    
    int from_id;
    int id;
    double distance;
    double theta;
    double phi;
        
    void set_theta(double);
    void set_phi(double);
    
    void theta_phi(const std::vector<std::vector<double>>&, const std::vector<double>&);
        
    //Constructor
    NbData(int, int, double);
    NbData(int, int, double, double, double);

};

int NbLess(const NbData&, const NbData&);

void calculate_theta_phi(std::vector<NbData>&,const std::vector<std::vector<double>>&, const std::vector<double>&);
void print_neighbours(const std::vector<NbData>&);

std::vector<NbData> getSANNs(int, std::vector<std::vector<double>>&);

void save_sanns_angles(std::string, std::string, bool, int, std::vector<NbData>&);
std::vector<NbData> read_sanns(std::string, std::string, bool, int);
void save_sanss_positions(std::string, int, std::vector<NbData>&, std::vector<std::vector<double>>&, std::vector<double>&);


void symmetrize_neighbours(std::string, std::string, int);
void collect_sanns_ids(std::string, std::string);


#endif /* neighbours_hpp */
