#ifndef distances_hpp
#define distances_hpp

#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

void read_positions(std::vector<std::vector<double>>&, std::vector<double>&, std::string, bool);
double distance(std::vector<double>&, std::vector<double>&, std::vector<double>&);

std::vector<std::vector<double>> distance_matrix(std::vector<std::vector<double>>&, std::vector<double>&);
void print_distance_matrix(const std::vector<std::vector<double>>&);
void write_distance_matrix(std::vector<std::vector<double>>&, std::string);
std::vector<std::vector<double>> read_distance_matrix(std::string, unsigned long);

std::vector<std::vector<double>> distance_matrix_from_positions(std::string, std::string);

#endif /* distances_hpp */
