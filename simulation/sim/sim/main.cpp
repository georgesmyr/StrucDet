#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <ctime> // for time()
#include <string>
#include "mt19937.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;

const string path {"/Users/georgesmyridis/Desktop/Trading/ModSim/StructDet/data/positions/fccs/fcc2/"};

// FCC PARAMETERS
const int NDIM = 3;
int Nx {8}, Ny {8}, Nz {8}; //change also POSITIONS size!!!
int N {0};

double RADIUS = 0.5;
double DIAMETER = 2 * RADIUS;
double LATTICE_SPACING = 4.0 * RADIUS;

//SYSTEM VARIABLES
double POSITIONS[4*7*7*7][4] = {{0}};
double BOX[4] = {0};

//POTENTIAL PARAMETERS
const double THRESH1 {DIAMETER};
const double LAMBDA {1.5};
const double THRESH2 {LAMBDA * THRESH1};
const double POTENTIAL {-1};

//MONTE CARLO SIMULATION PARAMETERS
int STEPS {300000};
double BETA {1};
const double DISPLACEMENT_FRACTION {0.025};
double DELTA_MAX {DISPLACEMENT_FRACTION * LATTICE_SPACING};
double NUMBER_DENSITY {1.01};
const unsigned int OUTPUT_STEPS {5000};


void generate_fcc(){
    /*
     Description:
     ------------
     This function creates a configuration of hard spheres on a face-centered cubic crystal.
     
     Parameteres:
     ------------
     Nx, Ny, Nz: The number of hard spheres in each direction.
     d: The diameter of the hard spheres.
     spcacing: The lattice spacing in units of the diameter.
     
     Output:
     -------
     xy.dat file in the directory (path). The first line has the total number of particles in the box. The
     second line has the coordinates of the box corner in x-direction, and the next for the y- and z- directions
     respectively. The rest of the lines have the x-, y-, and z- coordinates of the hard sheres, and their diameter.
    */
    cout << "Generating FCC..." << endl;
    ofstream outfile(path + "pos.dat");
    if (!outfile.is_open()) {
        cerr << "Error: Unable to open output file" << endl;
        return;
    }

    if (NDIM == 2){
        N = 2 * Nx * Ny;
        double Lx {Nx * LATTICE_SPACING};
        double Ly {Ny * LATTICE_SPACING};

        outfile << N << endl;
        outfile << Lx << endl;
        outfile << Ly << endl;
        outfile << 0 << endl;
            
        for (int i {0}; i < Nx; i++){
            for (int j {0}; j < Ny; j++){
                outfile << i * LATTICE_SPACING << " " << j * LATTICE_SPACING << " " << 0 << " " << DIAMETER << endl;
                outfile << (i + 0.5) * LATTICE_SPACING << " " << (j + 0.5) * LATTICE_SPACING << " " << 0 << " " << DIAMETER << endl;
            }
        }

    }
    else if(NDIM == 3){
        cout << NDIM << endl;
        N = 4 * (Nx - 1) * (Ny - 1) * (Nz - 1);
        cout << N << endl;
        double Lx {(Nx - 1) * LATTICE_SPACING}, Ly {(Ny - 1) * LATTICE_SPACING}, Lz {(Nz - 1) * LATTICE_SPACING};
        
        outfile << N << endl;
        outfile << Lx << endl;
        outfile << Ly << endl;
        outfile << Lz << endl;
        
        for(auto i = 0; i < Nx - 1; ++i){
            for(auto j = 0; j < Ny - 1; ++j){
                for(auto k = 0; k < Nz - 1; ++k){
                    outfile << i * LATTICE_SPACING << " " << j * LATTICE_SPACING << " " << k * LATTICE_SPACING << " " << DIAMETER << endl;
                    outfile << (i + 0.5) * LATTICE_SPACING << " " << (j + 0.5) * LATTICE_SPACING << " " << k * LATTICE_SPACING << " " << DIAMETER << endl;
                    outfile << (i + 0.5) * LATTICE_SPACING << " " << j * LATTICE_SPACING << " " << (k + 0.5) * LATTICE_SPACING << " " << DIAMETER << endl;
                    outfile << i * LATTICE_SPACING << " " << (j + 0.5) * LATTICE_SPACING << " " << (k + 0.5) * LATTICE_SPACING << " " << DIAMETER << endl;
                }
            }
        }
    }
    
    cout << "FCC generated." << endl;
    outfile.close();
}

void read_data(){
    /*
     DESCRIPTION: Reads all the initial positions from a dat file.
     */
    
    ifstream infile;
    infile.open(path + "pos.dat");
    if (not infile.is_open()){
        cerr << "Could not read the file." << endl;
        return;
    }
    infile >> BOX[3] >> BOX[0] >> BOX[1] >> BOX[2];
    for(int i = 0; i < N; ++i){
        for(int j = 0; j < 4; ++j){
            infile >> POSITIONS[i][j];
        }
    }
    N = BOX[3];
}


void write_data(unsigned int num){
    /*
     DESCRIPTION: This function writes the positions of the particles in a dat file named xy(num).dat.
     */

    ofstream outfile(path + "pos" + to_string(num) + ".dat");
    if (!outfile.is_open()){
        cerr << "Error: Unable to open output file." << endl;
        return;
    }
        outfile << BOX[3] << endl;
    for(int i = 0; i < 3; ++i){
        outfile << BOX[i] << endl;
    }
    for(auto i = 0; i < N; ++i){
        for(int j = 0; j < 4; ++j){
            outfile << POSITIONS[i][j] << " ";
        }
        outfile << endl;
    }
    outfile.close();
}


void print_box(){
    for(int i = 0; i < 4; ++i){
        cout << BOX[i] << " ";
    }
    cout << endl;
}


void print_positions(){
    for(int i = 0; i < N; ++i){
        cout << i << ") ";
        for(int j = 0; j < 4; ++j){
            cout << POSITIONS[i][j] << " ";
        }
        cout << endl;
    }
}


double distance(unsigned int index1, unsigned int index2){
    /*
     DESCRIPTION: Returns the distance between the centers of two particles, given the boundary conditions of the box.
     INPUT: index1, index2: The indices of the two particles.
     OUTPUT: The distance.
     */

    double minD {0}, dist2 {0}, factor {0};
    for (auto i {0}; i < NDIM; ++i){
        minD = POSITIONS[index1][i] - POSITIONS[index2][i];
        factor = (int)(2.0 * minD /  BOX[i]);
        minD -= factor * BOX[i];
        dist2 += minD * minD;
    }
    
    return sqrt(dist2);
}


double potential(unsigned int idx1, unsigned int idx2){
    /*
     DESCRIPTION: Returns the potential energy between the particles with idx1 and idx2.
     */
    double dist {distance(idx1, idx2)};
    if(NDIM == 2){
        if(dist < THRESH1){
            return INFINITY;
        }
        else if(dist < THRESH2){
            return POTENTIAL;
        }
        else{
            return 0;
        }
    }else if(NDIM == 3){
        if (dist < THRESH1){
            return INFINITY;
        }else{
            return 0;
        }
    }
}


double energy_contr(unsigned int idx){
    /*
     DESCRIPTION: Calculates the energy contribution for the particle idx.
     */
    double contribution {0};
    for(auto i = 0; i < N; ++i){
        if(i != idx){
            contribution += potential(idx, i);
        }
        if(contribution == INFINITY){
            return contribution;
        }
    }
    return contribution;
}

int check_overlaps(unsigned int idx){
    /*
     DESCRIPTION: Checks overlaps of particle idx with others.
     */
    
    double dist {0};
    for(auto i = 0; i < N; ++i){
        if(i != idx){
            dist = distance(idx, i);
            if(dist < THRESH1){
                return 1;
            }
        }
    }
    
    return 0;
}



int move_particle(){
    
    
    //CHOOSE RANDOM PARTICLE
    int idx = (int)(dsfmt_genrand() * N);
    double displacements[NDIM] = {0};
    
    //CALCULATE CURRENT ENERGY CONTRIBUTION FRO MPARTICLE idx
    double contr_before {energy_contr(idx)};
    
    //CALCULATE RANDOM DISPLACEMENTS AND CHANGE POSITION
    for(auto i {0}; i < NDIM; ++i){
        displacements[i] = (2 * dsfmt_genrand() - 1) * DELTA_MAX;
        POSITIONS[idx][i] += displacements[i];
    
        //Account for periodic boundary conditions
        if(POSITIONS[idx][i] < 0){
            POSITIONS[idx][i] += BOX[i];
        }else if(POSITIONS[idx][i] > BOX[i]){
            POSITIONS[idx][i] -= BOX[i];
        }
        assert(POSITIONS[idx][i] >= 0 and POSITIONS[idx][i] <= BOX[i]);
    }
    if(NDIM == 2){
        //CALCULATE CHANGE IN ENERGY
        double contr_after {energy_contr(idx)};
        double dE {contr_after - contr_before};
            
        //ACCEPT OR REJECT MOVE
        //Accept and return 1
        if(dE != INFINITY and (dE <= 0 or dsfmt_genrand() < exp(- BETA * dE))){
            return 1;
        }
    }else if(NDIM == 3){
        if(check_overlaps(idx) == 0){
            return 1;
        }
    }
    

    //Reject, reverse the displacements and return 0.
    for(int i = 0; i < NDIM; ++i){
        POSITIONS[idx][i] -= displacements[i];
        
        //Account for periodic boundary conditions
        if(POSITIONS[idx][i] < 0){
            POSITIONS[idx][i] += BOX[i];
        }else if(POSITIONS[idx][i] > BOX[i]){
            POSITIONS[idx][i] -= BOX[i];
        }
        assert(POSITIONS[idx][i] >= 0 and POSITIONS[idx][i] <= BOX[i]);
    }
    
    return 0;
}

void ensure_overlaps(){
    /*
     DESCRIPTION: The function checks if there are overlaps between the particles.
     */
    for(int i = 0; i < N - 1; ++i){
        for(int j = i + 1; j < N; ++j){
            assert(distance(i,j) > 2 * RADIUS);
        }
    }
}

void set_number_density(double number_density){
    /*
     DESCRIPTION: Set number density.
     */
    double volume {BOX[0] * BOX[1] * BOX[2]};
    cout << "Volume: " << volume << endl;
    double scale_factor = pow(N * pow(DIAMETER, 3) / (number_density * volume), 0.333);
    cout << "Scale factor: " << scale_factor << endl;
    for (int i {0}; i < BOX[3]; i++){
        for (int j{0}; j < NDIM; j++){
            POSITIONS[i][j] *= scale_factor;
        }
    }
    print_box();
    for (int i{0}; i < NDIM; i++){
        BOX[i] *= scale_factor;
    }
    print_box();
//    LATTICE_SPACING *= scale_factor;
//    DELTA_MAX = DISPLACEMENT_FRACTION * LATTICE_SPACING;
    
    cout << "Number density: " << N / (BOX[0] * BOX[1] * BOX[2]) << endl;
}

int main(int argc, const char * argv[]) {
    /*
     You can switch from hard spheres in 3D to square shoulder particles in 2D
     by changing the dimension from NDIM = 3 to NDIM = 2.
     IMPORTANT!: The number of particles changes, and thus you have to change
     the dimensions of the array POSITIONS accordingly.
     */
    
    dsfmt_seed(time(NULL));
    //GENERATE FCC
    generate_fcc();
    read_data();
    write_data(0);
    set_number_density(NUMBER_DENSITY);
    ensure_overlaps();
//    print_positions();
//    print_box();

    float accepted {0};
    for(auto step = 0; step < STEPS; ++step){
        accepted += move_particle();
        if(step % OUTPUT_STEPS == 0){
            ensure_overlaps();
            printf("Acceptance ratio: %f\n", accepted / OUTPUT_STEPS);
            write_data(step / OUTPUT_STEPS);
            accepted = 0;
        }
    }

//    print_positions();
    
    return 0;
}




