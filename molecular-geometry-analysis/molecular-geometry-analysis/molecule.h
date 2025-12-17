//
//  molecule.h
//  molecular-geometry-analysis
//
//  Created by Ethan Song on 12/16/25.
//

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace std;

class Molecule {
    
public:
    int natoms;
    int charge;
    int multiplicity;
    vector<int> zvals;
    vector<double> coords;
    string point_group;
    vector<double> distance_matrix;
    vector<double> angle_matrix;
    void print_geom();
    void compute_distances();
    void print_distance_matrix();
    void compute_angles();
    Molecule(const char* filename, int q, int s);
private:
    inline int idx(int i, int j) const noexcept {
        return i * natoms + j;
    }
    inline int idx3(int i, int j, int k) const noexcept {
        return (i * natoms + j) * natoms + k;
    }
};

void Molecule::print_geom() {
    cout << "Geometry\nNumber of atoms: " << natoms << endl;
    for (int i = 0; i < natoms; i++) {
        printf("%d %6.3f %6.3f %6.3f\n", zvals[i], coords[idx(i,0)], coords[idx(i,1)], coords[idx(i,2)]);
    }
}

void Molecule::compute_distances() {
    distance_matrix.resize(natoms * natoms);
    for (int i = 0; i < natoms; i++) {
        distance_matrix[idx(i,i)] = 0.;
        for (int j = i + 1; j < natoms; j++) {
            double dist = sqrt(pow(coords[idx(i,0)] - coords[idx(j,0)], 2) + pow(coords[idx(i,1)] - coords[idx(j,1)], 2) + pow(coords[idx(i,2)] - coords[idx(j,2)], 2));
            distance_matrix[idx(i,j)] = dist;
            distance_matrix[idx(j,i)] = dist;
        }
    }
}

void Molecule::print_distance_matrix() {
    cout << "Distance matrix\n";
    for (int i = 0; i < natoms; i++) {
        for (int j = 0; j < natoms; j++) {
            printf("%.3f ", distance_matrix[idx(i,j)]);
        }
        cout << endl;
    }
}

void Molecule::compute_angles() {
    angle_matrix.resize(pow(natoms,3));
}

Molecule::Molecule(const char* filename, int q, int s) {
    charge = q;
    multiplicity = s;
    ifstream input(filename);
    assert(input.good());
    input >> natoms;
    zvals.resize(natoms);
    coords.resize(natoms * natoms);
    for (int i = 0; i < natoms; i++) {
        input >> zvals[i] >> coords[idx(i,0)] >> coords[idx(i,1)] >> coords[idx(i,2)];
    }
    input.close();
}
