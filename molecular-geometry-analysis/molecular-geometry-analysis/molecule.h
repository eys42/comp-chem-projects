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

using namespace std;

class Molecule {
    
public:
    int natoms;
    int charge;
    int multiplicity;
    int *zvals;
    double **coords;
    string point_group;
    double **distance_matrix;
    void print_geom();
    void compute_distances();
    void print_distance_matrix();
    Molecule(const char* filename, int q, int s);
    ~Molecule();
};

void Molecule::print_geom() {
    cout << "Geometry\nNumber of atoms: " << natoms << endl;
    for (int i = 0; i < natoms; i++) {
        printf("%d %6.3f %6.3f %6.3f\n", zvals[i], coords[i][0], coords[i][1], coords[i][2]);
    }
}

void Molecule::compute_distances() {
    distance_matrix = new double*[natoms];
    for (int i = 0; i < natoms; i++) {
        distance_matrix[i] = new double[natoms];
    }
    for (int i = 0; i < natoms; i++) {
        distance_matrix[i][i] = 0.;
        for (int j = i + 1; j < natoms; j++) {
            double dist = sqrt(pow(coords[i][0] - coords[j][0], 2) + pow(coords[i][1] - coords[j][1], 2) + pow(coords[i][2] - coords[j][2], 2));
            distance_matrix[i][j] = dist;
            distance_matrix[j][i] = dist;
        }
    }
}

void Molecule::print_distance_matrix() {
    cout << "Distance matrix\n";
    for (int i = 0; i < natoms; i++) {
        for (int j = 0; j < natoms; j++) {
            printf("%.3f ", distance_matrix[i][j]);
        }
        cout << endl;
    }
}

Molecule::Molecule(const char* filename, int q, int s) {
    charge = q;
    multiplicity = s;
    ifstream input(filename);
    assert(input.good());
    input >> natoms;
    zvals = new int[natoms];
    coords = new double*[natoms];
    for (int i = 0; i < natoms; i++) {
        coords[i] = new double[3];
    }
    for (int i = 0; i < natoms; i++) {
        input >> zvals[i] >> coords[i][0] >> coords[i][1] >> coords[i][2];
    }
    input.close();
}

Molecule::~Molecule() {
    delete[] zvals;
    for (int i = 0; i < natoms; i++) {
        delete[] coords[i];
        delete[] distance_matrix[i];
    }
    delete[] coords; delete[] distance_matrix;
}
