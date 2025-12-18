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
#include <numbers>

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
    vector<double> eX, eY, eZ;
    void print_geom();
    void compute_distances();
    void print_distance_matrix();
    void compute_angles();
    void print_angle_matrix();
    Molecule(const char* filename, int q, int s);
private:
    inline int idx(int i, int j) const noexcept {
        return i * natoms + j;
    }
    inline int cidx(int i, int j) const noexcept {
        return i * 3 + j;
    }
    inline int idx3(int i, int j, int k) const noexcept {
        return (i * natoms + j) * natoms + k;
    }
    void compute_unit_vectors();
    inline double dot_product(int i1, int j1, int i2, int j2) const;
};

void Molecule::print_geom() {
    cout << "Geometry\nNumber of atoms: " << natoms << endl;
    for (int i = 0; i < natoms; i++) {
        printf("%d %6.3f %6.3f %6.3f\n", zvals[i], coords[cidx(i,0)], coords[cidx(i,1)], coords[cidx(i,2)]);
    }
}

void Molecule::compute_distances() {
    distance_matrix.resize(natoms * natoms);
    for (int i = 0; i < natoms; i++) {
        distance_matrix[idx(i,i)] = 0.;
        for (int j = i + 1; j < natoms; j++) {
            double dist = sqrt(pow(coords[cidx(i,0)] - coords[cidx(j,0)], 2) + pow(coords[cidx(i,1)] - coords[cidx(j,1)], 2) + pow(coords[cidx(i,2)] - coords[cidx(j,2)], 2));
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

void Molecule::compute_unit_vectors() {
    eX.resize(natoms * natoms);
    eY.resize(natoms * natoms);
    eZ.resize(natoms * natoms);
    for (int i = 0; i < natoms; i++) {
        eX[idx(i,i)] = 0.;
        eY[idx(i,i)] = 0.;
        eZ[idx(i,i)] = 0.;
        for (int j = i + 1; j < natoms; j++) {
            double x = (coords[cidx(j,0)] - coords[cidx(i,0)]) / distance_matrix[idx(i,j)];
            double y = (coords[cidx(j,1)] - coords[cidx(i,1)]) / distance_matrix[idx(i,j)];
            double z = (coords[cidx(j,2)] - coords[cidx(i,2)]) / distance_matrix[idx(i,j)];
            eX[idx(i,j)] = x;
            eX[idx(j,i)] = -x;
            eY[idx(i,j)] = y;
            eY[idx(j,i)] = -y;
            eZ[idx(i,j)] = z;
            eZ[idx(j,i)] = -z;
        }
    }
}

inline double Molecule::dot_product(int i1, int j1, int i2, int j2) const {
    return eX[idx(i1,j1)] * eX[idx(i2,j2)] + eY[idx(i1,j1)] * eY[idx(i2,j2)] + eZ[idx(i1,j1)] * eZ[idx(i2,j2)];
}

void Molecule::compute_angles() {
    angle_matrix.resize(natoms * natoms * natoms);
    compute_unit_vectors();
    for (int i = 0; i < natoms; i++) {
        for (int j = 0; j < natoms; j++) {
            angle_matrix[idx3(i,j,j)] = 0.;
            for (int k = j + 1; k < natoms; k++) {
                if (i != j && i != k) {
                    double angle = acos(dot_product(j, i, j, k)) * 180 / numbers::pi;
                    angle_matrix[idx3(i,j,k)] = angle;
                    angle_matrix[idx3(k,j,i)] = angle;
                }
            }
        }
    }
}

void Molecule::print_angle_matrix() {
    cout << "Angle matrix\n";
    for (int i = 0; i < natoms; i++) {
        cout << "i = " << i << ":" << endl;
        for (int j = 0; j < natoms; j++) {
            for (int k = 0; k < natoms; k++) {
                printf("%6.2f ", angle_matrix[idx3(i,j,k)]);
            }
            cout << endl;
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
    zvals.resize(natoms);
    coords.resize(natoms * 3);
    for (int i = 0; i < natoms; i++) {
        input >> zvals[i] >> coords[cidx(i,0)] >> coords[cidx(i,1)] >> coords[cidx(i,2)];
    }
    input.close();
}

