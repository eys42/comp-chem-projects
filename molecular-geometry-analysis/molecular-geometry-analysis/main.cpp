//
//  main.cpp
//  molecular-geometry-analysis
//
//  Created by Ethan Song on 12/16/25.
//

#include "molecule.h"


int main(int argc, const char * argv[]) {
    Molecule acetaldehyde("/Users/ethan/Documents/ProgrammingProjects/molecular-geometry-analysis/molecular-geometry-analysis/input/acetaldehyde.dat", 0, 1);
    acetaldehyde.print_geom();
    acetaldehyde.compute_distances();
    acetaldehyde.print_distance_matrix();
    acetaldehyde.compute_angles();
    acetaldehyde.print_angle_matrix();
    return EXIT_SUCCESS;
}
