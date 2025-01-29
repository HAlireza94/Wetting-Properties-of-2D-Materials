#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <string>
#include <iomanip>
//#include <bits/stdc++.h>
#include <set>

using namespace std;

struct AtomPosition {
    double x, y;
    AtomPosition(double x, double y) : x(x), y(y) {}

    bool operator<(const AtomPosition& other) const {
        return x < other.x || (x == other.x && y < other.y);
    }
};

double pbcWrap(double coord, double L) {
    if (coord < 0) {
        coord += L * ceil(abs(coord) / L);
    } else if (coord >= L) {
        coord -= L * floor(coord / L);
    }
    return coord;
}



int main() {

    double Lx = 63.25252;  // Desired length in x-direction in Angstroms
    double Ly = 49.19;     // Desired length in y-direction in Angstroms
    double a = 4.05;      // Lattice constant in Angstroms for FCC

    double distancenN = a * sqrt(2) / 2;   // Nearest-neighbor distance
    double distancebRows = a * sqrt(3) / 2;  // Distance between rows

    // Calculate the number of atoms and rows to fit within Lx and Ly
    float numberatomsX = static_cast<int>(Lx / distancenN);
    float numberRows = static_cast<int>(Ly / distancebRows);


    set<AtomPosition> uniquePositions;

    for (int i = 0; i < numberRows; i++) {
        for (int j = 0; j < numberatomsX; j++) {
            float x = j * distancenN + (i % 2) * (distancenN / 2);
            float y = i * distancebRows;

            x = pbcWrap(x, Lx);
            y = pbcWrap(y, Ly);

            uniquePositions.insert(AtomPosition(x, y));
        }
    }


    ofstream outFile("FCC.gro");
    ofstream outXYZ("FCC.xyz");
    int totalnumberAtoms = numberatomsX * numberRows;

    outFile << "Repulsive sheet" << '\n';
    outFile << to_string(totalnumberAtoms) << '\n';

    outXYZ << to_string(totalnumberAtoms) << '\n';
    outXYZ << "Repulsive sheet" << '\n';

    char atomName = 'X';
    float z = 0;





    int ii = 1;
    for (const auto& pos : uniquePositions) {

        outXYZ <<atomName << "     " << fixed << setprecision(6)<< pos.x << "     " << fixed << setprecision(6) << pos.y <<  "     " << fixed << setprecision(6) << z << "\n";
        outFile <<"    "<<ii<<"UFF      " << atomName << "     " << fixed << setprecision(6)<< pos.x/10 << "     " << fixed << setprecision(6) << pos.y/10 <<  "     " << fixed << setprecision(6) << z << "\n";
    }

    outFile <<"   "<< Lx/10 <<"   "<< Ly/10 <<"   "<< 10 << endl;
    outFile.close();
    outXYZ.close();

    return 0;
}
