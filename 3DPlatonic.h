//
// Created by arno on 3/10/19.
//

#ifndef ENGINE_3DPLATONIC_H
#define ENGINE_3DPLATONIC_H

#include "Figure.h"

Figure createCube();
Figure createTetrahedron();
Figure createOctahedron();
Figure createIcosahedron();
Figure createDodecahedron();
Figure createCylinder(int n, const double h, bool openBase);
Figure createCone(int n, double h);

Figure createSphere(int n);
Figure createTorus(double r, double R, int n, int m);

void generateFractal(Figure& fig, Figures3D& fractal,
                     int nr_iterations, const double& scale);

Figure createBuckeyBall();

Figure createMengerSponge(const int& nrIterations);

#endif //ENGINE_3DPLATONIC_H
