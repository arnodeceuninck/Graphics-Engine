//
// Created by arno on 3/10/19.
//

#include "3DPlatonic.h"
#include "TransMatrix.h"
#include <cmath>
#include <vector>

Figure removeMiddle(Figure figure, const int &iterations);

Figure createCube() {

    Figure cube;

    cube.points.emplace_back(Vector3D::point(1, -1, -1));
    cube.points.emplace_back(Vector3D::point(-1, 1, -1));
    cube.points.emplace_back(Vector3D::point(1, 1, 1));
    cube.points.emplace_back(Vector3D::point(-1, -1, 1));
    cube.points.emplace_back(Vector3D::point(1, 1, -1));
    cube.points.emplace_back(Vector3D::point(-1, -1, -1));
    cube.points.emplace_back(Vector3D::point(1, -1, 1));
    cube.points.emplace_back(Vector3D::point(-1, 1, 1));

    cube.faces.emplace_back(Face({1, 5, 3, 7}));
    cube.faces.emplace_back(Face({5, 2, 8, 3}));
    cube.faces.emplace_back(Face({2, 6, 4, 8}));
    cube.faces.emplace_back(Face({6, 1, 7, 4}));
    cube.faces.emplace_back(Face({7, 3, 8, 4}));
    cube.faces.emplace_back(Face({1, 6, 2, 5}));

    cube.decreaseFaces();

    return cube;

}

Figure createTetrahedron() {

    Figure tetrahedron;

    tetrahedron.points.emplace_back(Vector3D::point(1, -1, -1));
    tetrahedron.points.emplace_back(Vector3D::point(-1, 1, -1));
    tetrahedron.points.emplace_back(Vector3D::point(1, 1, 1));
    tetrahedron.points.emplace_back(Vector3D::point(-1, -1, 1));

    tetrahedron.faces.emplace_back(Face({1, 2, 3}));
    tetrahedron.faces.emplace_back(Face({2, 4, 3}));
    tetrahedron.faces.emplace_back(Face({1, 4, 2}));
    tetrahedron.faces.emplace_back(Face({1, 3, 4}));

    tetrahedron.decreaseFaces();

    return tetrahedron;

}

Figure createOctahedron() {

    Figure octahedron;

    octahedron.points.emplace_back(Vector3D::point(1, 0, 0));
    octahedron.points.emplace_back(Vector3D::point(0, 1, 0));
    octahedron.points.emplace_back(Vector3D::point(-1, 0, 0));
    octahedron.points.emplace_back(Vector3D::point(0, -1, 0));
    octahedron.points.emplace_back(Vector3D::point(0, 0, -1));
    octahedron.points.emplace_back(Vector3D::point(0, 0, 1));

    octahedron.faces.emplace_back(Face({1, 2, 6}));
    octahedron.faces.emplace_back(Face({2, 3, 6}));
    octahedron.faces.emplace_back(Face({3, 4, 6}));
    octahedron.faces.emplace_back(Face({4, 1, 6}));
    octahedron.faces.emplace_back(Face({2, 1, 5}));
    octahedron.faces.emplace_back(Face({3, 2, 5}));
    octahedron.faces.emplace_back(Face({4, 3, 5}));
    octahedron.faces.emplace_back(Face({1, 4, 5}));

    octahedron.decreaseFaces();

    return octahedron;

}

Figure createIcosahedron() {

    Figure icosahedron;

    icosahedron.points.emplace_back(Vector3D::point(0, 0, sqrt(5) / 2));

    for (int i = 1; i < 6; ++i) {
        icosahedron.points.emplace_back(Vector3D::point(cos((i - 1) * 2 * M_PI / 5), sin((i - 1) * 2 * M_PI / 5), 0.5));
    }
    for (int i = 6; i < 11; ++i) {
        icosahedron.points.emplace_back(
                Vector3D::point(cos(M_PI / 5 + (i - 6) * 2 * M_PI / 5), sin(M_PI / 5 + (i - 6) * 2 * M_PI / 5), -0.5));
    }

    icosahedron.points.emplace_back(Vector3D::point(0, 0, -sqrt(5) / 2));

    icosahedron.faces.emplace_back(Face({1, 2, 3}));
    icosahedron.faces.emplace_back(Face({1, 3, 4}));
    icosahedron.faces.emplace_back(Face({1, 4, 5}));
    icosahedron.faces.emplace_back(Face({1, 5, 6}));
    icosahedron.faces.emplace_back(Face({1, 6, 2}));
    icosahedron.faces.emplace_back(Face({2, 7, 3}));
    icosahedron.faces.emplace_back(Face({3, 7, 8}));
    icosahedron.faces.emplace_back(Face({3, 8, 4}));
    icosahedron.faces.emplace_back(Face({4, 8, 9}));
    icosahedron.faces.emplace_back(Face({4, 9, 5}));
    icosahedron.faces.emplace_back(Face({5, 9, 10}));
    icosahedron.faces.emplace_back(Face({5, 10, 6}));
    icosahedron.faces.emplace_back(Face({6, 10, 11}));
    icosahedron.faces.emplace_back(Face({6, 11, 2}));
    icosahedron.faces.emplace_back(Face({2, 11, 7}));
    icosahedron.faces.emplace_back(Face({12, 8, 7}));
    icosahedron.faces.emplace_back(Face({12, 9, 8}));
    icosahedron.faces.emplace_back(Face({12, 10, 9}));
    icosahedron.faces.emplace_back(Face({12, 11, 10}));
    icosahedron.faces.emplace_back(Face({12, 7, 11}));

    icosahedron.decreaseFaces();

    return icosahedron;
}

Figure createDodecahedron() {

    Figure icosahedron = createIcosahedron();
    Figure dodecahedron;

    for (const Face &triangle: icosahedron.faces) {

        double x = (icosahedron.points[triangle.point_indexes[0]].x +
                    icosahedron.points[triangle.point_indexes[1]].x +
                    icosahedron.points[triangle.point_indexes[2]].x) / 3;

        double y = (icosahedron.points[triangle.point_indexes[0]].y +
                    icosahedron.points[triangle.point_indexes[1]].y +
                    icosahedron.points[triangle.point_indexes[2]].y) / 3;

        double z = (icosahedron.points[triangle.point_indexes[0]].z +
                    icosahedron.points[triangle.point_indexes[1]].z +
                    icosahedron.points[triangle.point_indexes[2]].z) / 3;

        dodecahedron.points.emplace_back(Vector3D::point(x, y, z));

    }

    dodecahedron.faces = {};

    dodecahedron.faces.emplace_back(Face({1, 2, 3, 4, 5}));
    dodecahedron.faces.emplace_back(Face({1, 6, 7, 8, 2}));
    dodecahedron.faces.emplace_back(Face({2, 8, 9, 10, 3}));
    dodecahedron.faces.emplace_back(Face({3, 10, 11, 12, 4}));
    dodecahedron.faces.emplace_back(Face({4, 12, 13, 14, 5}));
    dodecahedron.faces.emplace_back(Face({5, 14, 15, 6, 1}));
    dodecahedron.faces.emplace_back(Face({20, 19, 18, 17, 16}));
    dodecahedron.faces.emplace_back(Face({20, 15, 14, 13, 19}));
    dodecahedron.faces.emplace_back(Face({19, 13, 12, 11, 18}));
    dodecahedron.faces.emplace_back(Face({18, 11, 10, 9, 17}));
    dodecahedron.faces.emplace_back(Face({17, 9, 8, 7, 16}));
    dodecahedron.faces.emplace_back(Face({16, 7, 6, 15, 20}));

    dodecahedron.decreaseFaces();

    return dodecahedron;

}

Figure createCylinder(int n, const double h, bool openBase) {

    Figure cylinder;

    for (int i = 0; i < n; ++i) {
        cylinder.points.push_back(Vector3D::point(cos(2 * i * M_PI / n), sin(2 * i * M_PI / n), h));
        cylinder.points.push_back(Vector3D::point(cos(2 * i * M_PI / n), sin(2 * i * M_PI / n), 0));
    }

    for (int j = 0; j < 2 * n; j += 2) {
        cylinder.faces.push_back(Face({j % (2 * n), (j + 1) % (2 * n), (j + 3) % (2 * n), (j + 2) % (2 * n)}));
    }

    if (!openBase) {
        std::vector<int> bottom;
        std::vector<int> top;

        for (int k = 0; k < 2 * n; k += 2) {
            top.push_back(k);
            bottom.push_back(2 * n - (k + 1));
        }

        cylinder.faces.emplace_back(bottom);
        cylinder.faces.emplace_back(top);
    }
    return cylinder;
}

Figure createCone(const int n, const double h) {
    Figure cone;

    for (int i = 0; i < n; ++i) {
        cone.points.push_back(Vector3D::point(cos(2 * i * M_PI / n), sin(2 * i * M_PI / n), 0));
    }

    cone.points.push_back(Vector3D::point(0, 0, h));

    for (int j = 0; j < n; ++j) {
        cone.faces.push_back(Face({n, j % n, (j + 1) % n}));
    }

    std::vector<int> bottom;

    for (int k = 0; k < n; ++k) {
        bottom.push_back(n - k - 1);
    }

    cone.faces.push_back(bottom);

    return cone;
}

Figure createSphere(const int n) {
    Figure sphere = createIcosahedron();
    for (int i = 0; i < n; ++i) {
        sphere.splitTriangles();
    }
    sphere.rescalePoints();
    return sphere;
}

int findIndex(int i, int j, int m) {
    return i * m + j;
}

Figure createTorus(double r, double R, int n, int m) {

    Figure torus;

    for (int i = 0; i < n; ++i) {

        double u = 2 * i * M_PI / n;

        for (int j = 0; j < m; j++) {

            double v = 2 * j * M_PI / m;

            torus.points.emplace_back(Vector3D::point((R + r * cos(v)) * cos(u),
                                                      (R + r * cos(v)) * sin(u),
                                                      r * sin(v)));

            torus.faces.emplace_back(Face({findIndex(i, j, m),
                                           findIndex((i + 1) % n, j, m),
                                           findIndex((i + 1) % n, (j + 1) % m, m),
                                           findIndex(i, (j + 1) % m, m)}));
        }
    }


    return torus;
}

void generateFractal(Figure &figr, Figures3D &fractal, int nr_iterations, const double &scale) {

    fractal.emplace_back(figr);

    while (nr_iterations > 0) {
        Figures3D newFractals;
        for (Figure fig: fractal) {
            for (int point = 0; point < fig.points.size(); ++point) {

                Figure fig2 = fig;
                fig2.applyTransformation(TransMatrix::scaleFigure(1 / scale));

                Vector3D pi = fig.points[point];
                Vector3D pib = fig2.points[point];

                fig2.applyTransformation(TransMatrix::translate(pi - pib));

                newFractals.emplace_back(fig2);
            }
        }
        fractal = newFractals;
        nr_iterations--;
    }
}

//Figure createBuckeyBall() {
//
//    Figure buckeyBall;
//
//    buckeyBall.points.emplace_back(Vector3D::point(0, 0, sqrt(5)/2));
//
//    for (int i = 1; i < 6; ++i) {
//        buckeyBall.points.emplace_back(Vector3D::point(cos((i-1)*2*M_PI/5), sin((i-1)*2*M_PI/5), 0.5));
//    }
//    for (int i = 6; i < 11; ++i) {
//        buckeyBall.points.emplace_back(Vector3D::point(cos(M_PI/5 + (i-6)*2*M_PI/5), sin(M_PI/5 + (i-6)*2*M_PI/5), -0.5));
//    }
//
//    buckeyBall.points.emplace_back(Vector3D::point(0, 0, -sqrt(5)/2));
//
//    buckeyBall.faces.emplace_back(Face({1, 2, 3}));
//    buckeyBall.faces.emplace_back(Face({1, 3, 4}));
//    buckeyBall.faces.emplace_back(Face({1, 4, 5}));
//    buckeyBall.faces.emplace_back(Face({1, 5, 6}));
//    buckeyBall.faces.emplace_back(Face({1, 6, 2}));
//    buckeyBall.faces.emplace_back(Face({2, 7, 3}));
//    buckeyBall.faces.emplace_back(Face({3, 7, 8}));
//    buckeyBall.faces.emplace_back(Face({3, 8, 4}));
//    buckeyBall.faces.emplace_back(Face({4, 8, 9}));
//    buckeyBall.faces.emplace_back(Face({4, 9, 5}));
//    buckeyBall.faces.emplace_back(Face({5, 9, 10}));
//    buckeyBall.faces.emplace_back(Face({5, 10, 6}));
//    buckeyBall.faces.emplace_back(Face({6, 10, 11}));
//    buckeyBall.faces.emplace_back(Face({6, 11, 2}));
//    buckeyBall.faces.emplace_back(Face({2, 11, 7}));
//    buckeyBall.faces.emplace_back(Face({12, 8, 7}));
//    buckeyBall.faces.emplace_back(Face({12, 9, 8}));
//    buckeyBall.faces.emplace_back(Face({12, 10, 9}));
//    buckeyBall.faces.emplace_back(Face({12, 11, 10}));
//    buckeyBall.faces.emplace_back(Face({12, 7, 11}));
//
//    buckeyBall.decreaseFaces();
//
//    return buckeyBall;
//
//}

Figure createBuckeyBall() {
    Figure buckeyBall = createIcosahedron();
    // Iedere driehoek opsplitsen in een zeshoek en drie driehoeken

    // Temp
    std::vector<Vector3D> newPoints;
//    std::vector<Face> newFaces;

//    std::vector<Vector3D> finalFaces;

    std::vector<Face> newFaces;

    for (int face = 0; face < buckeyBall.faces.size(); ++face) {
        Face currentFace = buckeyBall.faces[face];

        for (int i = 0; i < 3; ++i) {
//            newPoints.push_back((0.0 / 3) * (buckeyBall.points[currentFace.point_indexes[(i+1)%3]] +
//                                             (buckeyBall.points[currentFace.point_indexes[i%3]])));

//            newPoints.push_back(buckeyBall.points[currentFace.point_indexes[i%3]]);
            newPoints.push_back(buckeyBall.points[currentFace.point_indexes[i % 3]] +
                                // TODO: Werkt (a+b)/3 ook? (gebasseerd op formule gemiddelde)
                                (1.0 / 3) * (buckeyBall.points[currentFace.point_indexes[(i + 1) % 3]] -
                                             (buckeyBall.points[currentFace.point_indexes[i % 3]])));
            newPoints.push_back(buckeyBall.points[currentFace.point_indexes[i % 3]] +
                                (2.0 / 3) * (buckeyBall.points[currentFace.point_indexes[(i + 1) % 3]] -
                                             (buckeyBall.points[currentFace.point_indexes[i % 3]])));
        }

//        // Drie driehoeken
//        newFaces.push_back(Face({9*face+0, 9*face+1, 9*face+8}));
//        newFaces.push_back(Face({9*face+2, 9*face+3, 9*face+4}));
//        newFaces.push_back(Face({9*face+5, 9*face+6, 9*face+7}));
//
//        // En een zeshoek
//        newFaces.push_back(Face({9*face+1, 9*face+2, 9*face+4, 9*face+5, 9*face+7, 9*face+8}));

        std::vector<int> hexagon;
        for (int j = 6; j > 0; --j) {
            hexagon.push_back(newPoints.size() - j + 1); // +1 omdat we decreaseFaces gebruiken
        }
        newFaces.emplace_back(hexagon);

    }

    // Alle driehoeken met een gemeenschappelijk punt gaan samennemen en enkel het grondvlak overhouden
    // Vind eerst alle toppen. een

    buckeyBall.points = newPoints;
    buckeyBall.faces = newFaces;
//    buckeyBall.faces = {}; // TODO: remove
    // De zeshoeken
//    buckeyBall.faces.emplace_back(Face({1, 2, 3, 4, 5, 6}));
//    buckeyBall.faces.emplace_back(Face({6, 5, 9, 10, 11, 12}));
//    buckeyBall.faces.emplace_back(Face({12, 11, 15, 16, 17, 18}));
//    buckeyBall.faces.emplace_back(Face({18, 17, 21, 22, 23, 24}));
//    buckeyBall.faces.emplace_back(Face({24, 23, 27, 28, 2, 1}));
//    buckeyBall.faces.emplace_back(Face({31, 32, 33, 34, 4, 3}));
//    buckeyBall.faces.emplace_back(Face({34, 33, 39, 40, 41, 42}));
//    buckeyBall.faces.emplace_back(Face({42, 41, 45, 46, 10, 9}));
//    buckeyBall.faces.emplace_back(Face({46, 45, 51, 52, 53, 54}));
//    buckeyBall.faces.emplace_back(Face({54, 53, 57, 58, 16, 15}));
//    buckeyBall.faces.emplace_back(Face({58, 57, 63, 64, 65, 66}));
//    buckeyBall.faces.emplace_back(Face({66, 65, 69, 10, 22, 21}));
//    buckeyBall.faces.emplace_back(Face({10, 69, 75, 76, 77, 78}));
//    buckeyBall.faces.emplace_back(Face({78, 77, 81, 82, 28, 27}));
//    buckeyBall.faces.emplace_back(Face({82, 81, 87, 88, 32, 31}));
//    buckeyBall.faces.emplace_back(Face({91, 92, 40, 39, 95, 96}));
//    buckeyBall.faces.emplace_back(Face({108, 107, 52, 51, 92, 91}));
//    buckeyBall.faces.emplace_back(Face({103, 104, 64, 63, 107, 108}));
//    buckeyBall.faces.emplace_back(Face({107, 108, 76, 75, 104, 103}));
//    buckeyBall.faces.emplace_back(Face({96, 95, 88, 87, 108, 107}));

    // De vijfhoeken
    // TODO: fix vijfhoeken
    buckeyBall.faces.emplace_back(Face({1, 6, 12, 18, 24})); // 1
    buckeyBall.faces.emplace_back(Face({28, 82, 31, 3, 2})); // 2
    buckeyBall.faces.emplace_back(Face({9, 5, 4, 34, 42})); // 3
    buckeyBall.faces.emplace_back(Face({10, 46, 54, 15, 11})); // 4
    buckeyBall.faces.emplace_back(Face({17, 16, 58, 66, 21})); // 5
    buckeyBall.faces.emplace_back(Face({78, 27, 23, 22, 70})); // 6
    buckeyBall.faces.emplace_back(Face({33, 32, 117, 95, 39})); // 7
    buckeyBall.faces.emplace_back(Face({41, 40, 92, 51, 45})); // 8
    buckeyBall.faces.emplace_back(Face({53, 52, 107, 63, 57})); // 9
    buckeyBall.faces.emplace_back(Face({69, 65, 64, 104, 75})); // 10
    buckeyBall.faces.emplace_back(Face({110, 87, 81, 77, 76})); // 11 // getal 108 fout
    buckeyBall.faces.emplace_back(Face({96, 109, 103, 108, 91})); // 12

    buckeyBall.decreaseFaces();

    return buckeyBall;
}

Figure createMengerSponge(const int &nrIterations) {
    if (nrIterations == 0)
        return createCube();

    std::vector<Vector3D> cube;

    cube.emplace_back(Vector3D::point(1, -1, -1)); // 1
    cube.emplace_back(Vector3D::point(-1, 1, -1)); // 2
    cube.emplace_back(Vector3D::point(1, 1, 1)); // 3
    cube.emplace_back(Vector3D::point(-1, -1, 1)); // 4
    cube.emplace_back(Vector3D::point(1, 1, -1)); // 5
    cube.emplace_back(Vector3D::point(-1, -1, -1)); // 6
    cube.emplace_back(Vector3D::point(1, -1, 1)); // 7
    cube.emplace_back(Vector3D::point(-1, 1, 1)); // 8
    cube.emplace_back(Vector3D::point(1, 0, -1)); // tussen 1 en 5
    cube.emplace_back(Vector3D::point(-1, 0, -1)); // tussen 2 en 6
    cube.emplace_back(Vector3D::point(1, 0, 1));
    cube.emplace_back(Vector3D::point(-1, 0, 1));
    cube.emplace_back(Vector3D::point(0, 1, -1)); // tussen 2 en 5
    cube.emplace_back(Vector3D::point(0, -1, -1));
    cube.emplace_back(Vector3D::point(0, 1, 1));
    cube.emplace_back(Vector3D::point(0, -1, 1));
    cube.emplace_back(Vector3D::point(-1, -1, 0)); // tussen 3 en 5
    cube.emplace_back(Vector3D::point(1, -1, 0));
    cube.emplace_back(Vector3D::point(1, 1, 0));
    cube.emplace_back(Vector3D::point(-1, 1, 0));

    std::vector<Figure> sponge;

    for (int i = 0; i < cube.size(); ++i) {
        Figure fig = createMengerSponge(nrIterations - 1);
        fig.applyTransformation(TransMatrix::scaleFigure(1.0 / 3));
        fig.applyTransformation(TransMatrix::translate(Vector3D::vector((2.0 / 3) * cube[i].x,
                                                                        (2.0 / 3) * cube[i].y,
                                                                        (2.0 / 3) * cube[i].z)));
        sponge.emplace_back(fig);
    }

    Figure mengerSponge;

    for (int j = 0; j < sponge.size(); ++j) {
        mengerSponge.mergeFigure(sponge[j]);
    }

    return mengerSponge;
}

Figure removeMiddle(Figure figure, const int &iterations) {

    Figure newFigure;

    std::vector<Vector3D> newPoints;
    std::vector<Face> newFaces;

    for (int face = 0; face < figure.faces.size(); ++face) {
        Face currentFace = figure.faces[face];
    }

    figure.points = newPoints;
    figure.faces = newFaces;

    return Figure();
}

