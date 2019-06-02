//
// Created by arno on 26.02.19.
//

#ifndef ENGINE_FIGURE_H
#define ENGINE_FIGURE_H

#include "Face.h"
#include "vector3d.h"
#include "ColorD.h"
#include "Line2D.h"
#include <list>

class Figure;
typedef std::list<Figure> Figures3D;

class Figure {
public:
    std::vector<Vector3D> points;
    std::vector<Face> faces;

//    ColorD color;
    ColorD ambientReflection;
    ColorD diffuseReflection;
    ColorD specularReflection;

    double  reflectionCoeff;

    Lines2D doProjection();
    void applyTransformation(const Matrix &transformation);
    void decreaseFaces();

    void splitTriangles();

    void rescalePoints();

    void trinagulate();

    std::vector<Face> triangulate_face(const Face& face);

    void mergeFigure(Figure fig);

    /**
     *
     * @param resultingFigures
     * @param r The radius of the spheres and cilinders
     * @param n The number of side faces of the cilinders
     * @param m The number of iterations to generate the spheres with
     */
    void generateThickFigure(Figures3D& resultingFigures, const double r, const int n, const int m);

    bool pointInFace(int i);
};


#endif //ENGINE_FIGURE_H
