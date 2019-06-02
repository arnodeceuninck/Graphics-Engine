//
// Created by arno on 26.02.19.
//

#include "Figure.h"
#include "Point2D.h"
#include "3DPlatonic.h"
#include "TransMatrix.h"
#include <cmath>

Point2D do2DProjection(const Vector3D &point, const double d) {
    double x = d * point.x / -point.z;
    double y = d * point.y / -point.z;
    return Point2D(x, y, point.z);
}

void Figure::applyTransformation(const Matrix &transformation) {
    for(Vector3D& point: points){
        point *= transformation;
    }
}

Lines2D Figure::doProjection() {
    std::vector<Point2D> points2D;

    for(Vector3D vector3D: points){
        points2D.emplace_back(do2DProjection(vector3D, 1));
    }

    Lines2D lines2D;
    for(Face face: faces){
        for (int i = 0; i < face.point_indexes.size(); ++i) {
            lines2D.emplace_back(Line2D(points2D[face.point_indexes[i%face.point_indexes.size()]],
                    points2D[face.point_indexes[(i+1)%face.point_indexes.size()]], ambientReflection));
        }
    }

    return lines2D;
}

void Figure::decreaseFaces() {
    for(Face &face: faces){
        for(int &point: face.point_indexes){
            point -= 1;
        }
    }
}

void Figure::splitTriangles() {

    std::vector<Face> new_faces = {};

    int i = 0;

    for(Face face: faces){

        i++;

        const Vector3D& A = points[face.point_indexes[0]];
        const Vector3D& B = points[face.point_indexes[1]];
        const Vector3D& C = points[face.point_indexes[2]];

        Vector3D D = Vector3D::point((A.x + B.x)/2,
                                     (A.y + B.y)/2,
                                     (A.z + B.z)/2);
        Vector3D E = Vector3D::point((A.x + C.x)/2,
                                     (A.y + C.y)/2,
                                     (A.z + C.z)/2);
        Vector3D F = Vector3D::point((B.x + C.x)/2,
                                     (B.y + C.y)/2,
                                     (B.z + C.z)/2);

        points.push_back(D);
        points.push_back(E);
        points.push_back(F);

        int Ai = face.point_indexes[0];
        int Bi = face.point_indexes[1];
        int Ci = face.point_indexes[2];
        int Di = points.size()-3;
        int Ei = points.size()-2;
        int Fi = points.size()-1;

        new_faces.push_back(Face({Ai, Di, Ei}));
        new_faces.push_back(Face({Bi, Fi, Di}));
        new_faces.push_back(Face({Ci, Ei, Fi}));
        new_faces.push_back(Face({Di, Fi, Ei}));

    }

    faces = new_faces;
}

void Figure::rescalePoints() {
    for(Vector3D &point: points){
        point.normalise();
    }
}

std::vector<Face> Figure::triangulate_face(const Face &face) {
    int n = face.point_indexes.size();
    if(n<=3){
        return {face};
    }
    std::vector<Face> faces;
    for (int j = 1; j <= n-2 ; ++j) {
//        Fi = [P0, Pi, Pi+1]
        Face Fi = Face({face.point_indexes[0], face.point_indexes[j], face.point_indexes[j+1]});
        faces.emplace_back(Fi);
    }
    return faces;
}

void Figure::trinagulate() {
    std::vector<Face> newFaces;
    for(const Face& face: faces){
        std::vector<Face> triangulated = triangulate_face(face);
        for (int i = 0; i < triangulated.size(); ++i) {
            newFaces.emplace_back(triangulated[i]);
        }
    }
    faces = newFaces;
}

void Figure::mergeFigure(Figure fig) {

    int pointSize = points.size();

    for(Vector3D point: fig.points){
        points.emplace_back(point);
    }

    for(Face face: fig.faces){
        for(int &point: face.point_indexes){
            point += pointSize;
        }

        faces.emplace_back(face);
    }

}

void Figure::generateThickFigure(Figures3D &resultingFigures, const double r, const int n, const int m) {
    resultingFigures = {};

    // Generating all spheres
    int i = -1;
    for(Vector3D& point: points){
        i++;
//        if(!pointInFace(i)){
//            continue;
//        } else if(i > 0 and i < points.size()-1){
//            Vector3D p1p2 = points[i] - points[i-1];
//            Vector3D p2p3 = points[i+1] - points[i];
////            std::cout << "p1p2: " << p1p2 << std::endl;
////            std::cout << "p2p3; " << p2p3 << std::endl << std::endl;
//            if(p1p2.x == p2p3.x and p1p2.y == p1p2.y and p1p2.z == p1p2.z){
//                std::cout << "At least someonoe" << std::endl;
//                continue;
//            }
////            std::cout << "Compared neq: " << p1p2 << " " << p2p3 << std::endl;
//        }
        Figure sphere = createSphere(m);
        sphere.applyTransformation(TransMatrix::scaleFigure(r));
        sphere.applyTransformation(TransMatrix::translate(point));
        resultingFigures.push_back(sphere);
    }

    // Generating all cylinders
    for(Face& face: faces){
        for (int l = 0; l < face.point_indexes.size(); ++l) {

           Vector3D p1 = points[face.point_indexes[l]];
           Vector3D p2 = points[face.point_indexes[(l+1)%face.point_indexes.size()]];

           Vector3D p1p2 = p2 - p1;
           double h = p1p2.length() / r;
           Figure cylinder = createCylinder(n, h, false);

           cylinder.applyTransformation(TransMatrix::scaleFigure(r));

           Vector3D pointP1P2 = Vector3D::point(0, 0, 0) + p1p2;
           double theta, phi, rp;
           TransMatrix::toPolar(pointP1P2, theta, phi, rp);

//         cylinder.applyTransformation(TransMatrix::rotateZ(phi)); // TODO: check Z axis
//          cylinder.applyTransformation(TransMatrix::rotateX(theta)); // TODO: check X axis
           cylinder.applyTransformation(TransMatrix::rotateY(phi)); // TODO: check Z axis
           cylinder.applyTransformation(TransMatrix::rotateZ(theta)); // TODO: check X axis
           cylinder.applyTransformation(TransMatrix::translate(p1));

           resultingFigures.push_back(cylinder);
        }
    }
//    Figure figure = *this;
//    resultingFigures.push_back(figure); // TODO: remove this
}

bool Figure::pointInFace(int i) {
    for(Face& face: faces){
        for (int& j: face.point_indexes){
            if(i == j){
                return true;
            }
        }
    }
    return false;
}
