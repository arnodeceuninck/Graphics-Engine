//
// Created by arno on 2/24/19.
//

#include "TransMatrix.h"

Matrix TransMatrix::scaleFigure(const double &scale) {
    Matrix scale_matrix;
    for (int i = 1; i < 4; ++i) {
        scale_matrix(i, i) = scale;
    }
    return scale_matrix;
}

Matrix TransMatrix::rotateX(const double &angle) {
    Matrix rotate_x_matrix;
    rotate_x_matrix(2, 2) = std::cos(angle);
    rotate_x_matrix(2, 3) = std::sin(angle);
    rotate_x_matrix(3, 3) = std::cos(angle);
    rotate_x_matrix(3, 2) = -std::sin(angle);
    return rotate_x_matrix;
}

Matrix TransMatrix::rotateY(const double &angle) {
    Matrix rotate_y_matrix;
    rotate_y_matrix(1, 1) = std::cos(angle);
    rotate_y_matrix(3, 1) = std::sin(angle);
    rotate_y_matrix(3, 3) = std::cos(angle);
    rotate_y_matrix(1, 3) = -std::sin(angle);
    return rotate_y_matrix;
}

Matrix TransMatrix::rotateZ(const double &angle) {
    Matrix rotate_z_matrix;
    rotate_z_matrix(1, 1) = std::cos(angle);
    rotate_z_matrix(2, 1) = -std::sin(angle);
    rotate_z_matrix(2, 2) = std::cos(angle);
    rotate_z_matrix(1, 2) = std::sin(angle);
    return rotate_z_matrix;
}

Matrix TransMatrix::translate(const Vector3D &vector) {
    Matrix translate_matrix;
    translate_matrix(4, 1) = vector.x;
    translate_matrix(4, 2) = vector.y;
    translate_matrix(4, 3) = vector.z;
    return translate_matrix;
}

void TransMatrix::toPolar(
        const Vector3D &point,
        double &theta,
        double &phi,
        double &r){
    r = sqrt(point.x*point.x + point.y*point.y + point.z*point.z);
    theta = std::atan2(point.y, point.x);
    phi = std::acos(point.z / r);
}

Matrix TransMatrix::eyePointTrans(const Vector3D &eyepoint){
    double theta;
    double phi;
    double r;
    toPolar(eyepoint, theta, phi, r);
    Matrix eyePointTransMatrix = rotateZ(-1*M_PI / 2 - theta) * rotateX(-phi) * translate(Vector3D::point(0, 0, -r));
    return eyePointTransMatrix;
}

Matrix TransMatrix::invEyePointTrans(const Vector3D &eyepoint) {
    return Matrix::inv(eyePointTrans(eyepoint));
}
