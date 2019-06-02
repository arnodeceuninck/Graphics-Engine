//
// Created by arno on 3/4/19.
//

#ifndef ENGINE_MATRIX_H
#define ENGINE_MATRIX_H

#include "vector3d.h"
#include <cmath>

namespace TransMatrix{
    Matrix scaleFigure(const double &scale);

    Matrix rotateX(const double &angle);

    Matrix rotateY(const double &angle);

    Matrix rotateZ(const double &angle);

    Matrix translate(const Vector3D &vector);

    void toPolar(
            const Vector3D &point,
            double &theta,
            double &phi,
            double &r);

    Matrix invEyePointTrans(const Vector3D &eyepoint);
    Matrix eyePointTrans(const Vector3D &eyepoint);

}

#endif //ENGINE_MATRIX_H
