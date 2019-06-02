//
// Created by arno on 2/24/19.
//

#include "Point2D.h"

Point2D::Point2D(double x, double y) : x(x), y(y), z(0) {}

Point2D Point2D::move(const double &dx, const double &dy) {
    x += dx;
    y += dy;
    return *this;
}

Point2D::Point2D(const Point2D &point2D) : x(point2D.x), y(point2D.y), z(point2D.z){}

Point2D::Point2D(double x, double y, double z) : x(x), y(y), z(z){

}

