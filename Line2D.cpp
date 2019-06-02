//
// Created by arno on 2/24/19.
//

#include "Line2D.h"

Line2D::Line2D(const Point2D &p1, const Point2D &p2, const ColorD &colorD) : p1(p1), p2(p2), colorD(colorD) {}

void Line2D::setP1(const Point2D &p1) {
    Line2D::p1 = p1;
}

void Line2D::setP2(const Point2D &p2) {
    Line2D::p2 = p2;
}

void Line2D::scale(const double& factor) {
    p1.x *= factor;
    p1.y *= factor;
    p2.x *= factor;
    p2.y *= factor;
}

void Line2D::move(const double &dx, const double &dy) {
    p1.x += dx;
    p1.y += dy;
    p2.x += dx;
    p2.y += dy;
}

Line2D::Line2D(const Point2D &p1, const Point2D &p2) : p1(p1), p2(p2), colorD(ColorD(1, 1,1 )) {}
