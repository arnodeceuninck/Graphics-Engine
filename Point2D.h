//
// Created by arno on 2/24/19.
//

#ifndef ENGINE_POINT2D_H
#define ENGINE_POINT2D_H


class Point2D {
public:
    Point2D(double x, double y);

    double x;
    double y;

    Point2D(double x, double y, double z);
    double  z;

    Point2D move(const double& dx, const double& dy);

    Point2D(const Point2D& point2D);
};


#endif //ENGINE_POINT2D_H
