//
// Created by arno on 2/24/19.
//

#ifndef ENGINE_LINE2D_H
#define ENGINE_LINE2D_H


#include <list>
#include "Point2D.h"
#include "ColorD.h"

class Line2D;
typedef std::list<Line2D> Lines2D;

class Line2D {
public:

    void setP1(const Point2D &p1);

    void setP2(const Point2D &p2);

    void scale(const double& factor);

    void move(const double& dx, const double& dy);

    Line2D(const Point2D &p1, const Point2D &p2);

    Line2D(const Point2D &p1, const Point2D &p2, const ColorD &colorD);

    Point2D p1;
    Point2D p2;
    ColorD colorD;

};


#endif //ENGINE_LINE2D_H
