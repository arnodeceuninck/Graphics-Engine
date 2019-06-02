//
// Created by arno on 2/24/19.
//

#ifndef ENGINE_COLORD_H
#define ENGINE_COLORD_H

#include <vector>

class ColorD {
public:
    double red;
    double green;
    double blue;

    ColorD();

    ColorD(double red, double green, double blue);

    ColorD(std::vector<double> colors);
};


#endif //ENGINE_COLORD_H
