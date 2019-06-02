//
// Created by arno on 2/24/19.
//

#include "ColorD.h"

ColorD::ColorD(double red, double green, double blue) : red(red), green(green), blue(blue) {}

ColorD::ColorD() : red(0), green(0), blue(0) {}

ColorD::ColorD(std::vector<double> colors) : red(colors[0]), green(colors[1]), blue(colors[2]){}
