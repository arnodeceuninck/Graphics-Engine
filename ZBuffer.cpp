//
// Created by arno on 3/19/19.
//

#include <limits>
#include <iostream>
#include <fstream>
#include "ZBuffer.h"
#include "easy_image.h"

ZBuffer::ZBuffer(const int width, const int height) {

    for (int j = 0; j < width; ++j) {

        std::vector<double> row;
        for (int i = 0; i < height; ++i) {
            row.emplace_back(std::numeric_limits<double>::infinity());
        }

        buffer.emplace_back(row);
        minValue = std::numeric_limits<double >::infinity();
        maxValue = -std::numeric_limits<double >::infinity();
    }
}

bool ZBuffer::check(const int x, const int y, const double zvalue) {
    if (x >= getSize()  or y >= buffer[0].size() or x < 0 or y < 0){
//        std::cerr << "Warning: skipped request for ZBuffer " << x << " " << y << std::endl;
        return false;
    }

    // Check whether it's closer
    double a = buffer[x][y];
//    std::cout << x << " " << y << ' ' << zvalue << std::endl;
    if(zvalue < buffer[x][y]){
        buffer[x][y] = zvalue;
        minValue = minValue < zvalue ? minValue : zvalue;
        maxValue = maxValue > zvalue ? maxValue : zvalue;
        return true;
    }
    return false;
}

double ZBuffer::get(const int x, const int y) const {
    if (x >= getSize()  or y >= buffer[0].size() or x < 0 or y < 0){
//        std::cerr << "Warning: skipped request for ZBuffer " << x << " " << y << std::endl;
        return std::numeric_limits<double >::infinity();
    }
    return buffer[x][y];
}

ZBuffer::ZBuffer() {

}

void ZBuffer::generateImage(std::string filename) {
    int width = getSize();
    int height = buffer[0].size();
    // Scale all depths to become a number between 0.2 and 0.6;
    double range = maxValue - minValue;

    img::EasyImage image(width, height);
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            double colorFactor = ((buffer[x][y] - minValue)/range)*0.4 + 0.2 ;
//            if(colorFactor < 9999 and colorFactor > -9999) {
//                std::cout << colorFactor << std::endl;
//            }
            image(x, y) = img::Color(255*colorFactor,
                                     255*colorFactor,
                                     255*colorFactor);
        }
    }
 std::ofstream f_out(filename.c_str(),std::ios::trunc | std::ios::out | std::ios::binary);
                                        f_out << image;
}
