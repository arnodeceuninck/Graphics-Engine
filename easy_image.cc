/*
 * easy_image.cc
 * Copyright (C) 2011  Daniel van den Akker
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "easy_image.h"
#include "TransMatrix.h"
//#include "Light.h"
#include <algorithm>
#include <assert.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include <fstream>

#define le32toh(x) (x)

inline int roundInt(double d) {
    return static_cast<int>(std::round(d));
}

namespace {
    //structs borrowed from wikipedia's article on the BMP file format
    struct bmpfile_magic {
        uint8_t magic[2];
    };

    struct bmpfile_header {
        uint32_t file_size;
        uint16_t reserved_1;
        uint16_t reserved_2;
        uint32_t bmp_offset;
    };
    struct bmp_header {
        uint32_t header_size;
        int32_t width;
        int32_t height;
        uint16_t nplanes;
        uint16_t bits_per_pixel;
        uint32_t compress_type;
        uint32_t pixel_size;
        int32_t hres;
        int32_t vres;
        uint32_t ncolors;
        uint32_t nimpcolors;
    };

    //copy-pasted from lparser.cc to allow these classes to be used independently from each other
    class enable_exceptions {
    private:
        std::ios &ios;
        std::ios::iostate state;
    public:
        enable_exceptions(std::ios &an_ios, std::ios::iostate exceptions) :
                ios(an_ios) {
            state = ios.exceptions();
            ios.exceptions(exceptions);
        }

        ~enable_exceptions() {
            ios.exceptions(state);
        }
    };

    //helper function to convert a number (char, int, ...) to little endian
    //regardless of the endiannes of the system
    //more efficient machine-dependent functions exist, but this one is more portable
    template<typename T>
    T to_little_endian(T value) {
        //yes, unions must be used with caution, but this is a case in which a union is needed
        union {
            T t;
            uint8_t bytes[sizeof(T)];
        } temp_storage;

        for (uint8_t i = 0; i < sizeof(T); i++) {
            temp_storage.bytes[i] = value & 0xFF;
            value >>= 8;
        }
        return temp_storage.t;
    }

    template<typename T>
    T from_little_endian(T value) {
        //yes, unions must be used with caution, but this is a case in which a union is needed
        union {
            T t;
            uint8_t bytes[sizeof(T)];
        } temp_storage;
        temp_storage.t = value;
        T retVal = 0;

        for (uint8_t i = 0; i < sizeof(T); i++) {
            retVal = (retVal << 8) | temp_storage.bytes[sizeof(T) - i - 1];
        }
        return retVal;
    }

}

img::Color::Color() :
        blue(0), green(0), red(0) {
}

img::Color::Color(uint8_t r, uint8_t g, uint8_t b) :
        blue(b), green(g), red(r) {
}

img::Color::~Color() {
}

img::UnsupportedFileTypeException::UnsupportedFileTypeException(std::string const &msg) :
        message(msg) {
}

img::UnsupportedFileTypeException::UnsupportedFileTypeException(const UnsupportedFileTypeException &original)
        : std::exception(original), message(original.message) {
}

img::UnsupportedFileTypeException::~UnsupportedFileTypeException() throw() {
}

img::UnsupportedFileTypeException &
img::UnsupportedFileTypeException::operator=(UnsupportedFileTypeException const &original) {
    this->message = original.message;
    return *this;
}

const char *img::UnsupportedFileTypeException::what() const throw() {
    return message.c_str();
}

img::EasyImage::EasyImage() :
        width(0), height(0), bitmap(), texelmap() {
}

img::EasyImage::EasyImage(unsigned int _width, unsigned int _height, Color color) :
        width(_width), height(_height), bitmap(width * height, color), texelmap(width * height, Texel(-1, -1)) {
}

img::EasyImage::EasyImage(EasyImage const &img) :
        width(img.width), height(img.height), bitmap(img.bitmap), texelmap(img.texelmap) {
}

img::EasyImage::~EasyImage() {
    bitmap.clear();
    texelmap.clear();
}

img::EasyImage &img::EasyImage::operator=(img::EasyImage const &img) {
    width = img.width;
    height = img.height;
    bitmap.assign(img.bitmap.begin(), img.bitmap.end());
    texelmap.assign(img.texelmap.begin(), img.texelmap.end());
    return (*this);
}

unsigned int img::EasyImage::get_width() const {
    return width;
}

unsigned int img::EasyImage::get_height() const {
    return height;
}

void img::EasyImage::clear(Color color) {
    for (std::vector<Color>::iterator i = bitmap.begin(); i != bitmap.end(); i++) {
        *i = color;
    }
}

img::Color &img::EasyImage::operator()(unsigned int x, unsigned int y) {
    assert(x < this->width);
    assert(y < this->height);
    return bitmap.at(x * height + y);
}

img::Color const &img::EasyImage::operator()(unsigned int x, unsigned int y) const {
    assert(x < this->width);
    assert(y < this->height);
    return bitmap.at(x * height + y);
}

void img::EasyImage::draw_line(unsigned int x0, unsigned int y0, unsigned int x1, unsigned int y1, Color color) {
    assert(x0 < this->width && y0 < this->height);
    assert(x1 < this->width && y1 < this->height);
    if (x0 == x1 and y0 == y1) {
        return;
    }
    if (x0 == x1) {
        //special case for x0 == x1
        for (unsigned int i = std::min(y0, y1); i <= std::max(y0, y1); i++) {
            (*this)(x0, i) = color;
        }
    } else if (y0 == y1) {
        //special case for y0 == y1
        for (unsigned int i = std::min(x0, x1); i <= std::max(x0, x1); i++) {
            (*this)(i, y0) = color;
        }
    } else {
        if (x0 > x1) {
            //flip points if x1>x0: we want x0 to have the lowest value
            std::swap(x0, x1);
            std::swap(y0, y1);
        }
        double m = ((double) y1 - (double) y0) / ((double) x1 - (double) x0);
        if (-1.0 <= m && m <= 1.0) {
            for (unsigned int i = 0; i <= (x1 - x0); i++) {
                (*this)(x0 + i, (unsigned int) round(y0 + m * i)) = color;
            }
        } else if (m > 1.0) {
            for (unsigned int i = 0; i <= (y1 - y0); i++) {
                (*this)((unsigned int) round(x0 + (i / m)), y0 + i) = color;
            }
        } else if (m < -1.0) {
            for (unsigned int i = 0; i <= (y0 - y1); i++) {
                (*this)((unsigned int) round(x0 - (i / m)), y0 - i) = color;
            }
        }
    }
}

std::ostream &img::operator<<(std::ostream &out, EasyImage const &image) {

    //temporaryily enable exceptions on output stream
    enable_exceptions(out, std::ios::badbit | std::ios::failbit);
    //declare some struct-vars we're going to need:
    bmpfile_magic magic;
    bmpfile_header file_header;
    bmp_header header;
    uint8_t padding[] =
            {0, 0, 0, 0};
    //calculate the total size of the pixel data
    unsigned int line_width = image.get_width() * 3; //3 bytes per pixel
    unsigned int line_padding = 0;
    if (line_width % 4 != 0) {
        line_padding = 4 - (line_width % 4);
    }
    //lines must be aligned to a multiple of 4 bytes
    line_width += line_padding;
    unsigned int pixel_size = image.get_height() * line_width;

    //start filling the headers
    magic.magic[0] = 'B';
    magic.magic[1] = 'M';

    file_header.file_size = to_little_endian(pixel_size + sizeof(file_header) + sizeof(header) + sizeof(magic));
    file_header.bmp_offset = to_little_endian(sizeof(file_header) + sizeof(header) + sizeof(magic));
    file_header.reserved_1 = 0;
    file_header.reserved_2 = 0;
    header.header_size = to_little_endian(sizeof(header));
    header.width = to_little_endian(image.get_width());
    header.height = to_little_endian(image.get_height());
    header.nplanes = to_little_endian(1);
    header.bits_per_pixel = to_little_endian(24);//3bytes or 24 bits per pixel
    header.compress_type = 0; //no compression
    header.pixel_size = pixel_size;
    header.hres = to_little_endian(11811); //11811 pixels/meter or 300dpi
    header.vres = to_little_endian(11811); //11811 pixels/meter or 300dpi
    header.ncolors = 0; //no color palette
    header.nimpcolors = 0;//no important colors

    //okay that should be all the header stuff: let's write it to the stream
    out.write((char *) &magic, sizeof(magic));
    out.write((char *) &file_header, sizeof(file_header));
    out.write((char *) &header, sizeof(header));

    //okay let's write the pixels themselves:
    //they are arranged left->right, bottom->top, b,g,r
    for (unsigned int i = 0; i < image.get_height(); i++) {
        //loop over all lines
        for (unsigned int j = 0; j < image.get_width(); j++) {
            //loop over all pixels in a line
            //we cast &color to char*. since the color fields are ordered blue,green,red they should be written automatically
            //in the right order
            out.write((char *) &image(j, i), 3 * sizeof(uint8_t));
        }
        if (line_padding > 0)
            out.write((char *) padding, line_padding);
    }
    //okay we should be done
    return out;
}

std::istream &img::operator>>(std::istream &in, EasyImage &image) {
    enable_exceptions(in, std::ios::badbit | std::ios::failbit);
    //declare some struct-vars we're going to need
    bmpfile_magic magic;
    bmpfile_header file_header;
    bmp_header header;
    //a temp buffer for reading the padding at the end of each line
    uint8_t padding[] =
            {0, 0, 0, 0};

    //read the headers && do some sanity checks
    in.read((char *) &magic, sizeof(magic));
    if (magic.magic[0] != 'B' || magic.magic[1] != 'M')
        throw UnsupportedFileTypeException("Could not parse BMP File: invalid magic header");
    in.read((char *) &file_header, sizeof(file_header));
    in.read((char *) &header, sizeof(header));
    if (le32toh(header.pixel_size) + le32toh(file_header.bmp_offset) != le32toh(file_header.file_size))
        throw UnsupportedFileTypeException("Could not parse BMP File: file size mismatch");
    if (le32toh(header.header_size) != sizeof(header))
        throw UnsupportedFileTypeException("Could not parse BMP File: Unsupported BITMAPV5HEADER size");
    if (le32toh(header.compress_type) != 0)
        throw UnsupportedFileTypeException("Could not parse BMP File: Only uncompressed BMP files can be parsed");
    if (le32toh(header.nplanes) != 1)
        throw UnsupportedFileTypeException("Could not parse BMP File: Only one plane should exist in the BMP file");
    if (le32toh(header.bits_per_pixel) != 24)
        throw UnsupportedFileTypeException("Could not parse BMP File: Only 24bit/pixel BMP's are supported");
    //if height<0 -> read top to bottom instead of bottom to top
    bool invertedLines = from_little_endian(header.height) < 0;
    image.height = std::abs(from_little_endian(header.height));
    image.width = std::abs(from_little_endian(header.width));
    unsigned int line_padding = from_little_endian(header.pixel_size) / image.height - (3 * image.width);
    //re-initialize the image bitmap
    image.bitmap.clear();
    image.bitmap.assign(image.height * image.width, Color());
    //okay let's read the pixels themselves:
    //they are arranged left->right., bottom->top if height>0, top->bottom if height<0, b,g,r
    for (unsigned int i = 0; i < image.get_height(); i++) {
        //loop over all lines
        for (unsigned int j = 0; j < image.get_width(); j++) {
            //loop over all pixels in a line
            //we cast &color to char*. since the color fields are ordered blue,green,red, the data read should be written in the right variables
            if (invertedLines) {
                //store top-to-bottom
                in.read((char *) &image(j, image.height - 1 - i), 3 * sizeof(uint8_t));
            } else {
                //store bottom-to-top
                in.read((char *) &image(j, i), 3 * sizeof(uint8_t));
            }
        }
        if (line_padding > 0) {
            in.read((char *) padding, line_padding);
        }
    }
    //okay we're done
    return in;
}

void img::EasyImage::draw2DLines(Lines2D lines2D, const int size, bool zbuffered) {

    // Calculate xmin, xmax, ymin, ymax
    double xmin = lines2D.begin()->p1.x;
    double xmax = lines2D.begin()->p1.x;
    double ymin = lines2D.begin()->p1.y;
    double ymax = lines2D.begin()->p1.y;

    for (const Line2D &line : lines2D) {

        double p1x = line.p1.x;
        double p1y = line.p1.y;

        xmin = p1x < xmin ? p1x : xmin;
        ymin = p1y < ymin ? p1y : ymin;
        xmax = p1x > xmax ? p1x : xmax;
        ymax = p1y > ymax ? p1y : ymax;

        double p2x = line.p2.x;
        double p2y = line.p2.y;

        xmin = p2x < xmin ? p2x : xmin;
        ymin = p2y < ymin ? p2y : ymin;
        xmax = p2x > xmax ? p2x : xmax;
        ymax = p2y > ymax ? p2y : ymax;

    }

    double xrange = xmax - xmin;
    double yrange = ymax - ymin;

    double maxrange = xrange > yrange ? xrange : yrange;

    // Calculate the size of the image
    double image_x = size * xrange / maxrange;
    double image_y = size * yrange / maxrange;

    width = image_x;
    height = image_y;

    // Calculate scaling factor d
    double d = 0.95 * image_x / xrange;

    // Multiply the coordinates of all points with d
    for (Line2D &line2D : lines2D) {
        line2D.scale(d);
    }

    // Move the line drawing
    double DCx = d * (xmin + xmax) / 2;
    double DCy = d * (ymin + ymax) / 2;
    double dx = image_x / 2 - DCx;
    double dy = image_y / 2 - DCy;

//    std::cout << "Calculated variables xrange yrange image_x image_y d " << xrange << " " << yrange << " " << image_x
//              << " " << image_y << " " << d << std::endl;

    // Add (dx, dy) to all line coordinates
    for (Line2D &line2D : lines2D) {
        line2D.move(dx, dy);
    }

    int ymini = std::numeric_limits<int>::max();
    int ymaxi = -std::numeric_limits<int>::max();
    if (zbuffered) {

        ZBuffer zBuffer(width, height);

        for (const Line2D &line2D : lines2D) {

            ymini = line2D.p1.x < ymini ? line2D.p1.x : ymini;
            ymini = line2D.p2.x < ymini ? line2D.p2.x : ymini;
            ymaxi = line2D.p1.x > ymaxi ? line2D.p1.x : ymaxi;
            ymaxi = line2D.p2.x > ymaxi ? line2D.p2.x : ymaxi;


            this->draw_zbuf_line(zBuffer, roundInt(line2D.p1.x), roundInt(line2D.p1.y), line2D.p1.z,
                                 roundInt(line2D.p2.x), roundInt(line2D.p2.y), line2D.p2.z,
                                 img::Color(roundInt(line2D.colorD.red * 255),
                                            roundInt(line2D.colorD.green * 255),
                                            roundInt(line2D.colorD.blue * 255)));

        }
//        std::cout << "Calculated ymini ymaxi " << ymini << " " << ymaxi << std::endl;
    } else {

//        ZBuffer zBuffer(width, height);

        for (const Line2D &line2D : lines2D) {
            this->draw_line(roundInt(line2D.p1.x), roundInt(line2D.p1.y),
                            roundInt(line2D.p2.x), roundInt(line2D.p2.y),
                            img::Color(roundInt(line2D.colorD.red * 255),
                                       roundInt(line2D.colorD.green * 255),
                                       roundInt(line2D.colorD.blue * 255)));

        }
    }
}

void img::EasyImage::drawLSystem(const LParser::LSystem2D &l_system, const int &size, std::vector<double> color) {

    std::string str = replace_string(l_system.get_initiator(), l_system, l_system.get_nr_iterations());

    Lines2D lines;

    std::vector<std::vector<double>> round_brackets_stack = {};

    double current_angle = l_system.get_starting_angle() * M_PI / 180;
    Point2D current_position(0, 0);

    for (char c: str) {

        if (c == '+') {
            current_angle += l_system.get_angle() * M_PI / 180;
            continue;
        } else if (c == '-') {
            current_angle -= l_system.get_angle() * M_PI / 180;
            continue;
        } else if (c == '(') {
            round_brackets_stack.push_back({current_position.x, current_position.y, current_angle});
            continue;
        } else if (c == ')') {
            std::vector<double> old_location = round_brackets_stack.back();
            current_position = Point2D(old_location[0], old_location[1]);
            current_angle = old_location[2];
            round_brackets_stack.pop_back();
            continue;
        } else {
            Point2D old_position = current_position;
            current_position.move(std::cos(current_angle), std::sin(current_angle));

            if (l_system.draw(c)) {
                lines.push_back(Line2D(old_position, current_position, ColorD(color)));
            }
            continue;
        }
    }

    draw2DLines(lines, size, false);
    return;

}

std::string
img::EasyImage::replace_string(const std::string &str, const LParser::LSystem2D &l_system, const int &size) {

    if (size == 0) { return str; }

    std::string new_string = "";

    for (char c: str) {
        if (c == '+' or c == '-' or c == '(' or c == ')' or c == '^' or c == '&' or c == '/' or c == '\\' or c == '|') {
            new_string += c;
            continue;
        }
        new_string += l_system.get_replacement(c);
    }

    return replace_string(new_string, l_system, size - 1);
}

void img::EasyImage::drawFigures3D(const Figures3D &figures3D, const int size, bool zbuffered, Vector3D eyepoint) {

    Lines2D lines2D;

    for (Figure figure: figures3D) {
        figure.applyTransformation(TransMatrix::eyePointTrans(eyepoint));
        Lines2D lines = figure.doProjection();
        for (const Line2D &line: lines) {
            lines2D.push_back(line);
        }
    }

    draw2DLines(lines2D, size, zbuffered);
}

double zscore(int i, int a, double za, double zb) {
    if (za == 0 or zb == 0) {
    }
    return (i * 1.0 / a) / za + (1 - i * 1.0 / a) / zb;
}

void img::EasyImage::draw_zbuf_line(ZBuffer &zBuffer, unsigned int x0, unsigned int y0,
                                    double z0, unsigned int x1, unsigned int y1, double z1,
                                    const img::Color &color) {

    // TODO: swap van z waarden grondig controleren

    assert(x0 < this->width && y0 < this->height);
    assert(x1 < this->width && y1 < this->height);
    if (x0 == x1 and y0 == y1) {
        return;
    }
    if (x0 == x1) {
        //special case for x0 == x1
        for (unsigned int i = std::min(y0, y1); i <= std::max(y0, y1); i++) {
            if (zBuffer.check(x0, i, zscore(i - std::min(y0, y1), std::max(y0, y1) - std::min(y0, y1),
                                            y0 < y1 ? z0 : z1, y0 < y1 ? z1 : z0))) {
                (*this)(x0, i) = color;
            }
        }
    } else if (y0 == y1) {
        //special case for y0 == y1
        for (unsigned int i = std::min(x0, x1); i <= std::max(x0, x1); i++) {
            if (zBuffer.check(i, y0, zscore(i - std::min(x0, x1), std::max(x0, x1) - std::min(x0, x1),
                                            x0 < x1 ? z0 : z1, x0 < x1 ? z1 : z0))) {
                (*this)(i, y0) = color;
            }
        }
    } else {
        if (x0 > x1) {
            //flip points if x1>x0: we want x0 to have the lowest value
            std::swap(x0, x1);
            std::swap(y0, y1);
            std::swap(z0, z1);
        }
        double m = ((double) y1 - (double) y0) / ((double) x1 - (double) x0);
        if (-1.0 <= m && m <= 1.0) {
            for (unsigned int i = 0; i <= (x1 - x0); i++) {
                if (zBuffer.check(x0 + i, (unsigned int) round(y0 + m * i),
                                  zscore(i, x1 - x0, z0, z1))) {
                    (*this)(x0 + i, (unsigned int) round(y0 + m * i)) = color;
                }
            }
        } else if (m > 1.0) {
            for (unsigned int i = 0; i <= (y1 - y0); i++) {
                if (zBuffer.check((unsigned int) round(x0 + (i / m)), y0 + i,
                                  zscore(i, y1 - y0, z0, z1))) {
                    (*this)((unsigned int) round(x0 + (i / m)), y0 + i) = color;
                }
            }
        } else if (m < -1.0) {
            for (unsigned int i = 0; i <= (y0 - y1); i++) {
                if (zBuffer.check((unsigned int) round(x0 - (i / m)), y0 - i,
                                  zscore(i, y0 - y1, z0, z1))) {
                    (*this)((unsigned int) round(x0 - (i / m)), y0 - i) = color;
                }
            }
        }
    }

}

double saturate(double d) {
    if (d > 255) {
        return 255;
    } else {
        return d;
    }
}

double zCoord(const ZBuffer &zBuffer, int x, int y) {
    return 1.0 / zBuffer.get(x, y);
}

Vector3D pixel_to_eye(int xpixel, int ypixel, double dx, double dy, double d, ZBuffer &zBuffer) {
    double xproj = xpixel; // - dx;
    double yproj = ypixel; // - dy;

    double zEye = zCoord(zBuffer, xpixel, ypixel);

    double xEye = -zEye * (xproj - dx) / d; // todo: dx, dy: required or not?
    double yEye = -zEye * (yproj - dy) / d;
//    double xEye = -zEye * (xproj) / d; // todo: dx, dy: required or not?
//    double yEye = -zEye * (yproj) / d;

//    std::cout << xEye << " " << yEye << " " << zEye << std::endl; // x y en z zijn juist

    Vector3D point = Vector3D::point(xEye, yEye, zEye);

    return point;
};

double
surroundingzValues(double xcoord, double ycoord, const double &zcoord, const double &dx, const double &dy, double d,
                   const ZBuffer &zBuffer,
                   double &xpix, double &ypix) {
//    std::cout << "d value of light being: " << d << std::endl;
    xcoord = (xcoord * d) / (-zcoord) + dx;
    ycoord = (ycoord * d) / (-zcoord) + dy;

    xpix = xcoord;
    ypix = ycoord;

    double alphax = xcoord - floor(xcoord); // A = (x, y), B = (x+1, y), C = (x, y+1), D = (x+1, y+1)
    double alphay = ycoord - floor(ycoord);

//    int x = round(xcoord - alphax);
//    int y = round(ycoord - alphay);

    double zA = zBuffer.get(floor(xcoord), ceil(ycoord));
    double zB = zBuffer.get(ceil(xcoord), ceil(ycoord));
    double zC = zBuffer.get(floor(xcoord), floor(ycoord));
    double zD = zBuffer.get(ceil(xcoord), floor(ycoord));

//    std::cout << "Interpolating values: " << zA << " " << zB << " " << zC << " " << zD << std::endl;
//    std::cout << "With z-coordinates: " << 1.0/ zA << " " << 1.0/zB << " " << 1.0/zC << " " << 1.0/zD << std::endl;
//    std::cout << "Offsets: " << alphax << " " << alphay << std::endl;
//    double zA = zCoord(zBuffer, x, y);
//    double zB = zCoord(zBuffer, x + 1, y);
//    double zC = zCoord(zBuffer, x, y+1);
//    double zD = zCoord(zBuffer, x+1, y + 1);


    double zE = (1 - alphax) * zA + alphax * zB;
    double zF = (1 - alphax) * zC + alphax * zD;


    double zL = alphay * zE + (1 - alphay) * zF;


//    double xEye = -zEye*xproj/d;
//    double yEye = -zEye*yproj/d;

//    std::cout << xEye << " " << yEye << " " << zEye << std::endl; // x y en z zijn juist


//    Vector3D point = Vector3D::point(xEye, yEye, zEye);

    return zL;
}

img::Texel generateTexel(Vector3D P);

void img::EasyImage::diffuseColor(img::Color &color, int xpixel, int ypixel, const double &d, const double &dx,
                                  const double &dy, ZBuffer &zBuffer, const Vector3D &eyepoint, Lights3D &lights3D,
                                  const Vector3D &normaalVector, const img::Color &diffuseReflection,
                                  const img::Color &specularReflection, const double &reflectionCoeff, bool reflection,
                                  bool shadow, Matrix eyePointM, Matrix invEyepointMatrix, bool texture) {

//    std::cout << xpixel << " " << ypixel << std::endl;

    Vector3D point = pixel_to_eye(xpixel, ypixel, dx, dy, d, zBuffer);
    point *= invEyepointMatrix;
//    std::cout << point.x << " " << point.y << " " << point.z << std::endl; // lijken mij ook nog juist

    if (texture) {
        texelmap.at(xpixel * height + ypixel) = generateTexel(point);
    }

    for (Light *light: lights3D) {
        if (light->getType() == "POINT") {

            // Berekenen zbuffer waarde
            if (shadow) {
                Vector3D pointL = point * light->eye;

//                std::cout << light->dx << std::endl;
//                std::cout << "Check whether saved zL value is equal to calculated" << pointL.z << std::endl;
                double xcoord, ycoord;
                double zL = surroundingzValues(pointL.x, pointL.y, pointL.z, light->dx, light->dy, light->d,
                                               light->shadowMask, xcoord, ycoord);
//                std::cout << "Calculated:\t" << zL << std::endl;
//                std::cout << "Expected:\t" << 1.0/pointL.z << std::endl;
//                std::cout << "##############################" << std::endl;

//                int xpix = round(xcoord);
//                int ypix = round(ycoord);
//                double expectedZl = light->shadowMask.get(xpix, ypix);
                if (std::abs(zL - 1.0 / pointL.z) > std::pow(10, -4)) {
                    continue;
                }

//                std::cout << "Match found: " << zL << " " << 1.0/pointL.z << std::endl;
            }


//            std::cout << "YEEEETTT" << std::endl;

//            Vector3D lightPos = light->getLd();
//            Vector3D eLightPos = lightPos * TransMatrix::eyePointTrans(eyepoint);
            Vector3D ePoint = point * eyePointM; // TODO: Weird, want da zouden al eyepoint coords moeten zijn bij point

            Vector3D l = Vector3D::normalise(
                    light->eLightPos - ePoint); // de light->getLd() en point zijn allebei de orginele posities
            double alpha = normaalVector.x * l.x + normaalVector.y * l.y + normaalVector.z * l.z;
            if (alpha > 0) {
                color.red = saturate(color.red + diffuseReflection.red * light->diffuseLight.red * alpha);
                color.green = saturate(color.green + diffuseReflection.green * light->diffuseLight.green * alpha);
                color.blue = saturate(color.blue + diffuseReflection.blue * light->diffuseLight.blue * alpha);
            }

            if (reflection) {

                Vector3D r = Vector3D::normalise(2 * alpha * normaalVector - l);
//                Vector3D r = Vector3D::normalise(2 * alpha * normaalVector - eLightPos);
                Vector3D pointToEye = Vector3D::normalise(eyepoint * eyePointM - ePoint);

                double cosBeta = r.x * pointToEye.x + r.y * pointToEye.y + r.z * pointToEye.z;

//                std::cout << cosBeta << std::endl;

                if (cosBeta > 0) {
                    color.red = saturate(color.red + specularReflection.red * light->specularLight.red *
                                                     pow(cosBeta, reflectionCoeff));
                    color.green = saturate(color.green + specularReflection.green * light->specularLight.green *
                                                         pow(cosBeta, reflectionCoeff));
                    color.blue = saturate(color.blue + specularReflection.blue * light->specularLight.blue *
                                                       pow(cosBeta, reflectionCoeff));
                }
            }
        }
        if (light->getType() == "INFLIGHT" and reflection) {
//
//            Vector3D lightPos = light->getLd();
//            Vector3D eLightPos = lightPos; // * TransMatrix::eyePointTrans(eyepoint);
//            Vector3D ePoint = point * TransMatrix::eyePointTrans(eyepoint);
//
//            Vector3D l = Vector3D::normalise(eLightPos-ePoint); // de light->getLd() en point zijn allebei de orginele posities
//            double alpha = normaalVector.x * l.x + normaalVector.y * l.y + normaalVector.z * l.z;
//
//
////            Vector3D r = Vector3D::normalise(2 * alpha * normaalVector - l);
//            Vector3D r = Vector3D::normalise(2 * alpha * normaalVector - eLightPos);
//            Vector3D pointToEye = Vector3D::normalise(eyepoint * TransMatrix::eyePointTrans(eyepoint) - ePoint);
//
//            double cosBeta = r.x * pointToEye.x + r.y * pointToEye.y + r.z * pointToEye.z;
//
////                std::cout << cosBeta << std::endl;
//
//            if (cosBeta > 0) {
//                color.red = saturate(color.red + specularReflection.red * light->specularLight.red *
//                                                 pow(cosBeta, reflectionCoeff));
//                color.green = saturate(color.green + specularReflection.green * light->specularLight.green *
//                                                     pow(cosBeta, reflectionCoeff));
//                color.blue = saturate(color.blue + specularReflection.blue * light->specularLight.blue *
//                                                   pow(cosBeta, reflectionCoeff));
//            }


            // Alles in oorspronkelijke coordinaten
            Vector3D lightPos = light->getLd();
            Vector3D eLightPos = Vector3D::normalise(lightPos * TransMatrix::invEyePointTrans(eyepoint));
            Vector3D ePoint = point; // * TransMatrix::eyePointTrans(eyepoint);

//            Vector3D l = Vector3D::normalise(eLightPos-ePoint); // de light->getLd() en point zijn allebei de orginele posities
            Vector3D normaalVectorB = Vector3D::normalise(normaalVector * TransMatrix::invEyePointTrans(eyepoint));
            double alpha =
                    normaalVectorB.x * eLightPos.x + normaalVectorB.y * eLightPos.y + normaalVectorB.z * eLightPos.z;


//            Vector3D r = Vector3D::normalise(2 * alpha * normaalVector - l);
            Vector3D r = Vector3D::normalise(2 * alpha * normaalVector - eLightPos);
            Vector3D pointToEye = Vector3D::normalise(eyepoint - ePoint);

            double cosBeta = r.x * pointToEye.x + r.y * pointToEye.y + r.z * pointToEye.z;

//                std::cout << cosBeta << std::endl;

            if (cosBeta > 0) {
                color.red = saturate(color.red + specularReflection.red * light->specularLight.red *
                                                 pow(cosBeta, reflectionCoeff));
                color.green = saturate(color.green + specularReflection.green * light->specularLight.green *
                                                     pow(cosBeta, reflectionCoeff));
                color.blue = saturate(color.blue + specularReflection.blue * light->specularLight.blue *
                                                   pow(cosBeta, reflectionCoeff));
            }
        }
    }

    // Nu zou ik dus de oorspronkelijke x, y, z waarde moeten hebben
}


void
img::EasyImage::draw_zbuf_triag(ZBuffer &zBuffer, Vector3D const &A, Vector3D const &B, Vector3D const &C, double d,
                                double dx, double dy, img::Color ambientReflection, img::Color diffuseReflection,
                                img::Color specularReflection, double reflectionCoeff, Lights3D &lights,
                                const Vector3D &eyepoint, bool shadow, Matrix eyePointM, Matrix invEyepointMatrix,
                                bool texture) {

    // Berekenen van de normaalvector (mbv de eyepoint getransformeerde coordinaten A, B, C)
    Vector3D AB = Vector3D::normalise(B - A);
    Vector3D AC = Vector3D::normalise(C - A);
    Vector3D normaalVector = Vector3D::cross(AB, AC);
//    std::cout << normaalVector << std::endl;
    normaalVector = Vector3D::normalise(normaalVector);

    bool reflection = false;
    bool pointLight = false;

    // Calculate the color
    img::Color color(0, 0, 0);
    for (Light *light: lights) {
        color.red = saturate(color.red + ambientReflection.red * light->ambientLight.red);
        color.green = saturate(color.green + ambientReflection.green * light->ambientLight.green);
        color.blue = saturate(color.blue + ambientReflection.blue * light->ambientLight.blue);
        if (light->getType() == "INFLIGHT") {
            double alpha = light->alpha(normaalVector);
//            std::cout << alpha << std::endl;
            if (alpha > 0) {
                color.red = saturate(color.red + diffuseReflection.red * light->diffuseLight.red * alpha);
                color.green = saturate(color.green + diffuseReflection.green * light->diffuseLight.green * alpha);
                color.blue = saturate(color.blue + diffuseReflection.blue * light->diffuseLight.blue * alpha);
            }
        } else if (light->getType() == "POINT") {
            pointLight = true;
        }
        if (light->hasReflection()) {
            reflection = true;
        }
    }

    // Projecteren van de driehoek
    Point2D A2(d * A.x / (-A.z) + dx, d * A.y / (-A.z) + dy);
    Point2D B2(d * B.x / (-B.z) + dx, d * B.y / (-B.z) + dy);
    Point2D C2(d * C.x / (-C.z) + dx, d * C.y / (-C.z) + dy);


    // Bepalen welke pixels tot de driehoek behoren
    int ymax = roundInt(std::max(A2.y, std::max(B2.y, C2.y)) - 0.5); // aka yP id cursus
    int ymin = roundInt(std::min(A2.y, std::min(B2.y, C2.y)) + 0.5); // aka yQ in de cursus

    double xg = (A2.x + B2.x + C2.x) / 3;
    double yg = (A2.y + B2.y + C2.y) / 3;

    double zgr = 1 / (3 * A.z) + 1 / (3 * B.z) + 1 / (3 * C.z); // zgr = 1/zg


    // We gaan voor alle y1 waardes in Y
    for (int y = ymin; y <= ymax; ++y) {

        // Bepaal xL en xR
        // De xL en xR-waarde bepalen, eat aangeeft dat de pixels (xL, yI) tot (xR, yI) tot de driehoek behoren.

        // Om xL en xR te bepalen voor een enkele y-waarde yI gebruiken we 6 variabelen die we als volgt initialiseren
        double xLAB = std::numeric_limits<double>::infinity();
        double xLAC = std::numeric_limits<double>::infinity();
        double xLBC = std::numeric_limits<double>::infinity();

        double xRAB = -std::numeric_limits<double>::infinity();
        double xRAC = -std::numeric_limits<double>::infinity();
        double xRBC = -std::numeric_limits<double>::infinity();

        // We testen voor elk van de drie lijnstukken PQ gelijk aan AB, AC en BC of
        // TODO: En stel xL en xR hieraan gelijk: allebei dezelfde waarde
        // AB
        if ((y - A2.y) * (y - B2.y) <= 0 and A2.y != B2.y) {
            xLAB = B2.x + (A2.x - B2.x) * (y - B2.y) / (A2.y - B2.y);
            xRAB = xLAB;
        }
        // AC
        if ((y - A2.y) * (y - C2.y) <= 0 and A2.y != C2.y) {
            xLAC = C2.x + (A2.x - C2.x) * (y - C2.y) / (A2.y - C2.y);
            xRAC = xLAC;
        }
        // BC
        if ((y - B2.y) * (y - C2.y) <= 0 and C2.y != B2.y) {
            xLBC = C2.x + (B2.x - C2.x) * (y - C2.y) / (B2.y - C2.y);
            xRBC = xLBC;
        }

        if ((xLAB == std::numeric_limits<double>::infinity() and xLAC == xLAB and xLBC == xLAB) or
            (xRAB == -std::numeric_limits<double>::infinity() and xRAC == xRAB and xRBC == xRAB)) {
            return;
        }

        int xL = round(std::min(xLAB, std::min(xLAC, xLBC)) + 0.5);
        int xR = round(std::max(xRAB, std::max(xRAC, xRBC)) - 0.5);

        for (int x = xL; x <= xR; ++x) {

            Vector3D u = B - A;
            Vector3D v = C - A;
            Vector3D w = Vector3D::point(u.y * v.z - u.z * v.y,
                                         u.z * v.x - u.x * v.z,
                                         u.x * v.y - u.y * v.x);

            double k = w.x * A.x + w.y * A.y + w.z * A.z;

            double dzdx = w.x / (-d * k);
            double dzdy = w.y / (-d * k);

            double zr;
            if (shadow) {
                zr = zgr + (x - xg) * dzdx + (y - yg) * dzdy;
            } else {
                zr = 1.0001 * zgr + (x - xg) * dzdx + (y - yg) * dzdy; // zr = 1/z
            }

            if (zBuffer.check(x, y, zr)) {
                img::Color pixelColor = color;
                if (pointLight or reflection or texture)
                    diffuseColor(pixelColor, x, y, d, dx, dy, zBuffer, eyepoint, lights, normaalVector,
                                 diffuseReflection, specularReflection, reflectionCoeff, reflection, shadow, eyePointM,
                                 invEyepointMatrix, false);
//                if(texture){
//                    texelmap.at(x * height + y) = generateTexel();
//                }
                (*this)(x, y) = pixelColor;
            }
        }
    }
    // De 1/z waarde voor elke pixel berekenen
}

void calculateDDxDy(double &d, double &dx, double &dy, int &width, int &height, Figures3D &figures3D, const int size,
                    double &xmin, double &xmax, double &ymin, double &ymax) {
    Lines2D lines2D;

//    std::cout << "Calculating values for size " << size << std::endl;
    for (Figure &figure: figures3D) {
        Lines2D lines = figure.doProjection();
        for (const Line2D &line: lines) {
            lines2D.push_back(line);
        }
    }

//    img::EasyImage debugLines2D(size, size, img::Color(255, 255, 255));
//    debugLines2D.draw2DLines(lines2D, size, true);
//    std::string filename = "debug2DLines.bmp";
//    std::ofstream f_out(filename.c_str(),std::ios::trunc | std::ios::out | std::ios::binary);
//    f_out << debugLines2D;

    // Calculate xmin, xmax, ymin, ymax
    xmin = lines2D.begin()->p1.x;
    xmax = lines2D.begin()->p1.x;
    ymin = lines2D.begin()->p1.y;
    ymax = lines2D.begin()->p1.y;

    for (const Line2D &line : lines2D) {

        double p1x = line.p1.x;
        double p1y = line.p1.y;

        xmin = p1x < xmin ? p1x : xmin;
        ymin = p1y < ymin ? p1y : ymin;
        xmax = p1x > xmax ? p1x : xmax;
        ymax = p1y > ymax ? p1y : ymax;

        double p2x = line.p2.x;
        double p2y = line.p2.y;

        xmin = p2x < xmin ? p2x : xmin;
        ymin = p2y < ymin ? p2y : ymin;
        xmax = p2x > xmax ? p2x : xmax;
        ymax = p2y > ymax ? p2y : ymax;

    }

    double xrange = xmax - xmin;
    double yrange = ymax - ymin;

    double maxrange = xrange > yrange ? xrange : yrange;

    // Calculate the size of the image
    double image_x = size * xrange / maxrange;
    double image_y = size * yrange / maxrange;

    width = image_x;
    height = image_y;

//    std::cout << "The figure will be generated on an image with the following size: " << image_x << " " << image_y << std::endl;

    // Calculate scaling factor d
    d = 0.95 * image_x / xrange;

    //TODO: Niet doen bij shadowing
    // Multiply the coordinates of all points with d
    for (Line2D &line2D : lines2D) {
        line2D.scale(d);
    }

    // Move the line drawing
    double DCx = d * (xmin + xmax) / 2;
    double DCy = d * (ymin + ymax) / 2;
    dx = image_x / 2 - DCx;
    dy = image_y / 2 - DCy;

//    std::cout << "Calculated d value: " << d << std::endl;

//    std::cout << "Calculated variables xrange yrange image_x image_y d " << xrange << " " << yrange << " " << image_x
//              << " " << image_y << " " << d << std::endl;

}

void
makeZBuffer(ZBuffer &zBuffer, const double d, const double dx, const double dy, const Vector3D &A, const Vector3D &B,
            const Vector3D &C,
            bool shadow) {

//    std::cout << "Filling triangle in ZBuffer for A, B, C" << std::endl;
//    std::cout << "A = " << A << std::endl;
//    std::cout << "B = " << B << std::endl;
//    std::cout << "C = " << C << std::endl;

//    std::cout << "ZBuffer aan het berekenen met d waarde " << d << std::endl;
    // Projecteren van de driehoek
    Point2D A2(d * A.x / (-A.z) + dx, d * A.y / (-A.z) + dy);
    Point2D B2(d * B.x / (-B.z) + dx, d * B.y / (-B.z) + dy);
    Point2D C2(d * C.x / (-C.z) + dx, d * C.y / (-C.z) + dy);


    // Bepalen welke pixels tot de driehoek behoren
    int ymax = roundInt(std::max(A2.y, std::max(B2.y, C2.y)) - 0.5); // aka yP id cursus
    int ymin = roundInt(std::min(A2.y, std::min(B2.y, C2.y)) + 0.5); // aka yQ in de cursus

    double xg = (A2.x + B2.x + C2.x) / 3;
    double yg = (A2.y + B2.y + C2.y) / 3;

    double zgr = 1 / (3 * A.z) + 1 / (3 * B.z) + 1 / (3 * C.z); // zgr = 1/zg

//    std::cout << "Gravity point = " << xg << " " << yg << " " << 1.0/zgr << std::endl;
    // We gaan voor alle y1 waardes in Y
//    std::cout << "ymin ymax " << std::to_string(ymin) << " " << std::to_string(ymax) << std::endl;
    for (int y = ymin; y <= ymax; ++y) {

        // Bepaal xL en xR
        // De xL en xR-waarde bepalen, eat aangeeft dat de pixels (xL, yI) tot (xR, yI) tot de driehoek behoren.

        // Om xL en xR te bepalen voor een enkele y-waarde yI gebruiken we 6 variabelen die we als volgt initialiseren
        double xLAB = std::numeric_limits<double>::infinity();
        double xLAC = std::numeric_limits<double>::infinity();
        double xLBC = std::numeric_limits<double>::infinity();

        double xRAB = -std::numeric_limits<double>::infinity();
        double xRAC = -std::numeric_limits<double>::infinity();
        double xRBC = -std::numeric_limits<double>::infinity();

        // We testen voor elk van de drie lijnstukken PQ gelijk aan AB, AC en BC of
        // TODO: En stel xL en xR hieraan gelijk: allebei dezelfde waarde
        // AB
        if ((y - A2.y) * (y - B2.y) <= 0 and A2.y != B2.y) {
            xLAB = B2.x + (A2.x - B2.x) * (y - B2.y) / (A2.y - B2.y);
            xRAB = xLAB;
        }
        // AC
        if ((y - A2.y) * (y - C2.y) <= 0 and A2.y != C2.y) {
            xLAC = C2.x + (A2.x - C2.x) * (y - C2.y) / (A2.y - C2.y);
            xRAC = xLAC;
        }
        // BC
        if ((y - B2.y) * (y - C2.y) <= 0 and B2.y != C2.y) {
            xLBC = C2.x + (B2.x - C2.x) * (y - C2.y) / (B2.y - C2.y);
            xRBC = xLBC;
        }

        if ((xLAB == std::numeric_limits<double>::infinity() and xLAC == xLAB and xLBC == xLAB) or
            (xRAB == -std::numeric_limits<double>::infinity() and xRAC == xRAB and xRBC == xRAB)) {
//            std::cerr << "RIP" << y << std::endl;
//            return;
            continue;
        }

        int xL = round(std::min(xLAB, std::min(xLAC, xLBC)) + 0.5);
        int xR = round(std::max(xRAB, std::max(xRAC, xRBC)) - 0.5);

        for (int x = xL; x <= xR; ++x) {

            Vector3D u = B - A;
            Vector3D v = C - A;
            Vector3D w = Vector3D::point(u.y * v.z - u.z * v.y,
                                         u.z * v.x - u.x * v.z,
                                         u.x * v.y - u.y * v.x);

            double k = w.x * A.x + w.y * A.y + w.z * A.z;

            double dzdx = w.x / (-d * k);
            double dzdy = w.y / (-d * k);

            double zr;
            if (shadow) {
                zr = zgr + (x - xg) * dzdx + (y - yg) * dzdy;
            } else {
                zr = 1.0001 * zgr + (x - xg) * dzdx + (y - yg) * dzdy; // zr = 1/z
            }
//            std::cout << "zBuffer[" << x << "][" << y  << "] => " << zr << std::endl;
            zBuffer.check(x, y, zr);

        }
    }
    // De 1/z waarde voor elke pixel berekenen
}


void img::EasyImage::drawZbufTriangFig(Figures3D &figures3D, const int size, Lights3D lights, const Vector3D &eyepoint,
                                       bool shadows, bool background, std::string backgroundLocation, bool texture,
                                       std::string texturefile, bool ZbufferOutput) {

    for (Figure &figure: figures3D) {
        figure.trinagulate();
    }


    Matrix eyePointM = TransMatrix::eyePointTrans(eyepoint);
    Matrix invEyepointMatrix = Matrix::inv(eyePointM);

//    Matrix eyepointM = TransMatrix::eyePointTrans(eyepoint);

    if (shadows) {
//        std::cout << "Calculating shadowmasks" << std::endl;
        int i = 0;
        for (Light *light: lights) {
            Figures3D figuresCopy = figures3D;
            if (light->getType() == "POINT") {
                for (Figure &figure: figuresCopy) {
                    figure.applyTransformation(light->eye);
                }
//                light->eye = TransMatrix::eyePointTrans(light->getLd());
                int w, h;
                double xmin, xmax, ymin, ymax;
                calculateDDxDy(light->d, light->dx, light->dy, w, h, figuresCopy, light->shadowMask.getSize(),
                               xmin, xmax, ymin, ymax);
//                std::cout << "Calculated d value for light: " << light->d << std::endl;
                int j = 0;
                for (Figure &figure: figuresCopy) {
//                    figure.applyTransformation(light->eye);
//                    figure.applyTransformation(TransMatrix::scaleFigure(light->d));
                    for (Face &face: figure.faces) {
                        makeZBuffer(light->shadowMask, light->d, light->dx, light->dy,
                                    figure.points[face.point_indexes[0]],
                                    figure.points[face.point_indexes[1]],
                                    figure.points[face.point_indexes[2]], shadows);
//                        light->shadowMask.generateImage("LightZbufferPartial" + std::to_string(j) + ".bmp");
//                        std::cout << j << std::endl;
                        j++;
                    }
                }
//                for (Figure &figure: figures3D) {
//                    figure.applyTransformation(light->invEye);
//                }
            }
            if (ZbufferOutput) {
                std::string filename = "LightZbuffer" + std::to_string(i) + ".bmp";
                light->shadowMask.generateImage(filename);
            }
            i++;
        }
    }


//    std::cout << "Projecting figure" << std::endl;
    double d, dx, dy;
    int w, h;

    for (Figure &figure: figures3D) {
        figure.applyTransformation(eyePointM);
    }

    double xmin, xmax, ymin, ymax;
    calculateDDxDy(d, dx, dy, w, h, figures3D, size, xmin, xmax, ymin, ymax);
    width = w;
    height = h;

    if (texture) {
        setupTexture(figures3D, eyepoint, w, h, d);
    }

    if (background) {
        applyBackground(backgroundLocation);
    }

    ZBuffer zBuffer(width, height);

    for (Figure &figure : figures3D) {
        for (Face &face: figure.faces) {
            draw_zbuf_triag(zBuffer,
                            figure.points[face.point_indexes[0]],
                            figure.points[face.point_indexes[1]],
                            figure.points[face.point_indexes[2]],
                            d, dx, dy,
                            img::Color(figure.ambientReflection.red * 255, figure.ambientReflection.green * 255,
                                       figure.ambientReflection.blue * 255),
                            img::Color(figure.diffuseReflection.red * 255, figure.diffuseReflection.green * 255,
                                       figure.diffuseReflection.blue * 255),
                            img::Color(figure.specularReflection.red * 255, figure.specularReflection.green * 255,
                                       figure.specularReflection.blue * 255), figure.reflectionCoeff, lights, eyepoint,
                            shadows, eyePointM, invEyepointMatrix, false);
        }
    }

    if (texture) {
        img::EasyImage textureImage;
        std::ifstream fin(backgroundLocation);
        fin >> textureImage;
        fin.close();
        applyTexture(textureImage, 1);
    }

    if (ZbufferOutput) {
        zBuffer.generateImage("ZBuffer.bmp");
    }

}

void img::EasyImage::applyTexture(img::EasyImage texture, double reflection) {
    for (int l = 0; l < texelmap.size(); ++l) {
        if (texelmap[l].y != -1 and texelmap[l].x != -1) {
            img::Color texelColor = texture(texelmap[l].x, texelmap[l].y);
            img::Color figureColor = bitmap[l];
            img::Color newColor = img::Color(texelColor.red * reflection + figureColor.red * (1 - reflection),
                                             texelColor.green * reflection + figureColor.green * (1 - reflection),
                                             texelColor.blue * reflection + figureColor.blue * (1 - reflection));
            bitmap[l] = newColor;
        }
    }
}

img::Texel img::EasyImage::generateTexel(Vector3D P) {

//    Vector3D P; // het te projecteren punt
//    img::EasyImage texture;

    // Specifieren rechthoek die het oppervlak omvat
    // TODO: we willen alle geprojecteerde punten zien te omvatten, maar hoe is dit mogelijk? \
    //       * Bekijk de grootte van de texture en
    //         projecteer de figuur vanuit het eyepoint op een oppervlak van deze grootte \
    //         het vlak zal dus evenwijdig staan aan het vlak waarop we de afbeelding gaan projecteren \
    //         enige te bepalen factor is afstand tot het vlak. \
    //         zoek het meest linkse, rechtse, boven, onder geprojecteerde punt.\
    //         plaats de texture zodanig dat linkse, rechtse, boven, onder coordinaat hiermee overeen komen
    //       * Begin eerst met gewoon gelijk te stellen aan het vlak van waaruit we kijken. Daarna kunnen er nog dingen
    //         aangepast worden.

    Vector3D uvw = Vector3D::point(P.x - p.x, P.y - p.y, P.z - p.z) * Matrix::inv(abcMatrix);

    // Te gebruiken texel heeft het nummer (width * u, height * v)
//    img::Color texel;
    return img::Texel(textureWidth * uvw.x, textureHeight * uvw.y);
}

void img::EasyImage::setupTexture(Figures3D figures, Vector3D eyepoint, int pixelwidth, int pixelheight, double d) {

    // TODO: zoek meest linkse punt \
    //       zoek meest rechts gelegen punt \
    //       zoek meest boven punt \
    //       zoek meest onder punt \
    //       zoek meest vooraan gelegen punt \
    //       zoek meest vanacter gelegen punt \
    //       bereken het verschil in diepte \
    //       dit geeft een getal dat we uiteindelijk gaan moeten aftrekken van de berekende punten
    //       zoek de hoek linksonder
//    Vector3D p; // Een vector vanuit de oorsprong naar de linkerbenedenhoek van de rechthoek
    Vector3D a; // Start in p, eindigt in de rechterbenedenhoek
    Vector3D b; // Start in p, eindigt in linkerbovenhoek
    // TODO: vector p hier ook initialiseren
    p = Vector3D::vector(0, 0, -d);
    Vector3D upperCorner = Vector3D::vector(0, pixelheight, -d);
    Vector3D lowercorner = Vector3D::vector(pixelwidth, 0, -d);
    b = upperCorner - p;
    a = lowercorner - p;
    // Alle punten gelegen op dit vlak kunnen beschreven worden alsvolgt:
    // (x, y, z) = p + u*a + v*b ; Hierbij liggen u en v tussen 0 en 1
    // Texelnummer wordt bepaald door (width*u, height*v)

    // Her vlak ab wordt V genoemd, P' is de projectie van punt P op dit vlak
    // V (en dus ab) moet zodanig gekozen zijn dat deze alle geprojecteerde punten P' bevat

    Vector3D c = Vector3D::cross(a, b); // TODO: check of dit het juiste vectoriele product is

    // Bepaal coordinaten van P in assenstelsel met p als oorsprong en vectoren a, b, c als assenstelsen
    abcMatrix = Matrix();
    abcMatrix(1, 1) = a.x;
    abcMatrix(2, 1) = b.x;
    abcMatrix(3, 1) = c.x;
    abcMatrix(1, 1) = a.y;
    abcMatrix(2, 1) = b.y;
    abcMatrix(3, 1) = c.y;
    abcMatrix(1, 1) = a.z;
    abcMatrix(2, 1) = b.z;
    abcMatrix(3, 1) = c.z;

}

void img::EasyImage::applyBackground(std::string backgroundLocation) {

    img::EasyImage background;
    std::ifstream fin(backgroundLocation);
    fin >> background;
    fin.close();

    int backgroundWidth = background.get_width();
    int backgroundHeight = background.get_height();

    double scaleWidth = backgroundWidth * 1.0 / width;
    double scaleHeight = backgroundHeight * 1.0 / height;

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            (*this)(x, y) = background(std::floor(x * scaleWidth), std::floor(y * scaleHeight));
        }
    }

}

img::Texel::Texel(int x, int y) : x(x), y(y) {}

img::Texel::Texel() {}
