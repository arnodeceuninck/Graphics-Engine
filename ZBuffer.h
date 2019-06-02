//
// Created by arno on 3/19/19.
//

#ifndef ENGINE_ZBUFFER_H
#define ENGINE_ZBUFFER_H


#include <vector>

class ZBuffer {
public:
    //Constructor: maakt een Z-Buffer van de correcte
    //grootte aan en initialiseert alle velden op +inf
    ZBuffer(const int width, const int height);
    bool check(const int x, const int y, const double zvalue);
    double get(const int x, const int y) const;
    int getSize() const { return buffer.size(); };
    ZBuffer();
    void generateImage(std::string filename);
    double minValue;
    double maxValue;

private:
    std::vector<std::vector<double>> buffer;
};


#endif //ENGINE_ZBUFFER_H
