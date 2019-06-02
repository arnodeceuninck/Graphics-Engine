//
// Created by arno on 3/13/19.
//

#include "LSystem3D.h"
#include <cmath>

Figure draw3DLSystem(const LParser::LSystem3D &l_system, std::vector<double> color) {

    std::string str = replace_string(l_system.get_initiator(), l_system, l_system.get_nr_iterations());

    Figure lines;
    lines.ambientReflection = ColorD(color);

    Vector3D H = Vector3D::vector(1, 0, 0);
    Vector3D L = Vector3D::vector(0, 1, 0);
    Vector3D U = Vector3D::vector(0, 0, 1);

    std::vector<std::vector<Vector3D>> round_brackets_stack = {};

    Vector3D current_position = Vector3D::point(0, 0, 0);

    double angle = l_system.get_angle() * M_PI / 180;

    for(char c: str){

        if(c == '+'){
            Vector3D newH = H*cos(angle) + L * sin(angle);
            Vector3D newL = -H * sin(angle) + L*cos(angle);
            H = newH;
            L = newL;
            continue;
        } else if(c == '-'){
            Vector3D newH = H*cos(-angle) + L * sin(-angle);
            Vector3D newL = -H * sin(-angle) + L*cos(-angle);
            H = newH;
            L = newL;
            continue;
        } else if (c == '^') {
            Vector3D newH = H*cos(angle) + U*sin(angle);
            Vector3D newU = -H*sin(angle) + U*cos(angle);
            H = newH;
            U = newU;
        } else if (c == '&') {
            Vector3D newH = H*cos(-angle) + U*sin(-angle);
            Vector3D newU = -H*sin(-angle) + U*cos(-angle);
            H = newH;
            U = newU;
        } else if (c == '\\') {
            Vector3D newL = L*cos(angle) - U*sin(angle);
            Vector3D newU = L*sin(angle) + U*cos(angle);
            L = newL;
            U = newU;
        } else if (c == '/') {
            Vector3D newL = L*cos(-angle) - U*sin(-angle);
            Vector3D newU = L*sin(-angle) + U*cos(-angle);
            L = newL;
            U = newU;
        } else if (c == '|') {
            H = -H;
            L = -L;
        } else if(c=='(') {
            round_brackets_stack.push_back({current_position, H, L, U});
            continue;
        } else if(c == ')') {
            std::vector<Vector3D> old_location = round_brackets_stack.back();
            current_position = old_location[0];
            H = old_location[1];
            L = old_location[2];
            U = old_location[3];
            round_brackets_stack.pop_back();
            continue;
        } else {
            Vector3D old_position = current_position;
            current_position += H;
            if (l_system.draw(c)) {
                lines.points.push_back(old_position);
                lines.points.push_back(current_position);
                lines.faces.emplace_back(Face({static_cast<int>(lines.points.size() - 1),
                                               static_cast<int>(lines.points.size() - 2)}));
            }
            continue;
        }
    }

    return lines;
}

std::string replace_string(const std::string &str, const LParser::LSystem3D &l_system, const int &size) {

    if(size == 0){ return str; }

    std::string new_string = "";

    for (char c: str){
        if (c == '+' or c == '-' or c == '(' or c ==')' or c == '^' or c == '&' or c == '/' or c == '\\' or c == '|'){
            new_string += c;
            continue;
        }
        new_string += l_system.get_replacement(c);
    }

    return replace_string(new_string, l_system, size-1);
}