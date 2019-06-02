//
// Created by arno on 3/13/19.
//

#ifndef ENGINE_LSYSTEM3D_H
#define ENGINE_LSYSTEM3D_H

#include "Figure.h"
#include "l_parser.h"

/**
* \brief Draw a 3DL-System
*
* \param l_system 	a LParser::LSystem3D containing the L-system that's going to be drawn
* \param size   the amount of pixels on the largest side of the image
* \param color    A vector<double> containing three numbers in [0; 1], the RGB values of the lines
*
* \note an l_system contains the angel in degrees, all local variables here will use the angel in radians
*/
Figure draw3DLSystem(const LParser::LSystem3D &l_system, std::vector<double> color);

/**
 * \brief Generate an L-System full length string
 *
 * @param str The index string
 * @param l_system The l_system (required for the replacement rules)
 * @param size The nr of iterations for the replacement
 * @return
 */
std::string replace_string(const std::string& str, const LParser::LSystem3D& l_system, const int& size);

#endif //ENGINE_LSYSTEM3D_H
