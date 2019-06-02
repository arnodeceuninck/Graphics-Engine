//
// Created by arno on 26.02.19.
//

#ifndef ENGINE_FACE_H
#define ENGINE_FACE_H

#include <vector>

class Face {
public:
    // De indexen refereren naar
    // punten in de ‘points’ vector
    // van de Figure-klasse
    std::vector<int> point_indexes;

    Face(const std::vector<int> &point_indexes);
};


#endif //ENGINE_FACE_H
