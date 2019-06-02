cmake_minimum_required(VERSION 3.13)
project(Engine)

set(CMAKE_CXX_STANDARD 11)

add_executable(Engine
        engine.cc
        easy_image.cc
        easy_image.h
        ini_configuration.cc
        ini_configuration.h
        l_parser.cc
        l_parser.h
        ColorD.cpp
        ColorD.h
        Point2D.cpp
        Point2D.h
        Line2D.cpp
        Line2D.h
        TransMatrix.h
        Figure.cpp
        Figure.h
        vector3d.cc
        vector3d.h
        Face.cpp
        Face.h 3DPlatonic.cpp 3DPlatonic.h LSystem3D.cpp LSystem3D.h ZBuffer.cpp ZBuffer.h TransMatrix.cpp Light.cpp Light.h)
