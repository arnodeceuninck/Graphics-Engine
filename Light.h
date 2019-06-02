//
// Created by arno on 30.04.19.
//

#ifndef ENGINE_LIGHT_H
#define ENGINE_LIGHT_H


#include "ColorD.h"
#include "vector3d.h"
#include "ZBuffer.h"
#include <list>

class Light
{
public:
    Light();
    Light(std::vector<double> ambient, std::vector<double> diffuse, std::vector<double> specular);
    //de ambiente licht component
    ColorD ambientLight;
    //de diffuse licht component
    ColorD diffuseLight;
    //de diffuse licht component
    ColorD specularLight;
    virtual std::string getType();
    virtual double alpha(Vector3D direction);
    virtual Vector3D getLd();
    bool hasReflection();
    ZBuffer shadowMask;
    Matrix eye;
    Matrix invEye;
    Vector3D eLightPos;
    double d, dx, dy;
};

class InfLight: public Light
{
public:
    //de richting waarin het
    //licht schijnt
    Vector3D ldVector;

    InfLight(const std::vector<double> &ambient, const std::vector<double> &diffuse,
             const std::vector<double> &specular, const Vector3D &ldVector);

    virtual std::string getType();
    virtual double alpha(Vector3D direction);
    virtual Vector3D getLd();
};

class PointLight: public Light
{
public:
    //de locatie van de puntbron
    Vector3D location;

    PointLight(const std::vector<double> &ambient, const std::vector<double> &diffuse,
               const std::vector<double> &specular, const Vector3D &location, const Vector3D& eyepoint);
    PointLight(const std::vector<double> &ambient, const std::vector<double> &diffuse,
               const std::vector<double> &specular, const Vector3D &location, int maskSize, const Vector3D& eyepoint);

    virtual std::string getType();
    virtual Vector3D getLd();

};

typedef std::list<Light*> Lights3D;

#endif //ENGINE_LIGHT_H
