//
// Created by arno on 30.04.19.
//

#include "Light.h"
#include "TransMatrix.h"

Light::Light() {
    ambientLight = ColorD(1, 1, 1);
    diffuseLight = ColorD(0, 0, 0);
    specularLight = ColorD(0, 0, 0);
}

Light::Light(std::vector<double> ambient, std::vector<double> diffuse, std::vector<double> specular) {
    ambientLight = ColorD(ambient[0], ambient[1], ambient[2]);
    diffuseLight = ColorD(diffuse[0], diffuse[1], diffuse[2]);
    specularLight = ColorD(specular[0], specular[1], specular[2]);
}

std::string Light::getType() {
    return "NONE";
}

double Light::alpha(Vector3D direction) {
    return 0;
}

Vector3D Light::getLd() {
    return Vector3D();
}

bool Light::hasReflection() {
    return !(specularLight.red == 0 and specularLight.blue == 0 and specularLight.green == 0);
}

double InfLight::alpha(Vector3D direction) {
    return ldVector.x * direction.x + ldVector.y * direction.y + ldVector.z * direction.z;
}

std::string InfLight::getType() {
    return "INFLIGHT";
}

InfLight::InfLight(const std::vector<double> &ambient, const std::vector<double> &diffuse,
                   const std::vector<double> &specular, const Vector3D &ldVector) : Light(ambient, diffuse, specular),
                                                                                    ldVector(ldVector) {}

Vector3D InfLight::getLd() {
    return ldVector;
}

std::string PointLight::getType() {
    return "POINT";
}

PointLight::PointLight(const std::vector<double> &ambient, const std::vector<double> &diffuse,
                       const std::vector<double> &specular, const Vector3D &location, const Vector3D& eyepoint) : Light(ambient, diffuse,
                                                                                              specular),
                                                                                        location(location) {
    shadowMask = ZBuffer();
    eye = TransMatrix::eyePointTrans(location);
    invEye = Matrix::inv(eye);
    eLightPos = location * TransMatrix::eyePointTrans(eyepoint);
}

Vector3D PointLight::getLd() {
    return location;
}

PointLight::PointLight(const std::vector<double> &ambient, const std::vector<double> &diffuse,
                       const std::vector<double> &specular, const Vector3D &location, int maskSize, const Vector3D& eyepoint) : Light(ambient, diffuse,
                                                                                                            specular),
                                                                                                      location(location) {
    shadowMask = ZBuffer(maskSize, maskSize);
    eye = TransMatrix::eyePointTrans(location);
    invEye = Matrix::inv(eye);
    eLightPos = location * TransMatrix::eyePointTrans(eyepoint);
}
