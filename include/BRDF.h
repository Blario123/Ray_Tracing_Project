#ifndef BRDF_H
#define BRDF_H

#include "RayTracingProject.h"
#include "Vec3.h"

class BRDF {
public:
    BRDF(float r_, float g_, float b_) {
        pR = r_;
        pG = g_;
        pB = b_;
    };
    Radiance radiance() {
        return damping_factor * pi_reciprocal * Vec3(pR, pG, pB);
    }
private:
    float pR, pG, pB;
};

#endif // BRDF_H
