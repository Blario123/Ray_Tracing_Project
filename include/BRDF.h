#ifndef BRDF_H
#define BRDF_H

#include "RayTracingProject.h"
#include "Vec3.h"
#include <QString>

class BRDF {
public:
    BRDF() {
        pR = 0.0;
        pG = 0.0;
        pB = 0.0;
        pName = "";
    };

    BRDF(float r_, float g_, float b_, const QString name = nullptr) {
        pR = r_;
        pG = g_;
        pB = b_;
        pName = name;
    };

    BRDF &operator=(const BRDF &o) {
        pR = o.pR;
        pG = o.pG;
        pB = o.pB;
        return *this;
    };

    float r() {
        return pR;
    };

    float g() {
        return pG;
    };

    float b() {
        return pB;
    };

    Radiance radiance() {
        return damping_factor * pi_reciprocal * Vec3(pR, pG, pB);
    };

    const QString name() {
        return pName;
    };
private:
    float pR, pG, pB;
    QString pName;
};

inline BRDF red = {1.0f, 0.2f, 0.2f, "Red"};
inline BRDF blue = {0.2f, 0.2f, 1.0f, "Blue"};
inline BRDF white = {1.0f, 1.0f, 1.0f, "White"};
inline BRDF pink = {1.0f, 0.2f, 1.0f, "Pink"};
inline BRDF purple = {0.75f, 0.25f, 0.75f, "Purple"};
inline BRDF yellow = {1.0f, 1.0f, 0.0f, "Yellow"};
inline BRDF brown = {0.625f, 0.32f, 0.18f, "Brown"};
inline BRDF cyan = {0.0f, 1.0f, 0.76f, "Cyan"};

#endif // BRDF_H
