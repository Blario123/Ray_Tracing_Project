//
// Created by blario123 on 22/04/2022.
//

#ifndef RAYTRACINGPROJECT_H
#define RAYTRACINGPROJECT_H

#include <cmath>
#include "Vec3.h"

static const double pi = M_PI;
static const double pi_reciprocal = M_1_PI;

static double damping_factor = 0.5;

static Vec3 rotate_x_axis(const Vec3 &position, const double &angle) {
    return {position.x, position.y * cos(angle) - position.z * sin(angle), position.y * sin(angle) + position.z * cos(angle)};
}

static Vec3 rotate_y_axis(const Vec3 &position, const double &angle) {
    return {position.x * cos(angle) + position.z * sin(angle), position.y, -position.x * sin(angle) + position.z * cos(angle)};
}

static Vec3 rotate_z_axis(const Vec3 &position, const double &angle) {
    return {position.x * cos(angle) - position.y * sin(angle), position.x * sin(angle) + position.y * cos(angle), position.z};
}

static Vec3 move_torus_1(const Vec3 &position) {
    return rotate_y_axis(rotate_z_axis(position - Vec3(2.0, 0.0, 0.55), M_PI / 4.0), M_PI / 2.0);
}

static Vec3 move_cube(const Vec3 &position) {
    return position - Vec3(2.0, 0.0, 0.0);
}

static Vec3 move_tetra(const Vec3 &position) {
    return rotate_z_axis(rotate_x_axis(position - Vec3(2.0, 0.0, 0.17), M_PI / 3.0), M_PI / 8.0);
}

#endif //RAYTRACINGPROJECT_H
