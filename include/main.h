#ifndef MAIN_H
#define MAIN_H

#include "Ray.h"
#include "Image.h"
#include "Vec3.h"
#include "Observer.h"
#include "PhysicalObject.h"
#include "RayTracingProject.h"
#include "SceneRender.h"
#include "SDF.h"
#include "GUI.h"

Radiance red_BRDF(const Vec3 &position, const Vec3 &incident_light_vector, const Vec3 &outgoing_light_vector) {
	return damping_factor * pi_reciprocal * Radiance(1.0, 0.2, 0.2);
}

Radiance blue_BRDF(const Vec3 &position, const Vec3 &incident_light_vector, const Vec3 &outgoing_light_vector) {
	return damping_factor * pi_reciprocal * Radiance(0.2, 0.2, 1.0);
}

Radiance white_BRDF(const Vec3 &position, const Vec3 &incident_light_vector, const Vec3 &outgoing_light_vector) {
	return damping_factor * pi_reciprocal * Radiance(1.0, 1.0, 1.0);
}

Radiance pink_BRDF(const Vec3 &position, const Vec3 &incident_light_vector, const Vec3 &outgoing_light_vector) {
	return damping_factor * pi_reciprocal * Radiance(1.0, 0.2, 1.0);
}

Radiance purple_BRDF(const Vec3 &position, const Vec3 &incident_light_vector, const Vec3 &outgoing_light_vector) {
	return damping_factor * pi_reciprocal * Radiance(191.0, 64.0, 191.0) / 255.0;
}

Radiance yellow_BRDF(const Vec3 &position, const Vec3 &incident_light_vector, const Vec3 &outgoing_light_vector) {
    return damping_factor * pi_reciprocal * Radiance(1.0, 1.0, 0.0);
}

Radiance brown_BRDF(const Vec3 &position, const Vec3 &incident_light_vector, const Vec3 &outgoing_light_vector) {
    return damping_factor * pi_reciprocal * Radiance(160.0, 82.0, 45.0) / 255.0;
}

Radiance ceiling_light_emitted(const Vec3 &position) {
	return {1.0, 1.0, 1.0};
}

Vec3 rotate_x_axis(const Vec3 &position, const double &angle) {
    return {position.x, position.y * cos(angle) - position.z * sin(angle), position.y * sin(angle) + position.z * cos(angle)};
}

Vec3 rotate_y_axis(const Vec3 &position, const double &angle) {
    return {position.x * cos(angle) + position.z * sin(angle), position.y, -position.x * sin(angle) + position.z * cos(angle)};
}

Vec3 rotate_z_axis(const Vec3 &position, const double &angle) {
    return {position.x * cos(angle) - position.y * sin(angle), position.x * sin(angle) + position.y * cos(angle), position.z};
}


#endif //MAIN_H
