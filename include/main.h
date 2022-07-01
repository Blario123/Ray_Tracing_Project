#ifndef MAIN_H
#define MAIN_H

#include "Vec3.h"
#include "RayTracingProject.h"
#include "GUI.h"

Radiance ceiling_light_emitted(const Vec3 &position) {
	return {1.0, 1.0, 1.0};
}

#endif //MAIN_H
