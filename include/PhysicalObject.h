#ifndef PHYSICALOBJECT_H
#define PHYSICALOBJECT_H

#include "Ray.h"
#include "RayTracingProject.h"

class PhysicalObject {
public:
	virtual ~PhysicalObject() = default;
	virtual Radiance Light_Emitted(const Vec3 &position) {
        if (Light_Emitted_Fct_Pt == 0) {
            Radiance zero_vector(0.0, 0.0, 0.0);
            return zero_vector;
        } else {
            return (*Light_Emitted_Fct_Pt)(position);
        }
    }
	virtual Radiance BRDF(const Vec3 &position, const Vec3 &incident_light_vector, const Vec3 &outgoing_light_vector) {
        if (BRDF_Fct_Pt == 0) {
            Radiance zero_vector(0.0, 0.0, 0.0);
            return zero_vector;
        } else {
            return (*BRDF_Fct_Pt)(
                    position, incident_light_vector, outgoing_light_vector);
        }
    }
	virtual bool Intersection_Check(const Ray &incident_ray, double &distance) = 0;
	virtual Vec3 Orientated_Normal(const Vec3 &position, const Vec3 &direction) = 0;
	Radiance (*Light_Emitted_Fct_Pt)(const Vec3 &position) = nullptr;
	Radiance (*BRDF_Fct_Pt)(const Vec3 &position, const Vec3 &incident_light_vector, const Vec3 &outgoing_light_vector) = nullptr;
};

class Sphere : public PhysicalObject {
public:
	Sphere(const Vec3 &centre_, const double &radius_) : centre(centre_), radius(radius_) {

    }
	bool Intersection_Check(const Ray &light_ray, double &distance) override {
        Vec3 initial_position = light_ray.Get_Initial_Position();
        Vec3 direction_vector = light_ray.Get_Direction_Vector();
        double quadratic_coefficient = direction_vector.norm2();
        double linear_coefficient = 2.0 * (dot(initial_position, direction_vector) - dot(centre, direction_vector));
        double constant =
                initial_position.norm2() + centre.norm2() - 2.0 * dot(initial_position, centre) - radius * radius;
        double determinant = linear_coefficient * linear_coefficient - 4.0 * quadratic_coefficient * constant;
        if (determinant > 0.0) {
            const bool linear_coefficient_bigger_than_sqrt_determinant =
                    std::abs(linear_coefficient) > sqrt(determinant);
            if (std::abs(linear_coefficient) < 1.0e-8) {
                distance = sqrt(determinant) / (2.0 * quadratic_coefficient);
                return true;
            } else if (linear_coefficient > 0.0) {
                if (linear_coefficient_bigger_than_sqrt_determinant) {
                    return false;
                } else {
                    distance = (-linear_coefficient + sqrt(determinant)) / (2.0 * quadratic_coefficient);
                    return true;
                }
            } else {
                if (linear_coefficient_bigger_than_sqrt_determinant) {
                    distance = (-linear_coefficient - sqrt(determinant)) / (2.0 * quadratic_coefficient);
                    return true;
                } else {
                    distance = (-linear_coefficient + sqrt(determinant)) / (2.0 * quadratic_coefficient);
                    return true;
                }
            }
        }
        return false;
    }
	Vec3 Orientated_Normal(const Vec3 &position, const Vec3 &direction) override {
		Vec3 normal = position - centre;
		normal.normalise();
		if (dot(normal, direction) < 0.0) {
			normal = -normal;
		}
		return normal;
	}
private:
	Vec3 centre;
	double radius;
};

#endif //PHYSICALOBJECT_H
