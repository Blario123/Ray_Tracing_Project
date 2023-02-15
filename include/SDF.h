#ifndef SDF_H
#define SDF_H

#include "Vec3.h"

class SDF {
public:
	virtual ~SDF() {}
	virtual double SDF_Fct(const Vec3 &position) const = 0;
	virtual Vec3 Inverse_Transformation(const Vec3 &position) const {
		if (Inverse_Transformation_Fct_Pt == 0) {
			return position;
		} else {
			return Inverse_Transformation_Fct_Pt(position);
		}
	}
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
			return (*BRDF_Fct_Pt)(position, incident_light_vector, outgoing_light_vector);
		}
	}
	virtual Vec3 Outward_Normal(const Vec3 &position, const double &finite_difference_size) {
		Vec3 normal(0.0, 0.0, 0.0);
		normal.x = SDF_Fct(Vec3(position.x + finite_difference_size, position.y, position.z)) - SDF_Fct(Vec3(position.x - finite_difference_size, position.y, position.z));
		normal.y = SDF_Fct(Vec3(position.x, position.y + finite_difference_size, position.z)) - SDF_Fct(Vec3(position.x, position.y - finite_difference_size, position.z));
		normal.z = SDF_Fct(Vec3(position.x, position.y, position.z + finite_difference_size)) - SDF_Fct(Vec3(position.x, position.y, position.z - finite_difference_size));
		normal.normalise();
		return normal;
	}
	Vec3 (*Inverse_Transformation_Fct_Pt)(const Vec3 &position) = nullptr;
	Radiance (*Light_Emitted_Fct_Pt)(const Vec3 &position) = nullptr;
	Radiance (*BRDF_Fct_Pt)(const Vec3 &position, const Vec3 &incident_light_vector, const Vec3 &outgoing_light_vector) = nullptr;
	bool invert_SDF = false;
};

class SphereSDF : public SDF {
public:
	SphereSDF(const double &radius_, const bool &invert_SDF_) {
		radius = radius_;
		invert_SDF = invert_SDF_;
	}
	double SDF_Fct(const Vec3 &position) const override {
		double sdf_value = Inverse_Transformation(position).norm() - radius;
		if (invert_SDF) {
			return -sdf_value;
		} else {
			return sdf_value;
		}
	}
private:
	double radius;
};

class TorusSDF : public SDF {
public:
    TorusSDF(const double &larger_radius_, const double &smaller_radius_, const bool &invert_SDF_) : larger_radius(larger_radius_), smaller_radius(smaller_radius_) {
        invert_SDF = invert_SDF_;
    }
    double SDF_Fct(const Vec3 &position) const override {
        Vec3 transformed_position = Inverse_Transformation(position);
        double distance_to_inner_circle = sqrt(pow(transformed_position.x, 2) + pow(transformed_position.y, 2)) - larger_radius;
        double sdf_value = sqrt(pow(distance_to_inner_circle, 2) + pow(transformed_position.z, 2)) - smaller_radius;
        if (invert_SDF) {
            return -sdf_value;
        } else {
            return sdf_value;
        }
    }
private:
    double larger_radius;
    double smaller_radius;
};

class CuboidSDF : public SDF {
public:
    CuboidSDF(const Vec3 &vertex_, const bool &invert_SDF_) {
        vertex.x = std::abs(vertex_.x);
        vertex.y = std::abs(vertex_.y);
        vertex.z = std::abs(vertex_.z);
        invert_SDF = invert_SDF_;
    }
    double SDF_Fct(const Vec3 &position) const override {
        Vec3 transformed_position = Inverse_Transformation(position);
        Vec3 absolute_position(std::abs(transformed_position.x),std::abs(transformed_position.y),std::abs(transformed_position.z));
        Vec3 vertex_to_absolute_position = absolute_position - vertex;
        double first_term = Vec3(std::max(vertex_to_absolute_position.x, 0.0), std::max(vertex_to_absolute_position.y, 0.0), std::max(vertex_to_absolute_position.z, 0.0)).norm();
        double second_term = std::min(std::max(std::max(vertex_to_absolute_position.x, vertex_to_absolute_position.y), vertex_to_absolute_position.z), 0.0);
        double sdf_value = first_term + second_term;
        if (invert_SDF) {
            return -sdf_value;
        } else {
            return sdf_value;
        }
    }
private:
  Vec3 vertex;
};

class ReuleauxSDF : public SDF {
public:
    ReuleauxSDF(const double &side_length_, const bool &invert_SDF_) : side_length(side_length_) {
        invert_SDF = invert_SDF_;
    }
    double SDF_Fct(const Vec3 &position) const override {
        Vec3 transformed_position = Inverse_Transformation(position);
        double sdf_1 = (transformed_position - side_length * Vec3(0.5, 0.0, -0.5 / sqrt(2.0))).norm() - side_length;
        double sdf_2 = (transformed_position - side_length * Vec3(-0.5, 0.0, -0.5 / sqrt(2.0))).norm() - side_length;
        double sdf_3 =
                (transformed_position - side_length * Vec3(0.0, 0.5, 0.5 / sqrt(2.0)))
                        .norm() -
                side_length;
        double sdf_4 = (transformed_position - side_length * Vec3(0.0, -0.5, 0.5 / sqrt(2.0))).norm() - side_length;
        double sdf_value = std::max(std::max(sdf_1, sdf_2), std::max(sdf_3, sdf_4));
        if (invert_SDF) {
            return -sdf_value;
        } else {
            return sdf_value;
        }
    }
private:
    double side_length;
};

#endif //SDF_H
