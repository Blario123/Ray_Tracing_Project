#ifndef SDF_H
#define SDF_H

class SDF {
public:
	// Do nothing destructor
	virtual ~SDF() {}
	
	// Pure virtual function for the SDF to be defined in derived classes.
	virtual double SDF_Fct(const Vec3 &position) const = 0;
	
	// The function describing the inverse of the transformation of the position
	// of the current SDF-defined surface. The function is defined via a function
	// pointer.
	virtual Vec3 Inverse_Transformation(const Vec3 &position) const {
		// If there is no inverse transformation defined, don't apply the inverse
		// transformation to the position vector.
		if (Inverse_Transformation_Fct_Pt == 0) {
			return position;
		} else {
			// If there is an inverse transformation defined, apply it to the position
			// vector.
			return Inverse_Transformation_Fct_Pt(position);
		}
	}
	
	// This function returns the components of the light emitted from this
	// object's surface given a position on the surface and the direction of the
	// light.
	virtual Radiance Light_Emitted(const Vec3 &position) {
		if (Light_Emitted_Fct_Pt == 0) {
			// If the pointer to the light emitted is a null pointer, return the zero
			// vector
			Radiance zero_vector(0.0, 0.0, 0.0);
			return zero_vector;
		} else {
			// If Light_Emitted_Fct_Pt points to a function, evaluate this function
			return (*Light_Emitted_Fct_Pt)(position);
		}
	} // End of Light_Emitted
	
	
	// This function returns the Bidirection Reflectance Distribution Function
	// given a position, incident light vector and outgoing light vector.
	virtual Radiance BRDF(const Vec3 &position,
						  const Vec3 &incident_light_vector,
						  const Vec3 &outgoing_light_vector) {
		if (BRDF_Fct_Pt == 0) {
			// If the pointer to the BRDF is a null pointer, return the zero vector
			Radiance zero_vector(0.0, 0.0, 0.0);
			return zero_vector;
		} else {
			// If BRDF_Fct_Pt points to a function, evaluate this function
			return (*BRDF_Fct_Pt)(
					position, incident_light_vector, outgoing_light_vector);
		}
	} // End of BRDF
	
	
	// This function obtains the outward unit normal to the current instance of
	// the SDF-defined defined surface via the central finite difference method.
	virtual Vec3 Outward_Normal(const Vec3 &position,
								const double &finite_difference_size) {
		// The central finite difference of the SDF is taken in all three coordinate
		// directions in order to get the gradient of the SDF. This gradient is a
		// normal to any isosurface at 'position'. This gradient is then normalised
		// to get the unit normal, then the dot product is taken with 'direction'.
		// If the dot product is negative, we flip the signs of the normal and take
		// that as the normal.
		
		// Initialise storage
		Vec3 normal(0.0, 0.0, 0.0);
		
		// Find the non-unit normal, division by step size is unnecessary since the
		// normal will be normalised
		normal.x =
				SDF_Fct(
						Vec3(position.x + finite_difference_size, position.y, position.z)) -
				SDF_Fct(
						Vec3(position.x - finite_difference_size, position.y, position.z));
		
		normal.y =
				SDF_Fct(
						Vec3(position.x, position.y + finite_difference_size, position.z)) -
				SDF_Fct(
						Vec3(position.x, position.y - finite_difference_size, position.z));
		
		normal.z =
				SDF_Fct(
						Vec3(position.x, position.y, position.z + finite_difference_size)) -
				SDF_Fct(
						Vec3(position.x, position.y, position.z - finite_difference_size));
		
		// Find the unit normal
		normal.normalise();
		
		return normal;
	}
	
	// A function pointer to the inverse transformation of the SDF surface
	Vec3 (*Inverse_Transformation_Fct_Pt)(const Vec3 &position) = 0;
	
	// A function pointer to the light emitted
	Radiance (*Light_Emitted_Fct_Pt)(const Vec3 &position) = 0;
	
	// A function pointer to the BRDF
	Radiance (*BRDF_Fct_Pt)(const Vec3 &position,
							const Vec3 &incident_light_vector,
							const Vec3 &outgoing_light_vector) = 0;
	
	// In order to keep consistency in the distance of the global SDF, we must
	// enforce that light rays may only travel in the positive region of an SDF.
	// Therefore, we must be able to choose the inside and outside of an SDF which
	// we do by flipping the sign of an SDF defining a shape.
	bool invert_SDF = false;
};

class SphereSDF : public SDF {
public:
	// SphereSDF constructor: Set the radius and the inside/outside of the sphere
	SphereSDF(const double &radius_, const bool &invert_SDF_) {
		// Set the radius
		radius = radius_;
		
		// The SDF is positive outside of the sphere unless invert_SDF is true.
		invert_SDF = invert_SDF_;
	}
	
	// Define the SDF function of a sphere.
	double SDF_Fct(const Vec3 &position) const {
		// Apply the inverse transformation to the position vector in order to
		// transform the SDF. Calculate the SDF using this inversely transformed
		// position vector.
		double sdf_value = Inverse_Transformation(position).norm() - radius;
		
		// If we want to invert the SDF, we flip the negative region with the
		// positive region. This is done by flipping the sign of the SDF.
		if (invert_SDF) {
			return -sdf_value;
		} else {
			return sdf_value;
		}
	}

private:
	// Radius of the sphere
	double radius;
};

class TorusSDF : public SDF {
public:
  TorusSDF(const double &larger_radius_, const double &smaller_radius_, const bool &invert_SDF_) {
    larger_radius = larger_radius_;
    smaller_radius = smaller_radius_;
    invert_SDF = invert_SDF_;
  }

  double SDF_Fct(const Vec3 &position) const {
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

  double SDF_Fct(const Vec3 &position) const {
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

#endif //SDF_H
