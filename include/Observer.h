#ifndef OBSERVER_H
#define OBSERVER_H

class Observer {
public:
	Observer(const Vec3 &position_,
			 const Vec3 &camera_direction_,
			 const Vec3 &upward_direction_,
			 const double &horizontal_field_of_view_angle_,
			 const std::vector<unsigned> &resolution_) {
		// Define the member variables
		
		// Position of the camera
		position = position_;
		
		// Direction the camera is pointing in
		unit_camera_direction = camera_direction_ / camera_direction_.norm();
		
		// The upwards direction of the camera
		unit_upward_direction = upward_direction_ / upward_direction_.norm();
		
		// The vector pointing "right" from the point of view of the camera
		unit_sideways_direction =
				cross(unit_camera_direction, unit_upward_direction);
		
		// The horizontal viewing angle(radians)
		horizontal_field_of_view_angle = horizontal_field_of_view_angle_;
		
		// The resolution of the camera
		resolution = resolution_;
		
		// The vertical viewing angle
		vertical_field_of_view_angle =
				2.0 * atan((double(resolution[1]) / double(resolution[0])) *
						   tan(horizontal_field_of_view_angle / 2.0));
	}
	
	
	// When given a pixel, this function will return the normalised vector from
	// the camera to the given pixel on the "grid" in front of the camera
	virtual Ray Ray_To_Pixel_XY(const unsigned &x_pixel, const unsigned &y_pixel) {
		// Add the vector from the camera to the centre of the "tennis racket" with
		// the vector to the pixel specified from the centre of the "tennis racket"
		
		Vec3 direction_to_pixel = unit_camera_direction;
		
		// If one of the dimensions of the resolution is 1, the calculation of
		// 'direction_to_pixel' returns nan when we want 0. So we check the
		// dimensions of the resolution in turn.
		if (resolution[0] != 1) {
			direction_to_pixel +=
					tan(horizontal_field_of_view_angle / 2.0) *
					((2.0 * double(x_pixel)) / (double(resolution[0]) - 1.0) - 1.0) *
					unit_sideways_direction;
		}
		
		if (resolution[1] != 1) {
			direction_to_pixel -=
					tan(vertical_field_of_view_angle / 2.0) *
					((2.0 * double(y_pixel)) / (double(resolution[1]) - 1.0) - 1.0) *
					unit_upward_direction;
		}
		
		// Normalise this direction vector
		direction_to_pixel.normalise();
		
		return Ray(position, direction_to_pixel);
	} // End of Ray_To_Pixel_XY
	
	// Constant access function to the resolution of the camera
	std::vector<unsigned> Get_Resolution() const {
		// Just return the vector containing the resolution
		return resolution;
	} // End of Get_Resolution

protected:
	// The position of the observer
	Vec3 position;
	
	// The direction the observer is looking in (unit vector)
	Vec3 unit_camera_direction;
	
	// The vector defining the upward direction (unit vector)
	Vec3 unit_upward_direction;
	
	// The vector resulting from the cross product of the upward direction with
	// the direction of the observer
	Vec3 unit_sideways_direction;
	
	// The horizontal viewing angle
	double horizontal_field_of_view_angle;
	
	// The vertical viewing angle (calculated from the horizontal field of view
	// and the resolution in order to correspond to square pixels)
	double vertical_field_of_view_angle;
	
	// The resolution of the observer's camera (e.g 1920x1080 means that there are
	// 1920 columns and 1080 rows of pixels)
	std::vector<unsigned> resolution;
};

class DOFObserver : public Observer {
public:
	DOFObserver(const Vec3 &position_,
				const Vec3 &camera_direction_,
				const Vec3 &upward_direction_,
				const double &horizontal_field_of_view_angle_,
				const double &focal_distance_,
				const double &lens_radius_,
				const std::vector<unsigned> &resolution_) :
			Observer(position_,
					 camera_direction_,
					 upward_direction_,
					 horizontal_field_of_view_angle_,
					 resolution_) {
		// Define the member variables
		
		// Position of the camera
		position = position_;
		
		// Direction the camera is pointing in
		unit_camera_direction = camera_direction_ / camera_direction_.norm();
		
		// The upwards direction of the camera
		unit_upward_direction = upward_direction_ / upward_direction_.norm();
		
		// The vector pointing "right" from the point of view of the camera
		unit_sideways_direction =
				cross(unit_camera_direction, unit_upward_direction);
		
		// The horizontal viewing angle(radians)
		horizontal_field_of_view_angle = horizontal_field_of_view_angle_;
		
		// The resolution of the camera
		resolution = resolution_;
		
		// The vertical viewing angle
		vertical_field_of_view_angle =
				2.0 * atan((double(resolution[1]) / double(resolution[0])) *
						   tan(horizontal_field_of_view_angle / 2.0));
		
		// The distance from the observer to the plane of focus
		focal_distance = focal_distance_;
		
		// The radius of the thin lens
		lens_radius = lens_radius_;
		
		// Create a uniform random number generator
		// Get a seed with a random value
		std::random_device random_seed;
		
		// Create a uniform distribution between 0 and 1
		generator = std::default_random_engine(random_seed());
		distribution = std::uniform_real_distribution<double>(0.0, 1.0);
	}
	
	
	Ray Ray_To_Pixel_XY(const unsigned &x_pixel, const unsigned &y_pixel) {
		std::vector<double> random_coordinates = rnd_unit_disk();
		
		Vec3 pseudo_position =
				position +
				lens_radius * (random_coordinates[0] * unit_sideways_direction +
							   random_coordinates[1] * unit_upward_direction);
		
		// Add the vector from the camera to the centre of the "tennis racket" with
		// the vector to the pixel specified from the centre of the "tennis racket"
		
		Vec3 direction_to_pixel = unit_camera_direction;
		
		// If one of the dimensions of the resolution is 1, the calculation of
		// 'direction_to_pixel' returns nan when we want 0. So we check the
		// dimensions of the resolution in turn.
		if (resolution[0] != 1) {
			direction_to_pixel +=
					tan(horizontal_field_of_view_angle / 2.0) *
					((2.0 * double(x_pixel)) / (double(resolution[0]) - 1.0) - 1.0) *
					unit_sideways_direction;
		}
		
		if (resolution[1] != 1) {
			direction_to_pixel -=
					tan(vertical_field_of_view_angle / 2.0) *
					((2.0 * double(y_pixel)) / (double(resolution[1]) - 1.0) - 1.0) *
					unit_upward_direction;
		}
		
		// Normalise this direction vector
		direction_to_pixel.normalise();
		
		Vec3 direction_to_pixel_new =
				(position - pseudo_position) +
				(focal_distance / dot(direction_to_pixel, unit_camera_direction)) *
				direction_to_pixel;
		
		direction_to_pixel_new.normalise();
		
		return Ray(pseudo_position, direction_to_pixel_new);
	} // End of Ray_To_Pixel_XY


private:
	// Sample a random point on a unit disk and output the coordinates
	std::vector<double> rnd_unit_disk() {
		// Create a 2D coordinate with a magnitude greater than 1
		std::vector<double> coordinates{1.0, 1.0};
		
		// Continually create random coordinates until the magnitude is smaller than
		// 1
		while (coordinates[0] * coordinates[0] + coordinates[1] * coordinates[1] >
			   1.0) {
			// Set each component to a random number
			coordinates[0] = 2.0 * distribution(generator) - 1.0;
			coordinates[1] = 2.0 * distribution(generator) - 1.0;
		}
		
		// Return the coordinates when a suitable coordinate is found
		return coordinates;
	}
	
	// The distance from the observer to the plane of focus
	double focal_distance;
	
	// The radius of the thin lens
	double lens_radius;
	
	// Random engine and distribution for use in generating random vectors
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution;
};

#endif //OBSERVER_H
