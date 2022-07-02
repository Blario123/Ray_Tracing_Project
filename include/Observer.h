#ifndef OBSERVER_H
#define OBSERVER_H

#include <random>

#include "Ray.h"
#include "RayTracingProject.h"

class Observer {
public:
	Observer(const Vec3 &position_, const Vec3 &camera_direction_, const Vec3 &upward_direction_, const double &horizontal_field_of_view_angle_, const double &focal_distance_, const double &lens_radius_, const std::vector<int> &resolution_) :
            position(position_),
            unit_camera_direction(camera_direction_ / camera_direction_.norm()),
            unit_upward_direction(upward_direction_ / upward_direction_.norm()),
            unit_sideways_direction(cross(unit_camera_direction, unit_upward_direction)),
            horizontal_field_of_view_angle(horizontal_field_of_view_angle_),
            resolution(resolution_),
            vertical_field_of_view_angle(2.0 *  atan((double(resolution_[1]) / double (resolution_[0])) * tan(horizontal_field_of_view_angle / 2.0))),
            focal_distance(focal_distance_),
            lens_radius(lens_radius_) {
        std::random_device random_seed;

        generator = std::default_random_engine(random_seed());
        distribution = std::uniform_real_distribution<double>(0.0, 1.0);
	}
	Ray Ray_To_Pixel_XY(const unsigned &x_pixel, const unsigned &y_pixel, const unsigned &isample, const unsigned &number_of_random_samples, const unsigned &supersampling_value) {
        std::vector<double> random_coordinates = rnd_unit_disk(isample, number_of_random_samples);
        std::vector<unsigned> supersampling_resolution{resolution[0] * supersampling_value, resolution[1] * supersampling_value};
        Vec3 pseudo_position = position + lens_radius * (random_coordinates[0] * unit_sideways_direction + random_coordinates[1] * unit_upward_direction);
        Vec3 direction_to_pixel = unit_camera_direction;
        if (supersampling_resolution[0] != 1) {
            direction_to_pixel += tan(horizontal_field_of_view_angle / 2.0) * ((2.0 * double(x_pixel)) / (double(supersampling_resolution[0]) - 1.0) - 1.0) * unit_sideways_direction;
        }
        if (supersampling_resolution[1] != 1) {
            direction_to_pixel -= tan(vertical_field_of_view_angle / 2.0) * ((2.0 * double(y_pixel)) / (double(supersampling_resolution[1]) - 1.0) - 1.0) * unit_upward_direction;
        }
		direction_to_pixel.normalise();
		Vec3 direction_to_pixel_new = (position - pseudo_position) + (focal_distance / dot(direction_to_pixel, unit_camera_direction)) * direction_to_pixel;
        direction_to_pixel_new.normalise();
        return {pseudo_position, direction_to_pixel_new};
	}
	std::vector<int> Get_Resolution() const {
        return resolution;
    }
    virtual std::vector<double> rnd_unit_disk(const unsigned &isample, const unsigned &number_of_samples) {
        std::vector<double> coordinates{1.0, 1.0};
        while(coordinates[0] * coordinates[0] + coordinates[1] * coordinates[1] > 1.0) {
            coordinates[0] = 2.0 * distribution(generator) - 1.0;
            coordinates[1] = 2.0 * distribution(generator) - 1.0;
        }
        return coordinates;
    }
	Vec3 position;
	Vec3 unit_camera_direction;
	Vec3 unit_upward_direction;
	Vec3 unit_sideways_direction;
	double horizontal_field_of_view_angle;
	double vertical_field_of_view_angle;
	std::vector<int> resolution;
    double focal_distance;
    double lens_radius;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution;
};

class StratifiedObserver : public Observer {
public:
    StratifiedObserver(const Vec3 &position_, const Vec3 &camera_direction_, const Vec3 &upward_direction_, const double &horizontal_field_of_view_angle_, const double &focal_distance_, const double &lens_radius_, const std::vector<int> &resolution_) : Observer(position_, camera_direction_, upward_direction_, horizontal_field_of_view_angle_, focal_distance_, lens_radius_, resolution_) {

	}
	std::vector<double> rnd_unit_disk(const unsigned &isample, const unsigned &number_of_random_samples) {
		unsigned nside = int(sqrt(number_of_random_samples));
        double x = 2.0 * ((isample % nside) + distribution(generator)) / nside - 1.0;
        double y = 2.0 * (((isample / nside) % nside) + distribution(generator)) / nside - 1.0;
        double r = 0.0;
        double theta = 1.0;
		if (std::abs(x) > std::abs(y)) {
			r = x;
            theta = pi * (y / x) / 4.0;
		} else {
            r = y;
            theta = pi / 2.0 - (pi / 4.0) * (x / y);
		}
		return std::vector<double> {r * cos(theta), r * sin(theta)};
	}
};

#endif //OBSERVER_H
