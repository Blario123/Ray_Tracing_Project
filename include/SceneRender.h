#ifndef SCENERENDER_H
#define SCENERENDER_H

#include <chrono>
#include <memory>
#include <thread>

#include "Image.h"
#include "Observer.h"
#include "PhysicalObject.h"
#include "SDF.h"

class SceneRender {
public:
	SceneRender(Observer &observer_) {
		observer_pt = std::make_unique<Observer>(observer_);

		std::random_device random_seed;

		generator = std::default_random_engine(random_seed());
		distribution = std::uniform_real_distribution<double>(0.0, 1.0);
	}
	void Add_Object(std::unique_ptr<PhysicalObject> &object_upt) { object_pt_vector.push_back(std::move(object_upt)); }
	void Add_Object(std::unique_ptr<PhysicalObject> &&object_upt) { object_pt_vector.push_back(std::move(object_upt)); }
	void Add_Object(std::unique_ptr<SDF> &object_upt) { sdf_object_pt_vector.push_back(std::move(object_upt)); }
	void Add_Object(std::unique_ptr<SDF> &&object_upt) { sdf_object_pt_vector.push_back(std::move(object_upt)); }
	Vec3 Hemisphere_Vector_Generator(const Vec3 &normal) {
		Vec3 random_vector;
		while (true) {
			random_vector.x = 2.0 * distribution(generator) - 1.0;
			random_vector.y = 2.0 * distribution(generator) - 1.0;
			random_vector.z = 2.0 * distribution(generator) - 1.0;
			double modulus = random_vector.norm2();
			if (modulus <= 1.0 && modulus > 0.0) {
				random_vector.normalise();
				if (dot(random_vector, normal) < 0.0) {
					random_vector = -random_vector;
				}
				break;
			}
		}
		return random_vector;
	}
	double SDF_Fct(const Vec3 &position) const {
		double minimum_distance = 0.0;
		double current_distance = 0.0;
		if (!sdf_object_pt_vector.empty()) {
			minimum_distance = sdf_object_pt_vector[0]->SDF_Fct(position);
		}
		for (int i = 1; i < sdf_object_pt_vector.size(); i++) {
			current_distance = sdf_object_pt_vector[i]->SDF_Fct(position);
			if (current_distance < minimum_distance) {
				minimum_distance = current_distance;
			}
		}
		return minimum_distance;
	}
	double SDF_Fct(const Vec3 &position, int &sdf_intersection_index) const {
		double minimum_distance = 0.0;
		double current_distance = 0.0;
		if (!sdf_object_pt_vector.empty()) {
			minimum_distance = sdf_object_pt_vector[0]->SDF_Fct(position);
			sdf_intersection_index = 0;
		}
		for (int i = 1; i < sdf_object_pt_vector.size(); i++) {
			current_distance = sdf_object_pt_vector[i]->SDF_Fct(position);
			if (current_distance < minimum_distance) {
				minimum_distance = current_distance;
				sdf_intersection_index = i;
			}
		}
		return minimum_distance;
	}
	Vec3 First_Intersection_Point(const Ray &light_ray, int &object_intersection_index, int &sdf_intersection_index, const double &threshold_intersection_distance, const double &threshold_no_intersection_distance) {
		object_intersection_index = -1;
		sdf_intersection_index = -1;
		double safe_travel_distance = std::abs(SDF_Fct(light_ray.Get_Initial_Position()));
		double new_safe_travel_distance = safe_travel_distance;
		double smallest_distance = 0.0;
		double current_distance = 0.0;
		Vec3 intersection_point(0.0, 0.0, 0.0);
		for (int i = 0; i < object_pt_vector.size(); i++) {
			if (object_intersection_index != -1) {
				if (object_pt_vector[i]->Intersection_Check(light_ray,current_distance) &&
					current_distance < smallest_distance) {
					smallest_distance = current_distance;
					object_intersection_index = i;
				}
			} else {
				if (object_pt_vector[i]->Intersection_Check(light_ray,current_distance)) {
					smallest_distance = current_distance;
					object_intersection_index = i;
				}
			}
		}
		double distance_travelled = 0.0;
		if (!sdf_object_pt_vector.empty() && (safe_travel_distance < smallest_distance || object_intersection_index == -1)) {
			while (safe_travel_distance >= threshold_intersection_distance && safe_travel_distance <= threshold_no_intersection_distance && (distance_travelled < smallest_distance || object_intersection_index == -1)) {
				new_safe_travel_distance = std::abs(SDF_Fct(light_ray.Get_Initial_Position() + (distance_travelled + safe_travel_distance) * light_ray.Get_Direction_Vector()));
				distance_travelled += safe_travel_distance;
				safe_travel_distance = new_safe_travel_distance;
			}
			if (safe_travel_distance < threshold_intersection_distance) {
				intersection_point = light_ray.Get_Initial_Position() + distance_travelled * light_ray.Get_Direction_Vector();
				SDF_Fct(intersection_point, sdf_intersection_index);
				object_intersection_index = -1;
				return intersection_point;
			}
		}
		if (object_intersection_index != -1) {
			intersection_point = (light_ray.Get_Initial_Position() + smallest_distance * light_ray.Get_Direction_Vector());
			return intersection_point;
		}
		return {0.0, 0.0, 0.0};
	}
	Radiance Light_Out(const Ray &light_ray, unsigned bounces_remaining, const double &threshold_intersection_distance, const double &threshold_no_intersection_distance, const double &finite_difference_size) {
		if (bounces_remaining == 0) {
			return {0.0, 0.0, 0.0};
		}
		Radiance resulting_light(0.0, 0.0, 0.0);
		int object_intersection_index = -1;
		int sdf_intersection_index = -1;
		Vec3 normal(0.0, 0.0, 0.0);
		Radiance BRDF(0.0, 0.0, 0.0);
		Vec3 intersection_point = First_Intersection_Point(light_ray,object_intersection_index,sdf_intersection_index, threshold_intersection_distance, threshold_no_intersection_distance);
		if (object_intersection_index == -1 && sdf_intersection_index == -1) {
			return {0.0, 0.0, 0.0};
		}
		if (object_intersection_index == -1) {
			resulting_light = sdf_object_pt_vector[sdf_intersection_index]->Light_Emitted(intersection_point);
			normal = sdf_object_pt_vector[sdf_intersection_index]->Outward_Normal(intersection_point, finite_difference_size);
		} else {
			resulting_light = object_pt_vector[object_intersection_index]->Light_Emitted(intersection_point);
			normal = object_pt_vector[object_intersection_index]->Orientated_Normal(intersection_point, -light_ray.Get_Direction_Vector());
		}
		Vec3 new_direction = Hemisphere_Vector_Generator(normal);
		Vec3 new_position = intersection_point + 1.0e-8 * normal;
		Ray new_ray(new_position, new_direction);
		if (object_intersection_index == -1) {
			BRDF = sdf_object_pt_vector[sdf_intersection_index]->BRDF(intersection_point,light_ray.Get_Direction_Vector(),new_ray.Get_Direction_Vector());
		} else {
			BRDF = object_pt_vector[object_intersection_index]->BRDF(intersection_point,light_ray.Get_Direction_Vector(),new_ray.Get_Direction_Vector());
		}
		resulting_light += 2.0 * pi * BRDF * dot(new_direction, normal) * Light_Out(new_ray,bounces_remaining - 1,threshold_intersection_distance,threshold_no_intersection_distance,finite_difference_size);
		return resulting_light;
	}

	Image Render_Image(const unsigned &number_of_bounces, const unsigned &number_of_random_samples, const double &threshold_intersection_distance, const double &threshold_no_intersection_distance, const double &finite_difference_size, const unsigned &supersampling_value, const bool &silent = false) {
		std::vector<int> resolution = observer_pt->Get_Resolution();
		Image image(resolution[0], resolution[1]);
        Radiance pixel_radiance(0.0, 0.0, 0.0);
		unsigned tenth_percentiles = 0;
		double proportion_done = 0.0;
		for (int i = 0; i < resolution[0]; i++) {
			if (!silent) {
				proportion_done = double(i + 1) / double(resolution[0]);
				if (unsigned(10.0 * proportion_done) > tenth_percentiles) {
					tenth_percentiles = unsigned(10.0 * proportion_done);
					if (tenth_percentiles == 10) {
						std::cout << "The image has been rendered." << std::endl;
					} else {
                        std::cout << "Roughly " << 10 * tenth_percentiles << "\% of the pixels have been rendered." << std::endl;
					}
				}
			}
			for (int j = 0; j < resolution[1]; j++) {
			    pixel_radiance = Radiance(0.0, 0.0, 0.0);
                for(unsigned k = 0; k < number_of_random_samples / (supersampling_value * supersampling_value); k++) {
                    for(unsigned l = 0; l < supersampling_value * supersampling_value; l++) {
                        pixel_radiance += Light_Out(observer_pt->Ray_To_Pixel_XY(supersampling_value * i + l % supersampling_value, supersampling_value * j + l / supersampling_value, k / (supersampling_value * supersampling_value), number_of_random_samples / (supersampling_value * supersampling_value), supersampling_value), number_of_bounces, threshold_intersection_distance, threshold_no_intersection_distance, finite_difference_size);
                    }
                }
                pixel_radiance /= supersampling_value * supersampling_value;
                image(i, j) = pixel_radiance / number_of_random_samples;
            }
		}
		return image;
	}
	void Render_Image_Per_Thread(std::vector<Vec3> &partition, const unsigned &number_of_bounces, const unsigned &number_of_random_samples, const double &threshold_intersection_distance, const double &threshold_no_intersection_distance, const double &finite_difference_size, const unsigned &supersampling_value, const unsigned &thread_index, const unsigned &number_of_threads) {
		// Get the resolution of the image
		std::vector<int> resolution = observer_pt->Get_Resolution();
		unsigned number_of_pixels = resolution[0] * resolution[1];
		partition.resize(0);
		partition.reserve(number_of_pixels / number_of_threads);
		unsigned pixel_index = thread_index;
		Radiance pixel_radiance(0.0, 0.0, 0.0);
		unsigned pixel_index_i = 0;
		unsigned pixel_index_j = 0;
        unsigned progress = 0;
        auto start = std::chrono::steady_clock::now();
        auto previous = std::chrono::steady_clock::now();
        auto current = std::chrono::steady_clock::now();
        double slowest_interval = 0.0;
        double fastest_interval = 0.0;
        double current_interval = 0.0;
        double running_time = 0.0;
		while (pixel_index < number_of_pixels) {
            if(thread_index == 0) {
                if(pixel_index * 10 / number_of_pixels > progress) {
                    progress = pixel_index * 10 / number_of_pixels;
                    current = std::chrono::steady_clock::now();
                    current_interval = std::chrono::duration<double> (current - previous).count();
                    if(progress == 1) {
                        slowest_interval = current_interval;
                        fastest_interval = current_interval;
                    } else if (current_interval > slowest_interval) {
                        slowest_interval = current_interval;
                    } else if (current_interval < fastest_interval) {
                        fastest_interval = current_interval;
                    }
                    previous = current;
                    std::cout << "ETA: " << unsigned(fastest_interval * (10 - progress)) << "-" << unsigned(slowest_interval * double(10 - progress)) << " seconds." << std::endl;
                    std::cout << "The first thread is " << progress << "0\% done." << std::endl << std::endl;
                }
            }
			pixel_index_i = pixel_index % resolution[0];
			pixel_index_j = pixel_index / resolution[0];
			pixel_radiance = Radiance(0.0, 0.0, 0.0);
			for(unsigned i = 0; i < number_of_random_samples / (supersampling_value * supersampling_value); i++) {
                for(unsigned k = 0; k < supersampling_value * supersampling_value; k++) {
                    pixel_radiance += Light_Out(observer_pt->Ray_To_Pixel_XY(supersampling_value * pixel_index_i + k % supersampling_value, supersampling_value * pixel_index_j + k / supersampling_value, i / (supersampling_value * supersampling_value), number_of_random_samples / (supersampling_value * supersampling_value), supersampling_value),number_of_bounces, threshold_intersection_distance, threshold_no_intersection_distance,finite_difference_size);
                }
            }
            pixel_radiance /= supersampling_value * supersampling_value;
			partition.push_back(pixel_radiance / number_of_random_samples);
			pixel_index += number_of_threads;
		}
        if(thread_index == 0) {
            running_time = std::chrono::duration<double> (current - start).count();
            std::cout << "Finished in " << unsigned(running_time) << " seconds." << std::endl;
        }
	}
	Image Render_Image_Multithreaded(const unsigned &number_of_bounces, const unsigned &number_of_random_samples, const double &threshold_intersection_distance, const double &threshold_no_intersection_distance, const double &finite_difference_size, const unsigned &supersampling_value, const unsigned &number_of_threads = std::thread::hardware_concurrency()) {
		std::vector<std::thread> threads;
		std::vector<int> resolution = observer_pt->Get_Resolution();
		Image image(resolution[0], resolution[1]);
		std::vector<std::vector<Vec3>> partitions(number_of_threads);
		for (unsigned thread_index = 0; thread_index < number_of_threads;
			 thread_index++) {
			threads.emplace_back(&SceneRender::Render_Image_Per_Thread,this,std::ref(partitions[thread_index]),number_of_bounces,number_of_random_samples,threshold_intersection_distance,threshold_no_intersection_distance,finite_difference_size,supersampling_value,thread_index,number_of_threads);
		}
		for (std::thread &t: threads) {
			if (t.joinable()) {
				t.join();
			}
		}
		int pixel_index = 0;
		int pixel_index_i = 0;
		int pixel_index_j = 0;
		while (pixel_index < resolution[0] * resolution[1]) {
			pixel_index_i = pixel_index % resolution[0];
			pixel_index_j = pixel_index / resolution[0];
			image(pixel_index_i, pixel_index_j) = partitions[pixel_index % number_of_threads][pixel_index / number_of_threads];
			pixel_index += 1;
		}
		return image;
	}
	Radiance Light_Out_Russian(const Ray &light_ray, const double &threshold_intersection_distance, const double &threshold_no_intersection_distance, const double &finite_difference_size) {
		Radiance resulting_light(0.0, 0.0, 0.0);
		int object_intersection_index = -1;
		int sdf_intersection_index = -1;
		Vec3 normal(0.0, 0.0, 0.0);
		Radiance BRDF(0.0, 0.0, 0.0);
		Vec3 intersection_point = First_Intersection_Point(light_ray,object_intersection_index,sdf_intersection_index,threshold_intersection_distance,threshold_no_intersection_distance);
		if (object_intersection_index == -1 && sdf_intersection_index == -1) {
			return {0.0, 0.0, 0.0};
		}
		if (object_intersection_index == -1) {
			resulting_light = sdf_object_pt_vector[sdf_intersection_index]->Light_Emitted(intersection_point);
			normal = sdf_object_pt_vector[sdf_intersection_index]->Outward_Normal(intersection_point, finite_difference_size);
		} else {
			resulting_light = object_pt_vector[object_intersection_index]->Light_Emitted(intersection_point);
			normal = object_pt_vector[object_intersection_index]->Orientated_Normal(intersection_point, -light_ray.Get_Direction_Vector());
		}
		Vec3 new_direction = Hemisphere_Vector_Generator(normal);
		Vec3 new_position = intersection_point + 1.0e-8 * normal;
		Ray new_ray(new_position, new_direction);
		if (object_intersection_index == -1) {
            BRDF = sdf_object_pt_vector[sdf_intersection_index]->BRDF(intersection_point,light_ray.Get_Direction_Vector(),new_ray.Get_Direction_Vector());
		} else {
			BRDF = object_pt_vector[object_intersection_index]->BRDF(intersection_point,light_ray.Get_Direction_Vector(),new_ray.Get_Direction_Vector());
		}
		double dot_product = dot(new_direction, normal);
		if (distribution(generator) < dot_product) {
			resulting_light += 2.0 * pi * BRDF * Light_Out_Russian(new_ray, threshold_intersection_distance, threshold_no_intersection_distance, finite_difference_size);
		}
		return resulting_light;
	}
	
	Image Render_Image_Russian(const unsigned &number_of_random_samples, const double &threshold_intersection_distance, const double &threshold_no_intersection_distance, const double &finite_difference_size, const unsigned &supersampling_value, const bool &silent = false) {
		std::vector<int> resolution = observer_pt->Get_Resolution();
		Image image(resolution[0], resolution[1]);
        Radiance pixel_radiance(0.0, 0.0, 0.0);
		unsigned tenth_percentiles = 0;
		double proportion_done = 0.0;
		for (int i = 0; i < resolution[0]; i++) {
			if (!silent) {
				proportion_done = double(i + 1) / double(resolution[0]);
				if (unsigned(10.0 * proportion_done) > tenth_percentiles) {
					tenth_percentiles = unsigned(10.0 * proportion_done);
					if (tenth_percentiles == 10) {
						std::cout << "The image has been rendered." << std::endl;
					} else {
						std::cout << "Roughly " << 10 * tenth_percentiles << "\% of the pixels have been rendered." << std::endl;
					}
				}
			}
			for (int j = 0; j < resolution[1]; j++) {
			    pixel_radiance = Radiance(0.0, 0.0, 0.0);
                for(int k = 0; k < number_of_random_samples / (supersampling_value * supersampling_value); k++) {
                    for (int l = 0; l < supersampling_value * supersampling_value; l++) {
                        pixel_radiance += Light_Out_Russian(
                                observer_pt->Ray_To_Pixel_XY(supersampling_value * i + l % supersampling_value, supersampling_value * j + l / supersampling_value, k / (supersampling_value * supersampling_value), number_of_random_samples / (supersampling_value * supersampling_value), supersampling_value),
                                threshold_intersection_distance, threshold_no_intersection_distance,
                                finite_difference_size);
                    }
                }
                pixel_radiance /= supersampling_value * supersampling_value;
                image(i, j) = pixel_radiance / number_of_random_samples;
            }
		}
		return image;
	}
	void Render_Image_Russian_Per_Thread(std::vector<Vec3> &partition, const int &number_of_random_samples, const double &threshold_intersection_distance, const double &threshold_no_intersection_distance, const double &finite_difference_size, const int &supersampling_value, const int &thread_index, const int &number_of_threads) {
		std::vector<int> resolution = observer_pt->Get_Resolution();
		int number_of_pixels = resolution[0] * resolution[1];
		partition.resize(0);
		partition.reserve(number_of_pixels / number_of_threads);
		int pixel_index = thread_index;
		Radiance pixel_radiance(0.0, 0.0, 0.0);
		int pixel_index_i = 0;
		int pixel_index_j = 0;
        unsigned progress = 0;
        auto start = std::chrono::steady_clock::now();
        auto previous = std::chrono::steady_clock::now();
        auto current = std::chrono::steady_clock::now();
        double slowest_interval = 0.0;
        double fastest_interval = 0.0;
        double current_interval = 0.0;
        double running_time = 0.0;
		while (pixel_index < number_of_pixels) {
            if (thread_index == 0) {
                if (pixel_index * 10 / number_of_pixels > progress) {
                    progress = pixel_index * 10 / number_of_pixels;
                    current = std::chrono::steady_clock::now();
                    current_interval = std::chrono::duration<double>(current - previous).count();
                    if (progress == 1) {
                        slowest_interval = current_interval;
                        fastest_interval = current_interval;
                    } else if (current_interval > slowest_interval) {
                        slowest_interval = current_interval;
                    } else if (current_interval < fastest_interval) {
                        fastest_interval = current_interval;
                    }
                    previous = current;
                    std::cout << "ETA: " << unsigned(fastest_interval * (10 - progress)) << "-" << unsigned(slowest_interval * double(10 - progress)) << " seconds." << std::endl;
                    std::cout << "The first thread is " << progress << "0\% done." << std::endl;
                }
            }
			pixel_index_i = pixel_index % resolution[0];
			pixel_index_j = pixel_index / resolution[0];
			pixel_radiance.clear();
            for(int i = 0; i < number_of_random_samples / (supersampling_value * supersampling_value); i++) {
                for (int k = 0; k < supersampling_value * supersampling_value; k++) {
                    pixel_radiance += Light_Out_Russian(
                            observer_pt->Ray_To_Pixel_XY(supersampling_value * pixel_index_i + k % supersampling_value, supersampling_value * pixel_index_j + k / supersampling_value, i / (supersampling_value * supersampling_value), number_of_random_samples / (supersampling_value * supersampling_value), supersampling_value), threshold_intersection_distance, threshold_no_intersection_distance, finite_difference_size);
                }
            }
			partition.push_back(pixel_radiance / number_of_random_samples);
			pixel_index += number_of_threads;
		}
        if (thread_index == 0) {
            running_time = std::chrono::duration<double>(current - start).count();
            std::cout << "Finished in " << unsigned(running_time) << " seconds." << std::endl;
        }
	}
	Image Render_Image_Russian_Multithreaded(const unsigned &number_of_random_samples, const double &threshold_intersection_distance, const double &threshold_no_intersection_distance, const double &finite_difference_size, const unsigned &supersampling_value, const unsigned &number_of_threads = std::thread::hardware_concurrency()) {
		std::vector<std::thread> threads;
		std::vector<int> resolution = observer_pt->Get_Resolution();
		Image image(resolution[0], resolution[1]);
		std::vector<std::vector<Vec3>> partitions(number_of_threads);
		for (unsigned thread_index = 0; thread_index < number_of_threads; thread_index++) {
			threads.emplace_back(&SceneRender::Render_Image_Russian_Per_Thread, this, std::ref(partitions[thread_index]), number_of_random_samples, threshold_intersection_distance, threshold_no_intersection_distance, finite_difference_size, supersampling_value, thread_index, number_of_threads);
		}
		for (std::thread &t: threads) {
			if (t.joinable()) {
				t.join();
			}
		}
		int pixel_index = 0;
		int pixel_index_i = 0;
		int pixel_index_j = 0;
		while (pixel_index < resolution[0] * resolution[1]) {
			pixel_index_i = pixel_index % resolution[0];
			pixel_index_j = pixel_index / resolution[0];
			image(pixel_index_i, pixel_index_j) = partitions[pixel_index % number_of_threads][pixel_index / number_of_threads];
			pixel_index += 1;
		}
		return image;
	}
private:
	std::vector<std::unique_ptr<PhysicalObject>> object_pt_vector;
	std::vector<std::unique_ptr<SDF>> sdf_object_pt_vector;
	std::unique_ptr<Observer> observer_pt;
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution;
};

#endif //SCENERENDER_H
